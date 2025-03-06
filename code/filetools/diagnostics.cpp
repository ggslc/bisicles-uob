#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//===========================================================================
// diagnostics.cpp
// read in a bisicles plotfile 
// and an optional mask, and write out a bunch of diagnostics about the ice sheet. 
// These are
// 0. time
// 1. Volume of ice
// 2. Volume of ice above flotation
// dh/dt
// SMB
// BMB
// Volume of calved ice
// Calving flux  
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "LevelSigmaCS.H"
#include "IceConstants.H"
#include "computeSum.H"
#include "CellToEdge.H"
#include "amrIceF_F.H"
#include "AdvectPhysics.H"
#include "PatchGodunov.H"
#include "SigmaCSF_F.H"
#include "PiecewiseLinearFillPatch.H"
#include "FillFromReference.H"
#include "IceThicknessIBC.H"
#include "LevelDataIBC.H"
#include "DomainDiagnosticData.H"
#include <functional>
 

/**
   convert single level mask data, to which stores integer mask labels (as reals...) 
   between a_mask_label_start and a_mask_no, to multi-level real fractional coverage areas for each 
   distinct mask label. Resulting a_mask_fraction FABS will have a_mask_no components
*/
void maskToAMR(Vector<LevelData<FArrayBox>* >& a_mask_fraction,
	       const Vector<Real>& a_mask_fraction_dx,
	       const LevelData<FArrayBox>& a_mask,
	       int a_mask_id_start,
	       int a_n_mask_ids,
	       const Real& a_mask_dx)
{

  
  Real tol = 1.0e-10; // tolerance in Abs(x - y) < tol test
  
  for (int comp = 0; comp < a_n_mask_ids ; comp++)
    {
      Real mask_id = Real(comp + a_mask_id_start);
      
      // create field of mask areas on the mask level (all cells 1 or 0)
      const DisjointBoxLayout mask_grids = a_mask.disjointBoxLayout();
      LevelData<FArrayBox> mask_fraction(mask_grids,1,IntVect::Zero);
      for (DataIterator dit(mask_grids); dit.ok(); ++dit)
	{
	  for (BoxIterator bit(mask_grids[dit]); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      mask_fraction[dit](iv) = (Abs(a_mask[dit](iv) - mask_id) < tol)?1.0:0.0;  
	    }
	}
      // coarsen/refine the mask areas, store in component mask_no - a_mask_no_start)
      for (int lev = 0; lev < a_mask_fraction.size(); lev++)
	{
	  LevelData<FArrayBox>& dest = *a_mask_fraction[lev];
	  LevelData<FArrayBox> tmp(dest.disjointBoxLayout(),1,IntVect::Zero);
	  RealVect src_dx = RealVect::Unit * a_mask_dx;
	  RealVect dest_dx = RealVect::Unit * a_mask_fraction_dx[lev];
	  FillFromReference(tmp, mask_fraction, dest_dx, src_dx, false);
	  tmp.localCopyTo(Interval(0,0), dest, Interval(comp, comp));
	}
    }
}


struct NameUnitValue
{
  std::string name;
  std::string unit;
  Real value;
  NameUnitValue(std::string n, std::string u, Real v)
    :name(n), unit(u), value(v){}

};

std::ostream& json(std::ostream& o, NameUnitValue& a)
  {
    o << "\"" << a.name << "\":{\"unit\":\"" << a.unit << "\", \"value\":" << a.value << "}";
    return o;
  }

std::ostream& operator<<(std::ostream& o, NameUnitValue& a)
  {
    o << a.name << "," << a.unit << "," << a.value;
    return o;
  }


void reportConservationInside(Vector<NameUnitValue>& report,
			      function<bool(Real h, Real f, int mask)> inside,
			      const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
			      const Vector<LevelData<FluxBox>* >& fluxOfIce,
			      const Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
			      const Vector<LevelData<FArrayBox>* >& basalThicknessSource,
			      const Vector<LevelData<FArrayBox>* >& deltaThickness,
			      const Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux,
			      const Vector<LevelData<FArrayBox>* >& calvingFlux,
			      const Vector<LevelData<FArrayBox>* >& topography,
			      const Vector<LevelData<FArrayBox>* >& thickness,
			      const Vector<LevelData<FArrayBox>* >& iceFrac, 
			      const Vector<LevelData<FArrayBox>* >& sectorMaskFraction,
			      const Vector<Real>& dx, const Vector<int>& ratio, 
			      int finestLevel, int maskNo, int maskComp)
{

  CH_TIME("reportConservationInside");
  // just to make the args less reptitive...

  auto sumScalar = [inside,coords,topography, thickness, iceFrac, sectorMaskFraction, dx, ratio, finestLevel, maskNo, maskComp ]
    (const Vector<LevelData<FArrayBox>* >& scalar)
		   {
		     Real sumScalar;
		     MaskedIntegration::integrateScalarInside
		       (sumScalar, inside ,coords, scalar, topography, thickness,
			iceFrac,sectorMaskFraction, dx, ratio, finestLevel,  maskNo, maskComp);
		     return sumScalar;
		   };

  std::string dhunit("m3/a");
  Real SMB = sumScalar(surfaceThicknessSource);
  report.push_back(NameUnitValue("SMB",dhunit,SMB));

  Real BMB = sumScalar(basalThicknessSource);
  report.push_back(NameUnitValue("BMB",dhunit,BMB));
  
  Real dhdt = sumScalar(deltaThickness);
  report.push_back(NameUnitValue("dhdt",dhunit,dhdt));

  Real flxDivFile = sumScalar(divergenceThicknessFlux);
  report.push_back(NameUnitValue("flxDivFile",dhunit,flxDivFile));

  Real calving = sumScalar(calvingFlux);
  report.push_back(NameUnitValue("calving",dhunit,calving));

  //discharge, div(uh) reconstructed from uh
  Real sumDischarge, flxDivReconstr;
  MaskedIntegration::integrateDischargeInside(sumDischarge, flxDivReconstr,
			   inside, coords, fluxOfIce,topography,
			   thickness, iceFrac, dx, ratio, sectorMaskFraction, 
			   finestLevel, maskNo, maskComp);
  report.push_back(NameUnitValue("flxDivReconstr",dhunit,flxDivReconstr));
  report.push_back(NameUnitValue("discharge",dhunit,sumDischarge));

  //Conservation errors
  report.push_back(NameUnitValue("SMB+BMB-dhdt-calving-flxDivFile",dhunit,
				 SMB + BMB - dhdt - calving - flxDivFile));
  report.push_back(NameUnitValue("SMB+BMB-dhdt-calving-discharge",dhunit,
				 SMB + BMB - dhdt - calving - sumDischarge));
		   
  // volumes
  std::string vunit("m3");
  report.push_back(NameUnitValue("volume",vunit,sumScalar(thickness)));
  // volume above flotation
  {
    Vector<LevelData<FArrayBox>* > hab;
    for (int lev = 0; lev < coords.size(); lev++)
      {
	hab.push_back(const_cast<LevelData<FArrayBox>* >(&(coords[lev]->getThicknessOverFlotation())));
      }
    report.push_back(NameUnitValue("volumeAbove",vunit,sumScalar(hab)));
  }

  // Areas
  std::string aunit("m2");
  // ice covered are according to ice fraction
  report.push_back(NameUnitValue("fracArea",aunit,sumScalar(iceFrac)));
  
  // Area of the region specified by inside(h,f,mask) && (sector == maskNo).
  // Maybe put this elsehwere
  {
    Vector<LevelData<FArrayBox>* > one(thickness.size(),NULL);
    for (int lev = 0; lev < one.size(); lev++)
      {
	const DisjointBoxLayout& grids = (*thickness[lev]).disjointBoxLayout();
	one[lev] = new LevelData<FArrayBox>(grids, 1, IntVect::Zero);
	for (DataIterator dit(grids); dit.ok(); ++dit)
	  (*one[lev])[dit].setVal(1.0);
      }
    
    report.push_back(NameUnitValue("area",aunit,sumScalar(one)));
    for (int lev = 0; lev < one.size(); lev++)
      {
	if (one[lev])
	  {
	    delete one[lev]; one[lev] = NULL;
	  }
      } 
  }
}


void createDEM( Vector<LevelData<FArrayBox>* >& topography,  
		Vector<LevelData<FArrayBox>* >& thickness,
		Vector<LevelData<FArrayBox>* >& iceFrac,
		const Vector<std::string>& name, 
		const Vector<LevelData<FArrayBox>* >& data,
		const Vector<int>& ratio,
		const Vector<Real>& dx,
		Real hmin)
{
  CH_TIME("createDEM");
  int numLevels = data.size();

  for (int lev = 0; lev < numLevels; lev++)
    {

      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  (*topography[lev])[dit].setVal(0.0);
	  (*thickness[lev])[dit].setVal(0.0);
	}

      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "Z_base")
	    {
	      data[lev]->copyTo(Interval(j,j),*topography[lev],Interval(0,0));
	    }
	  else if (name[j] == "thickness")
	    {
	      data[lev]->copyTo(Interval(j,j),*thickness[lev],Interval(0,0));
	    }	  
	}

      //estimate frac
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  for (BoxIterator bit(grids[dit]); bit.ok(); ++bit)
	    {
	      (*iceFrac[lev])[dit](bit()) = ((*thickness[lev])[dit](bit()) > hmin)?1.0:0.0;
	    }
	}
      
      // replace estimated frac with actual frac, it it exists
      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "iceFrac")
	    {
	      data[lev]->copyTo(Interval(j,j),*iceFrac[lev],Interval(0,0));
	    }
	}
      

      if (lev > 0)
	{
	  
	  const DisjointBoxLayout& crseGrids = topography[lev-1]->disjointBoxLayout();
	  PiecewiseLinearFillPatch filler(grids , crseGrids, 1, 
					  crseGrids.physDomain(), ratio[lev-1], 2);
	  Real time_interp_coeff = 0.0;
	  filler.fillInterp(*topography[lev],*topography[lev-1] ,*topography[lev-1],
			    time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*thickness[lev],*thickness[lev-1] ,*thickness[lev-1],
			    time_interp_coeff,0, 0, 1);

	}
      thickness[lev] -> exchange();
      topography[lev] -> exchange();
    }
}

void createSigmaCS(Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		   Vector<LevelData<FArrayBox>* >& topography, 
		   Vector<LevelData<FArrayBox>* >& thickness,
		   Vector<LevelData<FArrayBox>* >& iceFrac,
		   Vector<Real>& dx, Vector<int>& ratio,
		   Real iceDensity, Real waterDensity, Real gravity)
{

  int numLevels = topography.size();

  IntVect sigmaCSGhost(2*IntVect::Unit);
       
  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& levelGrids = topography[lev]->disjointBoxLayout();
	   
      coords[lev] = RefCountedPtr<LevelSigmaCS> 
	(new LevelSigmaCS(levelGrids, RealVect::Unit*dx[lev], sigmaCSGhost));
      coords[lev]->setIceDensity(iceDensity);
      coords[lev]->setWaterDensity(waterDensity);
      coords[lev]->setGravity(gravity);
      // FAB - FAB copy copies ghosts 
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  coords[lev]->getTopography()[dit].copy( (*topography[lev])[dit]);
	  coords[lev]->getH()[dit].copy( (*thickness[lev])[dit]);
	}
      LevelSigmaCS* crseCoords = (lev > 0)?&(*coords[lev-1]):NULL;
      coords[lev]->recomputeGeometry(crseCoords, (lev > 0)?ratio[lev-1]:0);
	   
    }
       
}

void extractThicknessSource(Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
			    Vector<LevelData<FArrayBox>* >& basalThicknessSource,
			    Vector<LevelData<FArrayBox>* >& deltaThickness,
			    Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux,
			    Vector<LevelData<FArrayBox>* >& calvingFlux,
			    Vector<LevelData<FArrayBox>* >& melangeThickness,
			    Vector<LevelData<FArrayBox>* >& topography,
			    const Vector<Real>& dx, const Vector<int>& ratio, 
			    const Vector<std::string>& name, 
			    const Vector<LevelData<FArrayBox>* >& data)

{
  int numLevels = data.size();

  for (int lev = 0; lev < numLevels; lev++)
    {

      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
          (*deltaThickness[lev])[dit].setVal(0.0);
          (*surfaceThicknessSource[lev])[dit].setVal(0.0);
          (*basalThicknessSource[lev])[dit].setVal(0.0);
          (*divergenceThicknessFlux[lev])[dit].setVal(0.0);
          (*melangeThickness[lev])[dit].setVal(0.0);
          (*calvingFlux[lev])[dit].setVal(0.0);
	}

      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "dThickness/dt")
	    {
	      data[lev]->copyTo(Interval(j,j),*deltaThickness[lev],Interval(0,0));
	    }
	  else if (name[j] == "activeSurfaceThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*surfaceThicknessSource[lev],Interval(0,0));
	    }	  
	  else if (name[j] == "activeBasalThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*basalThicknessSource[lev],Interval(0,0));
	    }
	  else if (name[j] == "divergenceThicknessFlux")
	    {
	      data[lev]->copyTo(Interval(j,j),*divergenceThicknessFlux[lev],Interval(0,0));
	    }
	  else if (name[j] == "calvedThicknessSource" || name[j] == "calvingFlux")
	    {
	      data[lev]->copyTo(Interval(j,j),*calvingFlux[lev],Interval(0,0));
	    }
	  else if (name[j] == "calvedIceThickness" || name[j] == "melangeThickness")
	    {
	      data[lev]->copyTo(Interval(j,j),*melangeThickness[lev],Interval(0,0));
	    }
	  
	}
     

      if (lev > 0)
	{
	  
	  const DisjointBoxLayout& crseGrids = topography[lev-1]->disjointBoxLayout();
	  PiecewiseLinearFillPatch filler(grids , crseGrids, 1, 
					  crseGrids.physDomain(), ratio[lev-1], 1);
	  Real time_interp_coeff = 0.0;

	  filler.fillInterp(*deltaThickness[lev],*deltaThickness[lev-1],
			    *deltaThickness[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*surfaceThicknessSource[lev],*surfaceThicknessSource[lev-1],
			    *surfaceThicknessSource[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*basalThicknessSource[lev],*basalThicknessSource[lev-1],
			    *basalThicknessSource[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*divergenceThicknessFlux[lev],*divergenceThicknessFlux[lev-1],
			    *divergenceThicknessFlux[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*calvingFlux[lev],*calvingFlux[lev-1],
			    *calvingFlux[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*melangeThickness[lev],*melangeThickness[lev-1],
			    *melangeThickness[lev-1],time_interp_coeff,0, 0, 1);
	}
      deltaThickness[lev]->exchange();
      surfaceThicknessSource[lev]->exchange();
      basalThicknessSource[lev]->exchange();
      divergenceThicknessFlux[lev]->exchange();
      calvingFlux[lev]->exchange();
      melangeThickness[lev]->exchange();
    }
} 

void computeFlux(Vector<LevelData<FluxBox>* >& fluxOfIce, 
		 const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		 const Vector<LevelData<FArrayBox>* >& topography,
		 const Vector<LevelData<FArrayBox>* >& thickness,
		 const Vector<LevelData<FArrayBox>* >& iceFrac,
		 const Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
		 const Vector<LevelData<FArrayBox>* >& basalThicknessSource, 
		 const Vector<Real>& dx, Real dt,
		 const Vector<int>& ratio, 
		 const Vector<std::string>& name, 
		 const Vector<LevelData<FArrayBox>* >& data)
{

  
  int numLevels = data.size();
  Vector<LevelData<FArrayBox>* > ccVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcThck(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      int comp = 0;
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccVel[lev] = new LevelData<FArrayBox>(grids,SpaceDim,2*IntVect::Unit);
      fcVel[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      fcThck[lev] = new LevelData<FluxBox>(grids,1,IntVect::Zero);
      const LevelSigmaCS& levelCS = *coords[lev];

      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "xVel")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccVel[lev],Interval(0,0));
	      comp++;
	    }
	  else if (name[j] == "yVel")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccVel[lev],Interval(1,1));
	      comp++;
	    }
	}
      CH_assert(comp == 2);

      ccVel[lev]->exchange();
      if (lev > 0)
	{
	  const DisjointBoxLayout& crseGrids = topography[lev-1]->disjointBoxLayout();
	  PiecewiseLinearFillPatch velFiller(grids , crseGrids, ccVel[lev]->nComp(), 
					     crseGrids.physDomain(), ratio[lev-1], 2);
	  Real time_interp_coeff = 0.0;
	  velFiller.fillInterp(*ccVel[lev],*ccVel[lev-1] ,*ccVel[lev-1],
			       time_interp_coeff,0, 0, ccVel[lev]->nComp());
	}

      CellToEdge(*ccVel[lev], *fcVel[lev]);

      //modification to fluxes at the margins, that is where mask changes to open sea or land.
      for (DataIterator dit(grids); dit.ok(); ++dit)
      	{

      	  for (int dir = 0; dir < SpaceDim; ++dir)
      	      {
      		Box faceBox = grids[dit];
      		faceBox.surroundingNodes(dir);
      		FArrayBox& faceVel = (*fcVel[lev])[dit][dir];
      		Box grownFaceBox = faceBox;
      		CH_assert(faceVel.box().contains(grownFaceBox));
      		FArrayBox vface(faceBox,1);
      		FArrayBox faceVelCopy(faceVel.box(), 1); faceVelCopy.copy(faceVel);
      		const FArrayBox& cellVel = (*ccVel[lev])[dit];
      		const FArrayBox& usrf = levelCS.getSurfaceHeight()[dit];
      		const FArrayBox& thck = (*thickness[lev])[dit];
      		const FArrayBox& topg = (*topography[lev])[dit];

      		FORT_EXTRAPTOMARGIN(CHF_FRA1(faceVel,0),
                                    CHF_FRA1(vface,0),
      				    CHF_CONST_FRA1(faceVelCopy,0),
      				    CHF_CONST_FRA1(cellVel,dir),
      				    CHF_CONST_FRA1(usrf,0),
      				    CHF_CONST_FRA1(topg,0),
      				    CHF_CONST_FRA1(thck,0),
      				    CHF_CONST_INT(dir),
      				    CHF_BOX(faceBox));
      	      }
      	}

      // face centered thickness from PPM
      RealVect levelDx = RealVect::Unit * dx[lev];
      RefCountedPtr<LevelData<FArrayBox> > levelThck(thickness[0]); levelThck.neverDelete();
      RefCountedPtr<LevelData<FArrayBox> > levelTopg(topography[0]); levelTopg.neverDelete();
      LevelDataIBC thicknessIBC(levelThck,levelTopg,levelDx);
      AdvectPhysics advectPhys;
      advectPhys.setPhysIBC(&thicknessIBC);
      int normalPredOrder = 2;
      bool useFourthOrderSlopes = true;
      bool usePrimLimiting = true;
      bool useCharLimiting = false;
      bool useFlattening = false;
      bool useArtificialViscosity = false;
      Real artificialViscosity = 0.0;
      PatchGodunov patchGod;
      patchGod.define(grids.physDomain(), dx[lev], &advectPhys,
		      normalPredOrder,
		      useFourthOrderSlopes,
		      usePrimLimiting,
		      useCharLimiting,
		      useFlattening,
		      useArtificialViscosity,
		      artificialViscosity);
      AdvectPhysics* advectPhysPtr = dynamic_cast<AdvectPhysics*>(patchGod.getGodunovPhysicsPtr());
      CH_assert(advectPhysPtr != NULL);  
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  
	  FluxBox& fcvel = (*fcVel[lev])[dit];
	  FArrayBox& ccvel = (*ccVel[lev])[dit];
	  FluxBox& faceh = (*fcThck[lev])[dit];
	  FArrayBox& cch = (*thickness[lev])[dit];
	  
	  FArrayBox src(grids[dit],1); 
	  src.copy( (*surfaceThicknessSource[lev])[dit]) ;
	  src.plus( (*basalThicknessSource[lev])[dit]) ;
	  CH_assert(advectPhysPtr != NULL);
	  
	  advectPhysPtr->setVelocities(&ccvel,&fcvel);

	  patchGod.setCurrentTime(0.0);
	  patchGod.setCurrentBox(grids[dit]);
	  patchGod.computeWHalf(faceh, cch, src, dt, grids[dit]);

	}

      //work out uh
      for (DataIterator dit(grids);dit.ok();++dit)
    	{	  
	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      FArrayBox& flux = (*fluxOfIce[lev])[dit][dir];
	      const FArrayBox& vel = (*fcVel[lev])[dit][dir];
   
	      //CH_assert(vel.norm(0) < 1.0e+12);

	      flux.copy((*fcVel[lev])[dit][dir]);
	      flux.mult((*fcThck[lev])[dit][dir]);
	     
    	    } // end loop over direction
    	} // end loop over grids

      fluxOfIce[lev] -> exchange();

    } //end loop over levels


  for (int lev = 0; lev < numLevels; lev++)
    {
      if (ccVel[lev] != NULL)
	{
	  delete ccVel[lev];ccVel[lev]= NULL;
	}
      if (fcVel[lev] != NULL)
	{
	  delete fcVel[lev];fcVel[lev]= NULL;
	}
      if (fcThck[lev] != NULL)
	{
	  delete fcThck[lev];fcThck[lev]= NULL;
	}
      
    }


}

void stateDiagnostics(std::ostream& sout, bool append, std::string plot_file,
		      const Vector<LevelData<FArrayBox>* >& sectorMaskFraction,
		      int maskNoStart, int maskNoEnd,
		      const Vector<LevelData<FArrayBox>* >& data,
		      const Vector<DisjointBoxLayout >& grids,
		      const Vector<std::string>& name,
		      int numLevels,
		      Vector<Real>& dx, // cant be const due to createSigmaCS
		      Real dt, Real time,
		      Vector<int>& ratio,// cant be const due to createSigmaCS
		      Real iceDensity,
		      Real waterDensity,
		      Real gravity,
		      Real h_min, Real f_min)
{

  CH_TIME("stateDiagnostics");

  Vector<LevelData<FArrayBox>* > thickness(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > topography(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > deltaThickness(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > surfaceThicknessSource(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > basalThicknessSource(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > divergenceThicknessFlux(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > calvingFlux(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > melangeThickness(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > iceFrac(numLevels,NULL);
  Vector<LevelData<FluxBox>* > fluxOfIce(numLevels, NULL);
  

  for (int lev=0;lev<numLevels;++lev)
    {
      // Extra ghost cells are needed to calculate advection using PPM. 
      thickness[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
      topography[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
      deltaThickness[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      surfaceThicknessSource[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      basalThicknessSource[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      divergenceThicknessFlux[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      calvingFlux[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      melangeThickness[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      iceFrac[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      fluxOfIce[lev] = new LevelData<FluxBox>(grids[lev],1,IntVect::Zero);  
    }

    
  sout.setf(ios_base::scientific,ios_base::floatfield); 
  sout.precision(12);
  
  createDEM(topography, thickness, iceFrac, name, data, ratio, dx, h_min);
  
  Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
  createSigmaCS(coords,topography, thickness, iceFrac, 
		dx, ratio, iceDensity, waterDensity,gravity);
  
  extractThicknessSource(surfaceThicknessSource, basalThicknessSource, deltaThickness, 
			 divergenceThicknessFlux, calvingFlux, melangeThickness,
			 topography, dx, ratio, name, data);
  
  computeFlux(fluxOfIce,coords,topography,thickness,iceFrac,surfaceThicknessSource,basalThicknessSource,
	      dx,dt,ratio,name,data);
  
  
  // CSV style output of diagnostics
  if (!append)
    {
      sout << "csvheader,filename,time,maskNo,region,quantity,unit,value" << endl;
    }
  
  for (int maskNo = maskNoStart; maskNo <= maskNoEnd; ++maskNo)
    {
      int maskComp = maskNo - maskNoStart; // component of sectorMaskFraction FABs corresponding to maskNo
       
      typedef std::map<std::string, function<bool(Real, Real, int)> > MapSF;
      MapSF regions;
      
      auto entire = [h_min](Real h, Real f, int mask){ return true;} ;
      regions["entire"] = entire;
      
      auto grounded = [h_min](Real h, Real f, int mask){ return ((mask == GROUNDEDMASKVAL) && (h > h_min));} ;
      regions["grounded"] = grounded;
      
      auto floating = [h_min](Real h, Real f, int mask){ return ((mask == FLOATINGMASKVAL) && (h > h_min));} ;
      regions["floating"] = floating;
      
      auto ice = [h_min](Real h, Real f, int mask){ return (h > h_min);} ;
      regions["ice"] = ice;
      
      auto nonice = [h_min](Real h, Real f, int mask){ return (h <= h_min);} ;
      regions["nonice"] = nonice;
      
      
      for (MapSF::const_iterator mit = regions.begin(); mit != regions.end(); ++mit)
	{
	  Vector<NameUnitValue> report;
	 
	  reportConservationInside(report, mit->second,coords, fluxOfIce, surfaceThicknessSource,
				   basalThicknessSource, deltaThickness, divergenceThicknessFlux,
				   calvingFlux,topography, thickness, iceFrac,
				   sectorMaskFraction, dx, ratio, numLevels-1, maskNo, maskComp);
	  
	  for (int i = 0; i < report.size(); i++)
	    {
	      sout << "csvdata," << plot_file << "," << time << "," << maskNo << "," << mit->first << "," << report[i] << endl;
	    }
	}
    }
    
  for (int lev=0;lev<numLevels;++lev)
    {
      //if (sectorMask[lev] != NULL) delete sectorMask[lev];
      if (thickness[lev] != NULL) delete thickness[lev];
      if (topography[lev] != NULL) delete topography[lev];
      if (deltaThickness[lev] != NULL) delete deltaThickness[lev];
      if (surfaceThicknessSource[lev] != NULL) delete surfaceThicknessSource[lev];
      if (basalThicknessSource[lev] != NULL) delete basalThicknessSource[lev];
      if (divergenceThicknessFlux[lev] != NULL) delete divergenceThicknessFlux[lev];
      if (calvingFlux[lev] != NULL) delete calvingFlux[lev];
      if (melangeThickness[lev] != NULL) delete melangeThickness[lev];
      if (iceFrac[lev] != NULL) delete iceFrac[lev];
      if (fluxOfIce[lev] != NULL) delete fluxOfIce[lev];
    }
}


int main(int argc, char* argv[]) {

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif 

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif

    if(argc < 2) 
      {
	std::cerr << " usage: " << argv[0] << "plot_file=plot_file [-options] [key=value args]" << endl;
	std::cerr << " example: " << argv[0] << "plot_file=test.2d.hdf5 out_file=test.csv -append ice_density=918.0 water_density=1028.0"  << endl;
	std::cerr << "        computes diagnostic integrals from test.2d.hdf5, appends results to test.csv " << endl;
	std::cerr << " example: " << argv[0] << "plot_file=test.2d.hdf5 out_file=test.csv -append ice_density=918.0 water_density=1028.0 mask_file=mask2.d.hdf5 mask_no_start=0 mask_no_end=4"  << endl;
	std::cerr << "        computes diagnostic integrals from test.2d.hdf5, for the whole domain and each subdomain 1-4, appends results to test.csv, mask.2d.hdf5 indicates subdomains with a grid of values 1-4 " << endl;
	exit(1);
      }
  
    ParmParse pp(argc-1,argv+1,NULL,NULL);

    // the one mandatory arg : a plot file
    std::string plot_file; pp.get("plot_file",plot_file);

    // out_file : if not specified use stdio
    std::string out_file(""); pp.query("out_file", out_file);
    std::ofstream fout;
    bool append = pp.contains("append");
    if (out_file != "")
	{
	  fout.open(out_file, append?std::ofstream::app:std::ofstream::trunc);
	}
    std::ostream& sout = (out_file == "")?pout():fout;
      
    
    Real ice_density(917.0); pp.query("ice_density", ice_density);
    Real water_density(1028.0); pp.query("water_density", water_density);
    Real gravity(9.81); pp.query("gravity", gravity);

    std::string mask_file("");
    pp.query("mask_file",mask_file);
    bool maskFile = mask_file.size() > 0;

    int mask_no_start(0); pp.query("mask_no_start", mask_no_start);
    int mask_no_end(mask_no_start); pp.query("mask_no_end", mask_no_end);

    Real h_min(100.0); pp.query("h_min", h_min);
    Real f_min(1.0e-1); pp.query("f_min", f_min); 
    

    Box domainBox;
    Vector<std::string> name;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout > grids;
    Vector<int > ratio;
    int numLevels;

    Real dt ,crseDx, time;
   
    ReadAMRHierarchyHDF5(std::string(plot_file), grids, data, name , 
			 domainBox, crseDx, dt, time, ratio, numLevels);

    Vector<ProblemDomain> domain(numLevels,domainBox);
    Vector<RealVect> vdx(numLevels,RealVect::Unit*crseDx);
    Vector<Real> dx(numLevels,crseDx);
    for (int lev=1;lev<numLevels;++lev)
      {
	dx[lev] = dx[lev-1] / Real(ratio[lev-1]);
	vdx[lev] = vdx[lev-1] / Real(ratio[lev-1]);
	domain[lev] = domain[lev-1];
	domain[lev].refine(ratio[lev-1]);
      }

    Vector<LevelData<FArrayBox>* > sector_mask_fraction(numLevels);
    if (maskFile)
      {
	//load the sector mask, if it exists    
	Box mdomainBox;
	Vector<std::string> mname;
	Vector<LevelData<FArrayBox>* > mdata;
	Vector<DisjointBoxLayout > mgrids;
	Vector<int > mratio;
	int mnumLevels;
	Real mdt ,mcrseDx, mtime;
	ReadAMRHierarchyHDF5(std::string(mask_file), mgrids, mdata, mname , 
			     mdomainBox, mcrseDx, mdt, mtime, mratio, mnumLevels);

	if (mnumLevels > 1) MayDay::Error("mask files assumed to have one AMR level");

	// create storage for sector mask fractions - one FAB component for each mask label
	int n_mask_nos = 1 + mask_no_end - mask_no_start;
	
	for (int lev = 0; lev < data.size(); lev++)
	  {
	    sector_mask_fraction[lev] = new LevelData<FArrayBox>(grids[lev],n_mask_nos,IntVect::Unit);
	  }
	
	maskToAMR(sector_mask_fraction, dx, *mdata[0], mask_no_start, n_mask_nos, mcrseDx);

	//free heap
	for (int lev = 0; lev < mdata.size(); lev++)
	  {
	    if (mdata[lev]) delete mdata[lev];
	  }
	
      }
    

    
    
    stateDiagnostics(sout, append, plot_file, sector_mask_fraction, mask_no_start, mask_no_end,
		     data, grids, name, numLevels,  dx, dt, time,  ratio,
		     ice_density, water_density, gravity, h_min, f_min);


    // free heap
    for (int lev = 0; lev < data.size(); lev++)
      {
	if (sector_mask_fraction[lev]) delete sector_mask_fraction[lev];
	if (data[lev]) delete data[lev];
      }
   
    		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
