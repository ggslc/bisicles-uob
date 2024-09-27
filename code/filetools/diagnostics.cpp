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
// read in a bisicles plotfile (which must) 
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
#include <functional>


typedef BaseFab<int> IArrayBox; // just so it lines up neatly with FArrayBox :)

void FtoI(IArrayBox& dest, const FArrayBox& src)
{
  CH_assert(src.box().contains(dest.box()));
  CH_assert(src.nComp() >= dest.nComp());
  for (int comp = 0; comp < dest.nComp(); comp++)
    {
      for (BoxIterator bit(dest.box()); bit.ok(); ++bit)
	{
	  dest(bit(), comp) = src(bit(),comp);
	}
    }
}


void integrateScalarInside(Real& integral, 
		     function<bool(int mask, Real thickness)> inside,
		     const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		     const Vector<LevelData<FArrayBox>* >& integrand, 
		     const Vector<LevelData<FArrayBox>* >& topography,
		     const Vector<LevelData<FArrayBox>* >& thickness, 
		     const Vector<LevelData<IArrayBox>* >& sectorMask,
		     const Vector<Real>& dx, const Vector<int>& ratio, 
		     int maskNo)
{
  CH_TIME("integrateScalarInside");
  
  int numLevels = coords.size();
  Vector<LevelData<FArrayBox>* > maskedIntegrand(numLevels, NULL);
  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& grids = integrand[lev]->disjointBoxLayout();
      maskedIntegrand[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  const BaseFab<int>& mask = coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
	  const Box& b = grids[dit];
	  const FArrayBox& f = (*integrand[lev])[dit];
	  FArrayBox& g = (*maskedIntegrand[lev])[dit];
	  g.setVal(0.0);
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if (((*sectorMask[lev])[dit](iv) == maskNo) || maskNo == 0)
		{
		  if (inside(thck(iv), mask(iv)))
		    {
		      g(iv) = f(iv);
		    }
		}
	    } // bit
	} //dit
    } //lev

  integral = computeSum(maskedIntegrand, ratio, dx[0], Interval(0,0), 0);

  for (int lev = 0; lev < numLevels; lev++)
    {
      if (maskedIntegrand[lev] != NULL)
	{
	  delete maskedIntegrand[lev];maskedIntegrand[lev]= NULL;
	}  
    }
  
}

void integrateDischargeInside(Real& sumDischarge, Real& sumDivUH,
			      function<bool(int mask, Real thickness)> inside,
			      const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
			      const Vector<LevelData<FluxBox>* >& fluxOfIce,
			      const Vector<LevelData<FArrayBox>* >& topography,
			      const Vector<LevelData<FArrayBox>* >& thickness, 
			      const Vector<Real>& dx, const Vector<int>& ratio, 
			      const Vector<LevelData<IArrayBox>* >& sectorMask,  
			      int maskNo)
{

  CH_TIME("integrateDischargeInside");
  
  int numLevels = coords.size();
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDivUH(numLevels, NULL);
  for (int lev = 0; lev < numLevels; lev++)
    {
      int comp = 0;
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDivUH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
     
      //work out discharge across boundary and flux divergence inside region 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
	  FArrayBox& divuh = (*ccDivUH[lev])[dit];
	  divuh.setVal(0.0);
	  const BaseFab<int>& mask = coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  
    	  for (int dir =0; dir < SpaceDim; dir++)
    	    {
	      const FArrayBox& flux  = (*fluxOfIce[lev])[dit][dir];
    	      for (BoxIterator bit(b);bit.ok();++bit)
    		{
    		  const IntVect& iv = bit();
    		  if (((*sectorMask[lev])[dit](iv) == maskNo) || maskNo == 0)
		    {
		      if (inside(thck(iv), mask(iv)))
		       	{
			  //inside the sub-region - compute div(flux)
		       	  divuh(iv) += (flux(iv + BASISV(dir)) - flux(iv))/dx[lev];
		       	}
		      else
    			{
			  //outside - compute discharge/dx if neigbours are inside
			  if (inside(thck(iv + BASISV(dir)), mask(iv + BASISV(dir))))
			    {
			      discharge(iv) += -flux(iv + BASISV(dir)) / dx[lev];
    			    }
			  if (inside(thck(iv - BASISV(dir)), mask(iv - BASISV(dir))))
    			    {
    			      discharge(iv) += flux(iv) / dx[lev];
    			    }
    			}//end if inside else
    		    }//end if mask
    		} //end loop over cells
    	    } // end loop over direction
    	} // end loop over grids
    } //end loop over levels

  sumDischarge = computeSum(ccDischarge, ratio, dx[0], Interval(0,0), 0);
  sumDivUH = computeSum(ccDivUH, ratio, dx[0], Interval(0,0), 0);
  
  for (int lev = 0; lev < numLevels; lev++)
    {
      
      if (ccDischarge[lev] != NULL)
	{
	  delete ccDischarge[lev];ccDischarge[lev]= NULL;
	}
      if (ccDivUH[lev] != NULL)
	{
	  delete ccDivUH[lev];ccDivUH[lev]= NULL;
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


  void reportConservationInside(//std::map<std::string, Real>& values,
				Vector<NameUnitValue>& report,
			      function<bool(int mask, Real thickness)> inside,
			      const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
			      const Vector<LevelData<FluxBox>* >& fluxOfIce,
			      const Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
			      const Vector<LevelData<FArrayBox>* >& basalThicknessSource,
			      const Vector<LevelData<FArrayBox>* >& deltaThickness,
			      const Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux,
			      const Vector<LevelData<FArrayBox>* >& calvingFlux,
			      const Vector<LevelData<FArrayBox>* >& topography,
			      const Vector<LevelData<FArrayBox>* >& thickness, 
			      const Vector<LevelData<IArrayBox>* >& sectorMask,
			      const Vector<Real>& dx, const Vector<int>& ratio, 
			      int maskNo)
{

  CH_TIME("reportConservationInside");
  // just to make the args less reptitive...
  auto sumScalar = [inside,coords,topography, thickness, sectorMask, dx, ratio, maskNo ]
    (const Vector<LevelData<FArrayBox>* >& scalar)
		   {
		     Real sumScalar;
		     integrateScalarInside(sumScalar, inside ,coords, scalar, topography, thickness, 
					   sectorMask, dx, ratio, maskNo);
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
  integrateDischargeInside(sumDischarge, flxDivReconstr,
			   inside, coords, fluxOfIce,topography,
			   thickness, dx, ratio, sectorMask, maskNo);
  report.push_back(NameUnitValue("flxDivReconstr",dhunit,flxDivReconstr));
  report.push_back(NameUnitValue("discharge",dhunit,sumDischarge));

  //Conservation errors
  report.push_back(NameUnitValue("SMB+BMB-dhdt-calving-flxDivFile",dhunit,
				 SMB + BMB - dhdt - calving - flxDivFile));
  report.push_back(NameUnitValue("SMB+BMB-dhdt-calving-flxDivReconstr",dhunit,
				 SMB + BMB - dhdt - calving - flxDivReconstr));
		   
  // volumes
  std::string vunit("m3");
  report.push_back(NameUnitValue("volume",vunit,sumScalar(thickness)));
  // volume abive flotation
  {
    Vector<LevelData<FArrayBox>* > hab;
    for (int lev = 0; lev < coords.size(); lev++)
      {
	hab.push_back(const_cast<LevelData<FArrayBox>* >(&(coords[lev]->getThicknessOverFlotation())));
      }
    report.push_back(NameUnitValue("volumeAbove",vunit,sumScalar(hab)));
  }
  
  // Area calculation. Maybe put this elsehwere
  {
    Vector<LevelData<FArrayBox>* > one(thickness.size(),NULL);
    for (int lev = 0; lev < one.size(); lev++)
      {
	const DisjointBoxLayout& grids = (*thickness[lev]).disjointBoxLayout();
	one[lev] = new LevelData<FArrayBox>(grids, 1, IntVect::Zero);
	for (DataIterator dit(grids); dit.ok(); ++dit)
	  (*one[lev])[dit].setVal(1.0);
      }
    std::string aunit("m2");
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
		Vector<std::string>& name, 
		Vector<LevelData<FArrayBox>* >& data,
		Vector<int>& ratio,
		Vector<Real>& dx,
		Real mcrseDx)
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
      ccVel[lev] = new LevelData<FArrayBox>(grids,SpaceDim,IntVect::Unit);
      fcVel[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      fcThck[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      const LevelSigmaCS& levelCS = *coords[lev];

      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  (*ccVel[lev])[dit].setVal(0.0);
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      (*fcVel[lev])[dit][dir].setVal(0.0);
	      (*fluxOfIce[lev])[dit][dir].setVal(0.0);
	    }
	}

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
					     crseGrids.physDomain(), ratio[lev-1], 1);
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
   
	      CH_assert(vel.norm(0) < 1.0e+12);

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
   
    
    if(argc < 5) 
      { 
	std::cerr << " usage: " << argv[0] << " <plot file> <ice_density> <water_density> <gravity> [mask_file] [mask_no_start = 0] [mask_no_end = mask_no_start] " << std::endl; 
	exit(0); 
      }
    char* plotFile = argv[1];
    Real iceDensity = atof(argv[2]);
    Real waterDensity = atof(argv[3]);
    Real gravity = atof(argv[4]);
    char* maskFile = (argc > 5)?argv[5]:NULL;
    int maskNoStart = 0;
    if (maskFile && argc > 6)
      {
    	maskNoStart = atoi(argv[6]);
      }
    int maskNoEnd = maskNoStart;
    if (maskFile && argc > 7)
      {
    	maskNoEnd = atoi(argv[7]);
      }
    
    Real Hmin = 1.0e-10;

    Box domainBox;
    Vector<std::string> name;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout > grids;
    Vector<int > ratio;
    int numLevels;

    Real dt ,crseDx, time;
   
    ReadAMRHierarchyHDF5(std::string(plotFile), grids, data, name , 
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
    //load the sector mask, if it exists    
    Box mdomainBox;
    Vector<std::string> mname;
    Vector<LevelData<FArrayBox>* > mdata;
    Vector<DisjointBoxLayout > mgrids;
    Vector<int > mratio;
    int mnumLevels;
    Real mdt ,mcrseDx, mtime;

    if (maskFile)
      {   	
	ReadAMRHierarchyHDF5(std::string(maskFile), mgrids, mdata, mname , 
			     mdomainBox, mcrseDx, mdt, mtime, mratio, mnumLevels);
      }
   

    
    Vector<LevelData<IArrayBox>* > sectorMask(numLevels,NULL);
    for (int lev=0;lev<numLevels;++lev)
      {
	sectorMask[lev] = new LevelData<IArrayBox>(grids[lev],1,IntVect::Unit);

	for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
	  {
	    (*sectorMask[lev])[dit].setVal(0);
	  }

	if (maskFile)
	  {
	    // reading real values...
	    LevelData<FArrayBox> tmp(grids[lev],1,IntVect::Unit);
	    FillFromReference(tmp, *mdata[0], RealVect::Unit*dx[lev], RealVect::Unit*mcrseDx,true);
	    for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
	      {
		FtoI((*sectorMask[lev])[dit],tmp[dit]);
	      }
	  }

      }

    Vector<LevelData<FArrayBox>* > thickness(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > topography(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > deltaThickness(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > surfaceThicknessSource(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > basalThicknessSource(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > divergenceThicknessFlux(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > calvingFlux(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > melangeThickness(numLevels,NULL);
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
	fluxOfIce[lev] = new LevelData<FluxBox>(grids[lev],1,IntVect::Unit);

      }

    
    pout().setf(ios_base::scientific,ios_base::floatfield); 
    pout().precision(12);

    createDEM(topography, thickness, name, data, ratio, dx, mcrseDx);

    Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
    createSigmaCS(coords,topography, thickness, 
		  dx, ratio, iceDensity, waterDensity,gravity);

    extractThicknessSource(surfaceThicknessSource, basalThicknessSource, deltaThickness, 
		       divergenceThicknessFlux, calvingFlux, melangeThickness,
		       topography, dx, ratio, name, data);

    computeFlux(fluxOfIce,coords,topography,thickness,surfaceThicknessSource,basalThicknessSource,
		dx,dt,ratio,name,data);

    
    // CSV style output of diagnostics
    pout() << "csvheader,time,maskNo,region,quantity,unit,value" << endl;
    for (int maskNo = maskNoStart; maskNo <= maskNoEnd; ++maskNo)
      {
	
	typedef std::map<std::string, function<bool(int mask, Real thickness)> > MapSF;
	MapSF regions;

	auto entire = [Hmin](Real h, int mask){ return true;} ;
	regions["entire"] = entire;
	
	auto grounded = [Hmin](Real h, int mask){ return ((mask == GROUNDEDMASKVAL) && (h > Hmin));} ;
	regions["grounded"] = grounded;

	auto floating = [Hmin](Real h, int mask){ return ((mask == FLOATINGMASKVAL) && (h > Hmin));} ;
	regions["floating"] = floating;

	auto ice = [Hmin](Real h, int mask){ return (h > Hmin);} ;
	regions["ice"] = ice;

	auto nonice = [Hmin](Real h, int mask){ return (h < Hmin);} ;
	regions["nonice"] = nonice;

       
	for (MapSF::const_iterator mit = regions.begin(); mit != regions.end(); ++mit)
	  {
	    Vector<NameUnitValue> report;
	    reportConservationInside(report, mit->second,coords, fluxOfIce, surfaceThicknessSource,
				     basalThicknessSource, deltaThickness, divergenceThicknessFlux,
				     calvingFlux,topography, thickness, sectorMask, dx, ratio, maskNo);
	    for (int i = 0; i < report.size(); i++)
	      {
		pout() << "csvdata," << time << "," << maskNo << "," << mit->first << "," << report[i] << endl;
	      }
	  }
      }
    
    for (int lev=0;lev<numLevels;++lev)
      {
	if (sectorMask[lev] != NULL) delete sectorMask[lev];
	if (thickness[lev] != NULL) delete thickness[lev];
	if (topography[lev] != NULL) delete topography[lev];
	if (deltaThickness[lev] != NULL) delete deltaThickness[lev];
	if (surfaceThicknessSource[lev] != NULL) delete surfaceThicknessSource[lev];
	if (basalThicknessSource[lev] != NULL) delete basalThicknessSource[lev];
	if (divergenceThicknessFlux[lev] != NULL) delete divergenceThicknessFlux[lev];
	if (calvingFlux[lev] != NULL) delete calvingFlux[lev];
	if (melangeThickness[lev] != NULL) delete melangeThickness[lev];
	if (fluxOfIce[lev] != NULL) delete fluxOfIce[lev];
      }

		  
  }  // end nested scope
  CH_TIMER_REPORT();
  dumpmemoryatexit();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
