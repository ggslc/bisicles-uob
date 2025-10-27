#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ConstitutiveRelation.H"
#include "LevelMappedDerivatives.H"
#include "CellToEdge.H"
#include "ConstitutiveRelationF_F.H"
#include "IceConstants.H"
#include "IceThermodynamics.H"
#include "TensorCFInterp.H"
#include "ViscousTensorOpF_F.H"
#include "ViscousTensorOp.H"
#include "ParmParse.H"
#include "L1L2ConstitutiveRelation.H"
#include "NamespaceHeader.H"

void
ConstitutiveRelation::computeStrainRateInvariant(LevelData<FArrayBox>& a_epsilonSquared,
                                                 const LevelData<FArrayBox>& a_velocity,
                                                 const LevelData<FArrayBox>* a_crseVel,
                                                 int a_nRefCrse,
                                                 const LevelSigmaCS& a_coordSys,
                                                 const IntVect& a_ghostVect) const
{
  const DisjointBoxLayout& grids = a_epsilonSquared.getBoxes();
  LevelData<FArrayBox> derivs(grids, SpaceDim*SpaceDim, a_ghostVect);
  computeStrainRateInvariant(a_epsilonSquared,derivs,a_velocity,
                             a_crseVel, a_nRefCrse,
                             a_coordSys, a_ghostVect);

}

void
ConstitutiveRelation::computeStrainRateInvariantFace(LevelData<FluxBox>& a_epsilonSquared,
                                                     LevelData<FArrayBox>& a_velocity,
                                                     const LevelData<FArrayBox>* a_crseVel,
                                                     int a_nRefCrse,
                                                     const LevelSigmaCS& a_coordSys,
                                                     const IntVect& a_ghostVect) const
{
  
  const DisjointBoxLayout& grids = a_epsilonSquared.getBoxes();
  LevelData<FluxBox> derivs(grids, SpaceDim*SpaceDim, a_ghostVect);

  computeStrainRateInvariantFace(a_epsilonSquared, derivs, a_velocity,
                                 a_crseVel, a_nRefCrse, a_coordSys,
                                 a_ghostVect);

}



// extrapolate the velocity beyond the front: needed to avoid large artificial strains.
// From Sam K & Dan M
void ConstitutiveRelation::extendVelocity(LevelData<FArrayBox>& a_modifiedVel,
					  const LevelData<FArrayBox>& a_originalVel,
					  const LevelData<FArrayBox>& a_iceThickness) const
{
  
  DataIterator dit = a_modifiedVel.dataIterator();
  Real vel_eps = 1.0e-8;
  Real thickness_eps = 0.1;

  // turns out this might need a cornerCopier
  // put this into its own context because we're playing fast
  // and loose with const-ness 
  {
    const DisjointBoxLayout& grids = a_originalVel.getBoxes();
    const ProblemDomain& domain = grids.physDomain();
    LevelData<FArrayBox>& nonConstVel = *(const_cast<LevelData<FArrayBox>* >(&a_originalVel));
    CornerCopier cc(grids, grids, domain,
                    nonConstVel.ghostVect(), true);
    
    nonConstVel.exchange(cc);
  } // end cornerCopier scope

  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisVel = a_modifiedVel[dit];
      const FArrayBox& thisH = a_iceThickness[dit];
      thisVel.copy(a_originalVel[dit]);
      
      // loop over directions first because we want to be able to do ghost cells where possible
      for (int dir=0; dir<SpaceDim; dir++)
        {          
          Box velBox = a_modifiedVel[dit].box();
          velBox.grow(dir,-1);
          BoxIterator bit(velBox);
          SideIterator sit;                  
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              // only do this for openSea cells
              if (abs(a_originalVel[dit](iv,0)) < vel_eps)
                {
                  for (sit.begin(); sit.ok(); ++sit)
                    {
                      IntVect shiftVect = sign(sit())*BASISV(dir);
                      // are we at a calving front?
                      if (abs(thisH(iv+shiftVect,0)) > thickness_eps)                      

                        {
                          // linearly extrapolate if possible
                          if ((thisVel.box().contains((iv+2*shiftVect))) && abs(thisH(iv+2*shiftVect,0)) > thickness_eps)
                            {
                              thisVel(iv,0) = 2.0*thisVel(iv+shiftVect,0) - thisVel(iv+2*shiftVect,0);
                              thisVel(iv,1) = 2.0*thisVel(iv+shiftVect,1) - thisVel(iv+2*shiftVect,1);
                            }
                          else
                            {
                              // just copy if stencil doesn't exist
                              thisVel(iv,0) = thisVel(iv+shiftVect,0);
                              thisVel(iv,1) = thisVel(iv+shiftVect,1);
                            }
                        } // end if neigboring cell is floating
                    } // end loop over hi-lo
                } // end if this is an opensea cell
            } // end loop over cells in this box
        } // end loop over directions
    } // end loop over boxes
}



void
ConstitutiveRelation::computeStrainRateInvariant(LevelData<FArrayBox>& a_epsilonSquared,
						 LevelData<FArrayBox>& a_gradVelocity,
                                                 const LevelData<FArrayBox>& a_velocity,
                                                 const LevelData<FArrayBox>* a_crseVel,
                                                 int a_nRefCrse,
                                                 const LevelSigmaCS& a_coordSys,
                                                 const IntVect& a_ghostVect) const
{
  // first compute derivatives...
  // need all derivatives of all velocities
  Interval DerivDir(0,SpaceDim-1);
  Interval DerivComps(0,SpaceDim-1);

  CH_assert(a_gradVelocity.nComp() == SpaceDim*SpaceDim);

  bool mask = true; // switch to one-sided differences  next to fluid free cells, avoids spurious large gradient
  computeCCDerivatives(a_gradVelocity, a_velocity, a_coordSys,
                       DerivComps, DerivDir, a_ghostVect, mask);

 
  const DisjointBoxLayout& grids = a_epsilonSquared.getBoxes();
  
  DataIterator dit = grids.dataIterator();
  
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box derivBox = grids[dit];
      derivBox.grow(a_ghostVect);

      
      // now compute strainInvariant
      if (SpaceDim != 3)
        {
          // this is a bit hokey, just do it to ensure that
          // we get the components correct. Eventually hardwire
          // this into the fortran.
          int dudxComp = derivComponent(0,0);
          int dudyComp = derivComponent(1,0);
          int dvdxComp = derivComponent(0,1);
          int dvdyComp = derivComponent(1,1);
          
          // compute shallow-shelf approximation to epsilonSquared
          FORT_STRAININVARSSA(CHF_FRA1(a_epsilonSquared[dit],0),
                              CHF_FRA(a_gradVelocity[dit]),
                              CHF_INT(dudxComp),
                              CHF_INT(dudyComp),
                              CHF_INT(dvdxComp),
                              CHF_INT(dvdyComp),
                              CHF_BOX(derivBox));
          

        }
      else
        {
          MayDay::Error("3D cell-centered strain invariant not implemented yet");
        }
     
    } // end loop over grids 
}

/// compute face-centered epsilon^2 and grda(u) based on cell-centered velocity
void
ConstitutiveRelation::computeStrainRateInvariantFace(LevelData<FluxBox>& a_epsilonSquared,
						     LevelData<FluxBox>& a_gradVelocity,
                                                     LevelData<FArrayBox>& a_velocity,
                                                     const LevelData<FArrayBox>* a_crseVelPtr,
                                                     int a_nRefCrse,
                                                     const LevelSigmaCS& a_coordSys,
                                                     const IntVect& a_ghostVect) const
{

  // this looks really similar to the cell-centered version
  // first compute derivatives...
  // need all derivatives of all velocities
  Interval DerivDir(0,SpaceDim-1);
  Interval DerivComps(0,SpaceDim-1);

  CH_assert(a_gradVelocity.nComp() == SpaceDim*SpaceDim);  
  const DisjointBoxLayout& grids = a_epsilonSquared.getBoxes();
  
  //LevelData<FArrayBox> extVel(grids, SpaceDim, a_velocity.ghostVect());
  //extendVelocity(extVel ,a_velocity,a_coordSys.getH());
  
#define TENSORCF
#ifdef TENSORCF
  //calculate the ghost cell values of velocity and its tangential
  //gradients.

  LevelData<FArrayBox>& vel = a_velocity;
  LevelData<FArrayBox> gradvel(grids,SpaceDim*SpaceDim,IntVect::Unit);
  Real dx = a_coordSys.dx()[0];
  if (a_crseVelPtr != NULL)
    {
      //ghost cells of velocity and its transverse gradient components
      const DisjointBoxLayout& crseGrids = a_crseVelPtr->getBoxes();
      
      const ProblemDomain& pd = grids.physDomain();
      TensorCFInterp cf(grids,&crseGrids,dx,a_nRefCrse,SpaceDim,pd);
      cf.coarseFineInterp(vel,gradvel,*a_crseVelPtr);

      
    }
  
  //compute cell-centered gradients ;
  for (DataIterator dit(grids);dit.ok();++dit)
    {

      const Box& b = grids[dit];
      FArrayBox& gv = gradvel[dit];
      const FArrayBox& v = vel[dit];

      for(int derivDir = 0; derivDir < SpaceDim; derivDir++)
	{
	  for(int dir = 0; dir < SpaceDim; dir++)
	    {
	      int gradcomp = TensorCFInterp::gradIndex(dir,derivDir);
	      FORT_CELLGRADVTOP(CHF_FRA1(gv, gradcomp),
				CHF_CONST_FRA1(v, dir),
				CHF_BOX(b),
				CHF_CONST_REAL(dx),
				CHF_CONST_INT(derivDir));
	      
	    }
	}
    }
  gradvel.exchange();
  //compute face-centered gradients ;
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  FArrayBox& faceGrad = a_gradVelocity[dit][faceDir];
	  FArrayBox& gv = gradvel[dit];
	  FArrayBox& v = vel[dit];
	  Box faceBox(grids[dit]);
	  faceBox.surroundingNodes(faceDir);
	  FArrayBox faceDiv(faceBox,1); // we don't care about this

	  ViscousTensorOp::getFaceDivAndGrad(faceDiv,faceGrad,v,gv,grids.physDomain(),
					     faceBox,faceDir,dx);
	  
	}
    }

  
#else
  computeFCDerivatives(a_gradVelocity, a_velocity, a_coordSys,
                       DerivComps, DerivDir, a_ghostVect);
#endif

  // // replace with one-sided derivatives close to ice front
  //computeFCDerivatives(a_gradVelocity, a_velocity, a_coordSys,
  //                      DerivComps, DerivDir, a_ghostVect, true);
  
  
  // now compute strainInvariant
  
  DataIterator dit = grids.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      Box derivBox = grids[dit];
      derivBox.grow(a_ghostVect);
      
      FluxBox& epsSqr = a_epsilonSquared[dit];
      FluxBox& gradVel = a_gradVelocity[dit];
      
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          Box faceBox = derivBox;
          faceBox.surroundingNodes(faceDir);

            if (SpaceDim != 3)
              {
                // this is a bit hokey, just do it to ensure that
                // we get the components correct. Eventually hardwire
                // this into the fortran.
                int dudxComp = derivComponent(0,0);
                int dudyComp = derivComponent(1,0);
                int dvdxComp = derivComponent(0,1);
                int dvdyComp = derivComponent(1,1);
                
                // compute shallow-shelf approximation to epsilonSquared
                // note that both the cell-centered and face-centered versions
                // call the same fortran here.
                FORT_STRAININVARSSA(CHF_FRA1(epsSqr[faceDir],0),
                                    CHF_FRA(gradVel[faceDir]),
                                    CHF_INT(dudxComp),
                                    CHF_INT(dudyComp),
                                    CHF_INT(dvdxComp),
                                    CHF_INT(dvdyComp),
                                    CHF_BOX(faceBox));
                
		int dbg(0); dbg++;
		
              }
            else
              {
                MayDay::Error("3D face-centered strain invariant not implemented yet");
              }
          
        } // end loop over face directions
    } // end loop over grids
}
  

 


GlensFlowRelation::GlensFlowRelation() 
{
  // initialize parameter values
  setDefaultParameters();
}


GlensFlowRelation::~GlensFlowRelation()
{

}



void
GlensFlowRelation::computeMu(LevelData<FArrayBox>& a_mu,
                             const LevelData<FArrayBox>& a_vel, const Real& a_scale,
                             const LevelData<FArrayBox>* a_crseVel,
                             int a_nRefCrse,
                             const LevelData<FArrayBox>& a_A,
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect) const
{
  // first, compute epsilon^2
  const DisjointBoxLayout& grids = a_mu.getBoxes();
  LevelData<FArrayBox> epsSqr(grids, 1, a_ghostVect);
  computeStrainRateInvariant(epsSqr, a_vel, a_crseVel, a_nRefCrse, a_coordSys, a_ghostVect);

  // now compute mu
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisMu = a_mu[dit];
      const FArrayBox& thisA = a_A[dit];
      Box derivBox = grids[dit];
      derivBox.grow(a_ghostVect);
     
     if (thisA.nComp() == 1 && thisMu.nComp() == 1)
	{
	  // only one layer (or 3D)
	  computeMu0(thisMu, thisA, a_scale, derivBox);
	  FORT_COMPUTEGLENSMU(CHF_FRA1(thisMu,0),
			      CHF_CONST_FRA1(epsSqr[dit],0),
			      CHF_BOX(derivBox),
			      CHF_CONST_REAL(m_n),
			      CHF_CONST_REAL(m_epsSqr0),
			      CHF_CONST_REAL(m_delta)
			      );


	}
#if BISICLES_Z == BISICLES_LAYERED
     else if (thisA.nComp() > 1 && thisMu.nComp() == thisA.nComp())
	{
	  // multiple layers, compute mu in each layer
	  for (int layer = 0; layer < thisA.nComp(); ++layer)
	    {
	      FArrayBox layerMu;
	      layerMu.define(Interval(layer,layer), thisMu);
	      FArrayBox layerA(thisA.box(),1);
	      layerA.copy(thisA,layer,0);
	      computeMu0(layerMu, layerA, a_scale, derivBox);
	      FORT_COMPUTEGLENSMU(CHF_FRA1(layerMu,0),
				  CHF_CONST_FRA1(epsSqr[dit],0),
				  CHF_BOX(derivBox),
				  CHF_CONST_REAL(m_n),
				  CHF_CONST_REAL(m_epsSqr0),
				  CHF_CONST_REAL(m_delta)
				  );
	    }
	}
     else if (thisA.nComp() > 1 && thisMu.nComp() == 1)
	{
	  // multiple layers, compute vertical average of mu
	  // however, I don't think this will ever be needed
	  MayDay::Error("GlensFlowRelation::computeMu : vertical average not yet implemented");
	}
#endif
      else
	{
	  MayDay::Error("GlensFlowRelation::computeMu : a_mu.nComp() != 1 or a_A.nComp()");
	}
      

    } // end loop over grids 
   
}
  
void 
GlensFlowRelation::computeDissipation(LevelData<FArrayBox>& a_dissipation,
				      const LevelData<FArrayBox>& a_vel, 
				      const LevelData<FArrayBox>* a_crseVel,
				      int a_nRefCrse,
				      const LevelData<FArrayBox>& a_A,
				      const LevelSigmaCS& a_coordSys,
				      const ProblemDomain& a_domain,
				      const IntVect& a_ghostVect) const
{
  CH_assert(a_dissipation.nComp() == a_A.nComp());

  computeMu(a_dissipation , a_vel, 1.0,   a_crseVel, a_nRefCrse, a_A, a_coordSys, 
	    a_domain, a_ghostVect);
  LevelData<FArrayBox> epsSqr;
  epsSqr.define(a_vel);
  computeStrainRateInvariant(epsSqr, a_vel, a_crseVel,
			     a_nRefCrse,a_coordSys, a_ghostVect);
  
  DataIterator dit = a_dissipation.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int layer=0; layer < a_dissipation.nComp(); ++layer)
	{
	  a_dissipation[dit].mult(epsSqr[dit],0,layer,1);
	}
      a_dissipation[dit].mult(4.0);
    }
}


void
GlensFlowRelation::computeFaceMu(LevelData<FluxBox>& a_mu,
                                 LevelData<FArrayBox>& a_vel, const Real& a_scale,
                                 const LevelData<FArrayBox>* a_crseVel, 
                                 int a_nRefCrse,
                                 const LevelData<FluxBox>& a_A, 
                                 const LevelSigmaCS& a_coordSys,
				 const ProblemDomain& a_domain,
                                 const IntVect& a_ghostVect) const
{
  // average cell-centered theta to faces
  const DisjointBoxLayout& grids = a_mu.getBoxes();
#if 0
  LevelData<FluxBox> faceTheta(grids, 1, a_ghostVect);
  CellToEdge(a_A, faceTheta);
#endif
  LevelData<FluxBox> epsSqr(grids, 1, a_ghostVect);
  computeStrainRateInvariantFace(epsSqr, a_vel, a_crseVel, a_nRefCrse, a_coordSys,
                                 a_ghostVect);

  
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box derivBox = grids[dit];
      derivBox.grow(a_ghostVect);
      FluxBox& thisMu = a_mu[dit];
      const FluxBox& thisTheta = a_A[dit];
      FluxBox& thisEpsSqr = epsSqr[dit];

      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          Box faceDerivBox = derivBox;
          faceDerivBox.surroundingNodes(faceDir);

          // note that these are the same as for the cell-centered version...
          computeMu0(thisMu[faceDir], thisTheta[faceDir], a_scale, faceDerivBox);
          FORT_COMPUTEGLENSMU(CHF_FRA1(thisMu[faceDir],0),
                              CHF_CONST_FRA1(thisEpsSqr[faceDir],0),
                              CHF_BOX(faceDerivBox),
                              CHF_CONST_REAL(m_n),
                              CHF_CONST_REAL(m_epsSqr0),
			      CHF_CONST_REAL(m_delta)
                              );
        } // end loop over face directions
    } // end loop over boxes
  
}



/// creates a new copy of this {\tt ConstitutiveRelation} object.
ConstitutiveRelation* 
GlensFlowRelation::getNewConstitutiveRelation() const
{
  GlensFlowRelation* newPtr = new GlensFlowRelation;
  
  newPtr->setParameters(m_n,  m_epsSqr0, m_delta);
                                                
  return static_cast<ConstitutiveRelation*>(newPtr);
}



// set flow parameters based on (Pattyn, 2003)
void GlensFlowRelation::setDefaultParameters()
{
  /// Power law index = 3.0
  Real n = 3.0;
  /// default rate factor 
  // small number to ensure that viscosity remains finite 
  // Pattyn (2003) uses 1e-30, but that's likely way too small
  // so instead use 1e-12 to put us more or less in the right neighborhood
  Real epsSqr0 = 1.0e-12;
  /// another regularizaton parameter, used to keep viscosity > A^{-1/n}*delta at large strain rates.
  /// originally, not present, so default is 0.0 
  Real delta = 0.0;
  setParameters(n, epsSqr0, delta); 
}

// set flow parameters 
void
GlensFlowRelation::setParameters(Real a_n, Real a_epsSqr0, Real a_delta)
{
  
  /// Power law index
  m_n = a_n;
  // small number to ensure that viscosity remains finite (1e-30)
  m_epsSqr0 = a_epsSqr0;
  //keep viscosity > A^{-1/n}*delta
  m_delta = a_delta;
}


/// utility function to compute temperature-dependent part of Glens's
/// law \mu _0 = A ^{-1/n}
void
GlensFlowRelation::computeMu0(FArrayBox& a_mu0,
                              const FArrayBox& a_A,
			      const Real& a_scale,
                              const Box& a_box) const
{
  

  a_mu0.copy(a_A);
  /// scale A (from [P]^-n /[U] to [P]^-n/[U*] 
  a_mu0 *= a_scale;
  FORT_COMPUTEGLENSMU0(CHF_FRA1(a_mu0,0),
		       CHF_BOX(a_box),
		       CHF_CONST_REAL(m_n)
		       );                        

}



constMuRelation::constMuRelation()
{
  setDefaultParameters();
}

constMuRelation::~constMuRelation()
{

}


void 
constMuRelation::computeMu(LevelData<FArrayBox>& a_mu,
                           const LevelData<FArrayBox>& a_vel, const Real& a_scale,
                           const LevelData<FArrayBox>* a_crseVel,
                           int a_nRefCrse,
                           const LevelData<FArrayBox>& a_A,
                           const LevelSigmaCS& a_coordSys,
			   const ProblemDomain& a_domain,
                           const IntVect& a_ghostVect) const
{
  DataIterator dit = a_mu.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_mu[dit].setVal(m_mu);
    }
}


void 
constMuRelation::computeDissipation(LevelData<FArrayBox>& a_dissipation,
				      const LevelData<FArrayBox>& a_vel, 
				      const LevelData<FArrayBox>* a_crseVel,
				      int a_nRefCrse,
				      const LevelData<FArrayBox>& a_A,
				      const LevelSigmaCS& a_coordSys,
				      const ProblemDomain& a_domain,
				      const IntVect& a_ghostVect) const
{

  computeStrainRateInvariant(a_dissipation, a_vel, a_crseVel,
			     a_nRefCrse,a_coordSys, a_ghostVect);
  
  DataIterator dit = a_dissipation.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_dissipation[dit] *= 4.0 * m_mu ;
    }

}

void
constMuRelation::computeFaceMu(LevelData<FluxBox>& a_mu,
                               LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                               const LevelData<FArrayBox>* a_crseVel,
                               int a_nRefCrse,
                               const LevelData<FluxBox>& a_A, 
                               const LevelSigmaCS& a_coordSys,
			       const ProblemDomain& a_domain,
                               const IntVect& a_ghostVect) const
{
  DataIterator dit=a_mu.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          a_mu[dit][dir].setVal(m_mu/a_scale);
        }
    }
}



void
constMuRelation::setDefaultParameters()
{
  // default value is mu = 1.0
  m_mu = 1.0;
}




/// creates a new copy of this {\tt ConstitutiveRelation} object.
ConstitutiveRelation* 
constMuRelation::getNewConstitutiveRelation() const
{
  constMuRelation* newPtr = new constMuRelation;
  newPtr->m_mu = m_mu;
  

  return static_cast<ConstitutiveRelation*>(newPtr);

}


void ArrheniusRateFactor::setParameters(Real a_n,
					Real a_enhance,
					Real a_B0,
					Real a_theta_r,
					Real a_K,
					Real a_C,
					Real a_R,
					Real a_Q)
{
 
  ///power law exponent
  m_n = a_n;
  /// Pattyn's "enhancement factor" = 1.0
  m_enhance = a_enhance;
  /// flow rate factor (2.207 Pa a^(1/n) )
  m_B0 = a_B0;
  /// limit temperature in flow-rate factor (273.39 K)
  m_theta_r = a_theta_r;
  /// flow rate exponent (1.17)
  m_K = a_K;
  /// flow rate factor 0.16612 K^m_K
  m_C = a_C;
  /// universal gas constant 8.31 J/(mol K) )
  m_R = a_R;
  /// activation energy for creep (7.88e4 J/mol
  // (insert funny comment about lazy creeps here...)
  m_Q = a_Q; 

}

void ArrheniusRateFactor::setDefaultParameters(Real a_seconds_per_unit_time)
{

  /// power law exponent
  Real n = 3.0;
  /// Pattyns enhancement factor
  Real enhance = 1.0 ;
  /// flow rate factor (2.207 Pa a^(1/n) )
  Real B0 = 2.207 * std::pow(SECONDS_PER_TROPICAL_YEAR/a_seconds_per_unit_time, 1.0/n);
  /// limit temperature in flow-rate factor (273.39 K)
  Real theta_r =  273.39;
  /// flow rate exponent (1.17)
  Real K = 1.17;
  /// flow rate factor 0.16612 K^m_K
  Real C = 0.16612;
  /// universal gas constant 8.31 J/(mol K) )
  Real R = 8.31;
  /// activation energy for creep (7.88e4 J/mol
  // (insert funny comment about lazy creeps here...)
  Real Q =  7.88e4;
 

  setParameters(n, enhance, B0, theta_r, K,
                C, R, Q);
}

void ArrheniusRateFactor::computeA(FArrayBox& a_A, 
				   const FArrayBox& a_thetaStar, 
				   const FArrayBox& a_pressure,
				   const Box& a_box) const
{
  //CH_assert(a_thetaStar.max(a_box) < m_theta_r);
  FORT_COMPUTEARRHENIUSA(CHF_FRA1(a_A,0),
			 CHF_CONST_FRA1(a_thetaStar,0),
			 CHF_BOX(a_box),
			 CHF_CONST_REAL(m_n),
			 CHF_CONST_REAL(m_enhance),
			 CHF_CONST_REAL(m_B0),
			 CHF_CONST_REAL(m_theta_r),
			 CHF_CONST_REAL(m_K),
			 CHF_CONST_REAL(m_C),
			 CHF_CONST_REAL(m_R),
			 CHF_CONST_REAL(m_Q)
			 );                        

 

}

ArrheniusRateFactor::ArrheniusRateFactor(Real a_seconds_per_unit_time)
{
  setDefaultParameters(a_seconds_per_unit_time);
}


RateFactor* ArrheniusRateFactor::getNewRateFactor() const
{
  ArrheniusRateFactor* newPtr = new ArrheniusRateFactor(*this);
  //newPtr->setParameters(m_n, m_enhance, m_B0, m_theta_r, m_K,
  //		       m_C, m_R, m_Q);
  return static_cast<RateFactor*>(newPtr);
  
}

PatersonRateFactor::PatersonRateFactor(Real a_seconds_per_unit_time, ParmParse& a_pp)
{
  setDefaultParameters(a_seconds_per_unit_time);
  a_pp.query("A0",m_A0);
  Real A0_multiplier = 1.0;
  a_pp.query("A0_multiplier",A0_multiplier);
  m_A0 *= A0_multiplier;
}

void PatersonRateFactor::setDefaultParameters(Real a_seconds_per_unit_time)
{
  m_E  = 1.0;
  m_A0 = 3.5e-25 * a_seconds_per_unit_time;
  m_T0 = 263.0;
  m_R  = 8.314;
  m_Qm = 6.0e+4;
  m_Qp = 1.15e+5;
}

void PatersonRateFactor::setParameters
(Real a_E, Real a_A0, Real a_T0, Real a_R, 
 Real a_Qm, Real a_Qp)
{
  m_E  = a_E;
  m_A0 = a_A0; 
  m_T0 = a_T0;
  m_R  = a_R;
  m_Qm = a_Qm;
  m_Qp = a_Qp;

}
void PatersonRateFactor::computeA
(FArrayBox& a_A, 
 const FArrayBox& a_thetaPC,
 const FArrayBox& a_pressure,
 const Box& a_box) const
{
  FArrayBox theta0PC(a_box,1);
  theta0PC.copy(a_pressure);
  theta0PC *= IceThermodynamics::icepmeltfactor();
  theta0PC += m_T0;

  FORT_COMPUTEPATERSONA
    (CHF_FRA1(a_A,0),
     CHF_CONST_FRA1(a_thetaPC,0),
     CHF_CONST_FRA1(theta0PC,0),
     CHF_BOX(a_box),
     CHF_CONST_REAL(m_E),
     CHF_CONST_REAL(m_A0),
     CHF_CONST_REAL(m_R),
     CHF_CONST_REAL(m_Qm),
     CHF_CONST_REAL(m_Qp));
  
}



RateFactor* PatersonRateFactor::getNewRateFactor() const
{
  PatersonRateFactor* newPtr = new PatersonRateFactor(*this);
  //newPtr->setParameters(m_E,m_A0,m_T0,m_R,m_Qm,m_Qp);
  return static_cast<RateFactor*>(newPtr);
}


ConstitutiveRelation* ConstitutiveRelation::parse(const char* a_prefix)
{
  ConstitutiveRelation* constRelPtr = NULL;
  ParmParse pp(a_prefix);
  std::string constRelType;
  pp.get("constitutiveRelation",constRelType);

  if (constRelType == "constMu")
      {
        constMuRelation* newPtr = new constMuRelation;
        ParmParse crPP("constMu");
        Real muVal;
        crPP.get("mu", muVal);
        newPtr->setConstVal(muVal);
        constRelPtr = static_cast<ConstitutiveRelation*>(newPtr);
      }
    else if (constRelType == "GlensLaw")
      {
        GlensFlowRelation* ptr = new GlensFlowRelation;
	ParmParse glpp("GlensLaw");
	Real n = 3.0;
	glpp.query("n",n);
	Real epsSqr0 = 1.0e-12;
	glpp.query("epsSqr0",epsSqr0);
	Real delta = 0.0;
	glpp.query("delta",delta);
	ptr->setParameters(n, epsSqr0, delta);
	constRelPtr = static_cast<ConstitutiveRelation*>(ptr);
      }
    else if (constRelType == "L1L2")
      {
        L1L2ConstitutiveRelation* ptr = new L1L2ConstitutiveRelation;
        ptr->parseParameters();
        constRelPtr = static_cast<ConstitutiveRelation*>(ptr);
      }
    else 
      {
        MayDay::Error("bad Constitutive relation type");
      }
  
  return constRelPtr;
}



RateFactor* RateFactor::parse(const char* a_prefix, bool a_basal)
{

  /// Get time units needed in rate factor constructors
  Real seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
  ParmParse ppc("constants");
  ppc.query("seconds_per_unit_time",seconds_per_unit_time);

  
  RateFactor* rateFactorPtr = NULL;

  std::string rateFactorType = a_basal ? "" : "constRate"; // default rate factor for bulk ice is constRate
  std::string rateFactorKey = a_basal ? "basalRateFactor" : "rateFactor";

  ParmParse pp(a_prefix);
  pp.query(rateFactorKey.c_str(), rateFactorType);

  if (rateFactorType == "constRate")
      {
        ParmParse crPP("constRate");
        // default value for A
        Real A = 9.2e-18 * seconds_per_unit_time/SECONDS_PER_TROPICAL_YEAR;
        crPP.query("A", A);
        ConstantRateFactor* newPtr = new ConstantRateFactor(A);

        if (a_basal)
          {
            MayDay::Error("Only PatersonRate is supported as a basal rate factor.");
          }

        rateFactorPtr = static_cast<RateFactor*>(newPtr);

        ParmParse crPP2("ConstRate");
        if (crPP2.contains("A"))
          {
            MayDay::Error("With main.rateFactor = constRate, set options using e.g. constRate.A = X (lowercase 'c'). You have used ConstRate.A = X in the input file.");
          }
      }

  else if (rateFactorType == "ConstRate")
      {
        MayDay::Error("Use main.rateFactor = constRate (lowercase 'c'). You have used main.rateFactor = ConstRate in the input file.");
      }

  else if (rateFactorType == "arrheniusRate")
      {
        /// Currently only default parameters are supported (so this additional ParmParse call, and name checking, are unnecessary). Providing code for future use if needed.
        ParmParse arPP("ArrheniusRate");
        ArrheniusRateFactor* newPtr = new ArrheniusRateFactor(seconds_per_unit_time);

        if (a_basal)
          {
            MayDay::Error("Only PatersonRate is supported as a basal rate factor.");
          }

        rateFactorPtr = static_cast<RateFactor*>(newPtr);

        ///
        ParmParse arPP2("arrheniusRate");
        if (arPP2.contains("n") ||
            arPP2.contains("enhance") ||
            arPP2.contains("B0") ||
            arPP2.contains("theta_r") ||
            arPP2.contains("K") ||
            arPP2.contains("C") ||
            arPP2.contains("R") ||
            arPP2.contains("Q"))
          {
            MayDay::Error("With main.rateFactor = arrheniusRate, set options using e.g. ArrheniusRate.n = X (uppercase 'a'). You have used arrheniusRate.n = X in the input file.");
          }
      }

  else if (rateFactorType == "ArrheniusRate")
      {
        MayDay::Error("Use  main.rateFactor = arrheniusRate (lowercase 'a'). You have used main.rateFactor = ArrheniusRate in the input file.");
      }

  else if (rateFactorType == "patersonRate")
      {
        ParmParse prPP("PatersonRate");
        PatersonRateFactor* newPtr = new PatersonRateFactor(seconds_per_unit_time, prPP);

        if (a_basal)
          {
            // Always set A0 to 1.0 as we will get the non-temperature-dependent component from the basal friction law
            newPtr->setA0(1.0);
          }

        rateFactorPtr = static_cast<RateFactor*>(newPtr);

        ParmParse prPP2("patersonRate");
        if (prPP2.contains("E") ||
            prPP2.contains("A0") ||
            prPP2.contains("T0") ||
            prPP2.contains("R") ||
            prPP2.contains("Qm") ||
            prPP2.contains("Qp") ||
            prPP2.contains("A0_multiplier"))
          {
            MayDay::Error("With main.rateFactor = patersonRate, set options using e.g. PatersonRate.E = X (uppercase 'P'). You have used patersonRate.E = X in the input file.");
          }
      }

  else if (rateFactorType == "PatersonRate")
      {
        MayDay::Error("Use main.rateFactor = patersonRate (lowercase 'p'). You have used main.rateFactor = PatersonRate in the input file.");
      }

  return rateFactorPtr;
}

#include "NamespaceFooter.H"
