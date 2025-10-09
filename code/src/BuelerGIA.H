#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if FFTW_3

#ifndef _BUELERGIA_H_
#define _BUELERGIA_H_
#include "SurfaceFlux.H"
#include "RefCountedPtr.H"
#include <fftw3.h>

#include "NamespaceHeader.H"

/* Input file ParmParse structure:
 * topographyFlux.type = buelerGIA
 * topographyFlux.nlayers = 1
 * topographyFlux.flex = 1e23				// in N m
 * topographyFlux.visc = 1e18				// in Pa s
 * topographyFlux.pad = 2				// multiply domain for padding
 *
 * topographyFlux.layers = 2
 * topographyFlux.visc = 2e19 4e18			// Bottom layer to top layer
 * topographyFlux.thk  = 200				// in km
 *
 * topographyFlux.init = true
 * topographyFlux.init.file = /path/to/init.hdf5	// must match the domain
 * topographyFlux.init.name = udot0
 */

class BuelerGIAFlux : public SurfaceFlux
{
public:
  BuelerGIAFlux( Real a_iceDensity, Real a_mantleDensity, Real a_gravity, Real a_waterDensity );

  virtual ~BuelerGIAFlux();
 

  /// factory method
  /** return a pointerto a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  // The BISICLES flux interface function
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);


protected:
  const Real PI2=2.*3.14159267;
  const Real SECSPERYEAR=31536000.;

  int m_nlayers;
  Real m_flex, m_visc, m_thk;
  Vector<Real> m_viscvec, m_thkvec;
  Real m_lame1, m_lame2;
  Real m_dt;
  // physical constants
  Real m_iceDensity, m_gravity, m_mantleDensity, m_waterDensity;
  // Domain constants
  int m_Nx, m_Ny;
  Real m_Lx, m_Ly;
  IntVect m_domainOffset;
  bool m_inside_box;
  int m_gia_box_lox, m_gia_box_hix, m_gia_box_loy, m_gia_box_hiy;
  bool m_isDomainSet, m_init, m_oceanLoad, m_includeElas;
  bool m_ELRA;
  Real m_ELRAtau;
  int m_verbosity;
  Real m_updatedTime;
  int m_pad;
  // Constants for GIA step (computed once, stored) (size pad*Nx x pad*Ny)
  RefCountedPtr<LevelData<FArrayBox>> m_beta, m_elas, m_tau;
  // Padded quantities required during computation (size pad*Nx x pad*Ny)
  RefCountedPtr<LevelData<FArrayBox>> m_tafpadhat;   // Thickness above flotation and FFT'.
  RefCountedPtr<LevelData<FArrayBox>> m_udotpad, m_udotpadhat; // Surface velocity and FFT'd.
  // Pass-through arrays for the FFT and inverse FFT (size pad*Nxx pad*Ny)
  RefCountedPtr<LevelData<FArrayBox>> m_inpad, m_outpadhat;
  RefCountedPtr<LevelData<FArrayBox>> m_inpadhat, m_outpad;
  // FFTW transformations.
  fftw_plan fftfor, fftinv;

protected:
  // Check if updated velocities needed. For now update every step.
  bool updateCheck( Real t );
  // Update transformed velocity and uplift fields using Bueler, et al. 2007, eq 11. 
  void updateUdot( const AmrIceBase& a_amrIce, Real time );
  // Extract thickness above flotation from AmrIce, add ocean load if desired, 
  // and transform to FFT space.
  void computeAndTransformTAF( const AmrIceBase& a_amrIce );

  // Forward transform a_varin, using pass-through arrays.
  // a_varin is unpadded, a_varouthat is padded.
  void fftpadfor(LevelData<FArrayBox>& a_varin, LevelData<FArrayBox>& a_varouthat);
  // Inverse transform a_varinhat, normalize, and crop into a_varout
  void fftinvcrop(LevelData<FArrayBox>& a_varinhat, LevelData<FArrayBox>& a_varout);
  
public:
  /* set-up functions */
  // Set 1-layer viscosity
  void setViscosity( Real& a_visc, int a_nlayers );
  // Set 2-layer viscosity, from bottom to top, thickness in km
  void setViscosity( Vector<Real> a_viscvec, Real& a_thk, int a_nlayers );
  // Set N-layer viscosity (N>2)
  //void setViscosity( RealVect a_viscvec, RealVect a_thkvec ); 
  // Set the flexural rigidity (N m) of the lithosphere.
  void setFlexural( Real& a_flex ) {  m_flex = a_flex; };
  // Set the timestep of the GIA iteration. Not used currently.
  void setTimestep( Real& a_dt ) {m_dt = a_dt;};
  // Set domain size characteristics
  void setDomain( int a_Nx, int a_Ny, Real a_Lx, Real a_Ly, IntVect& a_offset,
                  int a_pad, bool a_inside_box, int a_gia_box_lox, int a_gia_box_hix,
                          int a_gia_box_loy, int a_gia_box_hiy );
  // Set the elastic properties.
  void setElastic( bool a_includeElas, Real& a_lame1, Real& a_lame2 );
  // Set ELRA
  void setELRA( bool a_ELRA, Real& a_ELRAtau );
  // Set Ocean Load bool
  void setOceanLoad( bool a_oceanLoad ) {m_oceanLoad = a_oceanLoad;};
  // Compute and store values for GIA steps 
  void precomputeGIAstep();
  void init( const AmrIceBase& a_amrIce );
  // Set the initial thickness above flotation (tof at equilibrium uplift)..
  void setInitialLoad( const AmrIceBase& a_amrIce ) ;
  // Set the initial uplift field.
  void setInitialUplift( LevelData<FArrayBox>& a_upl0 );
  // Set the initial velocity field (computes initial uplift with
  // computeUpliftFromVelocity in AmrIce.cpp.
  void setInitialVelocity( LevelData<FArrayBox>& a_upl0 );
  // Compute the initial uplift from velocity. Set velocity with BuelerGIAFlux::setInitialVelocity.
  // It works better to initialize from velocity than to input uplift field directly.
  void computeInitialUpliftFromVelocity( const AmrIceBase& a_amrIce ) ;

  /* Checkpoint restart functions */
  void setTAF0( LevelData<FArrayBox>& a_tafpadhat0 ) { a_tafpadhat0.copyTo(*m_tafpadhat0); };
  void setTAFold( LevelData<FArrayBox>& a_tafpadhatold ) { a_tafpadhatold.copyTo(*m_tafpadhatold); };
  void setUn( LevelData<FArrayBox>& a_upadhat ) { a_upadhat.copyTo(*m_upadhat); };
  void setUdot( LevelData<FArrayBox>& a_udot ) { a_udot.copyTo(*m_udot); };
  void setUpdatedTime( Real t ) { m_updatedTime = t; };
  void setInitIceRef0( bool a_init ) {m_init=a_init;};

  RefCountedPtr<LevelData<FArrayBox>> m_tafpadhat0, m_upadhat, m_tafpadhatold; // Surface velocity and FFT'd.
  //RefCountedPtr<LevelData<FArrayBox>> m_udot, m_taf, m_u, m_topo, m_load; // Surface velocity and FFT'd.
  RefCountedPtr<LevelData<FArrayBox>> m_udot, m_taf, m_u, m_topo, m_h, m_load; // Surface velocity and FFT'd.

  /* Checkpoint write functions */
  RefCountedPtr<LevelData<FArrayBox>> getTAF0() { return m_tafpadhat0; };
  RefCountedPtr<LevelData<FArrayBox>> getTAFold() { return m_tafpadhatold; };
  RefCountedPtr<LevelData<FArrayBox>> getUn() { return m_upadhat; };
  RefCountedPtr<LevelData<FArrayBox>> getUdot() { return m_udot; };
  RefCountedPtr<LevelData<FArrayBox>> getLoad() { return m_load; };
  Real getUpdatedTime() { return m_updatedTime; };
  

  bool isInitIceRef0() { return m_init; };
};

#include "NamespaceFooter.H"


#endif /*_BUELERGIA_H_*/
#else /*FFTW_3*/
#warning ('FFTW_3 != TRUE: BuelerGIA not compiled')
#endif /*FFTW_3*/

