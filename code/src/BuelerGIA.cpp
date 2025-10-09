#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/* c++ implementation of Bueler et al, 2007 GIA model
 * Author: Samuel B. Kachuck
 * Date: Jan 19, 2019
 *
 */

#if FFTW_3

#include "BuelerGIA.H"
#include <fftw3.h>
#include "BisiclesF_F.H"
#include "AmrIce.H"
#include "FillFromReference.H"
#include "RefCountedPtr.H"


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


BuelerGIAFlux::BuelerGIAFlux( Real a_iceDensity, Real a_mantleDensity, Real a_gravity, Real a_waterDensity)
  : m_flex(1e23), m_visc(1e21), m_dt(0.), m_Nx(0), m_Ny(0), m_Lx(0), m_Ly(0), m_nlayers(0), m_pad(1), 
    m_isDomainSet(false), m_updatedTime(0.), m_verbosity(5), m_init(false),
    m_iceDensity(a_iceDensity), m_mantleDensity(a_mantleDensity), m_gravity(a_gravity), m_waterDensity(a_waterDensity), m_inside_box(true), m_gia_box_lox(-9999999), m_gia_box_hix(9999999), m_gia_box_loy(-9999999), m_gia_box_hiy(9999999)
{
  // need to allocate pointers for LevelData
  RefCountedPtr<LevelData<FArrayBox> > betaTmpPtr(new LevelData<FArrayBox>());
  m_beta = betaTmpPtr;

  RefCountedPtr<LevelData<FArrayBox> > tauTmpPtr(new LevelData<FArrayBox>());
  m_tau = tauTmpPtr;

  RefCountedPtr<LevelData<FArrayBox> > tafTmpPtr(new LevelData<FArrayBox>());    
  m_taf = tafTmpPtr;

  RefCountedPtr<LevelData<FArrayBox> > topoTmpPtr(new LevelData<FArrayBox>());    
  m_topo = topoTmpPtr;

  RefCountedPtr<LevelData<FArrayBox> > hTmpPtr(new LevelData<FArrayBox>());    
  m_h = hTmpPtr;

  RefCountedPtr<LevelData<FArrayBox> > loadTmpPtr(new LevelData<FArrayBox>());    
  m_load = loadTmpPtr;

  RefCountedPtr<LevelData<FArrayBox> > udotTmpPtr(new LevelData<FArrayBox>());
  m_udot = udotTmpPtr; 

  RefCountedPtr<LevelData<FArrayBox> > uTmpPtr(new LevelData<FArrayBox>());
  m_u = uTmpPtr; 

  RefCountedPtr<LevelData<FArrayBox> > tafpadhat0TmpPtr(new LevelData<FArrayBox>());
  m_tafpadhat0 = tafpadhat0TmpPtr;
    
  RefCountedPtr<LevelData<FArrayBox> > tafpadhatTmpPtr(new LevelData<FArrayBox>());
  m_tafpadhat =  tafpadhatTmpPtr;
  
  RefCountedPtr<LevelData<FArrayBox> > udotpadTmpPtr(new LevelData<FArrayBox>());
  m_udotpad = udotpadTmpPtr; 
  
  RefCountedPtr<LevelData<FArrayBox> > udotpadhatTmpPtr(new LevelData<FArrayBox>());
  m_udotpadhat = udotpadhatTmpPtr;
  
  RefCountedPtr<LevelData<FArrayBox> > upadhatTmpPtr(new LevelData<FArrayBox>());
  m_upadhat = upadhatTmpPtr; 

  RefCountedPtr<LevelData<FArrayBox> > elasTmpPtr(new LevelData<FArrayBox>());
  m_elas = elasTmpPtr;

  RefCountedPtr<LevelData<FArrayBox> > tafpadhatoldTmpPtr(new LevelData<FArrayBox>());
  m_tafpadhatold =  tafpadhatoldTmpPtr;

  RefCountedPtr<LevelData<FArrayBox> > inpadTmpPtr(new LevelData<FArrayBox>());
  m_inpad = inpadTmpPtr; 

  RefCountedPtr<LevelData<FArrayBox> > outpadhatTmpPtr(new LevelData<FArrayBox>());
  m_outpadhat = outpadhatTmpPtr; 

  RefCountedPtr<LevelData<FArrayBox> > inpadhatTmpPtr(new LevelData<FArrayBox>());
  m_inpadhat = inpadhatTmpPtr; 

  RefCountedPtr<LevelData<FArrayBox> > outpadTmpPtr(new LevelData<FArrayBox>());
  m_outpad = outpadTmpPtr; 


  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::BuelerGIAFlux()" << endl;
  }
}

/// factory method
/** return a pointerto a new SurfaceFlux object
 */
SurfaceFlux* 
BuelerGIAFlux::new_surfaceFlux()
{
  BuelerGIAFlux* newPtr = new BuelerGIAFlux(m_iceDensity, m_mantleDensity, m_gravity, m_waterDensity);

  // Copy EVERYTHING over
  newPtr->m_Nx          = m_Nx; 
  newPtr->m_Ny          = m_Ny;
  newPtr->m_Lx          = m_Lx;
  newPtr->m_Ly          = m_Ly;
  newPtr->m_isDomainSet = m_isDomainSet;
  newPtr->m_domainOffset= m_domainOffset;
  newPtr->m_inside_box  = m_inside_box;
  newPtr->m_gia_box_lox = m_gia_box_lox;
  newPtr->m_gia_box_hix = m_gia_box_hix;
  newPtr->m_gia_box_loy = m_gia_box_loy;
  newPtr->m_gia_box_hiy = m_gia_box_hiy;
  newPtr->m_nlayers     = m_nlayers;
  newPtr->m_visc        = m_visc;
  newPtr->m_viscvec     = m_viscvec;
  newPtr->m_thk         = m_thk;
  newPtr->m_flex        = m_flex;
  newPtr->m_dt          = m_dt;
  newPtr->m_updatedTime = m_updatedTime;
  newPtr->m_pad         = m_pad;  
  newPtr->fftfor        = fftfor;
  newPtr->fftinv        = fftinv;

  newPtr->m_beta        = m_beta; 
  newPtr->m_tau         = m_tau;
  newPtr->m_taf         = m_taf;
  newPtr->m_topo        = m_topo;
  newPtr->m_h           = m_h;
  newPtr->m_load        = m_load; 
  newPtr->m_tafpadhat   = m_tafpadhat;
  newPtr->m_u           = m_u;
  newPtr->m_udot        = m_udot;
  newPtr->m_udotpad     = m_udotpad;
  newPtr->m_udotpadhat  = m_udotpadhat;
  newPtr->m_upadhat     = m_upadhat;
  newPtr->m_tafpadhat0  = m_tafpadhat0;

  newPtr->m_includeElas = m_includeElas;
  newPtr->m_lame1       = m_lame1;
  newPtr->m_lame2       = m_lame2;
  newPtr->m_elas        = m_elas;
  newPtr->m_tafpadhatold= m_tafpadhatold;

  newPtr->m_oceanLoad   = m_oceanLoad;

  newPtr->m_ELRA        = m_ELRA;
  newPtr->m_ELRAtau     = m_ELRAtau;

  newPtr->m_inpad       = m_inpad;
  newPtr->m_outpadhat   = m_outpadhat;
  newPtr->m_inpadhat    = m_inpadhat;
  newPtr->m_outpad      = m_outpad;

  newPtr->m_init        = m_init;

  return static_cast<SurfaceFlux*>(newPtr);
}

BuelerGIAFlux::~BuelerGIAFlux() {
  if (m_verbosity>3) {
    pout() << "BuelerGIAFlux::~BuelerGIAFlux" << endl;
  }
  //fftw_destroy_plan(fftfor);
  //fftw_destroy_plan(fftinv);
}

// Set domain size characteristics..
void
BuelerGIAFlux::setDomain( int a_Nx, int a_Ny, Real a_Lx, Real a_Ly, 
                          IntVect& a_offset, int a_pad,
                          bool a_inside_box,
                          int a_gia_box_lox, int a_gia_box_hix,
                          int a_gia_box_loy, int a_gia_box_hiy ) {
  m_Nx = a_Nx;
  m_Ny = a_Ny;
  m_Lx = a_Lx;
  m_Ly = a_Ly;
  m_domainOffset = a_offset;
  m_pad = a_pad;
  m_inside_box = a_inside_box;
  m_gia_box_lox = a_gia_box_lox;
  m_gia_box_hix = a_gia_box_hix;
  m_gia_box_loy = a_gia_box_loy;
  m_gia_box_hiy = a_gia_box_hiy;
  m_isDomainSet = true;
  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux domain set" << m_Nx << m_Ny << m_Lx << m_Ly << endl;
  }
}

// Set 1-layer viscosity, in Pa s
void 
BuelerGIAFlux::setViscosity( Real& a_visc, int a_nlayers ) {
  m_visc = a_visc;
  m_nlayers = a_nlayers;
}

// Set 2-layer viscosity, in Pa s, from bottom to top, thickness in km
void 
BuelerGIAFlux::setViscosity( Vector<Real> a_viscvec, Real& a_thk, int a_nlayers ) {
  m_viscvec = a_viscvec;
  m_visc = a_viscvec[0];
  m_thk = a_thk;
  m_nlayers = a_nlayers;
}

// Set N-layer viscosity (N>2)
//void 
//BuelerGIAFlux::setViscosity( RealVect a_viscvec, RealVect a_thkvec ) {
//	TO BE IMPLEMENTED
//}
//

// Set the elastic properties.
void
BuelerGIAFlux::setElastic( bool a_includeElas, Real& a_lame1, Real& a_lame2 ){
  m_includeElas = a_includeElas;
  m_lame1 = a_lame1;
  m_lame2 = a_lame2;
}

// Set ELRA variables
void
BuelerGIAFlux::setELRA( bool a_ELRA, Real& a_ELRAtau ){
  m_ELRA = a_ELRA;
  if ( m_ELRA ){ m_ELRAtau = a_ELRAtau; }
}

// Method called each timestep by SurfaceFlux.
void 
BuelerGIAFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt)
{
  // Get time and check if need to update
  // assume this is being called at the end of the timestep...
  Real time = a_amrIce.time() + a_dt;
  bool needToUpdate = updateCheck(time);

  if ( needToUpdate ) {  
    // If need to update, do so. 
    if (m_verbosity > 3 ) {
      pout() << "Last upd: " << m_updatedTime << endl;
      pout() << "AMR time: " << a_amrIce.time() << endl;
      pout() << "actualDt: " << time-m_updatedTime << endl;
      pout() << "new time: " << time << endl;
    }
    updateUdot(a_amrIce, time-m_updatedTime);
    m_updatedTime = time; 
  }
  RealVect dx = a_amrIce.dx(a_level);
  RealVect m_inputFluxDx = a_amrIce.dx(0);
  FillFromReference(a_flux,
                    *m_udot,
                    dx, m_inputFluxDx,
                    m_verbosity);
}

// Compute and store values for GIA steps.
void
BuelerGIAFlux::precomputeGIAstep() {

  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::precomputeGIAstep()" << endl;
  }

  IntVect loVect = IntVect::Zero;

  // Create the box for the unpadded, data arrays.
  IntVect hiVect(m_Nx-1, m_Ny-1);
  Box domBox(loVect, hiVect);
  domBox.shift(m_domainOffset);
  Vector<Box> thisVectBox(1);
  Vector<int> procAssign(1,0);
  thisVectBox[0] = domBox;
  DisjointBoxLayout dbl(thisVectBox, procAssign);

  // Create the box for padded FFT arrays.
  IntVect hiVectPad(m_pad*m_Nx-1, m_pad*m_Ny-1);
  Box domBoxPad(loVect, hiVectPad);
  domBoxPad.shift(m_domainOffset);
  Vector<Box> thisVectBoxPad(1);
  Vector<int> procAssignPad(1,0);
  thisVectBoxPad[0] = domBoxPad;
  DisjointBoxLayout dblPad(thisVectBoxPad, procAssignPad);

  // Resize the arrays
  (*m_taf).define(dbl,1); 
  (*m_topo).define(dbl,1); 
  (*m_h).define(dbl,1); 
  (*m_load).define(dbl,1); 
  (*m_udot).define(dbl,1);
  (*m_u).define(dbl,1);

  (*m_beta).define(dblPad,1);
  (*m_tau).define(dblPad,1);
  (*m_tafpadhat0).define(dblPad,1);          
  (*m_tafpadhat).define(dblPad,1); 
  (*m_udotpad).define(dblPad,1);
  (*m_udotpadhat).define(dblPad,1); 
  (*m_upadhat).define(dblPad,1);   

  (*m_inpad).define(dblPad,1);   
  (*m_outpadhat).define(dblPad,1);   
  (*m_inpadhat).define(dblPad,1);   
  (*m_outpad).define(dblPad,1);   

  (*m_elas).define(dblPad,1);
  (*m_tafpadhatold).define(dblPad,1); 

  DataIterator ditert(dbl);
  for (ditert.begin();ditert.ok();++ditert)
  {
    (*m_taf)[ditert].setVal(0);
    (*m_topo)[ditert].setVal(0);
    (*m_h)[ditert].setVal(0);
    (*m_load)[ditert].setVal(0);
    (*m_udot)[ditert].setVal(0);
    (*m_u)[ditert].setVal(0);
  }
 
  DataIterator diter(dblPad);
  for (diter.begin();diter.ok();++diter)
  {
    (*m_beta)[diter].setVal(0);
    (*m_tau)[diter].setVal(0);
    (*m_tafpadhat0)[diter].setVal(0);        
    (*m_tafpadhat)[diter].setVal(0); 
    (*m_udotpad)[diter].setVal(0); 
    (*m_udotpadhat)[diter].setVal(0); 
    (*m_upadhat)[diter].setVal(0); 

    (*m_inpad)[diter].setVal(0); 
    (*m_outpadhat)[diter].setVal(0); 
    (*m_inpadhat)[diter].setVal(0); 
    (*m_outpad)[diter].setVal(0); 

    (*m_elas)[diter].setVal(0);
    (*m_tafpadhatold)[diter].setVal(0); 

  }

  // We use real-to-real (discrete hartley) transformations.
  // Note: FFTW in column-major order by swapping order of Nx, Ny. 
  // Note: Inverse fft needs to be normalized by N.
  DataIterator dit = dblPad.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      fftfor = fftw_plan_r2r_2d(m_Ny*m_pad, m_Nx*m_pad, (*m_inpad)[dit].dataPtr(), (*m_outpadhat)[dit].dataPtr(), 
				     FFTW_DHT, FFTW_DHT, FFTW_ESTIMATE);
      fftinv = fftw_plan_r2r_2d(m_Ny*m_pad, m_Nx*m_pad, (*m_inpadhat)[dit].dataPtr(), (*m_outpad)[dit].dataPtr(), 
				     FFTW_DHT, FFTW_DHT, FFTW_ESTIMATE);
    }

  Real kx, ky, kij, tau, alpha_l;

  dit = dblPad.dataIterator();
  for (dit.begin();dit.ok();++dit) {
    // grab this FAB
    FArrayBox& m_betaFAB = (*m_beta)[dit];
    FArrayBox& m_tauFAB = (*m_tau)[dit];
    FArrayBox& m_elasFAB = (*m_elas)[dit];

    BoxIterator bit(m_betaFAB.box());
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      int ix = iv[0]-m_domainOffset[0];
      int iy = iv[1]-m_domainOffset[1];
      int comp = 0;

      // Form the wavenumber (m_N*m_pad is the number of points)
      //                     (m_N-1)/m_L is the spacing of the points.
      kx = PI2*min(ix,m_Nx*m_pad-ix)/m_Nx/m_pad*(m_Nx-1)/m_Lx;    // m^{-1}
      ky = PI2*min(iy,m_Ny*m_pad-iy)/m_Ny/m_pad*(m_Ny-1)/m_Ly;    // m^{-1}
      if (m_verbosity>5) {
	pout() << iv[0] << '\t' << iv[1] << '\t' << ix << '\t' << iy << '\t';
      };
      kij = sqrt(pow(kx,2)+pow(ky,2));                            // m^{-1}
      // The lithosphere filter.
      alpha_l = 1. + pow(kij,4)*m_flex/m_mantleDensity/m_gravity;
      // Two layer models have r != 1.
      Real r;
      if (m_nlayers == 2) {
        Real hk = m_thk*kij;
        Real c = std::cosh(hk);
        Real s = std::sinh(hk);
        Real vv = m_viscvec[1]/m_viscvec[0];
        Real vv2 = pow(vv, 2);
        r = 2*vv*c*s + (1-vv2)*pow(hk,2) + vv2*pow(s,2) + pow(c,2);
        r /= ((vv+pow(vv,-1))*c*s + (vv-pow(vv,-1))*hk + pow(s,2) + pow(c,2));
        if (vv < 1 and r<vv){
          r = vv;
        }
        if (vv > 1 and r>vv){
          r = vv;
        }
        if (isnan(r)) {
          r = vv;
        }
      }
      else {
        r = 1.;
      }
      // The explonential viscous relaxation time constant (s)
      tau = 2*m_visc*kij/m_gravity/m_mantleDensity/alpha_l*r;
      if (m_ELRA) {
        tau = m_ELRAtau*SECSPERYEAR;
      }
      m_tauFAB(iv,comp) = tau;
      // The Bueler et al., 2007 fields.  
      m_betaFAB(iv,comp) = m_mantleDensity*m_gravity + m_flex*pow(kij,4);               // Pa/m
      if (m_includeElas) {
        m_elasFAB(iv,comp) = -(1-pow(alpha_l,-1))*m_gravity*(1./m_lame2+1./(m_lame1+m_lame2))/(2*kij); // m/Pa
      }

      if (ix == 0 && iy == 0) {
        m_elasFAB(iv,comp) = 0.;
      }
      if (m_verbosity>5) {
        pout() << kij << '\t' << r << '\t' << alpha_l << '\t' << tau << '\t' << m_betaFAB(iv,comp) << endl;
      }
    }
  }

  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::precomputeGIAstep() successful" << endl;
  }
}

void BuelerGIAFlux::init( const AmrIceBase& a_amrIce )
{
  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::init()" << endl;
  }
  if (m_init) { 
    // If initializing with a velocity field, compute and initialize the committed uplift.
    // Also assumes that the topography is not in equlibrium to start, so keeps
    // the reference load at zero load.
    computeInitialUpliftFromVelocity( a_amrIce );
  }
  else {
    // Otherwise assume the initial topography, with ice, is in equilibrium and
    // establishes current ice thickness (and water load) as reference.
    setInitialLoad(a_amrIce);
  }
}

// Set the reference thickness over flotation (tof at equilibrium uplift).
void BuelerGIAFlux::setInitialLoad( const AmrIceBase& a_amrIce ) 
{
  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::setInitialLoad()" << endl;
  }
  computeAndTransformTAF(a_amrIce);
  m_tafpadhat->copyTo(*m_tafpadhat0);
  // Save previous load for elastic correct. Suboptimal memory.
  m_tafpadhat->copyTo(*m_tafpadhatold);
}

// Set initial uplift.
void BuelerGIAFlux::setInitialUplift( LevelData<FArrayBox>& a_upl0 )
{
   if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::setInitialUplift()" << endl;
  }
  a_upl0.copyTo(*m_u);

  fftpadfor(*m_u, *m_upadhat);
}

// Set the initial velocity field (computes initial uplift with
// computeUpliftFromVelocity in AmrIce.cpp.
void BuelerGIAFlux::setInitialVelocity( LevelData<FArrayBox>& a_udot0 )
{
   if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::setInitialVelocity()" << endl;
  }
  a_udot0.copyTo(*m_udot);

  fftpadfor(*m_udot, *m_udotpadhat);

  // a REALLY ineffient (but one time) mass conservation assurance.
  DataIterator dit = (*m_udotpadhat).dataIterator();

  for (dit.begin();dit.ok();++dit) {
    // grab the relevant FAB 
    FArrayBox& m_udothatFAB = (*m_udotpadhat)[dit];
    BoxIterator bit(m_udothatFAB.box());
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      int comp = 0;
      if (iv[0]-m_domainOffset[0] == 0 && iv[1]-m_domainOffset[1] == 0) {
        m_udothatFAB(iv,comp) = 0.; 
      }
    }
  }
}

// Extract thickness above flotation from AmrIce, add ocean load if desired, 
// and transform to FFT space.
void BuelerGIAFlux::computeAndTransformTAF( const AmrIceBase& a_amrIce )
{
  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::computeAndTransformTAF()" << endl;
  }
  // extract height above flotation for each level,
  // flatten it to a single level and compute response.
  int n = a_amrIce.finestLevel() + 1;
  Vector<LevelData<FArrayBox>* > data(n);
  Vector<RealVect> amrDx(n);

  // Thickness above flotation
  for (int lev=0; lev<n; lev++) {
    data[lev] = const_cast<LevelData<FArrayBox>* >(&(a_amrIce.geometry(lev)->getThicknessOverFlotation()));
    amrDx[lev] = a_amrIce.dx(lev);
  }
  //pout() << "BuelerGIAFlux::computeAndTransformTAF() found TAF" << endl; 
  RealVect m_destDx = a_amrIce.dx(0);
  flattenCellData(*m_taf, m_destDx,data,amrDx,m_verbosity); 

  //pout() << "Collecting thickness" << endl;
  Vector<LevelData<FArrayBox>* > datathk(n);
  Vector<RealVect> amrDxthk(n);
  // Thickness
  for (int lev=0; lev<n; lev++) {
    datathk[lev] = const_cast<LevelData<FArrayBox>* >(&(a_amrIce.geometry(lev)->getH()));
  //const LevelData<FArrayBox>& H = (&(a_amrIce.geometry(lev)->getH()));
    amrDxthk[lev] = a_amrIce.dx(lev);
  }
  //pout() << "Flattening thickness" << endl;
  flattenCellData(*m_h, m_destDx,datathk,amrDxthk,m_verbosity); 
  //pout() << "Thicknes collected" << endl;

  Vector<LevelData<FArrayBox>* > dataocean(n);
  Vector<RealVect> amrDxocean(n);
  // Ocean load
  for (int lev=0; lev<n; lev++) {
    dataocean[lev] = const_cast<LevelData<FArrayBox>* >(&(a_amrIce.geometry(lev)->getTopography()));
    amrDxocean[lev] = a_amrIce.dx(lev);
  }
  flattenCellData(*m_topo, m_destDx,dataocean,amrDxocean,m_verbosity); 

  DataIterator dit = (*m_taf).dataIterator();

  for (dit.begin();dit.ok();++dit) {
    // grab the relevant FAB
    FArrayBox& m_tafFAB = (*m_taf)[dit];
    FArrayBox& m_topoFAB = (*m_topo)[dit];
    FArrayBox& m_loadFAB = (*m_load)[dit];
    FArrayBox& m_hFAB = (*m_h)[dit];

    BoxIterator bit(m_loadFAB.box());
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      int comp = 0;
      bool inbox = (iv[0] > m_gia_box_lox) && (iv[0] < m_gia_box_hix) && (iv[1] > m_gia_box_loy) && (iv[1] < m_gia_box_hiy);
      inbox = m_inside_box ? inbox : !inbox;
      if ( inbox ) {
        if (m_oceanLoad && (m_topoFAB(iv,comp) <= 0) && (m_tafFAB(iv,comp) <= 0)) { 
        // Ocean-consistent load with no grounded ice (taf<=0) is water (topo<0)
          m_loadFAB(iv, comp) = -m_waterDensity*m_topoFAB(iv,comp);
        }
        else if (m_oceanLoad && (m_tafFAB(iv,comp) > 0)) {
        // Ocean-consstent load with grounded ice (taf>0) is total ice load (h)
          m_loadFAB(iv, comp) = m_iceDensity*m_hFAB(iv,comp);
        }
        else {
        // If not considering ocean load, using thickness above flotation.
          m_loadFAB(iv, comp) = m_iceDensity*m_tafFAB(iv,comp);
        }
      }
    //pout() << m_loadFAB(iv, comp) << endl;
    }
  }
   // Forward FFT the load (padding it with zeros first)
  fftpadfor(*m_load, *m_tafpadhat);
}

// Compute the initial uplift from velocity. Set velocity with BuelerGIAFlux::setInitialVelocity.
// It works better to initialize from velocity than to input uplift field directly.
void BuelerGIAFlux::computeInitialUpliftFromVelocity(  const AmrIceBase& a_amrIce )
{
  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::computeInitialUpliftFromVelocity()" << endl;
  }

  computeAndTransformTAF(a_amrIce);
  
  DataIterator dit = (*m_beta).dataIterator();

  for (dit.begin();dit.ok();++dit) {
    // grab the relevant FAB
    FArrayBox& m_betaFAB = (*m_beta)[dit];
    FArrayBox& m_tauFAB = (*m_tau)[dit];

    FArrayBox& m_tafhatFAB = (*m_tafpadhat)[dit];
    FArrayBox& m_udothatFAB = (*m_udotpadhat)[dit];
    FArrayBox& m_uhatFAB = (*m_upadhat)[dit];

    BoxIterator bit(m_betaFAB.box());
  
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      int comp = 0;

      // Bueler formula
      m_uhatFAB(iv,comp) = -(m_gravity*m_tafhatFAB(iv,comp))/m_betaFAB(iv,comp) - m_tauFAB(iv,comp)*m_udothatFAB(iv,comp)/SECSPERYEAR;
      // Mass conservation.
      if (iv[0]-m_domainOffset[0] == 0 && iv[1]-m_domainOffset[1] == 0) {
        m_uhatFAB(iv,comp) = 0.; 
      }
    }
  }
}

// Check if updated velocities needed. For now update every step.
bool BuelerGIAFlux::updateCheck(Real time) {
  return time > m_updatedTime;
}

// Update transformed velocity and uplift fields using Bueler, et al. 2007, eq 11. 
void 
BuelerGIAFlux::updateUdot( const AmrIceBase& a_amrIce, Real a_dt ) {
  if (m_verbosity>1) {
    pout() << "BuelerGIAFlux::updateUdot()" << endl;
  }

  CH_TIME("BuelerGIAFlux::updateUdot");  
  
  computeAndTransformTAF(a_amrIce);


  DataIterator dit = (*m_beta).dataIterator();


  for (dit.begin();dit.ok();++dit) { 
    // grab the relevant FAB
    FArrayBox& m_betaFAB = (*m_beta)[dit];
    FArrayBox& m_tauFAB = (*m_tau)[dit];

    FArrayBox& m_tafhatFAB = (*m_tafpadhat)[dit];
    FArrayBox& m_tafhat0FAB = (*m_tafpadhat0)[dit];
    FArrayBox& m_udothatFAB = (*m_udotpadhat)[dit];
    FArrayBox& m_uhatFAB = (*m_upadhat)[dit];

    FArrayBox& m_elasFAB = (*m_elas)[dit];
    FArrayBox& m_tafhatoldFAB = (*m_tafpadhatold)[dit];

    BoxIterator bit(m_betaFAB.box());
  
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      int comp = 0;
      // Compute the gamma array using the timestep taken.
      Real actualGamma = pow((m_betaFAB(iv,comp)*(m_tauFAB(iv,comp) + 0.5*a_dt*SECSPERYEAR)),-1);  // m/s/Pa

      // Compute stress change from initial t=0 state, Pa
      Real dL = (m_tafhatFAB(iv,comp)-m_tafhat0FAB(iv,comp))*m_gravity;

      // m/yr
      m_udothatFAB(iv,comp) = -actualGamma*(dL + m_betaFAB(iv,comp)*m_uhatFAB(iv,comp))*SECSPERYEAR;   
      // Mass conservation.
      if (iv[0]-m_domainOffset[0] == 0 && iv[1]-m_domainOffset[1] == 0) {
        m_udothatFAB(iv,comp) = 0.; 
      }
      // Update the viscously relaxed uplift stored here for approach to
      // equilibrium.
      m_uhatFAB(iv,comp) += m_udothatFAB(iv,comp)*a_dt;

      // If elastic uplift has been requested, it is added to the velocity
      // below, but does not affect the uplift stored herein for disequilbrium.
      if (m_includeElas) {
        m_udothatFAB(iv,comp) += m_elasFAB(iv,comp)*(m_tafhatFAB(iv,comp)-m_tafhatoldFAB(iv,comp))/a_dt;
      }
    }
  }
  // Save previous load for elastic correct. Suboptimal memory.
  m_tafpadhat->copyTo(*m_tafpadhatold);
  // Reverse FFT the velocity (and crop it)
  fftinvcrop(*m_udotpadhat, *m_udot);
}

// Forward transform a_varin, using pass-through arrays.
// a_varin is unpadded, a_varouthat is padded.
void 
BuelerGIAFlux::fftpadfor (LevelData<FArrayBox>& a_varin, LevelData<FArrayBox>& a_varouthat) {
  if (m_verbosity>3) {
    pout() << "BuelerGIAFlux::fftpadfor()" << endl;
  }
  // copy the load to the padded array.
  DataIterator dit = (a_varin).dataIterator();
  DataIterator ditpad = (*m_inpad).dataIterator();
  ditpad.begin();
  for (dit.begin();dit.ok();++dit) {
    // grab this FAB
    FArrayBox& m_varinFAB = a_varin[dit];
    FArrayBox& m_inpadFAB = (*m_inpad)[ditpad];
    // Fill FFT input with zeros (just in case)
    m_inpadFAB.setVal(0.);

    // Fill data into FFT input
    BoxIterator bit(m_varinFAB.box());
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      int comp = 0;
      m_inpadFAB(iv,comp) = m_varinFAB(iv,comp); 
     } // end box loop
    ++ditpad;
  } // end data iterator loop

  // Perform the FFT
  for (ditpad.begin();ditpad.ok();++ditpad){
    fftw_execute(fftfor);
  }

  // Copy into provided FFT LevelData
  m_outpadhat->copyTo(a_varouthat);  
}

// Inverse transform a_varinhat, normalize, and crop into a_varout
void 
BuelerGIAFlux::fftinvcrop (LevelData<FArrayBox>& a_varinhat, LevelData<FArrayBox>& a_varout) {
  if (m_verbosity>3) {
    pout() << "BuelerGIAFlux::fftinvcrop()" << endl;
  }
  DataIterator dit = (*m_inpadhat).dataIterator();

  a_varinhat.copyTo(*m_inpadhat);
  for (dit.begin();dit.ok();++dit) {
    fftw_execute(fftinv);
  }

  DataIterator ditcrop = a_varout.dataIterator();
  ditcrop.begin();

  // Unpad (crop) the FFT array into the provided LevelData
  for (dit.begin();dit.ok();++dit) {
    // grab this FAB
    FArrayBox& m_outpadFAB = (*m_outpad)[dit];
    FArrayBox& m_varoutFAB = a_varout[ditcrop];

    BoxIterator bit(m_varoutFAB.box());
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      int comp = 0;
      // Normalization for inverse transform, see FFTW docs
      m_varoutFAB(iv,comp) = m_outpadFAB(iv,comp)/(m_Nx*m_Ny*m_pad*m_pad);
    } // end box loop
    ++ditcrop;
  } // end data iterator loop  
}

#include "NamespaceFooter.H"
#endif /*FFTW_3*/

