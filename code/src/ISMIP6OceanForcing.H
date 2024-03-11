#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// concrete class encapsulating surface fluxes determined  
// by ISMIP6 ocean forcing.


// \todo replace the std::map<Real,string> mechanism with
// a suitable abstraction

#ifndef _ISMIP6OCEANFORCING_H_
#define _ISMIP6OCEANFORCING_H_

#include "SurfaceFlux.H"
#include "ParmParse.H"
#include <map>
#include "NamespaceHeader.H"


/**

   Concrete class encapsulating surface fluxes determined  
   by ISMIP6 ocean forcing and similar.

   It supports similar quadratic melt rate formula from two papers
   Jourdain et al 2020, https://tc.copernicus.org/articles/14/3111/2020/
   Burgard et al 2022 https://tc.copernicus.org/articles/16/4931/2022/t

   They require 3D thermal forcing data Tf(x,y,z) and (optionally) 2D corrections dT.
   It may also require basin maps.

   The class computes Tf = Tf(x,y,z = z_ice_shelf_lower_surface) for use im the
   formulae and <Tf>, which is the mean <Tf> over ice shelves in eaqch basin 
   
   Burgard e,g eq 14 does not include dT but does include salinity. The user needs
   to supply Tf(x,y,z) = sqrt(S_ocean)*(T_ocean - T_freeze) and set dT = 0. T_ocean
   is the far field ocean temperature, S_ocean the salinity.

   The thermal forcing is either 'local' or 'semilocal'

   local forcing:  melt = a(x,y) * Abs(Tf - dT) * (Tf + dT)
   semi-local forcing: melt = a(x,y) * Abs(<Tf> - dT) * (Tf + dT)

   a(x,y) can be a constant (e.g in the 'Antarctic slope' cases of Burgard), or a constant multiplied
   by the slope ('local slope') case of Burgard. 

   i.e a(x,y) = f or a(x,y) = f |grad(lower surface)|

   The user specifies factor (which should be negative) via ParmParse. 
   If it is not specified, the default is computed according
   to Jourdain 2020, with either the use-specified or defaul values of 'gamma0'


   Factor example: Burgard 2022 eq 14/18

   melt = K \frac{rhoo}{rhoi} (Co/Li)^2 \beta_S^2 \frac{g}{2|f|} \sin \theta |Tf|Tf

   melt is measured in kg/s. \sin \theta is the slope, either local, or a global constant. 
   Burgard gives a table 11 of K, we need 

   factor = - K \frac{\rho_o}{\rhoi_} (Co/Li)^2 \beta_S^2 \frac{g}{2|f|} * \rho_i * seconds_per_year

   Quadratic local antarctic slope: \sin \theta = 2.9e-3, median K = 12.0 e-5 -> factor = -0.0332
   Quadratic local local cavity slope: median K = 5.5e-5 -> factor = 5.246904 (slope computed)

   The melt rate can be applied as an anomaly, in which case the melt rate computed from the 
   current ice shelf geometry and forcing data of the start (start_year)
   is subtracted from the melt rate computed from current geoemtry and forcing data.


 */
class ISMIP6OceanForcing : public SurfaceFlux
{
public:

  
  ISMIP6OceanForcing(ParmParse& a_pp);

  /// destructor
  virtual ~ISMIP6OceanForcing() {};

  /// factory method
  /** return a pointer to a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** 
  */

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);


protected:

  std::map<int,std::string> m_year_file;
  int m_start_year, m_end_year;
  std::string m_name;
  std::string m_basin_file, m_basin_var_name;
  bool m_deltaT_specified;
  std::string m_deltaT_file,  m_deltaT_var_name;
  bool m_anomaly, m_local, m_local_slope;
  Real m_factor, m_basin_mean_min_thickness;
  bool m_tf_is_source;
  

  //- 1.0 * melt rate on a uniform mesh
  RefCountedPtr<LevelData<FArrayBox> > m_uniform_source;
  // mesh spacing of m_uniform_source
  RealVect m_dx;
  // timestamp for
  int m_uniform_source_year;

  // layer spacing of m_TF (between components)
  Real m_dz;
  int m_n_layer;

  // basin masks: mask(x,y,n) = basin_number(x,y) == n ? 1 : 0  
  RefCountedPtr<LevelData<FArrayBox> > m_basin_mask;
  int m_n_basin; // number of basins
  
  // deltaT(x,y)
  RefCountedPtr<LevelData<FArrayBox> > m_deltaT; // no-need to keep reading this
  
  /// compute source from thermal forcing
  void computeSource(LevelData<FArrayBox>& a_source,
		     LevelData<FArrayBox>& a_TFb,
		     LevelData<FArrayBox>& a_TFb_basin_mean,
		     LevelData<FArrayBox>& a_deltaT,
		     Real a_factor,
		     LevelData<FArrayBox>& a_slope);

  /// compute thermal forcing TFb at the ice shelf base from the 3D 'layered' field TF
  void computeTFb_and_slope(LevelData<FArrayBox>& a_TFb,
			    const LevelData<FArrayBox>& a_TF,
			    LevelData<FArrayBox>& a_slope,
			    const AmrIceBase& a_amrIce,
			    bool compute_slope);


  /// compute the basin means of TFb, store in TFb_mean
  void computeBasinMeans(LevelData<FArrayBox>&a_TFb_basin_mean,
			 LevelData<FArrayBox>&a_TFb,
			 LevelData<FArrayBox>&a_basin_mask,
			 const AmrIceBase& a_amrIce);
  
  /// update the uniform source
  void updateUniformSource(Real a_time, const AmrIceBase& a_amrIce);

  ///  Read a field TF(x,y,z) and compute source(x,y) 
  void readUniformSource ( RefCountedPtr<LevelData<FArrayBox> >& a_source,
			   RealVect& a_dx,
			   int a_year,
			   const AmrIceBase& a_amrIce);
    
   
};


#include "NamespaceFooter.H"

#endif
