#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DomainDiagnosticData.H"
#include "AmrIce.H"
#include "IceConstants.H"
#include "computeSum.H"
#include <functional>
#include "NamespaceHeader.H"



void  MaskedIntegration::integrateScalarInside
(Real& a_integral, 
 function<bool(Real h, Real f, int mask)> a_inside,
 const Vector<RefCountedPtr<LevelSigmaCS > >& a_coords,
 const Vector<LevelData<FArrayBox>* >& a_integrand, 
 const Vector<LevelData<FArrayBox>* >& a_topography,
 const Vector<LevelData<FArrayBox>* >& a_thickness,
 const Vector<LevelData<FArrayBox>* >& a_iceFrac,
 const Vector<LevelData<FArrayBox>* >& a_sectorMaskFraction,
 const Vector<Real>& a_dx, const Vector<int>& a_ratio,
 int a_finestLevel, 
 int a_maskNo, int a_maskComp)
{
  CH_TIME("MaskedIntegration::integrateScalarInside");
  
  int numLevels = a_finestLevel + 1;
  Vector<LevelData<FArrayBox>* > maskedIntegrand(numLevels, NULL);
  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& grids = a_integrand[lev]->disjointBoxLayout();
      maskedIntegrand[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  const BaseFab<int>& mask = a_coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*a_thickness[lev])[dit];
	  const FArrayBox& frac = (*a_iceFrac[lev])[dit];
	  const Box& b = grids[dit];
	  const FArrayBox& f = (*a_integrand[lev])[dit];
	  FArrayBox& g = (*maskedIntegrand[lev])[dit];
	  g.setVal(0.0);
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      // fractional area if cell within the sector mask...
	      //safe if sectorMaskFraction[lev] is NULL because maskNo == 0
	      Real sf = (a_maskNo == 0)?1.0:(*a_sectorMaskFraction[lev])[dit](iv, a_maskComp);
	      if (a_inside(thck(iv), frac(iv), mask(iv)))
		{
		  g(iv) = f(iv)*sf;
		}
	    } // bit
	} //dit
    } //lev

  a_integral = computeSum(maskedIntegrand, a_ratio, a_dx[0], Interval(0,0), 0);

  for (int lev = 0; lev < numLevels; lev++)
    {
      if (maskedIntegrand[lev] != NULL)
	{
	  delete maskedIntegrand[lev];maskedIntegrand[lev]= NULL;
	}  
    }
  
}

void MaskedIntegration::integrateDischargeInside
(Real& a_sumDischarge, Real& a_sumDivUH,
 function<bool(Real h, Real f, int mask)> a_inside,
 const Vector<RefCountedPtr<LevelSigmaCS > >& a_coords,
 const Vector<LevelData<FluxBox>* >& a_fluxOfIce,
 const Vector<LevelData<FArrayBox>* >& a_topography,
 const Vector<LevelData<FArrayBox>* >& a_thickness,
 const Vector<LevelData<FArrayBox>* >& a_iceFrac,
 const Vector<Real>& a_dx, const Vector<int>& a_ratio, 
 const Vector<LevelData<FArrayBox>* >& a_sectorMaskFraction,
 int a_finestLevel, 
 int a_maskNo, int a_maskComp)
{

  CH_TIME("integrateDischargeInside");
  
  int numLevels = a_finestLevel + 1;
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDivUH(numLevels, NULL);
  for (int lev = 0; lev < numLevels; lev++)
    {
      int comp = 0;
      const DisjointBoxLayout& grids = a_topography[lev]->disjointBoxLayout();
      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDivUH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
     
      //work out discharge across boundary and flux divergence inside region 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
	  FArrayBox& divuh = (*ccDivUH[lev])[dit];
	  divuh.setVal(0.0);
	  const BaseFab<int>& mask = a_coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*a_thickness[lev])[dit];
	  const FArrayBox& frac = (*a_iceFrac[lev])[dit];
    	  const Box& b = grids[dit];
	  
    	  for (int dir =0; dir < SpaceDim; dir++)
    	    {
	      const FArrayBox& flux  = (*a_fluxOfIce[lev])[dit][dir];
    	      for (BoxIterator bit(b);bit.ok();++bit)
    		{
    		  const IntVect& iv = bit();
		  //safe if sectorMaskFraction[lev] is NULL because maskNo == 0
    		  if (a_maskNo == 0 || ((*a_sectorMaskFraction[lev])[dit](iv, a_maskComp) > 0.5))
		    {
		      if (a_inside(thck(iv), frac(iv), mask(iv)))
		       	{
			  //inside the sub-region - compute div(flux)
		       	  divuh(iv) += (flux(iv + BASISV(dir)) - flux(iv))/a_dx[lev];
		       	}
		      else
    			{
			  //outside - compute discharge/dx if neigbours are inside
			  IntVect ivp = iv + BASISV(dir);
			  if (a_inside(thck(ivp), frac(ivp), mask(ivp)))
			    {
			      discharge(iv) += -flux(ivp) / a_dx[lev]; // ivp is the iv+ face
    			    }
			  IntVect ivm = iv - BASISV(dir);
			  if (a_inside(thck(ivm), frac(ivm), mask(ivm)))
    			    {
    			      discharge(iv) += flux(iv) / a_dx[lev]; // iv is the iv- face
    			    }
    			}//end if inside else
    		    }//end if mask
    		} //end loop over cells
    	    } // end loop over direction
    	} // end loop over grids
    } //end loop over levels

  a_sumDischarge = computeSum(ccDischarge, a_ratio, a_dx[0], Interval(0,0), 0);
  a_sumDivUH = computeSum(ccDivUH, a_ratio, a_dx[0], Interval(0,0), 0);
  
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


/// set up a struct of cf names to storage
void DomainDiagnosticData::setCFdata()
{

  cfDiagnostic cf_info;

  cf_info.short_name = CFIO_DIAGNOSTIC_TIME_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_TIME_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_TIME_NAME;
  cf_info.units = "yr";
  cf_info.data = &m_time;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_ICE_VOLUME_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_ICE_VOLUME_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_ICE_VOLUME_LONG_NAME;
  cf_info.units = "kg";
  cf_info.data = &m_ice_volume;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_ICE_VAF_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_ICE_VAF_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_ICE_VAF_LONG_NAME;
  cf_info.units = "kg";
  cf_info.data = &m_ice_vaf;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_LONG_NAME;
  cf_info.units = "m^2";
  cf_info.data = &m_ice_grounded_area;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_LONG_NAME;
  cf_info.units = "m^2";
  cf_info.data = &m_ice_floating_area;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_total_smb;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_total_bmb;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_floating_total_bmb;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_CALVING_FLUX_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_CALVING_FLUX_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_CALVING_FLUX_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_total_calving_flux;
  m_cf_stuff.push_back(cf_info);

}


void
DomainDiagnosticData::initDiagnostics
(AmrIce& a_amrIce,
 const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_refRatio, Real a_crseDx,
 Real a_time, int a_finestLevel)
{

  m_report_grounded_ice = false;
  m_report_area = false;
  m_report_total_flux = false;
  m_report_calving = false;

  ParmParse pp("amr");

  pp.query("report_sum_grounded_ice", m_report_grounded_ice);
  pp.query("report_ice_area", m_report_area);
  pp.query("report_total_flux", m_report_total_flux);
  pp.query("report_calving", m_report_calving);
  record(a_amrIce);
}


void DomainDiagnosticData::record(AmrIce& a_amrIce)
{

  Real dt = 1.2345678e+300;
  size_t ithis = m_time.size();
  size_t iprev = ithis - 1;
  if (ithis > 0)
    {
      dt = a_amrIce.time() - m_time[iprev]; 
    }
  Real rhoi = a_amrIce.m_iceDensity;


  const Vector<RefCountedPtr<LevelSigmaCS> >& coords = a_amrIce.amrGeometry();
  const Vector<int>& ratio = a_amrIce.refRatios();
  const Vector<Real>& dx = a_amrIce.amrDx();
  int n = a_amrIce.finestLevel() + 1;
  Vector<LevelData<FArrayBox>* > topg(n), thk(n), hab(n), sectorMaskFraction(n);
  for (int lev = 0; lev < n; lev++)
    {
      hab[lev] = (const_cast<LevelData<FArrayBox>* >(&(coords[lev]->getThicknessOverFlotation())));
      thk[lev] = (const_cast<LevelData<FArrayBox>* >(&(coords[lev]->getH())));
      topg[lev] = (const_cast<LevelData<FArrayBox>* >(&(coords[lev]->getTopography())));
    }
  Vector<LevelData<FArrayBox>* >& smb = a_amrIce.m_surfaceThicknessSource;
  Vector<LevelData<FArrayBox>* >& bmb = a_amrIce.m_basalThicknessSource;
  Vector<LevelData<FArrayBox>* >& cflux = a_amrIce.m_calvedIceThickness;
  Vector<LevelData<FArrayBox>* >& frac = a_amrIce.m_iceFrac;
  
  int maskNo = 0; int maskComp = 0; // might support this later
  Real h_min = TINY_THICKNESS;

  // subregions
  auto entire = [h_min](Real h, Real f, int mask){ return true;} ;
  auto grounded = [h_min](Real h, Real f, int mask){ return ((mask == GROUNDEDMASKVAL) && (h > h_min));} ;
  auto floating = [h_min](Real h, Real f, int mask){ return ((mask == FLOATINGMASKVAL) && (h > h_min));} ;
  auto ice = [h_min](Real h, Real f, int mask){ return (h > h_min);} ;

  // capture repetitive args
  auto sumScalar = [coords, topg, thk, hab, frac,
		    sectorMaskFraction, dx, ratio, n, maskNo, maskComp]
    (const Vector<LevelData<FArrayBox>* >& scalar, function<bool(Real h, Real f, int mask)> inside)
		   {
		     Real s;
		     MaskedIntegration::integrateScalarInside
		       (s, inside, coords, scalar, topg, thk,
			frac, sectorMaskFraction, dx, ratio, n-1, maskNo, maskComp);
		     return s;
		   };

  m_time.push_back(a_amrIce.time());
  m_ice_volume.push_back(rhoi*sumScalar(thk, ice));
  m_ice_vaf.push_back(rhoi*sumScalar(hab, ice));
  m_ice_grounded_area.push_back(sumScalar(frac, grounded));
  m_ice_floating_area.push_back(sumScalar(frac, floating));
  m_ice_total_smb.push_back(rhoi*sumScalar(smb, ice));
  m_ice_total_bmb.push_back(rhoi*sumScalar(bmb, ice));
  m_ice_floating_total_bmb.push_back(rhoi*sumScalar(bmb, floating));
  m_ice_total_calving_flux.push_back(rhoi*sumScalar(cflux, entire));

}

inline 
void last_to_first_resize(Vector<Real>& a_v)
{
  if (a_v.size() > 0)
    {
      a_v[0] = a_v[a_v.size()-1];
      a_v.resize(1);
    }
}


void DomainDiagnosticData::reset()
{

 if (m_time.size() > 0)
   {
     /// copy the last items from the previous set of records
     last_to_first_resize(m_time);
     last_to_first_resize(m_ice_volume);
     last_to_first_resize(m_ice_vaf);
     last_to_first_resize(m_ice_grounded_area);
     last_to_first_resize(m_ice_floating_area);
     last_to_first_resize(m_ice_total_smb);
     last_to_first_resize(m_ice_total_bmb);
     last_to_first_resize(m_ice_floating_total_bmb);
     last_to_first_resize(m_ice_total_calving_flux);
   }
 else
   {
     m_time.resize(0);
     m_ice_volume.resize(0);
     m_ice_vaf.resize(0);
     m_ice_grounded_area.resize(0);
     m_ice_floating_area.resize(0);
     m_ice_total_smb.resize(0);
     m_ice_total_bmb.resize(0);
     m_ice_floating_total_bmb.resize(0);
     m_ice_total_calving_flux.resize(0);
   }
}

#ifdef CH_USE_HDF5

void readV(HDF5Handle& a_handle, Vector<Real>& a_data,
	   const std::string& a_name)
{

  H5E_auto_t efunc; void* edata;
  // need these to turn auto error messaging off then back 
  
#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  hid_t set = H5Dopen(a_handle.groupID(), a_name.c_str());
  H5Eset_auto(efunc, edata);
#else
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  hid_t set = H5Dopen2(a_handle.groupID(), a_name.c_str(), H5P_DEFAULT);
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);
#endif
  if (set >= 0)
    {
      hid_t space = H5Dget_space(set);
      if (space >= 0)
	{
	  hsize_t dim,maxdim;
	  H5Sget_simple_extent_dims(space, &dim, &maxdim);
	  a_data.resize(dim);
	  readDataset(set, space, &a_data[0], 0,  a_data.size());
	  H5Sclose(space);
	}
      H5Dclose(set);
    }
}

void readAttr(hid_t loc_id, const std::string& a_name, std::string& a_value)
{
  herr_t ret = 0;
  char messg[1024];

#ifdef H516
  hid_t attr   = H5Aopen_name(loc_id, a_name.c_str());
#else
  hid_t attr   = H5Aopen_by_name(loc_id,".", a_name.c_str(), H5P_DEFAULT, H5P_DEFAULT );
#endif
  hid_t atype  = H5Aget_type(attr);
  hid_t aclass = H5Tget_class(atype);
  char* buf = NULL;  size_t size = 0;

  size = H5Tget_size(atype);
  buf = new char[size+1];
  ret = H5Aread(attr, atype, buf);
  if (ret < 0) 
    {
      sprintf(messg,"DomainDiagnosticData:: Problem reading attribute %s",a_name.c_str());
      MayDay::Warning(messg);
    }

  buf[size] = 0; // for some reason HDF5 is not null terminating strings correctly
  a_value = std::string(buf);

  delete[] buf;
  H5Tclose(atype);
  H5Aclose(attr);

}

void readStruct(HDF5Handle& a_handle, cfDiagnostic& a_cf_info)
{

  H5E_auto_t efunc; void* edata;
  // need these to turn auto error messaging off then back 

  Vector<Real>& data = *a_cf_info.data;
  std:: string data_name = a_cf_info.short_name;
  
#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  hid_t set = H5Dopen(a_handle.groupID(), data_name.c_str());
  H5Eset_auto(efunc, edata);
#else
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  hid_t set = H5Dopen2(a_handle.groupID(), data_name.c_str(), H5P_DEFAULT);
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);
#endif
  if (set >= 0)
    {
      hid_t space = H5Dget_space(set);
      if (space >= 0)
	{
	  hsize_t dim,maxdim;
	  H5Sget_simple_extent_dims(space, &dim, &maxdim);
	  data.resize(dim);
	  readDataset(set, space, &data[0], 0,  data.size());
	  H5Sclose(space);
	}

      H5Dclose(set);
    }
}

void writeV(HDF5Handle& a_handle, Vector<Real>& a_data,
	    const std::string& a_name)
{
  hid_t set,space,type;
  createDataset(set, space, a_handle, a_name, &a_data[0], a_data.size());
  if (procID() == uniqueProc(SerialTask::compute))
    {
      writeDataset(set, space, &a_data[0], 0,  a_data.size());
    }
  H5Sclose(space);
  H5Dclose(set);
}

void writeAttr(hid_t loc_id, const std::string& a_name, const std::string& a_value)
{
  H5E_auto_t efunc; void* edata;
#ifdef H516
  H5Eget_auto(&efunc, &edata);
#else
  H5Eget_auto2(H5E_DEFAULT, &efunc, &edata);
#endif
  herr_t ret = 0;
  char messg[1024];

  hid_t s_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(s_type, a_value.length());
  hid_t aid  = H5Screate(H5S_SCALAR);
#ifdef H516
  H5Eset_auto(NULL, NULL);
  hid_t attr = H5Acreate(loc_id, a_name.c_str(), s_type,
			 aid, H5P_DEFAULT);
#else
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  hid_t attr = H5Acreate2(loc_id, a_name.c_str(), s_type,
			  aid, H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (attr < 0)
    {
      H5Adelete(loc_id, a_name.c_str());
#ifdef H516
      attr = H5Acreate(loc_id, a_name.c_str(), s_type,
		       aid, H5P_DEFAULT);
#else
      attr = H5Acreate2(loc_id, a_name.c_str(), s_type,
			aid, H5P_DEFAULT,H5P_DEFAULT);
#endif
      if (attr < 0)
	{
	  sprintf(messg,"DomainDiagnosticData:: Problem creating attribute %s",a_name.c_str());
	  MayDay::Warning(messg);
	}
    }
#ifdef H516
  H5Eset_auto(efunc, edata);
#else
  H5Eset_auto2(H5E_DEFAULT, efunc, edata);
#endif
  char* tmp = (char*)a_value.c_str();
  ret = H5Awrite(attr, s_type, tmp);
  if (ret < 0) 
    {
      sprintf(messg,"DomainDiagnosticData:: Problem writing attribute %s",a_name.c_str());
      MayDay::Warning(messg);
    }

  H5Sclose(aid);
  H5Aclose(attr);
  H5Tclose(s_type);

}

void writeStruct(HDF5Handle& a_handle, const cfDiagnostic& a_cf_info)
{
  
  Vector<Real>& data = *a_cf_info.data;
  if (data.size() > 0)
    {
      hid_t set,space,type;
      std:: string data_name = a_cf_info.short_name;
      createDataset(set, space, a_handle, data_name, &data[0], data.size());
      if (procID() == uniqueProc(SerialTask::compute))
	{
	  writeDataset(set, space, &data[0], 0,  data.size());
	}
      writeAttr(set, "Short name", a_cf_info.short_name);
      writeAttr(set, "Long name", a_cf_info.long_name);
      writeAttr(set, "Units", a_cf_info.units);
      writeAttr(set, "Standard name", a_cf_info.cf_name);
      H5Sclose(space);
      H5Dclose(set);	
    }
}

void DomainDiagnosticData::write(HDF5Handle& a_handle)
{
      if (a_handle.pushGroup(HDF5_SUBGROUP_NAME) == 0)
	{     
	  for (int i = 0; i < m_cf_stuff.size(); ++i)
	    {
	      cfDiagnostic cf_info = m_cf_stuff[i];
	      writeStruct(a_handle,cf_info);
	    }
	  a_handle.popGroup();
	}
    
}

void DomainDiagnosticData::read(HDF5Handle& a_handle)
{
  if (a_handle.pushGroup(HDF5_SUBGROUP_NAME) == 0)
    {
      for (int i = 0; i < m_cf_stuff.size(); ++i)
	{
	  cfDiagnostic cf_info = m_cf_stuff[i];
	  readStruct(a_handle,cf_info);
	}
      a_handle.popGroup();
    }
}

DomainDiagnosticData::DomainDiagnosticData()
{
  setCFdata();
}

DomainDiagnosticData::DomainDiagnosticData(const DomainDiagnosticData& a)
{
  setCFdata();
}

DomainDiagnosticData&  DomainDiagnosticData::operator=(const DomainDiagnosticData& a)
{
  setCFdata();
  return *this;
}



#endif
#include "NamespaceFooter.H"
