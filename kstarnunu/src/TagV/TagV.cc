//
// $Id: TagV.cc 10002 2007-02-26 06:56:17Z katayama $
//
// $Log$
// Revision 1.14  2004/07/26 11:32:16  katayama
// new version from Kusaka san
//
// Revision 1.13  2003/04/08 03:36:04  sumisawa
// modify to fix OS dependence in the calculation of xi^2
//
// Revision 1.12  2002/06/04 13:48:43  sumisawa
// add functions which return modified chi2 and "setdefault" function
//
// Revision 1.11  2001/12/24 12:03:33  katayama
// gcc 3.0 and headers are cleaned up
//
// Revision 1.10  2001/12/23 09:58:26  katayama
// removed Strings.h
//
// Revision 1.9  2001/12/13 15:31:56  katayama
// MDST_OBS
//
// Revision 1.8  2001/12/09 15:45:11  katayama
// For gcc3.0
//
// Revision 1.7  2001/05/08 07:31:28  katayama
// Added float.h
//
//
#include <float.h>
#include <iostream>
#include <vector>
#include "belle.h"
#include MDST_H
#include "particle/Particle.h"
#include "particle/utility.h"
#include "kfitter/kvertexfitter.h"
// Add by T.N 2001/08/03 
extern double chisq2Confi(const int, const double);
#include "belleCLHEP/Vector/LorentzVector.h"
#include "belleCLHEP/Geometry/Point3D.h"

#include "hamlet/Hamlet.h"
#include "ip/IpProfile.h"
#include "tagv/TagV.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



// T.Nakadaira add temporary
static double
___trk_dzerr__(const Particle &p, const HepPoint3D &vertex_cp){
  double dzerr/*unused*/;
  HepSymMatrix err_from_cpvtx;
  
  Helix helix(calMdstChargedHelix(p));
  helix.pivot(vertex_cp); // move pivot to CP side vertex
  helix.x(0.0/*dPhi*/, err_from_cpvtx); // set err_from_cpvtx
  dzerr = sqrt(err_from_cpvtx[2][2]);
  return dzerr;
}

void TagV::push_back(const TagV::value_type& p)
{
  for (std::vector<Particle*>::iterator i = begin(), e = end();
       i != e; i++)
    if (*i == p) {
      dout(Debugout::ERR,"TagV") << "TagV::push_back: The same Particle* is already existing." << std::endl;
      return;
    }
  std::vector<Particle*>::push_back(p);
}

bool TagVB::isused(Particle* p)
{
  for (std::vector<Particle*>::iterator i = used_.begin(), e = used_.end();
       i != e; i++)
    if (*i == p) return 1;
  return 0;
}

//-------------------------------
int TagVB::chkKsDau(int id)
{
  int chk=0;
  double Mks=0.49767;

  // vee2 added (M.H, May 17, 2000)
  Mdst_vee2_Manager& ks2Mgr = Mdst_vee2_Manager::get_manager();
  if(ks2Mgr.count())
    {
      for (Mdst_vee2_Manager::iterator it = ks2Mgr.begin();
	   it != ks2Mgr.end(); it++) {
	Mdst_vee2& ks = *it;
	HepLorentzVector kshort(ks.px(),ks.py(),ks.pz(),ks.energy());
	if(kshort.m()>Mks - mks_veto_&&kshort.m()<Mks + mks_veto_ &&
	   (ks.chgd(0).get_ID()==id||
	    ks.chgd(1).get_ID()==id))
	  {
	    chk=1;
	    break;
	  }
      }
    }
  
  return chk;
}
//-------------------------------
void
TagVB::pivot(HepVector & helix, HepPoint3D & newPivot)
{

  m_alpha = 10000. / 2.99792458 / 15.0;
  M_PI2 = 2. * M_PI;
  M_PI4 = 4. * M_PI;
  M_PI8 = 8. * M_PI;

    double dr    = helix[0];
    double phi0  = helix[1];
    double kappa = helix[2];
    double dz    = helix[3];
    double tanl  = helix[4];
    double m_r;
    if (kappa != 0.0) {
      m_r = m_alpha / kappa;
    }
    else {
      m_r = (DBL_MAX);
    }

    double rdr = dr + m_r;
    double phi = fmod(phi0 + M_PI4, M_PI2);
    double csf0 = cos(phi);
    double snf0 = (1. - csf0) * (1. + csf0);
    snf0 = sqrt((snf0 > 0.) ? snf0 : 0.);
    if(phi > M_PI) snf0 = - snf0;

    double xc = rdr * csf0;
    double yc = rdr * snf0;
    double csf, snf;
    if(m_r != 0.0) {
      csf = (xc - newPivot.x()) / m_r;
      snf = (yc - newPivot.y()) / m_r;
      double anrm = sqrt(csf * csf + snf * snf);

      if(anrm != 0.0) {
        csf /= anrm;
        snf /= anrm;
        phi = atan2(snf, csf);
      } else {
        csf = 1.0;
        snf = 0.0;
        phi = 0.0;
      }
    } else {
      csf = 1.0;
      snf = 0.0;
      phi = 0.0;
    }
    double phid = fmod(phi - phi0 + M_PI8, M_PI2);
    if(phid > M_PI) phid = phid - M_PI2;
    double drp = (dr * csf0 + m_r * (csf0 - csf) - newPivot.x())
        * csf
        + (dr * snf0 + m_r * (snf0 - snf) - newPivot.y()) * snf;
    double dzp = dz - m_r * tanl * phid - newPivot.z();



    helix[0] = drp;
    helix[1] = fmod(phi + M_PI4, M_PI2);
    helix[2] = kappa;
    helix[3] = dzp;
    helix[4] = tanl;
}

//-------------------------------

int TagVK_impl::setdefault(class Particle& Bcp, 
		      const class HepPoint3D& Vcp,
		      const double tagvk_drcut,
		      const double tagvk_dzcut,
		      const double tagvk_errzcut,
		      const int event_by_event_ip){
  int status = 0;
  is_setdefault_called = true;

  Hamlet hmlt;
  hmlt.setBcp(Bcp, 0);
  
  if( abs(hmlt.flavor(Hamlet::HIGH_LEPTON)) == 1 ){
    const class HParticle &l(hmlt.get_plist("fbtg").front());
    this->tagl(l.mdstCharged().get_ID());
    this->useKeepTagl();
  }
  this->vcp(Vcp);
  // this->useKsVeto();
  this->dontUseKsVeto(); // modofied by T.H 2006/05/19
  this->dontUseKcRjct(); // added by T.H 2006/05/10
  this->dontUseTrkChi(); // added by T.H 2006/05/10
  this->dontUseTrkNdf(); // added by T.H 2006/05/10
  this->useTubeforFit(); // added by T.H 2006/05/19
  if(IpProfile::usable()){
//     HepPoint3D IpPosition(IpProfile::position(event_by_event_ip));
//     this->ip(IpPosition);
//     HepSymMatrix IpError(IpProfile::position_err_b_life_smeared(event_by_event_ip));
//     this->beam(IpError);
    this->ip(IpProfile::position(event_by_event_ip));
    this->beam(IpProfile::position_err_b_life_smeared(event_by_event_ip));
    // set ip and error information
    // but never use them for constraint
    this->dontUseIPforFit(); // added by T.H 2006/05/20
  }else{
    status = 1;
    this->dontUseIPforFit(); // added by T.H 2006/06/19
  }

  if(tagvk_drcut<=0.0){
    this->dontUseTrkDr();
  }else{
    this->useTrkDr();
    this->trk_dr(tagvk_drcut);
  }
  if(tagvk_dzcut<=0.0){
    this->dontUseTrkDz();
  }else{
    this->useTrkDz();
    this->trk_dz(tagvk_dzcut);
  }
  if(tagvk_errzcut<=0.0){
    this->dontUseTrkZerr();
  }else{
    this->useTrkZerr();
    this->trk_zerr(tagvk_errzcut);
  }

	use2005summer_ = false;

  return status;
}

int TagVK_impl::setdefault_2005summer(class Particle& Bcp, 
		      const class HepPoint3D& Vcp,
		      const double tagvk_drcut,
		      const double tagvk_dzcut,
		      const double tagvk_errzcut,
		      const int event_by_event_ip){
  int status = 0;
  is_setdefault_called = true;

  Hamlet hmlt;
  hmlt.setBcp(Bcp, 0);
  
  if( abs(hmlt.flavor(Hamlet::HIGH_LEPTON)) == 1 ){
    const class HParticle &l(hmlt.get_plist("fbtg").front());
    this->tagl(l.mdstCharged().get_ID());
    this->useKeepTagl();
  }
  this->vcp(Vcp);
  this->useKsVeto();
  this->dontUseKcRjct();
  this->dontUseTrkChi();
  this->dontUseTrkNdf();
  this->dontUseTubeforFit();
  if(IpProfile::usable()){
    this->ip(IpProfile::position(event_by_event_ip));
    this->beam(IpProfile::position_err_b_life_smeared(event_by_event_ip));
    this->useIPforFit(); // added by T.H 2006/05/20
  }else{
    status = 1;
  }

  if(tagvk_drcut<=0.0){
    this->dontUseTrkDr();
  }else{
    this->useTrkDr();
    this->trk_dr(tagvk_drcut);
  }
  if(tagvk_dzcut<=0.0){
    this->dontUseTrkDz();
  }else{
    this->useTrkDz();
    this->trk_dz(tagvk_dzcut);
  }
  if(tagvk_errzcut<=0.0){
    this->dontUseTrkZerr();
  }else{
    this->useTrkZerr();
    this->trk_zerr(tagvk_errzcut);
  }

	use2005summer_ = true;

  return status;
}

int TagVK_impl::selection()
{
  //  if (!ip_) return 1;

  plist_.clear();

  //TQK

  for (int j=0;j<20;j++) clist_[j]=0;

  int j=0;

  for (std::vector<Particle*>::iterator i = begin(), e = end();
       i != e; i++) {

    const Particle& p = **i;
    int mhyp = lund2mhyp(p.lund());

    if ((ip_ - p.x()).perp() <= max_perp_ &&
	p.mdstCharged().trk().mhyp(mhyp).nhits(3) >= min_nhits3_ && // Add by T.N. 2002/06
	p.mdstCharged().trk().mhyp(mhyp).nhits(4) >= min_nhits4_ ) {

      if(useKsVeto_==0||chkKsDau(p.mdstCharged().get_ID()) == 0){

	if(useKcRjct_==0||abs(p.lund()) !=321 ){

      //---------Add by TQK 2000/07/13
	  if(useTrkChi_==0||(p.mdstCharged().trk().mhyp(mhyp).chisq()/p.mdstCharged().trk().mhyp(mhyp).ndf()) < chi_cut_ ){

	    if(useTrkNdf_==0||p.mdstCharged().trk().mhyp(mhyp).ndf() > ndf_cut_ ){
	    
// 	      if(useTrkZerr_==0||
// 		 sqrt(abs(p.mdstCharged().trk().mhyp(mhyp).error(9)))
// 		< zerr_cut_ ){
//   	      T.Nakadaira change temporary
  	      if(useTrkZerr_==0||
		 ___trk_dzerr__(p, vcp_) < zerr_cut_){
	    
		//---------
		double vtagz = 0.0;
		double vtagr = 0.0;
		if(useTrkDz_==1 || useTrkDr_==1){
		  HepVector Helix(5);
		  for(int k=0;k<5;k++)
		    Helix[k] = p.mdstCharged().trk().mhyp(mhyp).helix(k);
		  
		  pivot(Helix, vcp_);
		  vtagz = Helix[3];
		  vtagr = Helix[0];
		}
		
		if((useTrkDz_==0||abs(vtagz) < dz_cut_) &&
		   (useTrkDr_==0||abs(vtagr) < dr_cut_)){

		  const Momentum& m = p.momentum();
		  const Ptype& t = p.pType();
		  
		  if(j<20){
		    clist_[j]=p.mdstCharged().get_ID();
		    j++;}
		  
		  plist_.push_back(kfitterparticle(m.p(), m.x(), m.dpx(),
						   t.charge(), t.mass()));
		  used_.push_back(*i);//Add TQK 20001122
		  //---
		}
	      }
	    }
	  }
	}
      }
    }
    //-----
    
  }
  
  return 0;
}

int TagVK_impl::selection_2005summer(void)
{
	return this->selection();
}

int TagVK_impl::selection_terse(void) // alternated by T.H 2006/05/10 (terse and improved)
{
  plist_.clear();

  //TQK

  for (int j=0;j<20;j++) clist_[j]=0;

  int j=0;

  for (std::vector<Particle*>::iterator i = begin(), e = end(); i != e; i++)
  {
    const Particle& p = **i;
    const int mhyp = lund2mhyp(p.lund());


    // if( (ip_ - p.x()).perp() > max_perp_ ) continue;

    if( p.mdstCharged().trk().mhyp(mhyp).nhits(3) < min_nhits3_ ) continue; // Add by T.N. 2002/06
    if( p.mdstCharged().trk().mhyp(mhyp).nhits(4) < min_nhits4_ ) continue;

    //---------Add by TQK 2000/07/13
    // useTrkChi_, useTrkNdf_, and useKcRjct_ are disabled in setdefault (T.H 2006/06/10)
    const double trkchi  = p.mdstCharged().trk().mhyp(mhyp).chisq();
    const int    trkndf  = p.mdstCharged().trk().mhyp(mhyp).ndf();
    const double trkRchi = trkchi / trkndf;
    if( useTrkChi_  &&  trkRchi >= chi_cut_                ) continue;
    if( useTrkNdf_  &&  trkndf <= ndf_cut_                 ) continue;

    if( useKcRjct_  &&  abs(p.lund())==321                 ) continue;


    // if not tagging lepton, further selection needed
    if( !(useKeepTagl_ && tagl_==(unsigned int)p.mdstCharged().get_ID()) )
    {
      // Tokyo cut (dzerr)
      const double trkzerr = ___trk_dzerr__(p, vcp_);
      if( useTrkZerr_ &&  trkzerr >= zerr_cut_               ) continue;


      // Tokyo cut (dz,dr)
      double vtagz = 0.0;
      double vtagr = 0.0;
      {
	HepVector Helix(5);
	for(int k=0; k<5; k++) Helix[k] = p.mdstCharged().trk().mhyp(mhyp).helix(k);

	pivot(Helix, vcp_);
	vtagz = Helix[3];
	vtagr = Helix[0];
      }

      if( useTrkDz_   &&  abs(vtagz) >= dz_cut_              ) continue;
      if( useTrkDr_   &&  abs(vtagr) >= dr_cut_              ) continue;

      // Ks veto
      if( useKsVeto_  &&  chkKsDau(p.mdstCharged().get_ID()) ) continue;
    }


    /*****************************/
    /* survived all selections!! */
    /* now add track to kfitter  */
    /*****************************/
    const Momentum& m = p.momentum();
    const Ptype& t = p.pType();

    if(j<20){
      clist_[j]=p.mdstCharged().get_ID();
      j++;
    }

    plist_.push_back(kfitterparticle(m.p(), m.x(), m.dpx(), t.charge(), t.mass()));
    used_.push_back(*i);//Add TQK 20001122
  }

  return 0;
}


// commented out by T.H 2006/05/09 -> use same ones in kvertexfitter
/*
// Add by T.N 2001/08/08
static double trackchisq(class kvertexfitter &kv, const int n,
			 const int usex, const int usey, const int usez){
  if(n<0||n>=(int)kv.tracks()){
    dout(Debugout::ERR,"TagV") << "Track # is out of range." << std::endl;
    return -1.0;
  }

  static const double hfact(1.e1);

  HepMatrix pbf = kv.track(n).getFitParameter(KF_BEFORE_FIT);
  HepMatrix paf = kv.track(n).getFitParameter(KF_AFTER_FIT);
  HepSymMatrix pbe = kv.track(n).getFitError(KF_BEFORE_FIT);
  
  int inv_error = 0;
  pbe.inverse(inv_error);
  while(inv_error && pbe.fast(6,6) < 1.e10){
    pbf *= hfact;
    paf *= hfact;
    pbe *= hfact * hfact;
    pbe.inverse(inv_error);
  }

  pbf[0][0] = 0.0;
  paf[0][0] = 0.0;
  pbf[1][0] = 0.0;
  paf[1][0] = 0.0;
  pbf[2][0] = 0.0;
  paf[2][0] = 0.0;
  if(!usex){
    pbf[3][0] = 0.0;
    paf[3][0] = 0.0;
  }
  if(!usey){
    pbf[4][0] = 0.0;
    paf[4][0] = 0.0;
  }
  if(!usez){
    pbf[5][0] = 0.0;
    paf[5][0] = 0.0;
  }
  HepMatrix dpa(pbf-paf);
  double chisq = (dpa.T()*pbe.inverse(inv_error)*dpa)[0][0];
  if(inv_error){
    dout(Debugout::ERR,"TagV") << "Can not calculate matrix inverse." << std::endl;
    return -1.0;
  }
  return chisq;
}

// Add by T.N 2001/08/08
static void chisq_wo_ip(class kvertexfitter &kv,
			double& chisq, unsigned int& ndf, double& cl,
			const int usex, const int usey, const int usez){
  chisq = 0.0;
  for(int itr=0;itr<(int)kv.tracks();++itr){
    const double i_chisq = trackchisq(kv, itr, usex, usey, usez);
    if(i_chisq<0.0){
      dout(Debugout::ERR,"TagV") << "Can not get chisq" << std::endl;
      chisq = -1.0;
      ndf = 0;
      cl = -1.0;
      return;
    }
    chisq += i_chisq;
  }
  ndf   = 2*kv.tracks();
  cl = chisq2Confi(ndf, chisq);
  return;
}
*/


int TagVK_impl::fit()
{
  // T.H
  if( !is_setdefault_called ){
    dout(Debugout::ERR,"TagV") << "Call TagVK::setdefault() prior to TagVK::fit()" << std::endl;
    exit(1);
  }

  // T.H
  if( useIP_ && useTube_ ){
    dout(Debugout::ERR,"TagV") << "TagVK::useIPforFit() and TagVK::useTubeforFit() cannot be set simultaneously" << std::endl;
    dout(Debugout::ERR,"TagV") << "Note that TagVK::setdefault() implicitly sets TagVK::useTube()" << std::endl;
    exit(1);
  }

	// T.H (2006/06/19)
	if( (useIP_ || useTube_) && !IpProfile::usable() )
	{
	  dout(Debugout::ERR,"TagV") << "TagVK: Requested IpProfile is not available for this run, return with error" << std::endl;
	  return 4;
	}

	if( use2005summer_ ){
		if (selection_2005summer()) return 3;
	} else {
		if (selection_terse()) return 3;
	}


  int size;
  int ntrk_cut(1);
  // if(useIP_) ntrk_cut = 0;
  if(useIP_||useTube_) ntrk_cut = 0; // T.H

  while ((size = plist_.size()) > ntrk_cut) {
    kvertexfitter fitter;
    if(useIP_) {
      fitter.initialVertex(ip_);
      fitter.beamProfile(beam_);
    }
    // T.H
    if(useTube_) {
      fitter.initialVertex(ip_);
      addTube2fit(fitter);
    }

    for (std::vector<kfitterparticle>::iterator i   = plist_.begin(),
	                                        end = plist_.end(); 
	 i != end; i++)
      fitter.add_track(*i);

    if (!fitter.start()) {
      chisq_ = fitter.get_chisq();

      // Modified by T.N 2001/08/16
      if ((!dontIterate_) &&
	  ((chisq_/fitter.dgf())  > max_chisq_)) { //use reduced chi**2

       	double max = 0.;
	int imax(0);

	for (int i = 0; i < size; i++) {
	  double chisq = fitter.chisq(i);

	  if(useKeepTagl_==0||clist_[i] != (int)tagl_ ){
	    
	    if (max < chisq) {
	      max = chisq;
	      imax = i;
	    }
	  }
	}
	
	plist_.erase(plist_.begin() + imax);
	used_.erase(used_.begin() + imax);//Add TQK 20001122

	for(int kk = imax;kk<19;kk++){
	  clist_[kk]=clist_[kk+1];
	}

      }
      else {

	//Add by TQK 2000/07/17
	Mdst_charged_Manager& ChgMgr = Mdst_charged_Manager::get_manager();
	if(useOneTrkIPcut_==1&&plist_.size()==1&&
 	   ((chisq_/fitter.dgf())>onetrk_ip_chi_&& //check and/or  
	    sqrt(sqr(ChgMgr((Panther_ID )clist_[0]).p(0))+	    
		 sqr(ChgMgr((Panther_ID )clist_[0]).p(1))+	    
		 sqr(ChgMgr((Panther_ID)clist_[0]).p(2)))<onetrk_ip_mom_))
	  return 3;
	//
	cl_ = fitter.get_CL();
	ntrk_ = plist_.size();
	//	vtx_ = fitter.get_vertex();
	vtx_ = fitter.vertex();
	ndf_ = fitter.dgf();
	errVtx_ = fitter.errVertex();

	// modified T.H 2006/05/09
	// // Add by T.N. 2001/08/03 --> modified 2001/08/08
	// chisq_wo_ip(fitter, chisq_woip_, ndf_woip_, cl_woip_, 0, 0, 1);
	chisq_woip_ = !use2005summer_ ? fitter.chisq_wo_ip() : fitter.chisq_wo_ip_2005summer();
	ndf_woip_   = !use2005summer_ ? fitter.dgf_wo_ip()   : fitter.dgf_wo_ip_2005summer();
	cl_woip_    = !use2005summer_ ? fitter.cl_wo_ip()    : fitter.cl_wo_ip_2005summer();
	return 0;
      }
    }
    else return 2;
  }

  return 1;
}


const bool
TagVK_impl::isTagLeptonVertex(void) const
{
  return (bool)(this->VertexTagLepton());
}

const Mdst_charged &
TagVK_impl::VertexTagLepton(void) const
{
	Mdst_charged_Manager &chgmgr = Mdst_charged_Manager::get_manager();

  if( !useKeepTagl_ ) return chgmgr.get_NULL();
  if( tagl_==0      ) return chgmgr.get_NULL(); // no tagging lepton
  if( cl_<0         ) return chgmgr.get_NULL(); // fit didn't converge

  // search for tagl in the remained (used) tracks
  // tagl would be dropped by KsVeto or SVDhits
  // i.e. any selection prior to KeepTagL
  bool found = false;
  for( std::vector<Particle*>::const_iterator u=used_.begin(); u!=used_.end(); u++ ){
    if( (unsigned int)((*u)->mdstCharged().get_ID())==tagl_ ){
      found = true;
      break;
    }
  }
  if( !found        ) return chgmgr.get_NULL(); // not found ...

  return (Mdst_charged_Manager::get_manager())[tagl_-1];
}


extern "C" int phy_vfit2_(int&, int&, int&, int[15], float[48], float[5],
			 float&, int&, float&, float[16], int&, float[3][3]);

int TagVP::selection()
{
  //  if (!ip_) return 0;

  int n = 0;

  for (std::vector<Particle*>::iterator i = begin(), e = end();
       i != e; i++) {

    const Particle& p = **i;
    int mhyp = lund2mhyp(p.lund());

    const Mdst_charged& charged = p.mdstCharged();

    if ((ip_ - p.x()).perp() <= max_perp_ &&
	charged.trk().mhyp(mhyp).nhits(3) >= min_nhits3_ &&  // Add by T.N. 2002/06
	charged.trk().mhyp(mhyp).nhits(4) >= min_nhits4_){
      
      if(useKsVeto_==0||chkKsDau(p.mdstCharged().get_ID()) == 0){
	
	if(useKcRjct_==0||abs(p.lund()) !=321 ){
	  
	  if(useKeepTagl_==0|| (unsigned)p.mdstCharged().get_ID() != tagl_ ){

	    //---------Add by TQK 2000/07/13
	    if(useTrkChi_==0||(p.mdstCharged().trk().mhyp(mhyp).chisq()/p.mdstCharged().trk().mhyp(mhyp).ndf()) < chi_cut_ ){

	      if(useTrkNdf_==0||p.mdstCharged().trk().mhyp(mhyp).ndf() > ndf_cut_ ){
		
// 		if(useTrkZerr_==0||
// 		   sqrt(abs(p.mdstCharged().trk().mhyp(mhyp).error(9)))
// 		  < zerr_cut_ ){
		
		if(useTrkZerr_==0||
		   ___trk_dzerr__(p, vcp_) < zerr_cut_){
	    
		  //---------
		  double vtagz = 0.0;
		  double vtagr = 0.0;
		  if(useTrkDz_==1 || useTrkDr_==1){
		    HepVector Helix(5);
		    for(int k=0;k<5;k++)
		      Helix[k] = p.mdstCharged().trk().mhyp(mhyp).helix(k);
		    pivot(Helix, vcp_);
		    vtagz = Helix[3];
		    vtagr = Helix[0];
		  }

		  if((useTrkDz_==0||abs(vtagz) < dz_cut_) &&
		     (useTrkDr_==0||abs(vtagr) < dr_cut_)){
		    
		    lstvtf_[n++] = charged.get_ID();
		    if (n == MAXTRK_) {
		      dout(Debugout::ERR,"TagV") << "TagVP: Too many tracks." << std::endl;
		      return 0;
		    }
		//
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    //
  }
  
  return n;
}

int TagVP::fit()
{
  int ntrkvf = selection();
  if (!ntrkvf) return 3;

  int icdc_only = 0;
  int modvtf = useIP_;
  float xvtf[48] = { 0 };
  float vprvtf[5] = {ip_.x(), ip_.y(), beam_(1,1), beam_(2,2), beam_(1,2)};

  float csqvtf;
  int ndfvtf;
  float cnfvtf;
  float csqevf[16];
  int iret;
  int ntrk_cut(1);
  if(useIP_) ntrk_cut = 0;
  float emat[3][3];

  while (ntrkvf > ntrk_cut) {
    phy_vfit2_(icdc_only, modvtf, ntrkvf, lstvtf_, xvtf, vprvtf, 
	      //          csqvtf, ndfvtf, cnfvtf, csqevf, iret);
	      csqvtf, ndfvtf, cnfvtf, csqevf, iret, emat);

    if (!iret) {
      //if (csqvtf > max_chisq_) {
      if (csqvtf/ndfvtf > max_chisq_) { //use reduced chi**2

       	float max = 0.;
	int imax(0);

	for (int i = 0; i < ntrkvf; i++) {
	  float chisq = csqevf[i];
	  if(useKeepTagl_==0||lstvtf_[i] != (int)tagl_ ){
	    if (max < chisq) {
	      max = chisq;
	      imax = i;
	    }
	  }
	}

	for (int i = 0, j = 0; i < ntrkvf; i++)
	  if (i != imax) lstvtf_[j++] = lstvtf_[i];
	ntrkvf--;
      }
      else {
	for (int i = 0; i < ntrkvf; i++) {
	  int id = lstvtf_[i];
	  for (std::vector<Particle*>::iterator i = begin(), e = end();
	       i != e; i++) {
	    if (id == (*i)->mdstCharged().get_ID())
	      used_.push_back(*i);
	  }
	}

	//Add by TQK 2000/07/17
	Mdst_charged_Manager& ChgMgr = Mdst_charged_Manager::get_manager();
	if(useOneTrkIPcut_==1&&ntrkvf==1&&
 	   ((csqvtf/ndfvtf)>onetrk_ip_chi_&& //check and/or  
	    sqrt(sqr(ChgMgr((Panther_ID )lstvtf_[0]).p(0))+	    
		 sqr(ChgMgr((Panther_ID )lstvtf_[0]).p(1))+	    
		 sqr(ChgMgr((Panther_ID)lstvtf_[0]).p(2)))<onetrk_ip_mom_))
	  return 3;
	//
	chisq_ = csqvtf;
	cl_ = cnfvtf;
	ntrk_ = ntrkvf;
	vtx_ = HepPoint3D(xvtf[0], xvtf[1], xvtf[2]);
	ndf_ = ndfvtf;

	errVtx_(1,1) = emat[0][0];
	errVtx_(1,2) = emat[0][1];
	errVtx_(1,3) = emat[0][2];
	errVtx_(2,2) = emat[1][1];
	errVtx_(2,3) = emat[1][2];
	errVtx_(3,3) = emat[2][2];
	// errVtx_(0,0) = emat[0][0];
	// errVtx_(0,1) = emat[0][1];
	// errVtx_(0,2) = emat[0][2];
	// errVtx_(1,1) = emat[1][1];
	// errVtx_(1,2) = emat[1][2];
	// errVtx_(2,2) = emat[2][2];
	return 0;
      }
    }
    else return 2;
  }

  return 1;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
