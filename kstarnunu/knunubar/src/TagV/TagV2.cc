//
// $Id: TagV2.cc 9944 2006-11-29 07:36:07Z katayama $
//
// $Log$
// Revision 1.1  2004/07/26 11:32:17  katayama
// new version from Kusaka san
//
//

#include "tagv/TagV.h"
#include <iostream>
#include <functional>
#include <algorithm>

#include "belle.h"
#include MDST_H
#include "particle/Particle.h"

#include "panther/panther.h"
#include BELLETDF_H

#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  void noticeTagV(bool useNewTagV, int expNo) {
    dout(Debugout::INFO,"TagV2")
      << "################ TagV notice ###################" << std::endl
      << __FILE__
      << ": Exp.No=" << expNo << (useNewTagV ? '>' : '<') << "30."
      << " I will use "
      << (useNewTagV ? "new" : "old")
      << " version TagV (version for SVD"
      << (useNewTagV ? '2' : '1')
      << ")." << std::endl
      << "################################################" << std::endl;
  }

  void showTagV2Version() {
    dout(Debugout::INFO,"TagV2")
      << "################ TagV notice ###################" << std::endl
      << "  July 19th, 2006 version." << std::endl
      << "################################################" << std::endl;
  }

  void showForceOldMessage() {
    dout(Debugout::INFO,"TagV2")
      << "################ TagV notice ###################" << std::endl
      << "WARNING: forceOldTagV() is called for SVD2 data!" << std::endl
      << "################################################" << std::endl;
  }

TagVK::TagVK()
  : impl_(), result_(), useNewTagV_(false)
{
  static int expNoBefore = -999;
  static bool firstTime = true;
  if(firstTime) {
    showTagV2Version();
    firstTime = false;
  }

  Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
  if(evtMgr.count() == 0){
    dout(Debugout::ERR,"TagV2")
      << "################ TagV error ####################" << std::endl
      << __FILE__
      << ": Failed to get Belle event manager."
      << " I will use old version TagV (version for SVD1)." << std::endl
      << "################################################" << std::endl;
    useNewTagV_ = false;
    return;
  }


  // Check event number and determine which method to use.
  const int expNo = evtMgr[0].ExpNo();
  useNewTagV_ = (expNo>30);

  if(expNo!=expNoBefore) {
    expNoBefore=expNo;
    noticeTagV(useNewTagV_, expNo);
  }
}


bool TagVK::lessThan2hitsIn234layer(const Particle* p) {
  const unsigned int hit_map =
    p->mdstCharged().trk().mhyp(lund2mhyp(p->lund())).hit_svd();
  const unsigned int nhits_234 =
    (hit_map&4||hit_map&8)
    + (hit_map&16||hit_map&32) + (hit_map&64||hit_map&128);
  
  return nhits_234<2;
}

bool TagVK::lessThan2hits(const Particle* p) {
  const unsigned int hit_map =
    p->mdstCharged().trk().mhyp(lund2mhyp(p->lund())).hit_svd();
  const unsigned int nhits =
    (hit_map&1||hit_map&2) + (hit_map&4||hit_map&8)
    + (hit_map&16||hit_map&32) + (hit_map&64||hit_map&128);
    
  return nhits<2;
}


int TagVK::setdefault(class Particle& Bcp,
		      const class HepPoint3D& Vcp)
{
  // if(useNewTagV_) {
  //   // Use no dr, dz, errz cut in TagV_impl
  //   return impl_.setdefault(Bcp, Vcp, -1, -1, -1);
  // }

  return impl_.setdefault(Bcp, Vcp);
}

int TagVK::setdefault(class Particle& Bcp,
		      const class HepPoint3D& Vcp,
		      const double tagvk_drcut,
		      const double tagvk_dzcut,
		      const double tagvk_errzcut,
		      const int event_by_event_ip)
{
  // if(useNewTagV_) {
  //   warnDiscardDefault();
  //   return impl_.setdefault(Bcp, Vcp, -1, -1, -1, event_by_event_ip);
  // }

  return impl_.setdefault(Bcp, Vcp,
			  tagvk_drcut, tagvk_dzcut, tagvk_errzcut,
			  event_by_event_ip);
}

int TagVK::setdefault_2005summer(class Particle& Bcp,
		      const class HepPoint3D& Vcp)
{
  return impl_.setdefault_2005summer(Bcp, Vcp);
}

int TagVK::setdefault_2005summer(class Particle& Bcp,
		      const class HepPoint3D& Vcp,
		      const double tagvk_drcut,
		      const double tagvk_dzcut,
		      const double tagvk_errzcut,
		      const int event_by_event_ip)
{
  return impl_.setdefault_2005summer(Bcp, Vcp,
			  tagvk_drcut, tagvk_dzcut, tagvk_errzcut,
			  event_by_event_ip);
}

void TagVK::push_back(std::vector<Particle*>::const_reference p) {
  particlesToUse_.push_back(p);
}


namespace {
  inline void copy_all(TagVK_impl& tagv, std::vector<Particle*>& vp) {
    for(std::vector<Particle*>::iterator it=vp.begin();it!=vp.end(); ++it)
      tagv.push_back(*it);
  }
}

int TagVK::fit() {
  if(!useNewTagV_) {
    // Old version fit
    result_ = impl_;
    copy_all(result_, particlesToUse_);
    return result_.fit();
  }


  //=========== Modified version for SVD2 ==============

  // Copy to temporary variable
  std::vector<Particle*> p_temp(particlesToUse_);

  {
    // Require 2 hits using hit pattern
    std::vector<Particle*>::iterator removed_end
      = std::remove_if(p_temp.begin(), p_temp.end(),
		       std::ptr_fun(lessThan2hits));
    p_temp.erase(removed_end, p_temp.end());
  }

  // Fit once
  result_ = impl_;
  copy_all(result_, p_temp);
  const int err = result_.fit();
  if(err) return err;

  // If the result is multi track vertex ...
  if(result_.ndf_woip()>2) return err;
  
  
  // ----- Fit again for single track vertex ------
  
  {
    // Require 2 hits in 234 layer using hit pattern
    std::vector<Particle*>::iterator removed_end
      = std::remove_if(p_temp.begin(), p_temp.end(),
  		       std::ptr_fun(lessThan2hitsIn234layer));
    p_temp.erase(removed_end, p_temp.end());
  }
  
  // Fit again
  result_ = impl_;
  copy_all(result_, p_temp);
  return result_.fit();
}

void TagVK::forceOldTagV() {
  static bool firstCall = true;
  if(firstCall && useNewTagV_) {
    showForceOldMessage();
    firstCall = false;
  }

  useNewTagV_ = false;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
