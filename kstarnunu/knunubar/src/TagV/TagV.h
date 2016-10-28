#ifndef __TAGV_H__
#define __TAGV_H__

//
// $Id: TagV.h 10002 2007-02-26 06:56:17Z katayama $
//
// $Log$
// Revision 1.7  2004/07/26 11:32:17  katayama
// new version from Kusaka san
//
// Revision 1.6  2002/07/11 17:44:46  sumisawa
// fixed bug in trk_dr() function.
//
// Revision 1.5  2002/06/04 13:48:44  sumisawa
// add functions which return modified chi2 and "setdefault" function
//
// Revision 1.4  2002/01/03 11:04:39  katayama
// Point3D and other header files are cleaned
//
// Revision 1.3  2000/07/20 12:04:31  katayama
// New (temporary???) version
//
// Revision 1.2  2000/05/18 21:23:31  katayama
// New version from Hazumi san
//
// Revision 1.1  2000/03/10 06:55:44  katayama
// New from Hazumi san
//
//

#include "belle.h"
#include <vector>

#ifndef CLHEP_POINT3D_H
#include "belleCLHEP/Geometry/Point3D.h"
#endif
#include "belleCLHEP/Matrix/SymMatrix.h"

#include "kfitter/kvertexfitter.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Particle;

class TagV : public std::vector<Particle*> {
protected:
  HepPoint3D ip_;
  HepSymMatrix beam_;
 
  HepPoint3D vtx_;
  double chisq_;
  double cl_;
  unsigned ndf_;
  unsigned ntrk_;
  //
  unsigned useIP_;
  unsigned useTube_; // added by T.H 2006/05/19
  unsigned use2005summer_; // added by T.H 2006/06/16
  //
  unsigned tagl_;
  HepSymMatrix errVtx_;
  unsigned useKsVeto_;
  unsigned useKcRjct_;
  unsigned useKeepTagl_;
  unsigned useTrkChi_;
  unsigned useTrkNdf_;
  unsigned useTrkZerr_;
  unsigned useTrkDz_;
  unsigned useTrkDr_; // Add by T.N. 2002/06
  unsigned useOneTrkIPcut_;

  HepPoint3D vcp_;

  // Add by T.N 2001/08/03
  double   chisq_woip_;
  double   cl_woip_;
  unsigned ndf_woip_;
  // Add by T.N 2001/08/16
  unsigned dontIterate_;

  double m_alpha,M_PI2,M_PI4,M_PI8;


public:
  //TagV() : ip_(0.,0.,0.), beam_(3,0), useIP_(0) {}
  TagV() :
    ip_(0.,0.,0.), beam_(3,0)
    // Add by T.H 2006/05/09 (to see if fit converges or not)
    ,cl_(-1)
    ,useIP_(0)
    ,useTube_(0)
    ,use2005summer_(0) // Add by T.H 2006/06/16
    // Add by T.H 2006/05/09 (tagl_=0 indicates no Hamlet info (or noKeepTagl)
    ,tagl_(0)
    ,errVtx_(3,0), useKsVeto_(0), useKcRjct_(0), useKeepTagl_(0), useTrkChi_(0), useTrkNdf_(0), useTrkZerr_(0), useTrkDz_(0), useTrkDr_(1), useOneTrkIPcut_(0), vcp_(0., 0., 0.)
    // Add by T.N 2001/08/03
    ,chisq_woip_(0.0), cl_woip_(0.0),ndf_woip_(0)
    // Add by T.N 2001/08/16
    ,dontIterate_(0)
  {}
  //Modified by TQK 00/07/17

  virtual ~TagV() {}
  virtual int fit() = 0;

  void push_back(const TagV::value_type& p);

  void ip(const class HepPoint3D* ip) { ip_ = *ip; }
  void ip(const class HepPoint3D& ip) { ip_ = ip; }
  void beam(const class HepSymMatrix* beam) { beam_ = *beam; useIP_ = 1;}
  void beam(const class HepSymMatrix& beam) { beam_ =  beam; useIP_ = 1;}
  void useIPforFit() {useIP_ = 1;}
  void dontUseIPforFit() {useIP_ = 0;}
  // added by T.H 2006/05/19
  void useTubeforFit() {useTube_ = 1;}
  void dontUseTubeforFit() {useTube_ = 0;}


  //Add by TQK 00/05/15
  void useKsVeto() {useKsVeto_ = 1;}
  void dontUseKsVeto() {useKsVeto_ = 0;}
  void useKcRjct() { useKcRjct_ = 1;}
  void dontUseKcRjct() {useKcRjct_ = 0;}
  void useKeepTagl() {useKeepTagl_ = 1;}
  void dontUseKeepTagl() {useKeepTagl_ = 0;}

  //Add by TQK 00/07/13
  void useTrkChi() { useTrkChi_ = 1;}
  void dontUseTrkChi() {useTrkChi_ = 0;}
  void useTrkNdf() { useTrkNdf_ = 1;}
  void dontUseTrkNdf() {useTrkNdf_ = 0;}
  void useTrkZerr() {useTrkZerr_ = 1;}
  void dontUseTrkZerr() {useTrkZerr_ = 0;}

  //Add by TQK 00/07/17
  void useTrkDz() {useTrkDz_ = 1;}
  void dontUseTrkDz() {useTrkDz_ = 0;}
 
  // Add by T.N. 2002/06
  void useTrkDr() {useTrkDr_ = 1;}
  void dontUseTrkDr() {useTrkDr_ = 0;}

  void useOneTrkIPcut() {useOneTrkIPcut_ = 1; useIP_= 1;}
  void dontUseOneTrkIPcut() {useOneTrkIPcut_ = 0; useIP_= 0;}

  void tagl(const int lchg) {tagl_ = lchg;}

  //  
  const HepPoint3D& ip() const { return ip_; }
  const HepSymMatrix& beam() const { return beam_; }

  const HepPoint3D& vtx() const { return vtx_; }
  double chisq() const { return chisq_; }
  double cl() const { return cl_; }
  double ndf() const { return ndf_; }
  double ntrk() const { return ntrk_; }
  // Add by T.N 2001/08/03
  // Now these 3 functions are not supported by TagVP.
  double chisq_woip() const { return chisq_woip_; }
  double cl_woip()    const { return cl_woip_; }
  double ndf_woip()   const { return ndf_woip_; }
  //
  const HepSymMatrix& errVtx() const { return errVtx_; }
  //
  // Add by T.N 2001/08/16
  void doIterate() { dontIterate_ = 0; };
  void dontIterate(){ dontIterate_ = 1; };
};

class TagVB : public TagV {
protected:
  std::vector<Particle*> used_;

  double mks_veto_;  
  double max_perp_;
  double max_chisq_;
  int min_nhits3_; // Add by T.N. 2002/06
  int min_nhits4_;
  double chi_cut_;
  int ndf_cut_;
  double zerr_cut_;
  double dz_cut_;
  double dr_cut_; // Add by T.N. 2002/06
  double onetrk_ip_chi_;
  double onetrk_ip_mom_;

public:
  //TagVB() : max_perp_(.1), max_chisq_(50), min_nhits4_(2) {}
  TagVB() : mks_veto_(0.015), max_perp_(.1), max_chisq_(20), min_nhits3_(1), min_nhits4_(2), chi_cut_(10.), ndf_cut_(20), zerr_cut_(.1), dz_cut_(2.), dr_cut_(.05), onetrk_ip_chi_(3.), onetrk_ip_mom_(.5) {}
  //Change threshold and add parameters by TQK 00/07/13

  virtual ~TagVB() {}

  void max_perp(const double a) { max_perp_ = a; }
  void max_chisq(const double a) { max_chisq_ = a; }
  void min_nhits3(const int a) { min_nhits3_ = a;} // Add by T.N. 2002/06
  void min_nhits4(const int a) { min_nhits4_ = a;}
  void mks_veto(const double a) { mks_veto_ = a; }
  void trk_chi(const double a) { chi_cut_ = a; }
  void trk_ndf(const int a) { ndf_cut_ = a; }
  void trk_zerr(const double a) { zerr_cut_ = a; }
  void trk_dz(const double a) { dz_cut_ = a; }
  void trk_dr(const double a) { dr_cut_ = a; }

  void vcp(const class HepPoint3D* vcp) { vcp_ = *vcp; }
  void vcp(const class HepPoint3D& vcp) { vcp_ = vcp; }

  void onetrk_ip_chi(const double a) { onetrk_ip_chi_ = a; }
  void onetrk_ip_mom(const double a) { onetrk_ip_mom_ = a; }

  const std::vector<Particle*>& used_particles() const { return used_; }
  bool isused(Particle*);
  //
  int chkKsDau(int);
//  HepVector & pivot(HepVector &, HepPoint3D &);
  void pivot(HepVector &, HepPoint3D &);
  //
};


/**
 * The implementation of TagVK.
 */
class TagVK_impl : public TagVB {
  std::vector<kfitterparticle> plist_;

  int clist_[20];

  bool is_setdefault_called;

public:
  TagVK_impl() : is_setdefault_called(false) {}
  virtual ~TagVK_impl() {}
  virtual int fit();
  int chgd(const int i) const { return clist_[i]; }

protected:
  virtual int selection();

// Add by T.H  2006/05/09
protected:
  virtual int selection_2005summer(void);
  virtual int selection_terse(void);
public:
	const bool isTagLeptonVertex(void) const;
	const Mdst_charged &VertexTagLepton(void) const;

public:
  // Add by T.N. 2002/06
  virtual int setdefault(class Particle& Bcp, 
			 const class HepPoint3D& Vcp,
			 const double tagvk_drcut = 0.05,
			 const double tagvk_dzcut = -0.18,
			 const double tagvk_errzcut = 0.05,
			 const int event_by_event_ip = 1); 
  // Add by T.H 2006/06/16 for backward compatibility
  virtual int setdefault_2005summer(class Particle& Bcp, 
			 const class HepPoint3D& Vcp,
			 const double tagvk_drcut = 0.05,
			 const double tagvk_dzcut = -0.18,
			 const double tagvk_errzcut = 0.05,
			 const int event_by_event_ip = 1); 
};


/**
 * The wrapper class of TagVK_impl
 * @note This class relies on implicit copy constructor of TagVK_impl.
 */
class TagVK {
private:
  TagVK_impl impl_;
  TagVK_impl result_;

  bool useNewTagV_;
  std::vector<Particle*> particlesToUse_;

  /**
   * @return True if nhits in 234 layer is less than 2. (Using hit pattern)
   */
  static bool lessThan2hitsIn234layer(const Particle* p);
  /**
   * @return True if nhits is less than 2. (Using hit pattern)
   */
  static bool lessThan2hits(const Particle* p);

public:
  // The new methods in Jul. 8th, 2004
  TagVK();
  /// Force to use old version (SVD1) TagV, even for SVD2 data.
  void forceOldTagV();
  int setdefault(class Particle& Bcp,
		 const class HepPoint3D& Vcp);
  int setdefault_2005summer(class Particle& Bcp,
		 const class HepPoint3D& Vcp);
  void push_back(std::vector<Particle*>::const_reference p);

  const std::vector<Particle*>& particlesToUse() const
    { return particlesToUse_; }
  std::vector<Particle*>& particlesToUse()
    { return particlesToUse_; }

  // Method from TagVK
  int setdefault(class Particle& Bcp,
		 const class HepPoint3D& Vcp,
		 const double tagvk_drcut/* = 0.05 */,
		 const double tagvk_dzcut = -0.18,
		 const double tagvk_errzcut = 0.05,
		 const int event_by_event_ip = 1);
  int setdefault_2005summer(class Particle& Bcp,
		 const class HepPoint3D& Vcp,
		 const double tagvk_drcut/* = 0.05 */,
		 const double tagvk_dzcut = -0.18,
		 const double tagvk_errzcut = 0.05,
		 const int event_by_event_ip = 1);
  int fit();

  // // Methods from TagVB
  // void max_perp(const double a)  { impl_.max_perp(a); }
  // void max_chisq(const double a) { impl_.max_chisq(a); }
  // void min_nhits3(const int a)   { impl_.min_nhits3(a); }
  // void min_nhits4(const int a)   { impl_.min_nhits4(a);}
  // void mks_veto(const double a)  { impl_.mks_veto(a); }
  // void trk_chi(const double a)   { impl_.trk_chi(a); }
  // void trk_ndf(const int a)      { impl_.trk_ndf(a); }
  void trk_zerr(const double a)  { impl_.trk_zerr(a); }
  void trk_dz(const double a)    { impl_.trk_dz(a); }
  void trk_dr(const double a)    { impl_.trk_dr(a); }

  void vcp(const class HepPoint3D* vcp) { impl_.vcp(vcp); }
  void vcp(const class HepPoint3D& vcp) { impl_.vcp(vcp); }

  void onetrk_ip_chi(const double a) { impl_.onetrk_ip_chi(a); }
  void onetrk_ip_mom(const double a) { impl_.onetrk_ip_mom(a); }

  // int chkKsDau(int);
  // void pivot(HepVector&, HepPoint3D&);


  // Methods from TagV
  void ip(const class HepPoint3D* ip) { impl_.ip(ip); }
  void ip(const class HepPoint3D& ip) { impl_.ip(ip); }
  void beam(const class HepSymMatrix* beam) { impl_.beam(beam);}
  void beam(const class HepSymMatrix& beam) { impl_.beam(beam);}
  void useIPforFit() { impl_.useIPforFit();}
  void dontUseIPforFit() { impl_.dontUseIPforFit();}
  // added by T.H 2006/05/19
  void useTubeforFit() { impl_.useTubeforFit();}
  void dontUseTubeforFit() { impl_.dontUseTubeforFit();}

  // modified by T.H 2006/05/19
  void useKsVeto()       { fprintf(stderr, "TagVK: useKsVeto is not recommended\nTagVK: setdefault() turns off this flag\n"); impl_.useKsVeto();}
  void dontUseKsVeto()   { impl_.dontUseKsVeto();}
  // commented out by T.H 2006/05/19
  // void useKcRjct()       { impl_.useKcRjct();}
  // void dontUseKcRjct()   { impl_.dontUseKcRjct();}
  void useKeepTagl()     { impl_.useKeepTagl();}
  void dontUseKeepTagl() { impl_.dontUseKeepTagl();}

  // void useTrkChi()          { impl_.useTrkChi(); }
  // void dontUseTrkChi()      { impl_.dontUseTrkChi(); }
  // void useTrkNdf()          { impl_.useTrkNdf(); }
  // void dontUseTrkNdf()      { impl_.dontUseTrkNdf(); }
  void useTrkZerr()         { impl_.useTrkZerr(); }
  void dontUseTrkZerr()     { impl_.dontUseTrkZerr(); }
  void useTrkDz()           { impl_.useTrkDz(); }
  void dontUseTrkDz()       { impl_.dontUseTrkDz(); }
  void useTrkDr()           { impl_.useTrkDr(); }
  void dontUseTrkDr()       { impl_.dontUseTrkDr(); }
  // void useOneTrkIPcut()     { impl_.useOneTrkIPcut(); }
  // void dontUseOneTrkIPcut() { impl_.dontUseOneTrkIPcut(); }

  void tagl(const int lchg) { impl_.tagl(lchg);}

  const HepPoint3D& ip() const { return impl_.ip(); }
  const HepSymMatrix& beam() const { return impl_.beam(); }


  // Accesses to the result.
  const std::vector<Particle*>& used_particles() const { return result_.used_particles(); }
  bool isused(Particle* p) { return result_.isused(p); };
  const HepPoint3D& vtx() const { return result_.vtx(); }
  double chisq() const          { return result_.chisq(); }
  double cl() const             { return result_.cl(); }
  double ndf() const            { return result_.ndf(); }
  double ntrk() const           { return result_.ntrk(); }
  double chisq_woip() const     { return result_.chisq_woip(); }
  double cl_woip()    const     { return result_.cl_woip(); }
  double ndf_woip()   const     { return result_.ndf_woip(); }
  const HepSymMatrix& errVtx() const { return result_.errVtx(); }
  const bool isTagLeptonVertex(void) const { return result_.isTagLeptonVertex();}      // T.H 206/05/09
  const Mdst_charged &VertexTagLepton(void) const { return result_.VertexTagLepton();} // T.H 206/05/09

  void doIterate() { impl_.doIterate(); }
  void dontIterate(){ impl_.doIterate(); }
};


class TagVP : public TagVB {
  static const int MAXTRK_ = 15;
  int lstvtf_[MAXTRK_];

public:
  TagVP() {}
  virtual ~TagVP() {}

  virtual int fit();
  int chgd(const int i) const { return lstvtf_[i]; }

protected:
  virtual int selection();
};

inline int lund2mhyp(const int lund)
{ 
  switch(std::abs(lund)) {
  case  11: return 0;
  case  13: return 1;
  case 321: return 3;
  default:  return 2;
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __TagV_H__ */

