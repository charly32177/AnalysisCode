#if !defined(FINDKS_H)
#define FINDKS_H

#include "belle.h"
#ifndef CLHEP_POINT3D_H
#include "belleCLHEP/Geometry/Point3D.h"
#endif
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Mdst_charged;
class Mdst_vee2;
class Mdst_vee_daughters;

class FindKs 
{
public:
  FindKs() {};
  virtual ~FindKs() {};

  void candidates(const Mdst_vee2& cand, HepPoint3D RunIP);
  int goodKs(); // 1: good
  double dr() { return _dr; }  //dr before vtx-fit
  double dphi() { return _dphi; }
  double zdist() { return _zdist; }
  double fl() {return _fl; }
  double pmag() {return _pmag; }
  double chisq() {return _chisq; }

private:
  HepPoint3D _IP;
  double _dr;
  double _dphi;
  double _fl;
  double _zdist;
  double _pmag;
  double _chisq;
  int _kind;

  // find dr of daughters if decays inside beampipe.
  double _impactParameter(const Mdst_vee_daughters& dau, 
			 HepPoint3D RunIP, 
			 HepPoint3D vtx,
			 int charge);
  double _impactParameter(const Mdst_charged& charged,
			  HepPoint3D RunIP);
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif
