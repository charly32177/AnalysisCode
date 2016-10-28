//*************************************************************
//makeKs change to Mdst_vee2, supported by b20000426 version
//
//*************************************************************
#include "particle/utility.h"
#include "./UserUtility.h"
#include "./exUserInfo.h"
#include "./gammac.h"
#include "mdst/mdst.h"
#include HEPEVT_H
#include MDST_H

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

/**************************************************************
Utility Funtion Source Code
***************************************************************/
/*
void
setUserInfo(Particle &p)
{
  p.userInfo(*(new UserInfo));
}

void
setUserInfo(std::vector<Particle> &p)
{
  for(int i=0;i<p.size();++i)setUserInfo(p[i]);
}
*/
/*
void setErr(Particle &gamm, HepPoint3D &vtx, HepSymMatrix &err){ 
 if(gamm.mdstGamma()){
      HepSymMatrix errG(3,0);
      errG[0][0] = gamm.mdstGamma().ecl().error(0);
      errG[1][1] = gamm.mdstGamma().ecl().error(2);
      errG[2][2] = gamm.mdstGamma().ecl().error(5);
      const HepSymMatrix errG1 = errG;
      GammaParticle test(double(gamm.mdstGamma().ecl().energy()),
                         double(gamm.mdstGamma().ecl().phi()),
                         double(gamm.mdstGamma().ecl().theta()),
			 errG1 );// 3x3 matrix

      test.vertex(vtx,err);
      gamm.momentum().momentum(test.momentumEnergy(),
                               test.errorMomentumEnergy());
      gamm.momentum().position(vtx,err);
 }else{
  cout << "This is not well reconstruted gamma, charge ="
       <<  gamm.charge() <<endl;
 }
}
*/
void 
makeKP(std::vector<Particle> &k_p, 
	std::vector<Particle> &k_m, 
	std::vector<Particle> &pi_p, 
	std::vector<Particle> &pi_m)
{
  Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
  
  //...Particle Type
  Ptype ptype_pion_plus("PI+");
  Ptype ptype_pion_minus("PI-");
  Ptype ptype_kaon_plus("K+");
  Ptype ptype_kaon_minus("K-");
 
  //...Fills pion and kaon lists with MDST_Charged Data Base
  for(std::vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); ++i){
//    if(!good_charged(*i))continue;
    if((*i).trk() && (*i).trk().quality()==0){
    if((*i).charge() > 0.){
      Particle tmp1(*i, ptype_kaon_plus);
      k_p.push_back(tmp1);
      Particle tmp2(*i, ptype_pion_plus);
      pi_p.push_back(tmp2);
    }else{
      Particle tmp1(*i, ptype_kaon_minus);
      k_m.push_back(tmp1);
      Particle tmp2(*i, ptype_pion_minus);
      pi_m.push_back(tmp2);
    }
   }
  }
}
/*
float helicity(Vector4& baseP4,Vector4& momP4,Vector4& dauP4){
        momP4.boost(-(baseP4.boostVector()));
        dauP4.boost(-(baseP4.boostVector()));
        dauP4.boost(-(momP4.boostVector()));
        Vector3 momP3 = momP4;
        Vector3 dauP3 = dauP4;
        if(momP3.mag()*dauP3.mag() >0)
        return (momP3*dauP3)/(momP3.mag()*dauP3.mag());
        else return -2;
}
*/

/********************************************************************
Do Mass fit
********************************************************************/
/*
void doKmFit(std::vector<Particle> &plist) {
  for(unsigned i=0;i<plist.size();++i)
    doKmFit(plist[i]);
}

unsigned doKmFit(Particle &p) {
  kmassfitter km;
  km.invariantMass(p.pType().mass());
  for(unsigned j=0;j<p.relation().nChildren();++j){
  km.addTrack(p.relation().child(j).momentum().p(),
              p.relation().child(j).momentum().x(),
              p.relation().child(j).momentum().dpx(),
              p.relation().child(j).pType().charge(),
              p.relation().child(j).pType().mass());
    
  }
  unsigned err = km.fit();

  setUserInfo(p);
  dynamic_cast<UserInfo&>(p.userInfo()).chisq(-7.);
  dynamic_cast<UserInfo&>(p.userInfo()).cl(-7.);
  dynamic_cast<UserInfo&>(p.userInfo()).dgf(-7.);

  if(err)return 0; //err!=0 -> failed return 0
  // kmassfitter always set dgf ==1
  dynamic_cast<UserInfo&>(p.userInfo()).chisq(km.chisq());
  dynamic_cast<UserInfo&>(p.userInfo()).cl(km.cl());
  dynamic_cast<UserInfo&>(p.userInfo()).dgf(km.dgf());
  return makeMother(km,p);
}
*/
/********************************************************************
Do kvertex fit
********************************************************************/
/*
void doKvFit(std::vector<Particle> &plist) {
  for(unsigned i=0;i<plist.size();++i)
    doKvFit(plist[i]);
}
unsigned doKvFit(Particle &p,HepPoint3D beam, HepSymMatrix errBeam) {
  kvertexfitter kv;
  kv.initialVertex(beam);
  kv.beamProfile(errBeam);

  for(unsigned j=0;j<p.relation().nChildren();++j){
    addTrack2fit(kv,p.relation().child(j));
  }
  unsigned err = kv.fit();

  setUserInfo(p);
  dynamic_cast<UserInfo&>(p.userInfo()).chisq(-7.);
  dynamic_cast<UserInfo&>(p.userInfo()).cl(-7.);
  dynamic_cast<UserInfo&>(p.userInfo()).dgf(-7.);

  if(err)return 0; //err!=0 -> failed return 0
  // kmassfitter always set dgf ==1
  dynamic_cast<UserInfo&>(p.userInfo()).chisq(kv.chisq());
  dynamic_cast<UserInfo&>(p.userInfo()).cl(kv.cl());
  dynamic_cast<UserInfo&>(p.userInfo()).dgf(kv.dgf());

  return makeMother(kv,p);
}

unsigned doKvFit(Particle &p) {
  kvertexfitter kv;
  for(unsigned j=0;j<p.relation().nChildren();++j){
    addTrack2fit(kv,p.relation().child(j));
  }
  unsigned err = kv.fit();

  setUserInfo(p);
  dynamic_cast<UserInfo&>(p.userInfo()).chisq(-7.);
  dynamic_cast<UserInfo&>(p.userInfo()).cl(-7.);
  dynamic_cast<UserInfo&>(p.userInfo()).dgf(-7.);

  if(err)return 0; //err!=0 -> failed return 0
  // kmassfitter always set dgf ==1
  dynamic_cast<UserInfo&>(p.userInfo()).chisq(kv.chisq());
  dynamic_cast<UserInfo&>(p.userInfo()).cl(kv.cl());
  dynamic_cast<UserInfo&>(p.userInfo()).dgf(kv.dgf());

  return makeMother(kv,p);
}
/*
/********************************************************************
Do vertex and mass fit
********************************************************************/
/*
void doKmvFit(std::vector<Particle> &plist) {
  for(unsigned i=0;i<plist.size();++i)
    doKmvFit(plist[i]);
}

unsigned doKmvFit(Particle &p) {
  kmassvertexfitter kmv;
  kmv.invariantMass(p.pType().mass());
  for(unsigned j=0;j<p.relation().nChildren();++j){
    addTrack2fit(kmv,p.relation().child(j));
  }
  unsigned err = kmv.fit();
  if(err)return 0;
  return makeMother(kmv,p);
}
unsigned
makeMother(kmassfitter &kv,
           Particle &mother)
{
  unsigned n = kv.tracks(); 
  kmakemother kmm;               // convert helix -> mom pos ??
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(kv.momentum(i),
                 kv.position(i),
                 kv.error(i),
                 mother.relation().child(i).pType().charge());
    mother.relation().child(i).momentum().momentumPosition(kv.momentum(i),
                                                                 kv.position(i),
                                                                 kv.error(i));

    //mother.relation().child(i).fittedMomentum().cl(kv.cl());
    //mother.relation().child(i).fittedMomentum().chisq(kv.chisq());
    //mother.relation().child(i).fittedMomentum().dof(kv.dgf());
  }
  for(unsigned i=0;i<n;++i){
    for(unsigned j=i+1;j<n;++j){
      HepMatrix tmp1(kv.correlation(j,i));
//      mother.relation().child(j).fittedMomentum().coMatrix(&mother.relation().child(i),tmp1);
    }
  }
  unsigned err = kmm.make();
  if(err != 0)return 0;
  mother.momentum().momentumPosition(kmm.momentum(),
				     kmm.position(),
				     kmm.error());
  return 1;
}
*/
int checkMultiUse(Particle& i, Particle& j){
// nFinalStateParticles()==1 => FinalStateParticles = itself  
  if(i.relation().nFinalStateParticles()>1){       
    for(int l=0;l<i.relation().nFinalStateParticles();l++){
      Particle Final(i.relation().finalStateParticle(l));
      if(checkMultiUse(Final,j)) return 1;   
    }
  }else if(j.relation().nFinalStateParticles()>1){
    for(int l=0;l<j.relation().nFinalStateParticles();l++){
      Particle Final(j.relation().finalStateParticle(l));
      if(checkMultiUse(i,Final)) return 1;  
    }  
  }else{
    if(i.charge() !=0 && j.charge()!=0 &&
       i.relation().mdstCharged() && 
       j.relation().mdstCharged() &&
       i.relation().mdstCharged().get_ID()==
       j.relation().mdstCharged().get_ID()) return 1;
    if(i.charge()==0 && j.charge()==0 &&
       i.relation().mdstGamma() &&
       j.relation().mdstGamma() &&
       i.relation().mdstGamma().get_ID()==
       j.relation().mdstGamma().get_ID()) return 1;
  }
  return 0;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
