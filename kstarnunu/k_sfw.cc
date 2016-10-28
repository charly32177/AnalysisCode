// -*- C++ -*-
//
// Package:     <package>
// Module:      k_sfw
// 
// Description: <one line class summary>
//
// Implimentation:
//     <Notes on implimentation>
//
// Author:      Hidekazu Kakuno
// Created:     ¿å  7·î  9 01:10:13 JST 2003
// $Id: k_sfw.cc,v 1.3 2003/07/25 14:53:14 kakuno Exp $
//
// Revision history
//
// $Log: k_sfw.cc,v $
// Revision 1.3  2003/07/25 14:53:14  kakuno
// Minor bug fix
//
// Revision 1.2  2003/07/17 18:11:54  kakuno
// k_sfw version 2.3: Modified for PDF(mm2) correction
//
// Revision 1.1.1.1  2003/07/17 17:54:38  kakuno
// k_sfw version 2.2:  Initial version
//
#include "belle.h"

#include "k_sfw.h"


#define USE_FINALSTATE_FOR_SIG 0
// 1 --> kk pipi kpi
// 0 --> k0k0 k0pi0

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//===============================
// coefficients
//===============================
        double k_sfw::alpha_k_sfw0[17] = {
                -0.006357,-0.0993837,0.0724986,0.024342,
                0.140341,0.00358143,-0.497933,0.0050519,
                -0.107423,0.0216523,0.0678463,0.0407537,
                -8.49392e-07,0.00140314,4.8663e-06,3.52144e-05,
                3.8196e-05};
 
 
        double k_sfw::alpha_k_sfw1[17] = {
                0.457001,-0.178582,0.597329,1.82596,
                0.354565,-3.3024,-3.51732,-1.00848,
                -0.286749,-0.96159,-0.744596,-0.52169,
                -21.7805,0.738892,-13.6898,-3.12317,
                7.22794};
 
 
        double k_sfw::alpha_k_sfw2[17] = {
                0.554658,-4.19462,-3.57367,0.0089497,
                0.357133,-3.82704,-4.49298,-1.62434,
                -0.327955,-0.940992,-0.519291,-0.970544,
                -28.3374,1.10484,-8.0394,-2.41847,
                5.81046};
 
 
        double k_sfw::alpha_k_sfw3[17] = {
                0.578035,-7.14771,-6.59138,-1.388,
                0.358921,-3.82085,-4.75458,-1.86224,
                -0.232916,-0.9284,-0.638403,-0.820278,
                -23.1408,4.21952,-2.85437,-1.71457,
                5.01739};
 
 
        double k_sfw::alpha_k_sfw4[17] = {
                0.655201,-7.4282,-6.98985,-0.973756,
                0.347536,-3.7924,-5.30328,-1.96806,
                -0.140798,-0.846041,-0.835415,-0.761708,
                -26.365,5.69024,3.80912,-4.25854,
                5.19989};
 
 
        double k_sfw::alpha_k_sfw5[17] = {
                0.775827,-6.79078,-6.45771,-0.452798,
                0.351063,-4.04164,-6.34687,-2.09724,
                0.156894,-0.834613,-0.414389,-1.16234,
                -25.0239,7.08732,12.0576,-3.87403,
                4.02663};
 
 
        double k_sfw::alpha_k_sfw6[17] = {
                1.14423,-3.87207,-3.6551,-0.0322557,
                0.485723,-4.19427,-6.96612,-1.87398,
                -0.0400051,-0.680864,-1.07484,-0.689534,
                -22.7616,14.1667,29.3356,-7.04065,
                6.08013};
//original ones
/*double k_sfw::alpha_k_sfw0[17] = {
  0.114787,2.70994,2.1902,3.38913,
  1.3674,-4.46261,-4.33278,-0.24901,
  -0.21691,-0.00708349,0.973717,-0.106363,
  -8.10913,7.91641,0.47269,-6.41955,
  2.05342};


double k_sfw::alpha_k_sfw1[17] = {
  0.0582606,1.48162,1.2074,2.80741,
  1.45865,-4.69118,-5.93734,-1.80851,
  -0.654154,-0.689644,-0.106414,-0.658869,
  -15.0721,4.25763,2.94593,-5.3332,
  5.79317};


double k_sfw::alpha_k_sfw2[17] = {
  0.0991841,0.960743,0.557976,2.77684,
  1.45846,-4.75921,-6.25443,-2.04521,
  -0.573428,-0.642085,-0.128891,-0.72723,
  -24.2345,6.7299,3.53798,-4.32955,
  3.27681};


double k_sfw::alpha_k_sfw3[17] = {
  0.15778,1.83783,1.41151,3.62139,
  1.41658,-4.90182,-6.67407,-2.32502,
  -0.57902,-0.632242,0.0280757,-0.709468,
  -26.7782,7.21913,4.43299,-2.10529,
  2.53483};


double k_sfw::alpha_k_sfw4[17] = {
  0.175397,2.93265,2.20164,4.00415,
  1.26478,-4.74294,-6.678,-2.38965,
  -0.412136,-0.691847,-0.0224391,-0.413954,
  -27.9105,6.41955,6.52874,-6.29493,
  -1.53581};


double k_sfw::alpha_k_sfw5[17] = {
  0.270484,2.87845,2.08277,2.94351,
  1.40778,-5.19078,-7.05213,-2.21526,
  -0.709095,-0.466759,-0.0264784,-0.151382,
  -23.9441,12.6333,6.40984,-5.92907,
  3.33395};


double k_sfw::alpha_k_sfw6[17] = {
  0.346219,-0.981985,-2.112,-2.04251,
  0.883324,-5.05439,-6.82521,-1.48502,
  -0.637536,-0.470887,-0.490566,-0.342834,
  11.1658,24.5873,19.364,-15.2887,
  7.5127};*/

// mm2_correction_sig = PDF_sig_data(mm2)/PDF_sig_mc(mm2)
double k_sfw::mm2_correction_sig[7] = {
  1.12284,0.997957, 0.993946, 1.0938, 1.01227, 0.890187, 0.886456};

// mm2_correction_bkg = PDF_bkg_data(mm2)/PDF_bkg_mc(mm2)
double k_sfw::mm2_correction_bkg[7] = {
  1., 1., 1., 1., 1., 1., 1.};
  
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif


//===============================
#include <fstream>
#include <iostream>
#include "CLHEP/Vector/LorentzVector.h"
#include "particle/Particle.h"
#include "ip/IpProfile.h"
#include "helix/Helix.h"
#include "brutus/brutus_f.h"
#include "panther/panther.h"
#include "fix_mdst/fix_mdst.h"
#include "kfitter/kvertexfitter.h"
#include "mdst/findKs.h"
#include "tuple/BelleTupleManager.h"
#include MDST_H

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

brutus_f*
k_sfw::bf = NULL;

BelleTuple*
k_sfw::nt = NULL;

double
k_sfw::sp[7][9];
double
k_sfw::bp[7][9];


void
k_sfw::initialize(brutus_f* brutus_fisher, const unsigned correct_mm2)
{
  extern BelleTupleManager *BASF_Histogram;
  BelleTupleManager *Tm = BASF_Histogram;
  bf = brutus_fisher;
  bf[0].histogram(17, "k_SuperFox-Wolfram0");
  bf[0].coeff(alpha_k_sfw0);
  bf[1].histogram(17, "k_SuperFox-Wolfram1");
  bf[1].coeff(alpha_k_sfw1);
  bf[2].histogram(17, "k_SuperFox-Wolfram2");
  bf[2].coeff(alpha_k_sfw2);
  bf[3].histogram(17, "k_SuperFox-Wolfram3");
  bf[3].coeff(alpha_k_sfw3);
  bf[4].histogram(17, "k_SuperFox-Wolfram4");
  bf[4].coeff(alpha_k_sfw4);
  bf[5].histogram(17, "k_SuperFox-Wolfram5");
  bf[5].coeff(alpha_k_sfw5);
  bf[6].histogram(17, "k_SuperFox-Wolfram6");
  bf[6].coeff(alpha_k_sfw6);
  nt = Tm->ntuple("k_SuperFox-Wolfram", "ksfw iksfw");
  std::ifstream fin("k_sfw.dat");
  for (int i = 0; i < 7; i++){
    for (int j = 0; j < 8; j++){
      fin >> sp[i][j];
    }
  }
  for (int i = 0; i < 7; i++){
    for (int j = 0; j < 8; j++){
      fin >> bp[i][j];
    }
  }
  for (int i = 0; i < 7; i++){
    sp[i][8] = (correct_mm2&1U) ? mm2_correction_sig[i] : 1.0;
    bp[i][8] = (correct_mm2&2U) ? mm2_correction_bkg[i] : 1.0;
  }
  if (correct_mm2&1U){
    std::cout << "k_sfw:: correct PDF(mm2) for signal" << std::endl;
  }
  if (correct_mm2&2U){
    std::cout << "k_sfw:: correct PDF(mm2) for background" << std::endl;
  }
  std::cout << "k_sfw:: initialized" << std::endl;
}


k_sfw::k_sfw(const Particle& b)
{
  // kinematically allowed maximum momentum for mbc>5.2
  // sqrt(5.29^2 - 5.2^2) ~ 1.0GeV
  // 10.58/4 + 1.0/2 = 3.145
  static const double P_MAX(3.2);

  //=====================
  // create particle list
  //=====================
  static const double ELER(3.5);
  static const double EHER(7.9965);
  static const double THETA(0.022);
  static const HepLorentzVector cm(-EHER*sin(THETA), 0., 
				   -EHER*cos(THETA)+ELER, EHER+ELER);
  static const Hep3Vector CMBoost(cm.boostVector());
  static const double PIMASS(0.139570);
  HepLorentzVector pb(b.p());
  pb.boost(CMBoost);
  const double Eb(fix_mdst::get_Benergy());
  const double Hso0_max(2.*(2.*Eb-pb.e())); // 2 * (Eb - delta_E)

  std::vector<Hep3Vector> p_cms_sig, p_cms_other;
  std::vector<int> Q_sig, Q_other;
  HepLorentzVector p_cms_miss(0.,0.,0.,2.*Eb);
  m_et = 0.;

  Mdst_vee2_Manager& veeMgr(Mdst_vee2_Manager::get_manager());
  Mdst_charged_Manager& chMgr(Mdst_charged_Manager::get_manager());
  Mdst_gamma_Manager& gMgr(Mdst_gamma_Manager::get_manager());
  const std::vector<Particle *> & final 
    = b.relation().finalStateParticles();

#if USE_FINALSTATE_FOR_SIG == 0
  // for signal
  for (int i = 0; i < b.nChildren(); i++){
    const Particle & dau(b.child(i));
    HepLorentzVector p_cms(dau.p());
    p_cms.boost(CMBoost);
    p_cms_sig.push_back((Hep3Vector)p_cms);
    Q_sig.push_back((int)dau.charge());
    p_cms_miss -= p_cms;
    m_et += p_cms.perp();
  }
#endif
  // ks
  FindKs fKs;
  std::vector<Mdst_charged *> ks_dau;
  for (std::vector<Mdst_vee2>::iterator i = veeMgr.begin();
       i != veeMgr.end(); i++){
    bool sig(false);
    for (std::vector<Particle *>::const_iterator j = final.begin();
	 j != final.end(); j++){
      if ((**j).mdstCharged() && 
	  ((**j).mdstCharged().get_ID() == (*i).chgd(0).get_ID() || 
	   (**j).mdstCharged().get_ID() == (*i).chgd(1).get_ID() )){
	sig = true;
	break;
      }
    }
    if (sig) continue;
    if ( (*i).kind() != 1 ) continue;
    fKs.candidates( *i, IpProfile::position() );
    if( !fKs.goodKs() ) continue;
    HepLorentzVector p_cms((*i).px(), (*i).py(), (*i).pz(), (*i).energy());
    p_cms.boost(CMBoost);
    if (p_cms.rho() > P_MAX) continue;
    p_cms_other.push_back((Hep3Vector)p_cms);
    Q_other.push_back(0);
    ks_dau.push_back(&(*i).chgd(0));
    ks_dau.push_back(&(*i).chgd(1));
    p_cms_miss -= p_cms;
    m_et += p_cms.perp();
  }

  // charged tracks
  HepVector a(5);
  HepSymMatrix Ea(5,0);
  for (std::vector<Mdst_charged>::const_iterator i = chMgr.begin();
       i != chMgr.end(); i++){
    bool ksdau(false);
    for (std::vector<Mdst_charged *>::const_iterator j = ks_dau.begin();
	 j != ks_dau.end(); j++){
      if ((**j).get_ID() == (*i).get_ID()){
	ksdau = true;
	break;
      }
    }
    if (ksdau) continue;
    bool sig(false);
    HepLorentzVector p_cms;
    for (std::vector<Particle *>::const_iterator j = final.begin();
	 j != final.end(); j++){
      if ((**j).mdstCharged() && (**j).mdstCharged().get_ID() == (*i).get_ID()){
	p_cms = (**j).p();
	sig = true;
      }
    }
    if (sig){
#if USE_FINALSTATE_FOR_SIG == 1
      p_cms.boost(CMBoost);
      p_cms_sig.push_back((Hep3Vector)p_cms);
      Q_sig.push_back((int)(*i).charge());
      p_cms_miss -= p_cms;
      m_et += p_cms.perp();
#endif
    }else{
      HepPoint3D p((*i).trk().mhyp(2).pivot_x(),
		   (*i).trk().mhyp(2).pivot_y(),
		   (*i).trk().mhyp(2).pivot_z());
      a[0] = (*i).trk().mhyp(2).helix(0);
      a[1] = (*i).trk().mhyp(2).helix(1);
      a[2] = (*i).trk().mhyp(2).helix(2);
      a[3] = (*i).trk().mhyp(2).helix(3);
      a[4] = (*i).trk().mhyp(2).helix(4);
      Ea[0][0] = (*i).trk().mhyp(2).error(0);
      Ea[1][0] = (*i).trk().mhyp(2).error(1);
      Ea[1][1] = (*i).trk().mhyp(2).error(2);
      Ea[2][0] = (*i).trk().mhyp(2).error(3);
      Ea[2][1] = (*i).trk().mhyp(2).error(4);
      Ea[2][2] = (*i).trk().mhyp(2).error(5);
      Ea[3][0] = (*i).trk().mhyp(2).error(6);
      Ea[3][1] = (*i).trk().mhyp(2).error(7);
      Ea[3][2] = (*i).trk().mhyp(2).error(8);
      Ea[3][3] = (*i).trk().mhyp(2).error(9);
      Ea[4][0] = (*i).trk().mhyp(2).error(10);
      Ea[4][1] = (*i).trk().mhyp(2).error(11);
      Ea[4][2] = (*i).trk().mhyp(2).error(12);
      Ea[4][3] = (*i).trk().mhyp(2).error(13);
      Ea[4][4] = (*i).trk().mhyp(2).error(14);
      Helix hel(p, a, Ea);
      hel.pivot(IpProfile::position());
      if (std::fabs(hel.dr()) > 5.0 || std::fabs(hel.dz()) > 10.0) continue;

      HepSymMatrix Err(7,0);
      HepPoint3D Pos;
      HepLorentzVector Mom(hel.momentum(0., PIMASS, Pos, Err));

      kvertexfitter kvf;
      kvf.initialVertex(IpProfile::position());
      kvf.beamProfile(IpProfile::position_err_b_life_smeared());
      kvf.addTrack(Mom, Pos, Err, (*i).charge(), PIMASS);
      
      p_cms = kvf.fit() ? Mom : kvf.momentum(0);
      //p_cms = Mom;
      p_cms.boost(CMBoost);
      if (p_cms.rho() > P_MAX) continue;
      p_cms_other.push_back((Hep3Vector)p_cms);
      Q_other.push_back((int)(*i).charge());
      p_cms_miss -= p_cms;
      m_et += p_cms.perp();
    }
  }

  //gammas
  for (std::vector<Mdst_gamma>::const_iterator i = gMgr.begin();
       i != gMgr.end(); i++){
    bool sig(false);
    HepLorentzVector p_cms;
    for (std::vector<Particle *>::const_iterator j = final.begin();
	 j != final.end(); j++){
      if ((**j).mdstGamma() && (**j).mdstGamma().get_ID() == (*i).get_ID()){
	p_cms = (**j).p();
	sig = true;
      }
    }
    if (sig){
#if USE_FINALSTATE_FOR_SIG == 1
      p_cms.boost(CMBoost);
      p_cms_sig.push_back((Hep3Vector)p_cms);
      Q_sig.push_back(0);
      p_cms_miss -= p_cms;
      m_et += p_cms.perp();
#endif
    }else{
      p_cms.setPx((*i).px());
      p_cms.setPy((*i).py());
      p_cms.setPz((*i).pz());
      p_cms.setE(p_cms.rho());
      if (p_cms.rho() < 0.05) continue;
      //if ((*i).ecl().match() || p_cms.rho() < 0.1) continue;
      p_cms.boost(CMBoost);
      if (p_cms.rho() > P_MAX) continue;
      p_cms_other.push_back((Hep3Vector)p_cms);
      Q_other.push_back(0);
      p_cms_miss -= p_cms;
      m_et += p_cms.perp();
    }
  }
  m_mm2 = p_cms_miss.e() > 0 
    ? p_cms_miss.mag2() 
    : -p_cms_miss.e()*p_cms_miss.e() - p_cms_miss.vect().mag2();

  //=======================
  // calculate discriminants
  //=======================
  std::vector<Hep3Vector>::iterator pi, pj;
  std::vector<int>::iterator Qi, Qj;

  // calculate Hso components
  for (int i = 0; i < 3; i++){
    for (int k = 0; k < 5; k++){
      m_Hso[i][k] = 0.;
    }
  }

  for (pi = p_cms_sig.begin(), Qi = Q_sig.begin();
       pi != p_cms_sig.end(); pi++, Qi++){
    const double pi_mag((*pi).mag());
    for (pj = p_cms_other.begin(), Qj = Q_other.begin();
	 pj != p_cms_other.end(); pj++, Qj++){
      const double pj_mag((*pj).mag());
      const double ij_cos((*pi)*(*pj)/pi_mag/pj_mag);
      const int c_or_n(0 == (*Qj) ? 1 : 0);  // 0: charged 1: neutral
      for (int k = 0; k < 5; k++){
	m_Hso[c_or_n][k] += ( k % 2 ) 
	  ? (*Qi)*(*Qj)*pj_mag*legendre(ij_cos, k)
	  : pj_mag*legendre(ij_cos, k);
      }
    }
    const double p_miss_mag(p_cms_miss.rho());
    const double i_miss_cos((*pi)*p_cms_miss.vect()/pi_mag/p_miss_mag);
    for (int k = 0; k < 5; k++){
      m_Hso[2][k] += ( k % 2 ) ? 0. : p_miss_mag*legendre(i_miss_cos, k);
    }
  }

  // add missing to the lists
  p_cms_other.push_back((Hep3Vector)p_cms_miss);
  Q_other.push_back(0);

  // calculate Hoo components
  for (int k = 0; k < 5; k++){
    m_Hoo[k] = 0.;
  }
  for (pi = p_cms_other.begin(), Qi = Q_other.begin();
       pi != p_cms_other.end(); pi++, Qi++){
    const double pi_mag((*pi).mag());
    for (pj = p_cms_other.begin(), Qj = Q_other.begin();
	 pj != pi; pj++, Qj++){
      //for (pj = p_cms_other.begin(), Qj = Q_other.begin();
      //pj != p_cms_other.end(); pj++, Qj++){
      const double pj_mag((*pj).mag());
      const double ij_cos((*pi)*(*pj)/pi_mag/pj_mag);
      for (int k = 0; k < 5; k++){
	m_Hoo[k] += ( k % 2 ) 
	  ? (*Qi)*(*Qj)*pi_mag*pj_mag*legendre(ij_cos, k)
	  : pi_mag*pj_mag*legendre(ij_cos, k);
      }
    }
  }

  // nomalize so that it does not dependent on delta_e
  for (int k = 0; k < 5; k++){
    for (int j = 0; j < ((k%2) ? 1 : 3); j++){
      m_Hso[j][k] /= Hso0_max;
    }
    m_Hoo[k] /= (Hso0_max*Hso0_max);
  }

  //=======================
  // fisher discriminant
  //=======================
  if (bf){
    int i = 0;
    bf[i_mm2()].put(i++,m_et);
    for (int k = 0; k < 5; k++){
      for (int j = 0; j < ((k%2) ? 1 : 3); j++){
	bf[i_mm2()].put(i++, m_Hso[j][k]);
      }
    }
    for (int k = 0; k < 5; k++){
      bf[i_mm2()].put(i++, m_Hoo[k]);
    }
    bf[i_mm2()].fill();
    m_fd = bf[i_mm2()].discriminant();
  }else{
    m_fd = 0.;
  }
  nt->column("ksfw", m_fd);
  nt->column("iksfw", i_mm2());
  nt->dumpData();
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

