//New version for B+->knunubar at 2015/09/07
// Hadronic Tag (full reconstruction)
// Charged and mixed B meson all included
// B+->h(*)nunubar
// By Bo-Yuan Yang
// Start date: 2015/05/14
// Selection: Rank1 B only for each event
// An utility: Can know wherher a B0bar meson 
////////////////////////////////////////////
// 2015 05 18 add : Thrust cos
// 2015 06 02 add : NB variable -> costheta between miss momentum and all event momentum direction, E_ecl in miss momemtum direction, remaining track before cut
// 2015 09 07 modify : fix the bug for Not-working issue of dup-track, cut on dr dz and cut on p_t.

#include "belle.h"
#include "shape.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include "mdst/findKs.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "kfitter/kvertexfitter.h"
#include "particle/Particle.h"
#include "particle/PID.h"
#include "particle/gammac.h"
#include "particle/utility.h"
#include "particle/combination.h"
#include "hamlet/Hamlet.h" //for qr
#include "kid/atc_pid.h"
#include "eid/eid.h"
#include "mdst/Muid_mdst.h"
#include "mdst/mdst.h"
#include "findLambda.h"
#include "benergy/BeamEnergy.h"
#include "ip/IpProfile.h"
#include "brutus/brutus_f.h"

#include "toolbox/FoxWolfr.h"
#include "toolbox/Thrust.h"
#include "tagv/TagV.h"
#include "ExKFitter.h"
#include "kid_eff_06.h"
#include "kfitter/kmassfitter.h"
#include "kfitter/kmassvertexfitter.h"
#include "kfitter/khelix2xyz.h"
#include "kfitter/kfitterparticle.h"
#include "kfitter/kfitter_ini.h"
#include "kfitter/kfitterbase.h"
#include "kfitter/kfittererror.h"
#include "kfitter/kmakemother.h"
#include "kchisq2confi.cc"


#include "helix/Helix.h"
#include "k_sfw.h"
#include "CLHEP/Vector/LorentzVector.h"  // /afs/afs11/belle/belle/b20090127_0910/src/util/belleCLHEP/Vector
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "userinfo.h"
#include "panther/panther.h"

//new including headfiles involving fullreconstruction
#include "hamlet/AnaBrecon.h"
#include "fullrecon/frec_util.h"
#include "tables/ekpfullrecon_panther.h"
#include BRECON_H
#include FULLRECON_H
#include HEPEVT_H
#include BELLETDF_H
#include MDST_H
#define THRESHOLD 0.25

////////////////////////////////////////////////
// mass list
const double protonmass   = 0.938272;
const double kaonmass     = 0.493677;
const double pionmass     = 0.139570;
const double muonmass     = 0.105658;
const double electronmass = 0.000510998;
//Others mass
const double Lambdamass = 1.115683;
////////////////////////////////////////////////

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
//utilities

int checkMultiUse(Particle& i, Particle& j);


class ana_fullreconk : public Module 
{
public:
     ana_fullreconk();
     ~ana_fullreconk() {}

     void init( int * );
     void term();
     void disp_stat( const char* ) {}
     void hist_def();
     void event( BelleEvent*, int* );
     void begin_run( BelleEvent*, int* );
     void end_run( BelleEvent*, int * );
     void other( int*, BelleEvent*, int* ) {}
     void shape(Particle &, float &, float &, float &, float &, float par_sfw[17],Vector4 &,HepPoint3D &,int &, HepPoint3D &);
     void GetImpactParameters( const Mdst_charged*, double* ,double* ,int );
     void Getb_tpls();

    double deltaZ(Particle, double&, double&, double&, double&);
    double GetCos(Particle a,Particle b);
    double pt(Particle);
    double pmag(Particle);
    void remove_dup_trk(std::vector<Particle>& , std::vector<Particle>& \
               ,std::vector<Particle>& , std::vector<Particle>& \
               ,std::vector<Particle>& , std::vector<Particle>& \
               ,std::vector<Particle>& , std::vector<Particle>& \
              ,std::vector<Particle>& , std::vector<Particle>& ,double , double );
    void remove_dup_trk(std::vector<Particle>& , double , double );

private:
     BelleTuple* TagB_tpl;
     BelleTuple* B_LX_tpl;
     BelleTuple* b_tpl;
	  int evtcount;
	 brutus_f Fisher_ksfw[7];
}; //end class definition

extern "C" Module_descr *mdcl_ana_fullreconk() 
{
    ana_fullreconk *module = new ana_fullreconk;
    Module_descr *dscr = new Module_descr( "ana_fullreconk", module );
    return dscr;
}//Module_descr *mdcl_ana_fullreconk()
ana_fullreconk::ana_fullreconk() {  return; }//ana_fullreconk::ana_fullreconk()
void ana_fullreconk::init(int *) 
{
	evtcount=0;
	Hamlet::init(); 
    return;
}

void ana_fullreconk::term (void)
{
   return;
}
// Set IP location
HepPoint3D IP( 0, 0, 0 );
HepSymMatrix IPerr( 3, 0 );
void ana_fullreconk::begin_run( BelleEvent* evptr, int *status ) 
{
  eid::init_data();             // available in the new lib (after b199907*)
  BeamEnergy::begin_run();      //beam energy
  IpProfile::begin_run();       // Get IP profile data from $BELLE_POSTGRES_SERVER
  IpProfile::dump();            // Dump IP profile data to STDOUT (optional)
  // Set IP and error
  IP    = IpProfile::position();
  IPerr = IpProfile::position_err();
  //tagging
  // To initialize LH tables by EvtGen MC
  Hamlet::begin_run(Hamlet::MULT_DIM_LH);
  return;
}//ana_fullreconk::begin_run
void ana_fullreconk::end_run(BelleEvent* evptr, int *status) { (void)evptr; (void)status; return; }

void ana_fullreconk::hist_def() 
{
    extern BelleTupleManager *BASF_Histogram;
    BelleTupleManager *tm = BASF_Histogram;

    TagB_tpl = BASF_Histogram->ntuple( "Rank1 tagB meson","ebeam Tag_RunID Tag_EvtID Tag_ExpID tagid best nbout csbest csnbout tagBmbc tagBde mcinfo nFS decmod ddecmod1 ddecmod2 ddecmod3 ddecmod4 tagBpc PX PY PZ PhEnergy BtagE Ek pxk pyk pzk massk sumecl hi1 mo1 gmo1 ggmo1 hindex mmissb pxmissb pymissb pzmissb emissb pmissb mmissCM pxmissCM pymissCM pzmissCM emissCM pmissCM kpCM kpxCM kpyCM kpzCM keCM mkcm keb kpb kpxb kpyb kpzb mmiss2b mmiss2CM cosmisslab cossklab bpd1 bpd2 bpd3 bpd4 bpd5 bpd6 bpd7 bpd8 bpd9 nbpd bnd1 bnd2 bnd3 bnd4 bnd5 bnd6 bnd7 bnd8 bnd9 nbnd charged chisqexk mmisslab mmiss2lab pxmisslab pymisslab pzmisslab emisslab pmisslab mklab kpxlab kpylab kpzlab kelab sum_eclinmissdir cosbmissandallcm anglebmissandallcm cosklab nb_trackremaining tnb_trackremaining tnb_pi0remaining DistTOOtherBdz imm2 mmm2 H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son H4son H0soc H1soc H2soc H3soc H4soc H0som H1som H2som H3som H4som et kdr kdz RatioToOtherType RaatioToSecondBest vchisq_sig DOF_sig pvalue_sig vchisq_tag DOF_tag pvalue_tag Distsign mB2 sbcm reenergy sblab qsquare cosbcm", 1 );


 //b_tpl = BASF_Histogram->ntuple("B->pi0", "pi0mass gamcos gam1e gam2e hi0gam mo0gam gmo0gam ggmo0gam hindexgam hi0gam mo0gam gmo0gam ggmo0gam coneangle_ecl",2);

 k_sfw::initialize(Fisher_ksfw);

    return;
}

/////////////////////////////////////////////////////////
void ana_fullreconk::event(BelleEvent* evptr, int* status)
{
    const HepPoint3D             &ip     = IpProfile::position();

  Belle_event_Manager& bevt_mgr 	= Belle_event_Manager::get_manager();
  Mdst_charged_Manager &charged_mag    	= Mdst_charged_Manager::get_manager();
  Mdst_vee2_Manager& Vee2Mgr            = Mdst_vee2_Manager::get_manager();
  Gen_hepevt_Manager&  gen_mgr          = Gen_hepevt_Manager::get_manager();
  Mdst_gamma_Manager& gamma_Mgr         = Mdst_gamma_Manager::get_manager();
  Mdst_ecl_Manager &ecl_mag             = Mdst_ecl_Manager::get_manager();
  Mdst_ecl_aux_Manager& mdst_Aux_Mgr    = Mdst_ecl_aux_Manager::get_manager();
  Mdst_ecl_trk_Manager &ecltrk_mag      = Mdst_ecl_trk_Manager::get_manager();

  Mdst_pi0_Manager &pi0_mag             = Mdst_pi0_Manager::get_manager();
  Mdst_event_add_Manager& mevtmgr       = Mdst_event_add_Manager::get_manager();
  Mdst_vee_daughters_Manager& veedmgr   = Mdst_vee_daughters_Manager::get_manager();
  Mdst_klm_mu_ex_Manager& klmmgr        = Mdst_klm_mu_ex_Manager::get_manager();

  double E_HER=BeamEnergy::E_HER();
  double E_LER=BeamEnergy::E_LER();
  double cross_angle=BeamEnergy::Cross_angle();//radian
  double ebeam = BeamEnergy::E_beam_corr();

   Mdst_ecl_aux_Manager &eclaux_mag     = Mdst_ecl_aux_Manager::get_manager();

  //Maximum allowed momentum
 // const double pmax = sqrt( ebeam*ebeam - Lambdamass*Lambdamass ); //used for scaled momentum
  static Vector4 P_BEAM ( -E_HER*sin(cross_angle), 0.0,  E_LER-E_HER*cos(cross_angle), E_HER+E_LER );
	HepLorentzVector boost_vector( -E_HER*sin(cross_angle), 0.0, E_LER-E_HER*cos(cross_angle), E_HER+E_LER );

  //series Run-Evt-Exp and event is One to One
  int Run=0, Evt=0, Exp=0, farm=0; 
  Run = bevt_mgr[0].RunNo();
  Evt = bevt_mgr[0].EvtNo() & 0x0FFFFFFF;
  Exp = bevt_mgr[0].ExpNo();
  farm= bevt_mgr[0].EvtNo() >> 28;

  // set IP and error
  int IPUsable = 0;
  if( IpProfile::usable() )
  {
	IP    = IpProfile::position(1);
        IPerr = IpProfile::position_err_b_life_smeared(1);
        IPUsable = 1;
  }
  else
  {
        IP    = HepPoint3D(0,0,0);
        IPerr = HepSymMatrix(3,0);
  }

        std::vector<Particle> AllTrkList;
        std::vector<Particle> gamma;
        std::vector<Particle> pi_0;
        std::vector<Particle> pi_plus;
        std::vector<Particle> pi_minus;
        std::vector<Particle> k_plus;
        std::vector<Particle> k_minus;
        std::vector<Particle> k_short;
        std::vector<Particle> p_plus;
        std::vector<Particle> p_minus;
        std::vector<Particle> mu_plus;
        std::vector<Particle> mu_minus;
        std::vector<Particle> e_plus;
        std::vector<Particle> e_minus;
        std::vector<Particle> D0_cand;
        std::vector<Particle> D0B_cand;
        std::vector<Particle> B_cand1;
        std::vector<Particle> B_cand2;
        std::vector<Particle> pi_plus_cand;
	std::vector<Particle> gamma1;
        std::vector<Particle> Lamda;
        std::vector<Particle> Lamdabar;
        std::vector<Particle> Sigma_plus;
        std::vector<Particle> D_minus_cand;

	std::vector<Particle> te_plus;
	std::vector<Particle> te_minus;
	std::vector<Particle> tmu_plus;
	std::vector<Particle> tmu_minus;
	std::vector<Particle> tpi_plus;
	std::vector<Particle> tpi_minus;
	std::vector<Particle> tk_plus;
	std::vector<Particle> tk_minus;
	std::vector<Particle> tp_plus;
	std::vector<Particle> tp_minus;

	std::vector<Particle> tpi_0;
	std::vector<Particle> tk_short;
	std::vector<Particle> tLamda;
	std::vector<Particle> tLamdabar;
	std::vector<Particle> tgamma;
	std::vector<Particle> tgamma1;

        Ptype ptype_gamma("GAMM");
        Ptype ptype_pi_0("PI0");
        Ptype ptype_pi_plus("PI+");
        Ptype ptype_pi_minus("PI-");
        Ptype ptype_k_plus("K+");
        Ptype ptype_k_minus("K-");
        Ptype ptype_k_short("K0S");
        Ptype ptype_p_plus("P+");
        Ptype ptype_p_minus("AP+");
        Ptype ptype_mu_plus("MU+");
        Ptype ptype_mu_minus("MU-");
        Ptype ptype_e_plus("E+");
        Ptype ptype_e_minus("E-");
        Ptype ptype_D0_cand("D0");
        Ptype ptype_D0B_cand("D0B");
        Ptype ptype_B_plus_cand("B+");
        Ptype ptype_B0_cand("B0");
        Ptype ptype_B_minus_cand("B-");
        Ptype ptype_B0B_cand("B0B");
        Ptype ptype_D_minus_cand(-411);

        atc_pid selkpi(3,1,5,3,2);
        atc_pid selppi(3,1,5,4,2);
        atc_pid selkp(3,1,5,3,4);
        atc_pid selpk(3,1,5,4,3);
        atc_pid selmuk(3,1,5,1,3);
        atc_pid selmupi(3,1,5,1,2);

 int MCstatus=0;
  for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin();
       i != gen_mgr.end(); i++)
        {
                MCstatus=1;
        }

////////////////////////////// full recontruction (hadronic B-tagging method) ////////////////////////////////////////////////

  // AnaBrecon class is used to create a Particle object from the entry in Brecon table
  AnaBrecon brecon;   //"hamlet/AnaBrecon.h"
  // Check if Brecon table exists and has at least one entry
  Brecon_header_Manager & bhMgr  = Brecon_header_Manager::get_manager();
  if ( bhMgr.count() <= 0 ) { *status=0; return; }


////////////////////////////////////////////////////////////////////////////////////////////
//                         B+ B- B0 B0bar Tagging (full reconstruction)                   //
////////////////////////////////////////////////////////////////////////////////////////////
        double this_type_rank = -1;
        double second_best_nbout = -1;
        double other_type_nbout = -1;
        double this_type_nbout = -1;
        double RatioToOtherType = 0;
        double RaatioToSecondBest = 0;
 
  // First Loop over candidates in the Ekpfullrecon panther table to find out B+ and B0
  Ekpfullrecon_Manager &frec_mgr = Ekpfullrecon_Manager::get_manager();
 for( Ekpfullrecon_Manager::iterator frecon_it = frec_mgr.begin(); frecon_it != frec_mgr.end(); frecon_it++ )
  {
         Ekpfullrecon &frecB = *frecon_it;
         // create a B particle
         Particle & Btag = const_cast<Particle &>( brecon.getParticle( (int)frecB.get_ID() ) );
        double bestCont  = frecB.cont_NBRank();
        double Contnbout = frecB.cont_NBout();
        if(abs(frecB.tag_id())==511 && bestCont==1)
        {
            this_type_rank =bestCont ;
            other_type_nbout = Contnbout ;
	  //  break;
        }
        if (bestCont!=1 && bestCont > this_type_rank && abs(frecB.tag_id())==521)     
	{
	    this_type_rank = bestCont ;
	    second_best_nbout = Contnbout ;
	   // break;
	}
  }
  for( Ekpfullrecon_Manager::iterator frecon_it = frec_mgr.begin(); frecon_it != frec_mgr.end(); frecon_it++ )
  {
        Ekpfullrecon &frecB = *frecon_it;
         // create a B particle
        double bestCont  = frecB.cont_NBRank();
        double Contnbout = frecB.cont_NBout();

        Particle & Btag = const_cast<Particle &>( brecon.getParticle( (int)frecB.get_ID() ) );	
	if(abs(frecB.tag_id())==521)
        {
            this_type_nbout = Contnbout ;
        }
        RatioToOtherType = other_type_nbout/ this_type_nbout;
        RaatioToSecondBest = second_best_nbout/ this_type_nbout;
        if(other_type_nbout==-1 || this_type_nbout==-1)
	  RatioToOtherType = -999;
        if(second_best_nbout==-1 || this_type_nbout==-1)
	  RaatioToSecondBest = -999;


	gamma.clear();
        gamma1.clear();
        pi_0.clear();
        pi_plus.clear();
        pi_minus.clear();
        k_plus.clear();
        k_minus.clear();
        k_short.clear();
        p_plus.clear();
        p_minus.clear();
        Lamda.clear();
        Lamdabar.clear();
        Sigma_plus.clear();
        mu_plus.clear();
        mu_minus.clear();
        e_plus.clear();
        e_minus.clear();
	D0_cand.clear();
	D_minus_cand.clear();

        te_plus.clear();
        te_minus.clear();
        tmu_plus.clear();
        tmu_minus.clear();
        tpi_plus.clear();
        tpi_minus.clear();
        tk_plus.clear();
        tk_minus.clear();
        tp_plus.clear();
        tp_minus.clear();
        tpi_0.clear();
	B_cand1.clear();

        tk_short.clear();
        tLamda.clear();
        tLamdabar.clear();
        tgamma.clear();
        tgamma1.clear();

    // fill Ks list from Mdst_vee2 bank
    for(std::vector<Mdst_vee2>::iterator i = Vee2Mgr.begin(); i != Vee2Mgr.end(); i++)
    {


        int flag = 0;
        //int nfinal = Btag.relation().nFinalStateParticles(); //the number of final particle used to tagB
        int nchild = Btag.nChildren();
        for( int k = 0; k < nchild; k++ )
        {
            Particle child = Btag.child(k);
            if ( child.mdstVee2().get_ID()==(*i).get_ID() )
            {
 
                flag = 1;
            }
            if (flag==1) continue;

            int nchildchild = Btag.child(k).nChildren();
            for(int l = 0; l <nchildchild; l++)
                {
                  Particle childchild = child.child(l);
                  if (childchild.mdstVee2().get_ID()==(*i).get_ID())
                  {
                    flag = 1;
                  }
                }
        }

        if (flag == 0)
        {
            if( (*i).kind() == 1 )    // 1 for Ks, 2 for Lambda, 3 for anti-Lambda , 4 for gamma -> ee
            {
                int goodVee = 0;
                FindKs findks;
                findks.candidates( (*i),ip);
                goodVee = findks.goodKs();
                if( goodVee==1)
                {
                Particle tmp(*i);
                k_short.push_back(tmp);
                }
            }

            if( (*i).kind() == 2)
            {
                int goodVee = 0;
                FindLambda findlambda;
                findlambda.candidates( (*i),ip);
                goodVee = findlambda.goodLambda();
                if(goodVee==1)
                {
                    Particle tmp(*i);
                    Lamda.push_back(tmp);
                }
            }
            if( (*i).kind() == 3)
            {
                int goodVee = 0;
                FindLambda findlambda;
                findlambda.candidates( (*i),ip);
                goodVee = findlambda.goodLambda();
                if(goodVee==1)
                {
                  Particle tmp(*i);  
                  Lamdabar.push_back(tmp);
                }
            }

            if( (*i).kind() == 4)
            {
                Particle tmp(*i);
                gamma.push_back(tmp);

            }
        }
    }
        for (std::vector<Particle>::iterator i = k_short.begin(); i != k_short.end(); i++)
        {
            if((*i).mass()>0.46&&(*i).mass()<0.53&&(*i).mdstVee2().chisq()<100)
             tk_short.push_back(*i);
            
        }
      	 for (std::vector<Particle>::iterator i = Lamda.begin();i !=Lamda.end();i++)
        {
            if((*i).mass()<1.05||(*i).mass()>1.18||(*i).mdstVee2().chisq()>100) continue;
            int flag = 0; 
            for(std::vector<Particle>::iterator j = tk_short.begin(); j!= tk_short.end(); j++)
            {
               if((*i).child(0).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                 ||(*i).child(0).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID()\
                 ||(*i).child(1).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                 ||(*i).child(1).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID())
              {flag = 1;break;}
            }
            if (flag==0)
            {
              tLamda.push_back(*i);
            }
        }
    
        for (std::vector<Particle>::iterator i = Lamdabar.begin();i !=Lamdabar.end();i++)
        {
            if((*i).mass()<1.05||(*i).mass()>1.18||(*i).mdstVee2().chisq()>100) continue;
            int flag = 0; 
            for(std::vector<Particle>::iterator j = tk_short.begin(); j!= tk_short.end(); j++)
            {
                 if((*i).child(0).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                 ||(*i).child(0).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID()\
                 ||(*i).child(1).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                 ||(*i).child(1).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID())

              {flag = 1;break;}
            }
            if (flag==0)
            {
              tLamdabar.push_back(*i);
            }
        }
        for (std::vector<Particle>::iterator i = gamma.begin(); i!= gamma.end(); i ++)
        {
            int flag =0;
            if ((*i).mass()>0.3) continue;
            for(std::vector<Particle>::iterator j = tk_short.begin(); j!= tk_short.end(); j++)
            {
               if((*i).child(0).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                ||(*i).child(0).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID()\
                ||(*i).child(1).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                ||(*i).child(1).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID())
		{flag=1;}
            }
            if (flag==1) continue;

            
             for(std::vector<Particle>::iterator j = tLamda.begin(); j !=tLamda.end(); j++)
             {
                 if((*i).child(0).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                 ||(*i).child(0).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID()\
                 ||(*i).child(1).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                 ||(*i).child(1).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID())

               {flag = 1;break;}
             }
            
            if (flag==1) continue;

            
             for(std::vector<Particle>::iterator j = tLamda.begin(); j !=tLamda.end(); j++)
             {
                 if((*i).child(0).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                 ||(*i).child(0).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID()\
                 ||(*i).child(1).mdstCharged().get_ID()==(*j).child(0).mdstCharged().get_ID()\
                 ||(*i).child(1).mdstCharged().get_ID()==(*j).child(1).mdstCharged().get_ID())
               {flag = 1;break;}
             }
            
            if(flag==1) continue;

            
              tgamma.push_back(*i);
            
        }
        for(std::vector<Mdst_gamma>::iterator it = gamma_Mgr.begin();it !=gamma_Mgr.end(); it++)
        {
            int flag=0;
            int nfinal = Btag.relation().nFinalStateParticles();
            for(int k=0; k < nfinal; k++)
            {
                Particle child = Btag.relation().finalStateParticle(k);
                if(child.mdstGamma())
                {
                    if(child.mdstGamma().get_ID()==(*it).get_ID())
                        flag=1;
                }
            }
            if (flag) continue;

            Particle tmp(*it);
            gamma1.push_back(tmp);

        }
	int nremainingtrack = -1;
        for(std::vector<Mdst_charged>::iterator i = charged_mag.begin(); i != charged_mag.end(); i++)
        {
	    nremainingtrack++;
            Mdst_charged& ch = *i;
            int flag = 0;
            int nfinal = Btag.relation().nFinalStateParticles(); // the number of final particle used to tagB
            for(int k = 0 ; k < nfinal; k++)
            {
            Particle child = Btag.relation().finalStateParticle(k);
            if ( child.mdstCharged() )
            {  
              if ( child.mdstCharged().get_ID() == ch.get_ID()) 
              {
                flag=1;break;
              }
            }
            }
	
            if(flag) continue;
/////////////////////////////////////////////////////////////////////
            int flag1 = 0;
            for(std::vector<Particle>::iterator j = tk_short.begin();j != tk_short.end();j++)
            {
                if((*j).child(0).mdstCharged()==(*i)||(*j).child(1).mdstCharged()==(*i))
		 flag1 =1;
            }
            if (flag1 ==1) continue;

	    for(std::vector<Particle>::iterator j = tLamda.begin();j != tLamda.end();j++)
            {
                if((*j).child(0).mdstCharged()==(*i)||(*j).child(1).mdstCharged()==(*i))
		 flag1 =1;
            }
            if (flag1 ==1) continue;
            for(std::vector<Particle>::iterator j = tLamdabar.begin();j != tLamdabar.end();j++)
            {
                if((*j).child(0).mdstCharged()==(*i)||(*j).child(1).mdstCharged()==(*i))
		 flag1 =1;
            }
            if (flag1 ==1) continue;
            for(std::vector<Particle>::iterator j = tgamma.begin();j != tgamma.end();j++)
            {
                if((*j).child(0).mdstCharged()==(*i)||(*j).child(1).mdstCharged()==(*i))
                flag1 =1;
            }
            if (flag1 ==1) continue;

////////////////////////////////////////////////////////////////////////
//                    Satrt Particle identification                   //
///////////////////////////////////////////////////////////////////////
                //eid
                eid sel_e(*i);
                float eid_prob = sel_e.prob(3,-1,5);
                //muid
                Muid_mdst muon( *i );
                int mu_sta = muon.Status();
                int outcome = muon.Outcome();
                int mu_level = muon.Level();
                int reject = muon.Prerejection();
                double mu_like = muon.Muon_likelihood();
                double chi2 = muon.Chi_2();

           // muon
            if( !reject && mu_like>0.95 )
            {
                if ((*i).charge()>0)
                {
                        Particle tmp(*i,ptype_mu_plus);
                        mu_plus.push_back(tmp);
                }
                else 
                {
                        Particle tmp(*i,ptype_mu_minus);
                        mu_minus.push_back(tmp);

                }
            }//end of muon



            // electron
            else if ( eid_prob>0.95 )
            {
                if ((*i).charge()>0)
                {
                  Particle tmp(*i, ptype_e_plus);
                  e_plus.push_back(tmp);
                } 
                if((*i).charge()<0) 
                {
                  Particle tmp(*i, ptype_e_minus);
                  e_minus.push_back(tmp);
                }
            }
            //end of electron


            // e, mu rejection
            else 
            {
                if(selkpi.prob(*i)>0.6)
                    //if(selkpi.prob(*i)>0.)
                {
                    if ((*i).charge()>0)
                    {
                        if (ch.trk().quality() == 0)  // Obtain "Good track"
                        {
                                Particle tmp(*i, ptype_k_plus);
                                k_plus.push_back(tmp);
                        }
                    }
                    else
                    {
                        if (ch.trk().quality() == 0)  // Obtain "Good track"
                        {
                                Particle tmp(*i, ptype_k_minus);
                                k_minus.push_back(tmp);
                        }
                    }
                }
                else if (selkpi.prob(*i)<0.4)   //if (selkpi.prob(*i)<1.0)
                {
                    if ((*i).charge()>0)
                    {
                        if (ch.trk().quality() == 0)  // Obtain "Good track"
                        {
                                Particle tmp(*i, ptype_pi_plus);
                                pi_plus.push_back(tmp);
                        }
                    }
                    else
                    {
                        if (ch.trk().quality() == 0)  // Obtain "Good track"
                        {
                                Particle tmp(*i, ptype_pi_minus);
                                pi_minus.push_back(tmp);
                        }
                    }
                }//end of selkpi
                else
                {
                    if ((*i).charge()>0)
                    {
                        if (ch.trk().quality() == 0)  // Obtain "Good track"
                        {
                                Particle tmp(*i, ptype_p_plus);
                                p_plus.push_back(tmp);
                        }

                        else
                        {
                            if (ch.trk().quality() == 0)  // Obtain "Good track"
                            {
                                    Particle tmp(*i, ptype_p_minus);
                                    p_minus.push_back(tmp);
                            }
                        }
                    }
                }
            }//end loop of e, mu rejection

        }//for(std::vector<Mdst_charged>::iterator i

    std::vector<int> flagep;
    std::vector<int> flagem;
    double alldr=2;
    double alldz=5;
    double allpt=0.1;

    for(std::vector<Particle>::iterator i = e_plus.begin(); i != e_plus.end(); i++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,0);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
        te_plus.push_back(*i);
    }
    for(std::vector<Particle>::iterator i = e_minus.begin();i != e_minus.end(); i++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,0);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
            te_minus.push_back(*i);
    }
    for(std::vector<Particle>::iterator i = mu_plus.begin();i != mu_plus.end(); i++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,1);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
            tmu_plus.push_back(*i);

    }
    for(std::vector<Particle>::iterator i = mu_minus.begin();i != mu_minus.end(); i++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,1);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
            tmu_minus.push_back(*i);
    }

    for(std::vector<Particle>::iterator i = pi_plus.begin(); i !=pi_plus.end(); i++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,2);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
        tpi_plus.push_back(*i);
        HepLorentzVector ch1_p((*i).p());
    }
    for(std::vector<Particle>::iterator i = pi_minus.begin(); i !=pi_minus.end(); i ++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,2);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
        tpi_minus.push_back(*i);
    }
    for(std::vector<Particle>::iterator i = k_plus.begin(); i != k_plus.end(); i++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,3);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
            tk_plus.push_back(*i);
    }
    for(std::vector<Particle>::iterator i = k_minus.begin(); i !=k_minus.end(); i ++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,3);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
            tk_minus.push_back(*i);
    }
	   // std::cout <<"sizepbefore="<< p_plus.size() <<std::endl;

    for(std::vector<Particle>::iterator i = p_plus.begin(); i != p_plus.end(); i++)
    {

        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,4);
           // std::cout <<"drdzptuncut="<<"\t"<< abs(dr) <<"\t"<<abs(dz)<< "\t"<<pt(*i)<<std::endl;

        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
        {
            tp_plus.push_back(*i);
           // std::cout <<"drdzptcut="<<"\t"<< abs(dr) <<"\t"<<abs(dz)<< "\t"<<pt(*i)<<std::endl;

	}
    }
           // std::cout <<"sizepafter="<< tp_plus.size() <<std::endl;

   // std::cout <<"sizebefore="<< tpi_plus.size()+tpi_minus.size()+te_plus.size()+te_minus.size()+tk_plus.size()+tk_minus.size()+tmu_plus.size()+tmu_minus.size()+tp_plus.size()+tp_minus.size() <<std::endl;
    remove_dup_trk(te_plus, te_minus,tmu_plus, tmu_minus,tpi_plus, tpi_minus,tk_plus, tk_minus,tp_plus,tp_minus,alldr,alldz);
  //  std::cout <<"sizeafter="<<  tpi_plus.size()+tpi_minus.size()+te_plus.size()+te_minus.size()+tk_plus.size()+tk_minus.size()+tmu_plus.size()+tmu_minus.size()+tp_plus.size()+tp_minus.size() <<std::endl;

//////////////////////////////////////////////////////////////
    // fill pi0 list from Mdst_pi0 bank
    for(std::vector<Mdst_pi0>::iterator i = pi0_mag.begin(); i != pi0_mag.end(); i++)
    {
        int flag = 0;
        int nchild = Btag.nChildren(); //the number of final particle used to tagB
        for( int k = 0; k < nchild; k++ )
        {
            Particle child = Btag.child(k);
                if ( child.mdstPi0().get_ID() ==(*i).get_ID() )
                {
                    flag = 1;
                }
                if (flag==1) continue;

                int nchildchild = Btag.child(k).nChildren();
              if (nchildchild>0)
              {
                for(int l = 0; l <nchildchild; l++)
                {
                    Particle childchild = child.child(l);
                    if (childchild.mdstPi0().get_ID()==(*i).get_ID())
                    {
                        flag = 1;
                    }

                    int nchildchildchild = Btag.child(k).child(l).nChildren();
                    if (nchildchildchild>0)
                    {
                        for(int m = 0 ; m<nchildchildchild ; m++)
                        {
                            Particle childchildchild = childchild.child(m);
                            if (childchildchild.mdstPi0().get_ID()==(*i).get_ID())
                            {
                                flag = 1;
                            }
                        }
                    }
                }
              }

        }
        if (flag == 0)
        {
            Particle tmp(*i);
            pi_0.push_back(tmp);
        }
    }
////////////////////

    for(std::vector<Particle>::iterator i = pi_0.begin(); i != pi_0.end(); i++)
    {
  // std::cout<<"start pi0"<<std::endl;
      int pi_n;
      double e1,e2,a,p1,p2,p,cos;
      e1 = (*i).child(0).e();
      e2 = (*i).child(1).e();
      a = abs(e1-e2)/(e1+e2);

      p1 = sqrt((*i).child(0).px()*(*i).child(0).px()+(*i).child(0).py()*(*i).child(0).py()+(*i).child(0).pz()*(*i).child(0).pz());
      p2 = sqrt((*i).child(1).px()*(*i).child(1).px()+(*i).child(1).py()*(*i).child(1).py()+(*i).child(1).pz()*(*i).child(1).pz());
      p  = sqrt((*i).child(0).px()*(*i).child(1).px()+(*i).child(0).py()*(*i).child(1).py()+(*i).child(0).pz()*(*i).child(1).pz());
      cos= p/p1*p2;
    
/////////////////////////////////
      if (e1>0.05&&e2>0.05&&(*i).p().mag()>0.1178&&(*i).p().mag()<0.1502&&a<0.9)
	{
          int flags1=-1;
          int flags2=-1; 
          int flags3=-1;


	 for(std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s !=eclaux_mag.end();s++ )
        {      

            Mdst_ecl_aux& ch = *s;
            int flag = 0;
            int nfinal = Btag.relation().nFinalStateParticles(); // the number of final particle used to tagB
            for(int k = 0 ; k < nfinal; k++)
            {

                Particle child = Btag.relation().finalStateParticle(k);
                if ( child.mdstCharged() )
                {
                    if ( child.mdstCharged().trk() == ch.trk())
                    {
                        flag=1;
                    }
                }
                if (child.mdstGamma())
                {
                    if (child.mdstGamma().ecl().get_ID()==ch.get_ID())
                    {
                        flag=1;
                    }
                }

            }
            if(flag) continue;

	  for(std::vector<Mdst_ecl>::iterator p = ecl_mag.begin(); p!=ecl_mag.end(); p++)
            {

              if ((*s).get_ID()==(*p).get_ID())
                  {
                  //forward e>0.1
                  if ((*p).theta() < 0.547 && (*p).energy() > 0.1 )
		      flags1=1;
		//barrel e>0.05
                  if ((*p).theta() > 0.562 && (*p).theta() < 2.246 && (*p).energy()>0.05 )
	              flags2=1;
		//backward e>0.15
                  if ((*p).theta() > 2.246 && (*p).energy() > 0.15)
                      flags3=1;
                 }

            }


	}//end ecl
	

		if(flags1)
		{
			Particle tmp(*i);
			tpi_0.push_back(tmp);
		}	
		if(flags2)
		{
			Particle tmp(*i);
			tpi_0.push_back(tmp);
		}
		if(flags3)
		{
			Particle tmp(*i);
			tpi_0.push_back(tmp);
		}

      }//end Gamma selection
//		   b_tpl->dumpData();
//                   *status = 1;
          // std::cout<<"dump!!"<<std::endl;
    }//end pi_0 

/////////////////////////////////////////////////////////////
//                      pi0 veto                           //
////////////////////////////////////////////////////////////
	int npi0 = -1;
            npi0 = tpi_0.size();
	if (npi0 != 0)
	continue;
//////////////////////////////////////////////////////////////
	int npip = -1;
	int npim = -1;
	int nremaining = -1;
	int npp = -1;
	int npm = -1;
	int nmup = -1;
	int nmum = -1;
	int nkp = -1;
	int nkm = -1;
	int nep = -1;
	int nem = -1;
	int n_tks = -1;
	int n_tlamda = -1;
	int n_tlamdabar = -1;
	int n_tgamma = -1;
	int n_tgamma1 = -1;
	int tnpip = -1;
	int tnpim = -1;
	int tnremaining = -1;
	int tnpp = -1;
	int tnpm = -1;
	int tnmup = -1;
	int tnmum = -1;
	int tnkp = -1;
	int tnkm = -1;
	int tnep = -1;
	int tnem = -1;

/*	 npip = pi_plus.size();
        npim = pi_minus.size();
        nep = e_plus.size(); 
        nem = e_minus.size();
        nkp = k_plus.size();
        nkm = k_minus.size();
        nmup = mu_plus.size();
        nmum = mu_minus.size();
        npm = p_plus.size();
        npp = p_minus.size();
*/
        n_tks = tk_short.size();
        n_tlamda = tLamda.size();
        n_tlamdabar = tLamdabar.size();
        n_tgamma = tgamma.size();
        n_tgamma1 = tgamma1.size();

	//ntot = npip+npim+nep+nem+nkp+nkm+nmup+nmum+npm+npp ;

//	if( frecB.tag_id() != -521)
//	continue;
	float eecl = 0;
	float ek = 0 ;
//        std::cout<<tk_plus.size()<<std::endl;
//        std::cout<<k_plus.size()<<"kplus"<<std::endl;

///////////////////////////////////k/////////////////////////////////////////////
 for(std::vector<Particle>::iterator i5 =tk_plus.begin(); i5 != tk_plus.end(); i5++)

{
//	HepLorentzVector pnu;
//	HepLorentzVector pantinu;
	std::cout<<"start combination"<<std::endl;
//        std::cout<<(*i5).momentum().p().mag()<<std::endl;

//	Particle nul(pnu , ptype_nu_l);
//        Particle antinul(pantinu , ptype_antinu_l);
        Particle Bsig((*i5).momentum().p(), ptype_B_plus_cand);
        Bsig.relation().append((*i5));
//        Bsig.relation().append(nul);
//        Bsig.relation().append(antinul);
        B_cand1.push_back(Bsig);
}

///////////////////////////////////k/////////////////////////////////////////////
 for(std::vector<Particle>::iterator i =B_cand1.begin(); i != B_cand1.end(); i++)
	
{

	tnpip = tpi_plus.size();
	tnpim = tpi_minus.size();
	tnep = te_plus.size();
	tnem = te_minus.size();
	tnkp = tk_plus.size();
	tnkm = tk_minus.size();
	tnmup = tmu_plus.size();
	tnmum = tmu_minus.size();
	tnpm = tp_plus.size();
	tnpp = tp_minus.size();
	

	npip = pi_plus.size();
        npim = pi_minus.size();
        nep = e_plus.size();
        nem = e_minus.size();
        nkp = k_plus.size();
        nkm = k_minus.size();
        nmup = mu_plus.size();
        nmum = mu_minus.size();
        npm = p_plus.size();
        npp = p_minus.size();

     	tnremaining = tnpip+tnpim+tnep+tnem+tnkp+tnkm+tnmup+tnmum+tnpm+tnpp - 1 ;
        nremaining = npip+npim+nep+nem+nkp+nkm+nmup+nmum+npm+npp - 1;

        TagB_tpl->column( "nb_trackremaining"    , nremaining         );//new
	TagB_tpl->column( "tnb_trackremaining"    , tnremaining         );//new
	TagB_tpl->column( "tnb_pi0remaining"    , npi0         );//new


//			std::cout<<"start tk"<<std::endl;
	double sum_ecl = 0;
        for(std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s !=eclaux_mag.end();s++ )
        {
        //               std::cout<<"start eclaux"<<std::endl;
            Mdst_ecl_aux& ch = *s;
            int flag = 0;
            int nfinal = Btag.relation().nFinalStateParticles(); // the number of final particle used to tagB
            for(int k = 0 ; k < nfinal; k++)
            {
                Particle child = Btag.relation().finalStateParticle(k);
                if ( child.mdstCharged() )
                {
                    if ( child.mdstCharged().trk() == ch.trk())
                    {
                        flag=1;
                    }
                }
                if (child.mdstGamma())
                {
                    if (child.mdstGamma().ecl().get_ID()==ch.get_ID())
                    {
                        flag=1;
                    }
                }
            }
            if((*i).child(0).mdstCharged().trk() == (*s).trk())
            {
            flag=1;
            }
            if(flag) continue;
            for(std::vector<Mdst_ecl>::iterator p = ecl_mag.begin(); p!=ecl_mag.end(); p++)
            {
		       std::cout<<"start ecl"<<std::endl;
              if ((*s).get_ID()==(*p).get_ID())
              	  {
                  //forward e>0.1
                  if ((*p).theta() < 0.547 && (*p).energy() > 0.1 )
                      sum_ecl = sum_ecl + (*p).energy();
                  //barrel e>0.05
                  if ((*p).theta() > 0.562 && (*p).theta() < 2.246 && (*p).energy()>0.05 )
                      sum_ecl = sum_ecl + (*p).energy();
                  //backward e>0.15
	          if ((*p).theta() > 2.246 && (*p).energy() > 0.15)
                      sum_ecl = sum_ecl + (*p).energy();
      		 }
            }
        }//ecl_aux end
          //             std::cout<<"start full"<<std::endl;

	 TagB_tpl->column("RatioToOtherType" , RatioToOtherType);
         TagB_tpl->column("RaatioToSecondBest" , RaatioToSecondBest);
	 //////////////////////K information///////////////////////
//	 TagB_tpl->column("kpx" , (*i).child(0).px());
//         TagB_tpl->column("kpy" , (*i).child(0).py());
//         TagB_tpl->column("kpz" , (*i).child(0).pz());
//         TagB_tpl->column("kmass" , (*i).child(0).p().mag());
	 TagB_tpl->column("sumecl",sum_ecl);
//	 TagB_tpl->column("kE" , (*i).child(0).e());
	//////
        double totalecl = 0;
        for(std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s !=eclaux_mag.end();s++ )
        {
            Mdst_ecl_aux& ch = *s;
            for(std::vector<Mdst_ecl>::iterator p = ecl_mag.begin(); p!=ecl_mag.end(); p++)
            {
              if ((*s).get_ID()==(*p).get_ID())
                  {
                  //forward e>0.1
                  if ((*p).theta() < 0.547 && (*p).energy() > 0.1 )
                      totalecl = totalecl + (*p).energy();
                  //barrel e>0.05
                  if ((*p).theta() > 0.562 && (*p).theta() < 2.246 && (*p).energy()>0.05 )
                      totalecl = totalecl + (*p).energy();
                  //backward e>0.15
                  if ((*p).theta() > 2.246 && (*p).energy() > 0.15)
                      totalecl = totalecl + (*p).energy();
                 }
            }
	}
	//////
	 double reenergy = P_BEAM.e() - totalecl;
         TagB_tpl->column( "reenergy"    , reenergy          );

////////////////////////Btag///////////////////////////////////////////////////////////
	 // Event ID
         TagB_tpl->column( "ebeam"    , ebeam          );
         TagB_tpl->column( "Tag_RunID", Run            );
         TagB_tpl->column( "Tag_EvtID", Evt            );
         TagB_tpl->column( "Tag_ExpID", Exp            );
         // TagB idhep
         TagB_tpl->column( "tagid"    , frecB.tag_id() );
         // Kinetic variable
         TagB_tpl->column( "tagBmbc"  , frecB.Mbc()    );
         TagB_tpl->column( "tagBde"   , frecB.DeltaE() );
         float BtagE = frecB.DeltaE() + ebeam ;
         TagB_tpl->column( "BtagE"  , BtagE   );
         // unique code for reconstructed D decay channels
         // for B -> D X decay modes only 1 entry is filled
         // for B -> D* ( -> D X) Y decay modes only 1 (D*) and 2 (D) entries are filled
         // for B -> D* ( -> D1 X) D2 decay modes only 1 (D*), 2 (D1) and 3 (D2) entries are filled
         // for B -> D*1 ( -> D1 X) D*2 ( -> D2 Y) decay modes all four entries are filled: D*1, D1, D*2, D2 
         TagB_tpl->column( "decmod"   , frecB.decay()       );
         TagB_tpl->column( "ddecmod1" , frecB.Ddec(0)       );
         TagB_tpl->column( "ddecmod2" , frecB.Ddec(1)       );
         TagB_tpl->column( "ddecmod3" , frecB.Ddec(2)       );
         TagB_tpl->column( "ddecmod4" , frecB.Ddec(3)       );
         TagB_tpl->column( "mcinfo"   , frecB.MCinfo()      );
         TagB_tpl->column( "nFS"      , frecB.nFS()         );
         // default selection (without continuum suppression applied)
         TagB_tpl->column( "best"     , frecB.NBRank()      );
         TagB_tpl->column( "nbout"    , frecB.NBout()       );
         // selection with continuum suppression
         TagB_tpl->column( "csbest"   , frecB.cont_NBRank() );
         TagB_tpl->column( "csnbout"  , frecB.cont_NBout()  );
	// std::cout<<"start miss"<<std::endl;
///////////////////////////////////////////////////////////////////////	
//                   Missing P, Mass in Lab frame                    //
///////////////////////////////////////////////////////////////////////
        HepLorentzVector tagB = Btag.p();
        HepLorentzVector Bsigplab (-tagB.px(),-tagB.py(),-tagB.pz(),tagB.e());
	HepLorentzVector tagBlab = Btag.p();
	HepLorentzVector BEAMLAB(E_HER*sin(cross_angle), 0.0, E_HER*cos(cross_angle)-E_LER, E_HER+E_LER);//BEAM
	HepLorentzVector k_plab((*i).child(0).px(),(*i).child(0).py(),(*i).child(0).pz(),(*i).child(0).e());//sigk
	HepLorentzVector tmissinglab = ( BEAMLAB - tagBlab - k_plab); //total miss
        HepLorentzVector qsqurlab = ( Bsigplab - k_plab); // signal side miss
        double Bsigmass2 = Bsigplab.mag()*Bsigplab.mag();//sigbmass^2
	
	double lab_total_miss = sqrt( tmissinglab.px()*tmissinglab.px() + tmissinglab.py()*tmissinglab.py() + tmissinglab.pz()*tmissinglab.pz() );
	TagB_tpl->column("tmmisslab",missinglab.mag());
	double tmmiss2lab = tmissinglab.mag()*tmissinglab.mag();
	double qsquare = qsqurlab.mag()*qsqurlab.mag();
	double sblab = qsquare / Bsigmass2;
        TagB_tpl->column("sblab",sblab);
        TagB_tpl->column("qsquare",qsquare);
       // TagB_tpl->column("mmiss2lab",mmiss2lab);
	TagB_tpl->column("tmmiss2lab",tmmiss2lab);
	TagB_tpl->column("tpxmisslab",tmissinglab.px());
	TagB_tpl->column("tpymisslab",tmissinglab.py());
	TagB_tpl->column("tpzmisslab",tmissinglab.pz());
	TagB_tpl->column("temisslab",tmissinglab.e());
	TagB_tpl->column("tpmisslab",lab_total_miss);
	TagB_tpl->column("kmalab",k_plab.mag());
	TagB_tpl->column("kpxlab",k_plab.px());
	TagB_tpl->column("kpylab",k_plab.py());
	TagB_tpl->column("kpzlab",k_plab.pz());
	TagB_tpl->column("kelab",k_plab.e());
	double lab_k_p = sqrt( k_plab.px()*k_plab.px() + k_plab.py()*k_plab.py() + k_plab.pz()*k_plab.pz() );
	TagB_tpl->column("kplab",lab_k_p);	
	//std::cout<<"p1"<<std::endl;
	//Lab frame miss cos
	double cosmisslab;
	double angelmisslab;
	tcosmisslabinz = ((tmissinglab.pz())/lab_total_miss);
	tangelmisslabinz = acos((tmissinglab.pz())/lab_total_miss);
	TagB_tpl->column("tcosmizl",tcosmisslabinz);
	//h(*) B+ cos lab
	HepLorentzVector Btagplab (Btag.px(),Btag.py(),Btag.pz(),Btag.e());
	HepLorentzVector all_BEAM_polar ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle),E_HER+E_LER);
	HepLorentzVector P_BEAM_polar = (all_BEAM_polar - Btagplab);//BEAM - BTAG = SINGAL
        double CM_2 = sqrt( P_BEAM_polar.px()*P_BEAM_polar.px() + P_BEAM_polar.py()*P_BEAM_polar.py() + P_BEAM_polar.pz()*P_BEAM_polar.pz() );

	double cossklab = (((*i).child(0).px()*P_BEAM_polar.px())+((*i).child(0).py()*P_BEAM_polar.py())+((*i).child(0).pz()*P_BEAM_polar.pz()))/(lab_k_p*CM_2);
	TagB_tpl->column( "cosklab"   , cosksblab);//new
/////////////////////////////////////////////////////////////////////////////////////
//                            Boost to upsilon 4S rest frame                       //
/////////////////////////////////////////////////////////////////////////////////////
  
         tagB.boost( boost_vector.boostVector() );
         double CM_Btag_p = sqrt( tagB.px()*tagB.px() + tagB.py()*tagB.py() + tagB.pz()*tagB.pz() );
	 double CM_Btag_pz = sqrt( tagB.pz()*tagB.pz() );
//	 double CM_Btag_ppz = tagB.pz()*tagB.px() + tagB.pz()*tagB.py() + tagB.pz()*tagB.pz() ;
	 double cosbcm = tagB.pz()/(CM_Btag_p);
         TagB_tpl->column( "cosbcm"   , cosbcm           );
         TagB_tpl->column( "tagBpcm"   , CM_Btag_p           );
         TagB_tpl->column( "tagbPXcm"   , tagB.px());
         TagB_tpl->column( "tagbPYcm"   , tagB.py());
         TagB_tpl->column( "tagbPZcm"   , tagB.pz());
	HepLorentzVector BEAM(E_HER*sin(cross_angle), 0.0, E_HER*cos(cross_angle)-E_LER, E_HER+E_LER);
        HepLorentzVector k_p((*i).child(0).px(),(*i).child(0).py(),(*i).child(0).pz(),(*i).child(0).e());
	BEAM.boost(boost_vector.boostVector());
        k_p.boost(boost_vector.boostVector() );
        HepLorentzVector totalmissincm = ( BEAM - tagB - k_p);
        HepLorentzVector Bsigpcm (-tagB.px(),-tagB.py(),-tagB.pz(),tagB.e());
	HepLorentzVector sigmisscm = (  B_pcm - k_p );
	double cm_totalmiss_p =sqrt(totalmissincm.px()*totalmissincm.px() + totalmissincm.py()*totalmissincm.py() + totalmissincm.pz()*totalmissincm.pz()) ;
        double bsigpcm =sqrt(Bsigpcm.px()*Bsigpcm.px() + Bsigpcm.py()*Bsigpcm.py() + Bsigpcm.pz()*Bsigpcm.pz()) ;
	double sigmisspcm = sqrt( sigmisscm.px()*sigmisscm.px() + sigmisscm.py()*sigmisscm.py() + sigmisscm.pz()*sigmisscm.pz() );
	double anglebsigandsigmiss = acos((Bsigpcm.px()*sigmisscm.px() + Bsigpcm.py()*sigmisscm.py() + Bsigpcm.pz()*sigmisscm.pz())/(bsigpcm*sigmisspcm));//new
        double cosbsigandsigmiss = (Bsigpcm.px()*sigmisscm.px() + Bsigpcm.py()*sigmisscm.py() + Bsigpcm.pz()*sigmisscm.pz())/(bsigpcm*sigmisspcm);//new

        TagB_tpl->column("angbssmi", anglebsigandsigmiss);//new
        TagB_tpl->column("cosbssmi",cosbsigandsigmiss);//new
        TagB_tpl->column("sigmissp",sigmisscm.mag());
	double sigmmiss2 = sigmisscm.mag()*sigmisscm.mag();
        //TagB_tpl->column("mmiss2CM",mmiss2CM);
//	double mb2 = B_p.mag()*B_p.mag();
        TagB_tpl->column("Bsigmas2",Bsigmass2);
	double sbcm = sigmmiss2 /Bsigmass2  ;
        TagB_tpl->column("sbcm",sbcm);
        TagB_tpl->column("simmis2CM",sigmmiss2);
        TagB_tpl->column("spxmisisC",sigmisscm.px());
        TagB_tpl->column("spymissC",sigmisscm.py());
        TagB_tpl->column("spzmissC",sigmisscm.pz());
        TagB_tpl->column("semissC",sigmisscm.e());
        TagB_tpl->column("spmissC",sigmisspcm);
	TagB_tpl->column("mkcm",k_p.mag());
	TagB_tpl->column("kpxCM",k_p.px());
	TagB_tpl->column("kpyCM",k_p.py());
	TagB_tpl->column("kpzCM",k_p.pz());
	TagB_tpl->column("keCM",k_p.e());
	double CM_k_p = sqrt( k_p.px()*k_p.px() + k_p.py()*k_p.py() + k_p.pz()*k_p.pz() );
	TagB_tpl->column("kpCM",CM_k_p);
	//HepLorentzVector B_p (-tagB.px(),-tagB.py(),-tagB.pz(),tagB.e());
        HepLorentzVector miss_b(sigmisscm);
	HepLorentzVector k_b(k_p);
        miss_b.boost(-B_p.boostVector());
	k_b.boost(-B_p.boostVector());
        TagB_tpl->column("mmissb",miss_b.mag());
        double mmiss2b = miss_b.mag()*miss_b.mag();
	double b_k_pmiss = sqrt( miss_b.px()*miss_b.px() + miss_b.py()*miss_b.py() + miss_b.pz()*miss_b.pz() );
        TagB_tpl->column("mmiss2b",mmiss2b);
        TagB_tpl->column("pxmissb",miss_b.px());
        TagB_tpl->column("pymissb",miss_b.py());
        TagB_tpl->column("pzmissb",miss_b.pz());
        TagB_tpl->column("emissb",miss_b.e());
        TagB_tpl->column("pmissb",b_k_pmiss);
        double b_k_p = sqrt( k_b.px()*k_b.px() + k_b.py()*k_b.py() + k_b.pz()*k_b.pz() );
        TagB_tpl->column("kpb", b_k_p);
        TagB_tpl->column("kpxb",k_b.px());
        TagB_tpl->column("kpyb",k_b.py());
        TagB_tpl->column("kpzb",k_b.pz());
	TagB_tpl->column("keb",k_b.e());

////////////////////////////////Exkfitter////////////////////////////////
/*        std::cout<<"start Exkfitter"<<std::endl;

	Mdst_charged Mdst_k =(*i).mdstCharged();
	ExKFitterParticle KF_1(Mdst_k, 3);
        ExKFitterVertex B_Vertex(IP,IPerr);
	ExKFitterParticle B;
        B.LinkParticle(&KF_1);
        B.LinkVertex(&B_Vertex);
        ExKFitterConstrain con;
        con.SetVertexConstrain();
        con.LinkParticle(&KF_1);
        con.LinkVertex(&B_Vertex);

        ExKFitter Core;
        Core.LinkConstrain(&con);
        int ret = Core.Minimize();
        float chisqExK = Core.Chisq();
        float dof_exk = Core.N_DegreeOfFreedom();
	HepPoint3D bvertex = B_Vertex.Vertex();
        
	if(ret==0)
        {
        B.Update();
        }


        TagB_tpl->column("chisqexk",chisqExK/dof_exk);
//                  TagB_tpl->column("sum_eclinmissdir",fabs(angelmisslab-(*p1).theta()));//new
*/
//////////////////////////////////dz///////////////////////////////////////////////////////
	HepPoint3D Btag_init;
	HepPoint3D bvertex;
        HepPoint3D bvertexsig;
        float vchisq_tag;
        int DOF_tag;
        float pvalue_tag;
       std::cout<<"start setx"<<std::endl;
//	Btag_init.setX(IP.x() + Btag.p().px()/(*i).p().rho());
//        Btag_init.setY(IP.y() + Btag.p().py()/(*i).p().rho());
//        Btag_init.setZ(IP.z() + Btag.p().pz()/(*i).p().rho());
	int IDp = 0;
        int ID = 0;
	int flagch = -1;
	int ncharge = 0;
	int vv = 0 ;
	int nfinal = Btag.relation().nFinalStateParticles(); // the number of final particle used to tagB
        for(int k = 0 ; k < nfinal; k++)
	  {
            Particle child = Btag.relation().finalStateParticle(k);
            if ( child.mdstCharged() )
              {
	        ncharge ++;
              }
          }
	if (ncharge>1)
	flagch = 1;
	std::cout<<"start flagch"<<std::endl;
	if(flagch)
	{	
	  kvertexfitter kvf;

//          ExKFitterVertex B_Vertex(IP,IPerr);
//          ExKFitterParticle B;
//          ExKFitterConstrain con;
//          B.LinkVertex(&B_Vertex);
//          con.SetVertexConstrain();
//          con.LinkVertex(&B_Vertex);
              //   std::vector<Mdst_charged> bfinalnumber;
	//	std::vector<ExKFitterParticle> btagf;
        //        std::vector<Particle> bfinal;
	  Mdst_charged bfinalnumber[nfinal];
	  ExKFitterParticle btagf[nfinal];
	  Particle bfinal[nfinal];
	  for(int k2 = 0 ; k2 < nfinal; k2++)
            {
              std::cout<<"start k2"<<std::endl;
	      int expid = 0;
              Particle bfinal_tmp = Btag.relation().finalStateParticle(k2);
	      bfinal[k2]=bfinal_tmp;
	      int nunp = -1 ;
              std::cout<<"start k3"<<std::endl;
	      if(!(bfinal[k2].mdstCharged()))
	      continue;
	      addTrack2fit(kvf,  bfinal[k2]);
	//	char bfinalnumber[32];
        //        sprintf (bfinalnumber,"%s%d","bfinal",k2+1);//bnd1,bnd2,bnd3....
  //             IDp = bfinal[k2].mdstCharged().get_ID();
//		Mdst_charged F_tmp = bfinal[k2].mdstCharged();
//		bfinalnumber[k2] = F_tmp;
	//	bfinalnumber[k2] = bfinal.mdstCharged();
/*		
		if(Btag.relation().finalStateParticle(k2).pType()==ptype_e_plus||Btag.relation().finalStateParticle(k2).pType()==ptype_e_minus)
		ID==0;
                else if(Btag.relation().finalStateParticle(k2).pType()==ptype_mu_plus||Btag.relation().finalStateParticle(k2).pType()==ptype_mu_minus)
                ID==1;
                else if(Btag.relation().finalStateParticle(k2).pType()==ptype_pi_plus||Btag.relation().finalStateParticle(k2).pType()==ptype_pi_minus)
                ID==2;
                else if(Btag.relation().finalStateParticle(k2).pType()==ptype_k_plus||Btag.relation().finalStateParticle(k2).pType()==ptype_k_minus)
                ID==3;
                else if(Btag.relation().finalStateParticle(k2).pType()==ptype_p_plus||Btag.relation().finalStateParticle(k2).pType()==ptype_p_minus)
                ID==4;

*/		
             //   char btagf[32];
 //		ExKFitterParticle btagf_tmp(F_tmp, IDp);
//		btagf[k2] = btagf_tmp;
//		vv++;	
                  //      std::cout<<IDp<<std::endl;

		//	std::cout<<ID<<std::endl;
                  //      std::cout<<"start k5"<<std::endl;
             //  sprintf (btagf,"%s%d","btagf",k2+1);//bnd1,bnd2,:bnd3....
//			for(int i =0 ; i<vv ; i++)
//			{	
        	//	B.LinkParticle(&btagf[k2]);
        	//	con.LinkParticle(&btagf[k2]);
		//	}

             std::cout<<"start k6"<<std::endl;
            }		
	  unsigned err = kvf.fit(); // do "fitting"
	  if(err == 0) // success.
          {
	    vchisq_tag = kvf.chisq();
	    DOF_tag = kvf.dgf();
	    pvalue_tag = Chisq2Confi(DOF_tag , vchisq_tag);
	  }          
	  else 
	  {
	    vchisq_tag = -1;
	    DOF_tag = -1;
            pvalue_tag=-1;
	  }
	    bvertex = kvf.get_vertex();
            std::cout<<"start k6.9"<<std::endl;
//	 ExKFitter Core;
//                Core.LinkConstrain(&con);
  //              int ret = Core.Minimize();
    //            float chisqExK = Core.Chisq();
      //          float dof_exk = Core.N_DegreeOfFreedom();
        //        bvertex = B_Vertex.Vertex();
	  std::cout<<"start k7"<<std::endl;

          //      if(ret==0)
            //    {
              //  B.Update();
               // }

	}	

        double dr,dz;
        const Mdst_charged* ch1 = &(*i).child(0).mdstCharged();
	Mdst_charged Mdst_k = (*i).child(0).mdstCharged();
	Particle sig_k = (*i).child(0);
/*        GetImpactParameters(ch1,&dr,&dz,3);
	ExKFitterParticle k1(Mdst_k,3);
        ExKFitterParticle B_sig;
        ExKFitterVertex sig_Vertex(IP,IPerr);
	ExKFitterParticle IPtube(IpProfile::ip_tube());
	ExKFitterConstrain con1;
	B_sig.LinkParticle(&IPtube);
        B_sig.LinkParticle(&k1);
        B_sig.LinkVertex(&sig_Vertex);
        con1.SetVertexConstrain();
	con1.LinkParticle(&IPtube);
	con1.LinkParticle(&k1);
        con1.LinkVertex(&sig_Vertex);
	ExKFitter Core;
	Core.LinkConstrain(&con1);
	int ret1 = Core.Minimize();

        float chisqExK = Core.Chisq();
        float dof_exk = Core.N_DegreeOfFreedom();
	HepPoint3D B_sigvertex = sig_Vertex.Vertex();
        if(ret1==0)
          {B_sig.Update();}*/
        kvertexfitter kvfsig;
        addTrack2fit(kvfsig,sig_k);
    	const int use_nevtdep_ip = 1;
        kvfsig.initialVertex(IpProfile::position(use_nevtdep_ip));
	addTrack2fit(kvfsig,IpProfile::ip_tube(use_nevtdep_ip));
        unsigned err = kvfsig.fit(); // do "fitting"
        float vchisq_sig;
	int DOF_sig;
	float pvalue_sig;
        if(err == 0) // success.
        { 
	  vchisq_sig = kvfsig.chisq();
	  DOF_sig = kvfsig.dgf();
          pvalue_sig = Chisq2Confi(DOF_sig , vchisq_sig);}
        else 
	{
	  vchisq_sig = -1;
	  DOF_sig = -1;
	  pvalue_sig =-1;
	}
        bvertexsig = kvfsig.get_vertex();
        dr = sqrt( pow(bvertexsig.x(),2) + pow(bvertexsig.y(),2) );	


	double DistTOOtherBdz =  abs(bvertexsig.z())  - abs(bvertex.z());
	TagB_tpl->column("DistTOOtherBdz",DistTOOtherBdz);
        TagB_tpl->column("kdz",bvertexsig.z());
        TagB_tpl->column("kdr", dr );
	TagB_tpl->column("vchisq_sig",vchisq_sig);
        TagB_tpl->column("DOF_sig",DOF_sig);
        TagB_tpl->column("pvalue_sig",pvalue_sig);
        TagB_tpl->column("vchisq_tag",vchisq_tag);
        TagB_tpl->column("DOF_tag",DOF_tag);
        TagB_tpl->column("pvalue_tag",pvalue_tag);
	double Distsign = pvalue_sig -pvalue_tag;
        TagB_tpl->column("Distsign",Distsign);
////////////////////////////////////////////////////////////////////////////////////////////
//                                ECL in miss direction                                   //
///////////////////////////////////////////////////////////////////////////////////////////
        HepLorentzVector EclVirtualVec  (1.0,1.0,1.0,1.0);
        double sum_eclinmissdir = 0;
        for(std::vector<Mdst_ecl_aux>::iterator s1 = eclaux_mag.begin(); s1 !=eclaux_mag.end();s1++ )
          {
          for(std::vector<Mdst_ecl>::iterator p1 = ecl_mag.begin(); p1!=ecl_mag.end(); p1++)
            {
              if ((*s1).get_ID()==(*p1).get_ID())
                {
                  //forward e>0.1
                  if ((*p1).theta() < 0.547 && (*p1).energy() > 0.1 && fabs(missinglab.angle(EclVirtualVec.vect()))<0.1 )	
                    sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
                  //barrel e>0.05
                  if ((*p1).theta() > 0.562 && (*p1).theta() < 2.246 && (*p1).energy()>0.05&& fabs(missinglab.angle(EclVirtualVec.vect()))< 0.1 )
                    sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
                  //backward e>0.15
                  if ((*p1).theta() > 2.246 && (*p1).energy() > 0.15&& fabs(missinglab.angle(EclVirtualVec.vect()))< 0.1)
                    sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
                }
		EclVirtualVec.setRho(1.0);EclVirtualVec.setTheta((*p1).theta());EclVirtualVec.setPhi((*p1).phi());
                        
//	                  b_tpl->column("coneangle_ecl",fabs(missinglab.angle(EclVirtualVec.vect())));//new
//			  b_tpl->dumpData();
//	                 *status = 1;
            }
          }//ecl_aux end
	TagB_tpl->column("sum_eclinmissdir",sum_eclinmissdir);//new
/////////////////////////////////////////////////////////////////////////////
//                           k_sfw variables                               //
////////////////////////////////////////////////////////////////////////////     
       k_sfw ksfw_obj((*i));
      const double ksfw(ksfw_obj.fd());
      const int iksfw(ksfw_obj.i_mm2());
      TagB_tpl->column("imm2",iksfw);
      double miss(ksfw_obj.mm2());
      TagB_tpl->column("mmm2",miss);
      double et(ksfw_obj.e_t());
      TagB_tpl->column("et",et);
      double H0oo(ksfw_obj.Hoo(0));
      double H1oo(ksfw_obj.Hoo(1));
      double H2oo(ksfw_obj.Hoo(2));
      double H3oo(ksfw_obj.Hoo(3));
      double H4oo(ksfw_obj.Hoo(4));
      TagB_tpl->column("H0oo",H0oo);
      TagB_tpl->column("H1oo",H1oo);
      TagB_tpl->column("H2oo",H2oo);
      TagB_tpl->column("H3oo",H3oo);
      TagB_tpl->column("H4oo",H4oo);
      double H0son(ksfw_obj.Hso_n(0));
      double H1son(ksfw_obj.Hso_n(1));
      double H2son(ksfw_obj.Hso_n(2));
      double H3son(ksfw_obj.Hso_n(3));
      double H4son(ksfw_obj.Hso_n(4));
      TagB_tpl->column("H0son",H0son);
      TagB_tpl->column("H1son",H1son);
      TagB_tpl->column("H2son",H2son);
      TagB_tpl->column("H3son",H3son);
      TagB_tpl->column("H4son",H4son);
      double H0soc(ksfw_obj.Hso_c(0));
      double H1soc(ksfw_obj.Hso_c(1));
      double H2soc(ksfw_obj.Hso_c(2));
      double H3soc(ksfw_obj.Hso_c(3));
      double H4soc(ksfw_obj.Hso_c(4));
      TagB_tpl->column("H0soc",H0soc);
      TagB_tpl->column("H1soc",H1soc);
      TagB_tpl->column("H2soc",H2soc);
      TagB_tpl->column("H3soc",H3soc);
      TagB_tpl->column("H4soc",H4soc);
      double H0som(ksfw_obj.Hso_m(0));
      double H1som(ksfw_obj.Hso_m(1));
      double H2som(ksfw_obj.Hso_m(2));
      double H3som(ksfw_obj.Hso_m(3));
      double H4som(ksfw_obj.Hso_m(4));
      TagB_tpl->column("H0som",H0som);
      TagB_tpl->column("H1som",H1som);
      TagB_tpl->column("H2som",H2som);
      TagB_tpl->column("H3som",H3som);
      TagB_tpl->column("H4som",H4som);
        //for shape variables
        Vector4 OtherBVector;
        float R2, spher, cos_thr, cos_thp, par_sfw[13];
	
	 int vertexflag;
              HepPoint3D overtex;
        std::cout<<"start shape"<<std::endl;

//        shape((*i), R2,  spher, cos_thr, cos_thp, par_sfw, OtherBVector, bvertex, vertexflag, overtex, Btag);
//        TagB_tpl->column("bvertexx",bvertex.x());
//        TagB_tpl->column("bvertexy",bvertex.y());
//       TagB_tpl->column("bvertexz",bvertex.z());
//        TagB_tpl->column("overtexx",overtex.x());
//        TagB_tpl->column("overtexy",overtex.y());
 //       TagB_tpl->column("overtexz",overtex.z());
//        TagB_tpl->column("exfit",ret);
//        TagB_tpl->column("chisqExK",chisqExK);

//        TagB_tpl->column("pmiss_tag",P_miss_tag.rho());
//        TagB_tpl->column("emiss_tag",P_miss_tag.e());

        // cosine of B and beam dir. the angle between B beam and P_beam but in rest frame or lab frame?
//        float cosb = Vector3(b_p.vect()).dot(P_BEAM.vect());
 //       cosb = cosb/(Vector3(b_p.vect()).mag()*P_BEAM.vect().mag());

        //fill ntuple of shape variables
//        TagB_tpl->column("vtflag",vertexflag);
//        if( !vertexflag)
//        {
//        double deltaZ=bvertex.z()-overtex.z();
//        TagB_tpl->column("deltaz",deltaZ);
//        }


//        TagB_tpl->column("R2", R2);
//        TagB_tpl->column("costhr", float(cos_thr));
//        TagB_tpl->column("costhp", float(cos_thp));
//        TagB_tpl->column("spher",  float(spher));
   //     b_tpl->column("cosb", cosb);


//	        std::cout<<"end Exkfitter"<<std::endl;


///////////////////////////////////////////////////////////////////////////////////////
//                      MC TRuth(I don't know it work or not)                        //
///////////////////////////////////////////////////////////////////////////////////////

        if(MCstatus == 1)
        {

            const Mdst_charged ch1 = (*i).child(0).mdstCharged();
	    Gen_hepevt Evtch1 = get_hepevt(ch1);
            TagB_tpl->column("hi1",Evtch1.idhep());
	
            Gen_hepevt EvtP1;
            Gen_hepevt EvtGP1;
            Gen_hepevt EvtGGP1;
   if(Evtch1.mo(0))
            {
                EvtP1 = gen_mgr[Evtch1.mo(0)-1];
                TagB_tpl->column("mo1",EvtP1.idhep());


                if(EvtP1.mo(0))
                {
                    EvtGP1 = gen_mgr[EvtP1.mo(0)-1];
                    TagB_tpl->column("gmo1",EvtGP1.idhep());


                    if(EvtGP1.mo(0))
                    {
                        EvtGGP1 = gen_mgr[EvtP1.mo(0)-1];
                        TagB_tpl->column("ggmo1",EvtGGP1.idhep());
                    }
                }
            }


            int hindex = 0 ;

            if(Evtch1.idhep()==321&&EvtP1.idhep()==521)
                {	
			hindex = 1;		
 		}

		else
			{
				hindex = 0;
			}			
	
	           TagB_tpl->column("hindex",hindex);

		}//truth event  end
//        std::cout<<"start truth event"<<std::endl;

/////////////////////////////////trace form top to down//////////////
	int nbpd=0,nbnd=0;
        for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin();i != gen_mgr.end(); i++)
		  {
                  if ((*i).idhep()==521||(*i).idhep()== 511)//511==B0
                  {
                    int da1=(*i).da(0),da2=(*i).da(1);

                    if ((*i).idhep() == 521)
			 TagB_tpl->column("charged",1);
                    else TagB_tpl->column("charged",0);
                    nbpd=(da2-da1+1);
                    TagB_tpl->column("nbpd",nbpd);
                      for(int start=0;start<(da2-da1+1);start++)
                      {
                        Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
                        char bpdnumber[32];
                        sprintf (bpdnumber,"%s%d","bpd",start+1);
                        TagB_tpl->column(bpdnumber,Evda.idhep());
                      }
                   }
                   else if((*i).idhep()== -521 || (*i).idhep()== -511)
                   {
                        int da1=(*i).da(0),da2=(*i).da(1);

                          if ((*i).idhep() == -521) 
			       TagB_tpl->column("charged",1);
                          else TagB_tpl->column("charged",0);
                        nbnd=(da2-da1+1);
                        TagB_tpl->column("nbnd",nbnd);
                         for(int start=0;start<(da2-da1+1);start++)
                         {
                           Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
                           char bndnumber[32];
                           sprintf (bndnumber,"%s%d","bnd",start+1);//bnd1,bnd2,bnd3....
                           TagB_tpl->column(bndnumber,Evda.idhep());
                         }
                   }
        	}

		  TagB_tpl->dumpData();
    	   	 *status = 1;
             }//k_plus end/////////////////////////////////////////////////////////////////////////
	  }//iterator frecon_it
	  return;
}// end   void ana_fullreconk::event( BelleEvent* evptr, int* status )

void ana_fullreconk::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13],Vector4 &otherB_P, HepPoint3D &bvertex, int & vertexflag, HepPoint3D & overtex)
{
   TagVK ver; // Vertex Reconstructor with kvertexfitter
   ver.setdefault(b,bvertex);
   //ver.useKsVeto();
   ver.dontUseIPforFit();
   ver.useTubeforFit();

  double E_HER=BeamEnergy::E_HER();
  double E_LER=BeamEnergy::E_LER();
  double cross_angle=BeamEnergy::Cross_angle();
//        std::cout<<"start shape1"<<std::endl;

  //Particle Type
  Ptype ptype_gamma("GAMM");
  static Vector4 P_BEAM ( -E_HER*sin(cross_angle), 0.0, E_LER-E_HER*cos(cross_angle),
                          E_HER+E_LER );

  std::vector<Particle> char_plus;
  std::vector<Particle> char_minus;
  std::vector<Particle> g; //gamma particle

  Ptype ptype_g("GAMM");
  Ptype ptype_char_plus("PI+");
  Ptype ptype_char_minus("PI-");

  std::vector<Particle> candiBfinalParticle;
  std::vector<Vector4> allFinal;
  std::vector<Vector4> candiBfinal;
  std::vector<Vector4> candiBfinalgamm;
  std::vector<Vector4> otherBfinal;
  std::vector<Vector4> Pmiss;
//fill all charged particle to charged bank
//into particle list
  Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
  for(std::vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); i++)
  {
    if ( (*i).charge() > 0.)
    {
      Particle tmp(*i, ptype_char_plus);
      char_plus.push_back(tmp);
    }
    else
    {
      Particle tmp(*i, ptype_char_minus);
      char_minus.push_back(tmp);
    }
  }
  //      std::cout<<"start shape2"<<std::endl;

//fill gamma list from Mdst_gamma
//into particle list
  Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
  for (std::vector<Mdst_gamma>::iterator i = gamma_mag.begin();
       i != gamma_mag.end(); i++)
  {
    Particle tmp(*i);
    g.push_back(tmp);
  }
//for super fox-wolfram
//cout << "no of B final particle = " << b.relation().nFinalStateParticles()<<endl;

  for( int i=0;i< b.relation().nFinalStateParticles();
       i++)
  {
    candiBfinalParticle.push_back(b.relation().finalStateParticle(i));
    //cout << " i = " << i << " mass = " <<  b.relation().finalStateParticle(i).mass() <<  endl;
    Vector4 temp_P4 = b.p();
    temp_P4.boost(P_BEAM.boostVector());
    candiBfinal.push_back(temp_P4);
  }        std::cout<<"start shape3"<<std::endl;

  //by one daughter from b (now omega)
  //cout << "omega final state no ="
  //     << b.relation().child(0).relation().nFinalStateParticles() << endl;
  for(int i=0; i<b.relation().child(0).relation().nFinalStateParticles();i++)
  {
    Vector4 temp_P4 =
      b.relation().child(0).relation().finalStateParticle(i).p();
    temp_P4.boost(P_BEAM.boostVector());
    candiBfinalgamm.push_back(temp_P4);
  }

  int numver=0;
  for( std::vector<Particle>::iterator it1=char_plus.begin();
       it1!=char_plus.end(); it1++ )
  {
    Vector4 temp_P4 = (*it1).p();
    temp_P4.boost(P_BEAM.boostVector());
    allFinal.push_back(temp_P4);
    int notBDauFlag=1;
    for(std::vector<Particle>::iterator j= candiBfinalParticle.begin();
        j!= candiBfinalParticle.end();j++)
      if(checkMultiUse(*it1,*j)){ notBDauFlag=0; break;}
      if(notBDauFlag)
      {
      otherBfinal.push_back(temp_P4);
      otherB_P=otherB_P+temp_P4;
      Particle& tmp=*it1;
      ver.push_back(&tmp);
      numver++;
        std::cout<<"start shape4"<<std::endl;

      }
  }
  for ( std::vector<Particle>::iterator it1=char_minus.begin();
        it1!=char_minus.end(); it1++ )
  {
    Vector4 temp_P4 = (*it1).p();
    temp_P4.boost(P_BEAM.boostVector());
    allFinal.push_back(temp_P4);
    int notBDauFlag=1;
    for(std::vector<Particle>::iterator j= candiBfinalParticle.begin();
        j!= candiBfinalParticle.end();j++)
      if(checkMultiUse(*it1,*j)){ notBDauFlag=0; break;}
      if(notBDauFlag)
      {
      otherBfinal.push_back(temp_P4);
      otherB_P=otherB_P+temp_P4;
      //cout <<" "<<(*it1).relation().mdstCharged().get_ID();
      Particle& tmp=*it1;
        std::cout<<"start shape5"<<std::endl;

      ver.push_back(&tmp);
      numver++;
      }
  }

  vertexflag=ver.fit();
  if(!vertexflag)
  {
    overtex = ver.vtx();
    //cout<<numver<<":"<<ver.vtx()<<" "<<ver.used_particles().size()<<endl;
  }
  for ( std::vector<Particle>::iterator it1=g.begin();
        it1!= g.end(); it1++ )
  {
    Vector4 temp_P4 = (*it1).p();
    temp_P4.boost(P_BEAM.boostVector());
    allFinal.push_back(temp_P4);
    int notBDauFlag=1;
    for(std::vector<Particle>::iterator j= candiBfinalParticle.begin();
        j!= candiBfinalParticle.end();j++)
      if(checkMultiUse(*it1,*j)){ notBDauFlag=0; break;}
      if(notBDauFlag)
      {
      otherBfinal.push_back(temp_P4);
      otherB_P=otherB_P+temp_P4;
      //cout <<" "<<(*it1).relation().mdstGamma().get_ID();
      }
        std::cout<<"start shape6"<<std::endl;

  }

//  float ebeam = Benergy();
      double ebeam = BeamEnergy::E_beam_corr();
  HepLorentzVector boost_vector(E_HER*sin(cross_angle), 0.0,-E_LER+E_HER*cos(cross_angle),E_HER+E_LER);
  Vector4 ExpectedOtherBVector(0.,0.,0.,ebeam);
  Vector4 OtherBVector;
  OtherBVector.boost(boost_vector.boostVector());
  Vector4 temp_Pmiss = ExpectedOtherBVector - OtherBVector;
  Pmiss.push_back(temp_Pmiss);

  //cout << endl;
  //cout << " no of ch+ = "<< char_plus.size();
  //cout << " no of ch- = "<< char_minus.size();
  //cout << " no of gam = "<< g.size() << endl;

  if((candiBfinal.size()+otherBfinal.size()) != allFinal.size())


  //Start to fill shape variables
  Vector3 candiBthr=thrust(candiBfinal);
  Vector3 otherBthr=thrust(otherBfinal);
  //Fox Wolfram Moment
        std::cout<<"start shape7"<<std::endl;

  SuperFoxWolfram foxWolfram;
  SuperFoxWolfram sfwSO;
  SuperFoxWolfram sfwS1O;
  SuperFoxWolfram sfwOO;
  SuperFoxWolfram sfwS3O;
  //SuperFoxWolfram sfwSS;
  foxWolfram.fill(allFinal,allFinal);
  sfwSO.fill(candiBfinal,otherBfinal);
  //sfwSS.fill(candiBfinal,candiBfinal);

  par_sfw[0] = (float)foxWolfram.R(2);
  par_sfw[1] = (float)sfwSO.R(1);
  par_sfw[2] = (float)sfwSO.R(2);
  par_sfw[3] = (float)sfwSO.R(3);
  par_sfw[4] = (float)sfwSO.R(4);

  sfwOO.fill(otherBfinal,otherBfinal);
  par_sfw[5] = (float)sfwOO.R(1);
  par_sfw[6] = (float)sfwOO.R(2);
  par_sfw[7] = (float)sfwOO.R(3);
  par_sfw[8] = (float)sfwOO.R(4);

  sfwS1O.fill(candiBfinalgamm,otherBfinal);
  par_sfw[9] = (float)sfwS1O.R(1);
  par_sfw[10] = (float)sfwS1O.R(2);
  par_sfw[11] = (float)sfwS1O.R(3);
  par_sfw[12] = (float)sfwS1O.R(4);

  //for shape variables
  float H0,R1,R3,R4;
  int ntrk=char_plus.size()+char_minus.size()+g.size(),trkc=0;
  //int ntrk=g_mag.size()+charged_mag.size(),trkcount=0;
  float ptrk[ntrk*3];
  int itrk[ntrk];
  //float ebeam = Benergy();
        std::cout<<"start shape8"<<std::endl;

  for (std::vector<Particle>::iterator i = g.begin();
       i != g.end(); i++)
  {
    Vector4 g_P= (*i).p();
    g_P.boost(P_BEAM.boostVector());
    //Vector3 g_boost = g_P;

    ptrk[trkc*3]=g_P.px();
    ptrk[trkc*3+1]=g_P.py();
    ptrk[trkc*3+2]=g_P.pz();
    trkc++;
  }

  for (std::vector<Particle>::iterator i = char_plus.begin();
       i != char_plus.end(); i++)
  {
    Vector4 charged_P = (*i).p();
    charged_P.boost(P_BEAM.boostVector());
    //Vector3 charged_boost = charged_P;

    ptrk[trkc*3]=charged_P.px();
    ptrk[trkc*3+1]=charged_P.py();
    ptrk[trkc*3+2]=charged_P.pz();
    trkc++;
  } //end for Particle_List

  for (std::vector<Particle>::iterator i = char_minus.begin();
       i != char_minus.end(); i++)
  {
    Vector4 charged_P = (*i).p();
    charged_P.boost(P_BEAM.boostVector());
    //Vector3 charged_boost = charged_P;

    ptrk[trkc*3]=charged_P.px();
    ptrk[trkc*3+1]=charged_P.py();
    ptrk[trkc*3+2]=charged_P.pz();
    trkc++;
  } //end for Particle_List
  //end of shape variables
        std::cout<<"start shape9"<<std::endl;

  fwjet2(ntrk,ptrk,ebeam,&H0,&R1,&R2,&R3,&R4);

  // variables for thrust

  std::vector<Vector4> ptl;
  std::vector<Vector4> ptlb;
  Vector3 Bthr;


  for (std::vector<Particle>::iterator j = g.begin();
       j != g.end(); j++)
  {
    int mask=0;
    for (int i=0; i<b.relation().nFinalStateParticles(); i++)
    {
      if (b.relation().finalStateParticle(i).pType() ==
          ptype_gamma)
      {
        if ( (*j).relation().mdstGamma().get_ID() ==
             b.relation().finalStateParticle(i).mdstGamma().get_ID() )
          mask=1;

      }
    }
      Vector4 g_P= (*j).p();
      g_P.boost(P_BEAM.boostVector());
      Vector4 tmpptl(g_P.px(),g_P.py(),g_P.pz(),
                     g_P.e());

        std::cout<<"start shape10"<<std::endl;

    if (mask==1)
    {
      ptlb.push_back(tmpptl);
    }
    else
    {
      ptl.push_back(tmpptl);
    }
  }


  for (std::vector<Particle>::iterator j = char_plus.begin();
       j != char_plus.end(); j++)
   {
    int mask=0;
    for (int i=0; i<b.relation().nFinalStateParticles(); i++)
    {
      if (b.relation().finalStateParticle(i).pType() !=
          ptype_gamma)
      {
        if ( (*j).relation().mdstCharged().get_ID() ==
             b.relation().finalStateParticle(i).mdstCharged().get_ID() )
          mask=1;
      }
    }

    Vector4 charged_P = (*j).p();
    charged_P.boost(P_BEAM.boostVector());
    Vector4 tmpptl(charged_P.px(),charged_P.py(),charged_P.pz(),
                   charged_P.e());

    if (mask==1)
    {
    ptlb.push_back(tmpptl);
    }
    else
    {
    ptl.push_back(tmpptl);
    }
  }
  for (std::vector<Particle>::iterator j = char_minus.begin();
       j != char_minus.end(); j++)
  {
    int mask=0;
    for (int i=0; i<b.relation().nFinalStateParticles(); i++)
    {
      if (b.relation().finalStateParticle(i).pType() !=
          ptype_gamma)
      {
        if ( (*j).relation().mdstCharged().get_ID() ==
             b.relation().finalStateParticle(i).mdstCharged().get_ID() )
          mask=1;
      }
    }
        std::cout<<"start shape11"<<std::endl;

    Vector4 charged_P = (*j).p();
    charged_P.boost(P_BEAM.boostVector());
    Vector4 tmpptl(charged_P.px(),charged_P.py(),charged_P.pz(),
                   charged_P.e());

    if (mask==1)
    {
    ptlb.push_back(tmpptl);
      //cout << "track " << j <<" is been masked" << endl;
    }
    else
    {

    ptl.push_back(tmpptl);
    }
  }
//end for Particle_List


//one of the B daughters, used to calculate the thrust angle

//  Vector4 ome_P = b.relation().child(0).p();
//  ome_P.boost(P_BEAM.boostVector());

  Vector4 b_P= b.p();
  b_P.boost(P_BEAM.boostVector());

  Vector3 thr_axis = thrust(ptl);
  Vector3 thr_axis_b = thrust(ptlb);

  cos_thr = (thr_axis.dot(thr_axis_b));
  cos_thr = cos_thr / (thr_axis.mag()*thr_axis_b.mag());

Vector3 tmpjet(b_P.px(),b_P.py(),b_P.pz());
//  Vector3 tmpjet(ome_P.px(),ome_P.py(),ome_P.pz());
  spher = Sper(ptl, tmpjet);
        std::cout<<"start shape12"<<std::endl;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
/*void ana_fullreconk::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13],Vector4 &otherB_P, HepPoint3D &bvertex, int & vertexflag, HepPoint3D & overtex,Particle &bta)
{

    const HepPoint3D             &ip     = IpProfile::position();

  Belle_event_Manager& bevt_mgr 	= Belle_event_Manager::get_manager();
  Mdst_charged_Manager &charged_mag    	= Mdst_charged_Manager::get_manager();
  Mdst_vee2_Manager& Vee2Mgr            = Mdst_vee2_Manager::get_manager();
  Gen_hepevt_Manager&  gen_mgr          = Gen_hepevt_Manager::get_manager();
  Mdst_gamma_Manager& gamma_Mgr         = Mdst_gamma_Manager::get_manager();
  Mdst_ecl_Manager &ecl_mag             = Mdst_ecl_Manager::get_manager();
  Mdst_ecl_aux_Manager& mdst_Aux_Mgr    = Mdst_ecl_aux_Manager::get_manager();
  Mdst_ecl_trk_Manager &ecltrk_mag      = Mdst_ecl_trk_Manager::get_manager();

  Mdst_pi0_Manager &pi0_mag             = Mdst_pi0_Manager::get_manager();
  Mdst_event_add_Manager& mevtmgr       = Mdst_event_add_Manager::get_manager();
  Mdst_vee_daughters_Manager& veedmgr   = Mdst_vee_daughters_Manager::get_manager();
  Mdst_klm_mu_ex_Manager& klmmgr        = Mdst_klm_mu_ex_Manager::get_manager();


TagVK ver; // Vertex Reconstructor with kvertexfitter
ver.setdefault(b,bvertex);
//ver.useKsVeto();
ver.dontUseIPforFit();
ver.useTubeforFit();
double E_HER=BeamEnergy::E_HER();
double E_LER=BeamEnergy::E_LER();
double cross_angle=BeamEnergy::Cross_angle();

static Vector4 P_BEAM ( -E_HER*sin(cross_angle), 0.0, E_LER-E_HER*cos(cross_angle),
							  E_HER+E_LER );
for(std::vector<Mdst_charged>::iterator i = charged_mag.begin(); i != charged_mag.end(); i++)
{
	Mdst_charged& ch = *i;
	int flag = 0;
	int nfinal = bta.relation().nFinalStateParticles(); // the number of final particle used to tagB
	for(int k = 0 ; k < nfinal; k++)
	{
		Particle child = bta.relation().finalStateParticle(k);
		if ( child.mdstCharged() )
		{
			if ( child.mdstCharged().get_ID() == ch.get_ID())
			{
				Particle& tmp=*i;
				ver.push_back(&tmp);

			}
		}
	}
	vertexflag=ver.fit();
	if(!vertexflag)
	{
		overtex = ver.vtx();
		//cout<<numver<<":"<<ver.vtx()<<" "<<ver.used_particles().size()<<endl;
	}
}




}*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//for int t--> 0:e 1:mu 2:pi 3:k 4:p 
void ana_fullreconk::GetImpactParameters(const Mdst_charged *charged, double *dr, double *dz, int t)
{
    *dr = 1000.0;
    *dz = 1000.0; 
    if(charged->trk())
    {
        Mdst_trk_fit& trkFit = charged->trk().mhyp(t);
        if( trkFit )
        {
            HepVector a(5,0);
            a[0] = trkFit.helix(0);
            a[1] = trkFit.helix(1);
            a[2] = trkFit.helix(2);
			a[3] = trkFit.helix(3);
            a[4] = trkFit.helix(4);
 
            HepPoint3D pivot(trkFit.pivot(0), trkFit.pivot(1), trkFit.pivot(2));
            Helix ltrk(pivot, a);
            ltrk.pivot(IP);
 
            *dr = ltrk.dr();
            *dz = ltrk.dz();
         }
    }
}
/////////////////////
double ana_fullreconk::pt(Particle p1)
{
    return sqrt(p1.px()*p1.px()+p1.py()*p1.py());
}
double ana_fullreconk::pmag(Particle p1)
{
    return sqrt(p1.px()*p1.px()+p1.py()*p1.py()+p1.pz()*p1.pz());
}
void ana_fullreconk::remove_dup_trk(
        std::vector<Particle> &de_plus, std::vector<Particle> &de_minus \
        ,std::vector<Particle> &dmu_plus, std::vector<Particle> &dmu_minus \
        ,std::vector<Particle> &dpi_plus, std::vector<Particle> &dpi_minus \
        ,std::vector<Particle> &dk_plus, std::vector<Particle> &dk_minus \
        ,std::vector<Particle> &dp_plus, std::vector<Particle> &dp_minus,double DRCUT, double DZCUT)
{   
    std::vector<Particle> pl; 
    pl.clear();
    pl.reserve(de_plus.size() + de_minus.size() + dmu_plus.size() + dmu_minus.size() +\
            dpi_plus.size() + dpi_minus.size() + dk_plus.size() + dk_minus.size() + dp_plus.size() + dp_minus.size());
    pl.insert(pl.end(), de_plus.begin(), de_plus.end());
    pl.insert(pl.end(),de_minus.begin(),de_minus.end());
    pl.insert(pl.end(),dmu_plus.begin(),dmu_plus.end());
    pl.insert(pl.end(),dmu_minus.begin(),dmu_minus.end());
    pl.insert(pl.end(),dpi_plus.begin(),dpi_plus.end());
    pl.insert(pl.end(),dpi_minus.begin(),dpi_minus.end());
    pl.insert(pl.end(),dk_plus.begin(),dk_plus.end());
    pl.insert(pl.end(),dk_minus.begin(),dk_minus.end());
    pl.insert(pl.end(),dp_plus.begin(),dp_plus.end());
    pl.insert(pl.end(),dp_minus.begin(),dp_minus.end());
    std::cout<< "plsize1 = "<<pl.size() << std::endl;
    int plb = pl.size();
    remove_dup_trk(pl, DRCUT, DZCUT);
    int pla = pl.size();
    std::cout<< "plsize2 = "<<pl.size() << std::endl;
    if(abs(pla - plb)!=0)
    std::cout<< " has dup track "<< std::endl;

    de_plus.clear();
    de_minus.clear();
    dmu_plus.clear();
    dmu_minus.clear();
    dpi_plus.clear();
    dpi_minus.clear();
    dk_plus.clear();
    dk_minus.clear();
    dp_plus.clear();
    dp_minus.clear();

    for (std::vector<Particle>::iterator i = pl.begin(); i!=pl.end(); i++)
    {

    std::cout <<"sizebeforeinfun="<< dpi_plus.size()+dpi_minus.size()+de_plus.size()+de_minus.size()+dk_plus.size()+dk_minus.size()+dmu_plus.size()+dmu_minus.size()+dp_plus.size()+dp_minus.size() <<std::endl;

	
        if( (*i).charge() >0)
	{
            if(abs((*i).lund())==11) de_plus.push_back(*i) ; 
            else if (abs((*i).lund())==13) dmu_plus.push_back(*i);
            else if (abs((*i).lund())==211) dpi_plus.push_back(*i);
            else if (abs((*i).lund())==321) dk_plus.push_back(*i);
            else if (abs((*i).lund())==2212) dp_plus.push_back(*i);
        }

        if( (*i).charge() <0)
	{
            if(abs((*i).lund())==11){ de_minus.push_back(*i);}
            else if (abs((*i).lund())==13) dmu_minus.push_back(*i);
            else if (abs((*i).lund())==211) dpi_minus.push_back(*i);
            else if (abs((*i).lund())==321) dk_minus.push_back(*i);
            else if (abs((*i).lund())==2212)dp_minus.push_back(*i);
        }
    std::cout <<"sizeafterinfun="<< dpi_plus.size()+dpi_minus.size()+de_plus.size()+de_minus.size()+dk_plus.size()+dk_minus.size()+dmu_plus.size()+dmu_minus.size()+dp_plus.size()+dp_minus.size() <<std::endl;

    }
}

void ana_fullreconk::remove_dup_trk(std::vector<Particle> &plist, double DRCUT, double DZCUT){
    // remove duplicate tracks
    // for low pt tracks with similar momentum
    if (plist.size() < 1) return;

    for (std::vector<Particle>::iterator i = plist.begin(); i != plist.end(); ++i)
    {
        Vector3 p3i = (*i).p().vect();
        if ( p3i.perp() > THRESHOLD) {
            continue; // only examine pt < 0.25 GeV
        }
        int mhyp1(2);
        if(abs((*i).lund()) == 11)mhyp1=0;
        else if(abs((*i).lund()) == 13)mhyp1=1;
        else if(abs((*i).lund()) == 321)mhyp1=3;
        else if(abs((*i).lund()) == 2212)mhyp1=4;

        for (std::vector<Particle>::iterator j = i+1; j != plist.end(); ++j)
        {
            Vector3 p3j = (*j).p().vect();
            if ( p3j.perp() > THRESHOLD) continue;
            int mhyp2(2);
            if(abs((*j).lund()) == 11)mhyp2=0;
            else if(abs((*j).lund()) == 13)mhyp2=1;
            else if(abs((*j).lund()) == 321)mhyp2=3;
            else if(abs((*j).lund()) == 2212)mhyp2=4;
            if (abs(p3i.mag() - p3j.mag()) > 0.1) continue;
            double cosij = p3i.dot(p3j)/(p3i.mag()*p3j.mag());
            if ( (abs( (*i).charge() - (*j).charge() ) < 0.1 && cosij > 0.95) ||
                    (abs( (*i).charge() + (*j).charge() ) < 0.1 && cosij < -0.95 ) )
		{
                double dr1, dr2, dz1, dz2;
                GetImpactParameters(&(*i).mdstCharged(), &dr1, &dz1, mhyp1);
                GetImpactParameters(&(*j).mdstCharged(), &dr2, &dz2, mhyp2);
                double dist1 = pow(dr1/DRCUT, 2) + pow(dz1/DZCUT, 2);
                double dist2 = pow(dr2/DRCUT, 2) + pow(dz2/DZCUT, 2);
                if (dist1 >= dist2)
		{
                   --(i = plist.erase(i));
                } 
		else 
		{
                   --(j = plist.erase(j));
                }
            }
        }
    }
}


///////////////////
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif	
