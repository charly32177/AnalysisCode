// Hadronic Tag (full reconstruction) + Lambda selection
// Charged and mixed B meson all included
// Suit for any MC sample and Data sample

// Time: 2015/02/11
// Selection: Rank1 B only for each event
// An utility: Can know wherher a B0bar meson 
//K mood//
#include "belle.h"
#include <iostream>
#include <cmath>
#include "mdst/findKs.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "kfitter/kvertexfitter.h"
#include "particle/Particle.h"
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

class ana_fullrecon : public Module 
{
public:
     ana_fullrecon();
     ~ana_fullrecon() {}

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
    void remove_dup_trk(std::vector<Particle> , std::vector<Particle> \
               ,std::vector<Particle> , std::vector<Particle> \
               ,std::vector<Particle> , std::vector<Particle> \
               ,std::vector<Particle> , std::vector<Particle> \
              ,std::vector<Particle> , std::vector<Particle> ,double , double );
    void remove_dup_trk(std::vector<Particle> , double , double );

private:
     BelleTuple* TagB_tpl;
     BelleTuple* B_LX_tpl;
     BelleTuple* b_tpl;

}; //end class definition

extern "C" Module_descr *mdcl_ana_fullrecon() 
{
    ana_fullrecon *module = new ana_fullrecon;
    Module_descr *dscr = new Module_descr( "ana_fullrecon", module );
    return dscr;
}//Module_descr *mdcl_ana_fullrecon()
ana_fullrecon::ana_fullrecon() {  return; }//ana_fullrecon::ana_fullrecon()
void ana_fullrecon::init(int *) { Hamlet::init(); return;}
// Set IP location
HepPoint3D IP( 0, 0, 0 );
HepSymMatrix IPerr( 3, 0 );
void ana_fullrecon::begin_run( BelleEvent* evptr, int *status ) 
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
}//ana_fullrecon::begin_run
void ana_fullrecon::end_run(BelleEvent* evptr, int *status) { (void)evptr; (void)status; return; }

void ana_fullrecon::hist_def() 
{
    extern BelleTupleManager *BASF_Histogram;
    BelleTupleManager *tm = BASF_Histogram;

    TagB_tpl = BASF_Histogram->ntuple( "Rank1 tagB meson","ebeam Tag_RunID Tag_EvtID Tag_ExpID tagid best nbout csbest csnbout tagBmbc tagBde mcinfo nFS decmod ddecmod1 ddecmod2 ddecmod3 ddecmod4 tagBpc PX PY PZ PhEnergy eecl_d BtagE Ek pxk pyk pzk massk sumecl eecl hi1 mo1 gmo1 ggmo1 hindex mmissb pxmissb pymissb pzmissb emissb pmissb mmissCM pxmissCM pymissCM pzmissCM emissCM pmissbCM kpCM kpxCM kpyCM kpzCM mmiss2b mmis2CM", 1 );


// b_tpl = BASF_Histogram->ntuple("B->pinunubar", "eecl mass",2);

    return;
}

/////////////////////////////////////////////////////////
void ana_fullrecon::event(BelleEvent* evptr, int* status)
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
  const double pmax = sqrt( ebeam*ebeam - Lambdamass*Lambdamass ); //used for scaled momentum
  static Vector4 P_BEAM ( -E_HER*sin(cross_angle), 0.0,  E_LER-E_HER*cos(cross_angle), E_HER+E_LER );
  HepLorentzVector boost_vector( -E_HER*sin(cross_angle), 0.0, E_LER-E_HER*cos(cross_angle), E_HER+E_LER );

        std::cout<<"1"<<std::endl;

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

        std::cout<<"2"<<std::endl;






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
//        Ptype ptype_pi_plus_cand("PICAN+");

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
  
  // First Loop over candidates in the Ekpfullrecon panther table to find out B+ and B0
  Ekpfullrecon_Manager &frec_mgr = Ekpfullrecon_Manager::get_manager();
  for( Ekpfullrecon_Manager::iterator frecon_it = frec_mgr.begin(); frecon_it != frec_mgr.end(); frecon_it++ )
  {
         Ekpfullrecon &frecB = *frecon_it;
         // create a B particle
         Particle & Btag = const_cast<Particle &>( brecon.getParticle( (int)frecB.get_ID() ) );
        int    bestCont  = frecB.cont_NBRank();
        if (bestCont!=1) continue;

         // Event ID
  /*       TagB_tpl->column( "ebeam"    , ebeam          );
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


         //Boost to upsilon 4S rest frame
         HepLorentzVector tagB = Btag.p();
         tagB.boost( boost_vector.boostVector() );
         double CM_Btag_p = sqrt( tagB.px()*tagB.px() + tagB.py()*tagB.py() + tagB.pz()*tagB.pz() );
         TagB_tpl->column( "tagBpc"   , CM_Btag_p           );
        
         TagB_tpl->column( "PX"   , tagB.px());
         TagB_tpl->column( "PY"   , tagB.py());
         TagB_tpl->column( "PZ"   , tagB.pz());
          */
//         TagB_tpl->dumpData();
//////////////////////////////////////////////////////////////////////////    
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


        tk_short.clear();
        tLamda.clear();
        tLamdabar.clear();
        tgamma.clear();
        tgamma1.clear();

    // fill Ks list from Mdst_vee2 bank
    for(std::vector<Mdst_vee2>::iterator i = Vee2Mgr.begin(); i != Vee2Mgr.end(); i++)
    {
        std::cout<<"3"<<std::endl;


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


















//////////////////////Satrt Particle identification ////////////////////


        for(std::vector<Mdst_charged>::iterator i = charged_mag.begin();  i != charged_mag.end(); i++)
        {
                Mdst_charged& ch = *i;

                //fullrecon
                int flag = 0;
                int nfinal = Btag.relation().nFinalStateParticles();//Btag paticle
                for(int k=0; k<nfinal; k++) 
		{
                  Particle child = Btag.relation().finalStateParticle(k);
                  if (child.mdstCharged())
			{
                 	 if (child.mdstCharged().get_ID()==ch.get_ID())
		 		{
					 flag = 1;break;
				}
			}
                }
                if (flag) continue;//if flag=1 jump to the top of this for loop
        std::cout<<"4"<<std::endl;




 int flag1 = 0;
            for(std::vector<Particle>::iterator j = tk_short.begin();j != tk_short.end();j++)
            {
                if((*j).child(0).mdstCharged()==(*i)||(*j).child(1).mdstCharged()==(*i)) flag1 =1;
            }
            if (flag1 ==1) continue;
            for(std::vector<Particle>::iterator j = tLamda.begin();j != tLamda.end();j++)
            {
                if((*j).child(0).mdstCharged()==(*i)||(*j).child(1).mdstCharged()==(*i)) flag1 =1;
            }
            if (flag1 ==1) continue;
            for(std::vector<Particle>::iterator j = tLamdabar.begin();j != tLamdabar.end();j++)
            {
                if((*j).child(0).mdstCharged()==(*i)||(*j).child(1).mdstCharged()==(*i)) flag1 =1;
            }
            if (flag1 ==1) continue;
            for(std::vector<Particle>::iterator j = tgamma.begin();j != tgamma.end();j++)
            {
                if((*j).child(0).mdstCharged()==(*i)||(*j).child(1).mdstCharged()==(*i)) {
                flag1 =1;
                }

            }
            if (flag1 ==1) continue;





                //AllTrkList.push_back
       /*       if((*i).charge()>0)
                {
                  Particle all(*i, ptype_pi_plus);
                  AllTrkList.push_back(all);
                }
                else
                {
                  Particle all(*i, ptype_pi_minus);
                  AllTrkList.push_back(all);
                }
       */
                //eid
                eid sel_e(*i);
                float eid_prob = sel_e.prob(3,-1,5);
                //muid
                Muid_mdst muon( *i );
                //int mu_sta = muon.Status();
                //int outcome = muon.Outcome();
                //int mu_level = muon.Level();
                int reject = muon.Prerejection();
                double mu_like = muon.Muon_likelihood();
                //double chi2 = muon.Chi_2();

                //put hadron
        
  //If particle in Btag,we would loss it.

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
    for(std::vector<Particle>::iterator i = p_plus.begin(); i != p_plus.end(); i++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,4);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
            tp_plus.push_back(*i);
    }
    for(std::vector<Particle>::iterator i = p_plus.begin(); i != p_plus.end(); i++)
    {
        double dr,dz;
        const Mdst_charged* ch1 = &(*i).mdstCharged();
        GetImpactParameters(ch1,&dr,&dz,4);
        if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
            tp_plus.push_back(*i);
    }
    remove_dup_trk(te_plus, te_minus,tmu_plus, tmu_minus,tpi_plus, tpi_minus,tk_plus, tk_minus,tp_plus,tp_minus,alldr,alldz);

                std::cout<<"5"<<std::endl;
                   //    float ecl;

	int npip = -1;
	int npim = -1;
	int ntot = -1;
        int npp = -1;
        int npm = -1;
        int nmup = -1;
        int nmum = -1;
        int nkp = -1;
        int nkm = -1;
        int nep = -1;
        int nem = -1;
 	npip = tpi_plus.size();
	npim = tpi_minus.size();
        nep = te_plus.size();
        nem = te_minus.size();
        nkp = tk_plus.size();
        nkm = tk_minus.size();
        nmup = tmu_plus.size();
        nmum = tmu_minus.size();
        npm = tp_plus.size();
        npp = tp_minus.size();

	ntot = npip+npim+nep+nem+nkp+nkm+nmup+nmum+npm+npp ;
        if( ntot!=1)
	continue;
//	if( frecB.tag_id() != -521)
//	continue;
	float eecl = 0;
	float epi = 0 ;

 for(std::vector<Particle>::iterator i =tk_plus.begin(); i != tk_plus.end(); i++)
	
{

	std::cout<<"start tpi"<<std::endl;
 	double sum_ecl = 0;
        for(std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s !=eclaux_mag.end();s++ )
        {
	        std::cout<<"start tpi1"<<std::endl;

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
            if((*i).mdstCharged().trk() == (*s).trk())
            {
            flag=1;
            }
            if(flag) continue;




	
            for(std::vector<Mdst_ecl>::iterator p = ecl_mag.begin(); p!=ecl_mag.end(); p++)
            {
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

		
        std::cout<<"start tpi2"<<std::endl;



            }


        }//ecl_aux end
 
      	

					

	
        std::cout<<"start tpi3"<<std::endl;

					
	
	TagB_tpl->column("pxk" , (*i).px());
        TagB_tpl->column("pyk" , (*i).py());
        TagB_tpl->column("pzk" , (*i).pz());
        TagB_tpl->column("massk" , (*i).p().mag());
	TagB_tpl->column("sumecl",sum_ecl);
//	epi=pi_plus_cand.e();
       	TagB_tpl->column("Ek" , (*i).e());
	eecl = sum_ecl - epi ;
        TagB_tpl->column("eecl" , eecl);

//	TagB_tpl->column("eecl" , eecl );
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

     //Boost to upsilon 4S rest frame
         HepLorentzVector tagB = Btag.p();
         tagB.boost( boost_vector.boostVector() );
         double CM_Btag_p = sqrt( tagB.px()*tagB.px() + tagB.py()*tagB.py() + tagB.pz()*tagB.pz() );
         TagB_tpl->column( "tagBpc"   , CM_Btag_p           );

         TagB_tpl->column( "PX"   , tagB.px());
         TagB_tpl->column( "PY"   , tagB.py());
         TagB_tpl->column( "PZ"   , tagB.pz());
/////////////////////////////////////////////////////////////////////////

	 HepLorentzVector BEAM(E_HER*sin(cross_angle), 0.0, E_HER*cos(cross_angle)-E_LER, E_HER+E_LER);
	 HepLorentzVector btag_p = Btag.p();
	 HepLorentzVector ch1_p((*i).p());
        ch1_p.boost(boost_vector.boostVector() );
        btag_p.boost( boost_vector.boostVector() );
        HepLorentzVector missingchildVec = ( BEAM - btag_p - ch1_p);
       
        TagB_tpl->column("mmissCM",missingchildVec.mag());
        TagB_tpl->column("mmiss2CM",missingchildVec.mag2());
        TagB_tpl->column("pxmissCM",missingchildVec.px());
        TagB_tpl->column("pymissCM",missingchildVec.py());
        TagB_tpl->column("pzmissCM",missingchildVec.pz());
        TagB_tpl->column("emissCM",missingchildVec.e());
        TagB_tpl->column("pmissCM",missingchildVec.vect().mag());





	 HepLorentzVector B_p (-btag_p.px(),-btag_p.py(),-btag_p.pz(),btag_p.e());
        HepLorentzVector miss_b(missingchildVec);
        miss_b.boost(-B_p.boostVector());
        TagB_tpl->column("mmissb",miss_b.mag());
        TagB_tpl->column("mmiss2b",miss_b.mag2());
        TagB_tpl->column("pxmissb",miss_b.px());
        TagB_tpl->column("pymissb",miss_b.py());
        TagB_tpl->column("pzmissb",miss_b.pz());
        TagB_tpl->column("emissb",miss_b.e());
        TagB_tpl->column("pmissb",miss_b.vect().mag());
        TagB_tpl->column("kpCM",ch1_p.vect().mag());
        TagB_tpl->column("kpxCM",ch1_p.px());
        TagB_tpl->column("kpyCM",ch1_p.py());
        TagB_tpl->column("kpzCM",ch1_p.pz());





////////////////////////////////////////////////////////////////////////

        if(MCstatus == 1)
        {
	        std::cout<<"start tpi4"<<std::endl;

            const Mdst_charged ch1 = (*i).mdstCharged();
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
         std::cout<<"before make love"<<std::endl;
	

	           TagB_tpl->column("hindex",hindex);

		}//truth event  end
 std::cout<<"make love0"<<std::endl;

  TagB_tpl->dumpData();
        *status = 1;
std::cout<<"make love1"<<std::endl;
}//pi_plus end 
  }//iterator frecon_it

  
////////////////////////////////////////
 /* // My addition to include photons
  for( std::vector<Mdst_gamma>::iterator it = gamma_Mgr.begin(); it != gamma_Mgr.end(); it++ )
{
    Mdst_gamma& Gamma = *it;
    Particle tmp(*it);

    int flag = 0;
    double deposit_energy = 0;
    for( std::vector<Mdst_ecl_aux>::iterator s = mdst_Aux_Mgr.begin(); s != mdst_Aux_Mgr.end(); s++ )
    {
            Mdst_ecl_aux& ch = *s;

            if( Gamma.ecl().get_ID() == ch.get_ID() )
                flag = 1;

            if( flag == 1 )
            {
                for( std::vector<Mdst_ecl>::iterator p = ecl_mag.begin(); p!=ecl_mag.end(); p++ )
                {
                    if( (*s).get_ID() == (*p).get_ID() )
                    {
                        //Forward endcap calorimeter
                        if( fabs((*p).theta())>0.2096 && fabs((*p).theta())<0.5473 && (*p).energy() > 0.1 )
                            deposit_energy = (*p).energy();
                        //Barral calorimeter
                        if( fabs((*p).theta())>0.5619 && fabs((*p).theta())<2.2466 && (*p).energy() > 0.05 )
                            deposit_energy = (*p).energy();
                        //Backward endcap calorimeter
                        if( fabs((*p).theta())>2.2951 && fabs((*p).theta())<2.7416 && (*p).energy() > 0.1 )
                            deposit_energy = (*p).energy();
                            
                    TagB_tpl->column( "PhEnergy"   , deposit_energy);         
 					TagB_tpl->dumpData();
                    }
                }
            }
    }
}
*/ // My addition to include photons
  
  return;
}// end   void ana_fullrecon::event( BelleEvent* evptr, int* status )

void ana_fullrecon::term() { return; }

//for int t--> 0:e 1:mu 2:pi 3:k 4:p 
void ana_fullrecon::GetImpactParameters(const Mdst_charged *charged, double *dr, double *dz, int t)
{
    *dr = 1000.0;
    *dz = 1000.0;
	std::cout<<"begin if"<< std::endl;
 
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
double ana_fullrecon::pt(Particle p1)
{
    return sqrt(p1.px()*p1.px()+p1.py()*p1.py());
}
double ana_fullrecon::pmag(Particle p1)
{
    return sqrt(p1.px()*p1.px()+p1.py()*p1.py()+p1.pz()*p1.pz());
}
void ana_fullrecon::remove_dup_trk(
        std::vector<Particle> e_plus, std::vector<Particle> e_minus \
        ,std::vector<Particle> mu_plus, std::vector<Particle> mu_minus \
        ,std::vector<Particle> pi_plus, std::vector<Particle> pi_minus \
        ,std::vector<Particle> k_plus, std::vector<Particle> k_minus \
        ,std::vector<Particle> p_plus, std::vector<Particle> p_minus,double DRCUT, double DZCUT)
{   
    std::cout<<"in trk"<<std::endl;
    std::vector<Particle> pl; 
    pl.clear();
    pl.reserve(e_plus.size() + e_minus.size() + mu_plus.size() + mu_minus.size() +\
            pi_plus.size() + pi_minus.size() + k_plus.size() + k_minus.size() + p_plus.size() + p_minus.size());
    pl.insert(pl.end(), e_plus.begin(), e_plus.end());
    pl.insert(pl.end(),e_minus.begin(),e_minus.end());
    pl.insert(pl.end(),mu_plus.begin(),mu_plus.end());
    pl.insert(pl.end(),mu_minus.begin(),mu_minus.end());
    pl.insert(pl.end(),pi_plus.begin(),pi_plus.end());
    pl.insert(pl.end(),pi_minus.begin(),pi_minus.end());
    pl.insert(pl.end(),k_plus.begin(),k_plus.end());
    pl.insert(pl.end(),k_minus.begin(),k_minus.end());
    pl.insert(pl.end(),p_plus.begin(),p_plus.end());
    pl.insert(pl.end(),p_minus.begin(),p_minus.end());
    std::cout<<"b4re"<<std::endl;
    remove_dup_trk(pl, DRCUT, DZCUT);
    std::cout<<"endre"<<std::endl;
    e_plus.clear();
    e_minus.clear();
    mu_plus.clear();
    mu_minus.clear();
    pi_plus.clear();
    pi_minus.clear();
    k_plus.clear();
    k_minus.clear();
    p_plus.clear();
    p_minus.clear();
    for (std::vector<Particle>::iterator i = pl.begin(); i!=pl.end(); i++)
    {
        if( (*i).charge() >0){
            if(abs((*i).lund()==11)) e_plus.push_back(*i);
            else if (abs((*i).lund()==13)) mu_plus.push_back(*i);
            else if (abs((*i).lund())==211) p_plus.push_back(*i);
            else if (abs((*i).lund())==321) k_plus.push_back(*i);
            else if (abs((*i).lund())==2212) p_plus.push_back(*i);
        }
        if( (*i).charge() <0){
            if(abs((*i).lund()==11)) e_minus.push_back(*i);
            else if (abs((*i).lund())==13) mu_minus.push_back(*i);
            else if (abs((*i).lund())==211) p_minus.push_back(*i);
            else if (abs((*i).lund())==321) k_minus.push_back(*i);
            else if (abs((*i).lund())==2212) p_minus.push_back(*i);
        }
    }
}
void ana_fullrecon::remove_dup_trk(std::vector<Particle> plist, double DRCUT, double DZCUT)
{
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
                    (abs( (*i).charge() + (*j).charge() ) < 0.1 && cosij < -0.95 )  )
	    {
                double dr1, dr2, dz1, dz2;
                GetImpactParameters(&(*i).mdstCharged(), &dr1, &dz1, mhyp1);
                GetImpactParameters(&(*j).mdstCharged(), &dr2, &dz2, mhyp2);
                double dist1 = pow(dr1/DRCUT, 2) + pow(dz1/DZCUT, 2);
                double dist2 = pow(dr2/DRCUT, 2) + pow(dz2/DZCUT, 2);
                if (dist1 >= dist2){
                    j = plist.erase(j)-1;
                } else {
                    i = plist.erase(i)-1;
                    break;
                }
            }
        }
    }
}


///////////////////
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif	
