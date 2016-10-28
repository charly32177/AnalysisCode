// First version for B+->kstarnunubar at 2016/08/12
// Hadronic Tag (full reconstruction)
// Charged B meson all included
// B+->h(*)nunubar
// By Bo-Yuan Yang
// Start date: 2016/08/12
// Selection: Rank1 B only for each event
// An utility: Can know wherher a B0bar meson
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
double fullreconeff = 0;
double sigselection = 0;
double belle = 0 ;
double event = 0;
//double fullreconeff = 0;
//double fullreconeff = 0;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
	//utilities
	//belle++;
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
			int toptodown( Gen_hepevt  , Gen_hepevt  , int & ,std::vector<Gen_hepevt> & );

			double deltaZ(Particle, double&, double&, double&, double&);
			double GetCos(Particle a,Particle b);
			double pt(Particle);
			double pmag(Particle);
			double pt(HepLorentzVector&);
			double pmag(HepLorentzVector&);
			double costheta(Particle , Particle);
			double costheta(HepLorentzVector& , HepLorentzVector& );
			double costheta(HepLorentzVector& , Particle);
			double costheta(Particle , HepLorentzVector&);

			void remove_dup_trk(std::vector<Particle>& , std::vector<Particle>& \
					,std::vector<Particle>& , std::vector<Particle>& \
					,std::vector<Particle>& , std::vector<Particle>& \
					,std::vector<Particle>& , std::vector<Particle>& \
					,std::vector<Particle>& , std::vector<Particle>& ,double , double );
			void remove_dup_trk(std::vector<Particle>& , double , double );
		private:
			BelleTuple* TagB_tpl;
			BelleTuple* TagB_tpl2;

			BelleTuple* B_LX_tpl;
			BelleTuple* b_tpl;
			BelleTuple* pi0_tpl;

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

		TagB_tpl = BASF_Histogram->ntuple( "Rank1 tagB meson","nb_trackremaining tnb_trackremaining tnb_pi0remaining RatioToOtherType RaatioToSecondBest sumecl reenergy ebeam Tag_RunID Tag_EvtID Tag_ExpID tagid tagBmbc tagBde BtagE decmod ddecmod1 ddecmod2 ddecmod3 ddecmod4 mcinfo nFS best nbout csbest csnbout tmmisslab tmmiss2lab tpxmisslab tpymisslab tpzmisslab temisslab tpmisslab sblab qsquare qpxlab qpylab qpzlab qelab spmisslab kmalab kpxlab kpylab kpzlab kelab kplab tcosmizl cosksblab cosbcm cosklab tagBpcm tagbPXcm tagbPYcm tagbPZcm angbssmi cosbssmi sigmissp Bsigmas2 sbcm simmis2CM spxmisisC spymissC spzmissC semissC spmissC mkcm kpxCM kpyCM kpzCM keCM kpCM tmmisscm tmmiss2cm tpxmisscm tpymisscm tpzmisscm temisscm tpmisscm smmissb smmiss2b spxmissb spymissb spzmissb semissb spmissb kpb kpxb kpyb kpzb keb tmmissb tmmiss2b tpxmissb tpymissb tpzmissb temissb tpmissb imm2 mmm2 H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son H4son H0soc H1soc H2soc H3soc H4soc H0som H1som H2som H3som H4som et mode hi0 hi10 hi11 mo0 gmo0 ggmo0 mo10 gmo10 ggmo10 mo11 gmo11 ggmo11 gggmo10 gggmo11 hindex bpd1 bpd2 bpd3 bpd4 bpd5 bpd6 bpd7 bpd8 bpd9 nbpd bnd1 bnd2 bnd3 bnd4 bnd5 bnd6 bnd7 bnd8 bnd9 nbnd charged sum_eclinmissdir DistTOOtherBdz kdz kdr vchisq_sig DOF_sig pvalue_sig vchisq_tag DOF_tag pvalue_tag Distsign tnb_gamma nbmulti1 ch0dr ch0dz n_kl_tot n_kl n_ks dpd1 dpd2 dpd3 dpd4 dpd5 dpd6 dpd7 dpd8 dpd9 dnd1 dnd2 dnd3 dnd4 dnd5 dnd6 dnd7 dnd8 dnd9 nb_pi015 kstc0px kstc0py kstc0pz kstc0mass kstc1px kstc1py kstc1pz kstc1mass kstarpx kstarpy kstarpz kstarmass mode check doublecheck hhindex kstcheck kstrank1 topcheck chisqexk krankc1", 1 );
		TagB_tpl2 = BASF_Histogram->ntuple( "Rank1 tagB meson2","nb_trackremaining tnb_trackremaining tnb_pi0remaining RatioToOtherType RaatioToSecondBest sumecl reenergy ebeam Tag_RunID Tag_EvtID Tag_ExpID tagid tagBmbc tagBde BtagE decmod ddecmod1 ddecmod2 ddecmod3 ddecmod4 mcinfo nFS best nbout csbest csnbout tmmisslab tmmiss2lab tpxmisslab tpymisslab tpzmisslab temisslab tpmisslab sblab qsquare qpxlab qpylab qpzlab qelab spmisslab kmalab kpxlab kpylab kpzlab kelab kplab tcosmizl cosksblab cosbcm cosklab tagBpcm tagbPXcm tagbPYcm tagbPZcm angbssmi cosbssmi sigmissp Bsigmas2 sbcm simmis2CM spxmisisC spymissC spzmissC semissC spmissC mkcm kpxCM kpyCM kpzCM keCM kpCM tmmisscm tmmiss2cm tpxmisscm tpymisscm tpzmisscm temisscm tpmisscm smmissb smmiss2b spxmissb spymissb spzmissb semissb spmissb kpb kpxb kpyb kpzb keb tmmissb tmmiss2b tpxmissb tpymissb tpzmissb temissb tpmissb imm2 mmm2 H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son H4son H0soc H1soc H2soc H3soc H4soc H0som H1som H2som H3som H4som et mode hi0 hi10 hi11 mo0 gmo0 ggmo0 mo10 gmo10 ggmo10 mo11 gmo11 ggmo11 gggmo10 gggmo11 hindex bpd1 bpd2 bpd3 bpd4 bpd5 bpd6 bpd7 bpd8 bpd9 nbpd bnd1 bnd2 bnd3 bnd4 bnd5 bnd6 bnd7 bnd8 bnd9 nbnd charged sum_eclinmissdir DistTOOtherBdz kdz kdr vchisq_sig DOF_sig pvalue_sig vchisq_tag DOF_tag pvalue_tag Distsign tnb_gamma nbmulti2 ch0dr ch0dz n_kl_tot n_kl n_ks dpd1 dpd2 dpd3 dpd4 dpd5 dpd6 dpd7 dpd8 dpd9 dnd1 dnd2 dnd3 dnd4 dnd5 dnd6 dnd7 dnd8 dnd9 nb_pi015 kstc0px kstc0py kstc0pz kstc0mass kstc1px kstc1py kstc1pz kstc1mass kstarpx kstarpy kstarpz kstarmass mode check doublecheck hhindex kstcheck kstrank2 topcheck chisqexk krankc2", 11 );


		b_tpl = BASF_Histogram->ntuple( "level1 tagB meson","gggmo10 gggmo11 nb_trackremaining tnb_trackremaining tnb_pi0remaining tnb_pi0remaining nbmulti nb_pi015 n_kl_tot n_kl n_ks kstc0px kstc0py kstc0pz kstc0mass kstc1px kstc1py kstc1pz kstc1mass kstarpx kstarpy kstarpz kstarmass reenergy ebeam Tag_RunID Tag_EvtID Tag_ExpID tagid tagBmbc tagBde BtagE decmod ddecmod1 ddecmod2 ddecmod3 ddecmod4 mcinfo nFS best nbout csbest csnbout mode hi0 hi10 hi11 mo0 gmo0 ggmo0 mo10 gmo10 ggmo10 mo11 gmo11 ggmo11 hindex", 2 );

		pi0_tpl = BASF_Histogram->ntuple("pi0mode","pi0mass pi0hindex a cos hi0 hi1 mo0 mo1 gmo0 gmo1 pi0p pi0e", 13);

		k_sfw::initialize(Fisher_ksfw);

		return;
	}

	/////////////////////////////////////////////////////////
	void ana_fullreconk::event(BelleEvent* evptr, int* status)
	{
		//event++;
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
		std::vector<Particle> tpi_015;
		std::vector<Particle> kstarp1;
		std::vector<Particle> kstarp2;

		std::vector<Particle> kstarm;

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
		Ptype ptype_kstarp("K*+");
		Ptype ptype_kstarm("K*-");

		atc_pid selkpi(3,1,5,3,2);
		atc_pid selppi(3,1,5,4,2);
		atc_pid selkp(3,1,5,3,4);
		atc_pid selpk(3,1,5,4,3);
		atc_pid selmuk(3,1,5,1,3);
		atc_pid selmupi(3,1,5,1,2);

		int MCstatus=0;
		for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin(); i != gen_mgr.end(); i++)
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
			if (bestCont != 1)
				continue;
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
			B_cand2.clear();

			kstarp1.clear();
			kstarp2.clear();
			
			kstarm.clear();
			tk_short.clear();
			tLamda.clear();
			tLamdabar.clear();
			tgamma.clear();
			tgamma1.clear();
			tpi_015.clear();
			fullreconeff++;

			////////////////////////////////////////////////////////////////////////
			//                        K_L identification                          //
			////////////////////////////////////////////////////////////////////////

			int n_kl_tot = 0;
			int n_kl = 0;
			int n_kc_tot = 0;
			int n_kc = 0;

			Mdst_klm_cluster_Manager &klmc_mag     = Mdst_klm_cluster_Manager::get_manager();
			Mdst_klong_Manager &klong_mag     = Mdst_klong_Manager::get_manager();
			for (Mdst_klong_Manager::iterator i = klong_mag.begin(); i != klong_mag.end(); ++i)
			{
				n_kl_tot++;
				if ((*i).klmc())
				{
					Mdst_klm_cluster& Clus = (*i).klmc();
					if (Clus.layers()>=2) n_kl++;
				}
			}
			n_kc_tot = klmc_mag.size();
			for (Mdst_klm_cluster_Manager::iterator i = klmc_mag.begin(); i != klmc_mag.end(); ++i)
			{
				if ((*i).layers()>=2) { n_kc++; }
			}

			////////////////////////////////////////////////////////////////////////
			//                        K_s identification                          //
			////////////////////////////////////////////////////////////////////////

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
				}

			}

			for (std::vector<Particle>::iterator i = k_short.begin(); i != k_short.end(); i++)
			{
	
				int flag = 0;
				int nfinal = Btag.relation().nFinalStateParticles(); // the number of final particle used to tagB
				for(int k = 0 ; k < nfinal; k++)
				{
					Particle child = Btag.relation().finalStateParticle(k);
					if ( child.mdstCharged() )
					{
						for(int k1 = 0 ; k1 < 2; k1++)
						{
							Particle kshortchild =(*i).child(k1);
							if ( child.mdstCharged().get_ID() == kshortchild.mdstCharged().get_ID())
							{
								flag=1;break;
							}
						}
					}
				}
				
				
				if(flag) continue;

				if((*i).mass()>0.482614&&(*i).mass()<0.552614&&(*i).mdstVee2().chisq()<100)
				{
					Particle tmp(*i);
					tk_short.push_back(tmp);
				}

			}
			int n_ks = tk_short.size();

			////////////////////////////////////////////////////////////////////////
			//                        Gamma identification                        //
			////////////////////////////////////////////////////////////////////////

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
				if(flag) continue;

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


				if(flags1||flags2||flags3)
				{
					Particle tmp(*it);
					gamma1.push_back(tmp);

				}
			}

			////////////////////////////////////////////////////////////////////////
			//                    Satrt Particle identification                   //
			////////////////////////////////////////////////////////////////////////

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
				for (std::vector<Particle>::iterator i2 = tk_short.begin(); i2 != tk_short.end(); i2++)
				{
					for(int k2 = 0 ; k2 < 2; k2++)
					{
						Particle childk_s = (*i2).child(k2);
						if ( childk_s.mdstCharged().get_ID() == ch.get_ID())
						{
							flag=1;break;
						}
					}
				}


				if(flag) continue;

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
				if( !reject && mu_like>0.9 )
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
				else if ( eid_prob>0.9 )
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
			//////////////////////////charged particle selection done/////////////////

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
				{
					tp_plus.push_back(*i);
				}
			}
			for(std::vector<Particle>::iterator i = p_minus.begin(); i != p_minus.end(); i++)
			{

				double dr,dz;
				const Mdst_charged* ch1 = &(*i).mdstCharged();
				GetImpactParameters(ch1,&dr,&dz,4);
				if (abs(dr)<alldr&&abs(dz)<alldz&&pt(*i)>allpt)
				{
					tp_minus.push_back(*i);
				}
			}


			remove_dup_trk(te_plus, te_minus,tmu_plus, tmu_minus,tpi_plus, tpi_minus,tk_plus, tk_minus,tp_plus,tp_minus,alldr,alldz);
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

				//				std::cout<<"start pi0"<<std::endl;
				int pi_n;
				double e1,e2,a,p1,p2,p,cos;
				int pi0hindex = -1;
				e1 = (*i).child(0).e();
				e2 = (*i).child(1).e();
				a = abs(e1-e2)/(e1+e2);
				cos= costheta((*i).child(0),(*i).child(1));
				Mdst_pi0 pi0yo = (*i).mdstPi0();
				double pi0mass = pi0yo.mass();
				if(MCstatus == 1)
				{
					const Mdst_gamma ch0 = (*i).child(0).mdstGamma();
					const Mdst_gamma ch1 = (*i).child(1).mdstGamma();
					Gen_hepevt Evtch0 = get_hepevt(ch0);
					pi0_tpl->column("hi0",Evtch0.idhep());
					Gen_hepevt Evtch1 = get_hepevt(ch1);
					pi0_tpl->column("hi1",Evtch1.idhep());
					Gen_hepevt EvtP0;
					Gen_hepevt EvtP1;
					Gen_hepevt EvtGP0;
					Gen_hepevt EvtGP1;
					if(Evtch0.mo(0))
					{
						EvtP0 = gen_mgr[Evtch0.mo(0)-1];
						pi0_tpl->column("mo0",EvtP0.idhep());

						if(EvtP0.mo(0))
						{
							EvtGP0 = gen_mgr[EvtP0.mo(0)-1];
							pi0_tpl->column("gmo0",EvtGP0.idhep());
						}

					}
					if(Evtch1.mo(0))
					{
						EvtP1 = gen_mgr[Evtch1.mo(0)-1];
						pi0_tpl->column("mo1",EvtP1.idhep());

						if(EvtP1.mo(0))
						{
							EvtGP1 = gen_mgr[EvtP1.mo(0)-1];
							pi0_tpl->column("gmo1",EvtGP1.idhep());
						}

					}
					if( Evtch0.idhep()==22&&Evtch1.idhep()==22&& EvtP1.idhep() ==111&&EvtP0.idhep()==111)
					{
						pi0hindex = 1;
					}
					else
					{
						pi0hindex = 0;
					}
				}
				//				pi0_tpl->column("pi0hindex",pi0hindex);
				//				pi0_tpl->column("pi0mass",pi0mass);
				//				pi0_tpl->column("a",a);
				//				pi0_tpl->column("cos",cos);
				//				pi0_tpl->column("pi0p",pmag(*i));
				//				pi0_tpl->column("pi0e",(*i).e());

				//				pi0_tpl->dumpData();

				/////////////////////////////////
				if (e1>0.05&&e2>0.05&&pi0mass>0.1178&&pi0mass<0.1502&&a<0.65)
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
			/*
					
					int flag = 0;
					int nfinal = Btag.relation().nFinalStateParticles(); // the number of final particle used to tagB
					for(int k = 0 ; k < nfinal; k++)
					{
						Particle child = Btag.relation().finalStateParticle(k);
						if ( child.mdstGamma() )
						{
							for(int k1 = 0 ; k1 < 2; k1++)
							{
								Particle pi0child =(*i).child(k1);
								if ( child.mdstGamma().get_ID() == pi0child.mdstGamma().get_ID())
								{
									flag=1;break;
								}
							}
						}
					}
					
					
					if(flag) continue;
*/
					
					
					if(flags1 || flags2 || flags3)
					{
						Particle tmp(*i);
						tpi_0.push_back(tmp);
						if (pmag(*i) > 1.5 )
						{
							Particle tmp15(*i);
							tpi_015.push_back(tmp15);
						}
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
			int npi015 = -1;
			npi015 = tpi_015.size();

			//	  if (npi0 != 0)
			//		  continue;
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

			n_tks = tk_short.size();
			n_tlamda = tLamda.size();
			n_tlamdabar = tLamdabar.size();
			n_tgamma = tgamma.size();
			n_tgamma1 = gamma1.size();

			float eecl = 0;
			float ek = 0 ;
			///////////////////////////////////k/////////////////////////////////////////////

			if (frecB.tag_id() == -521)
			{
				combination(kstarp1,ptype_kstarp,tk_plus,tpi_0);
				combination(kstarp2,ptype_kstarp,tk_short,tpi_plus);

			}
			///////////////////////////
			else if (frecB.tag_id() == 521)
			{
				combination(kstarp1,ptype_kstarm,tk_minus,tpi_0);
				combination(kstarp2,ptype_kstarm,tk_short,tpi_minus);
			}
/*
			for(std::vector<Particle>::iterator i =kstarp.begin(); i !=kstarp.end(); i++)
			{
				for(std::vector<Particle>::iterator j=i+1; j != kstarp.end(); j++)
				{
						double massd1;
						double massd2;
							massd1 = fabs( 0.8916 - (*i).p().mag() );
							massd2 = fabs( 0.8916 - (*j).p().mag() );
						if (massd1 > massd2)
						{
							--(i = kstarp.erase(i));
							break;
						}
						else if (massd1 < massd2)
						{
							--(j = kstarp.erase(j));
						}
						else
							--(j = kstarp.erase(j));
					}
				}
*/
			for(std::vector<Particle>::iterator i5 =kstarp1.begin(); i5 != kstarp1.end(); i5++)
			{
				if(frecB.tag_id() == -521)
				{
					Particle Bsig((*i5).momentum().p(), ptype_B_plus_cand);
					Bsig.relation().append((*i5));
					B_cand1.push_back(Bsig);
					
				}
				else
				{
					Particle Bsig((*i5).momentum().p(), ptype_B_minus_cand);
					Bsig.relation().append((*i5));
					B_cand1.push_back(Bsig);
				}
			}
			for(std::vector<Particle>::iterator i5 =kstarp2.begin(); i5 != kstarp2.end(); i5++)
			{
				if(frecB.tag_id() == -521)
				{
					Particle Bsig((*i5).momentum().p(), ptype_B_plus_cand);
					Bsig.relation().append((*i5));
					B_cand2.push_back(Bsig);
					
				}
				else
				{
					Particle Bsig((*i5).momentum().p(), ptype_B_minus_cand);
					Bsig.relation().append((*i5));
					B_cand2.push_back(Bsig);
				}
			}

//			double nbmulti = B_cand1.size();
//			for(std::vector<Particle>::iterator i =kstarp.begin(); i != kstarp.end(); i++)
//			{
//				std::cout<< "1" << std::endl;
//
//				tnpip = tpi_plus.size();
//				tnpim = tpi_minus.size();
//				tnep = te_plus.size();
//				tnem = te_minus.size();
//				tnkp = tk_plus.size();
//				tnkm = tk_minus.size();
//				tnmup = tmu_plus.size();
//				tnmum = tmu_minus.size();
//				tnpm = tp_plus.size();
//				tnpp = tp_minus.size();
//				npip = pi_plus.size();
//				npim = pi_minus.size();
//				nep = e_plus.size();
//				nem = e_minus.size();
//				nkp = k_plus.size();
//				nkm = k_minus.size();
//				nmup = mu_plus.size();
//				nmum = mu_minus.size();
//				npm = p_plus.size();
//				npp = p_minus.size();
//
//				tnremaining = tnpip+tnpim+tnep+tnem+tnkp+tnkm+tnmup+tnmum+tnpm+tnpp ;
//				nremaining = npip+npim+nep+nem+nkp+nkm+nmup+nmum+npm+npp ;
//				std::cout<< "2" << std::endl;
//				b_tpl->column( "nb_trackremaining"    , nremaining         );//new
//				b_tpl->column( "tnb_trackremaining"    , tnremaining         );//new
//				b_tpl->column( "tnb_pi0remaining"    , npi0         );//new
//				b_tpl->column( "tnb_gamma"    , n_tgamma1         );//new
//				b_tpl->column( "nbmulti"    , nbmulti         );//new
//				b_tpl->column( "nb_pi015"    , npi015         );//new
//
//				b_tpl->column( "n_kl_tot"    , n_kl_tot         );//new
//				b_tpl->column( "n_kl"    , n_kl         );//new
//				b_tpl->column( "n_ks"    , n_ks         );
//
//				b_tpl->column("kstc0px" , (*i).child(0).px());
//				b_tpl->column("kstc0py" , (*i).child(0).py());
//				b_tpl->column("kstc0pz" , (*i).child(0).pz());
//				b_tpl->column("kstc0mass" , (*i).child(0).p().mag());
//
//				b_tpl->column("kstc1px" , (*i).child(1).px());
//				b_tpl->column("kstc1py" , (*i).child(1).py());
//				b_tpl->column("kstc1pz" , (*i).child(1).pz());
//				b_tpl->column("kstc1mass" , (*i).child(1).p().mag());
//
//				b_tpl->column("kstarpx" , (*i).px());
//				b_tpl->column("kstarpy" , (*i).py());
//				b_tpl->column("kstarpz" , (*i).pz());
//				b_tpl->column("kstarmass" , (*i).p().mag());
//
//				//////
//				//		  double reenergy = P_BEAM.e() - totalecl;
//				//		  b_tpl->column( "reenergy"    , reenergy          );
//
//				////////////////////////Btag///////////////////////////////////////////////////////////
//				// Event ID
//				b_tpl->column( "ebeam"    , ebeam          );
//				b_tpl->column( "Tag_RunID", Run            );
//				b_tpl->column( "Tag_EvtID", Evt            );
//				b_tpl->column( "Tag_ExpID", Exp            );
//				// TagB idhep
//				b_tpl->column( "tagid"    , frecB.tag_id() );
//				// Kinetic variable
//				b_tpl->column( "tagBmbc"  , frecB.Mbc()    );
//				b_tpl->column( "tagBde"   , frecB.DeltaE() );
//				float BtagE = frecB.DeltaE() + ebeam ;
//				b_tpl->column( "BtagE"  , BtagE   );
//				// unique code for reconstructed D decay channels
//				// for B -> D X decay modes only 1 entry is filled
//				// for B -> D* ( -> D X) Y decay modes only 1 (D*) and 2 (D) entries are filled
//				// for B -> D* ( -> D1 X) D2 decay modes only 1 (D*), 2 (D1) and 3 (D2) entries are filled
//				// for B -> D*1 ( -> D1 X) D*2 ( -> D2 Y) decay modes all four entries are filled: D*1, D1, D*2, D2
//				b_tpl->column( "decmod"   , frecB.decay()       );
//				b_tpl->column( "ddecmod1" , frecB.Ddec(0)       );
//				b_tpl->column( "ddecmod2" , frecB.Ddec(1)       );
//				b_tpl->column( "ddecmod3" , frecB.Ddec(2)       );
//				b_tpl->column( "ddecmod4" , frecB.Ddec(3)       );
//				b_tpl->column( "mcinfo"   , frecB.MCinfo()      );
//				b_tpl->column( "nFS"      , frecB.nFS()         );
//				// default selection (without continuum suppression applied)
//				b_tpl->column( "best"     , frecB.NBRank()      );
//				b_tpl->column( "nbout"    , frecB.NBout()       );
//				// selection with continuum suppression
//				b_tpl->column( "csbest"   , frecB.cont_NBRank() );
//				b_tpl->column( "csnbout"  , frecB.cont_NBout()  );
//				///////////////////////////////////////////////////////////////////////////////////////
//				//                      MC TRuth(I don't know it work or not)                        //
//				///////////////////////////////////////////////////////////////////////////////////////
//				int mode = -1;
//				std::cout<< "3" << std::endl;
//				if ((*i).child(0).mdstCharged())
//				{
//					mode = 1;
//				}
//				else
//				{
//					mode = 2;
//				}
//				std::cout<< "mode = " << mode <<std::endl;
//				b_tpl->column( "mode"  , mode  );
//
//				if(MCstatus == 1)
//				{
//					if (mode == 1) // k+ pi0
//					{
//						const Mdst_charged ch0 = (*i).child(0).mdstCharged();
//						const Mdst_gamma ch10 = (*i).child(1).child(0).mdstGamma();
//						const Mdst_gamma ch11 = (*i).child(1).child(1).mdstGamma();
//
//						Gen_hepevt Evtch0 = get_hepevt(ch0);
//						Gen_hepevt Evtch10 = get_hepevt(ch10);
//						Gen_hepevt Evtch11 = get_hepevt(ch11);
//						
//						Gen_hepevt tEvtch0 ;
//						Gen_hepevt tEvtch10;
//						Gen_hepevt tEvtch11;
//						std::cout<< "debug0 ch0id = " << Evtch0.idhep() <<std::endl;
//						std::cout<< "debug0 ch10id = " << Evtch10.idhep() <<std::endl;
//						std::cout<< "debug0 ch11id = " << Evtch11.idhep() <<std::endl;
//						
//						b_tpl->column("hi0",Evtch0.idhep());
//						b_tpl->column("hi10",Evtch10.idhep());
//						b_tpl->column("hi11",Evtch11.idhep());
//
//						Gen_hepevt EvtP0;
//						Gen_hepevt EvtGP0;
//						Gen_hepevt EvtGGP0;
//
//						Gen_hepevt EvtP10;
//						Gen_hepevt EvtGP10;
//						Gen_hepevt EvtGGP10;
//
//						Gen_hepevt EvtP11;
//						Gen_hepevt EvtGP11;
//						Gen_hepevt EvtGGP11;
//
//						if(Evtch0.mo(0))
//						{
//							EvtP0 = gen_mgr[Evtch0.mo(0)-1];
//							b_tpl->column("mo0",EvtP0.idhep());
//
//							if(EvtP0.mo(0))
//							{
//								EvtGP0 = gen_mgr[EvtP0.mo(0)-1];
//								b_tpl->column("gmo0",EvtGP0.idhep());
//
//
//								if(EvtGP0.mo(0))
//								{
//									EvtGGP0 = gen_mgr[EvtGP0.mo(0)-1];
//									b_tpl->column("ggmo0",EvtGGP0.idhep());
//								}
//							}
//						}
//						if(Evtch10.mo(0))
//						{
//							EvtP10 = gen_mgr[Evtch10.mo(0)-1];
//							b_tpl->column("mo10",EvtP10.idhep());
//
//							if(EvtP10.mo(0))
//							{
//								EvtGP10 = gen_mgr[EvtP10.mo(0)-1];
//								b_tpl->column("gmo10",EvtGP10.idhep());
//
//
//								if(EvtGP10.mo(0))
//								{
//									EvtGGP10 = gen_mgr[EvtGP10.mo(0)-1];
//									b_tpl->column("ggmo10",EvtGGP10.idhep());
//								}
//							}
//						}
//						if(Evtch11.mo(0))
//						{
//							EvtP11 = gen_mgr[Evtch11.mo(0)-1];
//							b_tpl->column("mo11",EvtP11.idhep());
//
//							if(EvtP11.mo(0))
//							{
//								EvtGP11 = gen_mgr[EvtP11.mo(0)-1];
//								b_tpl->column("gmo11",EvtGP11.idhep());
//
//
//								if(EvtGP11.mo(0))
//								{
//									EvtGGP11 = gen_mgr[EvtGP11.mo(0)-1];
//									b_tpl->column("ggmo11",EvtGGP11.idhep());
//								}
//							}
//						}
//						int hindex = -1 ;
//
//						if( Evtch0.idhep()== 321 && Evtch10.idhep()==22 && Evtch11.idhep()==22 && EvtP0.idhep()==323 && EvtP10.idhep()==111 && EvtP11.idhep()==111 && EvtGP0.idhep()==521 && EvtGP10.idhep()==323 && EvtGP11.idhep()==323 && EvtGGP10.idhep()==521 && EvtGGP11.idhep()==521 )
//						{
//							hindex = 1;
//						}
//						else if( Evtch0.idhep()== -321 && Evtch10.idhep()==22 && Evtch11.idhep()==22 && EvtP0.idhep()==-323 && EvtP10.idhep()==111 && EvtP11.idhep()==111 && EvtGP0.idhep()==-521 && EvtGP10.idhep()==-323 && EvtGP11.idhep()==-323 && EvtGGP10.idhep()==-521 && EvtGGP11.idhep()==-521 )
//						{
//							hindex = 2;
//						}
//						else
//						{
//							hindex = 0;
//						}
//
//						b_tpl->column("hindex",hindex);
//
//					}
//					else if (mode == 2) // ks pi+
//					{
//						std::cout<< "4" << std::endl;
//						std::cout<< "id0 =" << (*i).child(1).lund() <<std::endl;
//						const Mdst_charged ch0 = (*i).child(1).mdstCharged();
//
//						std::cout<< "id1 =" << (*i).child(0).child(0).lund() <<std::endl;
//						std::cout<< "id2 =" << (*i).child(0).child(1).lund()<<std::endl;
//
//						const Mdst_charged ch10 = (*i).child(0).child(0).mdstCharged();
//						const Mdst_charged ch11 = (*i).child(0).child(1).mdstCharged();
//						std::cout<< "5" << std::endl;
//						Gen_hepevt Evtch0 = get_hepevt(ch0);
//						Gen_hepevt Evtch10 = get_hepevt(ch10);
//						Gen_hepevt Evtch11 = get_hepevt(ch11);
//						
//						Gen_hepevt tEvtch0 ;
//						Gen_hepevt tEvtch10;
//						Gen_hepevt tEvtch11;
//						
//						std::cout<< "6" << std::endl;
//						b_tpl->column("hi0",Evtch0.idhep());
//						b_tpl->column("hi10",Evtch10.idhep());
//						b_tpl->column("hi11",Evtch11.idhep());
//
//						Gen_hepevt EvtP0;
//						Gen_hepevt EvtGP0;
//						Gen_hepevt EvtGGP0;
//
//						Gen_hepevt EvtP10;
//						Gen_hepevt EvtGP10;
//						Gen_hepevt EvtGGP10;
//						Gen_hepevt EvtGGGP10;
//
//						Gen_hepevt EvtP11;
//						Gen_hepevt EvtGP11;
//						Gen_hepevt EvtGGP11;
//						Gen_hepevt EvtGGGP11;
//						std::cout<< "7" << std::endl;
//						if(Evtch0.mo(0))
//						{
//							EvtP0 = gen_mgr[Evtch0.mo(0)-1];
//							b_tpl->column("mo0",EvtP0.idhep());
//
//							if(EvtP0.mo(0))
//							{
//								EvtGP0 = gen_mgr[EvtP0.mo(0)-1];
//								b_tpl->column("gmo0",EvtGP0.idhep());
//
//
//								if(EvtGP0.mo(0))
//								{
//									EvtGGP0 = gen_mgr[EvtGP0.mo(0)-1];
//									b_tpl->column("ggmo0",EvtGGP0.idhep());
//									//								if(EvtGGP0.mo(0))
//									//								{
//									//									EvtGGGP0 = gen_mgr[EvtGGP0.mo(0)-1];
//									//									b_tpl->column("gggmo0",EvtGGGP0.idhep());
//									//								}
//								}
//							}
//						}
//						std::cout<< "8" << std::endl;
//						if(Evtch10.mo(0))
//						{
//							EvtP10 = gen_mgr[Evtch10.mo(0)-1];
//							b_tpl->column("mo10",EvtP10.idhep());
//
//							if(EvtP10.mo(0))
//							{
//								EvtGP10 = gen_mgr[EvtP10.mo(0)-1];
//								b_tpl->column("gmo10",EvtGP10.idhep());
//
//
//								if(EvtGP10.mo(0))
//								{
//									EvtGGP10 = gen_mgr[EvtGP10.mo(0)-1];
//									b_tpl->column("ggmo10",EvtGGP10.idhep());
//									if(EvtGGP10.mo(0))
//									{
//										EvtGGGP10 = gen_mgr[EvtGGP10.mo(0)-1];
//										b_tpl->column("gggmo10",EvtGGGP10.idhep());
//									}
//								}
//							}
//						}
//						std::cout<< "9" << std::endl;
//						if(Evtch11.mo(0))
//						{
//							EvtP11 = gen_mgr[Evtch11.mo(0)-1];
//							b_tpl->column("mo11",EvtP11.idhep());
//
//							if(EvtP11.mo(0))
//							{
//								EvtGP11 = gen_mgr[EvtP11.mo(0)-1];
//								b_tpl->column("gmo11",EvtGP11.idhep());
//
//
//								if(EvtGP11.mo(0))
//								{
//									EvtGGP11 = gen_mgr[EvtGP11.mo(0)-1];
//									b_tpl->column("ggmo11",EvtGGP11.idhep());
//									if(EvtGGP11.mo(0))
//									{
//										EvtGGGP11 = gen_mgr[EvtGGP11.mo(0)-1];
//										b_tpl->column("gggmo11",EvtGGGP11.idhep());
//
//									}
//								}
//							}
//						}
//						int hindex = -1 ;
//						std::cout<< "10" << std::endl;
//						if( Evtch0.idhep()== 211 && Evtch10.idhep()==211 && Evtch11.idhep()==-211 && EvtP0.idhep()==323 && EvtP10.idhep()==310 && EvtP11.idhep()==310 && EvtGP0.idhep()==521 && abs(EvtGP10.idhep())==311 && abs(EvtGP11.idhep())==311 && EvtGGP10.idhep()==323 && EvtGGP11.idhep()==323 && EvtGGGP10.idhep()==521 && EvtGGGP11.idhep()==521  )
//						{
//							hindex = 1;
//						}
//						else if( Evtch0.idhep()== -211 && Evtch10.idhep()==211 && Evtch11.idhep()==-211 && EvtP0.idhep()==-323 && EvtP10.idhep()==310 && EvtP11.idhep()==310 && EvtGP0.idhep()==-521 && abs(EvtGP10.idhep())==311 && abs(EvtGP11.idhep())==311 && EvtGGP10.idhep()==-323 && EvtGGP11.idhep()==-323 && EvtGGGP10.idhep()==-521 && EvtGGGP11.idhep()==-521 )
//						{
//							hindex = 2;
//						}
//						else
//						{
//							hindex = 0;
//						}
//
//						b_tpl->column("hindex",hindex);
//
//					}
//					std::cout<< "11" << std::endl;
//				}//truth event  end
//
//				/////////////////////////////////trace form top to down//////////////
//				/*
//					 int nbpd=0,nbnd=0;
//					 for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin();i != gen_mgr.end(); i++)
//					 {
//					 if ((*i).idhep()==521||(*i).idhep()== 511)//511==B0
//					 {
//					 int da1=(*i).da(0),da2=(*i).da(1);
//
//					 if ((*i).idhep() == 521)
//					 b_tpl->column("charged",1);
//					 else
//					 b_tpl->column("charged",0);
//
//					 nbpd=(da2-da1+1);
//					 b_tpl->column("nbpd",nbpd);
//
//					 for(int start=0;start<(da2-da1+1);start++)
//					 {
//					 Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
//					 char bpdnumber[32];
//					 sprintf (bpdnumber,"%s%d","bpd",start+1);
//					 b_tpl->column(bpdnumber,Evda.idhep());
//					 }
//
//					 Gen_hepevt Evda1=gen_mgr[(*i).da(0)-1];
//					 int dada1 = Evda1.da(0),dada2 = Evda1.da(1);
//
//
//					 if( Evda1.idhep()== -421 || Evda1.idhep()== 421)
//					 {
//
//					 for(int start=0;start<(dada2-dada1+1);start++)
//					 {
//					 Gen_hepevt Evdadada1=gen_mgr[Evda1.da(0)-1+start];
//					 char dndnumber[32];
//					 sprintf (dndnumber,"%s%d","dpd",start+1);//bnd1,bnd2,bnd3....
//					 b_tpl->column(dndnumber,Evdadada1.idhep());
//					 }
//					 }
//
//					 if( Evda1.idhep()== -423 || Evda1.idhep()== 423)
//					 {
//					 Gen_hepevt Evdada1=gen_mgr[Evda1.da(0)-1];
//					 int dadada1 = Evdada1.da(0),dadada2 = Evdada1.da(1);
//					 for(int start=0;start<(dadada2-dadada1+1);start++)
//					 {
//					 Gen_hepevt Evdadada1=gen_mgr[Evdada1.da(0)-1+start];
//					 char dndnumber[32];
//					 sprintf (dndnumber,"%s%d","dpd",start+1);//bnd1,bnd2,bnd3....
//					 b_tpl->column(dndnumber,Evdadada1.idhep());
//					 }
//					 }
//
//					 }
//					 else if((*i).idhep()== -521 || (*i).idhep()== -511)
//					 {
//					 int da1=(*i).da(0),da2=(*i).da(1);
//
//					 if ((*i).idhep() == -521)
//					 b_tpl->column("charged",1);
//					 else b_tpl->column("charged",0);
//					 nbnd=(da2-da1+1);
//					 b_tpl->column("nbnd",nbnd);
//					 for(int start=0;start<(da2-da1+1);start++)
//					 {
//					 Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
//					 char bndnumber[32];
//					 sprintf (bndnumber,"%s%d","bnd",start+1);//bnd1,bnd2,bnd3....
//					 b_tpl->column(bndnumber,Evda.idhep());
//					 }
//					 Gen_hepevt Evda1=gen_mgr[(*i).da(0)-1];
//				int dada1 = Evda1.da(0),dada2 = Evda1.da(1);
//
//				if( Evda1.idhep()== -421 || Evda1.idhep()== 421)
//				{
//					for(int start=0;start<(dada2-dada1+1);start++)
//					{
//						Gen_hepevt Evdadada1=gen_mgr[Evda1.da(0)-1+start];
//						char dndnumber[32];
//						sprintf (dndnumber,"%s%d","dnd",start+1);//bnd1,bnd2,bnd3....
//						b_tpl->column(dndnumber,Evdadada1.idhep());
//					}
//				}
//
//				if( Evda1.idhep()== -423 || Evda1.idhep()== 423)
//				{
//					Gen_hepevt Evdada1=gen_mgr[Evda1.da(0)-1];
//					int dadada1 = Evdada1.da(0),dadada2 = Evdada1.da(1);
//
//					for(int start=0;start<(dadada2-dadada1+1);start++)
//					{
//						Gen_hepevt Evdadada1=gen_mgr[Evdada1.da(0)-1+start];
//						char dndnumber[32];
//						sprintf (dndnumber,"%s%d","dnd",start+1);//bnd1,bnd2,bnd3....
//						b_tpl->column(dndnumber,Evdadada1.idhep());
//					}
//				}
//
//			}
//			}
//			*/
//
//
//				b_tpl->dumpData();
//				*status = 1;
//			}

			///////////////////////////////////k/////////////////////////////////////////////
			
			double nbmulti1 = B_cand1.size();
			double nbmulti2 = B_cand2.size();
			int kstrank1 = -999;
			int kstrank2 = -999;
			int krankc1 = -999;
			int krankc2 = -999;
			for(std::vector<Particle>::iterator i =B_cand1.begin(); i != B_cand1.end(); i++)
			{
				for(std::vector<Particle>::iterator j1=B_cand1.begin(); j1 != B_cand1.end(); j1++)
				{
					if (i == j1)
						continue;
					
					double massd1;
					double massd2;
					massd1 = fabs( 0.8916 - (*i).p().mag() );
					massd2 = fabs( 0.8916 - (*j1).p().mag() );
					if (massd1 > massd2)
					{
						kstrank1 = 0;
						break;
					}
					else if (massd1 < massd2)
					{
						kstrank1 = 1;
					}
					else
						kstrank1 = -1;
				}
				
				TagB_tpl->column( "kstrank1"  , kstrank1  );

				
				double mode;
				if ((*i).child(0).child(0).mdstCharged())
				{
					mode = 1;
				}
				else
				{
					mode = 2;
				}
				TagB_tpl->column( "mode"  , mode  );

				std::cout<< "mode = " << mode <<std::endl;
				Mdst_gamma Mdst_10 =(*i).child(0).child(1).child(0).mdstGamma();
				Mdst_gamma Mdst_11 =(*i).child(0).child(1).child(1).mdstGamma();
				Mdst_charged Mdst_0 =(*i).child(0).child(0).mdstCharged();
				
				ExKFitterParticle KF_10(Mdst_10);
				ExKFitterParticle KF_11(Mdst_11);
				ExKFitterParticle KF_0(Mdst_0,3);
				
				HepPoint3D kstar_init;
				
				kstar_init.setX(IP.x() + (*i).child(0).p().px()/(*i).child(0).p().rho());
				kstar_init.setY(IP.y() + (*i).child(0).p().py()/(*i).child(0).p().rho());
				kstar_init.setZ(IP.z() + (*i).child(0).p().pz()/(*i).child(0).p().rho());
				ExKFitterVertex kstar_Vertex(kstar_init);
				ExKFitterVertex B_Vertex(IP,IPerr);
				
				ExKFitterMass Pi0_Mass1(0.134976);
				ExKFitterMass Pi0_Mass2(0.134976);
				ExKFitterMass kstar_Mass(0.8916);
				
				ExKFitterParticle kstar;
				kstar.LinkParticle(&KF_10);
				kstar.LinkParticle(&KF_11);
				kstar.LinkParticle(&KF_0);
				kstar.LinkVertex(&kstar_Vertex);
				
				ExKFitterConstrain con1;
				con1.SetMassConstrain();
				con1.LinkParticle(&KF_10);
				con1.LinkParticle(&KF_11);
				con1.LinkVertex(&kstar_Vertex);
				con1.LinkMass(&Pi0_Mass1);
				
				ExKFitterConstrain con2;
				con2.SetMassConstrain();
				con2.LinkParticle(&KF_10);
				con2.LinkParticle(&KF_11);
				con2.LinkParticle(&KF_0);
				con2.LinkVertex(&kstar_Vertex);
				con2.LinkMass(&kstar_Mass);
				
				ExKFitterConstrain con3;
				con3.SetVertexConstrain();
				con3.LinkParticle(&kstar);
				con3.LinkVertex(&B_Vertex);
				
				ExKFitter Core;
				Core.LinkConstrain(&con1);
				Core.LinkConstrain(&con2);
				Core.LinkConstrain(&con3);
				
				int ret = Core.Minimize();
				float chisqExK = Core.Chisq();
				float dof_exk = Core.N_DegreeOfFreedom();
				HepPoint3D dvertex = kstar_Vertex.Vertex();
				HepPoint3D bvertex1 = B_Vertex.Vertex();
				
				if(ret==0)
				{
					kstar.Update();
				}
				
				TagB_tpl->column("chisqexk",chisqExK/dof_exk);
				
				for(std::vector<Particle>::iterator j1=B_cand1.begin(); j1 != B_cand1.end(); j1++)
				{

					Mdst_gamma jMdst_10 =(*j1).child(0).child(1).child(0).mdstGamma();
					Mdst_gamma jMdst_11 =(*j1).child(0).child(1).child(1).mdstGamma();
					Mdst_charged jMdst_0 =(*j1).child(0).child(0).mdstCharged();
					
					ExKFitterParticle jKF_10(jMdst_10);
					ExKFitterParticle jKF_11(jMdst_11);
					ExKFitterParticle jKF_0(jMdst_0,3);
					
					HepPoint3D jkstar_init;
					
					jkstar_init.setX(IP.x() + (*j1).child(0).p().px()/(*j1).child(0).p().rho());
					jkstar_init.setY(IP.y() + (*j1).child(0).p().py()/(*j1).child(0).p().rho());
					jkstar_init.setZ(IP.z() + (*j1).child(0).p().pz()/(*j1).child(0).p().rho());
					ExKFitterVertex jkstar_Vertex(jkstar_init);
					ExKFitterVertex jB_Vertex(IP,IPerr);
					
					ExKFitterMass jPi0_Mass1(0.134976);
					ExKFitterMass jPi0_Mass2(0.134976);
					ExKFitterMass jkstar_Mass(0.8916);
					
					ExKFitterParticle jkstar;
					jkstar.LinkParticle(&jKF_10);
					jkstar.LinkParticle(&jKF_11);
					jkstar.LinkParticle(&jKF_0);
					jkstar.LinkVertex(&jkstar_Vertex);
					
					ExKFitterConstrain jcon1;
					jcon1.SetMassConstrain();
					jcon1.LinkParticle(&jKF_10);
					jcon1.LinkParticle(&jKF_11);
					jcon1.LinkVertex(&jkstar_Vertex);
					jcon1.LinkMass(&jPi0_Mass1);
					
					ExKFitterConstrain jcon2;
					jcon2.SetMassConstrain();
					jcon2.LinkParticle(&jKF_10);
					jcon2.LinkParticle(&jKF_11);
					jcon2.LinkParticle(&jKF_0);
					jcon2.LinkVertex(&jkstar_Vertex);
					jcon2.LinkMass(&jkstar_Mass);
					
					ExKFitterConstrain jcon3;
					jcon3.SetVertexConstrain();
					jcon3.LinkParticle(&jkstar);
					jcon3.LinkVertex(&jB_Vertex);
					
					ExKFitter jCore;
					jCore.LinkConstrain(&jcon1);
					jCore.LinkConstrain(&jcon2);
					jCore.LinkConstrain(&jcon3);
					
					int jret = jCore.Minimize();
					float jchisqExK = jCore.Chisq();
					float jdof_exk = jCore.N_DegreeOfFreedom();
					HepPoint3D jdvertex = jkstar_Vertex.Vertex();
					HepPoint3D jbvertex1 = jB_Vertex.Vertex();
					if(jret==0)
					{
						jkstar.Update();
					}
					/*
					Mdst_gamma jMdst_10 =(*j1).child(0).child(1).child(0).mdstGamma();
					Mdst_gamma jMdst_11 =(*j1).child(0).child(1).child(1).mdstGamma();
					Mdst_charged jMdst_0 =(*j1).child(0).child(0).mdstCharged();
					
					ExKFitterParticle jKF_10(jMdst_10);
					ExKFitterParticle jKF_11(jMdst_11);
					ExKFitterParticle jKF_0(jMdst_0, 3);
					
					HepPoint3D jkstar_init;
					
					jkstar_init.setX(IP.x() + (*j1).child(0).child(1).p().px()/(*j1).child(0).child(1).p().rho());
					jkstar_init.setY(IP.y() + (*j1).child(0).child(1).p().py()/(*j1).child(0).child(1).p().rho());
					jkstar_init.setZ(IP.z() + (*j1).child(0).child(1).p().pz()/(*j1).child(0).child(1).p().rho());
					ExKFitterVertex jkstar_Vertex(jkstar_init);
					ExKFitterVertex jB_Vertex(IP,IPerr);
					
					ExKFitterMass jPi0_Mass(0.134976);
					ExKFitterMass jkstar_Mass(0.8916);
					
					ExKFitterParticle jkstar;
					jkstar.LinkParticle(&jKF_10);
					jkstar.LinkParticle(&jKF_11);
					jkstar.LinkVertex(&jkstar_Vertex);
					
					
					ExKFitterConstrain jcon1;
					jcon1.SetMassConstrain();
					jcon1.LinkParticle(&jKF_10);
					jcon1.LinkParticle(&jKF_11);
					jcon1.LinkVertex(&jkstar_Vertex);
					jcon1.LinkMass(&jPi0_Mass);
					
					ExKFitterParticle jB;
					jB.LinkParticle(&jkstar);
					jB.LinkParticle(&jKF_0);
					jB.LinkVertex(&jB_Vertex);
					
					ExKFitterConstrain jcon2;
					jcon2.SetMassConstrain();
					jcon2.LinkParticle(&jkstar);
					jcon2.LinkParticle(&jKF_0);
					jcon2.LinkVertex(&jB_Vertex);
					jcon2.LinkMass(&jkstar_Mass);
					
					ExKFitter jCore;
					jCore.LinkConstrain(&jcon1);
					jCore.LinkConstrain(&jcon2);
					int jret = jCore.Minimize();
					float jchisqExK = jCore.Chisq();
					float jdof_exk = jCore.N_DegreeOfFreedom();
					HepPoint3D jdvertex = jkstar_Vertex.Vertex();
					HepPoint3D jbvertex1 = jB_Vertex.Vertex();
					
					if(jret==0)
					{
						jB.Update();
					}
					*/
					double chi1 = chisqExK/dof_exk;
					double chi2 = jchisqExK/jdof_exk;
					
					if (i == j1)
						continue;
					
					if (chi1 > chi2)
					{
						krankc1 = 0;
						break;
					}
					else if (chi1 < chi2)
					{
						krankc1 = 1;
					}
					else
						krankc1 = -1;
				}
				
				TagB_tpl->column( "krankc1"  , krankc1  );

				
				sigselection++;
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
				
				tnremaining = tnpip+tnpim+tnep+tnem+tnkp+tnkm+tnmup+tnmum+tnpm+tnpp  ;
				nremaining = npip+npim+nep+nem+nkp+nkm+nmup+nmum+npm+npp ;
				
				TagB_tpl->column("kstc0px" , (*i).child(0).child(0).px());
				TagB_tpl->column("kstc0py" , (*i).child(0).child(0).py());
				TagB_tpl->column("kstc0pz" , (*i).child(0).child(0).pz());
				TagB_tpl->column("kstc0mass" , (*i).child(0).child(0).p().mag());
				
				TagB_tpl->column("kstc1px" , (*i).child(0).child(1).px());
				TagB_tpl->column("kstc1py" , (*i).child(0).child(1).py());
				TagB_tpl->column("kstc1pz" , (*i).child(0).child(1).pz());
				TagB_tpl->column("kstc1mass" , (*i).child(0).child(1).p().mag());
				
				TagB_tpl->column("kstarpx" , (*i).child(0).px());
				TagB_tpl->column("kstarpy" , (*i).child(0).py());
				TagB_tpl->column("kstarpz" , (*i).child(0).pz());
				TagB_tpl->column("kstarmass" , (*i).child(0).p().mag());
				
				TagB_tpl->column( "nb_trackremaining"    , nremaining         );//new
				TagB_tpl->column( "tnb_trackremaining"    , tnremaining         );//new
				TagB_tpl->column( "tnb_pi0remaining"    , npi0         );//new
				TagB_tpl->column( "tnb_gamma"    , n_tgamma1         );//new
				TagB_tpl->column( "nbmulti1"    , nbmulti1         );//new
				TagB_tpl->column( "nb_pi015"    , npi015         );//new
				
				TagB_tpl->column( "n_kl_tot"    , n_kl_tot         );//new
				TagB_tpl->column( "n_kl"    , n_kl         );//new
				TagB_tpl->column( "n_ks"    , n_ks         );//new
	
//				double sdr,sdz;
//				const Mdst_charged* sch1 = &(*i).child(0).mdstCharged();
//				GetImpactParameters(sch1,&sdr,&sdz,3);
//				TagB_tpl->column( "ch0dr"     , sdr         );//new
//				TagB_tpl->column( "ch0dz"    , sdz         );//new
				
				double sum_ecl = 0;
				for(std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s !=eclaux_mag.end();s++ )
				{
					//               std::cout<<"start eclaux"<<std::endl;
					Mdst_ecl_aux& ch = *s;
					int flag = 0;
					int nfinal = Btag.relation().nFinalStateParticles(); // the number of final particle used to tagB
					int signfinal = (*i).relation().nFinalStateParticles();
					
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
					
					
					for(int kk = 0 ; kk < signfinal ; kk ++)
					{
						Particle sigchild = (*i).relation().finalStateParticle(kk);
						if ( sigchild.mdstCharged() )
						{
							if ( sigchild.mdstCharged().trk() == ch.trk())
							{
								flag=1;
							}
						}
						if (sigchild.mdstGamma())
						{
							if (sigchild.mdstGamma().ecl().get_ID()==ch.get_ID())
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
				HepLorentzVector qlab = ( Bsigplab - k_plab); // signal side miss
				double Bsigmass2 = Bsigplab.mag()*Bsigplab.mag();//sigbmass^2
				
				double lab_total_miss = pmag(tmissinglab);
				double lab_sig_miss = pmag(qlab);
				double tmmiss2lab = tmissinglab.mag()*tmissinglab.mag();
				double qsquare = qlab.mag()*qlab.mag();
				double sblab = qsquare / Bsigmass2;
				double lab_k_p = pmag(k_plab);
				//////////////////////        total miss information        ////////////////////
				TagB_tpl->column("tmmisslab",tmissinglab.mag());
				TagB_tpl->column("tmmiss2lab",tmmiss2lab);
				TagB_tpl->column("tpxmisslab",tmissinglab.px());
				TagB_tpl->column("tpymisslab",tmissinglab.py());
				TagB_tpl->column("tpzmisslab",tmissinglab.pz());
				TagB_tpl->column("temisslab",tmissinglab.e());
				TagB_tpl->column("tpmisslab",lab_total_miss);
				//////////////////////        signal miss information        ////////////////////
				TagB_tpl->column("sblab",sblab);
				TagB_tpl->column("qsquare",qsquare);
				TagB_tpl->column("qpxlab",qlab.px());
				TagB_tpl->column("qpylab",qlab.py());
				TagB_tpl->column("qpzlab",qlab.pz());
				TagB_tpl->column("qelab",qlab.e());
				TagB_tpl->column("spmisslab",lab_sig_miss);
				//////////////////////        K information        ////////////////////
				
				TagB_tpl->column("kmalab",k_plab.mag());
				TagB_tpl->column("kpxlab",k_plab.px());
				TagB_tpl->column("kpylab",k_plab.py());
				TagB_tpl->column("kpzlab",k_plab.pz());
				TagB_tpl->column("kelab",k_plab.e());
				TagB_tpl->column("kplab",lab_k_p);
				//std::cout<<"p1"<<std::endl;
				//Lab frame miss cos
				double tcosmisslabinz;
				double tangelmisslabinz;
				tcosmisslabinz = ((tmissinglab.pz())/lab_total_miss);
				tangelmisslabinz = acos((tmissinglab.pz())/lab_total_miss);
				TagB_tpl->column("tcosmizl",tcosmisslabinz);
				//h(*) B+ cos lab
				HepLorentzVector Btagplab (Btag.px(),Btag.py(),Btag.pz(),Btag.e());
				HepLorentzVector all_BEAM_polar ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle),E_HER+E_LER);
				HepLorentzVector P_BEAM_polar = (all_BEAM_polar - Btagplab);//BEAM - BTAG = SINGAL
				double CM_2 = pmag(P_BEAM_polar);
				double kp2 = pmag((*i).child(0));
				double cosksblab = costheta((*i).child(0),P_BEAM_polar);
				double cosklab = (*i).child(0).pz()/kp2;
				TagB_tpl->column( "cosksblab"   , cosksblab);//new
				TagB_tpl->column( "cosklab"   , cosklab);//new
				
				///////////////////////////////////////////////////////////////////////////////////////
				//                            Boost to upsilon 4S rest frame                         //
				///////////////////////////////////////////////////////////////////////////////////////
				
				tagB.boost( boost_vector.boostVector() );//boost tagb
				double CM_Btag_p = sqrt( tagB.px()*tagB.px() + tagB.py()*tagB.py() + tagB.pz()*tagB.pz() );
				double CM_Btag_pz = sqrt( tagB.pz()*tagB.pz() );
				//	 double CM_Btag_ppz = tagB.pz()*tagB.px() + tagB.pz()*tagB.py() + tagB.pz()*tagB.pz() ;
				double cosbcm = tagB.pz()/(CM_Btag_p);
				TagB_tpl->column( "cosbcm"     , cosbcm );
				TagB_tpl->column( "tagBpcm"    , CM_Btag_p );
				TagB_tpl->column( "tagbPXcm"   , tagB.px());
				TagB_tpl->column( "tagbPYcm"   , tagB.py());
				TagB_tpl->column( "tagbPZcm"   , tagB.pz());
				HepLorentzVector BEAM(E_HER*sin(cross_angle), 0.0, E_HER*cos(cross_angle)-E_LER, E_HER+E_LER);
				HepLorentzVector k_p((*i).child(0).px(),(*i).child(0).py(),(*i).child(0).pz(),(*i).child(0).e());
				HepLorentzVector Bsigpcm (-tagB.px(),-tagB.py(),-tagB.pz(),tagB.e());
				BEAM.boost(boost_vector.boostVector());//boost beam
				k_p.boost(boost_vector.boostVector() );//boost signal k
				HepLorentzVector totalmissincm = ( BEAM - tagB - k_p); //total miss
				HepLorentzVector sigmisscm = (  Bsigpcm - k_p ); //signal miss
				double cm_totalmiss_p =pmag(totalmissincm);
				double bsigpcm =pmag(Bsigpcm) ;
				double sigmisspcm = pmag(sigmisscm);
				double anglebsigandsigmiss = acos((Bsigpcm.px()*sigmisscm.px() + Bsigpcm.py()*sigmisscm.py() + Bsigpcm.pz()*sigmisscm.pz())/(bsigpcm*sigmisspcm));//new
				double cosbsigandsigmiss = (Bsigpcm.px()*sigmisscm.px() + Bsigpcm.py()*sigmisscm.py() + Bsigpcm.pz()*sigmisscm.pz())/(bsigpcm*sigmisspcm);//new
				
				TagB_tpl->column("angbssmi", anglebsigandsigmiss);//new
				TagB_tpl->column("cosbssmi",cosbsigandsigmiss);//new
				/////////////////////        signal miss information        ////////////////////
				TagB_tpl->column("sigmissp",sigmisscm.mag());
				double sigmmiss2 = sigmisscm.mag()*sigmisscm.mag();
				TagB_tpl->column("Bsigmas2",Bsigmass2);
				double sbcm = sigmmiss2 /Bsigmass2  ;
				TagB_tpl->column("sbcm",sbcm);
				TagB_tpl->column("simmis2CM",sigmmiss2);
				TagB_tpl->column("spxmisisC",sigmisscm.px());
				TagB_tpl->column("spymissC",sigmisscm.py());
				TagB_tpl->column("spzmissC",sigmisscm.pz());
				TagB_tpl->column("semissC",sigmisscm.e());
				TagB_tpl->column("spmissC",sigmisspcm);
				//////////////////////        K information        ////////////////////
				
				TagB_tpl->column("mkcm",k_p.mag());
				TagB_tpl->column("kpxCM",k_p.px());
				TagB_tpl->column("kpyCM",k_p.py());
				TagB_tpl->column("kpzCM",k_p.pz());
				TagB_tpl->column("keCM",k_p.e());
				double CM_k_p = pmag(k_p);
				TagB_tpl->column("kpCM",CM_k_p);
				//////////////////////        total miss information        ////////////////////
				TagB_tpl->column("tmmisscm",totalmissincm.mag());
				TagB_tpl->column("tmmiss2cm",totalmissincm.mag()*totalmissincm.mag());
				TagB_tpl->column("tpxmisscm",totalmissincm.px());
				TagB_tpl->column("tpymisscm",totalmissincm.py());
				TagB_tpl->column("tpzmisscm",totalmissincm.pz());
				TagB_tpl->column("temisscm",totalmissincm.e());
				TagB_tpl->column("tpmisscm",cm_totalmiss_p);
				
				///////////////////////////////////////////////////////////////////////////////////////
				//                            Boost to upsilon B rest frame                         //
				///////////////////////////////////////////////////////////////////////////////////////
				HepLorentzVector sigmissb(sigmisscm);
				HepLorentzVector k_b(k_p);
				HepLorentzVector totalmissinb(totalmissincm);
				
				sigmissb.boost(-Bsigplab.boostVector());//signal miss
				totalmissinb.boost(-Bsigplab.boostVector());//total miss
				
				k_b.boost(-Bsigplab.boostVector());//signl K
				double sigmissm2b = sigmissb.mag()*sigmissb.mag();
				double sigmisspb = pmag(sigmissb);
				/////////////////////        signal miss information        ////////////////////
				TagB_tpl->column("smmissb",sigmissb.mag());
				TagB_tpl->column("smmiss2b",sigmissm2b);
				TagB_tpl->column("spxmissb",sigmissb.px());
				TagB_tpl->column("spymissb",sigmissb.py());
				TagB_tpl->column("spzmissb",sigmissb.pz());
				TagB_tpl->column("semissb",sigmissb.e());
				TagB_tpl->column("spmissb",sigmisspb);
				double b_k_p = pmag(k_b);
				//////////////////////        K information        ////////////////////
				TagB_tpl->column("kpb", b_k_p);
				TagB_tpl->column("kpxb",k_b.px());
				TagB_tpl->column("kpyb",k_b.py());
				TagB_tpl->column("kpzb",k_b.pz());
				TagB_tpl->column("keb",k_b.e());
				//////////////////////        total miss information        ////////////////////
				TagB_tpl->column("tmmissb",totalmissinb.mag());
				TagB_tpl->column("tmmiss2b",totalmissinb.mag()*totalmissinb.mag());
				TagB_tpl->column("tpxmissb",totalmissinb.px());
				TagB_tpl->column("tpymissb",totalmissinb.py());
				TagB_tpl->column("tpzmissb",totalmissinb.pz());
				TagB_tpl->column("temissb",totalmissinb.e());
				TagB_tpl->column("tpmissb",pmag(totalmissinb));
				
				//////////////////////////////////dz///////////////////////////////////////////////////////
				HepPoint3D Btag_init;
				HepPoint3D bvertex;
				HepPoint3D bvertexsig;
				float vchisq_tag;
				int DOF_tag;
				float pvalue_tag;
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
				if(flagch)
				{
					kvertexfitter kvf;
					Mdst_charged bfinalnumber[nfinal];
					ExKFitterParticle btagf[nfinal];
					Particle bfinal[nfinal];
					for(int k2 = 0 ; k2 < nfinal; k2++)
					{
						int expid = 0;
						Particle bfinal_tmp = Btag.relation().finalStateParticle(k2);
						bfinal[k2]=bfinal_tmp;
						int nunp = -1 ;
						if(!(bfinal[k2].mdstCharged()))
							continue;
						addTrack2fit(kvf,  bfinal[k2]);
						
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
					//  std::cout<<"start k6.9"<<std::endl;
				}
				double dr,dz;
				int signfinal = (*i).relation().nFinalStateParticles(); // the number of final particle used to tagB
				Mdst_charged sigbfinalnumber[signfinal];
				ExKFitterParticle bsisgf[signfinal];
				Particle bsigfinal[signfinal];
				kvertexfitter kvfsig;
				for(int k3 = 0 ; k3 < signfinal; k3++)
				{
					int expid = 0;
					Particle bsigfinal_tmp = (*i).relation().finalStateParticle(k3);
					bsigfinal[k3]=bsigfinal_tmp;
					int nunp = -1 ;
					if(!(bsigfinal[k3].mdstCharged()))
						continue;
					addTrack2fit(kvfsig,  bsigfinal[k3] );
					
				}

//				const Mdst_charged* ch1 = &(*i).child(0).mdstCharged();
//				Mdst_charged Mdst_k = (*i).child(0).mdstCharged();
//				Particle sig_k = (*i).child(0);
				
//				addTrack2fit(kvfsig,sig_k);
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
							if ((*p1).theta() < 0.547 && (*p1).energy() > 0.1 && fabs(tmissinglab.angle(EclVirtualVec.vect()))<0.1 )
								sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
							//barrel e>0.05
							if ((*p1).theta() > 0.562 && (*p1).theta() < 2.246 && (*p1).energy()>0.05&& fabs(tmissinglab.angle(EclVirtualVec.vect()))< 0.1 )
								sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
							//backward e>0.15
							if ((*p1).theta() > 2.246 && (*p1).energy() > 0.15&& fabs(tmissinglab.angle(EclVirtualVec.vect()))< 0.1)
								sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
						}
						EclVirtualVec.setRho(1.0);EclVirtualVec.setTheta((*p1).theta());EclVirtualVec.setPhi((*p1).phi());
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
				
				///////////////////////////////////////////////////////////////////////////////////////
				//                      MC TRuth(I don't know it work or not)                        //
				///////////////////////////////////////////////////////////////////////////////////////
				
				if(MCstatus == 1)
				{
					int check = -1;
					int hhindex= -1;
					int kstcheck = -1;
					///////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////
					
					
					if (frecB.tag_id() == -521)
					{
						int bnfinal = (*i).child(0).relation().nFinalStateParticles();
						for(int k = 0 ; k < bnfinal; k++)
						{
							Particle BB = (*i);
							Particle bchild = (*i).child(0).relation().finalStateParticle(k);
							if ( bchild.mdstCharged() )
							{
								const Mdst_charged c11h000 = bchild.mdstCharged();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( (Evt_decay.idhep() == 323  && Evt_current.idhep() == 521))
									{
										check = 1;
										break;
									}
									else
										check = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							else if ( bchild.mdstGamma() )
							{
								const Mdst_gamma c11h000 = bchild.mdstGamma();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( (Evt_decay.idhep() == 323  && Evt_current.idhep() == 521))
									{
										check = 1;
										break;
									}
									else
										check = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							if (check != 1)
								break;
						}//kend
						TagB_tpl->column("check" ,check);
						int doublecheck = -1 ;
						for (std::vector<Gen_hepevt>::iterator it3 = gen_mgr.begin();it3 != gen_mgr.end(); it3++)
						{
							if((*it3).idhep()== 521 )
							{
								Gen_hepevt d = gen_mgr[(*it3).da(0)-1];
								//std::cout << "Evt_dd =" << d.idhep() << std::endl;
								if ( d.idhep() == 323)
								{
									for(int k = 0 ; k < bnfinal; k++)
									{
										doublecheck = 0 ;
										Particle bchild = (*i).child(0).relation().finalStateParticle(k);
										if ( bchild.mdstCharged() )
										{
											const Mdst_charged dsub = bchild.mdstCharged();
											Gen_hepevt  Evt_Ddecay= get_hepevt(dsub);
											toptodown(Evt_Ddecay,d,doublecheck,gen_mgr);
										}
										else if ( bchild.mdstGamma() )
										{
											const Mdst_gamma dsub = bchild.mdstGamma();
											Gen_hepevt Evt_Ddecay= get_hepevt(dsub);
											toptodown(Evt_Ddecay,d,doublecheck,gen_mgr);
										}
									}
								}
								else
									break;
							}
						}
						if ( check == 1 && doublecheck == 1 )
							hhindex = 1;
						else
							hhindex =0;
						TagB_tpl->column("doublecheck" ,doublecheck);
						TagB_tpl->column("hhindex" ,hhindex);
						
					}
					else if (frecB.tag_id() == 521)
					{
						int bnfinal = (*i).child(0).relation().nFinalStateParticles();
						for(int k = 0 ; k < bnfinal; k++)
						{
							Particle BB = (*i);
							Particle bchild = (*i).child(0).relation().finalStateParticle(k);
							if ( bchild.mdstCharged() )
							{
								const Mdst_charged c11h000 = bchild.mdstCharged();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( (Evt_decay.idhep() == -323  && Evt_current.idhep() == -521))
									{
										check = 1;
										break;
									}
									else
										check = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							else if ( bchild.mdstGamma() )
							{
								const Mdst_gamma c11h000 = bchild.mdstGamma();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( (Evt_decay.idhep() == -323  && Evt_current.idhep() == -521))
									{
										check = 1;
										break;
									}
									else
										check = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
						}//kend
						TagB_tpl->column("check" ,check);
						int doublecheck = -1 ;
						for (std::vector<Gen_hepevt>::iterator it3 = gen_mgr.begin();it3 != gen_mgr.end(); it3++)
						{
							if((*it3).idhep()== -521 )
							{
								Gen_hepevt d = gen_mgr[(*it3).da(0)-1];
								//std::cout << "Evt_dd =" << d.idhep() << std::endl;
								if ( d.idhep() == -323)
								{
									for(int k = 0 ; k < bnfinal; k++)
									{
										doublecheck = 0 ;
										Particle bchild = (*i).child(0).relation().finalStateParticle(k);
										if ( bchild.mdstCharged() )
										{
											const Mdst_charged dsub = bchild.mdstCharged();
											Gen_hepevt  Evt_Ddecay= get_hepevt(dsub);
											toptodown(Evt_Ddecay,d,doublecheck,gen_mgr);
										}
										else if ( bchild.mdstGamma() )
										{
											const Mdst_gamma dsub = bchild.mdstGamma();
											Gen_hepevt Evt_Ddecay= get_hepevt(dsub);
											toptodown(Evt_Ddecay,d,doublecheck,gen_mgr);
										}
									}
								}
								else
									break;
							}
						}
						if ( check == 1 && doublecheck == 1 )
							hhindex = 2;
						else
							hhindex =0;
						TagB_tpl->column("doublecheck" ,doublecheck);
						TagB_tpl->column("hhindex" ,hhindex);
						
					}
					///////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////
						int bnfinal = (*i).child(0).relation().nFinalStateParticles();
						for(int k = 0 ; k < bnfinal; k++)
						{
							Particle BB = (*i);
							Particle bchild = (*i).child(0).relation().finalStateParticle(k);
							if ( bchild.mdstCharged() )
							{
								const Mdst_charged c11h000 = bchild.mdstCharged();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( Evt_decay.idhep() == 323)
									{
										kstcheck = 1;
										break;
									}
									else
										kstcheck = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							else if ( bchild.mdstGamma() )
							{
								const Mdst_gamma c11h000 = bchild.mdstGamma();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( Evt_decay.idhep()== 323 )
									{
										kstcheck = 1;
										break;
									}
									else
										kstcheck = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							if (kstcheck != 1)
								break;
						}//kend
					
					if (kstcheck != 1)
					{
						for(int k = 0 ; k < bnfinal; k++)
						{
							Particle BB = (*i);
							Particle bchild = (*i).child(0).relation().finalStateParticle(k);
							if ( bchild.mdstCharged() )
							{
								const Mdst_charged c11h000 = bchild.mdstCharged();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( Evt_decay.idhep() == -323)
									{
										kstcheck = -1;
										break;
									}
									else
										kstcheck = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							else if ( bchild.mdstGamma() )
							{
								const Mdst_gamma c11h000 = bchild.mdstGamma();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( Evt_decay.idhep()== -323 )
									{
										kstcheck = -1;
										break;
									}
									else
										kstcheck = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							if (kstcheck != -1)
								break;
						}//kend
					}
					
					TagB_tpl->column("kstcheck" ,kstcheck);
						
					
					
					///////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////
					
					if (mode == 1) // k+ pi0
					{
						const Mdst_charged ch0 = (*i).child(0).child(0).mdstCharged();
						const Mdst_gamma ch10 = (*i).child(0).child(1).child(0).mdstGamma();
						const Mdst_gamma ch11 = (*i).child(0).child(1).child(1).mdstGamma();
						
						Gen_hepevt Evtch0 = get_hepevt(ch0);
						Gen_hepevt Evtch10 = get_hepevt(ch10);
						Gen_hepevt Evtch11 = get_hepevt(ch11);
						
						Gen_hepevt tEvtch0;
						Gen_hepevt tEvtch10;
						Gen_hepevt tEvtch11;
						std::cout<< "debug0 ch0id = " << Evtch0.idhep() <<std::endl;
						std::cout<< "debug0 ch10id = " << Evtch10.idhep() <<std::endl;
						std::cout<< "debug0 ch11id = " << Evtch11.idhep() <<std::endl;
						
//						while( abs(Evtch0.idhep()) != 321 && Evtch0.mo(0))
//						{
//							tEvtch0 = gen_mgr[Evtch0.mo(0)-1];
//							std::cout<< "debug1 ch0id = " << tEvtch0.idhep() <<std::endl;
//							Evtch0 = tEvtch0 ;
//						}
/*						while( abs(Evtch10.idhep()) != 22  && Evtch10.mo(0))
						{
							tEvtch10 = gen_mgr[Evtch10.mo(0)-1];
							std::cout<< "debug1 ch10id = " << tEvtch10.idhep() <<std::endl;
							Evtch10 = tEvtch10 ;
						}
						while( abs(Evtch11.isthep()) != 22  && Evtch11.mo(0))
						{
							tEvtch11 = gen_mgr[Evtch11.mo(0)-1];
							std::cout<< "debug1 ch11id = " << tEvtch11.idhep() <<std::endl;
							Evtch11 = tEvtch11 ;
						}
*/
						TagB_tpl->column("hi0",Evtch0.idhep());
						TagB_tpl->column("hi10",Evtch10.idhep());
						TagB_tpl->column("hi11",Evtch11.idhep());
						
						Gen_hepevt EvtP0;
						Gen_hepevt EvtGP0;
						Gen_hepevt EvtGGP0;
						
						Gen_hepevt EvtP10;
						Gen_hepevt EvtGP10;
						Gen_hepevt EvtGGP10;
						
						Gen_hepevt EvtP11;
						Gen_hepevt EvtGP11;
						Gen_hepevt EvtGGP11;
						
						if(Evtch0.mo(0))
						{
							EvtP0 = gen_mgr[Evtch0.mo(0)-1];
							TagB_tpl->column("mo0",EvtP0.idhep());
							
							if(EvtP0.mo(0))
							{
								EvtGP0 = gen_mgr[EvtP0.mo(0)-1];
								TagB_tpl->column("gmo0",EvtGP0.idhep());
								
								
								if(EvtGP0.mo(0))
								{
									EvtGGP0 = gen_mgr[EvtGP0.mo(0)-1];
									TagB_tpl->column("ggmo0",EvtGGP0.idhep());
								}
							}
						}
						if(Evtch10.mo(0))
						{
							EvtP10 = gen_mgr[Evtch10.mo(0)-1];
							TagB_tpl->column("mo10",EvtP10.idhep());
							
							if(EvtP10.mo(0))
							{
								EvtGP10 = gen_mgr[EvtP10.mo(0)-1];
								TagB_tpl->column("gmo10",EvtGP10.idhep());
								
								
								if(EvtGP10.mo(0))
								{
									EvtGGP10 = gen_mgr[EvtGP10.mo(0)-1];
									TagB_tpl->column("ggmo10",EvtGGP10.idhep());
								}
							}
						}
						if(Evtch11.mo(0))
						{
							EvtP11 = gen_mgr[Evtch11.mo(0)-1];
							TagB_tpl->column("mo11",EvtP11.idhep());
							
							if(EvtP11.mo(0))
							{
								EvtGP11 = gen_mgr[EvtP11.mo(0)-1];
								TagB_tpl->column("gmo11",EvtGP11.idhep());
								
								
								if(EvtGP11.mo(0))
								{
									EvtGGP11 = gen_mgr[EvtGP11.mo(0)-1];
									TagB_tpl->column("ggmo11",EvtGGP11.idhep());
								}
							}
						}
						int hindex = -1 ;
						
						if( Evtch0.idhep()== 321 && Evtch10.idhep()==22 && Evtch11.idhep()==22 && EvtP0.idhep()==323 && EvtP10.idhep()==111 && EvtP11.idhep()==111 && EvtGP0.idhep()==521 && EvtGP10.idhep()==323 && EvtGP11.idhep()==323 && EvtGGP10.idhep()==521 && EvtGGP11.idhep()==521 )
						{
							hindex = 1;
						}
						else if( Evtch0.idhep()== -321 && Evtch10.idhep()==22 && Evtch11.idhep()==22 && EvtP0.idhep()==-323 && EvtP10.idhep()==111 && EvtP11.idhep()==111 && EvtGP0.idhep()==-521 && EvtGP10.idhep()==-323 && EvtGP11.idhep()==-323 && EvtGGP10.idhep()==-521 && EvtGGP11.idhep()==-521 )
						{
							hindex = 2;
						}
						else
						{
							hindex = 0;
						}
						
						TagB_tpl->column("hindex",hindex);
						
					}
					else if (mode == 2) // ks pi+
					{
						std::cout<< "4" << std::endl;
						std::cout<< "id0 =" << (*i).child(0).child(1).lund() <<std::endl;
						const Mdst_charged ch0 = (*i).child(0).child(1).mdstCharged();
						
						std::cout<< "id1 =" << (*i).child(0).child(0).child(0).lund() <<std::endl;
						std::cout<< "id2 =" << (*i).child(0).child(0).child(1).lund()<<std::endl;
						
						const Mdst_charged ch10 = (*i).child(0).child(0).child(0).mdstCharged();
						const Mdst_charged ch11 = (*i).child(0).child(0).child(1).mdstCharged();
						std::cout<< "5" << std::endl;
						Gen_hepevt Evtch0 = get_hepevt(ch0);
						Gen_hepevt Evtch10 = get_hepevt(ch10);
						Gen_hepevt Evtch11 = get_hepevt(ch11);
						
						Gen_hepevt tEvtch0 ;
						Gen_hepevt tEvtch10;
						Gen_hepevt tEvtch11;
						
/*						while( abs(Evtch0.idhep()) != 211 && Evtch0.mo(0) )
						{
							tEvtch0 = gen_mgr[Evtch0.mo(0)-1];
							Evtch0 = tEvtch0 ;
						}
						while( abs(Evtch10.isthep()) != 211 && Evtch10.mo(0))
						{
							tEvtch10 = gen_mgr[Evtch10.mo(0)-1];
							Evtch10 = tEvtch10 ;
						}
						while( abs(Evtch11.isthep()) != 211 && Evtch11.mo(0))
						{
							tEvtch11 = gen_mgr[Evtch11.mo(0)-1];
							Evtch11 = tEvtch11 ;
						}
*/						
						std::cout<< "6" << std::endl;
						TagB_tpl->column("hi0",Evtch0.idhep());
						TagB_tpl->column("hi10",Evtch10.idhep());
						TagB_tpl->column("hi11",Evtch11.idhep());
						
						Gen_hepevt EvtP0;
						Gen_hepevt EvtGP0;
						Gen_hepevt EvtGGP0;
						
						Gen_hepevt EvtP10;
						Gen_hepevt EvtGP10;
						Gen_hepevt EvtGGP10;
						Gen_hepevt EvtGGGP10;
						
						Gen_hepevt EvtP11;
						Gen_hepevt EvtGP11;
						Gen_hepevt EvtGGP11;
						Gen_hepevt EvtGGGP11;
						std::cout<< "7" << std::endl;
						if(Evtch0.mo(0))
						{
							EvtP0 = gen_mgr[Evtch0.mo(0)-1];
							TagB_tpl->column("mo0",EvtP0.idhep());
							
							if(EvtP0.mo(0))
							{
								EvtGP0 = gen_mgr[EvtP0.mo(0)-1];
								TagB_tpl->column("gmo0",EvtGP0.idhep());
								
								
								if(EvtGP0.mo(0))
								{
									EvtGGP0 = gen_mgr[EvtGP0.mo(0)-1];
									TagB_tpl->column("ggmo0",EvtGGP0.idhep());
									//								if(EvtGGP0.mo(0))
									//								{
									//									EvtGGGP0 = gen_mgr[EvtGGP0.mo(0)-1];
									//									TagB_tpl->column("gggmo0",EvtGGGP0.idhep());
									//								}
								}
							}
						}
						std::cout<< "8" << std::endl;
						if(Evtch10.mo(0))
						{
							EvtP10 = gen_mgr[Evtch10.mo(0)-1];
							TagB_tpl->column("mo10",EvtP10.idhep());
							
							if(EvtP10.mo(0))
							{
								EvtGP10 = gen_mgr[EvtP10.mo(0)-1];
								TagB_tpl->column("gmo10",EvtGP10.idhep());
								
								
								if(EvtGP10.mo(0))
								{
									EvtGGP10 = gen_mgr[EvtGP10.mo(0)-1];
									TagB_tpl->column("ggmo10",EvtGGP10.idhep());
									if(EvtGGP10.mo(0))
									{
										EvtGGGP10 = gen_mgr[EvtGGP10.mo(0)-1];
										TagB_tpl->column("gggmo10",EvtGGGP10.idhep());
									}
								}
							}
						}
						std::cout<< "9" << std::endl;
						if(Evtch11.mo(0))
						{
							EvtP11 = gen_mgr[Evtch11.mo(0)-1];
							TagB_tpl->column("mo11",EvtP11.idhep());
							
							if(EvtP11.mo(0))
							{
								EvtGP11 = gen_mgr[EvtP11.mo(0)-1];
								TagB_tpl->column("gmo11",EvtGP11.idhep());
								
								
								if(EvtGP11.mo(0))
								{
									EvtGGP11 = gen_mgr[EvtGP11.mo(0)-1];
									TagB_tpl->column("ggmo11",EvtGGP11.idhep());
									if(EvtGGP11.mo(0))
									{
										EvtGGGP11 = gen_mgr[EvtGGP11.mo(0)-1];
										TagB_tpl->column("gggmo11",EvtGGGP11.idhep());
										
									}
								}
							}
						}
						int hindex = -1 ;
						std::cout<< "10" << std::endl;
						if( Evtch0.idhep()== 211 && Evtch10.idhep()==211 && Evtch11.idhep()==-211 && EvtP0.idhep()==323 && EvtP10.idhep()==310 && EvtP11.idhep()==310 && EvtGP0.idhep()==521 && abs(EvtGP10.idhep())==311 && abs(EvtGP11.idhep())==311 && EvtGGP10.idhep()==323 && EvtGGP11.idhep()==323 && EvtGGGP10.idhep()==521 && EvtGGGP11.idhep()==521  )
						{
							hindex = 1;
						}
						else if( Evtch0.idhep()== -211 && Evtch10.idhep()==211 && Evtch11.idhep()==-211 && EvtP0.idhep()==-323 && EvtP10.idhep()==310 && EvtP11.idhep()==310 && EvtGP0.idhep()==-521 && abs(EvtGP10.idhep())==311 && abs(EvtGP11.idhep())==311 && EvtGGP10.idhep()==-323 && EvtGGP11.idhep()==-323 && EvtGGGP10.idhep()==-521 && EvtGGGP11.idhep()==-521 )
						{
							hindex = 2;
						}
						else
						{
							hindex = 0;
						}
						
						TagB_tpl->column("hindex",hindex);
						
					}
				}//truth event  end
				//        std::cout<<"start truth event"<<std::endl;
			
				/////////////////////////////////trace form top to down//////////////////////////////////
				int nbpd=0,nbnd=0,topcheck= 0;
				for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin();i != gen_mgr.end(); i++)
				{
					if (abs((*i).idhep())==521)
					{
						int da1=(*i).da(0),da2=(*i).da(1);
						
						if ((*i).idhep() == 521)
							TagB_tpl->column("charged",1);
						else if ((*i).idhep() == -521)
							TagB_tpl->column("charged",-1);
						else
							TagB_tpl->column("charged",0);
						
						nbpd=(da2-da1+1);
						TagB_tpl->column("nbpd",nbpd);
						Gen_hepevt kstar=gen_mgr[(*i).da(0)-1];
						Gen_hepevt nu=gen_mgr[(*i).da(0)];
						Gen_hepevt nubar=gen_mgr[(*i).da(0)+1];
						
						if ((*i).idhep() == 521)
						{
							for(int start=0;start<(da2-da1+1);start++)
							{
								Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
								char bpdnumber[32];
								sprintf (bpdnumber,"%s%d","bpd",start+1);
								TagB_tpl->column(bpdnumber,Evda.idhep());
							}

								if (kstar.idhep() == 323 && nu.idhep() == 18 && nubar.idhep() == -18)
									TagB_tpl->column("topcheck",1);
							
						}
						else if ((*i).idhep() == -521)
						{
							for(int start=0;start<(da2-da1+1);start++)
							{
								Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
								char bndnumber[32];
								sprintf (bndnumber,"%s%d","bnd",start+1);//bnd1,bnd2,bnd3....
								TagB_tpl->column(bndnumber,Evda.idhep());
							}

								if (kstar.idhep() == -323 && nu.idhep() == 18 && nubar.idhep() == -18)
									TagB_tpl->column("topcheck",-1);
							
						}
						
					}

				}
				
				TagB_tpl->dumpData();
				*status = 1;
			}//k_plus end/////////////////////////////////////////////////////////////////////////
			
			for(std::vector<Particle>::iterator i =B_cand2.begin(); i != B_cand2.end(); i++)
			{
				for(std::vector<Particle>::iterator j1=B_cand2.begin(); j1 != B_cand2.end(); j1++)
				{
					if (i == j1)
						continue;
					
					double massd1;
					double massd2;
					massd1 = fabs( 0.8916 - (*i).p().mag() );
					massd2 = fabs( 0.8916 - (*j1).p().mag() );
					if (massd1 > massd2)
					{
						kstrank2 = 0;
						break;
					}
					else if (massd1 < massd2)
					{
						kstrank2 = 1;
					}
					else
						kstrank2 = -1;
				}
				
				TagB_tpl2->column( "kstrank2"  , kstrank2  );
				
				double mode;
				if ((*i).child(0).child(0).mdstCharged())
				{
					mode = 1;
				}
				else
				{
					mode = 2;
				}
				TagB_tpl2->column( "mode"  , mode  );
				
				std::cout<< "mode = " << mode <<std::endl;
				// b_tpl->column( "mode"  , mode  );
				
				
				Mdst_charged Mdst_00 =(*i).child(0).child(0).child(0).mdstCharged();
				Mdst_charged Mdst_01 =(*i).child(0).child(0).child(1).mdstCharged();
				Mdst_charged Mdst_1 =(*i).child(0).child(1).mdstCharged();
				
				ExKFitterParticle KF_00(Mdst_00, 2);
				ExKFitterParticle KF_01(Mdst_01, 2);
				ExKFitterParticle KF_1(Mdst_1, 2);
				
				HepPoint3D ks_init;
				
				ks_init.setX(IP.x() + (*i).child(0).child(0).p().px()/(*i).child(0).child(0).p().rho());
				ks_init.setY(IP.y() + (*i).child(0).child(0).p().py()/(*i).child(0).child(0).p().rho());
				ks_init.setZ(IP.z() + (*i).child(0).child(0).p().pz()/(*i).child(0).child(0).p().rho());
				ExKFitterVertex ks_Vertex(ks_init);
				ExKFitterVertex B_Vertex(IP,IPerr);
				
				ExKFitterParticle ks;
				ks.LinkParticle(&KF_00);
				ks.LinkParticle(&KF_01);
				ks.LinkVertex(&ks_Vertex);
				ExKFitterConstrain con1;
				con1.SetVertexConstrain();
				con1.LinkParticle(&KF_00);
				con1.LinkParticle(&KF_01);
				con1.LinkVertex(&ks_Vertex);
				
				ExKFitterParticle B;
				B.LinkParticle(&ks);
				B.LinkParticle(&KF_1);
				B.LinkVertex(&B_Vertex);
				ExKFitterConstrain con2;
				con2.SetVertexConstrain();
				con2.LinkParticle(&ks);
				con2.LinkParticle(&KF_1);
				con2.LinkVertex(&B_Vertex);
				
				ExKFitter Core;
				Core.LinkConstrain(&con1);
				Core.LinkConstrain(&con2);
				int ret = Core.Minimize();
				float chisqExK = Core.Chisq();
				float dof_exk = Core.N_DegreeOfFreedom();
				HepPoint3D dvertex = ks_Vertex.Vertex();
				HepPoint3D bvertex1 = B_Vertex.Vertex();
				
				if(ret==0)
				{
					B.Update();
				}
				
				TagB_tpl2->column("chisqexk",chisqExK/dof_exk);
				for(std::vector<Particle>::iterator j1=B_cand2.begin(); j1 != B_cand2.end(); j1++)
				{
					if (i == j1)
						continue;
					
					double chi1;
					double chi2;
					Mdst_charged jMdst_00 =(*j1).child(0).child(0).child(0).mdstCharged();
					Mdst_charged jMdst_01 =(*j1).child(0).child(0).child(1).mdstCharged();
					Mdst_charged jMdst_1 =(*j1).child(0).child(1).mdstCharged();
					
					ExKFitterParticle jKF_00(Mdst_00, 2);
					ExKFitterParticle jKF_01(Mdst_01, 2);
					ExKFitterParticle jKF_1(Mdst_1, 2);
					
					HepPoint3D jks_init;
					
					jks_init.setX(IP.x() + (*j1).child(0).child(0).p().px()/(*j1).child(0).child(0).p().rho());
					jks_init.setY(IP.y() + (*j1).child(0).child(0).p().py()/(*j1).child(0).child(0).p().rho());
					jks_init.setZ(IP.z() + (*j1).child(0).child(0).p().pz()/(*j1).child(0).child(0).p().rho());
					ExKFitterVertex jks_Vertex(jks_init);
					ExKFitterVertex jB_Vertex(IP,IPerr);
					
					ExKFitterParticle jks;
					jks.LinkParticle(&jKF_00);
					jks.LinkParticle(&jKF_01);
					jks.LinkVertex(&jks_Vertex);
					ExKFitterConstrain jcon1;
					jcon1.SetVertexConstrain();
					jcon1.LinkParticle(&jKF_00);
					jcon1.LinkParticle(&jKF_01);
					jcon1.LinkVertex(&jks_Vertex);
					
					ExKFitterParticle jB;
					jB.LinkParticle(&jks);
					jB.LinkParticle(&jKF_1);
					jB.LinkVertex(&jB_Vertex);
					ExKFitterConstrain jcon2;
					jcon2.SetVertexConstrain();
					jcon2.LinkParticle(&jks);
					jcon2.LinkParticle(&jKF_1);
					jcon2.LinkVertex(&jB_Vertex);
					
					ExKFitter jCore;
					jCore.LinkConstrain(&jcon1);
					jCore.LinkConstrain(&jcon2);
					int jret = jCore.Minimize();
					float jchisqExK = jCore.Chisq();
					float jdof_exk = jCore.N_DegreeOfFreedom();
					HepPoint3D jdvertex = jks_Vertex.Vertex();
					HepPoint3D jbvertex1 = jB_Vertex.Vertex();
					
					if(jret==0)
					{
						jB.Update();
					}
					if (chi1 > chi2)
					{
						krankc2 = 0;
						break;
					}
					else if (chi1 < chi2)
					{
						krankc2 = 1;
					}
					else
						krankc2 = -1;
				}
				TagB_tpl2->column("krankc2",krankc2);

				
				
				sigselection++;
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
				
				tnremaining = tnpip+tnpim+tnep+tnem+tnkp+tnkm+tnmup+tnmum+tnpm+tnpp  ;
				nremaining = npip+npim+nep+nem+nkp+nkm+nmup+nmum+npm+npp ;
				
				TagB_tpl2->column("kstc0px" , (*i).child(0).child(0).px());
				TagB_tpl2->column("kstc0py" , (*i).child(0).child(0).py());
				TagB_tpl2->column("kstc0pz" , (*i).child(0).child(0).pz());
				TagB_tpl2->column("kstc0mass" , (*i).child(0).child(0).p().mag());
				
				TagB_tpl2->column("kstc1px" , (*i).child(0).child(1).px());
				TagB_tpl2->column("kstc1py" , (*i).child(0).child(1).py());
				TagB_tpl2->column("kstc1pz" , (*i).child(0).child(1).pz());
				TagB_tpl2->column("kstc1mass" , (*i).child(0).child(1).p().mag());
				
				TagB_tpl2->column("kstarpx" , (*i).child(0).px());
				TagB_tpl2->column("kstarpy" , (*i).child(0).py());
				TagB_tpl2->column("kstarpz" , (*i).child(0).pz());
				TagB_tpl2->column("kstarmass" , (*i).child(0).p().mag());
				
				TagB_tpl2->column( "nb_trackremaining"    , nremaining         );//new
				TagB_tpl2->column( "tnb_trackremaining"    , tnremaining         );//new
				TagB_tpl2->column( "tnb_pi0remaining"    , npi0         );//new
				TagB_tpl2->column( "tnb_gamma"    , n_tgamma1         );//new
				TagB_tpl2->column( "nbmulti2"    , nbmulti2        );//new
				TagB_tpl2->column( "nb_pi015"    , npi015         );//new
				
				TagB_tpl2->column( "n_kl_tot"    , n_kl_tot         );//new
				TagB_tpl2->column( "n_kl"    , n_kl         );//new
				TagB_tpl2->column( "n_ks"    , n_ks         );//new
				
				//				double sdr,sdz;
				//				const Mdst_charged* sch1 = &(*i).child(0).mdstCharged();
				//				GetImpactParameters(sch1,&sdr,&sdz,3);
				//				TagB_tpl2->column( "ch0dr"     , sdr         );//new
				//				TagB_tpl2->column( "ch0dz"    , sdz         );//new
				
				double sum_ecl = 0;
				for(std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s !=eclaux_mag.end();s++ )
				{
					//               std::cout<<"start eclaux"<<std::endl;
					Mdst_ecl_aux& ch = *s;
					int flag = 0;
					int nfinal = Btag.relation().nFinalStateParticles(); // the number of final particle used to tagB
					int signfinal = (*i).relation().nFinalStateParticles();
					
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
					
					
					for(int kk = 0 ; kk < signfinal ; kk ++)
					{
						Particle sigchild = (*i).relation().finalStateParticle(kk);
						if ( sigchild.mdstCharged() )
						{
							if ( sigchild.mdstCharged().trk() == ch.trk())
							{
								flag=1;
							}
						}
						if (sigchild.mdstGamma())
						{
							if (sigchild.mdstGamma().ecl().get_ID()==ch.get_ID())
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
				TagB_tpl2->column("RatioToOtherType" , RatioToOtherType);
				TagB_tpl2->column("RaatioToSecondBest" , RaatioToSecondBest);
				//////////////////////K information///////////////////////
				//	 TagB_tpl2->column("kpx" , (*i).child(0).px());
				//         TagB_tpl2->column("kpy" , (*i).child(0).py());
				//         TagB_tpl2->column("kpz" , (*i).child(0).pz());
				//         TagB_tpl2->column("kmass" , (*i).child(0).p().mag());
				TagB_tpl2->column("sumecl",sum_ecl);
				//	 TagB_tpl2->column("kE" , (*i).child(0).e());
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
				TagB_tpl2->column( "reenergy"    , reenergy          );
				
				////////////////////////Btag///////////////////////////////////////////////////////////
				// Event ID
				TagB_tpl2->column( "ebeam"    , ebeam          );
				TagB_tpl2->column( "Tag_RunID", Run            );
				TagB_tpl2->column( "Tag_EvtID", Evt            );
				TagB_tpl2->column( "Tag_ExpID", Exp            );
				// TagB idhep
				TagB_tpl2->column( "tagid"    , frecB.tag_id() );
				// Kinetic variable
				TagB_tpl2->column( "tagBmbc"  , frecB.Mbc()    );
				TagB_tpl2->column( "tagBde"   , frecB.DeltaE() );
				float BtagE = frecB.DeltaE() + ebeam ;
				TagB_tpl2->column( "BtagE"  , BtagE   );
				// unique code for reconstructed D decay channels
				// for B -> D X decay modes only 1 entry is filled
				// for B -> D* ( -> D X) Y decay modes only 1 (D*) and 2 (D) entries are filled
				// for B -> D* ( -> D1 X) D2 decay modes only 1 (D*), 2 (D1) and 3 (D2) entries are filled
				// for B -> D*1 ( -> D1 X) D*2 ( -> D2 Y) decay modes all four entries are filled: D*1, D1, D*2, D2
				TagB_tpl2->column( "decmod"   , frecB.decay()       );
				TagB_tpl2->column( "ddecmod1" , frecB.Ddec(0)       );
				TagB_tpl2->column( "ddecmod2" , frecB.Ddec(1)       );
				TagB_tpl2->column( "ddecmod3" , frecB.Ddec(2)       );
				TagB_tpl2->column( "ddecmod4" , frecB.Ddec(3)       );
				TagB_tpl2->column( "mcinfo"   , frecB.MCinfo()      );
				TagB_tpl2->column( "nFS"      , frecB.nFS()         );
				// default selection (without continuum suppression applied)
				TagB_tpl2->column( "best"     , frecB.NBRank()      );
				TagB_tpl2->column( "nbout"    , frecB.NBout()       );
				// selection with continuum suppression
				TagB_tpl2->column( "csbest"   , frecB.cont_NBRank() );
				TagB_tpl2->column( "csnbout"  , frecB.cont_NBout()  );
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
				HepLorentzVector qlab = ( Bsigplab - k_plab); // signal side miss
				double Bsigmass2 = Bsigplab.mag()*Bsigplab.mag();//sigbmass^2
				
				double lab_total_miss = pmag(tmissinglab);
				double lab_sig_miss = pmag(qlab);
				double tmmiss2lab = tmissinglab.mag()*tmissinglab.mag();
				double qsquare = qlab.mag()*qlab.mag();
				double sblab = qsquare / Bsigmass2;
				double lab_k_p = pmag(k_plab);
				//////////////////////        total miss information        ////////////////////
				TagB_tpl2->column("tmmisslab",tmissinglab.mag());
				TagB_tpl2->column("tmmiss2lab",tmmiss2lab);
				TagB_tpl2->column("tpxmisslab",tmissinglab.px());
				TagB_tpl2->column("tpymisslab",tmissinglab.py());
				TagB_tpl2->column("tpzmisslab",tmissinglab.pz());
				TagB_tpl2->column("temisslab",tmissinglab.e());
				TagB_tpl2->column("tpmisslab",lab_total_miss);
				//////////////////////        signal miss information        ////////////////////
				TagB_tpl2->column("sblab",sblab);
				TagB_tpl2->column("qsquare",qsquare);
				TagB_tpl2->column("qpxlab",qlab.px());
				TagB_tpl2->column("qpylab",qlab.py());
				TagB_tpl2->column("qpzlab",qlab.pz());
				TagB_tpl2->column("qelab",qlab.e());
				TagB_tpl2->column("spmisslab",lab_sig_miss);
				//////////////////////        K information        ////////////////////
				
				TagB_tpl2->column("kmalab",k_plab.mag());
				TagB_tpl2->column("kpxlab",k_plab.px());
				TagB_tpl2->column("kpylab",k_plab.py());
				TagB_tpl2->column("kpzlab",k_plab.pz());
				TagB_tpl2->column("kelab",k_plab.e());
				TagB_tpl2->column("kplab",lab_k_p);
				//std::cout<<"p1"<<std::endl;
				//Lab frame miss cos
				double tcosmisslabinz;
				double tangelmisslabinz;
				tcosmisslabinz = ((tmissinglab.pz())/lab_total_miss);
				tangelmisslabinz = acos((tmissinglab.pz())/lab_total_miss);
				TagB_tpl2->column("tcosmizl",tcosmisslabinz);
				//h(*) B+ cos lab
				HepLorentzVector Btagplab (Btag.px(),Btag.py(),Btag.pz(),Btag.e());
				HepLorentzVector all_BEAM_polar ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle),E_HER+E_LER);
				HepLorentzVector P_BEAM_polar = (all_BEAM_polar - Btagplab);//BEAM - BTAG = SINGAL
				double CM_2 = pmag(P_BEAM_polar);
				double kp2 = pmag((*i).child(0));
				double cosksblab = costheta((*i).child(0),P_BEAM_polar);
				double cosklab = (*i).child(0).pz()/kp2;
				TagB_tpl2->column( "cosksblab"   , cosksblab);//new
				TagB_tpl2->column( "cosklab"   , cosklab);//new
				
				///////////////////////////////////////////////////////////////////////////////////////
				//                            Boost to upsilon 4S rest frame                         //
				///////////////////////////////////////////////////////////////////////////////////////
				
				tagB.boost( boost_vector.boostVector() );//boost tagb
				double CM_Btag_p = sqrt( tagB.px()*tagB.px() + tagB.py()*tagB.py() + tagB.pz()*tagB.pz() );
				double CM_Btag_pz = sqrt( tagB.pz()*tagB.pz() );
				//	 double CM_Btag_ppz = tagB.pz()*tagB.px() + tagB.pz()*tagB.py() + tagB.pz()*tagB.pz() ;
				double cosbcm = tagB.pz()/(CM_Btag_p);
				TagB_tpl2->column( "cosbcm"     , cosbcm );
				TagB_tpl2->column( "tagBpcm"    , CM_Btag_p );
				TagB_tpl2->column( "tagbPXcm"   , tagB.px());
				TagB_tpl2->column( "tagbPYcm"   , tagB.py());
				TagB_tpl2->column( "tagbPZcm"   , tagB.pz());
				HepLorentzVector BEAM(E_HER*sin(cross_angle), 0.0, E_HER*cos(cross_angle)-E_LER, E_HER+E_LER);
				HepLorentzVector k_p((*i).child(0).px(),(*i).child(0).py(),(*i).child(0).pz(),(*i).child(0).e());
				HepLorentzVector Bsigpcm (-tagB.px(),-tagB.py(),-tagB.pz(),tagB.e());
				BEAM.boost(boost_vector.boostVector());//boost beam
				k_p.boost(boost_vector.boostVector() );//boost signal k
				HepLorentzVector totalmissincm = ( BEAM - tagB - k_p); //total miss
				HepLorentzVector sigmisscm = (  Bsigpcm - k_p ); //signal miss
				double cm_totalmiss_p =pmag(totalmissincm);
				double bsigpcm =pmag(Bsigpcm) ;
				double sigmisspcm = pmag(sigmisscm);
				double anglebsigandsigmiss = acos((Bsigpcm.px()*sigmisscm.px() + Bsigpcm.py()*sigmisscm.py() + Bsigpcm.pz()*sigmisscm.pz())/(bsigpcm*sigmisspcm));//new
				double cosbsigandsigmiss = (Bsigpcm.px()*sigmisscm.px() + Bsigpcm.py()*sigmisscm.py() + Bsigpcm.pz()*sigmisscm.pz())/(bsigpcm*sigmisspcm);//new
				
				TagB_tpl2->column("angbssmi", anglebsigandsigmiss);//new
				TagB_tpl2->column("cosbssmi",cosbsigandsigmiss);//new
				/////////////////////        signal miss information        ////////////////////
				TagB_tpl2->column("sigmissp",sigmisscm.mag());
				double sigmmiss2 = sigmisscm.mag()*sigmisscm.mag();
				TagB_tpl2->column("Bsigmas2",Bsigmass2);
				double sbcm = sigmmiss2 /Bsigmass2  ;
				TagB_tpl2->column("sbcm",sbcm);
				TagB_tpl2->column("simmis2CM",sigmmiss2);
				TagB_tpl2->column("spxmisisC",sigmisscm.px());
				TagB_tpl2->column("spymissC",sigmisscm.py());
				TagB_tpl2->column("spzmissC",sigmisscm.pz());
				TagB_tpl2->column("semissC",sigmisscm.e());
				TagB_tpl2->column("spmissC",sigmisspcm);
				//////////////////////        K information        ////////////////////
				
				TagB_tpl2->column("mkcm",k_p.mag());
				TagB_tpl2->column("kpxCM",k_p.px());
				TagB_tpl2->column("kpyCM",k_p.py());
				TagB_tpl2->column("kpzCM",k_p.pz());
				TagB_tpl2->column("keCM",k_p.e());
				double CM_k_p = pmag(k_p);
				TagB_tpl2->column("kpCM",CM_k_p);
				//////////////////////        total miss information        ////////////////////
				TagB_tpl2->column("tmmisscm",totalmissincm.mag());
				TagB_tpl2->column("tmmiss2cm",totalmissincm.mag()*totalmissincm.mag());
				TagB_tpl2->column("tpxmisscm",totalmissincm.px());
				TagB_tpl2->column("tpymisscm",totalmissincm.py());
				TagB_tpl2->column("tpzmisscm",totalmissincm.pz());
				TagB_tpl2->column("temisscm",totalmissincm.e());
				TagB_tpl2->column("tpmisscm",cm_totalmiss_p);
				
				///////////////////////////////////////////////////////////////////////////////////////
				//                            Boost to upsilon B rest frame                         //
				///////////////////////////////////////////////////////////////////////////////////////
				HepLorentzVector sigmissb(sigmisscm);
				HepLorentzVector k_b(k_p);
				HepLorentzVector totalmissinb(totalmissincm);
				
				sigmissb.boost(-Bsigplab.boostVector());//signal miss
				totalmissinb.boost(-Bsigplab.boostVector());//total miss
				
				k_b.boost(-Bsigplab.boostVector());//signl K
				double sigmissm2b = sigmissb.mag()*sigmissb.mag();
				double sigmisspb = pmag(sigmissb);
				/////////////////////        signal miss information        ////////////////////
				TagB_tpl2->column("smmissb",sigmissb.mag());
				TagB_tpl2->column("smmiss2b",sigmissm2b);
				TagB_tpl2->column("spxmissb",sigmissb.px());
				TagB_tpl2->column("spymissb",sigmissb.py());
				TagB_tpl2->column("spzmissb",sigmissb.pz());
				TagB_tpl2->column("semissb",sigmissb.e());
				TagB_tpl2->column("spmissb",sigmisspb);
				double b_k_p = pmag(k_b);
				//////////////////////        K information        ////////////////////
				TagB_tpl2->column("kpb", b_k_p);
				TagB_tpl2->column("kpxb",k_b.px());
				TagB_tpl2->column("kpyb",k_b.py());
				TagB_tpl2->column("kpzb",k_b.pz());
				TagB_tpl2->column("keb",k_b.e());
				//////////////////////        total miss information        ////////////////////
				TagB_tpl2->column("tmmissb",totalmissinb.mag());
				TagB_tpl2->column("tmmiss2b",totalmissinb.mag()*totalmissinb.mag());
				TagB_tpl2->column("tpxmissb",totalmissinb.px());
				TagB_tpl2->column("tpymissb",totalmissinb.py());
				TagB_tpl2->column("tpzmissb",totalmissinb.pz());
				TagB_tpl2->column("temissb",totalmissinb.e());
				TagB_tpl2->column("tpmissb",pmag(totalmissinb));
				
				//////////////////////////////////dz///////////////////////////////////////////////////////
				HepPoint3D Btag_init;
				HepPoint3D bvertex;
				HepPoint3D bvertexsig;
				float vchisq_tag;
				int DOF_tag;
				float pvalue_tag;
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
				if(flagch)
				{
					kvertexfitter kvf;
					Mdst_charged bfinalnumber[nfinal];
					ExKFitterParticle btagf[nfinal];
					Particle bfinal[nfinal];
					for(int k2 = 0 ; k2 < nfinal; k2++)
					{
						int expid = 0;
						Particle bfinal_tmp = Btag.relation().finalStateParticle(k2);
						bfinal[k2]=bfinal_tmp;
						int nunp = -1 ;
						if(!(bfinal[k2].mdstCharged()))
							continue;
						addTrack2fit(kvf,  bfinal[k2]);
						
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
					//  std::cout<<"start k6.9"<<std::endl;
				}
				double dr,dz;
				int signfinal = (*i).relation().nFinalStateParticles(); // the number of final particle used to tagB
				Mdst_charged sigbfinalnumber[signfinal];
				ExKFitterParticle bsisgf[signfinal];
				Particle bsigfinal[signfinal];
				kvertexfitter kvfsig;
				for(int k3 = 0 ; k3 < signfinal; k3++)
				{
					int expid = 0;
					Particle bsigfinal_tmp = (*i).relation().finalStateParticle(k3);
					bsigfinal[k3]=bsigfinal_tmp;
					int nunp = -1 ;
					if(!(bsigfinal[k3].mdstCharged()))
						continue;
					addTrack2fit(kvfsig,  bsigfinal[k3] );
					
				}
				
				//				const Mdst_charged* ch1 = &(*i).child(0).mdstCharged();
				//				Mdst_charged Mdst_k = (*i).child(0).mdstCharged();
				//				Particle sig_k = (*i).child(0);
				
				//				addTrack2fit(kvfsig,sig_k);
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
				TagB_tpl2->column("DistTOOtherBdz",DistTOOtherBdz);
				TagB_tpl2->column("kdz",bvertexsig.z());
				TagB_tpl2->column("kdr", dr );
				TagB_tpl2->column("vchisq_sig",vchisq_sig);
				TagB_tpl2->column("DOF_sig",DOF_sig);
				TagB_tpl2->column("pvalue_sig",pvalue_sig);
				TagB_tpl2->column("vchisq_tag",vchisq_tag);
				TagB_tpl2->column("DOF_tag",DOF_tag);
				TagB_tpl2->column("pvalue_tag",pvalue_tag);
				double Distsign = pvalue_sig -pvalue_tag;
				TagB_tpl2->column("Distsign",Distsign);
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
							if ((*p1).theta() < 0.547 && (*p1).energy() > 0.1 && fabs(tmissinglab.angle(EclVirtualVec.vect()))<0.1 )
								sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
							//barrel e>0.05
							if ((*p1).theta() > 0.562 && (*p1).theta() < 2.246 && (*p1).energy()>0.05&& fabs(tmissinglab.angle(EclVirtualVec.vect()))< 0.1 )
								sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
							//backward e>0.15
							if ((*p1).theta() > 2.246 && (*p1).energy() > 0.15&& fabs(tmissinglab.angle(EclVirtualVec.vect()))< 0.1)
								sum_eclinmissdir = sum_eclinmissdir + (*p1).energy();
						}
						EclVirtualVec.setRho(1.0);EclVirtualVec.setTheta((*p1).theta());EclVirtualVec.setPhi((*p1).phi());
					}
				}//ecl_aux end
				TagB_tpl2->column("sum_eclinmissdir",sum_eclinmissdir);//new
				/////////////////////////////////////////////////////////////////////////////
				//                           k_sfw variables                               //
				////////////////////////////////////////////////////////////////////////////
				k_sfw ksfw_obj((*i));
				const double ksfw(ksfw_obj.fd());
				const int iksfw(ksfw_obj.i_mm2());
				TagB_tpl2->column("imm2",iksfw);
				double miss(ksfw_obj.mm2());
				TagB_tpl2->column("mmm2",miss);
				double et(ksfw_obj.e_t());
				TagB_tpl2->column("et",et);
				double H0oo(ksfw_obj.Hoo(0));
				double H1oo(ksfw_obj.Hoo(1));
				double H2oo(ksfw_obj.Hoo(2));
				double H3oo(ksfw_obj.Hoo(3));
				double H4oo(ksfw_obj.Hoo(4));
				TagB_tpl2->column("H0oo",H0oo);
				TagB_tpl2->column("H1oo",H1oo);
				TagB_tpl2->column("H2oo",H2oo);
				TagB_tpl2->column("H3oo",H3oo);
				TagB_tpl2->column("H4oo",H4oo);
				double H0son(ksfw_obj.Hso_n(0));
				double H1son(ksfw_obj.Hso_n(1));
				double H2son(ksfw_obj.Hso_n(2));
				double H3son(ksfw_obj.Hso_n(3));
				double H4son(ksfw_obj.Hso_n(4));
				TagB_tpl2->column("H0son",H0son);
				TagB_tpl2->column("H1son",H1son);
				TagB_tpl2->column("H2son",H2son);
				TagB_tpl2->column("H3son",H3son);
				TagB_tpl2->column("H4son",H4son);
				double H0soc(ksfw_obj.Hso_c(0));
				double H1soc(ksfw_obj.Hso_c(1));
				double H2soc(ksfw_obj.Hso_c(2));
				double H3soc(ksfw_obj.Hso_c(3));
				double H4soc(ksfw_obj.Hso_c(4));
				TagB_tpl2->column("H0soc",H0soc);
				TagB_tpl2->column("H1soc",H1soc);
				TagB_tpl2->column("H2soc",H2soc);
				TagB_tpl2->column("H3soc",H3soc);
				TagB_tpl2->column("H4soc",H4soc);
				double H0som(ksfw_obj.Hso_m(0));
				double H1som(ksfw_obj.Hso_m(1));
				double H2som(ksfw_obj.Hso_m(2));
				double H3som(ksfw_obj.Hso_m(3));
				double H4som(ksfw_obj.Hso_m(4));
				TagB_tpl2->column("H0som",H0som);
				TagB_tpl2->column("H1som",H1som);
				TagB_tpl2->column("H2som",H2som);
				TagB_tpl2->column("H3som",H3som);
				TagB_tpl2->column("H4som",H4som);
				//for shape variables
				Vector4 OtherBVector;
				float R2, spher, cos_thr, cos_thp, par_sfw[13];
				
				int vertexflag;
				HepPoint3D overtex;
				
				///////////////////////////////////////////////////////////////////////////////////////
				//                      MC TRuth(I don't know it work or not)                        //
				///////////////////////////////////////////////////////////////////////////////////////
				
				if(MCstatus == 1)
				{
					int check = -1;
					int hhindex= -1;
					int kstcheck = -1;
					///////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////
					
					
					if (frecB.tag_id() == -521)
					{
						int bnfinal = (*i).child(0).relation().nFinalStateParticles();
						for(int k = 0 ; k < bnfinal; k++)
						{
							Particle BB = (*i);
							Particle bchild = (*i).child(0).relation().finalStateParticle(k);
							if ( bchild.mdstCharged() )
							{
								const Mdst_charged c11h000 = bchild.mdstCharged();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( (Evt_decay.idhep() == 323  && Evt_current.idhep() == 521))
									{
										check = 1;
										break;
									}
									else
										check = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							else if ( bchild.mdstGamma() )
							{
								const Mdst_gamma c11h000 = bchild.mdstGamma();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( (Evt_decay.idhep() == 323  && Evt_current.idhep() == 521))
									{
										check = 1;
										break;
									}
									else
										check = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							if (check != 1)
								break;
						}//kend
						TagB_tpl2->column("check" ,check);
						int doublecheck = -1 ;
						for (std::vector<Gen_hepevt>::iterator it3 = gen_mgr.begin();it3 != gen_mgr.end(); it3++)
						{
							if((*it3).idhep()== 521 )
							{
								Gen_hepevt d = gen_mgr[(*it3).da(0)-1];
								//std::cout << "Evt_dd =" << d.idhep() << std::endl;
								if ( d.idhep() == 323)
								{
									for(int k = 0 ; k < bnfinal; k++)
									{
										doublecheck = 0 ;
										Particle bchild = (*i).child(0).relation().finalStateParticle(k);
										if ( bchild.mdstCharged() )
										{
											const Mdst_charged dsub = bchild.mdstCharged();
											Gen_hepevt  Evt_Ddecay= get_hepevt(dsub);
											toptodown(Evt_Ddecay,d,doublecheck,gen_mgr);
										}
										else if ( bchild.mdstGamma() )
										{
											const Mdst_gamma dsub = bchild.mdstGamma();
											Gen_hepevt Evt_Ddecay= get_hepevt(dsub);
											toptodown(Evt_Ddecay,d,doublecheck,gen_mgr);
										}
									}
								}
								else
									break;
							}
						}
						if ( check == 1 && doublecheck == 1 )
							hhindex = 1;
						else
							hhindex =0;
						TagB_tpl2->column("doublecheck" ,doublecheck);
						TagB_tpl2->column("hhindex" ,hhindex);
						
					}
					else if (frecB.tag_id() == 521)
					{
						int bnfinal = (*i).child(0).relation().nFinalStateParticles();
						for(int k = 0 ; k < bnfinal; k++)
						{
							Particle BB = (*i);
							Particle bchild = (*i).child(0).relation().finalStateParticle(k);
							if ( bchild.mdstCharged() )
							{
								const Mdst_charged c11h000 = bchild.mdstCharged();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( (Evt_decay.idhep() == -323  && Evt_current.idhep() == -521))
									{
										check = 1;
										break;
									}
									else
										check = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							else if ( bchild.mdstGamma() )
							{
								const Mdst_gamma c11h000 = bchild.mdstGamma();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( (Evt_decay.idhep() == -323  && Evt_current.idhep() == -521))
									{
										check = 1;
										break;
									}
									else
										check = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
						}//kend
						TagB_tpl2->column("check" ,check);
						int doublecheck = -1 ;
						for (std::vector<Gen_hepevt>::iterator it3 = gen_mgr.begin();it3 != gen_mgr.end(); it3++)
						{
							if((*it3).idhep()== -521 )
							{
								Gen_hepevt d = gen_mgr[(*it3).da(0)-1];
								//std::cout << "Evt_dd =" << d.idhep() << std::endl;
								if ( d.idhep() == -323)
								{
									for(int k = 0 ; k < bnfinal; k++)
									{
										doublecheck = 0 ;
										Particle bchild = (*i).child(0).relation().finalStateParticle(k);
										if ( bchild.mdstCharged() )
										{
											const Mdst_charged dsub = bchild.mdstCharged();
											Gen_hepevt  Evt_Ddecay= get_hepevt(dsub);
											toptodown(Evt_Ddecay,d,doublecheck,gen_mgr);
										}
										else if ( bchild.mdstGamma() )
										{
											const Mdst_gamma dsub = bchild.mdstGamma();
											Gen_hepevt Evt_Ddecay= get_hepevt(dsub);
											toptodown(Evt_Ddecay,d,doublecheck,gen_mgr);
										}
									}
								}
								else
									break;
							}
						}
						if ( check == 1 && doublecheck == 1 )
							hhindex = 2;
						else
							hhindex =0;
						TagB_tpl2->column("doublecheck" ,doublecheck);
						TagB_tpl2->column("hhindex" ,hhindex);
						
					}
					///////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////
					int bnfinal = (*i).child(0).relation().nFinalStateParticles();
					for(int k = 0 ; k < bnfinal; k++)
					{
						Particle BB = (*i);
						Particle bchild = (*i).child(0).relation().finalStateParticle(k);
						if ( bchild.mdstCharged() )
						{
							const Mdst_charged c11h000 = bchild.mdstCharged();
							Gen_hepevt Evt_decay= get_hepevt(c11h000);
							while (Evt_decay.mo(0))
							{
								Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
								if ( Evt_decay.idhep() == 323)
								{
									kstcheck = 1;
									break;
								}
								else
									kstcheck = 0;
								
								Evt_decay = Evt_current;
							}
							
						}
						else if ( bchild.mdstGamma() )
						{
							const Mdst_gamma c11h000 = bchild.mdstGamma();
							Gen_hepevt Evt_decay= get_hepevt(c11h000);
							while (Evt_decay.mo(0))
							{
								Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
								if ( Evt_decay.idhep()== 323 )
								{
									kstcheck = 1;
									break;
								}
								else
									kstcheck = 0;
								
								Evt_decay = Evt_current;
							}
							
						}
						if (kstcheck != 1)
							break;
					}//kend
					
					if (kstcheck != 1)
					{
						for(int k = 0 ; k < bnfinal; k++)
						{
							Particle BB = (*i);
							Particle bchild = (*i).child(0).relation().finalStateParticle(k);
							if ( bchild.mdstCharged() )
							{
								const Mdst_charged c11h000 = bchild.mdstCharged();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( Evt_decay.idhep() == -323)
									{
										kstcheck = -1;
										break;
									}
									else
										kstcheck = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							else if ( bchild.mdstGamma() )
							{
								const Mdst_gamma c11h000 = bchild.mdstGamma();
								Gen_hepevt Evt_decay= get_hepevt(c11h000);
								while (Evt_decay.mo(0))
								{
									Gen_hepevt Evt_current = gen_mgr[Evt_decay.mo(0)-1];
									if ( Evt_decay.idhep()== -323 )
									{
										kstcheck = -1;
										break;
									}
									else
										kstcheck = 0;
									
									Evt_decay = Evt_current;
								}
								
							}
							if (kstcheck != -1)
								break;
						}//kend
					}
					
					TagB_tpl2->column("kstcheck" ,kstcheck);
					
					
					
					///////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////
					
					if (mode == 1) // k+ pi0
					{
						const Mdst_charged ch0 = (*i).child(0).child(0).mdstCharged();
						const Mdst_gamma ch10 = (*i).child(0).child(1).child(0).mdstGamma();
						const Mdst_gamma ch11 = (*i).child(0).child(1).child(1).mdstGamma();
						
						Gen_hepevt Evtch0 = get_hepevt(ch0);
						Gen_hepevt Evtch10 = get_hepevt(ch10);
						Gen_hepevt Evtch11 = get_hepevt(ch11);
						
						Gen_hepevt tEvtch0;
						Gen_hepevt tEvtch10;
						Gen_hepevt tEvtch11;
						std::cout<< "debug0 ch0id = " << Evtch0.idhep() <<std::endl;
						std::cout<< "debug0 ch10id = " << Evtch10.idhep() <<std::endl;
						std::cout<< "debug0 ch11id = " << Evtch11.idhep() <<std::endl;
						
						//						while( abs(Evtch0.idhep()) != 321 && Evtch0.mo(0))
						//						{
						//							tEvtch0 = gen_mgr[Evtch0.mo(0)-1];
						//							std::cout<< "debug1 ch0id = " << tEvtch0.idhep() <<std::endl;
						//							Evtch0 = tEvtch0 ;
						//						}
						/*						while( abs(Evtch10.idhep()) != 22  && Evtch10.mo(0))
						 {
							tEvtch10 = gen_mgr[Evtch10.mo(0)-1];
							std::cout<< "debug1 ch10id = " << tEvtch10.idhep() <<std::endl;
							Evtch10 = tEvtch10 ;
						 }
						 while( abs(Evtch11.isthep()) != 22  && Evtch11.mo(0))
						 {
							tEvtch11 = gen_mgr[Evtch11.mo(0)-1];
							std::cout<< "debug1 ch11id = " << tEvtch11.idhep() <<std::endl;
							Evtch11 = tEvtch11 ;
						 }
						 */
						TagB_tpl2->column("hi0",Evtch0.idhep());
						TagB_tpl2->column("hi10",Evtch10.idhep());
						TagB_tpl2->column("hi11",Evtch11.idhep());
						
						Gen_hepevt EvtP0;
						Gen_hepevt EvtGP0;
						Gen_hepevt EvtGGP0;
						
						Gen_hepevt EvtP10;
						Gen_hepevt EvtGP10;
						Gen_hepevt EvtGGP10;
						
						Gen_hepevt EvtP11;
						Gen_hepevt EvtGP11;
						Gen_hepevt EvtGGP11;
						
						if(Evtch0.mo(0))
						{
							EvtP0 = gen_mgr[Evtch0.mo(0)-1];
							TagB_tpl2->column("mo0",EvtP0.idhep());
							
							if(EvtP0.mo(0))
							{
								EvtGP0 = gen_mgr[EvtP0.mo(0)-1];
								TagB_tpl2->column("gmo0",EvtGP0.idhep());
								
								
								if(EvtGP0.mo(0))
								{
									EvtGGP0 = gen_mgr[EvtGP0.mo(0)-1];
									TagB_tpl2->column("ggmo0",EvtGGP0.idhep());
								}
							}
						}
						if(Evtch10.mo(0))
						{
							EvtP10 = gen_mgr[Evtch10.mo(0)-1];
							TagB_tpl2->column("mo10",EvtP10.idhep());
							
							if(EvtP10.mo(0))
							{
								EvtGP10 = gen_mgr[EvtP10.mo(0)-1];
								TagB_tpl2->column("gmo10",EvtGP10.idhep());
								
								
								if(EvtGP10.mo(0))
								{
									EvtGGP10 = gen_mgr[EvtGP10.mo(0)-1];
									TagB_tpl2->column("ggmo10",EvtGGP10.idhep());
								}
							}
						}
						if(Evtch11.mo(0))
						{
							EvtP11 = gen_mgr[Evtch11.mo(0)-1];
							TagB_tpl2->column("mo11",EvtP11.idhep());
							
							if(EvtP11.mo(0))
							{
								EvtGP11 = gen_mgr[EvtP11.mo(0)-1];
								TagB_tpl2->column("gmo11",EvtGP11.idhep());
								
								
								if(EvtGP11.mo(0))
								{
									EvtGGP11 = gen_mgr[EvtGP11.mo(0)-1];
									TagB_tpl2->column("ggmo11",EvtGGP11.idhep());
								}
							}
						}
						int hindex = -1 ;
						
						if( Evtch0.idhep()== 321 && Evtch10.idhep()==22 && Evtch11.idhep()==22 && EvtP0.idhep()==323 && EvtP10.idhep()==111 && EvtP11.idhep()==111 && EvtGP0.idhep()==521 && EvtGP10.idhep()==323 && EvtGP11.idhep()==323 && EvtGGP10.idhep()==521 && EvtGGP11.idhep()==521 )
						{
							hindex = 1;
						}
						else if( Evtch0.idhep()== -321 && Evtch10.idhep()==22 && Evtch11.idhep()==22 && EvtP0.idhep()==-323 && EvtP10.idhep()==111 && EvtP11.idhep()==111 && EvtGP0.idhep()==-521 && EvtGP10.idhep()==-323 && EvtGP11.idhep()==-323 && EvtGGP10.idhep()==-521 && EvtGGP11.idhep()==-521 )
						{
							hindex = 2;
						}
						else
						{
							hindex = 0;
						}
						
						TagB_tpl2->column("hindex",hindex);
						
					}
					else if (mode == 2) // ks pi+
					{
						std::cout<< "4" << std::endl;
						std::cout<< "id0 =" << (*i).child(0).child(1).lund() <<std::endl;
						const Mdst_charged ch0 = (*i).child(0).child(1).mdstCharged();
						
						std::cout<< "id1 =" << (*i).child(0).child(0).child(0).lund() <<std::endl;
						std::cout<< "id2 =" << (*i).child(0).child(0).child(1).lund()<<std::endl;
						
						const Mdst_charged ch10 = (*i).child(0).child(0).child(0).mdstCharged();
						const Mdst_charged ch11 = (*i).child(0).child(0).child(1).mdstCharged();
						std::cout<< "5" << std::endl;
						Gen_hepevt Evtch0 = get_hepevt(ch0);
						Gen_hepevt Evtch10 = get_hepevt(ch10);
						Gen_hepevt Evtch11 = get_hepevt(ch11);
						
						Gen_hepevt tEvtch0 ;
						Gen_hepevt tEvtch10;
						Gen_hepevt tEvtch11;
						
						/*						while( abs(Evtch0.idhep()) != 211 && Evtch0.mo(0) )
						 {
							tEvtch0 = gen_mgr[Evtch0.mo(0)-1];
							Evtch0 = tEvtch0 ;
						 }
						 while( abs(Evtch10.isthep()) != 211 && Evtch10.mo(0))
						 {
							tEvtch10 = gen_mgr[Evtch10.mo(0)-1];
							Evtch10 = tEvtch10 ;
						 }
						 while( abs(Evtch11.isthep()) != 211 && Evtch11.mo(0))
						 {
							tEvtch11 = gen_mgr[Evtch11.mo(0)-1];
							Evtch11 = tEvtch11 ;
						 }
						 */
						std::cout<< "6" << std::endl;
						TagB_tpl2->column("hi0",Evtch0.idhep());
						TagB_tpl2->column("hi10",Evtch10.idhep());
						TagB_tpl2->column("hi11",Evtch11.idhep());
						
						Gen_hepevt EvtP0;
						Gen_hepevt EvtGP0;
						Gen_hepevt EvtGGP0;
						
						Gen_hepevt EvtP10;
						Gen_hepevt EvtGP10;
						Gen_hepevt EvtGGP10;
						Gen_hepevt EvtGGGP10;
						
						Gen_hepevt EvtP11;
						Gen_hepevt EvtGP11;
						Gen_hepevt EvtGGP11;
						Gen_hepevt EvtGGGP11;
						std::cout<< "7" << std::endl;
						if(Evtch0.mo(0))
						{
							EvtP0 = gen_mgr[Evtch0.mo(0)-1];
							TagB_tpl2->column("mo0",EvtP0.idhep());
							
							if(EvtP0.mo(0))
							{
								EvtGP0 = gen_mgr[EvtP0.mo(0)-1];
								TagB_tpl2->column("gmo0",EvtGP0.idhep());
								
								
								if(EvtGP0.mo(0))
								{
									EvtGGP0 = gen_mgr[EvtGP0.mo(0)-1];
									TagB_tpl2->column("ggmo0",EvtGGP0.idhep());
									//								if(EvtGGP0.mo(0))
									//								{
									//									EvtGGGP0 = gen_mgr[EvtGGP0.mo(0)-1];
									//									TagB_tpl2->column("gggmo0",EvtGGGP0.idhep());
									//								}
								}
							}
						}
						std::cout<< "8" << std::endl;
						if(Evtch10.mo(0))
						{
							EvtP10 = gen_mgr[Evtch10.mo(0)-1];
							TagB_tpl2->column("mo10",EvtP10.idhep());
							
							if(EvtP10.mo(0))
							{
								EvtGP10 = gen_mgr[EvtP10.mo(0)-1];
								TagB_tpl2->column("gmo10",EvtGP10.idhep());
								
								
								if(EvtGP10.mo(0))
								{
									EvtGGP10 = gen_mgr[EvtGP10.mo(0)-1];
									TagB_tpl2->column("ggmo10",EvtGGP10.idhep());
									if(EvtGGP10.mo(0))
									{
										EvtGGGP10 = gen_mgr[EvtGGP10.mo(0)-1];
										TagB_tpl2->column("gggmo10",EvtGGGP10.idhep());
									}
								}
							}
						}
						std::cout<< "9" << std::endl;
						if(Evtch11.mo(0))
						{
							EvtP11 = gen_mgr[Evtch11.mo(0)-1];
							TagB_tpl2->column("mo11",EvtP11.idhep());
							
							if(EvtP11.mo(0))
							{
								EvtGP11 = gen_mgr[EvtP11.mo(0)-1];
								TagB_tpl2->column("gmo11",EvtGP11.idhep());
								
								
								if(EvtGP11.mo(0))
								{
									EvtGGP11 = gen_mgr[EvtGP11.mo(0)-1];
									TagB_tpl2->column("ggmo11",EvtGGP11.idhep());
									if(EvtGGP11.mo(0))
									{
										EvtGGGP11 = gen_mgr[EvtGGP11.mo(0)-1];
										TagB_tpl2->column("gggmo11",EvtGGGP11.idhep());
										
									}
								}
							}
						}
						int hindex = -1 ;
						std::cout<< "10" << std::endl;
						if( Evtch0.idhep()== 211 && Evtch10.idhep()==211 && Evtch11.idhep()==-211 && EvtP0.idhep()==323 && EvtP10.idhep()==310 && EvtP11.idhep()==310 && EvtGP0.idhep()==521 && abs(EvtGP10.idhep())==311 && abs(EvtGP11.idhep())==311 && EvtGGP10.idhep()==323 && EvtGGP11.idhep()==323 && EvtGGGP10.idhep()==521 && EvtGGGP11.idhep()==521  )
						{
							hindex = 1;
						}
						else if( Evtch0.idhep()== -211 && Evtch10.idhep()==211 && Evtch11.idhep()==-211 && EvtP0.idhep()==-323 && EvtP10.idhep()==310 && EvtP11.idhep()==310 && EvtGP0.idhep()==-521 && abs(EvtGP10.idhep())==311 && abs(EvtGP11.idhep())==311 && EvtGGP10.idhep()==-323 && EvtGGP11.idhep()==-323 && EvtGGGP10.idhep()==-521 && EvtGGGP11.idhep()==-521 )
						{
							hindex = 2;
						}
						else
						{
							hindex = 0;
						}
						
						TagB_tpl2->column("hindex",hindex);
						
					}
				}//truth event  end
				//        std::cout<<"start truth event"<<std::endl;
				
				/////////////////////////////////trace form top to down//////////////////////////////////
				int nbpd=0,nbnd=0,topcheck= 0;
				for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin();i != gen_mgr.end(); i++)
				{
					if (abs((*i).idhep())==521)
					{
						int da1=(*i).da(0),da2=(*i).da(1);
						
						if ((*i).idhep() == 521)
							TagB_tpl2->column("charged",1);
						else if ((*i).idhep() == -521)
							TagB_tpl2->column("charged",-1);
						else
							TagB_tpl2->column("charged",0);
						
						nbpd=(da2-da1+1);
						TagB_tpl2->column("nbpd",nbpd);
						Gen_hepevt kstar=gen_mgr[(*i).da(0)-1];
						Gen_hepevt nu=gen_mgr[(*i).da(0)];
						Gen_hepevt nubar=gen_mgr[(*i).da(0)+1];
						
						if ((*i).idhep() == 521)
						{
							for(int start=0;start<(da2-da1+1);start++)
							{
								Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
								char bpdnumber[32];
								sprintf (bpdnumber,"%s%d","bpd",start+1);
								TagB_tpl2->column(bpdnumber,Evda.idhep());
							}
							
							if (kstar.idhep() == 323 && nu.idhep() == 18 && nubar.idhep() == -18)
								TagB_tpl2->column("topcheck",1);
							
						}
						else if ((*i).idhep() == -521)
						{
							for(int start=0;start<(da2-da1+1);start++)
							{
								Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
								char bndnumber[32];
								sprintf (bndnumber,"%s%d","bnd",start+1);//bnd1,bnd2,bnd3....
								TagB_tpl2->column(bndnumber,Evda.idhep());
							}
							
							if (kstar.idhep() == -323 && nu.idhep() == 18 && nubar.idhep() == -18)
								TagB_tpl2->column("topcheck",-1);
							
						}
						
					}
					
				}
				
				TagB_tpl2->dumpData();
				*status = 1;
			}//k_plus end/////////////////////////////////////////////////////////////////////////
			
			
			
			
			
		//   std::cout <<"fullreconeff ="<< fullreconeff<<std::endl;
		//   std::cout <<"sigselection = "<< sigselection<<std::endl;
		


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
		}     //   std::cout<<"start shape3"<<std::endl;

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
				// std::cout<<"start shape4"<<std::endl;

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
				//  std::cout<<"start shape5"<<std::endl;

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
			// std::cout<<"start shape6"<<std::endl;

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
		// std::cout<<"start shape9"<<std::endl;

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

			// std::cout<<"start shape10"<<std::endl;

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
			// std::cout<<"start shape11"<<std::endl;

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
		// std::cout<<"start shape12"<<std::endl;

	}
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
	double ana_fullreconk::costheta(Particle p1 , Particle p2)
	{
		return ( p1.px()*p2.px()+p1.py()*p2.py()+p1.pz()*p2.pz() ) / (pmag(p1)*pmag(p2));
	}
	double ana_fullreconk::pt(Particle p1)
	{
		return sqrt(p1.px()*p1.px()+p1.py()*p1.py());
	}
	double ana_fullreconk::pmag(Particle p1)
	{
		return sqrt(p1.px()*p1.px()+p1.py()*p1.py()+p1.pz()*p1.pz());
	}
	double ana_fullreconk::pmag(HepLorentzVector& p1)
	{
		return sqrt(p1.px()*p1.px()+p1.py()*p1.py()+p1.pz()*p1.pz());
	}
	double ana_fullreconk::costheta(HepLorentzVector& p1 , Particle p2)
	{
		return ( p1.px()*p2.px()+p1.py()*p2.py()+p1.pz()*p2.pz() ) / (pmag(p1)*pmag(p2));
	}
	double ana_fullreconk::costheta(Particle p1 , HepLorentzVector& p2)
	{
		return ( p1.px()*p2.px()+p1.py()*p2.py()+p1.pz()*p2.pz() ) / (pmag(p1)*pmag(p2));
	}
	double ana_fullreconk::costheta(HepLorentzVector& p1 , HepLorentzVector& p2)
	{
		return ( p1.px()*p2.px()+p1.py()*p2.py()+p1.pz()*p2.pz() ) / (pmag(p1)*pmag(p2));
	}

	int ana_fullreconk::toptodown( Gen_hepevt Evt_Ddecay , Gen_hepevt x_d , int &doublecheck ,std::vector<Gen_hepevt> &gen_mgr )
	{
		if ( x_d.da(0) )
		{
			std::cout << "toptodown"  << std::endl;
			Gen_hepevt x_ddais=gen_mgr[x_d.da(0)-1];
			
			if (doublecheck == 1)
			{
				std::cout << "------------"  << std::endl;
				std::cout<< " x_d.idhep()) ==  " <<x_d.idhep()<<std::endl;
				
			}
			else
			{
				int da0 = x_d.da(0), da1 = x_d.da(1);
				for(int start=0;start<(da1-da0+1);start++)
				{
					Gen_hepevt x_dda=gen_mgr[x_d.da(0)-1+start];
					if (Evt_Ddecay == x_dda && Evt_Ddecay.idhep() == x_dda.idhep() )
					{
						doublecheck = 1;
					}
					toptodown(Evt_Ddecay , x_dda , doublecheck , gen_mgr );
				}
			}
		}
		else
		{
			if (x_d.isthep() > 0)
			{
				std::cout<< " x_d.idhep()) ==  " <<x_d.idhep()<<std::endl;
				
			}
		}
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
		// std::cout<< "plsize1 = "<<pl.size() << std::endl;
		int plb = pl.size();
		remove_dup_trk(pl, DRCUT, DZCUT);
		int pla = pl.size();
		// std::cout<< "plsize2 = "<<pl.size() << std::endl;
		//		if(abs(pla - plb)!=0)
		// std::cout<< " has dup track "<< std::endl;

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

			// std::cout <<"sizebeforeinfun="<< dpi_plus.size()+dpi_minus.size()+de_plus.size()+de_minus.size()+dk_plus.size()+dk_minus.size()+dmu_plus.size()+dmu_minus.size()+dp_plus.size()+dp_minus.size() <<std::endl;


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
			// std::cout <<"sizeafterinfun="<< dpi_plus.size()+dpi_minus.size()+de_plus.size()+de_minus.size()+dk_plus.size()+dk_minus.size()+dmu_plus.size()+dmu_minus.size()+dp_plus.size()+dp_minus.size() <<std::endl;

		}
	}

	void ana_fullreconk::remove_dup_trk(std::vector<Particle> &plist, double DRCUT, double DZCUT){
		// remove duplicate tracks
		// for low pt tracks with similar momentum
		if (plist.size() < 1) return;
		// std::cout <<"remove_dup_trk0"<<std::endl;

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
			// std::cout <<"remove_dup_trk1"<<std::endl;
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
					// std::cout <<"remove_dup_trk2"<<std::endl;

					double dr1, dr2, dz1, dz2;
					GetImpactParameters(&(*i).mdstCharged(), &dr1, &dz1, mhyp1);
					GetImpactParameters(&(*j).mdstCharged(), &dr2, &dz2, mhyp2);
					double dist1 = pow(dr1/DRCUT, 2) + pow(dz1/DZCUT, 2);
					double dist2 = pow(dr2/DRCUT, 2) + pow(dz2/DZCUT, 2);
					if (dist1 >= dist2)
					{
						// std::cout <<"remove_dup_trk3"<<std::endl;

						--(i = plist.erase(i));
						break;
					}
					else
					{
						// std::cout <<"remove_dup_trk4"<<std::endl;

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
