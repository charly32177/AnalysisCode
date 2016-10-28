//
// $Id: Beamenergy.cc 10143 2007-05-25 23:12:57Z katayama $
//
// $Log$
// Revision 1.3  2005/04/13 13:13:19  matumot
// - add function to call benergy for evtgen (from Ishikawa-san)
//
// Revision 1.2  2004/05/28 06:37:06  matumot
// updated to use run-independent beam energy
// from database
//
//
//======================================================================
//  File Name :  Beamenergy.cc
//----------------------------------------------------------------------
//  First Version 2004.01.07 by T.Shibata
//
//  Modified by T.Shibata 2004.03.10 
//  Modified by T.Shibata 2004.05.17 ; add Run-independent value
//
//======================================================================
#include "belle.h"
#include <iostream>
#include "basfmodule/Module.h"
#include "basfmodule/Module_descr.h"
#include "belleutil/belutil.h"

#include "pntdb/TPntFDDB.h"
#include "pntdb/TPntDB.h"

#include "benergy/BeamEnergy.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//=====================================================================
int      BeamEnergy::defaut_run_flag     = 0;
double   BeamEnergy::defaut_e_beam_corr  = 5.290000;
double   BeamEnergy::default_e_beam_orig = 5.290000;
double   BeamEnergy::default_e_beam2     = 5.290000;
double   BeamEnergy::default_e_beam_err  = 0;
double   BeamEnergy::default_e_her       = 7.998213;
double   BeamEnergy::default_e_ler       = 3.499218;
double   BeamEnergy::default_corss_angle = 0.022;
double   BeamEnergy::default_e_her_orig  = 7.998213;
double   BeamEnergy::default_e_ler_orig  = 3.499218;

double   BeamEnergy::default_2_e_her     = 7.996500;
double   BeamEnergy::default_2_e_ler     = 3.500000; 

bool     BeamEnergy::m_MC                = false;
int      BeamEnergy::m_exp               = -1;
int      BeamEnergy::m_run               = -1;
int      BeamEnergy::m_version           =  1;
int      BeamEnergy::m_flag_run          =  0; 
int      BeamEnergy::m_flag_db           =  1;
bool     BeamEnergy::m_valid             = false;
bool     BeamEnergy::m_usable            = false;

int      BeamEnergy::m_run_flag          = -1;
double   BeamEnergy::m_e_beam_corr       = -999;
double   BeamEnergy::m_e_beam_orig       = -999;
double   BeamEnergy::m_e_beam2           = -999;
double   BeamEnergy::m_e_beam_err        = -999; 
double   BeamEnergy::m_e_her             = -999;
double   BeamEnergy::m_e_ler             = -999;
double   BeamEnergy::m_corss_angle       = -999;
double   BeamEnergy::m_e_her_orig        = -999; 
double   BeamEnergy::m_e_ler_orig        = -999;
 
int      BeamEnergy::mt_expno[100]       = {-1};
int      BeamEnergy::mt_runno[100]       = {-1};
int      BeamEnergy::mt_run_flag[100]    = {-1};
double   BeamEnergy::mt_e_beam_corr[100] = {-1};
double   BeamEnergy::mt_e_beam_orig[100] = {-1};
double   BeamEnergy::mt_e_beam2[100]     = {-1};
double   BeamEnergy::mt_e_beam_err[100]  = {-1};
double   BeamEnergy::mt_e_her[100]       = {-1};
double   BeamEnergy::mt_e_ler[100]       = {-1};
double   BeamEnergy::mt_corss_angle[100] = {-1};
double   BeamEnergy::mt_e_her_orig[100]  = {-1};
double   BeamEnergy::mt_e_ler_orig[100]  = {-1};

//--------------------------------------------------------------------
void BeamEnergy::define_global(Module_descr* dscr){

  dscr->define_param("beam_enegy_version",
                     "beam_enegy_version", &BeamEnergy::m_version);
  dscr->define_param("flag_rundep",
                     "flag_run_dependent", &BeamEnergy::m_flag_run);
  dscr->define_param("flag_db",
                     "flag_read_database", &BeamEnergy::m_flag_db);
  return;
};
//-------------------------------------------------------------------
void BeamEnergy::begin_run( const int ver, const int flag_run_dep, 
                             const int flag_db ){

  m_valid  = false;
  m_usable = false;
  m_exp = -1;
  m_run = -1; 

  m_version  = ver;
  m_flag_run = flag_run_dep;
  m_flag_db  = flag_db;

  Belle_runhead_Manager& evtmgr = Belle_runhead_Manager::get_manager();
  Belle_runhead_Manager::const_iterator belleevt = evtmgr.begin();

  if(belleevt == evtmgr.end()){
     dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Error: Cannot read \"Belle_RunHead\"." << std::endl;
     set_default_value();

     return; 
  }
 
  m_MC = false;
  if( belleevt->ExpMC()==2 )  m_MC = true;

  // #EXP & #RUN
  int exp = belleevt->ExpNo(); 
  int run = belleevt->RunNo(); 
  if( m_exp==exp && m_run==run && m_valid ){
      set_default_value();
      m_valid  = true;
      return ;
  }
  m_exp = exp;
  m_run = run; 

  /* 
   -------------------------------------------------------------
    m_flag_run == 0 or >3  : Run Dependent  
    m_flag_run == 1        : Default Value 
    m_flag_run == 2        : Second Default Value 
    m_flag_run == 3        : Run Independent Value( Mean Value )
   -------------------------------------------------------------
  */

  if(m_flag_run==1){  
    set_default_value();  
    m_valid  = true;
    m_usable = true;
 
    return; 
  }  

  if(m_flag_run==2){
    set_2nd_default_value();
    m_valid  = true;
    m_usable = true;

    return;
  }

  /*
   -------------------------------------------------------------
    Run number ==  0 ; Return Run-Independent Value  
                       Database Version = 2
    Run number !=  0 ; Return Run-dependent Value
                       Database Version = 1
   -------------------------------------------------------------
  */

 if( m_run == 0 || m_flag_run==3 ){
  m_version=2;

  if( m_exp == mt_expno[0] ){
      if( mt_run_flag[0]!=-1 ){
	  set_matrix_value(0);
	  return;  
      }
      if( mt_run_flag[0]==-1 ){
	  dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] There is no Beam Energy data." << std::endl;
	  set_default_value();
	  return;
      }
  }

  if(GetBeamEnergy3(m_version)){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] There is no Beam Energy data." << std::endl;
    set_default_value();
    return ; 
  }else{
    if( mt_run_flag[0]==-1){    
       dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] There is no Beam Energy data." << std::endl;
       set_default_value();
       return ; 
    }
       set_matrix_value(0);
       return ; 
  }   
 }

 if( m_run != 0 ){

  m_version=1;

  /*
   --------------------------------------------------------
    Read Database Flag
    m_flag_db  = 0 ; Read Database at every run
    m_flag_db >= 1 ; Read Database at every 100 run  
   --------------------------------------------------------
  */

  if( m_flag_db==0 ){
    if(GetBeamEnergy2(m_version)){
      dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] There is no Beam Energy data." << std::endl;
      set_default_value();
      return ; 

    }else{

      Beam_energy_Manager& Beam_mgr = Beam_energy_Manager::get_manager();
      Beam_energy_Manager::iterator beit = Beam_mgr.begin();

      if( beit == Beam_mgr.end()){
	dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] There is no Beam Energy data." << std::endl;
        set_default_value(); 
        return; 
      }

      m_run_flag    = beit->run_flag();
      m_e_beam_corr = beit->E_beam();
      m_e_beam_orig = beit->E_beam_orig();
      m_e_beam2     = beit->E_beam2();
      m_e_beam_err  = beit->E_beam_err();
      m_e_her       = beit->E_HER();
      m_e_ler       = beit->E_LER();
      m_corss_angle = beit->CROSS();
      m_e_her_orig  = beit->E_HER_orig();
      m_e_ler_orig  = beit->E_LER_orig();
      m_valid  = true;
      m_usable = true;

      return ; 
    }
  }
 

 if( m_flag_db>0 ){

  if( m_exp == mt_expno[0] ){
      for( int i=0; i< 100; i++ ){
  	if( m_run == mt_runno[i] && mt_run_flag[i]!=-1 ){
	  set_matrix_value(i);
	  return;  
        }
	if( m_run == mt_runno[i] && mt_run_flag[i]==-1 ){
	  dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] There is no Beam Energy data." << std::endl;
	  set_default_value();
	  return;
        }
      }
  }    

  if(GetBeamEnergy(m_version)){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] There is no Beam Energy data." << std::endl;
      set_default_value();
      return ; 

  }else{

   if( mt_run_flag[0]==-1){    
       dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] There is no Beam Energy data." << std::endl;
       set_default_value();
       return ; 
   }
       set_matrix_value(0);
       return ; 
  }

  }
 }  

}
//--------------------------------------------------------------------
void BeamEnergy::set_default_value(void){

      m_run_flag    = defaut_run_flag;
      m_e_beam_corr = defaut_e_beam_corr;
      m_e_beam_orig = default_e_beam_orig;
      m_e_beam2     = default_e_beam2;
      m_e_beam_err  = default_e_beam_err;
      m_e_her       = default_e_her;
      m_e_ler       = default_e_ler;
      m_corss_angle = default_corss_angle;
      m_e_her_orig  = default_e_her_orig;
      m_e_ler_orig  = default_e_ler_orig;

}
//--------------------------------------------------------------------
void BeamEnergy::set_2nd_default_value(void){

    m_run_flag    = defaut_run_flag;
    m_e_beam_corr = defaut_e_beam_corr;
    m_e_beam_orig = default_e_beam_orig;
    m_e_beam2     = default_e_beam2;
    m_e_beam_err  = default_e_beam_err;
    m_e_her       = default_2_e_her;
    m_e_ler       = default_2_e_ler;
    m_corss_angle = default_corss_angle;
    m_e_her_orig  = default_e_her_orig;
    m_e_ler_orig  = default_e_ler_orig;

}
//--------------------------------------------------------------------
void BeamEnergy::set_matrix_value(int i){

      m_run_flag    = mt_run_flag[i];
      m_e_beam_corr = mt_e_beam_corr[i];
      m_e_beam_orig = mt_e_beam_orig[i];
      m_e_beam2     = mt_e_beam2[i];
      m_e_beam_err  = mt_e_beam_err[i];
      m_e_her       = mt_e_her[i];
      m_e_ler       = mt_e_ler[i];
      m_corss_angle = mt_corss_angle[i];
      m_e_her_orig  = mt_e_her_orig[i];
      m_e_ler_orig  = mt_e_ler_orig[i];

      m_valid  = true;
      m_usable = true;

return;
}
//--------------------------------------------------------------------
int BeamEnergy::GetBeamEnergy( const int ver ){

 Beam_energy_Manager& Beam_mgr = Beam_energy_Manager::get_manager();   
 Beam_mgr.remove();

 TPntFDDB fddb("rif");
 TPntDB db(fddb, "benergy");
 if(!db.IsOK()){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Error : Cannot access rif::beam_profile.\n";
    return 1;
 }

 db.quiet_mode(1);
 int cnt_err(0);
 for( int i=0; i<100; i++ ){ 
   int flg(1);
   int rt = db.Get(m_exp, m_run+i , ver );
   if( rt<0 ) {
     flg = 0;
     cnt_err++;
   }
   input_matrix( flg, Beam_mgr, m_exp, m_run, i );   
 }
 dout(Debugout::WARN,"Beamenergy")
   << " No beam energy information in " 
   << cnt_err << " run(s) out of 100." << std::endl;

//  for( int i=0; i<100; i++ ){ 
//  if(db.Get(m_exp, m_run+i , ver ) <= 0){
//    input_matrix( 0, Beam_mgr, m_exp, m_run, i );   
//    continue;
//  }
//    input_matrix( 1, Beam_mgr, m_exp, m_run, i );
//  }
 
 db.non_clear(BEAM_ENERGY);
 
 return 0;
}
//--------------------------------------------------------------------
int BeamEnergy::GetBeamEnergy2( const int ver ){

  Beam_energy_Manager& Beam_mgr = Beam_energy_Manager::get_manager();   
  Beam_mgr.remove();

  TPntFDDB fddb("rif");
  TPntDB db(fddb, "benergy");
  if(!db.IsOK()){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Error : Cannot access rif::beam_profile.\n";
    return 1;
  }

  if(db.Get(m_exp, m_run , ver ) <= 0){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Error = : Cannot get data from database correctly.\n";
    return 1;
  }
  db.non_clear(BEAM_ENERGY);

  return 0;
}
//--------------------------------------------------------------------
int BeamEnergy::GetBeamEnergy3( const int ver ){

  Beam_energy_Manager& Beam_mgr = Beam_energy_Manager::get_manager();
  Beam_mgr.remove();

  TPntFDDB fddb("rif");
  TPntDB db(fddb, "benergy");
  if(!db.IsOK()){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Error : Cannot access rif::beam_profile.\n";
    return 1;
  }

  if(db.Get( m_exp , 1 , ver ) <= 0){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Error = : Cannot get data from database correctly.\n";
    return 1;
  }
  input_matrix( 1, Beam_mgr, m_exp, m_run, 0 );

  db.non_clear(BEAM_ENERGY);

  return 0;
}
//--------------------------------------------------------------------
// run by run beam energy by A.I
int BeamEnergy::GetBeamEnergy4( const int ver,
				const int expno,
				const int runno ){

  Beam_energy_Manager& Beam_mgr = Beam_energy_Manager::get_manager();
  Beam_mgr.remove();

  TPntFDDB fddb("rif");
  TPntDB db(fddb, "benergy");
  if(!db.IsOK()){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Error : Cannot access rif::beam_profile.\n";
    return 1;
  }

  if(db.Get( expno , runno , ver ) <= 0){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Error = : Cannot get data from database correctly.\n";
    return 1;
  }
  input_matrix( 1, Beam_mgr, expno, runno, 0 );
  set_matrix_value(0);

  db.non_clear(BEAM_ENERGY);

  return 0;
}
//--------------------------------------------------------------------
void BeamEnergy::input_matrix( int flg, Beam_energy_Manager& Beam_mgr,
                               int expno, int runno, int i ){
  
 if( flg==0 ){
   mt_expno[i]        = expno;
   mt_runno[i]        = runno+i;
   mt_run_flag[i]     = -1; 
   mt_e_beam_corr[i]  = -1; 
   mt_e_beam_orig[i]  = -1; 
   mt_e_beam2[i]      = -1; 
   mt_e_beam_err[i]   = -1; 
   mt_e_her[i]        = -1; 
   mt_e_ler[i]        = -1; 
   mt_corss_angle[i]  = -1; 
   mt_e_her_orig[i]   = -1; 
   mt_e_ler_orig[i]   = -1;
 }

 if( flg==1 ){
  Beam_energy_Manager::iterator beit = Beam_mgr.begin();
 
  mt_expno[i]        = expno;
  mt_runno[i]        = runno+i;
  mt_run_flag[i]     = beit->run_flag(); 
  mt_e_beam_corr[i]  = beit->E_beam();
  mt_e_beam_orig[i]  = beit->E_beam_orig();
  mt_e_beam2[i]      = beit->E_beam2();
  mt_e_beam_err[i]   = beit->E_beam_err(); 
  mt_e_her[i]        = beit->E_HER(); 
  mt_e_ler[i]        = beit->E_LER();
  mt_corss_angle[i]  = beit->CROSS();
  mt_e_her_orig[i]   = beit->E_HER_orig();
  mt_e_ler_orig[i]   = beit->E_LER_orig();
 }

  Beam_mgr.remove();

  return;
}
//--------------------------------------------------------------------
void BeamEnergy::check_runno(){

  if(check_beginrun() == -1){
    dout(Debugout::ERR,"Beamenergy") << "[BeamEnergy] Warning : Calling BeamEnergy::begin_run() again."
	      << std::endl;
    begin_run();
  }
}
//--------------------------------------------------------------------
HepLorentzVector BeamEnergy::p_beam(void){

  // Energy information
  const double mass_e = 0.0;    //0.000510998902;
  const double her = m_e_her;
  const double ler = m_e_ler;

  HepLorentzVector p4_her(0.0, 0.0, +sqrt(her*her-mass_e*mass_e), her);
  p4_her.rotateY(m_corss_angle);
  HepLorentzVector p4_ler(0.0, 0.0, -sqrt(ler*ler-mass_e*mass_e), ler);

  return p4_her + p4_ler;
}
//---------------------------------------------------------------------
HepLorentzVector BeamEnergy::p_beam2(void){

  // Energy information
  const double mass_e = 0.0;    //0.000510998902;
  const double her = m_e_her*(m_e_beam2/m_e_beam_corr);
  const double ler = m_e_ler*(m_e_beam2/m_e_beam_corr);

  HepLorentzVector p4_her(0.0, 0.0, +sqrt(her*her-mass_e*mass_e), her);
  p4_her.rotateY(m_corss_angle);
  HepLorentzVector p4_ler(0.0, 0.0, -sqrt(ler*ler-mass_e*mass_e), ler);

  HepLorentzVector p_beam2 = p4_her + p4_ler;
  return p_beam2;
}
//--------------------------------------------------------------------
double BeamEnergy::Ecm(void){

  // Energy information
  const double mass_e = 0.0;      //0.000510998902;
  const double her = m_e_her;
  const double ler = m_e_ler;

  HepLorentzVector p4_her(0.0, 0.0, +sqrt(her*her-mass_e*mass_e), her);
  p4_her.rotateY(m_corss_angle);
  HepLorentzVector p4_ler(0.0, 0.0, -sqrt(ler*ler-mass_e*mass_e), ler);

  HepLorentzVector cm = p4_her+p4_ler;
  
  return cm.mag();
}
//--------------------------------------------------------------------
HepLorentzVector BeamEnergy::p_cm( HepLorentzVector p ){

  HepLorentzVector pbeam = p_beam();
  p.boost(-pbeam.boostVector());

  return p;
}
//--------------------------------------------------------------------
Hep3Vector BeamEnergy::CMBoost(void){

  // Energy information
  const double mass_e = 0.0;       //0.000510998902;
  const double her = m_e_her;
  const double ler = m_e_ler;

  HepLorentzVector p4_her(0.0, 0.0, +sqrt(her*her-mass_e*mass_e), her);
  p4_her.rotateY(m_corss_angle);
  HepLorentzVector p4_ler(0.0, 0.0, -sqrt(ler*ler-mass_e*mass_e), ler);

  HepLorentzVector cm = p4_her+p4_ler;

  return cm.boostVector();
}
//--------------------------------------------------------------------
void BeamEnergy::dump(void){

 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Beam Energy Summary of this events." << std::endl;
 if(m_usable){
    dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] available." << std::endl;
  }else{
    dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] !NOT! available." << std::endl;
  }

 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Exp No." << m_exp << ", Run No." << m_run << std::endl;
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] DB Version No. " << m_version 
      << " ( 1 = Run-dependent version ; 2 = Run-independent version )" << std::endl; 
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Run Dependent Flag is " << m_flag_run 
      << "( 0 = on; 1 >= off ) " << std::endl; 
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Beam Energy = " << m_e_beam_corr << std::endl;
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Origianl Beam Energy = " << m_e_beam_orig << std::endl;
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Beam Energy for on_resonance = " << m_e_beam2 << std::endl;
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Error od Beam Energy = " << m_e_beam_err << std::endl; 
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Energy of HER = " << m_e_her << std::endl; 
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Energy of LER = " << m_e_ler << std::endl; 
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Cross angle   = " << m_corss_angle << std::endl; 
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Origianl Energy of HER = " << m_e_her_orig << std::endl;  
 dout(Debugout::DUMP,"Beamenergy") << "[BeamEnergy] Origianl Energy of LER = " << m_e_ler_orig << std::endl;  

 return;
}
//--------------------------------------------------------------------
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
