//
// $Id: BeamEnergy.h 10013 2007-02-27 06:55:37Z katayama $
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
//==============================================================
//  File Name :  BeamEnergy.h
//--------------------------------------------------------------
//  First Version 2004.01.07 by T.Shibata
//  Mofidied : 2004.02.29 
//==============================================================
#ifndef __BEAMENERGY_H__
#define __BEAMENERGY_H__
//--------------------------------------
#include "belle.h"
#include "event/BelleEvent.h"
#include <belleCLHEP/Vector/LorentzVector.h>
#include <belleCLHEP/Vector/ThreeVector.h>
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h"
#include "panther/panther_manager.h"

#include BELLETDF_H
#include RUN_INFO_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
//==============================================================
class BeamEnergy{
public:
  // constructor
  BeamEnergy(){};

  // destructor
  virtual ~BeamEnergy(){};

  // register variables as module parameter
  static void define_global(Module_descr* dscr);

  // Set BeamEnergy
  static void begin_run( const int version  = m_version,
                         const int flag_run = m_flag_run,
                         const int flag_db  = m_flag_db );

  static const int& Run_flag() { return m_run_flag; }; 
  static const double& E_beam_corr() { return m_e_beam_corr; };
  static const double& E_beam_orig() { return m_e_beam_orig; };
  static const double& E_beam2() { return m_e_beam2; };
  static const double& E_beam_err() { return m_e_beam_err; };
  static const double& E_HER() { return m_e_her; };
  static const double& E_LER() { return m_e_ler; };
  static const double& Cross_angle() { return m_corss_angle; };
  static const double& E_HER_orig() { return m_e_her_orig; };
  static const double& E_LER_orig() { return m_e_ler_orig; };

  static HepLorentzVector p_beam(void);
  static HepLorentzVector p_beam2(void);
  static double Ecm(void);
  static Hep3Vector CMBoost(void);
  static HepLorentzVector p_cm( HepLorentzVector p );

  // Return Experimental number
  static int ExpNo() { return m_exp; };
  static int RunNo() { return m_run; };
  static int Version() { return m_version; };

  // Optical dump of data
  static void dump(void);

  // BeamEnergy is valid for run #m_run
  static bool valid() { return m_valid; };

  // BeamEnergy is usable for analysis
  static bool usable() { 
    check_runno();
    return m_usable;
  };

  // run by run beam energy called outside beam energy class A.I
  static int GetBeamEnergy4(const int version, 
			    const int expno, 
			    const int runno);



public:
  // Modifiers

  // set Exp# and run#
  static void ExpRunNo(const int exp, const int  run) {
    m_exp = exp;
    m_run = run;
    return;
  };

  // set Module parameters
  static int flag_rundep( const int flag_run ) { return m_flag_run = flag_run; };
  static int version(const int ver) { return m_version = ver;   };
  static int flag_database(const int flag_db) { return m_flag_db = flag_db;   };

  static bool valid(const bool vld) { return m_valid = vld;     };
  static bool usable(const bool usbl) { return m_usable = usbl; };

private:
  static void set_default_value(void);
  static void set_2nd_default_value(void);
  static void set_matrix_value(int i);

  static int GetBeamEnergy(const int version);
  static int GetBeamEnergy2(const int version);
  static int GetBeamEnergy3(const int version);
  static void input_matrix( int flg, Beam_energy_Manager& Beam_mgr,
                            int expno, int runno, int i );

  static void check_runno();

  static bool    m_MC;
  static int     m_exp;
  static int     m_run;
  static bool    m_valid;
  static bool    m_usable;

  static int m_version;
  static int m_flag_run;
  static int m_flag_db;

  static int     defaut_run_flag;
  static double  defaut_e_beam_corr;
  static double  default_e_beam_orig;   
  static double  default_e_beam2;
  static double  default_e_beam_err; 
  static double  default_e_her;
  static double  default_e_ler;     
  static double  default_corss_angle; 
  static double  default_e_her_orig;
  static double  default_e_ler_orig;

  static double  default_2_e_her;
  static double  default_2_e_ler;
 
  static int     m_run_flag;
  static double  m_e_beam_corr;
  static double  m_e_beam_orig;
  static double  m_e_beam2;
  static double  m_e_beam_err;
  static double  m_e_her;
  static double  m_e_ler;
  static double  m_corss_angle;
  static double  m_e_her_orig; 
  static double  m_e_ler_orig;

  static int     mt_expno[100];
  static int     mt_runno[100];
  static int     mt_run_flag[100];
  static double  mt_e_beam_corr[100];
  static double  mt_e_beam_orig[100];
  static double  mt_e_beam2[100];
  static double  mt_e_beam_err[100];
  static double  mt_e_her[100];
  static double  mt_e_ler[100];
  static double  mt_corss_angle[100];
  static double  mt_e_her_orig[100];
  static double  mt_e_ler_orig[100];

};
//==============================================================
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __BEAMENERGY_H__ */
