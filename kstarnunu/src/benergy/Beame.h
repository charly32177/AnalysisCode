//=================================================================
//  File Name : beame.h
//----------------------------------------------------------------
//  Creation ; 2004.01.10
//  Author   ; T.Shibata
//----------------------------------------------------------------
#ifndef BEAME_H
#define BEAME_H
//----------------------------------------------------------------

#include "belle.h"
#include "event/BelleEvent.h"
#include "panther/panther.h"

#include "basf/module.h"
#include "basf/module_descr.h"

#include "belleCLHEP/Vector/LorentzVector.h"
#include "belleCLHEP/Vector/ThreeVector.h"

#include "benergy/BeamEnergy.h"

#include BELLETDF_H
#include RUN_INFO_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
//============================================================
// Define the Module Class
//============================================================
class beame : public Module{
 
 public:
  beame();
  ~beame(){};
  void init(int*){}; 
  void hist_def(){};
  void begin_run( BelleEvent*,int* );
  void disp_stat(const char*){};
  void event(BelleEvent*, int*);
  void end_run(BelleEvent*, int*){};
  
  void get_event_id( int &no_exp, int &no_run, int &no_evt );

};
//============================================================
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif




