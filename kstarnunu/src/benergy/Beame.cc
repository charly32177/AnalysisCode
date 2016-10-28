//================================================================
//  File Name : beame.cc
//----------------------------------------------------------------
//  Creation ; 2004.01.10
//  Author   ; T.Shibata
//  Modified 2004.03.02 by T.Shibata
//----------------------------------------------------------------
#include "belle.h"
#include <iostream>
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h"
#include "panther/panther_manager.h"
#include "benergy/Beame.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//===============================================================
beame::beame( ){}
//----------------------------------------------------------------
void beame::begin_run( BelleEvent*, int* ) {

  BeamEnergy::begin_run();

}  
//----------------------------------------------------------------
void beame::event( BelleEvent*, int* status){

  *status = 1;

  int no_exp, no_run, no_evt;
  get_event_id(no_exp, no_run, no_evt);

  dout(Debugout::INFO,"Beame") << "  Exp ,Run , Event = " << no_exp << " "
            << no_run << "  " <<  no_evt << std::endl;


  HepLorentzVector pbeam = BeamEnergy::p_beam();
  double E_beam_corr = BeamEnergy::E_beam_corr();
  double ecm = BeamEnergy::Ecm();
  Hep3Vector cmboost = BeamEnergy::CMBoost();

  dout(Debugout::INFO,"Beame") << " p_beam      = " << pbeam       << std::endl;
  dout(Debugout::INFO,"Beame") << " E_beam_corr = " << E_beam_corr << std::endl;
  dout(Debugout::INFO,"Beame") << " ecm         = " << ecm         << std::endl;
  dout(Debugout::INFO,"Beame") << " cmboost     = " << cmboost     << std::endl;

  BeamEnergy::dump(); 

  return;

}
//================================================================
void beame::get_event_id(int &no_exp, int &no_run, int &no_evt )
{
  no_exp=-1, no_run=-1, no_evt=-1;

  belle_event *belle_event;
  belle_event = (struct belle_event*)BsGetEnt(BELLE_EVENT,1,BBS_No_Index);
  if( belle_event ){
    no_exp = belle_event->m_ExpNo;
    no_run = belle_event->m_RunNo;
    no_evt = belle_event->m_EvtNo & 0x0fffffff;
  }
  return;
}
//================================================================
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
