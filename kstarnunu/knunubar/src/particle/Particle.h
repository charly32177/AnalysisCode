//
// Particle.h
// 
// A class definition for <Particle> class 
// for the BELLE standard (hopefully!) Particle object liblary
// 
// <Particle> class supplies you interfaces to various particle information,
// such as momentum, pid, etc. through its private member objects,
// Momentum, PID, Ptype, and Relation.
//
// Filename : Particle.h
// Author : YOKOyama Masashi
// e-mail:yokoyama@belaxp1.phys.s.u-tokyo.ac.jp
// --> e-mail:jtanaka@belaxp1.phys.s.u-tokyo.ac.jp
// 
// $Id: Particle.h 10002 2007-02-26 06:56:17Z katayama $ 
//
// Revision history
// 
// $Log$
// Revision 1.31  2001/12/14 02:16:24  katayama
// MDST_OBS removed
//
// Revision 1.30  2001/12/13 15:31:54  katayama
// MDST_OBS
//
// Revision 1.29  2001/12/12 07:10:53  jtanaka
// (1) compatibility for obsoleted Mdst_vee/mdst_sim_xref table
// (2) compatibility for gcc3
//
// Revision 1.28  2000/05/16 14:46:27  jtanaka
// added constructor of Mdst_ecl.
//
// Revision 1.27  2000/04/27 12:57:58  jtanaka
// add new constructor Particle(Mdst_charged, string <-- "new", pivot).
//
// Revision 1.26  2000/04/14 12:40:57  jtanaka
// updated for new table "mdst_vee2".
//
// Revision 1.25  2000/03/08 07:56:20  jtanaka
// modify a constructor of the Mdst_charged and
// add new functions(please use the lastest kfitter:May 8,2000)
//
// Revision 1.24  2000/03/07 11:14:00  katayama
// compatibility with CC5.0
//
// Revision 1.23  2000/01/05 06:10:38  jtanaka
// major updates:please see BELLE whiteboard.
//
// Revision 1.22  1999/12/17 04:10:50  yiwasaki
// Particle::dump implemented
//
// Revision 1.21  1999/04/27 12:15:42  jtanaka
// added child() and mc().
//
// Revision 1.20  1999/04/09 15:01:08  jtanaka
// Added new member functions, "const", "&"...etc. Removed some member functions.
//
// Revision 1.19  1999/01/16 10:31:30  katayama
// clean up includes
//
// Revision 1.18  1999/01/09 12:50:59  jtanaka
// Bugs fix, added some members in Ptype Class, removed some member from Particle.h because of functions for setting member of other classes.
//
// Revision 1.17  1999/01/08 06:08:08  katayama
// Do not create PID objects as default
//
// Revision 1.16  1998/10/20 14:45:50  jtanaka
// added FittedMomentum.
//
// Revision 1.15  1998/10/15 09:52:14  jtanaka
// bug fix, and added deep_copy and deep_delete in Particle.h
//
// Revision 1.14  1998/10/12 19:25:45  jtanaka
// updated constructors Momentum.h, Particle.h, Relation.h. removed and added some functions in PID.h. added some static_PIDs to static_particle.h.
//
// Revision 1.13  1998/10/05 10:59:29  jtanaka
// added some comments for doc++ in some headers files and modified static objects in Relation.h, ParticleManager.h.
//
// Revision 1.12  1998/09/09 08:29:32  jtanaka
// added new constructor to Particle.
//
// Revision 1.11  1998/09/03 23:55:36  jtanaka
// We need to include "Particle.h" only. Comment out some members of Ptype.
//
// Revision 1.10  1998/07/23 12:51:41  katayama
// conform ANSI
//
// Revision 1.9  1998/07/22 13:04:38  jtanaka
// add some functions and modify const functions a little.
//
// Revision 1.8  1998/07/13 10:55:17  jtanaka
// modify a little about const functions and add vertex to Momentum
//
// Revision 1.7  1998/07/03 14:57:00  jtanaka
// add some constructors for Ptype.
//
// Revision 1.6  1998/07/03 12:10:13  jtanaka
// modified some parts
//
// Revision 1.5  1998/07/02 09:29:09  higuchit
// null flag -> `usable' flag
//
// Revision 1.4  1998/07/01 11:54:52  jtanaka
// add m_null.
//
// Revision 1.3  1998/06/19 09:09:25  jtanaka
// *** empty log message ***
//
// Revision 1.2  1998/06/15 07:30:16  yokoyamm
// First version..
//

#ifndef PARTICLE_CLASS_PARTICLE_H
#define PARTICLE_CLASS_PARTICLE_H

#include "belle.h"
#include <string>

#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Vector/LorentzVector.h"
#include "particle/Ptype.h"
#include "particle/PID.h"
#include "particle/Momentum.h"
//#include "particle/FittedMomentum.h"
#include "particle/constant.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Mdst_charged;
class Mdst_gamma;
class Mdst_trk;
class Mdst_vee2;
class Mdst_klong;
class Mdst_pi0;
class Mdst_ecl;
class Gen_hepevt;



//class PID;
class ParticleUserInfo;

//...wants to write [#include "particle/Relation.h"] in this line. But it is a compiler error. So, [#include "particle/Relation.h"] is written on the end of this file. 
class Relation;

/// Particle class supplies you interfaces to various particle information, such as momentum, pid, etc. through its private member objects, Momentum, PID, Ptype, and Relation.
class Particle {
  
public:
  /// Default constructor
  Particle();
  
  /// Copy constructor
  Particle(const Particle &);
  
  /// Constructor with momentum
  Particle(const HepLorentzVector &, const Ptype &);
  /// Constructors with Momentum
  Particle(const Momentum &, const Ptype &);
  
  /// Constructor with Mdst\_charged(by pType)
  Particle(const Mdst_charged &, const Ptype &,
	   const HepPoint3D & = HepPoint3D(0.,0.,0.));
  /// Constructor with Mdst\_charged(by name)
  Particle(const Mdst_charged &, const std::string &,
	   const HepPoint3D & = HepPoint3D(0.,0.,0.));
  /// Constructor with Mdst\_gamma
  Particle(const Mdst_gamma &);
  /// Constructor with Mdst\_vee2
  Particle(const Mdst_vee2 &, const bool makeRelation = true);
  /// Constructor with Mdst\_klong 
  Particle(const Mdst_klong &);
  /// Constructor with Mdst\_pi0 
  Particle(const Mdst_pi0 &, const bool makeRelation = true);
  /// Constructor with Mdst\_ecl
  Particle(const Mdst_ecl &);
  /// Constructor with Gen\_hepevt
  Particle(const Gen_hepevt &, const Ptype &);
  /// Constructor with Gen\_hepevt
  Particle(const Gen_hepevt &);

  /// Destructor
  virtual ~Particle();

public:// General interfaces
  /// returns object name.
  virtual const std::string & name(void) const;

  /// sets object name and returns it.
  virtual const std::string & name(const std::string &);

  /// returns class name.
  virtual std::string className(void){ return std::string("Particle"); }
  
  /// dumps debug information. Keywords are 'mass', 'momentum', 'position', 'recursive', and 'full'. 'recursive' is to dump children also. 'full' is equivalant to specifying all keywords.
  virtual void dump(const std::string & keyword = std::string("mass momentum"),
		    const std::string & prefix = std::string("")) const;

  /// create PID object
  virtual void createPID(void);

  /// append daughter
  virtual Particle &append_daughter(Particle &);

public:// Momentum interfaces
  /// returns a const reference to Momentum Class.
  virtual const Momentum & momentum(void) const;
  
  /// returns a reference to Momentum Class.
  virtual Momentum & momentum(void);
  
  /// sets a reference to Momentum Class and returns it.
  virtual const Momentum & momentum(const Momentum &);
  
  /// returns momentum.
  virtual const HepLorentzVector & p(void) const;

  /// returns 3momentum. (not reference)
  virtual const Hep3Vector p3(void) const;

  /// returns position.
  virtual const Hep3Vector & x(void) const;

  /// returns magnitude of momentum.
  virtual double ptot(void) const;

  /// returns x component of momentum.
  virtual double px(void) const;

  /// returns y component of momentum.
  virtual double py(void) const;

  /// returns z component of momentum.
  virtual double pz(void) const;

  /// returns energy.
  virtual double e(void) const;

  /// returns invariant mass.
  virtual double mass(void) const;

public:// PID interfaces
  /// returns a const reference to PID.
  virtual const PID & pId(void) const;

  /// returns a reference to PID.
  virtual PID & pId(void);

  /// sets a reference to PID and returns it.
  virtual const PID & pId(const PID &);
  
public:// Ptype interfaces
  /// returns a const reference to Ptype.
  virtual const Ptype & pType(void) const;

  /// returns a reference to Ptype.
  virtual Ptype & pType(void);
  
  /// sets a reference to Ptype and returns it.
  virtual const Ptype & pType(const Ptype &);
  
  /// returns charge. (in units of e)
  virtual double charge(void) const;
  
  /// returns LUND7 particle code.
  virtual int lund(void) const;
  
public:// Relation interfaces
  /// returns a const reference to Relation.
  virtual const Relation & relation(void) const;

  /// returns a reference to Relation.
  virtual Relation & relation(void);
  
  /// sets a reference to Relation and returns it.
  virtual const Relation & relation(const Relation &);
  
  /// returns a reference to mother.
  virtual const Particle & mother(void) const;
  
  /// returns a number of children.
  virtual unsigned nChildren(void) const;

  /// returns a reference to i'th child.
  virtual const Particle & child(unsigned i) const;

  /// returns a reference to i'th child.
  virtual Particle & child(unsigned i);
 
  /// returns a reference to MC particle.
  virtual const Particle & mc(void) const;

  /// returns a reference to MC particle.
  virtual Particle & mc(void);

  /// returns a reference to Mdst\_charged.
  virtual const Mdst_charged & mdstCharged(void) const;

  /// returns a reference to Mdst\_gamma.
  virtual const Mdst_gamma & mdstGamma(void) const;

  /// returns a reference to Mdst\_trk.
  virtual const Mdst_trk & mdstTrk(void) const;
  
  /// returns a reference to Mdst\_vee2.
  virtual const Mdst_vee2 & mdstVee2(void) const;

  /// returns a reference to Mdst\_klong.
  virtual const Mdst_klong & mdstKlong(void) const;

  /// returns a reference to Mdst\_pi0.
  virtual const Mdst_pi0 & mdstPi0(void) const;

  /// returns a reference to Mdst\_ecl.
  virtual const Mdst_ecl & mdstEcl(void) const;

  /// returns a reference to Gen\_hepevt.
  virtual const Gen_hepevt & genHepevt(void) const;

public:// FittedMomentum interfaces
  // returns a const reference to FittedMomentum.
  //virtual const FittedMomentum & fittedMomentum(void) const;

  // returns a reference to FittedMomentum.
  //virtual FittedMomentum & fittedMomentum(void);

  // sets a reference to FittedMomentum and returns it.
  //virtual const FittedMomentum & fittedMomentum(const FittedMomentum &);

public:// User Definition Object
  /// returns a pointer of "user definition object".
  virtual const ParticleUserInfo & userInfo(void) const;

  /// returns a pointer of "user definition object".
  virtual ParticleUserInfo & userInfo(void);

  /// sets and returns a pointer of "user definition object".
  virtual const ParticleUserInfo & userInfo(const ParticleUserInfo &);
  
public:// Operators
  /// copy operator
  Particle & operator = (const Particle &);

  /// bool operator : returns "usable" info.
  operator bool() const { return m_usable; }
  
  /// ! operator : returns "!usable" info.
  bool operator ! () const { return !m_usable; }

  /// returns "usable" info. (please use usable() not this func.)
  virtual bool isUsable(void){ return usable(); }

  /// returns "usable" info.
  virtual bool usable(void){ return m_usable; }

  /// sets "usable" info and returns it.
  virtual bool usable(const bool &);

  // sets "not usable" info and returns it.
  //virtual bool unusable(bool = UNUSABLE);  

  /// copies(all private member is created by "new") and returns it.
  virtual Particle deepCopy(void);

  /// deletes children's objects made by deepCopy().
  virtual void deepDelete(void);

private:
  bool m_usable; 

  std::string    m_name;
  Momentum *m_momentum;
  PID      *m_pId;
  Relation *m_relation;
  Ptype    *m_pType;
  //FittedMomentum *m_fittedMomentum;

  ParticleUserInfo *m_userInfo;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#include "particle/Relation.h"
#endif /* PARTICLE_CLASS_PARTICLE_H */
