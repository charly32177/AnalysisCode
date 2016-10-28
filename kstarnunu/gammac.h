#ifndef GAMMA_UTILITY_CLASS_H
#define GAMMA_UTILITY_CLASS_H
/*
   Files : gammac.h and gammac.cc

   UTILITY CLASS FOR GAMMA

   By jtanaka@belaxp1.phys.s.u-tokyo.ac.jp
   ---------------------------------------
   
   Note1 : This class recalculates gamma momentum and its error matrix given a vertex. 

   Note2 : I think "GammaParticle" is a strange name. 
           But it is not coflicted with classes provided by CsI groups in future.
*/
#include "belle.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Geometry/Point3D.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class GammaParticle
{
public:
  GammaParticle(const double energy,
                const double phi,
                const double theta,
                const double R,
                const HepSymMatrix & Eg);

  virtual ~GammaParticle();
  
  //Selectors
  const Hep3Vector       &  momentum(void){ return m_momentum; }
  const HepLorentzVector &  momentumEnergy(void){ return m_momentumEnergy; }
  const HepSymMatrix     &  errorMomentum(void){ return m_errorMomentum; }
  const HepSymMatrix     &  errorMomentumEnergy(void){ return m_errorMomentumEnergy; }

  const HepPoint3D   & vertex(void) { return m_vertex; }
  const HepSymMatrix & errorVertex(void){ return m_errorVertex; }

  const double energy(void){ return m_energy; }
  const double theta(void) { return m_theta; }
  const double phi(void)   { return m_phi; }
  const double R(void)     { return m_R; }
  const HepSymMatrix & errorGamma(void){ return m_errorGamma; }

  //Modifiers
  const HepPoint3D & vertex(const HepPoint3D &vertex);
  const HepPoint3D & vertex(const HepPoint3D &vertex, const HepSymMatrix &error);

  //Operators
  GammaParticle & operator = (const GammaParticle &);
  
private:
  //Mathmatical functions
  HepMatrix delPEdelG(void);
  HepMatrix delPEdelGV(void);
  
  //calculation new gamma
  void calcu(void);
  unsigned calcuVector(void);
  unsigned calcuError(void);
  
  //data
  double       m_energy, m_theta, m_phi;
  HepSymMatrix m_errorGamma;
  HepPoint3D   m_vertex;
  HepSymMatrix m_errorVertex;
  HepSymMatrix m_errorGammaVertex;
  Hep3Vector       m_momentum;
  HepSymMatrix     m_errorMomentum;
  HepLorentzVector m_momentumEnergy;
  HepSymMatrix     m_errorMomentumEnergy;

  double m_R; // cm

  unsigned m_ipFlag;
};
#endif /* GAMMA_UTILITY_CLASS_H */

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

