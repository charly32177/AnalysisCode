// -*- C++ -*-
//
// Package:     <package>
// Module:      k_sfw
// 
// Description: <one line class summary>
//
// Usage:
//
// Author:      Hidekazu Kakuno
// Created:     ²Ð  7·î  8 18:56:42 JST 2003
// $Id: k_sfw.h,v 1.2 2003/07/17 18:11:54 kakuno Exp $
//
// Revision history
//
// $Log: k_sfw.h,v $
// Revision 1.2  2003/07/17 18:11:54  kakuno
// k_sfw version 2.3: Modified for PDF(mm2) correction
//
// Revision 1.1.1.1  2003/07/17 17:54:38  kakuno
// k_sfw version 2.2:  Initial version
//
#include "belle.h"

#if !defined(PACKAGE_K_SFW_H_INCLUDED)
#define PACKAGE_K_SFW_H_INCLUDED

#include <cmath>

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


class brutus_f;
class Particle;
class BelleTuple;
class k_sfw
{
public:

  // initializer of brutus_f
  // correct_mm2 switch:
  //   0: no correction
  //   1: correction on signal PDF(mm2)
  //   2: correction on background PDF(mm2)
  //   3: correction on signal & background PDF(mm2)
  static void initialize(brutus_f* brutus_fisher,
			 const unsigned correct_mm2 = 0);

  // Constructor
  k_sfw(const Particle& b);

  // Destractor
  ~k_sfw(){}

  // returns fisher discriminant
  double fd() const;

  // returns likelihood ratio
  double lr() const;

  // returns charged so component
  double Hso_c(const int order) const;

  // returns neutral so component
  double Hso_n(const int order) const;

  // returns missing so component
  double Hso_m(const int order) const;

  // returns oo component
  double Hoo(const int order) const;

  // returns sum of transverse energy
  double e_t() const;

  // returns ID of e_miss-p_miss intervals
  double mm2() const;

  // returns ID of e_miss-p_miss intervals
  unsigned int i_mm2() const;

private:
  double legendre(const double z, const int i) const;

  inline double squ(const double a) const {return a*a;}

  // data members
  double m_Hso[3][5];
  double m_Hoo[5];
  double m_et;
  double m_fd;
  double m_mm2;
  static brutus_f* bf;
  static double alpha_k_sfw0[17];
  static double alpha_k_sfw1[17];
  static double alpha_k_sfw2[17];
  static double alpha_k_sfw3[17];
  static double alpha_k_sfw4[17];
  static double alpha_k_sfw5[17];
  static double alpha_k_sfw6[17];
  static BelleTuple *nt;
  static double mm2_correction_sig[7];
  static double mm2_correction_bkg[7];
  static double sp[7][9];
  static double bp[7][9];
};

//---------------------------------------
#ifdef  K_SFW_NO_INLINE
#define inline
#else
#undef  inline
#define K_SFW_INLINE_DEFINE_HERE
#endif /* K_SFW_NO_INLINE */

#ifdef    K_SFW_INLINE_DEFINE_HERE

inline
double
k_sfw::fd() const
{
  return m_fd;
}

inline
double
k_sfw::lr() const
{
  const int i = i_mm2();
  const double a(fd());
  double s = (a < sp[i][1])
    ? (sp[i][0]*sp[i][2]/(sp[i][2]+sp[i][3])/(std::sqrt(2.*M_PI)*sp[i][2])
       *std::exp(-0.5*squ((a-sp[i][1])/sp[i][2]))
       + std::fabs(sp[i][4])*std::exp(-0.5*squ((a-sp[i][5])/sp[i][6])))/sp[i][7]*sp[i][8]
    : (sp[i][0]*sp[i][3]/(sp[i][2]+sp[i][3])/(std::sqrt(2.*M_PI)*sp[i][3])
       *std::exp(-0.5*squ((a-sp[i][1])/sp[i][3]))
       + std::fabs(sp[i][4])*std::exp(-0.5*squ((a-sp[i][5])/sp[i][6])))/sp[i][7]*sp[i][8];
  double b = (a < bp[i][1])
    ? (bp[i][0]*bp[i][2]/(bp[i][2]+bp[i][3])/(std::sqrt(2.*M_PI)*bp[i][2])
       *std::exp(-0.5*squ((a-bp[i][1])/bp[i][2]))
       + std::fabs(bp[i][4])*std::exp(-0.5*squ((a-bp[i][5])/bp[i][6])))/bp[i][7]*bp[i][8]
    : (bp[i][0]*bp[i][3]/(bp[i][2]+bp[i][3])/(std::sqrt(2.*M_PI)*bp[i][3])
       *std::exp(-0.5*squ((a-bp[i][1])/bp[i][3]))
       + std::fabs(bp[i][4])*std::exp(-0.5*squ((a-bp[i][5])/bp[i][6])))/bp[i][7]*bp[i][8];
  if (s < 0) s = 0.;
  if (b < 0) b = 0.;
  return (s+b > 0) ? s/(s+b) : 0.5;
}

inline
unsigned int
k_sfw::i_mm2() const
{
  if      (m_mm2 < -0.5) return 0;
  else if (m_mm2 <  0.3) return 1;
  else if (m_mm2 <  1.0) return 2;
  else if (m_mm2 <  2.0) return 3;
  else if (m_mm2 <  3.5) return 4;
  else if (m_mm2 <  6.0) return 5;
  return 6;
}

inline
double
k_sfw::Hso_c(const int order) const
{
  return m_Hso[0][order];
}

inline
double
k_sfw::Hso_n(const int order) const
{
  return m_Hso[1][order];
}

inline
double
k_sfw::Hso_m(const int order) const
{
  return m_Hso[2][order];
}

inline
double
k_sfw::Hoo(const int order) const
{
  return m_Hoo[order];
}

inline
double
k_sfw::e_t() const
{
  return m_et;
}

inline
double
k_sfw::mm2() const
{
  return m_mm2;
}

inline
double
k_sfw::legendre(const double z, const int i) const{
  switch(i){
  case 0:
    return 1.;
  case 1:
    return z;
  case 2:
    return 1.5*z*z - 0.5;
  case 3:
    return z*(2.5*z*z - 1.5);
  case 4:
    return (4.735*z*z*z*z - 3.75*z*z + 0.375);
  default:
    return 0;
  }
}

#endif /* K_SFW_INLINE_DEFINE_HERE */

#undef inline

#endif /* PACKAGE_K_SFW_H_INCLUDED */



#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

