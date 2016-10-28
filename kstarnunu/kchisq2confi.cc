//
// $Id: kchisq2confi.cc 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.6  2002/03/27 23:41:57  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.5  2002/02/22 12:14:31  katayama
// math.h
//
// Revision 1.4  2000/03/07 17:09:00  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.3  1999/04/01 00:33:31  jtanaka
// add old interface.
//
// Revision 1.2  1999/03/29 05:39:13  jtanaka
// new class structure and new+old interfaces
//
// Revision 1.1  1997/10/04 05:30:05  katayama
// New from Tanaka san
//
//
//kchisq2confi of the KFitter
//
//jtanaka
//e-mail : jtanaka@hep.phys.s.u-tokyo.ac.jp
//


/* for erfc */
#include "belle.h"
#if defined(__sparc)
#  if defined(__EXTENSIONS__)
#    include <cmath>
#  else
#    define __EXTENSIONS__
#    include <cmath>
#    undef __EXTENSIONS__
#  endif
#elif defined(__GNUC__)
#  if defined(_XOPEN_SOURCE)
#    include <cmath>
#  else
#    define _XOPEN_SOURCE
#    include <cmath>
#    undef _XOPEN_SOURCE
#  endif
#endif

#if defined(__SUNPRO_CC)
// for erfc and other functions (lgamma and cbrt
#  include <math.h>
#endif

#include "kfitter/kfitter_ini.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//
//Orignal in CERNLIB and FSIM(by H.Ozaki)
//I change C++ from Fortran.
//I use probab.F of FSIM
//                        J.Tanaka
//version 1.0::1997/08/19
#if KF_WITH_OLD_INTERFACE
double chisq2Confi(const int n, const double chi2);

double 
Chisq2Confi(int n,double chi2)
{
  return chisq2Confi(n,chi2);
}
#endif //KF_WITH_OLD_INTERFACE

double 
chisq2Confi(const int n, const double chi2)
{
#define SRTOPI 0.7978846
#define UPL 170.0
#define ROOT2I 0.70710678

  double prob = 0.0;
  double sum,term;
  int m;
  int i,k;
  double temp_i,temp_n;
  double srty;

  if((n <= 0)||(chi2 < 0.0)){
    return prob;
  }
  if(n > 60){
    temp_n = (double)n;
    srty = std::sqrt(chi2) - std::sqrt(temp_n-0.5);
    if (srty < 12.0){
      prob = 0.5*erfc(srty);
      return prob;
    }
    return prob;
  }
  if(chi2 > UPL){
    return prob;
  }
  sum = exp( -0.5 * chi2 );
  term = sum;
  m = (int)floor(n/2);

  if( 2*m == n ){
    if( m == 1 ){
      prob = sum;
      return prob;
    }else{
      for(i=2;i<m+1;i++){
	temp_i = (double)i;
	term = 0.5*chi2*term/(temp_i-1.0);
	sum = sum + term;
      }
      prob = sum;
      return prob;
    }
  }else{
    srty = std::sqrt(chi2);
    prob = erfc(ROOT2I*srty);
    if(n == 1){
      return prob;
    }
    if(n == 3){
      prob = SRTOPI*srty*sum + prob;
      return prob;
    }else{
      k = m - 1;
      for(i=1;i<k+1;i++){
	temp_i = (double)i;
	term = term*chi2/(2.0*temp_i + 1.0);
	sum = sum + term;
      }
      prob = SRTOPI*srty*sum + prob;
      return prob;
    }
  }
}



#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
