#include "belle.h"
#include <math.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <vector>
#include "toolbox/FuncPtr.h"

#ifndef SHAPE_H
#define SHAPE_H

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class SuperFoxWolfram{
  private:
    double sum[5];
  public:
   SuperFoxWolfram();
   ~SuperFoxWolfram();
   void fill(std::vector<Vector4>&, std::vector<Vector4>&); 
   double R(int i);
};


void fwjet2(int NTrk, float *PTrk, float EBeam, 
           float *H0, float *R1, float *R2, float *R3, float *R4);

int thrust(int ntrk, double ptrk[][3], double* thr, double tvec[3], int itrk[]);

int spherp(int ntrk, double ptrk[][3], double *sper, double jet[3]);

float Sper(std::vector<Vector4>& ptl, Vector3& Bthr);

Vector3 thrust(std::vector<Vector4>& ptl);

/*class Nine_Cone{
  private:
    double m_cone[18];
  public:
    template <class Iterator,class Function>
    Nine_Cone(std::vector<Iterator>& it,Hep3Vector axis ,Function func){
 
      for(int i=0;i<18;i++) m_cone[i]=0.;
   
      for (std::vector<Iterator>::iterator go = it.begin(); go != it.end(); go++){
	        Iterator& p = *go;
		HepLorentzVector p4 = func(p);
		Hep3Vector p3 = (Hep3Vector)func(p);
		float costh = p3*axis;
		if(p4.e()<12. && p3.mag()*axis.mag()>0)
		  costh /= (p3.mag()*axis.mag());
		else
		  continue;
		double th= acos(costh) * 180./ 3.1415926;
		int level=(int)(th/10.);
		m_cone[level]+=p4.e();
      }

    };
    ~Nine_Cone(){};
    double cone(int i){
      if(i<1||i>9) return -1.;  
      else return (m_cone[i-1]+m_cone[18-i]);
    };
    double back_cone(int i){
      if(i<1||i>9) return -1.;  
      else return m_cone[18-i];
    };
    double front_cone(int i){    
      if(i<1||i>9) return -1.;      
      else return m_cone[i-1];
    };
};*/


#endif

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

