/*
  	File : findLambda.h
  	Description : select good lambda from MDST_Vee2
	Version : 0.2

  	Author : Kai-Feng Chen, NTU
  	Date : Mar/30,2001
*/
#include "belle.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/Vector.h"

#include "helix/Helix.h"

#include "panther/panther.h"
#include MDST_H

#define MIN(m,n) ((m)<(n)?(m):(n))

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Class definition
class FindLambda
{
public:
    FindLambda() {};
    ~FindLambda() {};

    void candidates ( const Mdst_vee2& Cand, HepPoint3D runIP ) {
    	m_IP = runIP;
	m_kind = Cand.kind();
	m_chisq = Cand.chisq();
	
	HepPoint3D Vtx(Cand.vx(),Cand.vy(),Cand.vz());
	HepPoint3D VtxmIP = Vtx - m_IP;
	    
	m_fl = sqrt(VtxmIP.x()*VtxmIP.x()+VtxmIP.y()*VtxmIP.y());
	
	double dot =  VtxmIP.x()*Cand.px() + VtxmIP.y()*Cand.py();
    	double atot = m_fl*sqrt(Cand.px()*Cand.px() + 
	    	    	      Cand.py()*Cand.py());
			      
    	m_dphi = ((atot==0.0) ? 0.0:acos(dot/atot));
	
	double dr1,dr2,dz1,dz2;
	
	if (m_kind==2) {
	    GetImpactParameters(&Cand.chgd(0),&dr1,&dz1,4);
	    GetImpactParameters(&Cand.chgd(1),&dr2,&dz2,2);
	}else {
	    GetImpactParameters(&Cand.chgd(1),&dr1,&dz1,4);
	    GetImpactParameters(&Cand.chgd(0),&dr2,&dz2,2);
	}
	
	m_dr = MIN(fabs(dr1),fabs(dr2));
	m_dz = MIN(fabs(dz1),fabs(dz2));

    	m_pmag = sqrt(	Cand.px()*Cand.px()+
	    	    	Cand.py()*Cand.py()+
			Cand.pz()*Cand.pz());  

    	m_zdist=Cand.z_dist();
    }
  
    HepPoint3D IP()    { return m_IP;  	}
    int kind()      { return m_kind;	}
    double dr()     { return m_dr;  	}
    double dz()     { return m_dz;  	} 
    double dphi()   { return m_dphi; 	}
    double zdist()  { return m_zdist; 	}
    double fl()     { return m_fl;  	}
    double pmag()   { return m_pmag; 	}
    double chisq()  { return m_chisq; 	}
    
    void GetImpactParameters(Mdst_charged *charged, double *dr, double *dz, int t) { 	
	
	*dr = 1000.0;
	*dz = 1000.0;
	
	if(charged->trk()){	
	    Mdst_trk_fit& trkFit = charged->trk().mhyp(t);
    	    if(trkFit){      
      	    	HepVector a(5,0);    
      	    	a[0] = trkFit.helix(0);
      	    	a[1] = trkFit.helix(1);
      	    	a[2] = trkFit.helix(2);
      	    	a[3] = trkFit.helix(3);
      	    	a[4] = trkFit.helix(4);

      	    	HepPoint3D pivot(trkFit.pivot(0), trkFit.pivot(1), trkFit.pivot(2));
      	    	Helix ltrk(pivot, a);         
      	    	ltrk.pivot(m_IP);          
      	    	
		*dr = ltrk.dr();
		*dz = ltrk.dz();		      
    	    }
    	}
    }
    
    int goodLambda() {
    	if (m_kind==2 || m_kind==3) {    
    	if (m_pmag>=1.5) { 
	    if (m_zdist< 7.7 && m_dr>0.018 && m_dphi<0.07 && m_fl>.35) return 2;
	    if (m_zdist<12.9 && m_dr>0.008 && m_dphi<0.09 && m_fl>.22) return 1;	    
	}
	if (m_pmag>=0.5 && m_pmag<1.5) {
	    if  (m_zdist<2.1 && m_dr>0.033 && m_dphi<0.1 && m_fl>.24) return 2;
	    if  (m_zdist<9.8 && m_dr>0.010 && m_dphi<0.18 && m_fl>.16) return 1;
	}
	if (m_pmag<0.5) { 
	    if  (m_zdist<1.9 && m_dr>0.059 && m_dphi<0.6 && m_fl>.17) return 2;
	    if  (m_zdist<2.4 && m_dr>0.027 && m_dphi<1.2 && m_fl>.11) return 1;	    
	}
	}
	return 0;
    }   

private: 
    HepPoint3D m_IP;
    int m_kind;
    double m_dr;
    double m_dz;
    double m_dphi;
    double m_fl;
    double m_zdist;
    double m_pmag;
    double m_chisq;
}; 


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

