#include "shape.h"
#include <math.h>
#include <iostream.h>
#define LGNDR2(X) ((3*X*X-1)/2)
#define LGNDR3(X) ((5*X*X*X-3*X)/2)
#define LGNDR4(X) ((35*X*X*X*X-30*X*X+3)/8)


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

Vector3 thrust(std::vector<Vector4>& ptl)
{
  int ntrk = ptl.size();
  double p_trk[ntrk][3],thr, tvect[3];
  int i_trk[50];

  for(int i=0;i<ntrk;i++)
  {
    p_trk[i][0]=ptl[i].px();
    p_trk[i][1]=ptl[i].py();
    p_trk[i][2]=ptl[i].pz();
  }
  
  thrust(ntrk, p_trk, &thr, tvect, i_trk);
  Vector3 thr_axis(tvect[0],tvect[1],tvect[2]);
  if(abs(thr_axis.mag())>1.5)
  {
    cout << "size "<<thr_axis.mag()<<endl;
    Vector3 temp(0.,0.,0.);
    return temp;    
  }
  else
  { 
    return thr_axis; 
  }
}

float Sper(std::vector<Vector4>& ptl, Vector3& Bthr)
{
 double B_tvec[3]={Bthr.x(),Bthr.y(),Bthr.z()};
 int ntrk = ptl.size();
 if(ntrk<1) return -7.;
 double p_trk[ntrk][3];
 double sper=-2.;
 for(int i=0;i<ntrk;i++)
 {
      p_trk[i][0] = ptl[i].px();
      p_trk[i][1] = ptl[i].py();
      p_trk[i][2] = ptl[i].pz();
 }
 spherp(ntrk, p_trk, &sper, B_tvec);
 if(sper >=0. && sper<=1.)
   return sper;
 else
   return -7.;
}

SuperFoxWolfram::SuperFoxWolfram(){
        for (int i=0; i<5; i++) sum[i] = 0;

}
SuperFoxWolfram::~SuperFoxWolfram(){}
void SuperFoxWolfram::fill( std::vector<Vector4>& plist,std::vector<Vector4>& qlist){
  if(plist.size()*qlist.size()>0){
  for ( std::vector<Vector4>::iterator it1=plist.begin();
        it1!=plist.end(); it1++ ){
  for ( std::vector<Vector4>::iterator it2=qlist.begin();
        it2!=qlist.end(); it2++ ){
	
        Vector4&  p = *it1;
        Vector4&  q = *it2;
        Vector3 pvec = p;
        Vector3 qvec = q;
	
        double mag = pvec.mag() * qvec.mag();
        double costh = pvec.dot(qvec) / mag;
        double cost2 = costh * costh;
        sum[0] += mag;
        sum[1] += mag * costh;
        sum[2] += mag * LGNDR2(costh);
        sum[3] += mag * LGNDR3(costh);
        sum[4] += mag * LGNDR4(costh);
	}
	}
   }
   else for(int i=0;i<5;i++) sum[i]=0.;
}
double SuperFoxWolfram::R(int i){
    if( i < 0 || i > 4 || sum[0] == 0. ) return -7.;
    else {
        double sfwMoment = sum[i]/sum[0];	
	if(abs(sfwMoment)>100.) sfwMoment=-7.;
	return sfwMoment;
    }
}


/*.......................................................................
 . FWJET2 - Subroutine to compute event sphericity 
 . 
 . Inputs    : NTRK -- number of tracks in momentum array 
 .           : PTRK -- array dimensioned (3,NTRK) containing 3-momenta 
 .           :         of tracks in the event 
 .           : EBEAM -- beam energy of the event 
 . Outputs   : H0 -- 0th moment of the event 
 .           : R1 -- H1/H0 of the event 
 .           : R2 -- H2/H0 of the event 
 .           : R3 -- H3/H0 of the event 
 .           : R4 -- H4/H0 of the event 
 . 
 . Calls     : LEGN(X) -- INTERNAL functions for Legendre polynomials 
 . Called    : Anyone who wants the Fox/Wolfram moments of the event 
 . Author    : F. Morrow  09/02/89  14.47.58 
 .             Transfer to C by kfjack 2/9,1999
 . 
 . Detailed description 
 .     FWJET2 determines the Fox/Wolfram moments of the event, up to 
 .     order 5. 
 .......................................................................
*/

void fwjet2(int NTrk, float *PTrk, float EBeam, 
           float *H0, float *R1, float *R2, float *R3, float *R4)
{
	int i, j;
	float sss ,cxang, ppari, pparj, delh0;
	float h, h1, h2, h3, h4;	

	if (NTrk >= 2) {
	
    		if (EBeam < 4.0 || EBeam > 6.5) EBeam = 5.27;
		sss = EBeam*EBeam*4;
		
		h=h1=h2=h3=h4=0;
    		
    		for (i=0;i<NTrk;i++) {
			ppari = sqrt(PTrk[i*3]*PTrk[i*3]+
			             PTrk[i*3+1]*PTrk[i*3+1]+
				     PTrk[i*3+2]*PTrk[i*3+2]);
		
		for (j=i;j<NTrk;j++) {
		
		if (i!=j) {
			pparj = sqrt(PTrk[j*3]*PTrk[j*3]+
			             PTrk[j*3+1]*PTrk[j*3+1]+
				     PTrk[j*3+2]*PTrk[j*3+2]);

			cxang = (PTrk[i*3]*PTrk[j*3]+ 
		        	 PTrk[i*3+1]*PTrk[j*3+1]+
				 PTrk[i*3+2]*PTrk[j*3+2])/(ppari * pparj);
	
			delh0 = 2*ppari*pparj;
		
			h += delh0;
			h1 += delh0 * cxang;
			h2 += delh0 * LGNDR2(cxang);
			h3 += delh0 * LGNDR3(cxang);
			h4 += delh0 * LGNDR4(cxang);
	    	} else {
			delh0 = ppari*ppari;
			h += delh0;
			h1 += delh0;
			h2 += delh0;
			h3 += delh0;
			h4 += delh0;
	    	}
		}
		}
    	}

    	if (h > 1e-7f) {
		*H0 = h / sss;
		*R1 = h1 / h;
		*R2 = h2 / h;
		*R3 = h3 / h;
		*R4 = h4 / h;
    	}else {
		*H0=*R1=*R2=*R3=*R4=-1.0;
	}
}
/*

  int  thrust(int n_trk, double Ptrk[][3], doubel thr, double thr_vec[3], int Itrk[]);

C. THRUST - Routine to determine the THRUST of an event
C.
C. Inputs    : NTRK -- number of tracks in momentum array
C.           : PTRK -- array dimensioned (3,NTRK) containing 3-momenta
C.           :         of tracks in the event
C.           :
C. Outputs   : THR -- the thrust of the event
C.           : TVEC -- the thrust axis of the event
C.           : ITRK -- 0,1 -- set to "1" if the particle is assigned to
C.           :         the jet lying along the given thrust axis (TVEC)
C.
C. Called    : Anyone wishing to compute the thrust of an event
C. Author    : F. Morrow  10/02/89  11.46.36
   C Version : Yuan Chao  04/20/99
C.
C. Detailed description
C.      Routine that finds the THRUST axis, as defined by the
C.      principal axis method, by running through all combinations
C.      of pairs of particles in an event.  The plane containing each
C.      of the pair of particles is determined and used to partition
C.      the event into two "jets".  The two particles are added to
C.      either "jet" and the thrust for each partioning is determined.
C.      The combination resulting in the maximum thrust determines the
C.      THRUST of the event.
C.      This algorithm is approximate and is formerly SUBROUTINE THRST2.
C.......................................................................
*/
/* Subroutine */ 
int thrust(int ntrk, double ptrk[][3], double* thr, double tvec[3], int itrk[]){

    /* Local variables */
    int nmax;
    double pcos, pdot;
    int npls;
    double psqi, psqj, pprp[3], ptst;
    double psqij, p1, p2, p3, pt[50][3], pjet[3], denmax;
    double rnumer, psqmax, pi1, pi2, pi3, pj1, pj2, pj3, alp, apt, psq;
    double pij1, pij2, pij3;
    int iup1, iup2, iup3;

    /* Function Body */
    *thr = -1.;
    tvec[0] = 0.;
    tvec[1] = 0.;
    tvec[2] = 0.;

    for (int i = 1; i < ntrk; ++i) 
    {
	itrk[i] = 0;
    }

    if (ntrk <= 50) {
	nmax = ntrk;
	npls = nmax;
	pt[npls][0] = 0.;
	pt[npls][1] = 0.;
	pt[npls][2] = 0.;

	for (int i = 0; i < nmax; ++i) {
	    pt[i][0] = ptrk[i][0];
	    pt[i][1] = ptrk[i][1];
	    pt[i][2] = ptrk[i][2];
	    pt[npls][0] += pt[i][0];
	    pt[npls][1] += pt[i][1];
	    pt[npls][2] += pt[i][2];	    
	}

	ptst = sqrt( pt[npls][0] * pt[npls][0] + pt[npls][1] * pt[npls][1] + pt[npls][2] * pt[npls][2]);

	if (ptst >= 1e-7) {
	    pt[nmax][0] = -pt[nmax][0] / (float)2.;
	    pt[nmax][1] = -pt[nmax][1] / (float)2.;
	    pt[nmax][2] = -pt[nmax][2] / (float)2.;
	}
    }

    if (nmax > 2 && nmax <= 51) {
	psqmax = 0.;

/* ---THE FOLLOWING TWO DO-LOOPS RUN THROUGH ALL COMBINATIONS */
/*   OF PAIRS OF PARTICLES */

	iup1 = 0;
	for (int i = 0; i <= nmax; ++i) {
	    if (i == nmax) {
		iup1 = 1;
	    }
	    
	    iup2 = 0;
	    
	    for (int j = 0; j <= nmax; ++j) {
		if (i != j) {
		    if (j == nmax) {
			iup2 = 1;
		    }
		    p1 = p2 = p3 = 0.; //---GET NORMAL TO PLANE

		    pprp[0] = pt[i][1] * pt[j][2] - pt[i][2] * pt[j][1];
		    pprp[1] = pt[i][2] * pt[j][0] - pt[i][0] * pt[j][2];
		    pprp[2] = pt[i][0] * pt[j][1] - pt[i][1] * pt[j][0];
			    
		    if (iup1 != 1 && iup2 != 1){
			pdot = pt[nmax][0] * pprp[0] + pt[nmax][1] *
				 pprp[1] + pt[nmax][2] * pprp[2];
		    }

		    if (pdot >= 0. || iup1 == 1 || iup2 == 1) {
			/* ---SELECT THE RIGHT SIDE OF THE PLANE */
			iup3 = 0;

			for (int k = 0; k <= nmax; ++k) {
			    if (k != i && k != j) {
				pdot = pt[k][0] * pprp[0] + pt[k][1]
					 * pprp[1] + pt[k][2] * pprp[2];
					 
				if (pdot >= 0.) {
				    if (k == nmax) {
					iup3 = 1;
				    }
				    p1 += pt[k][0];
				    p2 += pt[k][1];
				    p3 += pt[k][2];
				}
			    }
			}

			pi1 = p1 + pt[i][0];
			pi2 = p2 + pt[i][1];
			pi3 = p3 + pt[i][2];
			pj1 = p1 + pt[j][0];
			pj2 = p2 + pt[j][1];
			pj3 = p3 + pt[j][2];
			pij1 = pi1 + pt[j][0];
			pij2 = pi2 + pt[j][1];
			pij3 = pi3 + pt[j][2];

			psq = p1 * p1 + p2 * p2 + p3 * p3;
			psqi = pi1 * pi1 + pi2 * pi2 + pi3 * pi3;
			psqj = pj1 * pj1 + pj2 * pj2 + pj3 * pj3;
			psqij = pij1 * pij1 + pij2 * pij2 + pij3 * pij3;
			
			/* ---GET MAXIMUM MOMENTUM, ASSIGN PLANAR PARTICLES */
			if (psq > psqmax) {
			    if (iup3 == 1) {
				psqmax = psq;
				pjet[0] = p1;
				pjet[1] = p2;
				pjet[2] = p3;
			    }
			}
			if (psqi > psqmax) {
			    if (iup3 == 1 || iup1 == 1) {
				psqmax = psqi;
				pjet[0] = pi1;
				pjet[1] = pi2;
				pjet[2] = pi3;
			    }
			}
			if (psqj > psqmax) {
			    if (iup3 == 1 || iup2 == 1) {
				psqmax = psqj;
				pjet[0] = pj1;
				pjet[1] = pj2;
				pjet[2] = pj3;
			    }
			}
			if (psqij > psqmax) {
			    psqmax = psqij;
			    pjet[0] = pij1;
			    pjet[1] = pij2;
			    pjet[2] = pij3;
			}
		    }
		}
	    }
	}

/* ---EFFECT FINAL THRUST RESULTS */
	apt = sqrt(psqmax);
	tvec[0] = pjet[0] / apt;
	tvec[1] = pjet[1] / apt;
	tvec[2] = pjet[2] / apt;
	denmax = 0.;
	rnumer = 0.;
	for (int i= 0; i< ntrk; ++i) {
	    alp = sqrt(pt[i][0]*pt[i][0]+pt[i][1]*pt[i][1]+pt[i][2]*pt[i][2]);
	    
	    pcos = (pjet[0] * pt[i][0] + pjet[1] * pt[i][1] + pjet[2] * pt[i][2]) / apt;
		    
	    if (pcos > 0.) {
		itrk[i] = 1;
	    }
	    rnumer += fabs(pcos);
	    denmax += alp;
	}
	*thr = rnumer / denmax;

    } else if (nmax == 2) {
	*thr = (float)1.;

	p1 = sqrt( pt[0][0]* pt[0][0] + pt[0][1] * pt[0][1] + pt[0][2] * pt[0][2]);
	p2 = sqrt( pt[1][0]* pt[1][0] + pt[1][1] * pt[1][1] + pt[1][2] * pt[1][2]);
	if (p1 > p2) {
	    tvec[0] = pt[0][0] / p1;
	    tvec[1] = pt[0][1] / p1;
	    tvec[2] = pt[0][2] / p1;
	    itrk[0] = 1;
	} else {
	    tvec[0] = pt[1][0] / p2;
	    tvec[1] = pt[1][1] / p2;
	    tvec[2] = pt[1][2] / p2;
	    itrk[0] = 1;
	}
    }
/*
cout << "----------------------------------------" <<endl;

    for(int i=0;i<ntrk;i++){
cout << "Input " << i << " = (" << ptrk[i][0] <<" , " 
                                << ptrk[i][1] <<" , " 
                                << ptrk[i][2] <<" )" <<endl;
   } 

cout << "p thrst = (" << tvec[0] <<" , " 
                      << tvec[1] <<" , " 
                      << tvec[2] <<" )" <<endl<<endl;
*/
  
    return 0;
} /* thrust_ */


/* ....................................................................... */
/* . */
/* . SPHERP - Subroutine to compute event spherocity */
/* . */
/* . Inputs    : NTRK -- number of tracks in momentum array */
/* .           : PTRK -- array dimensioned (3,NTRK) containing 3-momenta */
/* .           :         of the tracks in the event */
/* .           : JET  --  direction be used */
/* .           :          as a given axis to calculate Sper */
/* . Outputs   : SPHER -- SPHEROCITY of the event */
/* . */
/* . Called    : Anyone who wants the spherocity of a set of tracks. */
/* . Author    : F. Morrow  08/02/89   9.01.34 */
/* .           : Modified by C.H. Wang 04/12/1996 */
// . C Version : Yuan Chao 04/20/99
/* . */
/* . Detailed description */
/* .     SPHERP determines the spherocity of an event using */
/*              ( SUM ( ABS(PT) ) ) */
/*        S = ------------------------- */
/*              ( SUM ( ABS(P) ) ) */
/* WHERE PT IS THE MOMENTUM COMPONENT PERPENDICULAR TO JET(3). */
/*     THE SUM IN THE NUMERATOR DOES NOT INCLUDE PARTICLES FALL WITHIN */
/*     45 DEGREE FORWARD AND BACKWARD CONES W.R.T. JET DIRECTION. */
/*    THE SUM IN THE DENOMINATOR INCLUDE ALL PARTICLES EXCEPT PARTICLE WITH*/
/* .    DIRECTION JET. */

/*      TAKING THE AXIS TO BE GIVEN BY A PARTICLE DIRECTION, SAY PJ, */
/*   WE HAVE (FOR PARTICLE K) THE TRANSVERSE MOMENTUM (PTK) GIVEN BY: */
/*        PTK**2 = ( PK*SINE(IK) )**2 */
/*               = PK**2 * (1 - COS(IK)**2) */
/*               = PK**2 - PK**2 * ( (PJ.DOT.PK)/(PJ*PK) )**2 */
/*               = PK**2 - ( (PJ.DOT.PK)/PJ )**2 */
/* ....................................................................... */


/* Subroutine */ 
int spherp(int ntrk, double ptrk[][3], double *sper, double jet[3])
{

    /* Local variables */
    double pmod, jetp;
    double spmod, ptmin, pl, pt1, cut, psq, xpt1;


    /* Function Body */
    if (ntrk <= 1) {
        *sper = -1;
	return -1;
    }

    ptmin = 1e5;
    spmod = 0.;
    pt1 = 0.;
    cut = sqrt(2.) / 2.;
    jetp = sqrt(jet[0] * jet[0] + jet[1] * jet[1] + jet[2] * jet[2]);

    for (int i = 0; i < ntrk; ++i) {
	psq = ptrk[i][0] * ptrk[i][0] + ptrk[i][1] * ptrk[i][1] + ptrk[i][2] * ptrk[i][2];
	pmod = sqrt(psq);
	pl = (jet[0] * ptrk[i][0] + jet[1] * ptrk[i][1] + jet[2]
		 * ptrk[i][2]) / (pmod * jetp);

	spmod += pmod;	
/* ........................................................... */
/* ....... SELECT PARTOCLE OUTSIDE THE CONE .................. */
/* ........................................................... */
	if (fabs(pl) < cut) {
	    xpt1 = psq * (1 - pl * pl);
	    pt1 += sqrt((fabs(xpt1)));

	}
    }

//    spmod -= jetp;
    *sper = pt1 / spmod;

} /* spherp_ */


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
