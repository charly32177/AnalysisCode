/*
	File : ExKFitter.cc (The main source file)
	Description :	Experimental Extended Kinematic Fitter
		A set of c++ class developed for general vertex/mass 
		constrain fitting with the supports of multi-constrains 
		and simultaneous minimization.
	        	      
	Author : Kai-Feng Chen, National Taiwan University
	Version : 1.32
*/	

#include "ExKFitter.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

double ExKF_MASS_ARRAY[5] = {0.000511,0.105658,0.139570,0.493677,0.938272}; 
	
ExKFitterParticle::ExKFitterParticle() 
{
	m_Momentum = HepLorentzVector(0.,0.,0.,0.);
	m_Position = HepPoint3D(0.,0.,0.);
	m_IsErrAvailable = 0;
	m_IsPositionAvailable = 0;
	m_IsVirual = 0;
	m_Charge = 0.;
	m_Mass = 0.;
	m_Bfield = ExKF_DEF_B_FIELD;
	m_R = 0.;

	m_IsHelixAvailable = 0;
}
	
ExKFitterParticle::ExKFitterParticle(HepLorentzVector Momentum,
        		  	     HepPoint3D       Position,
			  	     HepSymMatrix     ErrMomentumPosition,
				     double	      Charge) 
{
	m_Momentum = Momentum;
	m_Position = Position;
	m_ErrMomentumPosition = ErrMomentumPosition;
	m_IsErrAvailable = 1;
	m_IsPositionAvailable = 1;
	m_IsVirual = 0;
	m_Charge = Charge;
	m_Mass = Momentum.mag();
	m_Bfield = ExKF_DEF_B_FIELD;
	m_R = 0.;

	m_IsHelixAvailable = 0;
}	

ExKFitterParticle::ExKFitterParticle(Mdst_charged& charged,
			  	     int 	   pid)
{
	HepPoint3D pivot; 
	HepVector helix(5);
	HepSymMatrix error(5,0);

	GetHelixParameters(charged,pid,pivot,helix,error);

	Helix helixclass(pivot, helix, error);

	m_Momentum = helixclass.momentum(0.0, ExKF_MASS_ARRAY[pid], 
	        			 m_Position, m_ErrMomentumPosition);
	              
	m_Pivot = pivot;
	m_Helix = helix;
	m_ErrHelix = error;   
	m_IsHelixAvailable = 1;       
	              
	m_IsErrAvailable = 1;
	m_IsPositionAvailable = 1;
	m_IsVirual = 0;
	m_Charge = charged.charge();  
	m_Mass = ExKF_MASS_ARRAY[pid];
	m_Bfield = ExKF_DEF_B_FIELD;
	m_R = 0.;
}		

ExKFitterParticle::ExKFitterParticle(Mdst_vee2&	Vee2,
				     double	Charge)
{
	HepPoint3D pivot;
	HepVector helix(5);
	HepSymMatrix error(5,0);

	int massHyp1 = 2, massHyp2 = 2;
	switch(Vee2.kind()) {
		case 0: massHyp1 = 0; massHyp2 = 0; break;
		case 1: massHyp1 = 2; massHyp2 = 2; break;
		case 2: massHyp1 = 4; massHyp2 = 2; break;
		case 3: massHyp1 = 2; massHyp2 = 4; break;
	}
	int pid = (Charge>0 ? massHyp1:massHyp2);

	if (Vee2.daut())
		GetHelixParameters(Vee2,Charge,pivot,helix,error);
	else
		GetHelixParameters((Charge>0 ? Vee2.chgd(0):Vee2.chgd(1)),pid,pivot,helix,error);

	Helix helixclass(pivot, helix, error);

	m_Momentum = helixclass.momentum(0.0, ExKF_MASS_ARRAY[pid], 
	        			 m_Position, m_ErrMomentumPosition);
	              
	m_Pivot = pivot;
	m_Helix = helix;
	m_ErrHelix = error;   
	m_IsHelixAvailable = 1;       
	              
	m_IsErrAvailable = 1;
	m_IsPositionAvailable = 1;
	m_IsVirual = 0;
	m_Charge = Charge;    
	m_Mass = ExKF_MASS_ARRAY[pid];
	m_Bfield = ExKF_DEF_B_FIELD;
	m_R = 0.;
}

ExKFitterParticle::ExKFitterParticle(Mdst_gamma&   gamma,
			  	     HepPoint3D    Vertex,
			  	     HepSymMatrix  ErrVertex)
{
	Mdst_ecl& ecl = gamma.ecl(); 

	HepSymMatrix error(3,1);
	error[0][0] = ecl.error(0);
	error[1][0] = ecl.error(1);
	error[1][1] = ecl.error(2);
	error[2][0] = ecl.error(3);
	error[2][1] = ecl.error(4);
	error[2][2] = ecl.error(5);	

	GammaParticle gamma_part(ecl.energy(),ecl.phi(),ecl.theta(),ecl.r(),error);

	gamma_part.vertex(Vertex,ErrVertex); 

	m_Momentum = gamma_part.momentumEnergy();
	m_Position = gamma_part.vertex();

	m_ErrMomentumPosition = HepSymMatrix(7,0);

	m_ErrMomentumPosition.sub(1,gamma_part.errorMomentumEnergy());
	m_ErrMomentumPosition.sub(5,ErrVertex); 

	m_IsErrAvailable = 1;
	m_IsPositionAvailable = 1;
	m_IsVirual = 0;
	m_Charge = 0.;        
	m_Mass = 0.;
	m_Bfield = ExKF_DEF_B_FIELD;
	m_R = ecl.r();

	m_IsHelixAvailable = 0;
}			

ExKFitterParticle::ExKFitterParticle(Mdst_gamma& gamma)
{
	Mdst_ecl& ecl = gamma.ecl(); 

	HepPoint3D Vertex(0.,0.,0.);
	HepSymMatrix ErrVertex(3,0);

	HepSymMatrix error(3,1);
	error[0][0] = ecl.error(0);
	error[1][0] = ecl.error(1);
	error[1][1] = ecl.error(2);
	error[2][0] = ecl.error(3);
	error[2][1] = ecl.error(4);
	error[2][2] = ecl.error(5);	

	GammaParticle gamma_part(ecl.energy(),ecl.phi(),ecl.theta(),ecl.r(),error);

	gamma_part.vertex(Vertex,ErrVertex); 

	m_Momentum = gamma_part.momentumEnergy();
	m_Position = gamma_part.vertex();

	m_ErrMomentumPosition = HepSymMatrix(7,0);

	m_ErrMomentumPosition.sub(1,gamma_part.errorMomentumEnergy());
	m_ErrMomentumPosition.sub(5,ErrVertex); 

	m_IsErrAvailable = 1;
	m_IsPositionAvailable = 0;
	m_IsVirual = 0;
	m_Charge = 0.;        
	m_Mass = 0.;
	m_Bfield = ExKF_DEF_B_FIELD;
	m_R = ecl.r();

	m_IsHelixAvailable = 0;
}			

ExKFitterParticle & ExKFitterParticle::operator = (const ExKFitterParticle &p) // Assigment
{ 
	m_Momentum = p.m_Momentum;
	m_Position = p.m_Position;
	m_ErrMomentumPosition = p.m_ErrMomentumPosition;
	m_IsErrAvailable = p.m_IsErrAvailable;
	m_IsVirual = p.m_IsVirual;  
	m_Charge = p.m_Charge;    
	m_Bfield = p.m_Bfield;
	m_Pivot = p.m_Pivot;
	m_Helix = p.m_Helix;
	m_ErrHelix = p.m_ErrHelix;
	m_IsHelixAvailable = p.m_IsHelixAvailable;
	m_Vertex = p.m_Vertex;        
	m_ParticleList = p.m_ParticleList;
	m_IsPositionAvailable = p.m_IsPositionAvailable;
	m_R = p.m_R;  
	return *this; 
}	     		     		  
	
void ExKFitterParticle::GetHelixParameters( Mdst_charged& charged,
				 	    int 	  pid,
			 	 	    HepPoint3D&   pivot,
			 	 	    HepVector& 	  helix,
        		 	 	    HepSymMatrix& error)
{
	Mdst_trk_fit& trk_fit = charged.trk().mhyp(pid); 

	pivot = HepPoint3D(trk_fit.pivot_x(),
	        	   trk_fit.pivot_y(),
	        	   trk_fit.pivot_z());
	               
	for(int i=0;i<5;i++) helix[i]=trk_fit.helix(i);

	error[0][0] = trk_fit.error(0);
	error[1][0] = trk_fit.error(1); 
	error[1][1] = trk_fit.error(2);
	error[2][0] = trk_fit.error(3);
	error[2][1] = trk_fit.error(4);
	error[2][2] = trk_fit.error(5);
	error[3][0] = trk_fit.error(6);
	error[3][1] = trk_fit.error(7);
	error[3][2] = trk_fit.error(8);
	error[3][3] = trk_fit.error(9);
	error[4][0] = trk_fit.error(10);
	error[4][1] = trk_fit.error(11);
	error[4][2] = trk_fit.error(12);
	error[4][3] = trk_fit.error(13);
	error[4][4] = trk_fit.error(14); 
}					    

void ExKFitterParticle::GetHelixParameters( Mdst_vee2&	  Vee2,
			 		    double	  charge,
					    HepPoint3D&   pivot,
					    HepVector& 	  helix,
        				    HepSymMatrix& error) 
{
	Mdst_vee_daughters& vee2dau = Vee2.daut();
	              
	pivot = HepPoint3D(Vee2.vx(),Vee2.vy(),Vee2.vz());

	if (charge>0) {
					     
		for(int i=0;i<5;i++) helix[i]=vee2dau.helix_p(i);

		error[0][0] = vee2dau.error_p(0);
		error[1][0] = vee2dau.error_p(1); 
		error[1][1] = vee2dau.error_p(2);
		error[2][0] = vee2dau.error_p(3);
		error[2][1] = vee2dau.error_p(4);
		error[2][2] = vee2dau.error_p(5);
		error[3][0] = vee2dau.error_p(6);
		error[3][1] = vee2dau.error_p(7);
		error[3][2] = vee2dau.error_p(8);
		error[3][3] = vee2dau.error_p(9);
		error[4][0] = vee2dau.error_p(10);
		error[4][1] = vee2dau.error_p(11);
		error[4][2] = vee2dau.error_p(12);
		error[4][3] = vee2dau.error_p(13);
		error[4][4] = vee2dau.error_p(14); 
	}else {
		
		for(int i=0;i<5;i++) helix[i]=vee2dau.helix_m(i);

		error[0][0] = vee2dau.error_m(0);
		error[1][0] = vee2dau.error_m(1); 
		error[1][1] = vee2dau.error_m(2);
		error[2][0] = vee2dau.error_m(3);
		error[2][1] = vee2dau.error_m(4);
		error[2][2] = vee2dau.error_m(5);
		error[3][0] = vee2dau.error_m(6);
		error[3][1] = vee2dau.error_m(7);
		error[3][2] = vee2dau.error_m(8);
		error[3][3] = vee2dau.error_m(9);
		error[4][0] = vee2dau.error_m(10);
		error[4][1] = vee2dau.error_m(11);
		error[4][2] = vee2dau.error_m(12);
		error[4][3] = vee2dau.error_m(13);
		error[4][4] = vee2dau.error_m(14); 
	}
}

double ExKFitterParticle::TrackChisq( HepLorentzVector 	Momentum,
        	   		      HepPoint3D	Position)
{  
	if (!m_IsHelixAvailable) return 0.;   
	              
	Helix hel(Position,Momentum.vect(),m_Charge);
	hel.pivot(m_Pivot);

	HepVector diff = m_Helix - hel.a();

	int ifail;
	HepSymMatrix InvEa = m_ErrHelix.inverse(ifail);
	if (ifail!=0) {
		fprintf(stderr,"[ExKFitterParticle] Error : Cannot get inverse matrix in function TrackChisq().\n");
		return -1.;
	}

	return (diff.T() * InvEa * diff)(1,1);
}

double ExKFitterParticle::TrackXi()
{  
	if (!m_IsHelixAvailable) return 0.;   

	Helix hel(m_Pivot, m_Helix, m_ErrHelix);

	HepLorentzVector momentum;
	HepPoint3D position;
	HepSymMatrix error;   

	momentum = hel.momentum(0.0, m_Mass,position, error);
	        		      
	position.setZ(m_Position.z());

	return TrackChisq(momentum,position);
}

double ExKFitterParticle::TrackZeta()
{  
	if (!m_IsHelixAvailable) return 0.;   

	Helix hel(m_Pivot, m_Helix, m_ErrHelix);

	HepLorentzVector momentum;
	HepPoint3D position;
	HepSymMatrix error;   

	momentum = hel.momentum(0.0, m_Mass, 
	        		position, error);
	        					      
	position.setZ(m_Position.z());        
	momentum.setPz(m_Momentum.pz());

	double scale; 
	scale = sqrt(pow(m_Momentum.px(),2)+pow(m_Momentum.py(),2))/
	  	sqrt(pow(  momentum.px(),2)+pow(  momentum.py(),2));
	momentum.setPx(momentum.px()*scale);
	momentum.setPy(momentum.py()*scale);  

	return TrackChisq(momentum,position);
}

void ExKFitterParticle::LinkParticle(ExKFitterParticle* Particle) 
{
	int flag = 0;
	for(unsigned int j=0;j<m_ParticleList.size();j++) {  
		if (Particle==m_ParticleList[j]) flag = 1;
	}
	if (!flag) {
		m_ParticleList.push_back(Particle);
			    
		if (Particle->m_IsPositionAvailable) m_IsPositionAvailable = 1;
		m_IsVirual = 1;
	}
}

void ExKFitterParticle::LinkParticleRecursively(ExKFitterParticle* Particle)
{
	if (Particle->m_IsVirual) {
		for(unsigned int j=0;j<Particle->m_ParticleList.size();j++) 
			LinkParticleRecursively(Particle->m_ParticleList[j]);     
		return;
	}

	int flag = 0;
	for(unsigned int j=0;j<m_ParticleList.size();j++) {  
		if (Particle==m_ParticleList[j]) flag = 1;
	}
	if (!flag) {
		m_ParticleList.push_back(Particle);
			    
		if (Particle->m_IsPositionAvailable) m_IsPositionAvailable = 1;
		m_IsVirual = 1;
	}
}

void ExKFitterParticle::LinkVertex(ExKFitterVertex* Vertex) 
{
	m_Vertex = Vertex;
	m_IsVirual = 1;       
	if (Vertex->m_IsErrAvailable) m_IsPositionAvailable = 1;
}

int ExKFitterParticle::CheckRelation(ExKFitterParticle* Particle) 
{
	if (Particle==this) return 1;
	if (!m_IsVirual) return 0;

	for(unsigned int i=0;i<m_ParticleList.size();i++) { 
		ExKFitterParticle *par = m_ParticleList[i];
			    
		if (par->CheckRelation(Particle)) return 1; 
	}
	return 0;
}

int ExKFitterParticle::CheckRelation(ExKFitterVertex* Vertex) 
{
	if (!m_IsVirual) return 0;
	if (Vertex==m_Vertex) return 1;

	for(unsigned int i=0;i<m_ParticleList.size();i++) { 
		ExKFitterParticle *par = m_ParticleList[i];
			    
		if (par->CheckRelation(Vertex)) return 1; 
	}
	return 0;
}

HepMatrix ExKFitterParticle::RelationMatrix(ExKFitterParticle* Particle) 
{
	HepMatrix J(7,7,0); // px,py,pz,x,y,z,m       

	if (Particle==this) { // = ident      
		for(int i=1;i<=7;i++) J(i,i) = 1.;
		return J;
	}
	        	      
	if (!m_IsVirual) {    // = null       
		return J;
	}     
	        	      
	for(unsigned int i=0;i<m_ParticleList.size();i++) { 
		ExKFitterParticle *par = m_ParticleList[i];
			    
		if (par->CheckRelation(Particle)) {
			    
		    if (par->m_IsPositionAvailable) {
				
		 	 double px = par->m_Momentum.px();
		 	 double py = par->m_Momentum.py();
		 	 double pz = par->m_Momentum.pz();
		 	 double m  = par->m_Mass;
		 	 double e  = sqrt(px*px + py*py + pz*pz + m*m);
		 	 double z  = par->m_Position.z();
		 	 double vz = m_Vertex->m_Vertex.z();
		 	     		 
		 	 double al = par->alpha();
		 	     	 
		 	 //double PX  = \sum px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz))
		 	     	 
		 	 double dPX_dpx = cos(al/pz*(z-vz));
		 	 double dPX_dpy = sin(al/pz*(z-vz));		 
		 	 double dPX_dpz = + px * sin(al/pz*(z-vz)) * al/pow(pz,2)*(z-vz) 
		 	   		  - py * cos(al/pz*(z-vz)) * al/pow(pz,2)*(z-vz);
		 	     		 
		 	 double dPX_dx = 0.;
		 	 double dPX_dy = 0.;
		 	 double dPX_dz = - px * sin(al/pz*(z-vz)) * al/pz 
		 	   		 + py * cos(al/pz*(z-vz)) * al/pz;
		 	 double dPX_dm = 0.;							 
		 	     	 
		 	 J(1,1) = dPX_dpx;	 // del px
		 	 J(1,2) = dPX_dpy;	 // del py
		 	 J(1,3) = dPX_dpz;	 // del pz
		 	 J(1,4) = dPX_dx;	 // del x
		 	 J(1,5) = dPX_dy;	 // del y
		 	 J(1,6) = dPX_dz;	 // del z
		 	 J(1,7) = dPX_dm;	 // del m
		 	     		 
		 	 //double PY  = \sum py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz))
		 	     		 
		 	 double dPY_dpx = - sin(al/pz*(z-vz));
		 	 double dPY_dpy =   cos(al/pz*(z-vz));  	
		 	 double dPY_dpz = + py * sin(al/pz*(z-vz)) * al/pow(pz,2)*(z-vz)	
		 	   + px * cos(al/pz*(z-vz)) * al/pow(pz,2)*(z-vz);
		 	     		 
		 	 double dPY_dx = 0.;
		 	 double dPY_dy = 0.;
		 	 double dPY_dz = - py * sin(al/pz*(z-vz)) * al/pz 
		 	   - px * cos(al/pz*(z-vz)) * al/pz;
		 	 double dPY_dm = 0.;		 
		 	     				 
		 	 J(2,1) = dPY_dpx;	 // del px
		 	 J(2,2) = dPY_dpy;	 // del py
		 	 J(2,3) = dPY_dpz;	 // del pz
		 	 J(2,4) = dPY_dx;	 // del x
		 	 J(2,5) = dPY_dy;	 // del y
		 	 J(2,6) = dPY_dz;	 // del z
		 	 J(2,7) = dPY_dm;	 // del m	 
		 	     		 
		 	 //double PZ  = \sum pz
		 	     		 
		 	 J(3,1) = 0.;	 // del px
		 	 J(3,2) = 0.;	 // del py
		 	 J(3,3) = 1.;	 // del pz
		 	 J(3,4) = 0.;	 // del x
		 	 J(3,5) = 0.;	 // del y
		 	 J(3,6) = 0.;	 // del z
		 	 J(3,7) = 0.;	 // del m
		 	     		 
		 	 //double X  = vx
		 	 //double Y  = vy
		 	 //double Z  = vz
		 	 
		 	 //double M  = sqrt(E^2 - PX^2 -PY^2 -PZ^2)			 
		 	 //double PX = \sum px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz))
		 	 //double PY = \sum py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz))
		 	 //double PZ = \sum pz
		 	 //double E  = \sum sqrt(px*px + py*py + pz*pz + m*m)
		 	     		 
		 	 double PX = m_Momentum.px();
		 	 double PY = m_Momentum.py();
		 	 double PZ = m_Momentum.pz();
		 	 double M  = m_Mass;
		 	 double E  = sqrt(PX*PX + PY*PY + PZ*PZ + M*M);
		 	     		 
		 	 double dM_dE  =   E/M;
		 	 double dM_dPX = -PX/M;
		 	 double dM_dPY = -PY/M;
		 	 double dM_dPZ = -PZ/M;
		 	     		 
		 	 double dE_dpx = px/e;
		 	 double dE_dpy = py/e;
		 	 double dE_dpz = pz/e;  		  
		 	 double dE_dx  = 0.;
		 	 double dE_dy  = 0.;
		 	 double dE_dz  = 0.;
		 	 double dE_dm  = m/e;
		 	     		 
		 	 J(7,1) = dM_dE*dE_dpx + dM_dPX*dPX_dpx + dM_dPY*dPY_dpx;		 // del px
		 	 J(7,2) = dM_dE*dE_dpy + dM_dPX*dPX_dpy + dM_dPY*dPY_dpy;		 // del py
		 	 J(7,3) = dM_dE*dE_dpz + dM_dPX*dPX_dpz + dM_dPY*dPY_dpz + dM_dPZ;	 // del pz
		 	 J(7,4) = dM_dE*dE_dx + dM_dPX*dPX_dx + dM_dPY*dPY_dx;  		 // del x
		 	 J(7,5) = dM_dE*dE_dy + dM_dPX*dPX_dy + dM_dPY*dPY_dy;  		 // del y
		 	 J(7,6) = dM_dE*dE_dz + dM_dPX*dPX_dz + dM_dPY*dPY_dz;  		 // del z
		 	 J(7,7) = dM_dE*dE_dm + dM_dPX*dPX_dm + dM_dPY*dPY_dm;  		 // del m
		 	     		 
		 	 // Others are null matrix
		 	     	 
		 	 J = J*par->RelationMatrix(Particle);	 
				    
		    }else {
				
			// Only used for photons 
			// x,y,z,R are constants
			    	    
			double px = par->m_Momentum.px();
			double py = par->m_Momentum.py();
			double pz = par->m_Momentum.pz();
			double p  = sqrt(px*px + py*py + pz*pz);
			double m  = par->m_Mass;
			double e  = sqrt(px*px + py*py + pz*pz + m*m);
			double x  = par->m_Position.x();
			double y  = par->m_Position.y();
			double z  = par->m_Position.z();
			double vx = m_Vertex->m_Vertex.x();
			double vy = m_Vertex->m_Vertex.y();
			double vz = m_Vertex->m_Vertex.z();
			double R  = par->m_R;
			    		
			double Rprimesq = pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)+pow(R,2) +
			  		  2*R/p*(px*(x-vx) + py*(y-vy) + pz*(z-vz));
			double Rprime	= sqrt(Rprimesq);
			    		
			//double PX = \sum px*R/Rprime + (x-vx)/Rprime * p
			    		
			double dRprime_dpx = R/Rprime*((x-vx)/p - px*(px*(x-vx) + py*(y-vy) + pz*(z-vz))/pow(p,3));
			double dRprime_dpy = R/Rprime*((y-vy)/p - py*(px*(x-vx) + py*(y-vy) + pz*(z-vz))/pow(p,3));
			double dRprime_dpz = R/Rprime*((z-vz)/p - pz*(px*(x-vx) + py*(y-vy) + pz*(z-vz))/pow(p,3));
			    		
			double dPX_dpx = R/Rprime - px*R/Rprimesq*dRprime_dpx + (x-vx)/Rprime*px/p - (x-vx)*p/Rprimesq*dRprime_dpx;
			double dPX_dpy =	  - px*R/Rprimesq*dRprime_dpy + (x-vx)/Rprime*py/p - (x-vx)*p/Rprimesq*dRprime_dpy;
			double dPX_dpz =	  - px*R/Rprimesq*dRprime_dpz + (x-vx)/Rprime*pz/p - (x-vx)*p/Rprimesq*dRprime_dpz;
			double dPX_dx = 0.;
			double dPX_dy = 0.;
			double dPX_dz = 0.;
			double dPX_dm = 0.;							
			    	
			J(1,1) = dPX_dpx;	// del px
			J(1,2) = dPX_dpy;	// del py
			J(1,3) = dPX_dpz;	// del pz
			J(1,4) = dPX_dx;	// del x
			J(1,5) = dPX_dy;	// del y
			J(1,6) = dPX_dz;	// del z
			J(1,7) = dPX_dm;	// del m
			    		
			//double PY = \sum py*R/Rprime + (y-vy)/Rprime * p
			    		
			double dPY_dpx =	  - py*R/Rprimesq*dRprime_dpx + (y-vy)/Rprime*px/p - (y-vy)*p/Rprimesq*dRprime_dpx;
			double dPY_dpy = R/Rprime - py*R/Rprimesq*dRprime_dpy + (y-vy)/Rprime*py/p - (y-vy)*p/Rprimesq*dRprime_dpy;
			double dPY_dpz =	  - py*R/Rprimesq*dRprime_dpz + (y-vy)/Rprime*pz/p - (y-vy)*p/Rprimesq*dRprime_dpz;
			double dPY_dx = 0.;
			double dPY_dy = 0.;
			double dPY_dz = 0.;
			double dPY_dm = 0.;	
			    				
			J(2,1) = dPY_dpx;	// del px
			J(2,2) = dPY_dpy;	// del py
			J(2,3) = dPY_dpz;	// del pz
			J(2,4) = dPY_dx;	// del x
			J(2,5) = dPY_dy;	// del y
			J(2,6) = dPY_dz;	// del z
			J(2,7) = dPY_dm;	// del m	
			    		
			//double PZ = \sum pz*R/Rprime + (z-vz)/Rprime * p
			    		
			double dPZ_dpx =	  - pz*R/Rprimesq*dRprime_dpx + (z-vz)/Rprime*px/p - (z-vz)*p/Rprimesq*dRprime_dpx;
			double dPZ_dpy =	  - pz*R/Rprimesq*dRprime_dpy + (z-vz)/Rprime*py/p - (z-vz)*p/Rprimesq*dRprime_dpy;
			double dPZ_dpz = R/Rprime - pz*R/Rprimesq*dRprime_dpz + (z-vz)/Rprime*pz/p - (z-vz)*p/Rprimesq*dRprime_dpz;
			double dPZ_dx = 0.;
			double dPZ_dy = 0.;
			double dPZ_dz = 0.;
			double dPZ_dm = 0.;
			    		
			J(3,1) = dPZ_dpx;	// del px
			J(3,2) = dPZ_dpy;	// del py
			J(3,3) = dPZ_dpz;	// del pz
			J(3,4) = dPZ_dx;	// del x
			J(3,5) = dPZ_dy;	// del y
			J(3,6) = dPZ_dz;	// del z
			J(3,7) = dPZ_dm;	// del m
			    		
			//double X  = vx
			//double Y  = vy
			//double Z  = vz
			    		
			//double M  = sqrt(E^2 - PX^2 -PY^2 -PZ^2)			
			//double PX = \sum px*R/Rprime + (x-vx)/Rprime * p
			//double PY = \sum py*R/Rprime + (y-vy)/Rprime * p
			//double PZ = \sum pz*R/Rprime + (z-vz)/Rprime * p
			//double E  = \sum sqrt(px*px + py*py + pz*pz + m*m)
			    		
			double PX = m_Momentum.px();
			double PY = m_Momentum.py();
			double PZ = m_Momentum.pz();
			double M  = m_Mass;
			double E  = sqrt(PX*PX + PY*PY + PZ*PZ + M*M);
			    		
			double dM_dE  =   E/M;
			double dM_dPX = -PX/M;
			double dM_dPY = -PY/M;
			double dM_dPZ = -PZ/M;
			    		
			double dE_dpx = px/e;
			double dE_dpy = py/e;
			double dE_dpz = pz/e;			 
			double dE_dx  = 0.;
			double dE_dy  = 0.;
			double dE_dz  = 0.;
			double dE_dm  = m/e;
			    		
			J(7,1) = dM_dE*dE_dpx + dM_dPX*dPX_dpx + dM_dPY*dPY_dpx;		// del px
			J(7,2) = dM_dE*dE_dpy + dM_dPX*dPX_dpy + dM_dPY*dPY_dpy;		// del py
			J(7,3) = dM_dE*dE_dpz + dM_dPX*dPX_dpz + dM_dPY*dPY_dpz + dM_dPZ;	// del pz
			J(7,4) = dM_dE*dE_dx + dM_dPX*dPX_dx + dM_dPY*dPY_dx;			// del x
			J(7,5) = dM_dE*dE_dy + dM_dPX*dPX_dy + dM_dPY*dPY_dy;			// del y
			J(7,6) = dM_dE*dE_dz + dM_dPX*dPX_dz + dM_dPY*dPY_dz;			// del z
			J(7,7) = dM_dE*dE_dm + dM_dPX*dPX_dm + dM_dPY*dPY_dm;			// del m
			    		
			// Others are null matrix
			    	
			J = J*par->RelationMatrix(Particle);	
		    }				    
		}
	}
	return J;
}

HepMatrix ExKFitterParticle::RelationMatrix(ExKFitterVertex* Vertex) 
{
	HepMatrix J(7,3,0);

	if (!m_IsVirual) return J;
	if (Vertex==m_Vertex) {
		
		for(unsigned int i=0;i<m_ParticleList.size();i++) {
				    
		    ExKFitterParticle *par = m_ParticleList[i];
				    
		    if (par->m_IsPositionAvailable) {

			double px = par->m_Momentum.px();
			double py = par->m_Momentum.py();
			double pz = par->m_Momentum.pz();
			double z  = par->m_Position.z();
			double vz = m_Vertex->m_Vertex.z();
			    	
			double al = par->alpha();
			    	
			//double PX  = \sum px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz))
			    		
			double dPX_dvz = + px * sin(al/pz*(z-vz)) * al/pz 
			  		 - py * cos(al/pz*(z-vz)) * al/pz;					
			    		
			J(1,1) += 0.;		// del vx
			J(1,2) += 0.;		// del vy
			J(1,3) += dPX_dvz;	// del vz
			    		
			//double PY  = \sum py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz))
			    		
			double dPY_dvz = + py * sin(al/pz*(z-vz)) * al/pz
			  		 + px * cos(al/pz*(z-vz)) * al/pz;
			    		
			J(2,1) += 0.;		// del vx
			J(2,2) += 0.;		// del vy
			J(2,3) += dPY_dvz;	// del vz
			    		
			//double PZ  = \sum pz
			    		
			//double X  = vx			    								
			J(4,1) = 1.;	// del vx
			J(4,2) = 0.;	// del vy
			J(4,3) = 0.;	// del vz
			    		
			//double Y  = vy			 	
			J(5,1) = 0.;	// del vx
			J(5,2) = 1.;	// del vy
			J(5,3) = 0.;	// del vz
			    		
			//double Z  = vz			    		
			J(6,1) = 0.;	// del vx
			J(6,2) = 0.;	// del vy
			J(6,3) = 1.;	// del vz
			    		
			//double M  = sqrt(E^2 - PX^2 -PY^2 -PZ^2)			
			//double PX = \sum px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz))
			//double PY = \sum py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz))
			//double PZ = \sum pz
			//double E  = \sum sqrt(px*px + py*py + pz*pz + m*m)
			    		
			double PX = m_Momentum.px();
			double PY = m_Momentum.py();
			double M  = m_Mass;
			    		
			double dM_dPX = -PX/M;
			double dM_dPY = -PY/M;
			    		
			J(7,1) += 0.;	  // del vx
			J(7,2) += 0.;	  // del vy
			J(7,3) += dM_dPX*dPX_dvz + dM_dPY*dPY_dvz;	 // del vz
		    }else {
				    
			double px = par->m_Momentum.px();
			double py = par->m_Momentum.py();
			double pz = par->m_Momentum.pz();
			double p  = sqrt(px*px + py*py + pz*pz);
			double x  = par->m_Position.x();
			double y  = par->m_Position.y();
			double z  = par->m_Position.z();
			double vx = m_Vertex->m_Vertex.x();
			double vy = m_Vertex->m_Vertex.y();
			double vz = m_Vertex->m_Vertex.z();
			double R  = par->m_R;
			    	
			double Rprimesq = pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)+pow(R,2) +
			  	  	  2*R/p*(px*(x-vx) + py*(y-vy) + pz*(z-vz));
			double Rprime	= sqrt(Rprimesq);
			    		
			double dRprime_dvx = -((x-vx) + R*px/p)/Rprime;
			double dRprime_dvy = -((y-vy) + R*py/p)/Rprime;
			double dRprime_dvz = -((z-vz) + R*pz/p)/Rprime;
			    		
			//double PX = \sum px*R/Rprime + (x-vx)/Rprime * p	
			    										
			double dPX_dvx = -p/Rprime - (px*R + (x-vx)*p)/Rprimesq*dRprime_dvx;
			double dPX_dvy =	   - (px*R + (x-vx)*p)/Rprimesq*dRprime_dvy;
			double dPX_dvz =	   - (px*R + (x-vx)*p)/Rprimesq*dRprime_dvz;				
			    		
			J(1,1) += dPX_dvx;	// del vx
			J(1,2) += dPX_dvy;	// del vy
			J(1,3) += dPX_dvz;	// del vz
			    		
			//double PY = \sum py*R/Rprime + (y-vy)/Rprime * p		     
			double dPY_dvx =	   - (py*R + (y-vy)*p)/Rprimesq*dRprime_dvx;
			double dPY_dvy = -p/Rprime - (py*R + (y-vy)*p)/Rprimesq*dRprime_dvy;
			double dPY_dvz =	   - (py*R + (y-vy)*p)/Rprimesq*dRprime_dvz;
			    		
			J(2,1) += dPY_dvx;	// del vx
			J(2,2) += dPY_dvy;	// del vy
			J(2,3) += dPY_dvz;	// del vz
			    		
			//double PZ = \sum pz*R/Rprime + (z-vz)/Rprime * p
			double dPZ_dvx =	   - (pz*R + (z-vz)*p)/Rprimesq*dRprime_dvx;
			double dPZ_dvy =	   - (pz*R + (z-vz)*p)/Rprimesq*dRprime_dvy;
			double dPZ_dvz = -p/Rprime - (pz*R + (z-vz)*p)/Rprimesq*dRprime_dvz;
			    		
			J(3,1) += dPZ_dvx;	// del vx
			J(3,2) += dPZ_dvy;	// del vy
			J(3,3) += dPZ_dvz;	// del vz
			    		
			//double X  = vx			    								
			J(4,1) = 1.;	// del vx
			J(4,2) = 0.;	// del vy
			J(4,3) = 0.;	// del vz
			    		
			//double Y  = vy			    		
			J(5,1) = 0.;	// del vx
			J(5,2) = 1.;	// del vy
			J(5,3) = 0.;	// del vz
			    		
			//double Z  = vz			    		
			J(6,1) = 0.;	// del vx
			J(6,2) = 0.;	// del vy
			J(6,3) = 1.;	// del vz
			    		
			//double M  = sqrt(E^2 - PX^2 -PY^2 -PZ^2)			
			//double PX = \sum px*R/Rprime + (x-vx)/Rprime * p
			//double PY = \sum py*R/Rprime + (y-vy)/Rprime * p
			//double PZ = \sum pz*R/Rprime + (z-vz)/Rprime * p
			//double E  = \sum sqrt(px*px + py*py + pz*pz + m*m)
			    		
			double PX = m_Momentum.px();
			double PY = m_Momentum.py();
			double PZ = m_Momentum.pz();
			double M  = m_Mass;
			    		
			double dM_dPX = -PX/M;
			double dM_dPY = -PY/M;
			double dM_dPZ = -PZ/M;
			    		
			J(7,1) += dM_dPX*dPX_dvx;	// del vx
			J(7,2) += dM_dPY*dPY_dvy;	// del vy
			J(7,3) += dM_dPZ*dPZ_dvz;	// del vz
		    }
		}
			    
		return J;
	}

	for(unsigned int i=0;i<m_ParticleList.size();i++) { 
		ExKFitterParticle *par = m_ParticleList[i];
							    
		if (par->CheckRelation(Vertex)) return RelationMatrix(par) * par->RelationMatrix(Vertex); 
	}
	return J;
}

void ExKFitterParticle::Update()
{
	if (!m_IsVirual) {    
		m_Momentum.setE(sqrt(m_Momentum.rho()*m_Momentum.rho() + m_Mass*m_Mass));
		return;
	}

	if (m_Vertex == NULL || m_ParticleList.size()==0) {
		fprintf(stderr,"[ExKFitterParticle] Error : 4-vector cannot be updated - lack of daughters/vertex information.\n");
		return; 
	}
	// if (par->m_IsPositionAvailable) {
	//    double PX  = \sum px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz))
	//    double PY  = \sum py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz))
	//    double PZ  = \sum pz
	// }else {
	//    double PX  = \sum px*R/Rprime + (x-vx)/Rprime * p
	//    double PY  = \sum py*R/Rprime + (y-vy)/Rprime * p
	//    double PZ  = \sum pz*R/Rprime + (z-vz)/Rprime * p
	// }

	// double X  = vx
	// double Y  = vy
	// double Z  = vz

	double PX,PY,PZ,E;
	PX = PY = PZ = E = 0.;
	double CHARGE = 0.;
	
	while (1) {
	    int flag = 0;   
	    for(unsigned int j=0;j<m_ParticleList.size();j++) {
		if (!m_ParticleList[j]->m_IsPositionAvailable &&
		     m_ParticleList[j]->m_IsVirual) {
			ExKFitterParticle *sub = m_ParticleList[j];
						    
			m_ParticleList.erase(m_ParticleList.begin()+j);
			LinkParticleRecursively(sub);
						    
			flag = 1;
			break;
		}
	    }
	    if (!flag) break;
	}

	for(unsigned int i=0;i<m_ParticleList.size();i++) {
		ExKFitterParticle *par = m_ParticleList[i];				
			    
		par->Update();				
			    
		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double p  = sqrt(px*px + py*py + pz*pz);
		double e  = sqrt(px*px + py*py + pz*pz + par->m_Mass*par->m_Mass);
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double z  = par->m_Position.z();
		double vx = m_Vertex->m_Vertex.x();
		double vy = m_Vertex->m_Vertex.y();
		double vz = m_Vertex->m_Vertex.z(); 
		double R  = par->m_R;
			    
		double al = par->alpha();
			    
		double Rprimesq = pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)+pow(R,2) +
		  		  2*R/p*(px*(x-vx) + py*(y-vy) + pz*(z-vz));
		double Rprime	= sqrt(Rprimesq);
			    
		if (par->m_IsPositionAvailable) {
			PX += px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz));
			PY += py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz));
			PZ += pz; 
			E  += e;
		}else {
			PX += px*R/Rprime + (x-vx)/Rprime * p;
			PY += py*R/Rprime + (y-vy)/Rprime * p;
			PZ += pz*R/Rprime + (z-vz)/Rprime * p;
			E  += e;
		}
			    
		CHARGE += par->m_Charge;	    
	}

	m_Momentum.setPx(PX);
	m_Momentum.setPy(PY);
	m_Momentum.setPz(PZ);
	m_Momentum.setE(E);

	m_Position = m_Vertex->m_Vertex;
	m_Charge = CHARGE;
	m_Mass = m_Momentum.mag();
}

void ExKFitterParticle::Bfield(double d)
{
	if (!m_IsVirual) {    
		m_Bfield = d;
		return;
	}

	for(unsigned int i=0;i<m_ParticleList.size();i++) { 
		ExKFitterParticle *par = m_ParticleList[i];
			    
		par->Bfield(d); 	    
	}
}

ExKFitterVertex::ExKFitterVertex() 
{
	m_Vertex = HepPoint3D(0.,0.,0.);
	m_ErrVertex = HepSymMatrix(3,0);
	m_IsErrAvailable = 0;
}
	
ExKFitterVertex::ExKFitterVertex(HepPoint3D   Vertex) 
{
	m_Vertex    = Vertex;
	m_ErrVertex = HepSymMatrix(3,0);
	m_IsErrAvailable = 0;
}
	
ExKFitterVertex::ExKFitterVertex(HepPoint3D   Vertex,
			         HepSymMatrix ErrVertex) 
{
	m_Vertex    = Vertex;
	m_ErrVertex = ErrVertex;
	m_IsErrAvailable = 1;
}

ExKFitterVertex & ExKFitterVertex::operator = (const ExKFitterVertex &p)
{
	m_Vertex = p.m_Vertex;
	m_ErrVertex = p.m_ErrVertex;
	m_IsErrAvailable = p.m_IsErrAvailable; return *this;
}
	
ExKFitterMass::ExKFitterMass()
{
	m_Mass = 0.;
	m_ErrMass = 0.;
	m_IsErrAvailable = 0;
}

ExKFitterMass::ExKFitterMass(double Mass)
{
	m_Mass = Mass;
	m_ErrMass = 0.;
	m_IsErrAvailable = 0;
}

ExKFitterMass::ExKFitterMass(double Mass, double ErrMass)
{
	m_Mass = Mass;
	m_ErrMass = ErrMass;
	m_IsErrAvailable = 1;
}

ExKFitterMass & ExKFitterMass::operator = (const ExKFitterMass &p)
{
	m_Mass = p.m_Mass;
	m_ErrMass = p.m_ErrMass;
	m_IsErrAvailable = p.m_IsErrAvailable; 
	return *this;
}

ExKFitterConstrain & ExKFitterConstrain::operator = (const ExKFitterConstrain &p)
{
	m_Vertex = p.m_Vertex;
	m_Mass = p.m_Mass;
	m_ErrorFlag = p.m_ErrorFlag;  
	m_ParticleList = p.m_ParticleList; 
	return *this;
}
	
void ExKFitterConstrain::LinkParticle(ExKFitterParticle* Particle) 
{
	int flag = 0;
	for(unsigned int j=0;j<m_ParticleList.size();j++) { 
		if (Particle==m_ParticleList[j]) flag = 1;
	}
	if (!flag) m_ParticleList.push_back(Particle);
}

void ExKFitterConstrain::LinkParticleRecursively(ExKFitterParticle* Particle)
{
	if (Particle->m_IsVirual) {
		for(unsigned int j=0;j<Particle->m_ParticleList.size();j++) 
			LinkParticleRecursively(Particle->m_ParticleList[j]);	  
		return;
	}

	int flag = 0;
	for(unsigned int j=0;j<m_ParticleList.size();j++) { 
		if (Particle==m_ParticleList[j]) flag = 1;
	}
	if (!flag) m_ParticleList.push_back(Particle);
}
	
HepVector ExKFitterConstrain::Hvc() // = constrain vector H(a,v) = d
{
	HepVector d(m_ParticleList.size()*2);

	for(unsigned int i=0;i<m_ParticleList.size();i++) {
			    
		ExKFitterParticle *par = m_ParticleList[i];
			    
		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double z  = par->m_Position.z();
		double vx = m_Vertex->m_Vertex.x();
		double vy = m_Vertex->m_Vertex.y();
		double vz = m_Vertex->m_Vertex.z();
		double pt = sqrt(px*px+py*py);
			    
		double W = px*(vx-x) + py*(vy-y);
			    
		if (par->m_Charge!=0.) {
			double al = par->alpha(); 
			double W_ov_ptsq = W/pow(pt,2); 	  

			if (al*W_ov_ptsq<-1. || al*W_ov_ptsq>1.) {m_ErrorFlag = ExKF_MATHERROR;break;}

			double asin_alxW_ov_ptsq = asin(al*W_ov_ptsq);

			d(i*2+1) = px*(vy-y) - py*(vx-x) - al/2. *(pow(vx-x,2) + pow(vy-y,2)); 
			d(i*2+2) = (vz-z)*pt - pz*pt/al * asin_alxW_ov_ptsq;
		}else {
			d(i*2+1) = px*(vy-y) - py*(vx-x); 
			d(i*2+2) = (vz-z)*pt - pz/pt * W; 
		}
	}

	return d;
}

HepVector ExKFitterConstrain::Hmc() // = constrain vector H(a,v) = d
{
	HepVector d(1);

	// if (par->m_IsPositionAvailable) {
	//    double PX  = \sum px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz))
	//    double PY  = \sum py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz))
	//    double PZ  = \sum pz
	// }else {
	//    double PX  = \sum px*R/Rprime + (x-vx)/Rprime * p;
	//    double PY  = \sum py*R/Rprime + (y-vy)/Rprime * p;
	//    double PZ  = \sum pz*R/Rprime + (z-vz)/Rprime * p;
	// }

	double PX,PY,PZ,E;
	PX = PY = PZ = E = 0.;

	for(unsigned int i=0;i<m_ParticleList.size();i++) {
		
		ExKFitterParticle *par = m_ParticleList[i];
			    
		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double p  = sqrt(px*px + py*py + pz*pz);
		double e  = sqrt(px*px + py*py + pz*pz + par->m_Mass*par->m_Mass);
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double z  = par->m_Position.z();
		double vx = m_Vertex->m_Vertex.x();
		double vy = m_Vertex->m_Vertex.y();
		double vz = m_Vertex->m_Vertex.z();
		double R  = par->m_R;	    
			    
		double al = par->alpha();
			    
		double Rprimesq = pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)+pow(R,2) +
				  2*R/p*(px*(x-vx) + py*(y-vy) + pz*(z-vz));
		double Rprime	= sqrt(Rprimesq);
			    
		if (par->m_IsPositionAvailable) {
			PX += px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz));
			PY += py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz));
			PZ += pz; 
			E  += e;
		}else {
			PX += px*R/Rprime + (x-vx)/Rprime * p;
			PY += py*R/Rprime + (y-vy)/Rprime * p;
			PZ += pz*R/Rprime + (z-vz)/Rprime * p;
			E  += e;
		}   
	}

	d(1) = E*E - PX*PX - PY*PY - PZ*PZ - m_Mass->m_Mass*m_Mass->m_Mass;

	return d;
}
	
HepMatrix ExKFitterConstrain::delHvc_delParticle(ExKFitterParticle* Particle) 
{
	HepMatrix J(m_ParticleList.size()*2,7,0);

	unsigned int index;
	for(index=0;index<m_ParticleList.size();index++) {
		
		if (!m_ParticleList[index]->CheckRelation(Particle)) continue;

		HepMatrix K(2,7,0);

		ExKFitterParticle* par = m_ParticleList[index];

		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double z  = par->m_Position.z();
		double vx = m_Vertex->m_Vertex.x();
		double vy = m_Vertex->m_Vertex.y();
		double vz = m_Vertex->m_Vertex.z();
		double pt = sqrt(px*px+py*py);

		double W = px*(vx-x) + py*(vy-y);

		double dpt_dpx = px/pt;
		double dpt_dpy = py/pt;

		double dW_dpx = (vx-x);
		double dW_dpy = (vy-y);
		double dW_dx  = -px;
		double dW_dy  = -py;

		if (par->m_Charge!=0.) {
			          
			double al = par->alpha(); 
			double W_ov_ptsq = W/pow(pt,2);

			if (al*W_ov_ptsq<-1. || al*W_ov_ptsq>1.) {m_ErrorFlag = ExKF_MATHERROR;return J;}

			double asin_alxW_ov_ptsq = asin(al*W_ov_ptsq);  	   

			double dW_ov_ptsq_dpx = dW_dpx/pow(pt,2) + W/pow(pt,4)*(-2.*px);
			double dW_ov_ptsq_dpy = dW_dpy/pow(pt,2) + W/pow(pt,4)*(-2.*py);       
			double dW_ov_ptsq_dx  = dW_dx/pow(pt,2);
			double dW_ov_ptsq_dy  = dW_dy/pow(pt,2);

			double dasin_alxW_ov_ptsq_dpx = al*dW_ov_ptsq_dpx / sqrt(1.-pow(al*W_ov_ptsq,2)) ;
			double dasin_alxW_ov_ptsq_dpy = al*dW_ov_ptsq_dpy / sqrt(1.-pow(al*W_ov_ptsq,2)) ;
			double dasin_alxW_ov_ptsq_dx  = al*dW_ov_ptsq_dx  / sqrt(1.-pow(al*W_ov_ptsq,2)) ;
			double dasin_alxW_ov_ptsq_dy  = al*dW_ov_ptsq_dy  / sqrt(1.-pow(al*W_ov_ptsq,2)) ;

			//H(1) = px*(vy-y) - py*(vx-x) - al/2. *(pow(vx-x,2) + pow(vy-y,2)); 

			K(1,1)  = (vy-y);	  	  // del px
			K(1,2)  = -(vx-x);		  // del py
			K(1,3)  = 0.;			  // del pz
			K(1,4)  =   py + al*(vx-x);	  // del x
			K(1,5)  = - px + al*(vy-y);	  // del y
			K(1,6)  = 0.;			  // del z
			K(1,7)  = 0.;			  // del m		  

			//H(2) = (vz-z)*pt - pz*pt/al * asin_alxW_ov_ptsq;

			K(2,1)  = (vz-z)*dpt_dpx - pz*dpt_dpx/al * asin_alxW_ov_ptsq - pz*pt/al * dasin_alxW_ov_ptsq_dpx; // del px			      
			K(2,2)  = (vz-z)*dpt_dpy - pz*dpt_dpy/al * asin_alxW_ov_ptsq - pz*pt/al * dasin_alxW_ov_ptsq_dpy; // del py
			K(2,3)  = - pt/al * asin_alxW_ov_ptsq;  							  // del pz
			K(2,4)  = - pz*pt/al * dasin_alxW_ov_ptsq_dx;							  // del x
			K(2,5)  = - pz*pt/al * dasin_alxW_ov_ptsq_dy;							  // del y
			K(2,6)  = - pt; 										  // del z
			K(2,7)  = 0.;											  // del m
		}else {
			
			//H(1) = px*(vy-y) - py*(vx-x); 

			K(1,1)  = (vy-y); 	  // del px
			K(1,2)  = -(vx-x);	  // del py
			K(1,3)  = 0.;		  // del pz
			K(1,4)  =   py;   	  // del x
			K(1,5)  = - px;   	  // del y
			K(1,6)  = 0.;		  // del z
			K(1,7)  = 0.;		  // del m
			          
			//H(2) = (vz-z)*pt - pz/pt * W;

			K(2,1)  = (vz-z)*dpt_dpx + pz * W / pow(pt,2) * dpt_dpx - pz/pt * dW_dpx;  // del px
			K(2,2)  = (vz-z)*dpt_dpy + pz * W / pow(pt,2) * dpt_dpy - pz/pt * dW_dpy;  // del py
			K(2,3)  = - W/pt;							   // del pz
			K(2,4)  = - pz/pt * dW_dx;						   // del x
			K(2,5)  = - pz/pt * dW_dy;						   // del y
			K(2,6)  = - pt; 							   // del z
			K(2,7)  = 0.;								   // del m				   
		}

		K = K*par->RelationMatrix(Particle);

		J.sub(index*2+1,1,K);

	}

	return J;     
}  

HepMatrix ExKFitterConstrain::delHmc_delParticle(ExKFitterParticle* Particle) 
{	
	double PX,PY,PZ,E;
	PX = PY = PZ = E = 0.;

	for(unsigned int i=0;i<m_ParticleList.size();i++) {
		
		ExKFitterParticle *par = m_ParticleList[i];
			    
		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double p  = sqrt(px*px + py*py + pz*pz);
		double e  = sqrt(px*px + py*py + pz*pz + par->m_Mass*par->m_Mass);
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double z  = par->m_Position.z();
		double vx = m_Vertex->m_Vertex.x();
		double vy = m_Vertex->m_Vertex.y();
		double vz = m_Vertex->m_Vertex.z(); 
		double R  = par->m_R;
			    
		double al = par->alpha();
			    
		double Rprimesq = pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)+pow(R,2) +
		  		  2*R/p*(px*(x-vx) + py*(y-vy) + pz*(z-vz));
		double Rprime	= sqrt(Rprimesq);
			    
		if (par->m_IsPositionAvailable) {
			PX += px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz));
			PY += py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz));
			PZ += pz; 
			E  += e;
		}else {
			PX += px*R/Rprime + (x-vx)/Rprime * p;
			PY += py*R/Rprime + (y-vy)/Rprime * p;
			PZ += pz*R/Rprime + (z-vz)/Rprime * p;
			E  += e;
		}
	}

	HepMatrix J(1,7,0);

	unsigned int index;
	for(index=0;index<m_ParticleList.size();index++) {
		
		if (!m_ParticleList[index]->CheckRelation(Particle)) continue;  		    
			    
		HepMatrix K(1,7,0);

		ExKFitterParticle* par = m_ParticleList[index];

		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double p  = sqrt(px*px + py*py + pz*pz);
		double m  = par->m_Mass;
		double e  = sqrt(px*px + py*py + pz*pz + m*m);
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double z  = par->m_Position.z();
		double vx = m_Vertex->m_Vertex.x();
		double vy = m_Vertex->m_Vertex.y();
		double vz = m_Vertex->m_Vertex.z();
		double R  = par->m_R;

		double al = par->alpha();	
			    
		if (par->m_IsPositionAvailable) {   
			          
			//double PX  = \sum px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz))
			//double PY  = \sum py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz))
			//double PZ  = \sum pz
			//double E   = \sum sqrt(px*px + py*py + pz*pz + m*m)
			          
			double dE_dpx = px/e;
			double dE_dpy = py/e;
			double dE_dpz = pz/e;			   
			double dE_dx  = 0.;
			double dE_dy  = 0.;
			double dE_dz  = 0.;
			double dE_dm  = m/e;
			          
			double dPX_dpx = cos(al/pz*(z-vz));
			double dPX_dpy = sin(al/pz*(z-vz));
			double dPX_dpz = + px*sin(al/pz*(z-vz)) * al*(z-vz)/pow(pz,2)	  
			  		 - py*cos(al/pz*(z-vz)) * al*(z-vz)/pow(pz,2); 			   
			double dPX_dx  = 0.;
			double dPX_dy  = 0.;
			double dPX_dz  = - px*sin(al/pz*(z-vz)) * al/pz 
			  		 + py*cos(al/pz*(z-vz)) * al/pz;
			        		   
			double dPY_dpx = - sin(al/pz*(z-vz));
			double dPY_dpy = + cos(al/pz*(z-vz));
			double dPY_dpz = + py*sin(al/pz*(z-vz)) * al*(z-vz)/pow(pz,2)	  
			  		 + px*cos(al/pz*(z-vz)) * al*(z-vz)/pow(pz,2); 			   
			double dPY_dx  = 0.;
			double dPY_dy  = 0.;
			double dPY_dz  = - py*sin(al/pz*(z-vz)) * al/pz 
			  		 - px*cos(al/pz*(z-vz)) * al/pz;
			        		   
			double dPZ_dpx = 0.;
			double dPZ_dpy = 0.;
			double dPZ_dpz = 1.;				   
			double dPZ_dx  = 0.;
			double dPZ_dy  = 0.;
			double dPZ_dz  = 0.;					  

			//H(1) = E*E - PX*PX - PY*PY - PZ*PZ - Mass*Mass;
			          
			K(1,1) = 2*E*dE_dpx - 2*PX*dPX_dpx - 2*PY*dPY_dpx - 2*PZ*dPZ_dpx;    // del px
			K(1,2) = 2*E*dE_dpy - 2*PX*dPX_dpy - 2*PY*dPY_dpy - 2*PZ*dPZ_dpy;    // del py
			K(1,3) = 2*E*dE_dpz - 2*PX*dPX_dpz - 2*PY*dPY_dpz - 2*PZ*dPZ_dpz;    // del pz
			K(1,4) = 2*E*dE_dx - 2*PX*dPX_dx - 2*PY*dPY_dx - 2*PZ*dPZ_dx;	     // del x
			K(1,5) = 2*E*dE_dy - 2*PX*dPX_dy - 2*PY*dPY_dy - 2*PZ*dPZ_dy;	     // del y
			K(1,6) = 2*E*dE_dz - 2*PX*dPX_dz - 2*PY*dPY_dz - 2*PZ*dPZ_dz;	     // del z 
			K(1,7) = 2*E*dE_dm;						     // del m	
		}else {
			          
			//double PX  = \sum px*R/Rprime + (x-vx)/Rprime * p;
			//double PY  = \sum py*R/Rprime + (y-vy)/Rprime * p;
			//double PZ  = \sum pz*R/Rprime + (z-vz)/Rprime * p;
			//double E   = \sum sqrt(px*px + py*py + pz*pz + m*m)
			        	  
			double Rprimesq = pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)+pow(R,2) +
					  2*R/p*(px*(x-vx) + py*(y-vy) + pz*(z-vz));
			double Rprime	= sqrt(Rprimesq);
			        	  
			double dRprime_dpx = R/Rprime*((x-vx)/p - px*(px*(x-vx) + py*(y-vy) + pz*(z-vz))/pow(p,3));
			double dRprime_dpy = R/Rprime*((y-vy)/p - py*(px*(x-vx) + py*(y-vy) + pz*(z-vz))/pow(p,3));
			double dRprime_dpz = R/Rprime*((z-vz)/p - pz*(px*(x-vx) + py*(y-vy) + pz*(z-vz))/pow(p,3));			  
			        			  
			double dE_dpx = px/e;
			double dE_dpy = py/e;
			double dE_dpz = pz/e;			   
			double dE_dx  = 0.;
			double dE_dy  = 0.;
			double dE_dz  = 0.;
			double dE_dm  = m/e;						  
			        	  
			double dPX_dpx = R/Rprime - px*R/Rprimesq*dRprime_dpx + (x-vx)/Rprime*px/p - (x-vx)*p/Rprimesq*dRprime_dpx;
			double dPX_dpy =	  - px*R/Rprimesq*dRprime_dpy + (x-vx)/Rprime*py/p - (x-vx)*p/Rprimesq*dRprime_dpy;
			double dPX_dpz =	  - px*R/Rprimesq*dRprime_dpz + (x-vx)/Rprime*pz/p - (x-vx)*p/Rprimesq*dRprime_dpz;
			double dPX_dx = 0.;
			double dPX_dy = 0.;
			double dPX_dz = 0.;
			        	  
			double dPY_dpx =	  - py*R/Rprimesq*dRprime_dpx + (y-vy)/Rprime*px/p - (y-vy)*p/Rprimesq*dRprime_dpx;
			double dPY_dpy = R/Rprime - py*R/Rprimesq*dRprime_dpy + (y-vy)/Rprime*py/p - (y-vy)*p/Rprimesq*dRprime_dpy;
			double dPY_dpz =	  - py*R/Rprimesq*dRprime_dpz + (y-vy)/Rprime*pz/p - (y-vy)*p/Rprimesq*dRprime_dpz;
			double dPY_dx = 0.;
			double dPY_dy = 0.;
			double dPY_dz = 0.;
			        	  
			double dPZ_dpx =	  - pz*R/Rprimesq*dRprime_dpx + (z-vz)/Rprime*px/p - (z-vz)*p/Rprimesq*dRprime_dpx;
			double dPZ_dpy =	  - pz*R/Rprimesq*dRprime_dpy + (z-vz)/Rprime*py/p - (z-vz)*p/Rprimesq*dRprime_dpy;
			double dPZ_dpz = R/Rprime - pz*R/Rprimesq*dRprime_dpz + (z-vz)/Rprime*pz/p - (z-vz)*p/Rprimesq*dRprime_dpz;
			double dPZ_dx = 0.;
			double dPZ_dy = 0.;
			double dPZ_dz = 0.;

			//H(1) = E*E - PX*PX - PY*PY - PZ*PZ - Mass*Mass;
			          
			K(1,1) = 2*E*dE_dpx - 2*PX*dPX_dpx - 2*PY*dPY_dpx - 2*PZ*dPZ_dpx;     // del px
			K(1,2) = 2*E*dE_dpy - 2*PX*dPX_dpy - 2*PY*dPY_dpy - 2*PZ*dPZ_dpy;     // del py
			K(1,3) = 2*E*dE_dpz - 2*PX*dPX_dpz - 2*PY*dPY_dpz - 2*PZ*dPZ_dpz;     // del pz
			K(1,4) = 2*E*dE_dx - 2*PX*dPX_dx - 2*PY*dPY_dx - 2*PZ*dPZ_dx;	      // del x
			K(1,5) = 2*E*dE_dy - 2*PX*dPX_dy - 2*PY*dPY_dy - 2*PZ*dPZ_dy;	      // del y
			K(1,6) = 2*E*dE_dz - 2*PX*dPX_dz - 2*PY*dPY_dz - 2*PZ*dPZ_dz;	      // del z 
			K(1,7) = 2*E*dE_dm;						      // del m   
		}
			    
		K = K*par->RelationMatrix(Particle);	    
			    
		J = J + K;
	}
	        	      
	return J;     
}

HepMatrix ExKFitterConstrain::delHvc_delVertex(ExKFitterVertex* Vertex) 
{
	HepMatrix J(m_ParticleList.size()*2,3,0);

	if (Vertex!=m_Vertex) {
		
		unsigned int index;
		for(index=0;index<m_ParticleList.size();index++) {
			
			if (!m_ParticleList[index]->CheckRelation(Vertex)) continue;
			        	  
			HepMatrix K(m_ParticleList.size()*2,3,0);
			ExKFitterParticle* par = m_ParticleList[index];
			        	  
			K = delHvc_delParticle(par)*par->RelationMatrix(Vertex);
			        	  
			J = J + K;		  
		}
		return J;
	}

	if (m_FollowVertex) return J;

	unsigned int index; 
	for(index=0;index<m_ParticleList.size();index++) {
			
		HepMatrix K(2,3,0);    

		ExKFitterParticle* par = m_ParticleList[index];

		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double vx = Vertex->m_Vertex.x();
		double vy = Vertex->m_Vertex.y();
		double pt = sqrt(px*px+py*py);

		double W = px*(vx-x) + py*(vy-y);

		double dW_dvx  =  px;
		double dW_dvy  =  py;

		if (par->m_Charge!=0.) {
			
			double al = par->alpha();
			double W_ov_ptsq = W/pow(pt,2); 	  

			double dW_ov_ptsq_dvx  = dW_dvx/pow(pt,2);
			double dW_ov_ptsq_dvy  = dW_dvy/pow(pt,2);

			double dasin_alxW_ov_ptsq_dvx  = al*dW_ov_ptsq_dvx  / sqrt(1.-pow(al*W_ov_ptsq,2)) ;
			double dasin_alxW_ov_ptsq_dvy  = al*dW_ov_ptsq_dvy  / sqrt(1.-pow(al*W_ov_ptsq,2)) ;

			//H(1) = px*(vy-y) - py*(vx-x) - al/2. *(pow(vx-x,2) + pow(vy-y,2)); 
			K(1,1)  = - py - al*(vx-x);	  // del vx
			K(1,2)  =   px - al*(vy-y);	  // del vy
			K(1,3)  = 0.;			  // del vz		  

			//H(2) = (vz-z)*pt - pz*pt/al * asin_alxW_ov_ptsq;
			K(2,1)  = - pz*pt/al * dasin_alxW_ov_ptsq_dvx;  // del vx
			K(2,2)  = - pz*pt/al * dasin_alxW_ov_ptsq_dvy;  // del vy
			K(2,3)  = pt;					  // del vz

		}else {
			
			//H(1) = px*(vy-y) - py*(vx-x); 
			K(1,1)  = - py;   // del vx
			K(1,2)  =   px;   // del vy
			K(1,3)  = 0.;		  // del vz
			          
			//H(2) = (vz-z)*pt - pz/pt * W;
			K(2,1)  = - pz/pt * dW_dvx;	  // del vx
			K(2,2)  = - pz/pt * dW_dvy;	  // del vy
			K(2,3)  = pt;			  // del vz
		}   

		J.sub(index*2+1,1,K);

	}
	return J;
}

HepMatrix ExKFitterConstrain::delHmc_delVertex(ExKFitterVertex* Vertex) 
{	
	HepMatrix J(1,3,0);

	if (Vertex!=m_Vertex) {
		
		unsigned int index; 
		for(index=0;index<m_ParticleList.size();index++) {
			
			if (!m_ParticleList[index]->CheckRelation(Vertex)) continue;
			          
			HepMatrix K(1,3,0);
			ExKFitterParticle* par = m_ParticleList[index];
			          
			K = delHmc_delParticle(par)*par->RelationMatrix(Vertex);
			        	  
			J = J + K;
			          
		}
		return J;
	}

	double PX,PY,PZ,E;
	PX = PY = PZ = E = 0.;

	for(unsigned int i=0;i<m_ParticleList.size();i++) {
		ExKFitterParticle *par = m_ParticleList[i];
			    
		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double p  = sqrt(px*px + py*py + pz*pz);
		double e  = sqrt(px*px + py*py + pz*pz + par->m_Mass*par->m_Mass);
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double z  = par->m_Position.z();
		double vx = m_Vertex->m_Vertex.x();
		double vy = m_Vertex->m_Vertex.y();
		double vz = m_Vertex->m_Vertex.z();
		double R  = par->m_R;	    
			    
		double al = par->alpha();
			    
		double Rprimesq = pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)+pow(R,2) +
		  		  2*R/p*(px*(x-vx) + py*(y-vy) + pz*(z-vz));
		double Rprime	= sqrt(Rprimesq);
			    
		if (par->m_IsPositionAvailable) {
			PX += px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz));
			PY += py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz));
			PZ += pz; 
			E  += e;
		}else {
			PX += px*R/Rprime + (x-vx)/Rprime * p;
			PY += py*R/Rprime + (y-vy)/Rprime * p;
			PZ += pz*R/Rprime + (z-vz)/Rprime * p;
			E  += e;
		}	    
	}	      

	unsigned int index; 
	for(index=0;index<m_ParticleList.size();index++) {	      
			    
		HepMatrix K(1,3,0);

		ExKFitterParticle* par = m_ParticleList[index];

		double px = par->m_Momentum.px();
		double py = par->m_Momentum.py();
		double pz = par->m_Momentum.pz();
		double p  = sqrt(px*px + py*py + pz*pz);
		double x  = par->m_Position.x();
		double y  = par->m_Position.y();
		double z  = par->m_Position.z();
		double vx = m_Vertex->m_Vertex.x();
		double vy = m_Vertex->m_Vertex.y();
		double vz = m_Vertex->m_Vertex.z();
		double R  = par->m_R;	    

		double al = par->alpha();
			    
		if (par->m_IsPositionAvailable) {
			          
			//double PX  = \sum px*cos(al/pz*(z-vz)) + py*sin(al/pz*(z-vz))
			//double PY  = \sum py*cos(al/pz*(z-vz)) - px*sin(al/pz*(z-vz))
			//double PZ  = \sum pz
			//double E   = \sum sqrt(px*px + py*py + pz*pz + m*m)
			        				   
			double dE_dvx  = 0.;
			double dE_dvy  = 0.;
			double dE_dvz  = 0.;
			        				   
			double dPX_dvx  = 0.;
			double dPX_dvy  = 0.;
			double dPX_dvz  = + px*sin(al/pz*(z-vz)) * al/pz 
			  		  - py*cos(al/pz*(z-vz)) * al/pz;
			           
			double dPY_dvx  = 0.;
			double dPY_dvy  = 0.;
			double dPY_dvz  = + py*sin(al/pz*(z-vz)) * al/pz 
			  	  	  + px*cos(al/pz*(z-vz)) * al/pz;
			        						   
			double dPZ_dvx  = 0.;
			double dPZ_dvy  = 0.;
			double dPZ_dvz  = 0.;						  

			//H(1) = E*E - PX*PX - PY*PY - PZ*PZ - Mass*Mass;			          
			K(1,1) = 2*E*dE_dvx - 2*PX*dPX_dvx - 2*PY*dPY_dvx - 2*PZ*dPZ_dvx; // del vx
			K(1,2) = 2*E*dE_dvy - 2*PX*dPX_dvy - 2*PY*dPY_dvy - 2*PZ*dPZ_dvy; // del vy
			K(1,3) = 2*E*dE_dvz - 2*PX*dPX_dvz - 2*PY*dPY_dvz - 2*PZ*dPZ_dvz; // del vz
		}else {
			          
			//double PX  = \sum px*R/Rprime + (x-vx)/Rprime * p;
			//double PY  = \sum py*R/Rprime + (y-vy)/Rprime * p;
			//double PZ  = \sum pz*R/Rprime + (z-vz)/Rprime * p;
			//double E   = \sum sqrt(px*px + py*py + pz*pz + m*m)
			        	  
			double Rprimesq = pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)+pow(R,2) +
			  		  2*R/p*(px*(x-vx) + py*(y-vy) + pz*(z-vz));
			double Rprime	= sqrt(Rprimesq);
			        	  
			double dRprime_dvx = -((x-vx) + R*px/p)/Rprime;
			double dRprime_dvy = -((y-vy) + R*py/p)/Rprime;
			double dRprime_dvz = -((z-vz) + R*pz/p)/Rprime;
			        				   
			double dE_dvx  = 0.;
			double dE_dvy  = 0.;
			double dE_dvz  = 0.;
			        	  
			double dPX_dvx = -p/Rprime - (px*R + (x-vx)*p)/Rprimesq*dRprime_dvx;
			double dPX_dvy =	   - (px*R + (x-vx)*p)/Rprimesq*dRprime_dvy;
			double dPX_dvz =	   - (px*R + (x-vx)*p)/Rprimesq*dRprime_dvz;				  
			        			       
			double dPY_dvx =	   - (py*R + (y-vy)*p)/Rprimesq*dRprime_dvx;
			double dPY_dvy = -p/Rprime - (py*R + (y-vy)*p)/Rprimesq*dRprime_dvy;
			double dPY_dvz =	   - (py*R + (y-vy)*p)/Rprimesq*dRprime_dvz;

			double dPZ_dvx =      	   - (pz*R + (z-vz)*p)/Rprimesq*dRprime_dvx;
			double dPZ_dvy =     	   - (pz*R + (z-vz)*p)/Rprimesq*dRprime_dvy;
			double dPZ_dvz = -p/Rprime - (pz*R + (z-vz)*p)/Rprimesq*dRprime_dvz;					    

			//H(1) = E*E - PX*PX - PY*PY - PZ*PZ - Mass*Mass;			          
			K(1,1) = 2*E*dE_dvx - 2*PX*dPX_dvx - 2*PY*dPY_dvx - 2*PZ*dPZ_dvx; // del vx
			K(1,2) = 2*E*dE_dvy - 2*PX*dPX_dvy - 2*PY*dPY_dvy - 2*PZ*dPZ_dvy; // del vy
			K(1,3) = 2*E*dE_dvz - 2*PX*dPX_dvz - 2*PY*dPY_dvz - 2*PZ*dPZ_dvz; // del vz
		}	     
			    
		J = J + K;
	}
	        	      
	return J;     
}

HepMatrix ExKFitterConstrain::delHvc_delMass(ExKFitterMass* /*Mass*/)
{
	HepMatrix J(m_ParticleList.size()*2,1,0);

	return J;
}

HepMatrix ExKFitterConstrain::delHmc_delMass(ExKFitterMass* Mass)
{
	HepMatrix J(1,1,0);

	if (Mass != m_Mass) return J;

	//H(1) = E*E - PX*PX - PY*PY - PZ*PZ - Mass*Mass;     
	J(1,1) = -2*Mass->m_Mass;     // del Mass

	return J;
}

void ExKFitterConstrain::Update()
{
	for(unsigned int i=0;i<m_ParticleList.size();i++) { 
		ExKFitterParticle *par = m_ParticleList[i];
			    
		par->Update();  	    
	}
}

void ExKFitterConstrain::Bfield(double d)
{
	for(unsigned int i=0;i<m_ParticleList.size();i++) { 
		ExKFitterParticle *par = m_ParticleList[i];
			    
		par->Bfield(d); 	    
	}
}

void ExKFitter::LinkConstrain(ExKFitterConstrain* Constrain) 
{
	int flag = 0;
	for(unsigned int j=0;j<m_ConstrainList.size();j++) {
		if (Constrain==m_ConstrainList[j]) flag = 1;
	}
	if (!flag) m_ConstrainList.push_back(Constrain);
}

ExKFitter & ExKFitter::operator = (const ExKFitter &p)
{
	m_ParticleList = p.m_ParticleList;
	m_VertexList = p.m_VertexList;
	m_MassList = p.m_MassList;    
	m_ConstrainList = p.m_ConstrainList;

	m_Max_Iteractions = p.m_Max_Iteractions;
	m_Verbose = p.m_Verbose;
	m_ErrorFlag = p.m_ErrorFlag;  
	m_Chisq = p.m_Chisq;
	m_N_Constraints = p.m_N_Constraints;  
	m_N_DegreeOfFreedom = p.m_N_DegreeOfFreedom;
	return *this;
}

void ExKFitter::LinkParticleRecursively(ExKFitterParticle* Particle)
{
	if (Particle->m_IsVirual) {
		for(unsigned int j=0;j<Particle->m_ParticleList.size();j++) 
			LinkParticleRecursively(Particle->m_ParticleList[j]);     
		return;
	}

	int flag = 0;
	for(unsigned int j=0;j<m_ParticleList.size();j++) { 
		if (Particle==m_ParticleList[j]) flag = 1;
	}
	if (!flag) m_ParticleList.push_back(Particle);
}

void ExKFitter::LinkExternalParticle(ExKFitterParticle* Particle)
{
	int flag = 0;
	for(unsigned int j=0;j<m_ExternalParticleList.size();j++) {
		if (Particle==m_ExternalParticleList[j]) flag = 1;
	}
	if (!flag) m_ExternalParticleList.push_back(Particle);
}

void ExKFitter::LinkVertex(ExKFitterVertex* Vertex) 
{
	int flag = 0;
	for(unsigned int j=0;j<m_VertexList.size();j++) {
		if (Vertex==m_VertexList[j]) flag = 1;
	}
	if (!flag) m_VertexList.push_back(Vertex);
}

void ExKFitter::LinkMass(ExKFitterMass* Mass) 
{
	int flag = 0;
	for(unsigned int j=0;j<m_MassList.size();j++) {
		if (Mass==m_MassList[j]) flag = 1;
	}
	if (!flag) m_MassList.push_back(Mass);
}

int ExKFitter::Minimize()
{
	if (PreMinimize()!=ExKF_NOERROR) return m_ErrorFlag; 

	if (m_N_Constraints-m_N_DegreeOfFreedom > 0) // N(unknown) = N(Constrains) - N(d.o.f.)
		return Minimize_w_UnknownPar();
	else 
		return Minimize_wo_UnknownPar();	     
}
	
int ExKFitter::PreMinimize()	
{
	for(unsigned int it=0;it<m_ConstrainList.size();it++) { 
		
		ExKFitterConstrain* con = m_ConstrainList[it];
			    
		if (con->m_ConstrainType == ExKF_MASSCONSTRAIN) {	    
			        	  
		    while (1) {
		    	int flag = 0;
		    	for(unsigned int i=0;i<con->m_ParticleList.size();i++) {
		    	    if (!con->m_ParticleList[i]->m_IsPositionAvailable &&
		    	    	 con->m_ParticleList[i]->m_IsVirual) {
		    	    	    ExKFitterParticle *par = con->m_ParticleList[i];
		    	    					
		    	    	    con->m_ParticleList.erase(con->m_ParticleList.begin()+i);
		    	    	    con->LinkParticleRecursively(par);
		    	    					
		    	    	    m_ExternalParticleList.push_back(par);
		    	    					
		    	    	    flag = 1;
		    	    	    break;
		    	    }
		    	}
		    	if (!flag) break;
		    }
		}
	}

	for(unsigned int it=0;it<m_ConstrainList.size();it++) {
		
		ExKFitterConstrain* con = m_ConstrainList[it];
			    
		for(unsigned int i=0;i<con->m_ParticleList.size();i++) {
		    ExKFitterParticle *par = con->m_ParticleList[i];
		  	  		      
		    while (1) {
		    	int flag = 0;	
		    	for(unsigned int j=0;j<par->m_ParticleList.size();j++) {
		    	    if (!par->m_ParticleList[j]->m_IsPositionAvailable &&
		    		 par->m_ParticleList[j]->m_IsVirual) {
		    		    ExKFitterParticle *sub = par->m_ParticleList[j];
		    						
		    		    par->m_ParticleList.erase(par->m_ParticleList.begin()+j);
		    		    par->LinkParticleRecursively(sub);
		    						
		    		    m_ExternalParticleList.push_back(sub);
		    						
		    		    flag = 1;
		    		    break;
		    	    }
		    	}
		    	if (!flag) break;
		    }
		}
	}

	for(unsigned int it=0;it<m_ConstrainList.size();it++) {
		
		ExKFitterConstrain* con = m_ConstrainList[it];
			    
		con->Update();
			    
		for(unsigned int i=0;i<con->m_ParticleList.size();i++) 
		  	LinkParticleRecursively(con->m_ParticleList[i]);
				    
		if (con->m_Vertex!=NULL && !con->m_FollowVertex) LinkVertex(con->m_Vertex);
		if (con->m_Mass!=NULL) LinkMass(con->m_Mass);
	}     

	for(unsigned int it=0;it<m_ConstrainList.size();it++) {
		
		ExKFitterConstrain* con = m_ConstrainList[it];
			    
		if (!con->CheckAvailavility()) {
			if (m_Verbose) fprintf(stderr,"[ExKFitter] Error : One of the constrains is not available.\n");
			m_ErrorFlag = ExKF_BADCONSTRAINS;
			return m_ErrorFlag;
		}
	}     

	int N_UnknownPar = 0; 
	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		if (par->m_IsPositionAvailable) {
			if (!par->m_IsErrAvailable) N_UnknownPar += 6;
		}else {
			if (!par->m_IsErrAvailable) N_UnknownPar += 3;
		}
	}     
	for(unsigned int it=0;it<m_VertexList.size();it++) {
		ExKFitterVertex* ver = m_VertexList[it];
			    
		if (!ver->m_IsErrAvailable) N_UnknownPar += 3;
	}

	m_N_Constraints = 0;
	for(unsigned int it=0;it<m_ConstrainList.size();it++) {       
		ExKFitterConstrain* con = m_ConstrainList[it];
			    
		m_N_Constraints+= con->N_Constraints();
	}     

	m_N_DegreeOfFreedom = m_N_Constraints - N_UnknownPar;

	if (m_N_DegreeOfFreedom<1) {
		m_ErrorFlag = ExKF_BADCONSTRAINS;
		return m_ErrorFlag;
	}

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		for(unsigned int it1=0;it1<m_ParticleList.size();it1++) {
			if (par != m_ParticleList[it1] &&
			    par->m_Momentum == m_ParticleList[it1]->m_Momentum &&
			    par->m_Position == m_ParticleList[it1]->m_Position) {			  
			  	if (m_Verbose) fprintf(stderr,"[ExKFitter] Warning : Two distinct tracks have the same initial condition.\n");
			}
		}
	}

	for(unsigned int it=0;it<m_VertexList.size();it++) {
		ExKFitterVertex* ver = m_VertexList[it];
			    
		for(unsigned int it1=0;it1<m_VertexList.size();it1++) {
			if (ver != m_VertexList[it1] &&
			    ver->m_Vertex == m_VertexList[it1]->m_Vertex) {				  
				if (m_Verbose) fprintf(stderr,"[ExKFitter] Warning : Two distinct vertexes have the same initial condition.\n");
			}
		}
	}     

	if (m_Max_Iteractions<0) 
	  	m_Max_Iteractions = ExKF_DEF_MAX_ITERACTIONS +
	    		(ExKF_DEF_MAX_ITERACTIONS/2)*m_ConstrainList.size();

	return ExKF_NOERROR;
}	
	
int ExKFitter::Minimize_w_UnknownPar()
{

	//Create related matrices
	int D_col;
	int E_col;
	int a_row, v_row;

	D_col = 0;
	E_col = 0;
	a_row = v_row = 0;

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		if (par->m_IsPositionAvailable) {
			if (par->m_IsErrAvailable) {D_col+=6;a_row+=6;}
			else			   {E_col+=6;v_row+=6;}
		}else {
			if (par->m_IsErrAvailable) {D_col+=3;a_row+=3;}
			else			   {E_col+=3;v_row+=3;}
		}
	}

	for(unsigned int it=0;it<m_VertexList.size();it++) {
		ExKFitterVertex* ver = m_VertexList[it];
			    
		if (ver->m_IsErrAvailable) {D_col+=3;a_row+=3;}
		else			   {E_col+=3;v_row+=3;}
	}     

	for(unsigned int it=0;it<m_MassList.size();it++) {
		ExKFitterMass* mas = m_MassList[it];
			    
		if (mas->m_IsErrAvailable) {D_col+=1;a_row+=1;}
	}
	              
	double chisq_keep   = HUGE;
	double chisq_keep_v = HUGE;
	m_Chisq = -1.0;

	int ifail = 0;
	int fit_updated = 0;

	HepVector a0(a_row),a(a_row);
	HepVector v0(v_row),v(v_row);
	HepSymMatrix Va0(a_row,0);

	int ai = 0;
	int vi = 0;

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		if (par->m_IsPositionAvailable) {
		    if (par->m_IsErrAvailable) {
			a0(ai+1) = par->m_Momentum.px();
			a0(ai+2) = par->m_Momentum.py();
			a0(ai+3) = par->m_Momentum.pz();
			a0(ai+4) = par->m_Position.x();
			a0(ai+5) = par->m_Position.y();
			a0(ai+6) = par->m_Position.z();
					
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
				Va0(ai+i,ai+j)     = par->m_ErrMomentumPosition(i,j);
				Va0(ai+3+i,ai+3+j) = par->m_ErrMomentumPosition(4+i,4+j);
			}
			}				
			for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
			  	Va0(ai+i,ai+3+j)   = par->m_ErrMomentumPosition(i,4+j);
			}
			}  
								
			ai += 6;
		    }else {
			v0(vi+1) = par->m_Momentum.px();
			v0(vi+2) = par->m_Momentum.py();
			v0(vi+3) = par->m_Momentum.pz();
			v0(vi+4) = par->m_Position.x();
			v0(vi+5) = par->m_Position.y();
			v0(vi+6) = par->m_Position.z();
			vi += 6;		
		    }
		}else {
		    if (par->m_IsErrAvailable) {
			a0(ai+1) = par->m_Momentum.px();
			a0(ai+2) = par->m_Momentum.py();
			a0(ai+3) = par->m_Momentum.pz();
					
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
			  	Va0(ai+i,ai+j)     = par->m_ErrMomentumPosition(i,j);
			}
			}							
			ai += 3;
		    }else {
			v0(vi+1) = par->m_Momentum.px();
			v0(vi+2) = par->m_Momentum.py();
			v0(vi+3) = par->m_Momentum.pz();
			vi += 3;		
		    }
		}
	}

	for(unsigned int it=0;it<m_VertexList.size();it++) {
		ExKFitterVertex* ver = m_VertexList[it];
			    
		if (ver->m_IsErrAvailable) {
			a0(ai+1) = ver->m_Vertex.x();
			a0(ai+2) = ver->m_Vertex.y();
			a0(ai+3) = ver->m_Vertex.z();
			        	  
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
			  	Va0(ai+i,ai+j)     = ver->m_ErrVertex(i,j);
			}
			}			  
			        	  
			ai += 3;
		} else {
			v0(vi+1) = ver->m_Vertex.x();
			v0(vi+2) = ver->m_Vertex.y();
			v0(vi+3) = ver->m_Vertex.z();
			vi += 3;
		}
	}

	for(unsigned int it=0;it<m_MassList.size();it++) {
		ExKFitterMass* mas = m_MassList[it];
			    
		if (mas->m_IsErrAvailable) {
			a0(ai+1) 	= mas->m_Mass;
			Va0(ai+1,ai+1) 	= mas->m_ErrMass;	  
			ai += 1;
		}
	}

	a = a0;
	v = v0;

	HepMatrix D(m_N_Constraints,D_col,0);
	HepMatrix E(m_N_Constraints,E_col,0);	      

	HepVector d(m_N_Constraints);

	HepVector da0(a_row); 
	HepVector dv0(v_row);

	HepMatrix VD(m_N_Constraints,m_N_Constraints);     
	HepMatrix VE(E_col,E_col);

	HepMatrix Lambda0(m_N_Constraints,1);
	HepMatrix Lambda(m_N_Constraints,1);

	HepVector v_keep	      = v;
	HepMatrix VE_keep	      = VE;
	HepMatrix VD_keep	      = VD;
	HepMatrix Lambda0_keep        = Lambda0;
	HepMatrix E_keep	      = E;
	HepMatrix D_keep	      = D;
	HepVector a_keep	      = a;

	HepVector v_keep_v	      = v;
	HepMatrix VE_keep_v	      = VE;
	HepMatrix VD_keep_v	      = VD;
	HepMatrix Lambda0_keep_v      = Lambda0;
	HepMatrix E_keep_v	      = E;
	HepMatrix D_keep_v	      = D;    

	//=============================================================       
	// start of minizing loop

	int it_a,it_v;
	for(it_a=0;it_a<m_Max_Iteractions;it_a++) {   // Global minizing loop
	chisq_keep_v = HUGE;				    // reset chisq_keep_v
	for(it_v=0;it_v<m_Max_Iteractions;it_v++) { // 'v' minizing loop    
		
		v0 = v; 		  // set initial vertex to be expanding point
		m_ErrorFlag = ExKF_NOERROR;	  // reset error flag		  

		int ci = 0;		  
		for(unsigned int index=0;index<m_ConstrainList.size();index++) {
			
			ExKFitterConstrain* con = m_ConstrainList[index];
				
			int di = 0;
			int ei = 0;
				
			d.sub(ci+1,con->H());
				
			for(unsigned int it=0;it<m_ParticleList.size();it++) { 
				ExKFitterParticle* par = m_ParticleList[it];					      
			
				if (par->m_IsPositionAvailable) {
				    if (par->m_IsErrAvailable) {			    
				  	D.sub(ci+1,di+1,con->delH_delParticle(par).sub(1,con->N_Constraints(),1,6));
				  	di+=6;
				    } else {
				  	E.sub(ci+1,ei+1,con->delH_delParticle(par).sub(1,con->N_Constraints(),1,6));
				  	ei+=6;
				    }
				}else {
				    if (par->m_IsErrAvailable) {			    
				  	D.sub(ci+1,di+1,con->delH_delParticle(par).sub(1,con->N_Constraints(),1,3));
				  	di+=3;
				    } else {
				  	E.sub(ci+1,ei+1,con->delH_delParticle(par).sub(1,con->N_Constraints(),1,3));
				  	ei+=3;
				    }
				}
			}
				
			for(unsigned int it=0;it<m_VertexList.size();it++) { 
				ExKFitterVertex* ver = m_VertexList[it];
				              
				if (ver->m_IsErrAvailable) {
					D.sub(ci+1,di+1,con->delH_delVertex(ver));
					di+=3;
				} else {
					E.sub(ci+1,ei+1,con->delH_delVertex(ver));
					ei+=3;
				}
			}		
				
			for(unsigned int it=0;it<m_MassList.size();it++) { 
				ExKFitterMass* mas = m_MassList[it];
			
				if (mas->m_IsErrAvailable) {
					D.sub(ci+1,di+1,con->delH_delMass(mas));
					di+=1;
				}
			}				
				
			if (con->m_ErrorFlag!=ExKF_NOERROR) m_ErrorFlag = con->m_ErrorFlag;
				
			ci+= con->N_Constraints();
		} 

		da0 = a0 - a;
		dv0 = v0 - v;		  

		VD = (D*Va0*D.T()).inverse(ifail);	  
		if (ifail!=0) {
			if (m_Verbose) fprintf(stderr,"[ExKFitter] Warning : Cannot get inverse matrix (VD).\n");
			m_ErrorFlag = ExKF_MATRIXINVERSION;
		}else {
			
			VE = (E.T()*VD*E).inverse(ifail);
			if (ifail!=0) {
				if (m_Verbose) fprintf(stderr,"[ExKFitter] Warning : Cannot get inverse matrix (VE).\n");
				m_ErrorFlag = ExKF_MATRIXINVERSION;
			}else {
				
				Lambda0 = VD*(D*da0 + d);
				m_Chisq = (Lambda0.T() * (D*da0 + E*dv0 + d))(1,1);		      
				v = v - VE * E.T() * Lambda0;				      
			}
		}

		if (m_Chisq < chisq_keep_v && m_ErrorFlag==ExKF_NOERROR && m_Chisq>=0.) { 
			chisq_keep_v	= m_Chisq;
			v_keep_v	= v;
			VE_keep_v	= VE;
			VD_keep_v	= VD;
			Lambda0_keep_v  = Lambda0;
			E_keep_v	= E;
			D_keep_v	= D;
			fit_updated	= 1;
		}else {
			m_Chisq = chisq_keep_v;
			v	= v_keep_v;
			VE	= VE_keep_v;
			VD	= VD_keep_v;
			Lambda0 = Lambda0_keep_v;
			E	= E_keep_v;
			D	= D_keep_v;
			break;
		}

		vi = 0;   
		for(unsigned int it=0;it<m_ParticleList.size();it++) {
			ExKFitterParticle* par = m_ParticleList[it];
				
			if (par->m_IsPositionAvailable) {
			    if (!par->m_IsErrAvailable) {
				par->m_Momentum.setPx(v(vi+1));
				par->m_Momentum.setPy(v(vi+2));
				par->m_Momentum.setPz(v(vi+3));
				par->m_Position.setX(v(vi+4));
				par->m_Position.setY(v(vi+5));
				par->m_Position.setZ(v(vi+6));
					    
				vi += 6;
			    }
			}else {
			    if (!par->m_IsErrAvailable) {
				par->m_Momentum.setPx(v(vi+1));
				par->m_Momentum.setPy(v(vi+2));
				par->m_Momentum.setPz(v(vi+3));
					    
				vi += 3;
			    }
			}
		}
		for(unsigned int it=0;it<m_VertexList.size();it++) {
			ExKFitterVertex* ver = m_VertexList[it];
				
			if (!ver->m_IsErrAvailable) {
				ver->m_Vertex.setX(v(vi+1));
				ver->m_Vertex.setY(v(vi+2));
				ver->m_Vertex.setZ(v(vi+3));
				              
				vi += 3;
			}
		}
		Update();

	} // end of 'v' minizing loop

	if (m_ErrorFlag==ExKF_NOERROR) {
		Lambda  = Lambda0 - VD * E * VE * E.T() * Lambda0;
		a = a0 - Va0 * D.T() * Lambda;
	}	    

	if ((it_a == 0 || m_Chisq < chisq_keep) && m_ErrorFlag==ExKF_NOERROR && m_Chisq>=0.) {
		chisq_keep	  = m_Chisq;
		v_keep  	  = v;
		VE_keep 	  = VE;
		VD_keep 	  = VD;
		Lambda0_keep	  = Lambda0;
		E_keep  	  = E;
		D_keep  	  = D;
		a_keep  	  = a;    
		fit_updated	  = 1;    
	}else {
		m_Chisq   = chisq_keep;
		v	  = v_keep;
		VE	  = VE_keep;
		VD	  = VD_keep;
		Lambda0   = Lambda0_keep;
		E	  = E_keep;
		D	  = D_keep;
		a = a_keep;
		break;
	}

	ai = 0;     
	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
		          
		if (par->m_IsPositionAvailable) {
		    if (par->m_IsErrAvailable) {
		  	  par->m_Momentum.setPx(a(ai+1));
		  	  par->m_Momentum.setPy(a(ai+2));
		  	  par->m_Momentum.setPz(a(ai+3));
		  	  par->m_Position.setX(a(ai+4));
		  	  par->m_Position.setY(a(ai+5));
		  	  par->m_Position.setZ(a(ai+6));
		  	  	
		  	  ai += 6;
		    }
		}else {
		    if (par->m_IsErrAvailable) {
		  	  par->m_Momentum.setPx(a(ai+1));
		  	  par->m_Momentum.setPy(a(ai+2));
		  	  par->m_Momentum.setPz(a(ai+3));
		  	  	
		  	  ai += 3;
		    }
		}
	}
	for(unsigned int it=0;it<m_VertexList.size();it++) {
		ExKFitterVertex* ver = m_VertexList[it];
		          
		if (ver->m_IsErrAvailable) {
			ver->m_Vertex.setX(a(ai+1));
			ver->m_Vertex.setY(a(ai+2));
			ver->m_Vertex.setZ(a(ai+3));
					
			ai += 3;
		}
	}
	for(unsigned int it=0;it<m_MassList.size();it++) {
		ExKFitterMass* mas = m_MassList[it];
		          
		if (mas->m_IsErrAvailable) {
			mas->m_Mass = a(ai+1);
					
			ai += 1;
		}
	}

	Update();

	if (it_a==m_Max_Iteractions-1)
		if (m_Verbose) fprintf(stderr,"[ExKFitter] Warning : Reach the maximum iteractions.\n");

	} // end of 'a' minizing loop 

	if (fit_updated == 0) { // fit failed.
		if (m_Verbose) fprintf(stderr,"[ExKFitter] Error : Fit failed, return chi^2 = -1.\n");
		m_Chisq = -1.;
		if (m_ErrorFlag==0) {
			m_ErrorFlag=ExKF_BADINITIAL;				  
			if (m_Verbose) fprintf(stderr,"[ExKFitter] Error : Possible due to bad initial condition.\n");
		}
			    
		return m_ErrorFlag;
	}

	HepMatrix VDtu = VD - VD * E * VE * E.T() * VD;
	HepMatrix Va = Va0 - Va0 * D.T() * VDtu * D * Va0; // final error matrix

	ai = 0;       
	vi = 0;

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		if (par->m_IsPositionAvailable) {
		    if (par->m_IsErrAvailable) {
				
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
				par->m_ErrMomentumPosition(i,j)     = Va(ai+i,ai+j);
				par->m_ErrMomentumPosition(4+i,4+j) = Va(ai+3+i,ai+3+j);
			}
			}				
			for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
				par->m_ErrMomentumPosition(i,4+j)   = Va(ai+i,ai+3+j);
			}
			}							
							
			ai += 6;
		    }else {
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
				 par->m_ErrMomentumPosition(i,j)     = VE(vi+i,vi+j);
				 par->m_ErrMomentumPosition(4+i,4+j) = VE(vi+3+i,vi+3+j);
			}
			}				 
			for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
				 par->m_ErrMomentumPosition(i,4+j)   = VE(vi+i,vi+3+j);
			}
			}		 
			 		 
			par->m_IsErrAvailable = 1;											 
			vi += 6;
		    }
		}else {
		    if (par->m_IsErrAvailable) {
					
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
				par->m_ErrMomentumPosition(i,j)     = Va(ai+i,ai+j);
			}
			}											
								
			ai += 3;
		    }else {
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
			  	par->m_ErrMomentumPosition(i,j)     = VE(vi+i,vi+j);
			}
			}					
						
			par->m_IsErrAvailable = 1;											
			vi += 3;
		    }
		}
	}

	for(unsigned int it=0;it<m_VertexList.size();it++) { 
		ExKFitterVertex* ver = m_VertexList[it];
			    
		if (ver->m_IsErrAvailable) {
			          
			for(int i=1;i<=3;i++) 
			for(int j=1;j<=3;j++)
				ver->m_ErrVertex(i,j) = Va(ai+i,ai+j);
			ai += 3;
		}else {
			for(int i=1;i<=3;i++) 
			for(int j=1;j<=3;j++)
				ver->m_ErrVertex(i,j) = VE(vi+i,vi+j);	  
			        			  
			ver->m_IsErrAvailable = 1;						  
			vi += 3;
		}
	}

	for(unsigned int it=0;it<m_MassList.size();it++) { 
		ExKFitterMass* mas = m_MassList[it];
			    
		if (mas->m_IsErrAvailable) {		    
			mas->m_ErrMass = Va(ai+1,ai+1);
			ai += 1;
		}
	}

	return ExKF_NOERROR;
}

int ExKFitter::Minimize_wo_UnknownPar()
{			
	//Create related matrices
	int D_col;
	int a_row;

	D_col = 0;
	a_row = 0;

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		if (par->m_IsPositionAvailable) {
			if (par->m_IsErrAvailable) {D_col+=6;a_row+=6;}
		}else {
			if (par->m_IsErrAvailable) {D_col+=3;a_row+=3;}
		}
	}

	for(unsigned int it=0;it<m_VertexList.size();it++) {
		ExKFitterVertex* ver = m_VertexList[it];
			    
		if (ver->m_IsErrAvailable) {D_col+=3;a_row+=3;}
	}     

	for(unsigned int it=0;it<m_MassList.size();it++) {
		ExKFitterMass* mas = m_MassList[it];
			    
		if (mas->m_IsErrAvailable) {D_col+=1;a_row+=1;}
	}
	              
	double chisq_keep   = HUGE;
	m_Chisq = -1.0;

	int ifail = 0;
	int fit_updated = 0;	      

	HepVector a0(a_row),a(a_row);
	HepSymMatrix Va0(a_row,0);

	int ai = 0;

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		if (par->m_IsPositionAvailable) {
		    if (par->m_IsErrAvailable) {
			a0(ai+1) = par->m_Momentum.px();
			a0(ai+2) = par->m_Momentum.py();
			a0(ai+3) = par->m_Momentum.pz();
			a0(ai+4) = par->m_Position.x();
			a0(ai+5) = par->m_Position.y();
			a0(ai+6) = par->m_Position.z();
			    		
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
			  	Va0(ai+i,ai+j)     = par->m_ErrMomentumPosition(i,j);
			  	Va0(ai+3+i,ai+3+j) = par->m_ErrMomentumPosition(4+i,4+j);
			}
			}				
			for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
				Va0(ai+i,ai+3+j)   = par->m_ErrMomentumPosition(i,4+j);
			}
			}  
			    					
			ai += 6;
		    }
		}else {
		    if (par->m_IsErrAvailable) {
			a0(ai+1) = par->m_Momentum.px();
			a0(ai+2) = par->m_Momentum.py();
			a0(ai+3) = par->m_Momentum.pz();
			    		
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
				Va0(ai+i,ai+j)     = par->m_ErrMomentumPosition(i,j);
			}
			}				
			    					
			ai += 3;
		    }
		}
	}

	for(unsigned int it=0;it<m_VertexList.size();it++) { 
		ExKFitterVertex* ver = m_VertexList[it];
			    
		if (ver->m_IsErrAvailable) {
			a0(ai+1) = ver->m_Vertex.x();
			a0(ai+2) = ver->m_Vertex.y();
			a0(ai+3) = ver->m_Vertex.z();
			        	  
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
			  	Va0(ai+i,ai+j)     = ver->m_ErrVertex(i,j);
			}
			}			  
			        	  
			ai += 3;
		}
	}

	for(unsigned int it=0;it<m_MassList.size();it++) {
		ExKFitterMass* mas = m_MassList[it];
			    
		if (mas->m_IsErrAvailable) {
			a0(ai+1) 	= mas->m_Mass;
			Va0(ai+1,ai+1)  = mas->m_ErrMass;	  
			ai += 1;
		}
	}

	a = a0;

	HepMatrix D(m_N_Constraints,D_col,0);	      

	HepVector d(m_N_Constraints);

	HepVector da0(a_row); 

	HepMatrix VD(m_N_Constraints,m_N_Constraints);     

	HepMatrix Lambda(m_N_Constraints,1);

	HepMatrix VD_keep	= VD;
	HepMatrix Lambda_keep	= Lambda;
	HepMatrix D_keep	= D;
	HepVector a_keep	= a;    

	//=============================================================       
	// start of minizing loop

	int it_a;
	for(it_a=0;it_a<m_Max_Iteractions;it_a++) {	      // Global minizing loop
		
		m_ErrorFlag = ExKF_NOERROR; // reset error flag

		int ci = 0;	    
		for(unsigned int index=0;index<m_ConstrainList.size();index++) { 
			
			ExKFitterConstrain* con = m_ConstrainList[index];
			          
			int di = 0;
			          
			d.sub(ci+1,con->H());
			          
			for(unsigned int it=0;it<m_ParticleList.size();it++) { 
				ExKFitterParticle* par = m_ParticleList[it];						
					
				if (par->m_IsPositionAvailable) {
				    if (par->m_IsErrAvailable) {		      
					D.sub(ci+1,di+1,con->delH_delParticle(par).sub(1,con->N_Constraints(),1,6));
					di+=6;
				    }
				}else {
				    if (par->m_IsErrAvailable) {		      
					D.sub(ci+1,di+1,con->delH_delParticle(par).sub(1,con->N_Constraints(),1,3));
					di+=3;
				    }
				}
			}
			          
			for(unsigned int it=0;it<m_VertexList.size();it++) { 
				ExKFitterVertex* ver = m_VertexList[it];
						
				if (ver->m_IsErrAvailable) {
					D.sub(ci+1,di+1,con->delH_delVertex(ver));
					di+=3;
				}
			}	  
			          
			for(unsigned int it=0;it<m_MassList.size();it++) { 
				ExKFitterMass* mas = m_MassList[it];
					
				if (mas->m_IsErrAvailable) {
					D.sub(ci+1,di+1,con->delH_delMass(mas));
					di+=1;
				}
			}
			          
			if (con->m_ErrorFlag!=ExKF_NOERROR) m_ErrorFlag = con->m_ErrorFlag;
			          
			ci+= con->N_Constraints();
		}   

		da0 = a0 - a;	    

		VD = (D*Va0*D.T()).inverse(ifail);  
		if (ifail!=0) {
			if (m_Verbose) fprintf(stderr,"[ExKFitter] Warning : Cannot get inverse matrix (VD).\n");
			m_ErrorFlag = ExKF_MATRIXINVERSION;
		}

		if (m_ErrorFlag==ExKF_NOERROR) {
			Lambda  = VD*(D*da0 + d);
			m_Chisq = (Lambda.T() * (D*Va0*D.T()) * Lambda)(1,1);	  
			          
			a = a0 - Va0 * D.T() * Lambda;
		}	    

		if ((it_a == 0 || m_Chisq < chisq_keep) && m_ErrorFlag==ExKF_NOERROR && m_Chisq>=0.) {
			chisq_keep	  = m_Chisq;
			VD_keep 	  = VD;
			Lambda_keep	  = Lambda;
			D_keep  	  = D;
			a_keep  	  = a;    
			fit_updated	  = 1;    
		}else {
			m_Chisq   = chisq_keep;
			VD	  = VD_keep;
			Lambda    = Lambda_keep;
			D	  = D_keep;
			a = a_keep;
			break;
		}

		ai = 0;     
		for(unsigned int it=0;it<m_ParticleList.size();it++) {
			ExKFitterParticle* par = m_ParticleList[it];
			          
			if (par->m_IsPositionAvailable) {
			    if (par->m_IsErrAvailable) {
				par->m_Momentum.setPx(a(ai+1));
				par->m_Momentum.setPy(a(ai+2));
				par->m_Momentum.setPz(a(ai+3));
				par->m_Position.setX(a(ai+4));
				par->m_Position.setY(a(ai+5));
				par->m_Position.setZ(a(ai+6));
			
				ai += 6;
			    }
			}else {
			    if (par->m_IsErrAvailable) {
				par->m_Momentum.setPx(a(ai+1));
				par->m_Momentum.setPy(a(ai+2));
				par->m_Momentum.setPz(a(ai+3));
			
				ai += 3;
			    }
			}
		}
		for(unsigned int it=0;it<m_VertexList.size();it++) {
			ExKFitterVertex* ver = m_VertexList[it];
			          
			if (ver->m_IsErrAvailable) {
				ver->m_Vertex.setX(a(ai+1));
				ver->m_Vertex.setY(a(ai+2));
				ver->m_Vertex.setZ(a(ai+3));
						
				ai += 3;
			}
		}
		for(unsigned int it=0;it<m_MassList.size();it++) {
			ExKFitterMass* mas = m_MassList[it];
			          
			if (mas->m_IsErrAvailable) {
				mas->m_Mass = a(ai+1);
						
				ai += 1;
			}
		}

		Update();
			    
		if (it_a==m_Max_Iteractions-1)
			if (m_Verbose) fprintf(stderr,"[ExKFitter] Warning : Reach the maximum iteractions.\n");  

	} // end of 'a' minizing loop 

	if (fit_updated == 0) { // fit failed.
		if (m_Verbose) fprintf(stderr,"[ExKFitter] Error : Fit failed, return chi^2 = -1.\n");
		m_Chisq = -1.;
		if (m_ErrorFlag==0) m_ErrorFlag = ExKF_BADINITIAL;
			    
		return m_ErrorFlag;
	}

	HepMatrix Va = Va0 - Va0 * D.T() * VD * D * Va0; // final error matrix

	ai = 0;       

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		if (par->m_IsPositionAvailable) {
		    if (par->m_IsErrAvailable) {
		          
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
				par->m_ErrMomentumPosition(i,j)     = Va(ai+i,ai+j);
				par->m_ErrMomentumPosition(4+i,4+j) = Va(ai+3+i,ai+3+j);
			}
			}				
			for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
				par->m_ErrMomentumPosition(i,4+j)   = Va(ai+i,ai+3+j);
			}
			}							
							
			ai += 6;
		    }
		}else {
		    if (par->m_IsErrAvailable) {
		          
			for(int i=1;i<=3;i++) {
			for(int j=i;j<=3;j++) {
				par->m_ErrMomentumPosition(i,j)     = Va(ai+i,ai+j);
			}
			}											
							
			ai += 3;
		    }
		}
	}

	for(unsigned int it=0;it<m_VertexList.size();it++) { 
		ExKFitterVertex* ver = m_VertexList[it];
			    
		if (ver->m_IsErrAvailable) {
			          
			for(int i=1;i<=3;i++) 
			for(int j=1;j<=3;j++)
			    	ver->m_ErrVertex(i,j) = Va(ai+i,ai+j);
			ai += 3;
		}
	}

	for(unsigned int it=0;it<m_MassList.size();it++) { 
		ExKFitterMass* mas = m_MassList[it];
			    
		if (mas->m_IsErrAvailable) {		    
			mas->m_ErrMass = Va(ai+1,ai+1);
			ai += 1;
		}
	}

	return ExKF_NOERROR;
}

double ExKFitter::EffectiveXi(ExKFitterConstrain* Master)
{
	if (m_Chisq<0.) return -1.;
	if (Master->m_ConstrainType != ExKF_VERTEXCONSTRAIN) return -1.;

	double Xi_sum = 0.;

	for(unsigned int it=0;it<m_ParticleList.size();it++) { 
		ExKFitterParticle* par = m_ParticleList[it];
			    
		if (par->m_Charge==0.) continue; // Only take charged tracks
			    
		int flag = 0;
		for(unsigned int i=0;i<m_ConstrainList.size();i++) { 
			ExKFitterConstrain* con = m_ConstrainList[i];		  
			if ( con->m_ConstrainType == ExKF_VERTEXCONSTRAIN &&
			     !con->m_FollowVertex) {
				for(unsigned int j=0;j<con->m_ParticleList.size();j++) 
				  	if (con->m_ParticleList[j] == par) flag = 1;
			}
		}
		if (!flag) continue; // Only take the tracks related to one of the vertex constraints.
			    
		flag = 0;
		for(unsigned int j=0;j<Master->m_ParticleList.size();j++) 
			if (Master->m_ParticleList[j]->CheckRelation(par)) flag = 1;
				    
		if (!flag) continue; // Only if this track related to the master constraint.	    
			    
		int Norm = 1;
		for(unsigned int i=0;i<m_ConstrainList.size();i++) { 
			ExKFitterConstrain* con = m_ConstrainList[i];	  
			        		  
			if ( con->m_ConstrainType == ExKF_VERTEXCONSTRAIN &&
			     !con->m_FollowVertex) {
			  	for(unsigned int j=0;j<con->m_ParticleList.size();j++)
				if (con->m_ParticleList[j]->CheckRelation(par)) {				      
					Norm *= con->m_ParticleList.size();									    
					break;
				}
			}
		}
			    
		Xi_sum += par->TrackXi()/(double)Norm;
	}

	return Xi_sum/2.;
}

int ExKFitter::N_EffectiveDOF(ExKFitterConstrain* Master)
{
	if (m_Chisq<0.) return 0;
	if (Master->m_ConstrainType != ExKF_VERTEXCONSTRAIN) return 0;

	return Master->m_ParticleList.size()*2;
}

int ExKFitter::N_VertexingTracks()
{
	int ntrks = 0;

	for(unsigned int it=0;it<m_ParticleList.size();it++) { 
		ExKFitterParticle* par = m_ParticleList[it];
			    
		int flag = 0;
		for(unsigned int i=0;i<m_ConstrainList.size();i++) { 
			ExKFitterConstrain* con = m_ConstrainList[i];		  
			if (con->m_ConstrainType == ExKF_VERTEXCONSTRAIN) {
				for(unsigned int j=0;j<con->m_ParticleList.size();j++) 
					if (con->m_ParticleList[j] == par) flag = 1;
			}
		}
			    
		if (flag && par->m_Charge!=0.) ntrks++;
	}

	return ntrks;
}

double ExKFitter::Xi()
{
	double Xi_sum = 0.;
	int ntrks = 0;

	if (m_Chisq<0.) return -1.;
	if (m_Verbose) fprintf(stderr,"[ExKFitter] Warning : Usual Xi() is called. Please check your definition of xi.\n");

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		int flag = 0;
		for(unsigned int i=0;i<m_ConstrainList.size();i++) { 
			ExKFitterConstrain* con = m_ConstrainList[i];		  
			if (con->m_ConstrainType == ExKF_VERTEXCONSTRAIN) {
				for(unsigned int j=0;j<con->m_ParticleList.size();j++)
					if (con->m_ParticleList[j] == par) flag = 1;
			}
		}
			    
		if (flag && par->m_Charge!=0.) {
			Xi_sum += par->TrackXi();
			ntrks++;
		}
	}

	// Use the normalization of 2*ntrks instead of N(d.o.f.)
	return Xi_sum/(double)ntrks/2.;
}

double ExKFitter::Zeta()
{
	double Zeta_sum = 0.;
	int ntrks = 0;

	if (m_Chisq<0.) return -1.;

	for(unsigned int it=0;it<m_ParticleList.size();it++) {
		ExKFitterParticle* par = m_ParticleList[it];
			    
		int flag = 0;
		for(unsigned int i=0;i<m_ConstrainList.size();i++) {
			ExKFitterConstrain* con = m_ConstrainList[i];		  
			if (con->m_ConstrainType == ExKF_VERTEXCONSTRAIN) {
				for(unsigned int j=0;j<con->m_ParticleList.size();j++)
					if (con->m_ParticleList[j] == par) flag = 1;
			}
		}
			    
		if (flag && par->m_Charge!=0.) {
			Zeta_sum += par->TrackZeta();
			ntrks++;
		}
	}

	//3 degree of freedom per Ks = 1.5 per track
	return Zeta_sum/((double)ntrks/2.*3.);
}


void ExKFitter::Update()
{
	for(unsigned int it=0;it<m_ConstrainList.size();it++) 
		m_ConstrainList[it]->Update();
	              
	for(unsigned int it=0;it<m_ExternalParticleList.size();it++) 
		m_ExternalParticleList[it]->Update();	    
}

void ExKFitter::Bfield(double d)
{
	for(unsigned int it=0;it<m_ConstrainList.size();it++) 
		m_ConstrainList[it]->Bfield(d);
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

