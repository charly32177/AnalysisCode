/*
	File : ExKFitter.h (The main header file)
	Description : 	Experimental Extended Kinematic Fitter
		A set of c++ class developed for general vertex/mass 
		constrain fitting with the supports of multi-constraints 
		and simultaneous minimization.
			
	Author : Kai-Feng Chen, National Taiwan University	      
	Version : 1.31
	Usage (An example, vertex+mass+IP constraints for B+->D0 pi+):
	    
	Mdst_charged Mdst_D_K;
	Mdst_charged Mdst_D_Pi;
	Mdst_charged Mdst_Pi;
	    
	ExKFitterParticle D_K(Mdst_D_K,3);
	ExKFitterParticle D_Pi(Mdst_D_Pi,2);
	ExKFitterParticle Pi(Mdst_Pi,2);
		
	// Remember to set this Initial_D_Vertex
	// for the case of multi-vertex simultaneous fit
	ExKFitterVertex D_Vertex(Initial_D_Vertex); 
		
	// Assign IP and IPerr to B_Vertex = IP constrained fit 
	ExKFitterVertex B_Vertex(IP,IPerr);	    
		
	ExKFitterParticle D;
	D.LinkParticle(&D_K);
	D.LinkParticle(&D_Pi);
	D.LinkVertex(&D_Vertex);			      
	              
	ExKFitterMass D_Mass(1.8645);				      
	              
	ExKFitterConstrain con1;
	con1.SetVertexConstrain();
	con1.LinkParticle(&D_K);
	con1.LinkParticle(&D_Pi);
	con1.LinkVertex(&D_Vertex);
	              
	ExKFitterConstrain con2;
	con2.SetVertexConstrain();
	con2.LinkParticle(&D);   
	con2.LinkParticle(&Pi); 	 
	con2.LinkVertex(&B_Vertex);
	              
	ExKFitterConstrain con3;
	con3.SetMassConstrain();
	con3.LinkParticle(&D_K);
	con3.LinkParticle(&D_Pi);
	con3.LinkVertex(&D_Vertex);
	con3.LinkMass(&D_Mass);
	              
	ExKFitter Core; 
	Core.LinkConstrain(&con1);
	Core.LinkConstrain(&con2);
	Core.LinkConstrain(&con3);
	int ret = Core.Minimize();
	        			      
	return of Minimize() = 0   -> Success.
			     < 0   -> Failed. (the chi^2 will be -1. also.)   
	        			   
	The results will direct store back into the linked 'ExKFitterParticle' class.
	B_Vertex.Vertex()     -> The results of B vertex
	B_Vertex.ErrVertex()  -> The error matrix of the vertex
	              
	For the B 4-momentum after D-mass constraining, it can be calculated in this way:
	              
	ExKFitterParticle B;
	B.LinkParticle(&D);
	B.LinkParticle(&Pi);
	B.LinkVertex(&B_Vertex);
	B.Update();
	        		      
	HepLorentzVector B_P4 = B.Momentum();

	Date :	Apr/03,2004 - First released (developing) version
		Apr/06,2004 - Bug fixed in chi^2 return of Minimize_wo_UnknownPar()
		Apr/07,2004 - Add assignment operator in all classes.
		Apr/11,2004 - Bug fixed in delHvc_delVertex(), wrong assigned matrix size.
		Apr/15,2004 - Change the 'xi' definition from chi_z^2/d.o.f -> chi_z^2/(n*ntrks)
			    - Some of the names of class members were changed.
		Apr/20,2004 - Add 'zeta' for Ks-B vertexing as well
		Apr/23,2004 - Add int ExKFitter::N_VertexingTracks() and only count real used 
			      vertexing tracks for xi/zeta calculation.
		Apr/25,2004 - Try to use adaptive no. of iteractions.
		May/06,2004 - Add the 'followVertex' (not-to-fit) option in the ExKFitterConstrain 
			      class, and change the normalization of Xi.
		Jun/03,2004 - Error messages updated.
		Jun/24,2004 - Better treatment for photons (remove the position
			      from the fit.) Auto replace the unsupportted struture.
		Jul/01,2004 - C++ corrections by Ushiroda-san (Thanks a lot!)
		Jul/02,2004 - Remove Ushiroda-san's mark
*/

#include <belle.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "helix/Helix.h"

#include "particle/Particle.h"  
#include "particle/gammac.h"

#include "panther/panther.h"
#include "panther/panther_manager.h"
#include MDST_H

#define ExKF_DEF_MAX_ITERACTIONS 	20
#define ExKF_DEF_B_FIELD		1.5

#define ExKF_VERTEXCONSTRAIN		1
#define ExKF_MASSCONSTRAIN		2

#define ExKF_NOERROR 		 0
#define ExKF_MATHERROR		-1
#define ExKF_MATRIXINVERSION	-2
#define ExKF_BADINITIAL		-3
#define ExKF_BADCONSTRAINS	-4
#define ExKF_UNDEFINED		-9

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class ExKFitterParticle;
class ExKFitterVertex;
class ExKFitterMass;
class ExKFitterConstrain; 
class ExKFitter;

class ExKFitterParticle {
public:
	HepLorentzVector m_Momentum;
	HepPoint3D	 m_Position;
	HepSymMatrix	 m_ErrMomentumPosition;
	int	       m_IsErrAvailable;
	int	       m_IsVirual;    
	int	       m_IsPositionAvailable;
	double m_Charge;
	double m_Mass;
	double m_Bfield;
	double m_R;

	HepPoint3D    m_Pivot;
	HepVector     m_Helix;
	HepSymMatrix  m_ErrHelix;
	int	      m_IsHelixAvailable;

	std::vector<ExKFitterParticle*> m_ParticleList; 
	ExKFitterVertex* m_Vertex;  

	int CheckRelation(ExKFitterParticle* Particle);
	int CheckRelation(ExKFitterVertex* Vertex);

	HepMatrix RelationMatrix(ExKFitterParticle* Particle);
	HepMatrix RelationMatrix(ExKFitterVertex* Vertex);

	void LinkParticle(ExKFitterParticle* Particle);
	void LinkParticleRecursively(ExKFitterParticle* Particle);
	void LinkVertex(ExKFitterVertex* Vertex);

	void Update();
	void Bfield(double d);

	// Constructors	              
	ExKFitterParticle();  
	ExKFitterParticle(HepLorentzVector Momentum,
	        	  HepPoint3D	   Position,
	        	  HepSymMatrix     ErrMomentumPosition,
	        	  double	   Charge);
	ExKFitterParticle(Mdst_charged& charged,
	        	  int		pid);   
	ExKFitterParticle(Mdst_vee2&  	Vee2,
	        	  double      	Charge);  		
	ExKFitterParticle(Mdst_gamma& 	gamma,
	        	  HepPoint3D	Vertex,
	        	  HepSymMatrix  ErrVertex);
	ExKFitterParticle(Mdst_gamma&	gamma);			      
	        		
	ExKFitterParticle & operator = (const ExKFitterParticle &p);		     
	        					
	double alpha() {return - 0.00299792458 * m_Bfield * m_Charge;}

	void GetHelixParameters( Mdst_charged&	charged,
	        		 int	      	pid,
	        		 HepPoint3D&	pivot,
	        		 HepVector&   	helix,
	        		 HepSymMatrix&  error); 
	void GetHelixParameters( Mdst_vee2&    	Vee2,
	        		 double        	charge,
	        		 HepPoint3D&   	pivot,
	        		 HepVector&    	helix,
	        		 HepSymMatrix& 	error);        
	        		       
	double TrackChisq( HepLorentzVector  Momentum,
	        	   HepPoint3D	   Position);  
	double TrackXi(); 
	double TrackZeta();	      
	              
	// Direct member access references    

	std::vector<ExKFitterParticle*>& ParticleList() {return m_ParticleList;}      
	ExKFitterVertex* Vertex() {return m_Vertex;}

	HepLorentzVector& Momentum()	       {return m_Momentum;}
	HepPoint3D&	Position()	       {return m_Position;}
	HepSymMatrix&	ErrMomentumPosition()  {return m_ErrMomentumPosition;}
	int&		IsErrAvailable()       {return m_IsErrAvailable;} 
	int&		IsVirual()	       {return m_IsVirual;} 
	int&		IsPositionAvailable()  {return m_IsPositionAvailable;} 
	double& 	Charge()	       {return m_Charge;}
	double& 	Mass()  	       {return m_Mass;}
	double& 	Bfield()	       {return m_Bfield;}
	double& 	R()		       {return m_R;}

	void Momentum(HepLorentzVector m)      {m_Momentum = m;}
	void Position(HepPoint3D p)	       {m_Position = p;}
	void ErrMomentumPosition(HepSymMatrix e) {m_ErrMomentumPosition = e;}
	void IsErrAvailable(int n)	       {m_IsErrAvailable = n;} 
	void IsVirual(int n)		       {m_IsVirual = n;} 
	void IsPositionAvailable(int n)        {m_IsPositionAvailable = n;} 
	void Charge(double d)		       {m_Charge = d;} 
	void Mass(double d)		       {m_Mass = d;}  
	void R(double d)		       {m_R = d;}		      
};

class ExKFitterVertex {
public:
	HepPoint3D   m_Vertex;
	HepSymMatrix m_ErrVertex;
	int	     m_IsErrAvailable;

	ExKFitterVertex();
	ExKFitterVertex(HepPoint3D   Vertex);
	ExKFitterVertex(HepPoint3D   Vertex,
	        	HepSymMatrix ErrVertex);
	        	      
	ExKFitterVertex & operator = (const ExKFitterVertex &p);	      
	        	      
	HepPoint3D&   Vertex()         {return m_Vertex;}
	HepSymMatrix& ErrVertex()      {return m_ErrVertex;}
	int&	      IsErrAvailable() {return m_IsErrAvailable;} 

	void Vertex(HepPoint3D v)      {m_Vertex = v;}
	void ErrVertex(HepSymMatrix e) {m_ErrVertex = e;}
	void IsErrAvailable(int n)     {m_IsErrAvailable = n;}        
};

class ExKFitterMass {
public:
	double m_Mass;
	double m_ErrMass;
	int    m_IsErrAvailable;

	ExKFitterMass();
	ExKFitterMass(double Mass);
	ExKFitterMass(double Mass, double ErrMass);

	ExKFitterMass & operator = (const ExKFitterMass &p);

	double& Mass()  	   {return m_Mass;}
	double& ErrMass()	   {return m_ErrMass;}
	int&	IsErrAvailable()   {return m_IsErrAvailable;} 

	void Mass(double m)	   {m_Mass = m;}
	void ErrMass(double e)     {m_ErrMass = e;}
	void IsErrAvailable(int n) {m_IsErrAvailable = n;} 
};

class ExKFitterConstrain {
public:
	std::vector<ExKFitterParticle*> m_ParticleList; 
	ExKFitterVertex* m_Vertex;  
	ExKFitterMass* m_Mass;
	int m_ErrorFlag;
	int m_ConstrainType;
	int m_FollowVertex;

	ExKFitterConstrain() {
		m_Vertex    = NULL;
		m_Mass      = NULL;
		m_ErrorFlag = ExKF_NOERROR;
		m_ConstrainType = ExKF_VERTEXCONSTRAIN; // It will perform a vertex constrained fit by default.
		m_FollowVertex  = 0;
	}

	ExKFitterConstrain & operator = (const ExKFitterConstrain &p);
	              
	void LinkParticle(ExKFitterParticle* Particle);
	void LinkParticleRecursively(ExKFitterParticle* Particle);
	void LinkVertex(ExKFitterVertex* Vertex) {m_Vertex = Vertex;}
	void LinkMass(ExKFitterMass* Mass) {m_Mass = Mass;}
	              
	void SetVertexConstrain()  {m_ConstrainType = ExKF_VERTEXCONSTRAIN;}
	void SetMassConstrain()    {m_ConstrainType = ExKF_MASSCONSTRAIN;}
	void ConstrainType(int c)  {m_ConstrainType = c;}
	void FollowVertex(int f=1) {m_FollowVertex  = f;}

	HepVector H() {        // = constrain vector H(a,v) = d
		if (m_ConstrainType == ExKF_VERTEXCONSTRAIN)
			return Hvc();
		else if (m_ConstrainType == ExKF_MASSCONSTRAIN)
			return Hmc();
		else {
			fprintf(stderr,"[ExKFitterConstrain] Error : No such ConstrainType defined.\n");
			return HepVector();
		}
	}
	HepVector Hvc();
	HepVector Hmc();

	HepMatrix delH_delParticle(ExKFitterParticle* Particle) { 
		if (m_ConstrainType == ExKF_VERTEXCONSTRAIN)
			return delHvc_delParticle(Particle);
		else if (m_ConstrainType == ExKF_MASSCONSTRAIN)
			return delHmc_delParticle(Particle);
		else {
			fprintf(stderr,"[ExKFitterConstrain] Error : No such ConstrainType defined.\n");
			return HepMatrix();
		}
	}
	HepMatrix delH_delVertex(ExKFitterVertex* Vertex) {   
		if (m_ConstrainType == ExKF_VERTEXCONSTRAIN)
			return delHvc_delVertex(Vertex);
		else if (m_ConstrainType == ExKF_MASSCONSTRAIN)
			return delHmc_delVertex(Vertex);  
		else {
			fprintf(stderr,"[ExKFitterConstrain] Error : No such ConstrainType defined.\n");
			return HepMatrix();
		}	    
	}
	HepMatrix delH_delMass(ExKFitterMass* Mass) { 
		if (m_ConstrainType == ExKF_VERTEXCONSTRAIN)
			return delHvc_delMass(Mass);
		else if (m_ConstrainType == ExKF_MASSCONSTRAIN)
			return delHmc_delMass(Mass);  
		else { 
			fprintf(stderr,"[ExKFitterConstrain] Error : No such ConstrainType defined.\n");
			return HepMatrix();
		}
	}

	HepMatrix delHvc_delParticle(ExKFitterParticle* Particle); 
	HepMatrix delHvc_delVertex(ExKFitterVertex* Vertex);   
	HepMatrix delHvc_delMass(ExKFitterMass* Mass);         
	HepMatrix delHmc_delParticle(ExKFitterParticle* Particle);
	HepMatrix delHmc_delVertex(ExKFitterVertex* Vertex); 
	HepMatrix delHmc_delMass(ExKFitterMass* Mass);         
	              
	int N_Constraints() {
		if (m_ConstrainType == ExKF_VERTEXCONSTRAIN)
			return 2*m_ParticleList.size();
		else if (m_ConstrainType == ExKF_MASSCONSTRAIN)
			return 1;
		else{ 
			fprintf(stderr,"[ExKFitterConstrain] Error : No such ConstrainType defined.\n");
			return -1;
		}
	}

	int CheckAvailavility() {
		if (m_ConstrainType == ExKF_VERTEXCONSTRAIN &&
		    m_ParticleList.size()>=2 && 
		    m_Vertex != NULL) return 1;
				
		if (m_ConstrainType == ExKF_VERTEXCONSTRAIN &&
		    m_ParticleList.size()>=1 && 
		    m_Vertex != NULL &&
		    (m_Vertex->m_IsErrAvailable || m_FollowVertex)) return 1;  
				
		if (m_ConstrainType == ExKF_MASSCONSTRAIN &&
		    m_ParticleList.size()>=2 && 
		    m_Vertex != NULL) return 1;       
				
		return 0;
	}

	void Update();
	void Bfield(double d);
	        		      
	std::vector<ExKFitterParticle*>& ParticleList() {return m_ParticleList;}      
	ExKFitterVertex* Vertex() {return m_Vertex;}
	ExKFitterMass* Mass() {return m_Mass;}

	int ConstrainType() {return m_ConstrainType;}
	int ErrorFlag() {return m_ErrorFlag;}
};

class ExKFitter {	
public:
	std::vector<ExKFitterParticle*> m_ParticleList; 
	std::vector<ExKFitterVertex*> m_VertexList;
	std::vector<ExKFitterMass*> m_MassList;  
	std::vector<ExKFitterConstrain*> m_ConstrainList;
	std::vector<ExKFitterParticle*> m_ExternalParticleList; 

	int m_Max_Iteractions;
	int m_Verbose;
	        	      
	int m_ErrorFlag;
	double m_Chisq;

	int m_N_Constraints;
	int m_N_DegreeOfFreedom;

	ExKFitter() {
		m_ErrorFlag = ExKF_NOERROR; 
		m_Chisq = -1.;
		m_N_Constraints = m_N_DegreeOfFreedom = 0;
		m_Max_Iteractions = -1; //Adaptive += (ExKF_DEF_MAX_ITERACTIONS/2) * m_ConstrainList.size();
		m_Verbose = 1;
	}

	ExKFitter & operator = (const ExKFitter &p);

	void LinkParticleRecursively(ExKFitterParticle* Particle);
	void LinkExternalParticle(ExKFitterParticle* Particle);
	void LinkVertex(ExKFitterVertex* Vertex);
	void LinkMass(ExKFitterMass* Mass);
	void LinkConstrain(ExKFitterConstrain* Constrain);

	int PreMinimize();
	int N_Constraints() {return m_N_Constraints;}
	int N_DegreeOfFreedom() {return m_N_DegreeOfFreedom;}
	int ErrorFlag() {return m_ErrorFlag;}
	              
	int Minimize();
	int Minimize_w_UnknownPar();
	int Minimize_wo_UnknownPar();

	void Update();
	void Bfield(double d);

	double Chisq() {return m_Chisq;}
	              
	double EffectiveXi(ExKFitterConstrain* Master);  
	int N_EffectiveDOF(ExKFitterConstrain* Master);       
	              
	double Xi();
	double Zeta();
	int N_VertexingTracks();

	std::vector<ExKFitterParticle*>& ParticleList() {return m_ParticleList;}      
	std::vector<ExKFitterVertex*>& VertexList() {return m_VertexList;}
	std::vector<ExKFitterMass*>& MassList() {return m_MassList;}
	std::vector<ExKFitterConstrain*>& ConstrainList() {return m_ConstrainList;}
	std::vector<ExKFitterParticle*>& ExternalParticleList() {return m_ExternalParticleList;}

	int& MaxIteractions()	   {return m_Max_Iteractions;}  	    
	int& Verbose()  	   {return m_Verbose;}
	        	      
	void MaxIteractions(int n) {m_Max_Iteractions = n;}
	void Verbose(int n)	   {m_Verbose = n;}

};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif


