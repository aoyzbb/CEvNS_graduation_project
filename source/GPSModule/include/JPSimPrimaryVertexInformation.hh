#ifndef JPSimPrimaryVertexInformation_h
#define JPSimPrimaryVertexInformation_h 1

#include "globals.hh"
#include "G4VUserPrimaryVertexInformation.hh"

class JPSimPrimaryVertexInformation :
	public G4VUserPrimaryVertexInformation
{
public:
	JPSimPrimaryVertexInformation() : 
		RadZ(0), RadA(0) {}
	JPSimPrimaryVertexInformation(int z, int a) :
		RadZ(z), RadA(a) {}

	~JPSimPrimaryVertexInformation() {}

	void Print() const 
	{ G4cout<<"Z: "<<RadZ<<", A: "<<RadA<<G4endl;}

	int GetZ() {return RadZ;}
	int GetA() {return RadA;}


private:
	int RadZ;
	int RadA;








};






#endif
