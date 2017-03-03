#include "DrapingMesh.h"
#include <stdio.h>

#if defined _MSC_VER && _MSC_VER >= 1200

#define _sprintf sprintf_s

#else

#define _sprintf std::sprintf

#endif 

class Fabric
{
public:
	Fabric(){fabricID= numOfFabrics; numOfFabrics++;};
	Fabric(std::string fileName){initilizeFromFile(fileName);fabricID= numOfFabrics; numOfFabrics++;};
	~Fabric(){};
	void initilizeFromFile(std::string fileName);
	void advanceStep(double dt);
	void findContacts();	// pass a kinetic object to the function
	void print();	// outputs vtu file
	void print(int step);	// outputs vtu file
	static void printPVDfile(int numSteps);
	void initializeTest2();
private:
	static int numOfFabrics;
	int fabricID;
	int stepID;
	semiDiscreteMesh mesh;
	double rho;		// material density
	double T11;		// material tensile stiffness in direction 1
	double T22;		// material tensile stiffness in direction 2
	double ct1;		// cti, are constants to compute the shear stifness
	double ct2;		// 
	double ct3;		// 
	double B_k;		// Bending stifness 

	vector3D k1;	// initial warp direction
	vector3D k2;	// initial weft direction

	double n1;		// threads per m
	double n2;		
};

