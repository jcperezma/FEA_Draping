#pragma once
#include "SemiDiscreteTriangle.h"
#include <vector>
#include <sstream> 
#include <fstream>
#include <string>
#include <iostream>
#include <memory>
#include <algorithm>
#include <omp.h>
// mesh that has semi discrete triangles,
// it will be used to solve the explicit ecuation of motion for 
// anisotropic fabrics. 

struct BC{
	double value;
	double ID;
};

class semiDiscreteMesh
{
public:
	semiDiscreteMesh(){};
	~semiDiscreteMesh(){};

	// set variables
	void setCoords(std::vector<vector3D> Coords){
		coords=Coords; int nDOFs = Coords.size()*3;
		F_int.resize(nDOFs); F_ext.resize(nDOFs);
		u.resize(nDOFs); v.resize(nDOFs); a.resize(nDOFs);
		M_D.resize(nDOFs); F_11.resize(nDOFs); F_22.resize(nDOFs); F_s.resize(nDOFs);
	}
	void setElements(std::vector<SDTriangle> Tris){tris=Tris;};
	void setBCs(std::vector<BC> newBCS){BCs= newBCS;}

	// at initializiation compute mass matrix and a and b for each triangle, and find neighbouring triangles
	void computeMassMatrix(double rho);
	void findAandB(vector3D k1, vector3D k2 );
	void findNeighbours();

	// 
	void computeBendingVariables();
	void computeInternalForces(double C1, double C2, double ct1, double ct2, double ct3,double n1, double n2, double B_k);
	void updatePositions(double dt);

	void initializeTest1();
	void initializeTest2(double rho, vector3D k_1, vector3D k_2);
	void printVTUfile(std::string fileName);

private:
	std::vector<SDTriangle> tris;
	std::vector<vector3D> coords;	// vector with the nodal coordinates "x"
	std::vector<double> F_int;		// vector of internal forces for each node {Fx1, Fy1, Fz1, Fx2, Fy2,Fz2...}
	std::vector<double> F_11;		// vector of internal forces caused by extension in direction 11
	std::vector<double> F_22;		// vector of internal forces caused by extension in direction 22
	std::vector<double> F_s;		// vector of internal forces caused by shear
	std::vector<double> F_ext;		// vector of external forces for each node, same format as above
	std::vector<double> u;			// vector of nodal displacements
	std::vector<double> v;			// vector of nodal velocities
	std::vector<double> a;			// vector of nodal accelerations
	std::vector<double> M_D;		// Lumped mass matrix
	std::vector<BC> BCs;			// vector with values for the velocity BCs

};

