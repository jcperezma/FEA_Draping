#pragma once
#include <vector>
#include "vectorOps.h"

#define sidesPerElement 3 

struct side{
	int node1;
	int node2;
	int triID;	// triangle are not aware of their IDs
	int sideID; // side is not aware of which side it is
	int oposingNode;

	bool operator<(side const & sideB){
		 if(node1==sideB.node1){
			return node2<sideB.node2;
		} else{
			return node1<sideB.node1;
		 }
	}

	bool operator>(side const & sideB){
		if(node1==sideB.node1){
			return node2>sideB.node2;
		} else{
			return node1>sideB.node1;
		 }
	}

	bool isEqual(side const & sideB){
		return (node1==sideB.node1) && (node2==sideB.node2);
	}
};

class SDTriangle
{
public:
	SDTriangle(){B11.resize(9);B22.resize(9);Bgamma.resize(9);epsilon11=0,epsilon22=0;gamma=0;Xi1=0;Xi2=0;
	neighborNodes.resize(3);std::fill(neighborNodes.begin(),neighborNodes.end(),-1);neighborSideID.resize(3); neighborTriID.resize(3);
	cosBeta.resize(3), d.resize(3); Bb1.resize(18);Bb2.resize(18);}
	~SDTriangle(){}

	void setNodes(std::vector<int> & nodes){this->nodes =nodes;};	// done
	void computeMassMatrix(std::vector<vector3D>const &coords, std::vector<double> & M, double rho); // done
	void computeInternalForce(double C11, double C22 , std::vector<double> & F_int, std::vector<vector3D> & coords  );
	void computeInternalForce(double C11, double C22 ,double ct1, double ct2, double ct3, double n1, double n2, std::vector<double> & F_int, std::vector<vector3D> & coords , std::vector<double> & displacements, std::vector<SDTriangle> & tris, double B_k, std::vector<double> & F_11, std::vector<double> & F_22, std::vector<double> & F_s );
	void findaAndB(vector2D const &h1, vector2D const &h2, std::vector<vector3D>& coords );
	void computeBvalues(std::vector<vector3D> & coords);
	void updateStrain(std::vector<double> &displacements);
	void computeBendingVars(std::vector<vector3D> & coords);
	void computeBbs(std::vector<vector3D> & coords, std::vector<SDTriangle> & tris, vector3D & h1, vector3D & h2 );

	void setNeighborNode(int sideID, int nodeID){neighborNodes[sideID] = nodeID;}
	void setNeighboringSide(int sideIDA, int sideIDB){neighborSideID[sideIDA] = sideIDB;}
	void setNeighboringTri(int sideID, int triID){neighborTriID[sideID] = triID;}

	std::vector<int> getNodes(){return nodes;};
	side getSide(int sideID,int triID);
	double getEpsilon11()const{return epsilon11;};
	double getEpsilon22()const{return epsilon22;};
	double getGamma()const{return gamma;};
	double getXi1()const{return Xi1 ;};
	double getXi2()const{return Xi2;};
	vector3D getk_1(std::vector<vector3D> & coords)const;
	vector3D getk_2(std::vector<vector3D> & coords)const;

private:
	double computeArea(std::vector<vector3D>const &coords);			// done
	double a;			// to help reconstruct k_1
	double b;			// to help reconstruct k_2
	double gamma;		// shear deformation
	double epsilon11;	// elongation in direction 1
	double epsilon22;	// elongation in direction 2
	double Xi1 ;		// curvature in direction 1
	double Xi2 ;		// curvature in direction 2
	double nCER;		// number of unitary cells
	std::vector<double> B11;
	std::vector<double> B22;
	std::vector<double> Bgamma;
	std::vector<double> Bb1;
	std::vector<double> Bb2;
	std::vector<int> nodes;
	// for bending
	std::vector<double> cosBeta;
	std::vector<double> d;			// distance between each vertex to opposing side
	std::vector<int> neighborNodes; // nodes of the neighbouring elements, if its -1 it does not have neighbors
	std::vector<int> neighborSideID; // local sideID within neighbouring triangle
	std::vector<int> neighborTriID; // ID of neighbouring triangle
	vector3D e_z;	// vector normal to triangle
};


