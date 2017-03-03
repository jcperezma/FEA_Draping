#include "SemiDiscreteTriangle.h"
double const EPSILON = std::numeric_limits<double>::epsilon()*1e8;

// returns the area of the triangle as half the norm of the cross product bewtween its consittuent vectors
double SDTriangle::computeArea(std::vector<vector3D>const &coords){
	vector3D pointA = coords[nodes[0]]; 
	vector3D pointB = coords[nodes[1]];	
	vector3D pointC = coords[nodes[2]];	

	return 0.5*norm(cross(pointB-pointA,pointC-pointA) ); // 0.5 * |(B-A)x(C-A)|
};
	
void SDTriangle::computeMassMatrix(std::vector<vector3D>const &coords, std::vector<double> & M, double rho){
	double mass = rho*computeArea(coords);	

	for (int ii = 0; ii < nodes.size(); ii++)
	{
		for (int dim = 0; dim < 3; dim++)
		{
			int M_index = nodes[ii]*3+dim;
			M[M_index] += 1.0/3.0*mass;
		}
	}

}

double computeC_gamma(double gamma, double ct1, double ct2, double cte3){
	double gamma3 = gamma*gamma*gamma;
	double gamma5 = gamma3 *gamma*gamma;
	return ct1*gamma+ ct2*gamma3+ cte3*gamma5;
}

vector3D SDTriangle::getk_1(std::vector<vector3D> & coords)const{
	// need to find k1, k_1, k2 and k_2
	// also g_1 and g_2

	vector3D g_1 = coords[nodes[1]]-coords[nodes[0]];
	vector3D g_2 = coords[nodes[2]]-coords[nodes[0]];

	vector3D k_1 = g_1 - a * g_2;

	return k_1;
}

vector3D SDTriangle::getk_2(std::vector<vector3D> & coords)const{
	// need to find k1, k_1, k2 and k_2
	// also g_1 and g_2

	vector3D g_1 = coords[nodes[1]]-coords[nodes[0]];
	vector3D g_2 = coords[nodes[2]]-coords[nodes[0]];

	vector3D k_2 = -b*g_1 + g_2;

	return k_2;
}

void SDTriangle::computeBvalues(std::vector<vector3D> & coords){
	// need to find k1, k_1, k2 and k_2
	// also g_1 and g_2

	vector3D g_1 = coords[nodes[1]]-coords[nodes[0]];
	vector3D g_2 = coords[nodes[2]]-coords[nodes[0]];
	
	vector3D k_1 = g_1 - a * g_2;
	vector3D k_2 = -b*g_1 + g_2;
	vector3D k_3 = cross(k_1,k_2);


	vector3D D = cross(k_1,k_3);

	vector3D K1 = cross(k_2,k_3)*(1/dot(k_1,cross(k_2,k_3)));
	vector3D K2 = cross(k_3,k_1)*(1/dot(k_2,cross(k_3,k_1)));
	vector3D K3 = cross(k_1,k_2)*(1/dot(k_3,cross(k_1,k_2)));

	double n1 =10; // I dont know what this number might be, but is related to the area of the triangle, and unit cell density
	double n2 =10;
	double  nCER =10;

	double mag_k_1 = norm(k_1);
	double mag_k_1sq = mag_k_1 *mag_k_1; 
	double mag_k_2 = norm(k_2);
	double mag_k_2sq = mag_k_2 *mag_k_2; 
	double mag_k1 = norm(K1);
	double mag_k2 = norm(K2);

	// B11, will help update the strain in the 1 direction
	// First node
	B11[0] = (a-1)*k_1.x/mag_k_1sq; B11[1] = (a-1)*k_1.y/mag_k_1sq; B11[2] = (a-1)*k_1.z/mag_k_1sq; 
	// second node
	B11[3] = k_1.x/mag_k_1sq;		  B11[4] = k_1.y/mag_k_1sq;		B11[5] = k_1.z/mag_k_1sq; 
	// third node
	B11[6] = -a*k_1.x/mag_k_1sq;	  B11[7] = -a*k_1.y/mag_k_1sq;	B11[8] = -a*k_1.z/mag_k_1sq; 

	// B22
	// First node
	B22[0] = (b-1)*k_2.x/mag_k_2sq; B22[1] = (b-1)*k_2.y/mag_k_2sq; B22[2] = (b-1)*k_2.z/mag_k_2sq; 
	// second node
	B22[3] = -b*k_2.x/mag_k_2sq;		  B22[4] = -b*k_2.y/mag_k_2sq;		B22[5] = -b*k_2.z/mag_k_2sq; 
	// third node
	B22[6] = k_2.x/mag_k_2sq;	  B22[7] = k_2.y/mag_k_2sq;	B22[8] = k_2.z/mag_k_2sq; 

	std::vector<double> B3(9), B4(9);

	
		// B3
	// First node
	B3[0] = (b-1)*k_1.x; B3[1] = (b-1)*k_1.y; B3[2] = (b-1)*k_1.z; 
	// second node
	B3[3] = -b*k_1.x;		  B3[4] = -b*k_1.y;		B3[5] = -b*k_1.z; 
	// third node
	B3[6] = k_1.x;	  B3[7] = k_1.y;	B3[8] = k_1.z; 


	// B4 for each node
		// First node
	B4[0] = (a-1)*k_2.x; B4[1] = (a-1)*k_2.y; B4[2] = (a-1)*k_2.z; 
	// second node
	B4[3] = k_2.x;		  B4[4] = k_2.y;		B4[5] = k_2.z; 
	// third node
	B4[6] = -a*k_2.x;	  B4[7] = -a*k_2.y;	B4[8] = -a*k_2.z; 

		double part1 = dot(K1/mag_k_1,K2/mag_k2);
	double part2 = dot(K2/mag_k_2,K1/mag_k1);
	// B_gamma
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		for (int j = 0; j < 3; j++)
		{
			int localIndex = ii*3+j;
			Bgamma[localIndex] += B11[localIndex]* part1 + B22[localIndex] *part2 + B3[localIndex] * mag_k2/mag_k_1 + B4[localIndex] * mag_k1/mag_k_2 ;
		}
	}
	

}

void SDTriangle::computeBbs(std::vector<vector3D> & coords, std::vector<SDTriangle> & tris, vector3D & h1, vector3D & h2 ){

	// I'm going to do the operations one by one instead of matricially, this will reduce the number of computations tremendously;
	// eq 3.3.16-19 on hamila's thesis
	vector3D R11 = 1/d[0] *e_z;
	vector3D R12 = -cosBeta[0]/d[1] *e_z;
	vector3D R13 = -cosBeta[1]/d[2] *e_z;

	vector3D R21 = -cosBeta[2]/d[0] *e_z;
	vector3D R22 = 1/d[1] *e_z;
	vector3D R23 = -cosBeta[0]/d[2] *e_z;

	vector3D R31 = -cosBeta[1]/d[0] *e_z;
	vector3D R32 = -cosBeta[0]/d[1] *e_z;
	vector3D R33 = 1/d[2] *e_z;

	double d4=0, d5=0, d6 =0;

	// side 0
	vector3D R41(0,0,0);
	vector3D R42(0,0,0);
	vector3D R44(0,0,0);
	int localSideID = 0;
	// see hamila's thesis for more info about the gammas and betas and so forth, page 57
	if (neighborNodes[localSideID]!=-1) // theres is a neighboring triangle
	{
		SDTriangle const & neighborTri = tris[neighborTriID[localSideID]];
		int index = neighborSideID[localSideID]; // index of the "lower" vertex of the neighboring triangle in the neighbouring side, lower in CCW winding sense
		double gamma2 = neighborTri.cosBeta[index];
		double gamma1 = neighborTri.cosBeta[(index+1)%sidesPerElement];
			   d4 = neighborTri.d[(index+2)%sidesPerElement];
		double q2 = neighborTri.d[(index)];
		double q1 = neighborTri.d[(index+1)%sidesPerElement];
		vector3D e_z4=  neighborTri.e_z;
		R41 = -gamma2/q1 * e_z4;
		R42 = -gamma1/q2 * e_z4;;
		R44 = 1.0/d4 *e_z4;;
	}
	

	// side 1
	vector3D R52(0,0,0);
	vector3D R53(0,0,0);
	vector3D R55(0,0,0);
	localSideID = 1;
	if (neighborNodes[localSideID]!=-1) // theres is a neighboring triangle
	{
		SDTriangle const & neighborTri = tris[neighborTriID[localSideID]];
		int index = neighborSideID[localSideID]; // index of the "lower" vertex of the neighboring triangle in the neighbouring side, lower in CCW winding sense
		double psi3 = neighborTri.cosBeta[index];
		double psi2 = neighborTri.cosBeta[(index+1)%sidesPerElement];
		         d5 = neighborTri.d[(index+2)%sidesPerElement];
		double r3 = neighborTri.d[(index)];
		double r2 = neighborTri.d[(index+1)%sidesPerElement];
		vector3D e_z5=  neighborTri.e_z;
	
		R52 = -psi3/r2 * e_z5;
		R52 = -psi2/r3 * e_z5;
		R55 = 1.0/d5 *e_z5;
	}

	//side 2
	vector3D R61(0,0,0);
	vector3D R63(0,0,0);
	vector3D R66(0,0,0);
	localSideID = 2;
	if (neighborNodes[localSideID]!=-1) // theres is a neighboring triangle
	{
		SDTriangle const & neighborTri = tris[neighborTriID[localSideID]];
		int index = neighborSideID[localSideID]; // index of the "lower" vertex of the neighboring triangle in the neighbouring side, lower in CCW winding sense
		double phi1 = neighborTri.cosBeta[index];
		double phi3 = neighborTri.cosBeta[(index+1)%sidesPerElement];
		         d6 = neighborTri.d[(index+2)%sidesPerElement];
		double p1 = neighborTri.d[(index)];
		double p3 = neighborTri.d[(index+1)%sidesPerElement];
		vector3D e_z6=  neighborTri.e_z;
	
		R61 = -phi3/p1 * e_z6;
		R63 = -phi1/p3 * e_z6;
		R66 = 1.0/d6 *e_z6;
	}

	double P11 = d[0]/(d[0]+d5);
	double P15 = d[0]/(d[0]+d5);
	double P22 = d[1]/(d[1]+d6);
	double P26 = d[1]/(d[1]+d6);
	double P33 = d[2]/(d[2]+d4);
	double P34 = d[2]/(d[2]+d4);



	std::vector<vector3D> PR1(6); // Row 1 of matrix PR 
	PR1[0] = P11*R11;	
	PR1[1] = P11*R12; + P15 *R52;
	PR1[2] = P11*R13; + P15 *R53;	
	PR1[3] = 0;
	PR1[4] = P15*R55;

	std::vector<vector3D> PR2(6); // Row 2 of matrix PR 
	PR2[0] = P22*R21 + P26 *R61;	
	PR2[1] = P22*R22;
	PR2[2] = P22*R23 + P26 *R63;	
	PR2[5] = P26*R66;

	std::vector<vector3D> PR3(6); // Row 3 of matrix PR 
	PR3[0] = P33*R31 + P34 *R41;	
	PR3[1] = P33*R32 + P34 *R42;
	PR3[2] = P33*R33 ;	
	PR3[3] = P34*R44;


	// find normals
	vector3D n1 = cross(coords[nodes[1]]-coords[nodes[0]],e_z);
	n1 = n1/norm(n1);

	vector3D n2 = cross(coords[nodes[2]]-coords[nodes[1]],e_z);
	n2 = n2/norm(n2);

	vector3D n3 = cross(coords[nodes[0]]-coords[nodes[2]],e_z);
	n3 = n3/norm(n3);

	// for direction 1
	std::vector<double> D(3);

	D[0] = -2.0/d[0] * dot(h1,n1)*dot(h1,n1);
	D[1] = -2.0/d[1] * dot(h1,n2)*dot(h1,n2);
	D[2] = -2.0/d[2] * dot(h1,n3)*dot(h1,n3);

	for (int i = 0; i < 6; i++)
	{
		Bb1[i*3]   = D[0] * PR1[i].x + D[1] * PR2[i].x + D[2] * PR3[i].x;
		Bb1[i*3+1] = D[0] * PR1[i].y + D[1] * PR2[i].y + D[2] * PR3[i].y;
		Bb1[i*3+2] = D[0] * PR1[i].z + D[1] * PR2[i].z + D[2] * PR3[i].z;
	}

	// for direction 2
	D[0] = -2.0/d[0] * dot(h2,n1)*dot(h2,n1);
	D[1] = -2.0/d[1] * dot(h2,n2)*dot(h2,n2);
	D[2] = -2.0/d[2] * dot(h2,n3)*dot(h2,n3);
	for (int i = 0; i < 6; i++)
	{
		Bb2[i*3]   = D[0] * PR1[i].x + D[1] * PR2[i].x + D[2] * PR3[i].x;
		Bb2[i*3+1] = D[0] * PR1[i].y + D[1] * PR2[i].y + D[2] * PR3[i].y;
		Bb2[i*3+2] = D[0] * PR1[i].z + D[1] * PR2[i].z + D[2] * PR3[i].z;
	}

}

void SDTriangle::computeInternalForce(double C1, double C2 , double ct1, double ct2, double ct3, double n1, double n2, std::vector<double> & F_int, std::vector<vector3D> & coords, std::vector<double> & displacements, std::vector<SDTriangle> & tris , double B_k, std::vector<double> & F_11, std::vector<double> & F_22, std::vector<double> & F_s ){
	// need to find k1, k_1, k2 and k_2
	// also g_1 and g_2

	vector3D g_1 = coords[nodes[1]]-coords[nodes[0]];
	vector3D g_2 = coords[nodes[2]]-coords[nodes[0]];

	vector3D k_1 = g_1 - a * g_2;
	vector3D k_2 = -b*g_1 + g_2;
	vector3D k_3 = cross(k_1,k_2);

	vector3D K1 = cross(k_2,k_3)*(1/dot(k_1,cross(k_2,k_3)));
	vector3D K2 = cross(k_3,k_1)*(1/dot(k_2,cross(k_3,k_1)));
	vector3D K3 = cross(k_1,k_2)*(1/dot(k_3,cross(k_1,k_2)));
	
	nCER =	0.0001;

	double mag_k_1 = norm(k_1);
	double mag_k_2 = norm(k_2);
	double mag_k1 = norm(K1);
	double mag_k2 = norm(K2);
	double mag_g_1 = norm(g_1);
	double mag_g_2 = norm(g_2);


	// compute the B matrices for tension, shear and bending

	// Tensile B
	// B11, will help update the strain in the 1 direction
	// First node
	B11[0] = (a-1)*k_1.x; B11[1] = (a-1)*k_1.y; B11[2] = (a-1)*k_1.z; 
	// second node
	B11[3] = k_1.x;		  B11[4] = k_1.y;		B11[5] = k_1.z; 
	// third node
	B11[6] = -a*k_1.x;	  B11[7] = -a*k_1.y;	B11[8] = -a*k_1.z; 

	// B22
	// First node
	B22[0] = (b-1)*k_2.x; B22[1] = (b-1)*k_2.y; B22[2] = (b-1)*k_2.z; 
	// second node
	B22[3] = -b*k_2.x;		  B22[4] = -b*k_2.y;		B22[5] = -b*k_2.z; 
	// third node
	B22[6] = k_2.x;	  B22[7] = k_2.y;	B22[8] = k_2.z; 

	
	// Shear B
	std::vector<double> B3(9), B4(9);
	// B3
	// First node
	B3[0] = (b-1)*k_1.x;	B3[1] = (b-1)*k_1.y;	B3[2] = (b-1)*k_1.z; 
	// second node
	B3[3] = -b*k_1.x;		B3[4] = -b*k_1.y;		B3[5] = -b*k_1.z; 
	// third node
	B3[6] = k_1.x;			B3[7] = k_1.y;			B3[8] = k_1.z; 


	// B4 for each node
		// First node
	B4[0] = (a-1)*k_2.x;	B4[1] = (a-1)*k_2.y;	B4[2] = (a-1)*k_2.z; 
	// second node
	B4[3] = k_2.x;			B4[4] = k_2.y;			B4[5] = k_2.z; 
	// third node
	B4[6] = -a*k_2.x;		B4[7] = -a*k_2.y;		B4[8] = -a*k_2.z; 

	double part1 = dot(K1/mag_k_1,K2/mag_k2);
	double part2 = dot(K2/mag_k_2,K1/mag_k1);
	// B_gamma
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		for (int j = 0; j < 3; j++)
		{
			int localIndex = ii*3+j;
			Bgamma[localIndex] = B11[localIndex]* part1 + B22[localIndex] *part2 + B3[localIndex] * mag_k2/mag_k_1 + B4[localIndex] * mag_k1/mag_k_2 ;
		}
	}
	
	// bending B

	vector3D h1 =  k_1 /mag_k_1;
	vector3D h2 =  k_2 /mag_k_2;

	computeBbs(coords,tris,h1,h2);

	// find curvature in both directions


	// first nodes that are part of the triangle
	for (int i = 0; i < nodes.size(); i++)
	{
		int globalID = nodes[i]*3;
		for (int j = 0; j < 3; j++)
		{
			Xi1 += Bb1[i*3+j] *  displacements[globalID+j];
			Xi2 += Bb2[i*3+j] *  displacements[globalID+j];
		}
	}

	// fnow nodes that are neighboring the triangle
	for (int i = 0; i < nodes.size(); i++)
	{
		int globalID = neighborNodes[i]*3;
		if(globalID>=0){
			for (int j = 0; j < 3; j++)
			{
				Xi1 += Bb1[9+i*3+j] *  displacements[globalID+j];
				Xi2 += Bb2[9+i*3+j] *  displacements[globalID+j];
			}
		}
	}


	// update shear deformation and deformation in the direction of the threads
	for (int i = 0; i < nodes.size(); i++)
	{
		int globalID = nodes[i]*3;
		for (int j = 0; j < 3; j++)
		{
			epsilon11+= B11[i*3+j]*displacements[globalID+j];
			epsilon22+= B22[i*3+j]*displacements[globalID+j];
			gamma+= Bgamma[i*3+j]*displacements[globalID+j];
		}
	}

	// if gamma is very high, wrinkle will form, add a force on all triangles normal to the triangle

	if(abs(gamma)>0.4){
		vector3D N = cross(g_1,g_2);
		N = N/ (norm(N)*1000000);
		for (int ii = 0; ii < nodes.size(); ii++)//nodes.size()
		{
			int index = nodes[ii]*3;
			F_int[index] += N.x;
			F_int[index+1] += N.y;
			F_int[index+2] += N.z;
		}

	}

	// Add forces to global force vector

	// tensile forces direction 1
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		int index = nodes[ii]*3;
		int localIndex = ii*3;
		for (int j = 0; j < 3; j++)
		{
			F_int[index+j] += n2 / mag_k_1 * C1*epsilon11 * B11[localIndex+j];
			F_11[index+j] += n2 / mag_k_1 * C1*epsilon11 * B11[localIndex+j];
		}
	}

	// tensile forces direction 2
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		int index = nodes[ii]*3;
		int localIndex = ii*3;
		for (int j = 0; j < 3; j++)
		{
			F_int[index+j] += n1 / mag_k_2 * C2*epsilon22 * B22[localIndex+j];
			F_22[index+j] += n1 / mag_k_2 * C2*epsilon22 * B22[localIndex+j];
		}
	}

	// shear forces
	double C_gamma = computeC_gamma(gamma, ct1,ct2,ct3);
	//C_gamma =0.001; // Nm
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		int index = nodes[ii]*3;
		int localIndex = ii*3;
		for (int j = 0; j < 3; j++)
		{
			F_s[index+j] += nCER * C_gamma *Bgamma[localIndex+j]; 
			F_int[index+j] += nCER * C_gamma *Bgamma[localIndex+j];   
		}
	}

	// Bending forces direction 1
	double M11 = B_k *Xi1;
	// for nodes within the element
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		int index = nodes[ii]*3;
		int localIndex = ii*3;
		for (int j = 0; j < 3; j++)
		{
			//F_int[index+j] += n1 /2.0 * mag_g_1 * M11 *Bb1[localIndex+j];   
		}
	}
	// for nodes of neighboring elements
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		int index = neighborNodes[ii]*3;
		if(index>=0){
			int localIndex = 9+ii*3;
			for (int j = 0; j < 3; j++)
			{
				//F_int[index+j] += n1 /2.0 * mag_g_1 * M11 *Bb1[localIndex+j];   
			}
		}
	}

	//For the direction 2
	double M22 = B_k *Xi2;
	// for nodes within the element
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		int index = nodes[ii]*3;
		int localIndex = ii*3;
		for (int j = 0; j < 3; j++)
		{
			//F_int[index+j] += n2 /2.0 * mag_g_2 * M22 *Bb2[localIndex+j];   
		}
	}
	// for nodes of neighboring elements
	for (int ii = 0; ii < nodes.size(); ii++)
	{
		int index = neighborNodes[ii]*3;
		if(index>=0){// if there is a neighbor
			int localIndex = 9+ii*3;
			for (int j = 0; j < 3; j++)
			{
				//F_int[index+j] += n2 /2.0 * mag_g_2 * M22 *Bb2[localIndex+j];   
			}
		}
	}



	double mag_k_1sq = mag_k_1 *mag_k_1; 

	for (int i = 0; i < B11.size(); i++) B11[i] /= mag_k_1sq;

	double mag_k_2sq = mag_k_2 *mag_k_2; 
	for (int i = 0; i < B22.size(); i++) B22[i] /= mag_k_2sq;


};

	// check if the multiplication of the z components of the vectors AB and AC is positive,  if it is
	// that means that the vector h "passes first" through vertex A;
bool testPointAndVector(vector3D const & A, vector3D const &B, vector3D const  &C , vector2D const & h){

	vector3D p1 = B-A;
	vector3D p2 = C-A;
	double z1 = p1.x*h.y-p1.y*h.x;
	double z2 = h.x*p2.y-h.y*p2.x;
	return (z1*z2>=-1e-7);

}

#define perp(u,v)  ((u).x * (v).y - (u).y * (v).x)  // perp product  (2D)

// finds intersaction of ray A h on segment BC
double findIntersection(vector3D const & A, vector3D const & B, vector3D  const & C, vector2D const & h ){

	vector3D    u(h.x,h.y,0.0); //S1.P1 - S1.P0;
    vector3D    v = C-B;//S2.P1 - S2.P0;
    vector3D    w = A - B;
    float     D = perp(u,v);

    
    // the segments are skew and may intersect in a point
    // get the intersect parameter for S1
    double     tI = perp(u,w) / D;

	//
	tI = tI>1.0?1:tI;
	tI = tI<0.0?0:tI;
	return tI;

    

}

// finds the values a and b that hold the information for the material coordinates
void SDTriangle::findaAndB(vector2D const &h1, vector2D const &h2, std::vector<vector3D>& coords ){

	// kinda complicated logic, but the idea is that the direction of the vectors h1, h2 qont affect the equations later, 
	// they are based on the fact that the origin of both the material coordinates and element coordinates is the same node A,
	// if one of the h_i enters through the origin node A, and if the vector is not colineal with the edges of the triangle, 
	// A cannot be the origin and the nodes indexing is rotates counterclockwise, if the vector is colineal then it can still be the origiin, 
	// and in that case either a or b are 0.
	
	bool repeat = true;
	while(repeat){
		vector3D A = coords[nodes[0]];
		vector3D B = coords[nodes[1]];
		vector3D C = coords[nodes[2]];
		repeat  = false;
		bool aFound=false;
		bool bFound=false;
		if (testPointAndVector(A,B,C,h1) ){
			// check if the vector A->h1 is colineal with A->B
			double result = norm(cross(vector3D(h1.x,h1.y,0),B-A));
			if(norm(cross(vector3D(h1.x,h1.y,0),B-A))<1e-6 ){
				a=0;
				aFound=true;
			} // check if the vector A->h1 is colineal with A->C
			else if(norm(cross(vector3D(h1.x,h1.y,0),C-A))<1e-6){
				int temp = nodes[0];
				nodes[0] = nodes[1];
				nodes[1] = nodes[2];
				nodes[2] = temp;
				repeat = true;
				continue;
			}else{
			// If true need to rotate nodes
			int temp = nodes[0];
			nodes[0] = nodes[2];
			nodes[2] = nodes[1];
			nodes[1] = temp;
			repeat = true;
			continue;
			}
		} 
		if (!aFound&&testPointAndVector(B,C,A,h1) ){a = findIntersection(B,A,C,h1);aFound=true;}
		if (!bFound&&testPointAndVector(C,A,B,h1) ){
			int temp = nodes[0];
				nodes[0] = nodes[1];
				nodes[1] = nodes[2];
				nodes[2] = temp;
				repeat = true;
				continue;
		}
	
		
		if (testPointAndVector(A,B,C,h2) ){
			// check if the vector A->h2 is colineal with A->B
			double result = norm(cross(vector3D(h2.x,h2.y,0),B-A));
			if(norm(cross(vector3D(h2.x,h2.y,0),B-A))<1e-6 ){
				// h2 needs to be parallel to k_2, rotate nodes CCW
				int temp = nodes[0];
				nodes[0] = nodes[1];
				nodes[1] = nodes[2];
				nodes[2] = temp;
				repeat = true;
				continue;
			} // check if the vector A->h2 is colineal with A->C
			else if(norm(cross(vector3D(h2.x,h2.y,0),C-A))<1e-6){
				b=0;
				bFound=true;
			}else{
			// If true need to rotate nodes
			int temp = nodes[0];
				nodes[0] = nodes[1];
				nodes[1] = nodes[2];
				nodes[2] = temp;
				repeat = true;
				continue;
			}
		} 
		if (!aFound&&testPointAndVector(B,C,A,h2) ){
			int temp = nodes[0];
				nodes[0] = nodes[1];
				nodes[1] = nodes[2];
				nodes[2] = temp;
				repeat = true;
				continue;
		}
		if (!bFound&&testPointAndVector(C,A,B,h2) ){b = findIntersection(C,A,B,h2);}
	
	}
	
	// the input files have a precision of about 1e6
	a = a>0.999?1:a;
	a = a<0.001?0:a;

	b = b>0.999?1:b;
	b = b<0.001?0:b;

};

void SDTriangle::updateStrain(std::vector<double> &displacements){
	
	for (int i = 0; i < nodes.size(); i++)
	{
		int globalID = nodes[i]*3;
		for (int j = 0; j < 3; j++)
		{
			epsilon11+= B11[i*3+j]*displacements[globalID+j];
			epsilon22+= B22[i*3+j]*displacements[globalID+j];
			gamma+= Bgamma[i*3+j]*displacements[globalID+j];
		}
	}

};

template <typename T>
inline T min(T a, T b ) {
	if (a<b) return a;
	else return b;
} ;

template <typename T>
inline T max(T a, T b ) {
	if (a>b) return a;
	else return b;
} ;

side SDTriangle::getSide(int sideID,int triID){
	// a triangle has three sides obviously, 
	// the side ID are given by the node order
	// side 0 node 0-1, side 1 node 1-2 etc..
	// the nodes of sides are in a sorted order, that is
	// side.node1 < side.node2

	int node1 = nodes[sideID];
	int node2 = nodes[(sideID+1)%3];

	side aSide;
	aSide.node1 = min(node1,node2);
	aSide.node2 = max(node1,node2);
	aSide.triID = triID;
	aSide.sideID= sideID;
	aSide.oposingNode = nodes[(sideID+2)%3];
	return aSide;
}

void SDTriangle::computeBendingVars(std::vector<vector3D> & coords){
	
	vector3D pointA = coords[nodes[0]];
	vector3D pointB = coords[nodes[1]];
	vector3D pointC = coords[nodes[2]];

	vector3D side0 = pointB -pointA;
	vector3D side1 = pointC-pointB;
	vector3D side2 = pointA-pointC;

	// vector normal to the triangle
	e_z = cross(side0,-1*side2);
	e_z = e_z/norm(e_z); 

	// triangle heights from each corner
	d[0] = norm( cross( -1*side0,side2  )  ) / norm(side1);
	d[1] = norm( cross( side0,-1*side1  )  ) / norm(side2);
	d[2] = norm( cross( side1,-1*side2  )  ) / norm(side0);

	// cosine of angles of each of the triangle corners
	cosBeta[0] =  dot(side0,-1*side2)/ (norm(side0)*norm(side2)) ;
	cosBeta[1] =  dot(side1,-1*side0)/ (norm(side1)*norm(side0)) ;
	cosBeta[2] =  dot(side2,-1*side1)/ (norm(side2)*norm(side1)) ;

}