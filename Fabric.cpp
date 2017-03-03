#include "Fabric.h"

int Fabric::numOfFabrics =0;

std::vector<vector3D> readNodesCoords(std::ifstream & tempFile){
	std::vector<vector3D>coords;
	std::string line;
	// read number of nodes

	int NumNodes;
	std::getline(tempFile,line);

	{ // have not found a better way to flush the stream
	std::stringstream ss;
	ss<<line;
	ss>>NumNodes;
	}
	

	//read node coords
	for (int i = 0; i < NumNodes; i++)
	{
		vector3D p;
		getline (tempFile,line);
		{ // have not found a better way to flush the stream
		std::stringstream ss;
		ss<<line;
		ss>>p.x;
		ss>>p.y;
		ss>>p.z;
		}
		coords.push_back(p);
	}
	return coords;
}

std::vector<SDTriangle> readElements(std::ifstream & tempFile){
	std::vector<SDTriangle> elements;
	std::string line;
	std::getline (tempFile,line);

	int NumElements;
	{ // have not found a better way to flush the stream
		std::stringstream ss;
		ss<<line;
		ss>>NumElements;
	}
	elements.resize(NumElements);
	// read element topology
	for (int i = 0; i < NumElements; i++)
	{
		std::vector<int> nodeID(3);
		int node;
		int elementType;
		getline (tempFile,line);
			{ // have not found a better way to flush the stream
			std::stringstream ss;
			ss<<line;
				for (int j = 0; j < 3; j++)
				{
					ss>>node;
					nodeID[j]=node;
				}
			}
		SDTriangle tri;
		tri.setNodes(nodeID);
		elements[i]=tri;	
	}

	return elements;
}

std::vector<BC> readBC(std::ifstream & tempFile ){
	std::vector<double>BC_value;
	std::vector<int>BC_ID;

	std::string line;
	// read displacement constrains
	std::getline (tempFile,line);
		int NumUConstrains;
	{ // have not found a better way to flush the stream
		std::stringstream ss;
		ss<<line;
		ss>>NumUConstrains;
		std::stringstream().swap(ss);// flush stream
	}
	//read indices of the constrains
	for (int i = 0; i < NumUConstrains; i++)
	{
		int U_ID;
		getline (tempFile,line);
		{ // have not found a better way to flush the stream
			std::stringstream ss;
			ss<<line;
			ss>>U_ID;
		}
		BC_ID.push_back(U_ID);
	}

	// read displacement constrain values
	for (int i = 0; i < NumUConstrains; i++)
	{
		double U_value;
		std::getline (tempFile,line);
		{ // have not found a better way to flush the stream
		std::stringstream ss;
		ss<<line;
		ss>>U_value;
		}
		BC_value.push_back(U_value);
	}

	std::vector<BC> BCs;

	for (int i = 0; i < BC_ID.size(); i++)
	{
		BC aBC;
		aBC.ID = BC_ID[i];
		aBC.value = BC_value[i];
		BCs.push_back(aBC);
	}
	return BCs;
}


void Fabric::initilizeFromFile(std::string fileName){
	// 1. read data form file
	
	std::string line;

	//open file
	std::ifstream tempFile (fileName);
	bool didntFail = tempFile.is_open();


	// read simulation parameters
	getline (tempFile,line);
	{ // have not found a better way to flush the stream
	std::stringstream ss;
	ss<<line;
	// Tensile properties of the fabric
	ss>>T11; ss>>T22; 
	// Shear properties of the fabric
	ss>>ct1; ss>>ct2; ss>>ct3;
	// bending stiffness
	ss>>B_k;
	}

	


	// read k1 and k2
	getline (tempFile,line);
	{ // have not found a better way to flush the stream
	std::stringstream ss;
	ss<<line;
	// Tensile properties of the fabric
	ss>>k1.x; ss>>k1.y; ss>>k1.z; 
	// Shear properties of the fabric
	ss>>k2.x; ss>>k2.y; ss>>k2.z; 
	}

	// read k1 and k2
	getline (tempFile,line);
	{ // have not found a better way to flush the stream
	std::stringstream ss;
	ss<<line;
	// surface density of the fabric
	ss>>rho;
	// threads per m
	ss>>n1; ss>>n2; 
	}

	// read mesh topology
	std::vector<vector3D> coords = readNodesCoords(tempFile); 
	mesh.setCoords(coords);

	// read element connectivity
	std::vector<SDTriangle> elements = readElements(tempFile);
	mesh.setElements(elements);

	// read BC for the velocity
	std::vector<BC> BCs =  readBC(tempFile);
	mesh.setBCs(BCs);

	tempFile.close();
	
	// Compute M_D matrix
	mesh.computeMassMatrix(rho);

	// find A and B for every triangle
	mesh.findAandB(k1,k2);

	mesh.findNeighbours();

}

void Fabric::advanceStep(double dt){
	mesh.computeBendingVariables();
	mesh.computeInternalForces(T11,T22,ct1,ct2,ct3,n1,n2,B_k);
	mesh.updatePositions(dt);

}

void Fabric::print(){
	char filename[32];
	_sprintf(filename, "fabric_%2d.vtu", fabricID);
	mesh.printVTUfile(filename);
};	// outputs vtu file

void Fabric::print(int step){
	char filename[32];
	_sprintf(filename, "fabric_%2dstep%05d.vtu", fabricID,step);
	mesh.printVTUfile(filename);
};	// outputs vtu file


void prindPVDDataset (std::ofstream &out, int step, std::string filename ){
	out << "<DataSet timestep=\"" << step << "\"  group = \"\" part=\"0\" file= \""<< filename << "\"/>" <<std::endl;
}

void  printPVDfileFooter(std::ofstream & out){
	// footer
	out << "  </Collection>" << std::endl;
	out << "</VTKFile>";
	out.close();
}
void printPVDfileHeading(std::ofstream &out){
	// header
	out << "<VTKFile type=\"" << "Collection" << "\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
	out << "  <Collection>" << std::endl;
}

void Fabric::printPVDfile(int numSteps){
	// print a PVD file for all the Fabric objects in the simulation

	for (int i = 0; i < numOfFabrics; i++)
	{
		char pvdfilename[32];
		_sprintf(pvdfilename, "fabric_%02d.pvd", i);
		std::ofstream ss;
		ss.open(pvdfilename);
		// print header
		printPVDfileHeading(ss);
		for (int step = 0; step < numSteps; step++)
		{
			char filename[32];
			_sprintf(filename, "fabric_%2dstep%05d.vtu", i,step);
			prindPVDDataset(ss,step,filename);
		}
		printPVDfileFooter(ss);
	}

}

void Fabric::initializeTest2(){
	rho = 0.6;
	k1.x=1; k1.y=0;k1.z=0;
	k2.x=0; k2.y=1;k2.z=0;
	T11 = 2600; T22 = 2600;
	ct1=0.37044; ct2= 0.84105; ct3 = 1.03113;
	n1 = 100; n2 = 100;
	mesh.initializeTest2(rho, k1,k2);


}