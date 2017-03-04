#pragma once
#include "DrapingMesh.h"



void semiDiscreteMesh::findAandB(vector3D k1, vector3D k2 ){
	for (auto & tri : tris)
	{
		tri.findaAndB(k1,k2,coords);
	}
}
void semiDiscreteMesh::computeMassMatrix(double rho){
	// the mass matrix is lumped
	// it is computed at each element and added to the global lumped matrix
	for (auto & tri : tris)
	{
		tri.computeMassMatrix(coords,M_D, rho);
	}

}
void semiDiscreteMesh::computeInternalForces(double C11, double C22, double ct1, double ct2, double ct3, double n1, double n2, double B_k){
	for (int i = 0; i < F_int.size(); i++){ F_int[i]=0; F_11[i]=0; F_22[i]=0; F_s[i]=0;}

	// all variables are shared
#pragma omp parallel for 
	for (int i =0 ; i< tris.size();i++ )//for (auto & triangle : tris)
	{
		tris[i].computeInternalForce(C11,C22,ct1,ct2,ct3,n1, n2,F_int,coords,u,tris, B_k,F_11,F_22,F_s);
		//triangle.computeInternalForce(C11,C22,ct1,ct2,ct3,n1, n2,F_int,coords,u,tris, B_k,F_11,F_22,F_s);
	}

}

void semiDiscreteMesh::updatePositions(double dt){
	double alpha =0.8;
	// TODO apply boundary conditions

	// update acceleration
	for (int i = 0; i < u.size(); i++)
	{
		// new acceleration
		//a[i] = 1.0/M_D[i]*(F_ext[i]-F_int[i]- alpha*M_D[i]*v[i]);
		a[i] = 1.0/M_D[i]*(F_ext[i]-F_int[i]);
	}
	// new velocity
	for (int i = 0; i < v.size(); i++)
	{
		v[i] = v[i] + dt * a[i];
	}
	// TODO set BC
		
	for (auto const & bc : BCs)
	{
		v[bc.ID] = bc.value;
	}

	std::vector<double> u_old=u;
	std::vector<double> du(u.size());
	// find displacements
	for (int i = 0; i < u.size(); i++)
	{
		u[i] = dt * v[i];
	}

	// compute du
	for (int i = 0; i < u.size(); i++)
	{
		du[i] = u[i]-u_old[i] ;
	}

	for (int i = 0; i < coords.size(); i++)
	{
		coords[i].x += u[3*i] ;
		coords[i].y += u[3*i+1] ;
		coords[i].z += u[3*i+2] ;
	}




	/*
	// update strains 
	for (auto & tri: tris)
	{
		tri.updateStrain(u);
	}
*/

}

void semiDiscreteMesh::initializeTest1(){
	coords.push_back(vector3D(0,0,0));
	coords.push_back(vector3D(0.1,0,0));
	coords.push_back(vector3D(0.1,0.1,0));
	coords.push_back(vector3D(0,0.1,0));
	coords.push_back(vector3D(0.05,0.05,0));

	tris.resize(4);
	std::vector<int> nodes(3);
	nodes[0] = 0; nodes[1] = 1; nodes[2] = 4;
	tris[0].setNodes(nodes);

	nodes[0] = 1; nodes[1] = 2; nodes[2] = 4;
	tris[1].setNodes(nodes);

	nodes[0] = 2; nodes[1] = 3; nodes[2] = 4;
	tris[2].setNodes(nodes);

	nodes[0] = 3; nodes[1] = 0; nodes[2] = 4;
	tris[3].setNodes(nodes);

	// find a and b values based on the orientation of the fibers
	for (auto& tri : tris)
	{
		tri.findaAndB(vector2D(1,0),vector2D(0,1),coords);
	}



	
}

void semiDiscreteMesh::initializeTest2(double rho, vector3D k_1, vector3D k_2){
	coords.push_back(vector3D(0,0,0));
	coords.push_back(vector3D(0.1,0,0));
	coords.push_back(vector3D(0.1,0.1,0));
	coords.push_back(vector3D(0,0.1,0));


	tris.resize(2);
	std::vector<int> nodes(3);
	nodes[0] = 0; nodes[1] = 1; nodes[2] = 3;
	tris[0].setNodes(nodes);

	nodes[0] = 1; nodes[1] = 2; nodes[2] = 3;
	tris[1].setNodes(nodes);



	// find a and b values based on the orientation of the fibers
	for (auto& tri : tris)
	{
		tri.findaAndB(vector2D(1,0),vector2D(0,1),coords);
	}

	BC bc;
	bc.ID = 2*3;bc.value = 6e-4;
	BCs.push_back(bc);
	bc.ID = 2*3+1;bc.value = 0;
	BCs.push_back(bc);
	bc.ID = 3*3;bc.value = 6e-4;
	BCs.push_back(bc);
	bc.ID = 3*3+1;bc.value = 0;
	BCs.push_back(bc);
	int nDOFs = coords.size()*3;
	F_int.resize(nDOFs); F_ext.resize(nDOFs);
	u.resize(nDOFs); v.resize(nDOFs); a.resize(nDOFs);
	M_D.resize(nDOFs);

		// Compute M_D matrix
	computeMassMatrix(rho);

	// find A and B for every triangle
	findAandB(k_1,k_2);
	
}


void writeVtkVectorArray(std::ostream &out, const std::string &name, const std::vector<vector3D> &data) {
	out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size(); c++) {
		
		out <<"			 "<< data[c].x << " " << data[c].y <<" " << data[c].z <<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void writeVtkVectorArray(std::ostream &out, const std::string &name, const std::vector<double> &data) {
	out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size()/3; c++) {
		
		out <<"			 "<< data[c*3] << " " << data[c*3+1] <<" " << data[c*3+2] <<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void writeVtkEpsilon11(std::ostream &out, const std::vector<SDTriangle> &data) {
	out << "        <DataArray type=\"Float64\" Name=\"" << "epsilon11" << "\" NumberOfComponents=\"1\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size(); c++) {
		
		out <<"			 "<< data[c].getEpsilon11()<<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void writeVtkEpsilon22(std::ostream &out, const std::vector<SDTriangle> &data) {
	out << "        <DataArray type=\"Float64\" Name=\"" << "epsilon22" << "\" NumberOfComponents=\"1\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size(); c++) {
		
		out <<"			 "<< data[c].getEpsilon22()<<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void writeVtkGamma(std::ostream &out, const std::vector<SDTriangle> &data) {
	out << "        <DataArray type=\"Float64\" Name=\"" << "gamma" << "\" NumberOfComponents=\"1\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size(); c++) {
		
		out <<"			 "<< data[c].getGamma()<<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void writeVtkXi1(std::ostream &out, const std::vector<SDTriangle> &data) {
	out << "        <DataArray type=\"Float64\" Name=\"" << "Xi1" << "\" NumberOfComponents=\"1\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size(); c++) {
		
		out <<"			 "<< data[c].getXi1()<<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void writeVtkXi2(std::ostream &out, const std::vector<SDTriangle> &data) {
	out << "        <DataArray type=\"Float64\" Name=\"" << "Xi2" << "\" NumberOfComponents=\"1\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size(); c++) {
		
		out <<"			 "<< data[c].getXi2()<<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void writeVtkk_1(std::ostream &out, const std::vector<SDTriangle> &data, std::vector<vector3D> & coords) {
	out << "        <DataArray type=\"Float64\" Name=\"" << "k_1" << "\" NumberOfComponents=\"3\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size(); c++) {
		vector3D k_1 = data[c].getk_1(coords);
		out <<"			 "<<  k_1.x << " " << k_1.y <<" " << k_1.z <<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void writeVtkk_2(std::ostream &out, const std::vector<SDTriangle> &data, std::vector<vector3D> & coords) {
	out << "        <DataArray type=\"Float64\" Name=\"" << "k_2" << "\" NumberOfComponents=\"3\" Format=\"ascii\">"<<std::endl;
	for (int c = 0; c < data.size(); c++) {
		vector3D k_2 = data[c].getk_2(coords);
		out <<"			 "<<  k_2.x << " " << k_2.y <<" " << k_2.z <<std::endl;
	}
	out << "		 </DataArray>" << std::endl;
}

void semiDiscreteMesh::printVTUfile(std::string fileName){
	std::ofstream out(fileName);
	// header
	out << "<VTKFile type=\"" << "UnstructuredGrid" << "\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
	
	out << "  <UnstructuredGrid>" << std::endl;
	out << "    <Piece NumberOfPoints=\"" << coords.size() << "\" NumberOfCells=\"" << tris.size() << "\" >" << std::endl;
	
	out << "      <Points>" << std::endl;
	writeVtkVectorArray(out, "Positions", coords);
	out << "      </Points>" <<std::endl;
	
	out << "      <Cells>" << std::endl; {
		out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl; {
			for (auto &elem : tris) {
				out<<"			";
				for (auto const &index : elem.getNodes()) {
					out << " " << index;
				}
				out<<std::endl;
			}
		} out << "		</DataArray>" << std::endl;

		out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"; {
			int endIndex = 0;
			for (auto const &elem : tris) {
				endIndex += 3;
				out << " " << endIndex;
			}
		} out << " </DataArray>" << std::endl;

		out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"; {
			for (auto &elem : tris) {
				out << " " << 5;
			}
		} out << " </DataArray>" << std::endl;
	} out << "      </Cells>" << std::endl;

	// write values
	out << "      <PointData Scalars= \"scalars\" \>" << std::endl;

	writeVtkVectorArray(out,"u",u);
	writeVtkVectorArray(out,"v",v);
	writeVtkVectorArray(out,"a",a);
	writeVtkVectorArray(out,"f",F_int);
	writeVtkVectorArray(out,"F_11",F_11);
	writeVtkVectorArray(out,"F_22",F_22);
	writeVtkVectorArray(out,"F_s",F_s);

	out << "      </PointData>" << std::endl;

	out << "      <CellData Scalars= \"scalars\" \>" << std::endl;
	writeVtkEpsilon11(out,tris);
	writeVtkEpsilon22(out,tris);
	writeVtkGamma(out,tris);
	writeVtkk_1(out,tris,coords);
	writeVtkk_2(out,tris,coords);
	writeVtkXi1(out,tris);
	writeVtkXi2(out,tris);

	out << "      </CellData>" << std::endl;

	//writeVtkCellData(out,coords.size(),45);
	//writeVtkVectorArray(out, "Velocity", time.velocities);
	//writeVtkVectorArray(out, "Vorticity", time.vorticities);
	

	out << "    </Piece>" << std::endl;
	out << "  </UnstructuredGrid>" << std::endl;
	out << "</VTKFile>" << std::endl;

}

void semiDiscreteMesh::findNeighbours(){
	// for each triangle finds the opposing nodes of the sides that are shared with other triangles

	// first fill a vector with all the sides of the triangles
	std::vector<side> sides(tris.size()*sidesPerElement);

	for (int i = 0; i < tris.size(); i++)
	{
		for (int j = 0; j < sidesPerElement; j++)
		{
			sides[i*sidesPerElement+j] = tris[i].getSide(j,i);
		}
	}

	// sort the array, the contiguous repeated sides will be neighbors, each side can be repeated at most two times

	std::sort(sides.begin(),sides.end());

	// loop through the sides array to fiind repeated sides
	for (int i = 0; i < sides.size()-1; i++)
	{
		if (sides[i].isEqual(sides[i+1]))
		{
			// these are sides of neighbouring triangles
			// set neighbour node of triangle i
			int nodeID = sides[i+1].oposingNode;
			int sideIDA = sides[i].sideID;
			int sideIDB = sides[i+1].sideID;
			tris[sides[i].triID].setNeighborNode(sideIDA,nodeID);
			tris[sides[i].triID].setNeighboringSide(sideIDA,sideIDB);
			tris[sides[i].triID].setNeighboringTri(sideIDA,sides[i+1].triID);
			// set neighbour node of triangle i+1
			nodeID = sides[i].oposingNode;
			tris[sides[i+1].triID].setNeighborNode(sideIDB,nodeID);
			tris[sides[i+1].triID].setNeighboringSide(sideIDB,sideIDA);
			tris[sides[i+1].triID].setNeighboringTri(sideIDB,sides[i].triID);
		}
	}


}

void semiDiscreteMesh::computeBendingVariables(){

#pragma omp parallel for 
	for (int i =0; i<tris.size(); i++)	//for (auto & tri :tris)
	{
		tris[i].computeBendingVars(coords);
		//tri.computeBendingVars(coords);
	}
}