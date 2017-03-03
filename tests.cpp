#include "tests.h"


void singleTriangle(){
	std::vector<vector3D> coords;

	coords.push_back(vector3D(0,0,0));
	coords.push_back(vector3D(0.5,0,0));
	coords.push_back(vector3D(0.25,0.3,0));

	SDTriangle tri;

	std::vector<int> nodes;
	nodes.push_back(0);
	nodes.push_back(1);
	nodes.push_back(2);

	tri.setNodes(nodes);
	vector3D k1(1,0,0);
	vector3D k2(0,1,0);

	vector3D h1 = k1 / norm(k1);
	vector3D h2 = k2 / norm(k1);

	tri.findaAndB(h1,h2,coords);

	std::vector<double> M_D(9);
	std::vector<double> F(9);
	tri.computeMassMatrix(coords,M_D,10);
	//tri.computeInternalForce(46000,48000,0.37044,0.84105,1.0311,F,coords);
}

void test(){

	semiDiscreteMesh mesh;
	mesh.initializeTest1();


}