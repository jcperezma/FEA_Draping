#include "DrapingSimulation.h"




void DrapingSimulation::initilizeFabricFromFile(std::string fileName){
	Fabric newFabric;
	newFabric.initilizeFromFile(fileName);
	fabrics.push_back(newFabric);

}

