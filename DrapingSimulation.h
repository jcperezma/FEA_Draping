#pragma once
#include "Fabric.h"
#include <vector>
#include <sstream> 
#include <fstream>
#include <string>
#include <iostream>
#include <memory>

class DrapingSimulation
{
public:
	DrapingSimulation(){};
	~DrapingSimulation(){};
	void initilizeFabricFromFile(std::string fileName);

private:
	std::vector<Fabric> fabrics;
	double dt;	

};

