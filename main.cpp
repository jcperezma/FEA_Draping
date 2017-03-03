#pragma once
#include "Fabric.h"
#include "tests.h"


int main(){


	Fabric glass("compressionTestmeshOUT.txt");

	//Fabric glass;
	//glass.initializeTest2();
	//
	int initialStep = 0;
	int finalStep = 1e6;
	int numFrames =1000;
	int frameStep =(finalStep-initialStep)/numFrames;

	int snapID=0;
	for (int i = 0; i < 1e6; i++)
	{
		glass.advanceStep(1e-6);
		std::cout<<i<<std::endl;

		if (i>=initialStep && i<finalStep && i%frameStep==0)
		{
			glass.print(snapID);
			snapID++;
		}
	}
	
	Fabric::printPVDfile(snapID);
	test();


	
	

return 0;
}