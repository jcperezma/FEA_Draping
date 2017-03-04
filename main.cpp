#pragma once
#include "Fabric.h"
#include "tests.h"
#include <omp.h>

#define NUM_OMP_THREADS 1

int main(){

	int numMaxThreads = omp_get_max_threads();
	int thradsToUse = NUM_OMP_THREADS;
	omp_set_num_threads(thradsToUse);
	std::cout<<"Max num of threads "<< numMaxThreads <<", threads being used "<< thradsToUse<<std::endl;
	Fabric glass("compressionTestmeshOUT.txt");

	//Fabric glass;
	//glass.initializeTest2();
	//
	/*
	int initialStep = 0;
	int finalStep = 1e6;
	int numFrames =1000;
	int frameStep =(finalStep-initialStep)/numFrames;


	double timeStart = omp_get_wtime();
	int snapID=0;
	for (int i = 0; i < 1000; i++)
	{
		glass.advanceStep(1e-6);
		std::cout<<i<<std::endl;

		if (i>=initialStep && i<finalStep && i%frameStep==0)
		{
			//glass.print(snapID);
			snapID++;
		}
	}
	
	double timeElapsed = omp_get_wtime()-timeStart;
	std::cout<<" Time elapsed: " << timeElapsed<< std::endl;
	//Fabric::printPVDfile(snapID);
	glass.print();

	*/

	double timeStart = omp_get_wtime();
	test();
	double timeElapsed = omp_get_wtime()-timeStart;
	std::cout<<" Time elapsed: " << timeElapsed<< std::endl;

	
	

return 0;
}