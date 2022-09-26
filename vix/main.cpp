//
//  main.cpp
//  agosto22
//
//  Created by Antonio Santos on 29/08/2022.
//
//
#include <stdio.h>
#include "nprnd.hpp"
#include "runfunc.hpp"

int main(int argc,const char* argv[])
{
	printf("Start .... \n");
	double r = 0.02;
	double tau = 1.0/12.0;
	//
	int grid[5] = {16,32,64,128,256};
	//
	for(int i=0;i<5;i++){
	   int ngrid = grid[i];
	   int nhc = ngrid;
	   int nhp = ngrid;
	   //
	   run_in_cpu(0.75, 2.0, nhc, 0.75, 2.0, nhp, 10.0, 47.5, 128, r, tau);
	}
	printf("Done ... \n");
	//
	return 0;
}
