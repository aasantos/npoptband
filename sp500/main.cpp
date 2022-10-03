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
	//
	double r = 0.02;
	double tau = 1.0/12.0;
	//
	int ngrid = atoi(argv[1]);
	int nhc = ngrid;
	int nhp = ngrid;
	run_in_cpu(0.25, 1.25, nhc, 0.25, 1.25, nhp, 24.0, 28.0, 128, r, tau);
	//
	printf("Done ... \n");
	//
	return 0;
}
