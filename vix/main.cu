#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/io.cuh"
#include "../include/vector.cuh"
#include "../include/ls.cuh"
#include "../include/qp.cuh"
#include "../include/rndestim.cuh"
#include "../include/kernelfunc.cuh"
#include "../include/runfunc.cuh"

int main(int argc,char *argv[])
{
std::cout<<"Starting estimation ...\n";
//
	//int ngrid = 32;
	//int nhc = ngrid;
	//int nhp = ngrid;
	//float hcmin = 0.75;
	//float hcmax = 2.0;
	//float hpmin = 0.75;
	//float hpmax = 2.0;
	//int nx = 128;
	//float xmin = 10.0;
	//float xmax = 47.5;
	//int cpugpu = 0;
//
int ngrid;
int nhc;
int nhp;
float hcmin;
float hcmax;
float hpmin;
float hpmax;
int nx;
float xmin;
float xmax;
int cpugpu;
//
if(argc == 2){
	ngrid = atoi(argv[1]);
  	nhc = ngrid;
  	nhp = ngrid;
  	std::cout<<"hc min: ";
  	std::cin>>hcmin;
  	std::cout<<"hc max: ";
  	std::cin>>hcmax;
  	std::cout<<"hp min: ";
  	std::cin>>hpmin;
  	std::cout<<"hp max: ";
  	std::cin>>hpmax;
  	std::cout<<"nx (int): ";
  	std::cin>>nx;
  	std::cout<<"xmin: ";
  	std::cin>>xmin;
  	std::cout<<"xmax: ";
  	std::cin>>xmax;
  	std::cout<<"cpu/gpu (0/1/2): ";
  	std::cin>>cpugpu;
}else{
	if(argc == 4){
		ngrid = atoi(argv[1]);
  		nhc = ngrid;
  		nhp = ngrid;
  		hcmin = atof(argv[2]);
  		hcmax = atof(argv[3]);
  		std::cout<<"hp min: ";
  		std::cin>>hpmin;
  		std::cout<<"hp max: ";
  		std::cin>>hpmax;
  		std::cout<<"nx (int): ";
  		std::cin>>nx;
  		std::cout<<"xmin: ";
  		std::cin>>xmin;
  		std::cout<<"xmax: ";
  		std::cin>>xmax;
  		std::cout<<"cpu/gpu (0/1): ";
  		std::cin>>cpugpu;
	}
	if(argc == 6){
		ngrid = atoi(argv[1]);
  		nhc = ngrid;
  		nhp = ngrid;
  		hcmin = atof(argv[2]);
  		hcmax = atof(argv[3]);
  		hpmin = atof(argv[4]);
  		hpmax = atof(argv[5]);
  		std::cout<<"nx (int): ";
  		std::cin>>nx;
  		std::cout<<"xmin: ";
  		std::cin>>xmin;
  		std::cout<<"xmax: ";
  		std::cin>>xmax;
  		std::cout<<"cpu/gpu (0/1/2): ";
  		std::cin>>cpugpu;
	}
	if(argc == 10){
		ngrid = atoi(argv[1]);
  		nhc = ngrid;
  		nhp = ngrid;
  		hcmin = atof(argv[2]);
  		hcmax = atof(argv[3]);
  		hpmin = atof(argv[4]);
  		hpmax = atof(argv[5]);
  		nx = atoi(argv[6]);
  		xmin = atof(argv[7]);
  		xmax = atof(argv[8]);
  		cpugpu = atoi(argv[9]);
	}
	else{
  		std::cout<<"grid size: ";
  		std::cin>>ngrid;
  		nhc = ngrid;
  		nhp = ngrid;
  		std::cout<<"hc min: ";
  		std::cin>>hcmin;
  		std::cout<<"hc max: ";
  		std::cin>>hcmax;
  		std::cout<<"hp min: ";
  		std::cin>>hpmin;
  		std::cout<<"hp max: ";
  		std::cin>>hpmax;
  		std::cout<<"nx (int): ";
  		std::cin>>nx;
  		std::cout<<"xmin: ";
  		std::cin>>xmin;
  		std::cout<<"xmax: ";
  		std::cin>>xmax;
  		std::cout<<"cpu/gpu (0/1/2): ";
  		std::cin>>cpugpu;
	}
}
//
if(cpugpu == 0)
{
  run_in_cpu(hcmin,hcmax,nhc,hpmin,hpmax,nhp,xmin,xmax,nx);
}
if(cpugpu == 1){
  run_in_gpu(hcmin,hcmax,nhc,hpmin,hpmax,nhp,xmin,xmax,nx);
}
if(cpugpu == 2){
  run_in_gpu_shared(hcmin,hcmax,nhc,hpmin,hpmax,nhp,xmin,xmax,nx);
}
//
//
return 0;
}
