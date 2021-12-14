# gpu_band (still in construction)
Parallel computing with GPU for bandwidth selection applied to nonparametric estimation of risk-neutral densities associated with option prices.

Implements MATLAB code, C code for sequential computing, and CUDA code for parallel computing.
The repository includes 3 main directories:

matlab

It includes code and data to perform bandwidth selection and risk-neutral density estimation 
using MATLAB code. There are two data sets, VIX and S&P500, for which to scripts vix.m and
sp500.m, can be used to perform the calculations. For more information on how to run the 
scripts, within the directory exists a file README.md, and the script are anotated to
facilitate eventual parameters tuning.

vix 

It includes code and data o perform bandwidth selection and risk-neutral density estimation 
using using C and CUDA, the first to a sequential implementation, whereas the second implements a sequential.
The code needs to be compiled through the compiler nvcc, and the choice of sequential (CPU) or 
parallel (GPU) is made through an input in running the application. 

Example:
- compile the code
     nvcc main.cu -o vix
          ./vix
          
          
     
- running the application
   - iterative introduction of the inputs
          
          ./vix

      
   - inputs given when the application is deployed
        ngrid number of elements in hc and hp grid
        hcmin lower value for hc
        hcmax upper vallue for hc
        hpmin lower valor for hp
        hpmax upper value for hp
        nx number of elements of the x grid
        xmin lower value for x
        xmax upper value for x
        cpugpu a flag for sequential or parallel -- 0 for cpu; 1 for gpu
        
          example for cpu
          ./vix 256 0.75 2.0 0.75 2.0 128 10.0 47.5 0
          example for gpu
          ./vix 256 0.75 2.0 0.75 2.0 128 10.0 47.5 1
        
        
sp500

The structure is in all identical to the vix, serves only to separate data sets.

Example:
- compile the code
     nvcc main.cu -o sp500
        
- running the application
   - iterative introduction of the inputs
       ./sp500
       
   -inputs given when the application is deployed
    
         example for cpu
        ./sp500 256 0.25 1.25 0.25 1.25 128 24.0 28.0 0
        
        example for gpu
        ./vix 256 0.25 1.25 0.25 1.25 128 24.0 28.0 1
        
        



