# npoptband (under construction)
Parallel computing with GPU for bandwidth selection applied to nonparametric estimation of risk-neutral densities associated with option prices.

Implements MATLAB code, C code for sequential computing, and CUDA code for parallel computing.

matlab

It includes code and data to perform bandwidth selection and risk-neutral density estimation 
using MATLAB code. There are two data sets, VIX and S&P500, for which two scripts vix.m and
sp500.m, can be used to perform the calculations. For more information on how to run the 
scripts, within the directory exists a file README.md, and the script are anotated to
facilitate eventual parameters tuning.

vix 

It includes code and data o perform bandwidth selection and risk-neutral density estimation 
using using C and CUDA, the first to a sequential implementation, whereas the second implements a sequential.
The code needs to be compiled through the compiler nvcc, and the choice of sequential (CPU) or 
parallel (GPU) is made through an input in running the application. 
        
        
sp500
