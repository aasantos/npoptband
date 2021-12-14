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



