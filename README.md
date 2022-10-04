# npoptband (under construction)

Parallel computing with GPU for bandwidth selection applied to nonparametric estimation of risk-neutral densities associated with option prices.
Implements MATLAB code, C code for sequential computing, and CUDA code for parallel computing.

# matlab

It includes code and data to perform bandwidth selection and risk-neutral density estimation using MATLAB code. There are two data sets, VIX and S&P500, for which two scripts vix.m and sp500.m, can be used to perform the calculations. For more information on how to run the scripts, within the directory exists a file README.md, and the script is anotated to facilitate eventual parameter tuning.
Sequential and parallel computing: Four directories are presented to test sequential and parallel for two datasets VIX and S&P500, vix, vixgpu, sp500, sp500gpu. Each directory contains the dataset, callprice, callstrike, callopenint (serve as weights), putprice, putstrike, putopenint (serve as weights), and strike, representing a series of unique strikes needed to implement the leave_k_i_estimate.

# vix

Sequential implementation for VIX. Applied to the sample included, the grid [0.75 2.0] was assumed for h_c and h_p. The underlying [10.0 47.5] range within a grid of 128 values was considered. These values can be changed in the main.cpp file. The app should normally run by changing these values accordingly and replacing the sample files with the same others using the same name. Future modifications will be implemented to allow a more straightforward way of defining the parameters

The code was compiled using

clang++ -std=c++11 -O2 -I"../include" runfunc.cpp main.cpp -o vix

Running the app for a h_c and h_p grid of 256

./vix 256

# sp500

For the S&P500, the same applies as with VIX, with grids [0.25 1.25] for h_c and h_p, and for the underlying [24.0 28.0]

The code was compiled using

clang++ -std=c++11 -O2 -I"../include" runfunc.cpp main.cpp -o sp500

Running the app for a h_c and h_p grid of 256

./sp500 256  (using the current dataset it can take several days to run)

# vixgpu and sp500gpu

The same datasets were used and the same grids. Running the code assume the same form, and the compiling was
performed as

nvcc -I"../include" main.cu -o vix          (in the vixgpu directory)

nvcc -I"../include" main.cu -o sp500       (in the sp500gpu directory)

running:

./vix 256    (using the current dataset, running time depending on the GPU, with a NVIDIA Tesla V100, 371 sec.)

./sp500 256  (using the current dataset, running time depending on the GPU, with a NVIDIA Tesla V100, 1370 sec.)


