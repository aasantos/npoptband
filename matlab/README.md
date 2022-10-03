# Bandwidth selection for NP estimation of RND -- MATLAB example
Estimate RND functions using option prices
Data available: sample of data from VIX and S&P500 options contracts, from April 16 to April 20, 2018, 
with maturity date at May 18, 2018. The MATLAB data format files are vixData.mat and sp500Data.mat
There is a series of functions used to implement the RND estimation and respective analysis:

areadensity.m <br />
entropy.m <br />
minMatrix.m <br />
callputfunctionf.m <br />
npcallputoptimLG.m <br />
optimalbandwidth.m <br />

#

Two scripts (functions) are defined applied to do the analysis for VIX and S&P500

vix.m <br />
sp500.m <br />

By running the scripts (functions), it loads the data, calculates the optimal
bandwidth values, and estimate the respective RNDs. The results are
compared with a series of graphic analysis.
