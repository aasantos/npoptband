%%
clc;clear;
load vixData.mat callstrike callprice callopenint ...
                 putstrike putprice putopenint
%%
%
index = callstrike > 9 & callstrike < 51;
callstrike = callstrike(index);
callprice = callprice(index);
callopenint = callopenint(index);
index = putstrike > 9 & putstrike < 51;
putstrike = putstrike(index);
putprice = putprice(index);
putopenint = putopenint(index);
%
%
%% Optimal bandwith choice (It may take several hours to run)
%
%  It uses the function optimalbandwidth function, which has as input a 
%  structure. The function outputs a vector with the optimal values for the
%  h_c and h_p, and the respective matrices used to calculate such values.
%
S.callprice = callprice;
S.callstrike = callstrike;
S.callopenint = callopenint;
%
S.putprice = putprice;
S.putstrike = putstrike;
S.putopenint = putopenint;
%
S.r = 0.02;
S.tau = 1/12;
%
S.hcmin = 0.75;
S.hcmax = 2.0;
S.hpmin = 0.75;
S.hpmax = 2.0;
S.ngrid = 32;
%
tic
[hoptim,cvmat,areamat,diffmat,entropymat,matcrit] = optimalbandwidth(S);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
