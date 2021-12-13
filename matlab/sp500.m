%%
clc;clear;
load sp500Data.mat callstrike callprice callopenint ...
                 putstrike putprice putopenint
%%
%
index = callstrike > 23 & callstrike < 29;
callstrike = callstrike(index);
callprice = callprice(index);
callopenint = callopenint(index);
index = putstrike > 23 & putstrike < 29;
putstrike = putstrike(index);
putprice = putprice(index);
putopenint = putopenint(index);
%%
r = 0.02;
tau = 1/12;
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
S.r = r ;
S.tau = tau;
%
S.hcmin = 0.15;
S.hcmax = 0.9;
S.hpmin = 0.15;
S.hpmax = 0.9;
S.ngrid = 64;
%
tic
[hoptim,cvmat,areamat,diffmat,entropymat,matcrit] = optimalbandwidth(S);
toc
%%%%%%
