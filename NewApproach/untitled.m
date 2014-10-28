clc; clear;

fitLog = @(x) cLib_profile(10.^x);

minn = [-2 -5 -5 -5 -5 -4 1];
maxx = [5 5 5 5 5 0 3];

params = minn + ((maxx - minn) .* rand(size(maxx)));


A = [0 0 -1  0 1 0 0;
     0 1  0 -1 0 0 0];
b = [0; 0];

options = optimset('Display','iter');

x = patternsearch(fitLog,params,A,b,[],[],minn,maxx,[],options)

fitter(10.^x, 1)