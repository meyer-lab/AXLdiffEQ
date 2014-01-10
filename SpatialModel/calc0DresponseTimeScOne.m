clc; clear;
params = getOptimParams(1);
TPS = linspace(0,10,1000);

    
YY = cLib_profile (TPS, params, 359.46, 0.0581, 0.25, 0);
plot(TPS,YY/YY(1));
