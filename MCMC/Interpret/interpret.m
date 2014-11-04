clc; clear;

params = flipud(csvread('8LV41B.csv'));
sshape = @(x,xx) cos(xx*pi/3).^x(1);
xx = linspace(0,1,50);
shapeParam = 60;
GasConc = 16;

GassF = @(x,xx) (xx < x)*GasConc;

D = ones(1,14);
D(2) = 1;

for ii = 1:size(params,1)
    
    B = cLib_diff_profile_pYavg ([0 240], 10.^params(ii,:), GassF(0.05,xx), D, 1, 1, 0);
    outter(ii) = B(2) / B(1);
    
    plot(outter);
    drawnow;
end