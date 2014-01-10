% Diff profile test

% This script generated Figure 2D

clc; clear;

xx = linspace(0,1,100);
sshape = @(x) cos((10^x(1))*xx).^x(3);

GasConc = -1.204;
Gass = @(x) 2*sshape(x)/mean(sshape(x) .* xx)*(10^x(2));


meanify = @(x) 2*x'*xx'/length(xx);




t = 30;
frac = 0;

D = ones(1,14)*0.1;
D(1) = 0;
D(3:7) = 0;
D(8:end) = 0;

ffSpat = @(x) cLib_diff_profile (t, 359.46, 0.0581, Gass(x), D, 1, 1);

A = linspace(-1.5,1.5,20);
shapeParam = 60;



for ii = 1:length(A)
    B(ii,:) = meanify(ffSpat([A(ii) GasConc shapeParam]));
    ii
end



%%


for ii = 1:size(B,2)
    
    C(:,ii) = B(:,ii) / B(1,ii);
    
end