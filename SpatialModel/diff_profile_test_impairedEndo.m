% Diff profile test

% This script generated Figure XX

clc; clear;

xx = linspace(0,1,10);
sshape = @(x) cos((10^x(1))*xx).^x(3);
Gass = @(x) 2*sshape(x)/mean(sshape(x) .* xx)*(10^x(2));

t = 30;
frac = 0;

D = ones(1,14)*0.1;
D(1) = 0;
D(3:7) = 0;
D(8:end) = 0;


ffnonSpat = @(x) cLib_diff_profile_pYavg (t, 359.46, 0.0581, Gass([-100 -1.204 2]), ones(1,14)*100, x, 1, frac);

A = logspace(-2,0,100);
B = ones(length(A),1);

for ii = 1:length(A)

    B(ii) = ffnonSpat(A(ii));
    
    semilogx(A,B/B(end),'r-');
    drawnow;
end
