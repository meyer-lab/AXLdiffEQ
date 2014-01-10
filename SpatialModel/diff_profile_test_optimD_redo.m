% Diff profile test
% This script generated Figure 2B
clc; clear;

xx = linspace(0,1,80);
sshape = @(x) cos(xx*pi/3).^x(1);

GasConc = -1.204;
Gass = @(x) 2*sshape(x)/mean(sshape(x).*xx)*(10^x(2));

t = 30;
frac = 0;

D = ones(1,14)*1;
D(1) = 0;
D(3:7) = 0;
D(8:end) = 0;

ffSpat = @(x,y) cLib_diff_profile_pYavg (t, 359.46, 0.0581, Gass(x), D*y, 1, 1, frac);
ffnonSpat = @(x) cLib_diff_profile_pYavg (t, 359.46, 0.0581, Gass(x), zeros(1,14), 1, 1, frac);

A = logspace(-2,3,20);

B = zeros(length(A),1);
B2 = B;
D = B;
C = B;

for ii = 1:length(A)
    D(ii) = ffSpat([A(ii) GasConc],0.1);
    B(ii) = ffSpat([A(ii) GasConc],1);
    C(ii) = ffSpat([A(ii) GasConc],10);
    B2(ii) = ffnonSpat([A(ii) GasConc]);
    
    hold off;
    semilogx(A,B./B(1),'r-');
    hold on;
    semilogx(A,D./D(1),'g-');
    semilogx(A,C./C(1),'k-');
    semilogx(A,B2./B(1),'b-');
    drawnow;
end