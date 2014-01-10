% Diff profile test

% This script generated Figure 2B
clc; clear;

xx = linspace(0,1,90);
sshape = @(x) cos(xx*pi/3).^x(1);

GasConc = -1.204;
Gass = @(x) 2*sshape(x)/mean(sshape(x).*xx)*(10^x(2));

t = 30;
frac = 0;

D = ones(1,14)*1;
D(1) = 0;
D(3:7) = 0;
D(8:end) = 0;

ffSpat = @(x,y) cLib_diff_profile_pY (t, 359.46, 0.0581, Gass(x), D*y, 1, 1, frac);
ffnonSpat = @(x) cLib_diff_profile_pY (t, 359.46, 0.0581, Gass(x), zeros(1,14), 1, 1, frac);

A = logspace(-2,3,20);

B = zeros(length(A),1);
B2 = B;
D = B;
C = B;

for ii = 1:length(A)
    temp = ffSpat([A(ii) GasConc],0.1);
    C(ii) = temp(1);
    temp = ffSpat([A(ii) GasConc],1);
    D(ii) = temp(1);
    temp = ffSpat([A(ii) GasConc],10);
    B(ii) = temp(1);
    temp = ffnonSpat([A(ii) GasConc]);
    B2(ii) = temp(1);
    
    hold off;
    semilogx(A,B./B(1),'r-');
    hold on;
    semilogx(A,C./B(1),'b-');
    semilogx(A,D./B(1),'b-');
    semilogx(A,B2./B(1),'b-');
    drawnow;
end
