% Diff profile test

% This script generated Figure 2B

clc; clear;

xx = linspace(0,1,101);
sshape = @(x) cos((10^x(1))*xx).^x(3);

GasConc = -1.204;
Gass = @(x) 2*sshape(x)/mean(sshape(x) .* xx)*(10^x(2));

t = 30;
frac = 0;

D = ones(1,14)*0.1;
D(1) = 0;
D(3:7) = 0;
D(8:end) = 0;

ffSpat = @(x) cLib_diff_profile_pY (t, 359.46, 0.0581, Gass(x), D, 1, 1, frac);
ffnonSpat = @(x) cLib_diff_profile_pY (t, 359.46, 0.0581, Gass(x), zeros(1,14), 1, 1, frac);
ffnonSpat2 = @(x) cLib_diff_profile_pY (t, 359.46, 0.0581, Gass([-100 x(2) x(3)]), ones(1,14)*100, 1, 1, frac);

A = linspace(-1.5,1.5,40);
shapeParam = 60;

B = zeros(length(A),1);
B2 = B;
B3 = B;

for ii = 1:length(A)
    temp = ffSpat([A(ii) GasConc shapeParam]);
    B(ii) = temp(1);
    temp = ffnonSpat([A(ii) GasConc shapeParam]);
    B2(ii) = temp(1);
    
    hold off;
    plot(A,B./B(1),'r-');
    hold on;
    plot(A,B2./B(1),'b-');
    title(['Shape parameter ' mat2str(shapeParam)]);
    drawnow;
end
