% Diff profile test

clc; clear;

xx = linspace(0,1,50);

sshape = @(x) cos((10^x(1))*xx + x(3)).^x(4);

GasConc = -1.204;
Gass = @(x) 2*sshape(x)/mean(sshape(x) .* xx)*(10^x(2));


t = linspace(0,30,200);



frac = 0;
shapeParam = 60;
offset = 0;


D = ones(1,14)*0.1;
D(1) = 0;
D(3:7) = 0;
D(8:end) = 0;

D2 = ones(1,14)*100;

% A549 values
ffSpat = @(x) cLib_diff_profile_pYavg (t, 359.46, 0.0581, Gass(x), D, 1, 1, frac);

ffnonSpat = @(x) cLib_diff_profile_pYavg (t, 359.46, 0.0581, Gass(x), D2, 1, 1, frac);


hold off;
spat = ffSpat([3 GasConc offset shapeParam]);
spat = spat / spat(1);
plot(t, spat, 'r');

nonspat = ffnonSpat([-60 GasConc offset shapeParam]);
nonspat = nonspat / nonspat(1);
hold on;
plot(t, nonspat, 'b');
axis([min(t) max(t) 0 max(spat)]);