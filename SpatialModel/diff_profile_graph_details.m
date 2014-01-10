% Diff profile test

clc; clear;

sigmas = 1.5; % Above 1 weird effects


xx = linspace(0,1,100);
Gass = @(x) 2*(cos((10^x(1))*xx + x(3)) + 1)/mean((cos((10^x(1))*xx + x(3)) + 1) .* xx)*(10^x(2));
vv = [sigmas -0.699 0];
Gasv = Gass(vv);




D = 0*ones(1,14);
D(1) = 0;

t = 30;


out = cLib_diff_profile_pYavg (t, 7781, 0.0002, Gass(vv), D);

details = cLib_diff_profile (t, 7781, 0.0002, Gass(vv), D);


