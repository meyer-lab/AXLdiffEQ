% Diff profile test

clc; clear;

xx = linspace(0,1,100);
Gass = @(x) 2*(cos((10^x(1))*xx + x(3)) + 1)/mean((cos((10^x(1))*xx + x(3)) + 1) .* xx)*(10^x(2));
D = 0*ones(1,14);
D(1) = 0;
D(2) = 0.1;


sigmas = logspace(-2,1.5,10); % Above 1 weird effects

for ii = 1:length(sigmas)
    vv = [sigmas(ii) 0 0];
    t = 10;
    
    
    out(ii) = cLib_diff_profile_pYavg (t, 7781, 0.0002, Gass(vv), D);
    
    %details = mean(cLib_diff_profile (t, 7781, 0.0002, Gass(vv), D));
    %totAXL(ii) = sum(details) + details(8) + details(9);
    
    plot(out);
    drawnow;
end
