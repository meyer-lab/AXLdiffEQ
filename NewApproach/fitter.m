function error = fitter(paramIn, plott)

meas = [417618, 1011357, 817580, 850500, 1110500, 880520];
meas = meas / mean(meas);

params = [0.6 paramIn(1:5) 0.03 0.3 0.058 0.0022 0.1 0.5 623];

names = {'B2','U2','xFwd1','xRev1','xFwd3','xRev3','AXLint1','AXLint2',...
	'kRec','kDeg','fElse','endoArea','internalV'};

try
    outter = foldActivate(params, paramIn(6), paramIn(7));
    
    outter = outter / mean(outter);
    error = sqrt(sum((outter - meas).^2));
catch
    error = 1E6;
end

if plott == 1
    plot(meas,'g');
    hold on;
    plot(outter,'r');
end