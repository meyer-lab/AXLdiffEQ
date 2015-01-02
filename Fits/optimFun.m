function outter = optimFun(params, xx)

pYbefore = cLib_profile ([0 0.01], params, 1.25*10^(xx(2)), 0);

params(3) = params(3).*(10.^xx(1));
params(5) = params(5).*(10.^xx(1));

pY = cLib_profile (linspace(0,10,10), params, 1.25*10^(xx(2)), 0);
pY2 = cLib_profile (linspace(0,10,10), params, 0.25*10^(xx(2)), 0);

if numel(pY) > 0 && numel(pY2) > 0
    outter = abs(pY(1)/pYbefore(1) - 1) * 10;
    
    pY = pY / pY(1);
    pY2 = pY2 / pY2(1);
    
    outter = outter + sum(abs(pY - linspace(1,3,length(pY))));
    outter = outter + sum(abs(pY2 - linspace(1,1,length(pY2))));
else
    outter = 1E6;
end
