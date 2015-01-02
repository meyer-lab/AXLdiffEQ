function outter = ensembleOpt (pps, params)

pYs = zeros(size(params,1),1);

for ii = 1:size(params,1)
    pY = cLib_profile ([0 10], params(ii,:).*(10.^pps'), 10, 0);

    if ~isempty(pY)
        pYs(ii) = pY(2) / pY(1);
    end
end

outter = sum((pYs - 5).^2);
