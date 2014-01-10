function fit_full_sensitivity_analysis(ii)

params = log10(getOptimParams(1));

names = {'B1', 'B2', 'U1', 'U2', 'xFwd1', 'xRev1', 'xFwd3', 'xRev3',...
    'AXLint1', 'AXLint2', 'scaleA', 'kRec', 'kDeg', 'fElse', 'fD2',...
    'Gas1', 'Gas2', 'AXL1', 'AXL2'};

minn = log10([1.2, 1E-6, 0.042, 1E-5,... % B1, B2, U1, U2
        1E-5,1E-5,1E-5,1E-5,... % xFwd1, xRev1, xFwd3, xRev3
1E-6, 1E-3,... % AXLint1, AXLint2
1E-6,... % scaleA
1E-3, 1E-4,1E-5,1E-1,... % kRec, kDeg, fElse, fD2
1E-6,1E-6,1,1]); % Gas1, Gas2, AXL1, AXL2

maxx = log10([1.2, 0.6, 0.042, 1E5... % B1, B2, U1, U2
1E5,1E5,1E5,1E5,... % xFwd1, xRev1, xFwd3, xRev3
1, 100,... % AXLint1, AXLint2
1,... % scaleA
0.1,0.1,1,1,...
100, 100, 1E5, 1E5]); % Gas1, Gas2, AXL1, AXL2

N = 25;

vv = linspace(max(params(ii)-1,minn(ii)),min(params(ii)+1,maxx(ii)),N);

for jj = 1:N
    A = zeros(1,19);
    b = vv(jj);
    A(ii) = 1;

    params2 = params;
    params2(ii) = vv(jj);

    paramsIDX{jj} = params2;
    AIDX{jj} = A;
    bIDX{jj} = b;
    fitIDX(jj) = cLib(params2);
end

fitIDXlocal = fitIDX;
fitIDXglobal = fitIDX;

Fopts = optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',1e6);
Dopts = psoptimset('Vectorized','on','TimeLimit',1600);

for jj = 1:N
    [~,fitIDXlocal(jj)] = fmincon(@cLib,paramsIDX{jj},[],[],AIDX{jj},bIDX{jj},minn,maxx,[],Fopts);
    
    [~,fitIDXglobal(jj)] = patternsearch(@cLib,paramsIDX{jj},[],[],AIDX{jj},bIDX{jj},minn,maxx,[],Dopts);
    
    
    save(mat2str(ii));
    
end





