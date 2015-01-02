clc; clear;
rng('shuffle');

symbols = ['a':'z' 'A':'Z' '0':'9'];
randName = symbols (randi(numel(symbols),[1 8]));

minn = log10([1E-5,1E-5,... % B2, U2
    1E-5,1E-5,1E-5,1E-5,1E-3,1E-3,... % 'xFwd1','xRev1','xFwd3','xRev3',AXLint1','AXLint2
    1E-3,1E-4,1E-1,10,1E-6, 1]); % 'kRec','kDeg','fElse','AXL2','Gas1'

maxx = log10([1,1E5,... % B2, U2
    1,1E5,1,1E5, 1, 1,... % 'xFwd1','xRev1','xFwd3','xRev3', AXLint1, AXLint2
    1E-1,1E-1,1,1E5,1E-2, 10]); % kRec, kDeg, fElse, AXL2, Gas

%%

Dopts = psoptimset('TimeLimit',60*60,'MaxFunEvals',1E10,'MaxIter',1E10,'Display','diagnose');

for ii = 1:200000000
    [paramOpt(ii,:),fitIDXglobal(ii)] = patternsearch(@cLibA549,minn + (rand(size(minn)) .* (maxx - minn)),...
        [],[],[],[],minn,maxx,[],Dopts);
    
    save(randName, 'paramOpt', 'fitIDXglobal', 'ii');
end

%%

paramOpt(fitIDXglobal > 23, :) = [];
fitIDXglobal(fitIDXglobal > 23) = [];

plot(paramOpt');
legend(arrayfun(@num2str, fitIDXglobal, 'unif', 0));
hold on;
plot(minn, 'k');
plot(maxx, 'k');
