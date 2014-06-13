function A549sensFull(IDX)
clc;
rng('shuffle');
slices = 25;

symbols = ['a':'z' 'A':'Z' '0':'9'];
randName = symbols (randi(numel(symbols),[1 8]));

minn = log10([0.0006,0.0006,1E-5,1E-5,... % 'B1','B2','U1','U2'
    1E-5,1E-5,1E-5,1E-5,... % 'xFwd1','xRev1','xFwd3','xRev3'
    1E-6,1E-6,1E-8,... % 'AXLint1','AXLint2','scaleA'
    1E-3,1E-4,1E-5,1E-1,1E-6,1]); % 'kRec','kDeg','fElse','fD2','Gas1','AXL2'
% 'AXLint2','scaleA' from 1E-3 and 1E-6

maxx = log10([6,6,1E5,1E5,... % 'B1','B2','U1','U2'
    1E5,1E5,1E5,1E5,... % 'xFwd1','xRev1','xFwd3','xRev3'
    1,100,1,... % 'AXLint1','AXLint2','scaleA'
    1,1,1,1,1,1E5]); % 'kRec','kDeg','fElse','fD2','Gas1','AXL2'
% Bumped up kDeg and kRec from 0.1

Dopts = psoptimset('Vectorized','on','TimeLimit',60*60,'CompletePoll','on','CompleteSearch','on',...
        'MaxFunEvals',1E10,'MaxIter',1E10,'Display','diagnose');

vv = linspace(minn(IDX),maxx(IDX),slices);

A = zeros(2,17);
A(1,3) = 1;
A(1,1) = -1;
A(2,4) = -1;
A(2,2) = 1;
b = [1; 0];

for ii = 1:slices
    minn2 = minn; minn2(IDX) = vv(ii);
    maxx2 = maxx; maxx2(IDX) = vv(ii); 
    
    params = minn2 + (rand(size(minn2)) .* (maxx2 - minn2));
    
    [paramOpt{ii},fitIDXglobal(ii)] = patternsearch(@cLibA549,params,A,b,[],[],minn2,maxx2,[],Dopts);

    save([randName '-' mat2str(IDX)]);
end