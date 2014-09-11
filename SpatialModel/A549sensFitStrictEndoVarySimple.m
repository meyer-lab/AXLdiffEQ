%function A549sensFitStrictEndoVarySimple(IDX)
clc;
rng('shuffle');
slices = 80;

symbols = ['a':'z' 'A':'Z' '0':'9'];
randName = symbols (randi(numel(symbols),[1 8]));

minn = log10([1.2,0.01,0.042,1E-5,... % 'B1','B2','U1','U2'
    1E-5,1E-5,1E-5,1E-5,... % 'xFwd1','xRev1','xFwd3','xRev3'
    1E-6,1E-3,1E-6,... % 'AXLint1','AXLint2','scaleA'
    1E-3,1E-4,1E-5,1,1E-6,1,0.01,10]); % 'kRec','kDeg','fElse','fD2','Gas1','AXL2'

maxx = log10([1.2,1.2,0.042,1E5,... % 'B1','B2','U1','U2'
    1E5,1E5,1E5,1E5,... % 'xFwd1','xRev1','xFwd3','xRev3'
    1,1,1,... % 'AXLint1','AXLint2','scaleA'
    0.1,0.1,1,1,0.1,1E5,100,1E4]); % 'kRec','kDeg','fElse','fD2','Gas1','AXL2'

Dopts = psoptimset('Vectorized','on','TimeLimit',8*60*60,'CompletePoll','on','CompleteSearch','on',...
        'MaxFunEvals',1E10,'MaxIter',1E10,'Display','diagnose');

Gopts = gaoptimset('Display','iter','Vectorized','on','PopulationSize',1e4);
    
A = zeros(2,length(minn));
A(1,3) = 1;
A(1,1) = -1;
A(2,4) = -1;
A(2,2) = 1;
b = [1; 0];

Aeq = zeros(1,length(minn));
A(1,5) = 1;
A(1,7) = -1;
beq = 0;

for ii = 1:slices
    
    params = minn + (rand(size(minn)) .* (maxx - minn));
    
    try
        %[paramOpt{ii},fitIDXglobal(ii)] = patternsearch(@cLibA549endoVary,params,A,b,Aeq,beq,minn,maxx,[],Dopts);
        ga(@cLibA549endoVary,length(minn),[],[],[],[],minn,maxx,[],Gopts);
    catch
        fitIDXglobal(ii) = 1E6;
    end

    save(['./StrictSensEndoVary/' randName '-' mat2str(IDX)]);
end