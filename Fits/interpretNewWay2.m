clc; clear;

load widertEXKsu
fitStruct2 = fitStruct;
load widerY4h4NgrH
fitStruct = [fitStruct, fitStruct2];

names ={'B2', 'U2', 'xFwd1','xRev1','xFwd3','xRev3', 'AXLint1','AXLint2',...
    'kRec','kDeg','fElse','AXL2','Gas1','picker'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

cutoff = min(fitt)+2;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

%IDXrem = params(:,7) < params(:,6);

%fitt(IDXrem) = [];
%params(IDXrem,:) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);

clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

params = [ones(size(params,1),1)*0.06, 10.^params];
params(:,end) = params(:,end) > 0.5;

% % % % % % % % % % % % % % %

tps = linspace(0,10,100);
stimProfile = @(x, pp) cLib_profile (tps, x, 12.5, 1);

xxx = 1000;

for ii = 1:size(params,1)
    temp = stimProfile(params(ii,:));
    
    if ~isempty(temp)
        B(ii,1:100) = temp / temp(1);
    end
    
    params2 = params(ii,:);
    params2(3) = params2(3)*xxx;
    params2(5) = params2(5)*xxx;

    temp = stimProfile(params2);
    
    if ~isempty(temp)
        C(ii,1:100) = temp / temp(1);
    end
end

subplot(1,3,1)
plot(tps,B');
axis([0 10 0 2])
subplot(1,3,2)
plot(tps,C');
axis([0 10 0 2])
subplot(1,3,3)
plot(tps,D');
axis([0 10 0 2])