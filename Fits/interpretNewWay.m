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

IDXrem = params(:,7) < params(:,6);

fitt(IDXrem) = [];
params(IDXrem,:) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);

clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

params = [ones(size(params,1),1)*0.06, 10.^params];
params(:,end) = params(:,end) > 0.05;

% % % % % % % % % % % % % % %

xxx = logspace(-2.0,2.0,20);

stimProfile = @(x, pp) cLib_profile (240, x, 10, 1);

outter = zeros(size(params,1),9,length(xxx));

%IDXX = [1:10 12:13];

IDXX = 3:4;

for ii = 1:size(params,1)
    B = stimProfile(params(ii,:));
    
    for jj = 1:length(IDXX)
        parfor xx = 1:length(xxx)
            params2 = params(ii,:);
            params2(IDXX(jj)) = params2(IDXX(jj))*xxx(xx);
            params2(IDXX(jj)+2) = params2(IDXX(jj)+2)*xxx(xx);
            %params2(IDXX(jj)+4) = params2(IDXX(jj)+4)*xxx(xx);
            
            C = stimProfile(params2);

            if ~isempty(C) && ~isempty(B)
                outter(ii,jj,xx) = C / B;
            else
                outter(ii,jj,xx) = NaN;
            end
        end
        
        
    end 
    
    if mod(ii,10) == 0
        for jj = 1:length(IDXX)
            subplot(3,4,jj);
            loglog(xxx,squeeze(outter(:,jj,:))');
            title(names(IDXX(jj)));
            axis([min(xxx) max(xxx) 0.01 100]);
        end
        drawnow;
    end
end
