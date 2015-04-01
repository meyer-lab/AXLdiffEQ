clc; clear;

load WithpYnewBalance_3YD;

names = {'U2','xFwd1','xRev3','AXLint1','AXLint2','kRec','kDeg','fElse','AXL','Gas'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

cutoff = min(fitt)+6;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);

clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

params = 10.^params;
params(:,end) = params(:,end) > 0.05;

% % % % % % % % % % % % % % %

xxx = logspace(-2.0,2.0,20);

stimProfile = @(x, pp) cLib_profile (240, x, 10, 1);

outter = zeros(size(params,1),9,length(xxx));

IDXX = [1:7 9:10];

for ii = 1:size(params,1)
    B = stimProfile(params(ii,:));
    
    for jj = 1:length(IDXX)
        for xx = 1:length(xxx)
            params2 = params(ii,:);
            params2(IDXX(jj)) = params2(IDXX(jj))*xxx(xx);
            
            C = stimProfile(params2);

            if ~isempty(C) && ~isempty(B)
                outter(ii,jj,xx) = C / B;
            else
                outter(ii,jj,xx) = NaN;
            end
        end
        
        
    end 
    
    for jj = 1:length(IDXX)
        subplot(3,4,jj);
        loglog(xxx,squeeze(outter(ii,jj,:))','Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
        hold on;
        title(names(IDXX(jj)));
        axis([min(xxx) max(xxx) 0.01 100]);
    end
    drawnow;
end
