clc; clear;

load 77j87EXi
fitStruct3 = fitStruct;
load SwiderY4h4NgrH
fitStruct2 = fitStruct;
load widerY4h4NgrH
fitStruct = [fitStruct, fitStruct2, fitStruct3];


names ={'U2', 'xFwd1','xRev1','xFwd3','xRev3', 'AXLint1','AXLint2',...
    'kRec','kDeg','fElse','AXL2','Gas1','picker'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

params(fitt > min(fitt)+1,:) = [];
fitt(fitt > min(fitt)+1) = [];

params = [ones(size(params,1),1)*0.06, 10.^params];
params(:,end) = params(:,end) > 0.05;

%%
clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

% for ii = 1:13
%     subplot(3,5,ii);
%     
%     plot(params(:,ii), fitt,'o');
%     title(names(ii));
%     axis([minn(ii) maxx(ii) min(fitt) max(fitt)]);
% end
%%

% xx = linspace(0,1,50);
% GasConc = 10;
% 
% GassF = @(x, au) (xx < x)*GasConc; 
% 
% D = zeros(1,14);
% D(2) = 1;
% D(8:end) = 0.1;
% 
% %figure(2);
% 
% params = [ones(size(params,1),1)*0.06, 10.^params];
% params(:,end) = params(:,end) > 0.05;
% 
% %parfor_progress(size(params,1));
% 
% for ii = 1:size(params,1)
%     B = cLib_diff_profile_pYavg (linspace(0,120,100), params(ii,:), GassF(0.04, params(ii,end-1)), D, 1);
%     %parfor_progress;
%     
%     try
%         B(end) / B(1)
%         outter(ii) = B(end) / B(1); %#ok<SAGROW>
%         
%         plot(sort(outter));
%         drawnow;
%     catch
%         outter(ii) = NaN;
%     end
%     
% end

%%

Gass = logspace(log10(0.25),log10(64),20);

for ii = 1:size(params,1)
    for jj = 1:length(Gass)
        stimProfile = cLib_profile ([0 10 240], params(ii,:), Gass(jj), 0);

        if numel(stimProfile) == 0
            outHour(ii,jj) = nan;
            outFour(ii,jj) = nan;
            continue;
        end

        outHour(ii,jj) = stimProfile(2) / stimProfile(1);
        outFour(ii,jj) = stimProfile(3) / stimProfile(1);
    end
end

semilogx(Gass,outHour','b');
hold on;
semilogx(Gass,outFour','k');
axis([min(Gass) max(Gass) 0 max(1,max(max([outHour outFour])))]);



