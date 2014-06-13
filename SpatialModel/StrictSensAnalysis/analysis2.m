clc; clear;

names = {'B1','B2','U1','U2','xFwd1','xRev1','xFwd3','xRev3',...
    'AXLint1','AXLint2','scaleA','kRec','kDeg','fElse','fD2',...
    'Gas1','AXL2'};

names{18} = 'KD1';
names{19} = 'KD2';
names{20} = 'K1';
names{21} = 'K2';
names{22} = 'KD_2-KD_1';
names{23} = 'K_2-K_1';

filess = dir('*.mat');

IDX = 1;

for ii = 1:length(filess)
    
    load(filess(ii).name,'fitIDXglobal');
    load(filess(ii).name,'paramOpt');
    
    for jj = 1:length(fitIDXglobal)
        if (fitIDXglobal(jj) < 1E4)
            
            fitt(IDX) = fitIDXglobal(jj);
            gParamOpt(IDX,:) = paramOpt{jj};
            
            IDX = IDX + 1;
        end
    end
end

gParamOpt = [gParamOpt gParamOpt(:,3)-gParamOpt(:,1) gParamOpt(:,4)-gParamOpt(:,2) ...
    gParamOpt(:,6)-gParamOpt(:,5) gParamOpt(:,8)-gParamOpt(:,7)];

gParamOpt = [gParamOpt gParamOpt(:,19)-gParamOpt(:,18) gParamOpt(:,21)-gParamOpt(:,20)];

iRemove = gParamOpt(:,10) < -3;
gParamOpt(iRemove,:) = [];
fitt(iRemove) = [];

iRemove = fitt > 200;
gParamOpt(iRemove,:) = [];
fitt(iRemove) = [];

iRemove = gParamOpt(:,11) > -4;
gParamOpt(iRemove,:) = [];
fitt(iRemove) = [];

% 
% iRemove = gParamOpt(:,7) > 1;
% gParamOpt(iRemove,:) = [];
% fitt(iRemove) = [];
% 
% 
% iRemove = gParamOpt(:,19) < gParamOpt(:,18);
% gParamOpt(iRemove,:) = [];
% fitt(iRemove) = [];




for varID = 1:23
    subplot(6,4,varID);

    plot(gParamOpt(:,varID),fitt,'.r');
    axis([min(gParamOpt(:,varID)) max(gParamOpt(:,varID)) 68 100]);
    title(names(varID)); 
end

ggParamOpt = gParamOpt(fitt < 80,:);
gFit = fitt(fitt < 80);


[LL, AIC, AICc, BIC, pValue1, pValue2] = pwGoodnessOfFit(min(fitt), 61, 15);

csvwrite('processed.csv',[ggParamOpt gFit']);
