clc; clear;

namesF = {'new_v4d_noBT','new_v4e_noBT','new_cochran_noBT','new_cochran2_noBT','new_cochran2_noAsig_noBT'};
p = [0 0 0 0 0 0 0 0 0 0 0];
N = 104; % 

for ii = 1:length(namesF)
    load(namesF{ii});
    
    chiSqs{ii} = chiSq;
    
    paramss{ii} = params;
    p(ii) = p(ii) + size(params,2);

    [chis(ii), aIDX(ii)] = min(chiSq);  %#ok<*SAGROW>
    
    
    BestPars{ii} = params(aIDX(ii),:);

    [LL(ii), AIC(ii), AICc(ii), BIC(ii), pValue1(ii), pValue2(ii)] = pwGoodnessOfFit(chis(ii), N, p(ii));
end

