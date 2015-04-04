clc; clear;

[fitt, ~, names, params, cutoff, ] = loadParam;

xxx = logspace(-2.0,2.0,20);

stimProfile = @(x, pp) cLib_profile (240, x, 10);

outter = zeros(size(params,1),9,length(xxx));

IDXX = [1:7 9:10];

for ii = 1:size(params,1)
    [B, BT] = stimProfile(params(ii,:));
    B = B / BT;
    
    for jj = 1:length(IDXX)
        for xx = 1:length(xxx)
            params2 = params(ii,:);
            params2(IDXX(jj)) = params2(IDXX(jj))*xxx(xx);
            
            [C, CT] = stimProfile(params2);
            C = C / CT;

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
