clc; clear;

load 77j87EXi

names ={'B2', 'U2', 'xFwd1','xRev1','xFwd3','xRev3', 'AXLint1','AXLint2',...
    'kRec','kDeg','fElse','AXL2','Gas1','picker'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

params(fitt > min(fitt)+3,:) = [];
fitt(fitt > min(fitt)+3) = [];

%%
clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

% for ii = 1:13
%     subplot(3,5,ii);
%     
%     plot(params(:,ii), fitt,'.');
%     title(names(ii));
%     axis([minn(ii) maxx(ii) min(fitt) max(fitt)]);
% end

params = [ones(size(params,1),1)*0.06, 10.^params];
params(:,end) = params(:,end) > 0.05;

%%

xxx = logspace(-2.0,2.0,100);

stimProfile = @(x, pp) cLib_profile ([0 240], x, 10, pp);

outter = zeros(size(params,1),9,length(xxx));

for ii = 1:size(params,1)
    
    for jj = 1:1
        for xx = 1:length(xxx)
            params2 = params(ii,:);
            params2(jj) = params2(jj)*xxx(xx);
            params2(jj+2) = params2(jj+2)*xxx(xx);
            params2(jj+4) = params2(jj+4)*xxx(xx);
            

            B = stimProfile(params(ii,:), 1);
            C = stimProfile(params2, 1);
            %D = stimProfile(params2, 3);
            

            if ~isempty(C) && ~isempty(B)
                %B = B(2) / B(1);
                %C = C(2) / C(1);
                %D = D(2) / D(1);
                
                outter(ii,jj,xx) = C(2) / B(2);
            else
                outter(ii,jj,xx) = NaN;
            end
        end
        
        
    end 
    
    for jj = 1:1
        %subplot(3,3,jj);
        loglog(xxx,squeeze(outter(:,jj,:))');
        title(names(jj));
        
    end
    drawnow;
end



