
clc; clear;

names = {'B1','B2','U1','U2','xFwd1','xRev1','xFwd3','xRev3',...
    'AXLint1','AXLint2','scaleA','kRec','kDeg','fElse','fD2',...
    'Gas1','AXL2'};

IDX = 1;

for ii = 1:17
    subplot(4,5,ii);
    
    filess = dir(['*-' mat2str(ii) '.mat']);
    
    
    for jj = 1:length(filess)
        load(filess(jj).name,'fitIDXglobal');
        load(filess(jj).name,'vv');
        load(filess(jj).name,'paramOpt');
        
        if (length(fitIDXglobal) < 25)
            fitIDXglobal(25) = 0;
        end
        
        fitIDXglobal(fitIDXglobal == 10E5) = 10E5;
        fitIDXglobal(fitIDXglobal == 0) = 10E5;
        
        
        for xx = 1:length(fitIDXglobal)
            if fitIDXglobal(xx) < 80
                gParamOpt(IDX,:) = paramOpt{xx};
                IDX = IDX + 1;
            end
        end
        
        if (jj == 1)
            vBest = fitIDXglobal;
        else
            vBest = min(vBest,fitIDXglobal);
        end
        
        
        plot(vv,vBest);
        axis([min(vv) max(vv) 70 120]);
        title(names(ii));
    end
end

gParamOpt(gParamOpt(:,10) < -4,:) = [];