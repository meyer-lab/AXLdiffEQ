function interpretVar()
clc; clear;
    params2 = U87SsensFullWider;
    
    paramGrid (params2);
end

function paramGrid (params2)
    
    subplot(1,3,1);
    tps = linspace(0,10,100);
    
    % % % % The real data
    pyK = [0.911, 0.904, 0.935, 1.08, 1.25];
    pyKE = [0.134, 0.0991, 0.152, 0.135, 0.130];
    kTps = [0, 0.5, 1, 5, 10];
    pY = [1.00, 1.93, 1.35, 1.74, 0.89, 0.00; 1.00, 1.33, 1.05, 1.10, 1.61, 0.00]';
    pYE = [0.19, 0.17, 0.41, 0.30, 0.44, 0.00; 0.19, 0.15, 0.12, 0.28, 0.20, 0.00]';
    
    GasConc = [64, 16, 4, 1, 0.25, 0.1];
    GasConc = fliplr(GasConc);
    
    tot = [1.00, 0.75, 0.72, 0.87, 0.83, 0.81; 0.83, 0.85, 0.76, 0.68, 0.75, 0.75]';
    totE = [0.12, 0.05, 0.14, 0.05, 0.04, 0.08; 0.05, 0.05, 0.09, 0.11, 0.03, 0.04]';
    % % % End of real data
    

    for ii = 1:size(params2,1)
            [stimProfile2, ~, ~, convF] = cLib_profile (tps, params2(ii,:), 1.25);

            plot(tps,stimProfile2*convF(3)/9,'k');
            axis([min(tps) max(tps) 0 2]);
            hold on;
    end
    
    errorbar(kTps, pyK, pyKE,'k','LineWidth',2);
    
    ylabel('AXL pan-pY');
    xlabel('Time (min)');
    fixFig(gca);
    
    Gass = [0 logspace(log10(0.25),log10(64),10)];

    for ii = 1:size(params2,1)
        for jj = 1:length(Gass)
                [stimProfile, stimProfileTot, ~, convF] = ...
                    cLib_profile ([0 60 240], params2(ii,:), Gass(jj));
              
                outFig{1}(ii,jj) = stimProfile(2)*convF(1);
                outFig{2}(ii,jj) = stimProfile(3)*convF(1);

                outFig{3}(ii,jj) = stimProfileTot(2); 
                outFig{4}(ii,jj) = stimProfileTot(3);
        end
    end
    
    median(outFig{3}(:,1)./3000)
    
    Gass(1) = 0.1;

    ylab = {'AXL pan-pY \div Total (RU)', 'Total AXL (RU)'};
    conv = [8, 3000, 0.25];
    
    for ii = 1:size(outFig{4},1)
        subplot(1,3,2);
        semilogx(Gass,outFig{1}(ii,:)' / conv(1),'b');
        hold on;
        semilogx(Gass,outFig{2}(ii,:)' / conv(1),'k');
        axis([min(Gass) max(Gass)+0.01 0 1.5]);
        fixFig(gca);
        set(gca,'Xtick',Gass);
        set(gca,'XtickLabel',[0.0 0.25 1.0 4.0 16 64]);
        ylabel(ylab{1});
        xlabel('Exogenous {\Delta}Gla Gas6 [nM]');
    end
    
    for ii = 1:size(outFig{4},1)
        subplot(1,3,3);
        semilogx(Gass,outFig{3}(ii,:)' / mean(outFig{3}(ii,:)) / 1.2,'b');
        hold on;
        semilogx(Gass,outFig{4}(ii,:)' / mean(outFig{4}(ii,:)) / 1.2,'k');
        axis([min(Gass) max(Gass)+0.01 0 1.5]);
        fixFig(gca);
        set(gca,'Xtick',[0.1 0.25 1.0 4.0 16 64]);
        set(gca,'XtickLabel',[0.0 0.25 1.0 4.0 16 64]);
        ylabel(ylab{2});
        xlabel('Exogenous {\Delta}Gla Gas6 [nM]');
    end
    
    subplot(1,3,2);
    Y = pY ./ tot / 2;
    YE = sqrt((pYE ./ pY).^2 + (totE./tot).^2) .* Y;
    
    errorbar(GasConc(1:5), Y(1:5,1), YE(1:5,1),'b','LineWidth',2);
    errorbar(GasConc(1:5), Y(1:5,2), YE(1:5,2),'k','LineWidth',2);
    fixFig(gca);

    subplot(1,3,3);
    errorbar(GasConc, tot(:,1), totE(:,1),'b','LineWidth',2);
    errorbar(GasConc, tot(:,2), totE(:,2),'k','LineWidth',2);
    fixFig(gca);
    
    set(gcf, 'Position', [100 100 1200 1000])
    
    export_fig('U87','-pdf');
end

function fixFig (gcc)
    set(gcc,'FontName','Helvetica Neue');
    set(gcc,'FontSize',12);
    axis(gcc,'square');
end

function [fitt, params, names, params2, cutoff, minn, maxx, minD2] = loadParam () %#ok<STOUT>

    load fitCombined;

    names = {'U2','xFwd1','xRev4','int1','int2','kRec','kDeg','fElse','AXL','Gas','pD'};

    fitt = [];
    params = [];

    for ii = 1:length(fitStruct) %#ok<USENS>
        fitt = [fitt, fitStruct{ii}.fitIDXglobal]; %#ok<AGROW>
        params = [params; fitStruct{ii}.paramOpt]; %#ok<AGROW>
    end

    minD2(1) = min(fitt(params(:,end) > 0.5));
    minD2(2) = min(fitt(params(:,end) < 0.5));

    cutoff = min(fitt)+3;

    params(fitt > cutoff,:) = [];
    fitt(fitt > cutoff) = [];

    [fitt, IDX] = sort(fitt);
    params = params(IDX,:);
    params = flipud(params);
    fitt = fliplr(fitt);

    params2 = 10.^params;
    params2(:,end) = params2(:,end) > 0.05;

end

% Performs the fitting procedure.
function params2 = U87SsensFullWider ()
    [~, ~, ~, params2] = loadParam ();

    Dopts = psoptimset('TimeLimit',45*60,'Display','off',...
        'CompletePoll','on','CompleteSearch','on','Vectorized','on');

    parfor ii = 1:size(params2,1)
        [paramOpt(ii,:)] = ...
            patternsearch(@(x) U87cLibLoc(x, params2(ii,:), 2),...
            [3 -2],[],[],[],[],[1 -4],[5 1],[],Dopts);
    end
    
    params2(:,9:10) = 10.^paramOpt;
end

% Calls the C routine to calculate error from parameter sets.
function outter = U87cLibLoc (in, pFit, N)
    in2 = ones(numel(in)/N,1)*pFit;
    
    in2(:,9:10) = 10.^in;

    if ~libisloaded('libOptimize')
        loadlibrary('libOptimize.dylib');
    end

    outP = libpointer('doublePtr',zeros(1,numel(in)/N));

    calllib('libOptimize','U87pyEntryVec',libpointer('doublePtr',in2'),outP,numel(in)/N);

    outter = outP.Value';
end

% This calculates the response to stimulation for the uniform case
function [pY, tot, surf, convF, species] = cLib_profile (tps, params, GasStim)

    if ~libisloaded('libOptimize')
        loadlibrary('libOptimize.dylib')
    end

    out1 = libpointer('doublePtr',1:length(tps));
    out2 = libpointer('doublePtr',1:length(tps));
    out3 = libpointer('doublePtr',1:length(tps));
    out4 = libpointer('doublePtr',1:(length(tps)*13));
    convFac = libpointer('doublePtr',1:3);

    x = calllib('libOptimize','calcProfileMatlab',out1,out2,out3,out4,...
        libpointer('doublePtr',params),libpointer('doublePtr',tps),uint32(length(tps)),GasStim,convFac);

    if x == 0
        pY = out1.Value;
        tot = out2.Value;
        surf = out3.Value;
        convF = convFac.Value;
        species = out4.Value;
    else
        pY = [];
        tot = [];
        surf = [];
        convF = [];
        species = [];
    end

end