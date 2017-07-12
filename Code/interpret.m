function interpret()
clc; clear;
    [~, params, names, params2, ~, minn, maxx, minD2] = loadParam;
    
    figure(1);
    %paramGrid (params, names, minn, maxx, params2, minD2);
    %figure(3)
    varParams (params2, names);
end

function paramGrid (params, names, minn, maxx, params2, minD2)

names

    names{1} = 'k_{rb2} [Log_{10}(min^{-1})]';
    names{2} = 'k_{f} [Log_{10}(cell{\cdot}receptors^{-1}{\cdot}min^{-1})]';
    names{3} = 'k_{r4} [Log_{10}(min^{-1})]';
    names{4} = 'k_{int,1} [Log_{10}(min^{-1})]';
    names{5} = 'k_{int,2} [Log_{10}(min^{-1})]';
    names{6} = 'k_{rec} [Log_{10}(min^{-1})]';
    names{7} = 'k_{deg} [Log_{10}(min^{-1})]';
    names{8} = 'f_{else}';
    names{9} = 'AXL Expression [Log_{10}(min^{-1})]';
    names{10} = 'Autocrine Gas6 [Log_{10}(nM)]';

    plotD = @(idx1, idx2) plotData(names, minn, maxx, params, idx1, idx2);
    
    subplot(4,4,3);
    plotD(3, 2);
    subplot(4,4,4);
    plotD(1, 4);
    subplot(4,4,7);
    plotD(5, 6);
    subplot(4,4,8);
    plotD(7, 8);
    subplot(4,4,11);
    plotD(9, 10);
    
    
    
    subplot(4,4,12);
    ax = bar(minD2,0.4,'k');
    set(gca,'XtickLabel',{'Active','Inactive'});
    ylabel('Sum of Squared Error');
    fixFig(gca);
    
    
    
    
    subplot(4,4,13);
    tps = linspace(0,10,100);
    
    
    
    % % % % The real data
    pyK = [4.1, 3.6, 5.3, 11.3, 11.6];
    pyKE = [1.0, 1.3, 1.1, 0.7, 0.4];
    kTps = [0, 0.5, 1, 5, 10];
    
    pY = [10.75952427, 8.305264139; 7.390399159, 7.056438019; 7.144036441, 7.680851079;
        4.570826833, 8.184089069; 6.107714557, 7.204021903; 7.535575387, 7.535575387];
    
    pYE = [1.74431775, 2.100723242; 1.267611, 1.260108508; 0.898008437, 1.680415875;
        1.521844479, 0.871927763; 0.932623012, 0.563182873; 0.812417951, 0.812417951];
    
    GasConc = [64, 16, 4, 1, 0.25, 0.1];
    
    
    tot = [3443.11, 3219.69; 3143.41, 3353.82; 3018.88, 3611.82; 2608.88, 3448.21;
        2690.24, 3168.14; 2672.00, 2672.00];
    
    totE = [174.38, 132.10; 189.03, 129.93; 245.75, 225.42; 154.89, 203.72; 128.72, 187.34; 82.62, 82.62];
    
    surf = [0.206, 0.239; 0.274, 0.316; 0.281, 0.251; 0.220, 0.302; 0.256, 0.281; 0.257, 0.337];
    surfE = [0.043, 0.015; 0.047, 0.037; 0.032, 0.030; 0.025, 0.036; 0.044, 0.035; 0.030, 0.023];
    
    % % % End of real data
    

    for ii = 1:10%size(params,1)
            [stimProfile2, ~, ~, convF] = cLib_profile (tps, params2(ii,:), 1.25);

            plot(tps,stimProfile2*convF(3),'k');
            axis([min(tps) max(tps) 0 14]);
            hold on;
            fixFig(gca);
    end
    
    errorbar(kTps, pyK, pyKE,'k','LineWidth',2);
    
    ylabel('AXL pan-pY');
    xlabel('Time (min)');
    
    Gass = [0 logspace(log10(0.25),log10(64),10)];

    for ii = 1:10%size(params2,1)
        for jj = 1:length(Gass)
                [stimProfile, stimProfileTot, stimProfileSurf, convF] = ...
                    cLib_profile ([0 60 240], params2(ii,:), Gass(jj));

                outFig{1}(ii,jj) = stimProfile(2)*convF(1); %#ok<AGROW>
                outFig{2}(ii,jj) = stimProfile(3)*convF(1); 

                outFig{3}(ii,jj) = stimProfileTot(2); 
                outFig{4}(ii,jj) = stimProfileTot(3); 

                outFig{5}(ii,jj) = stimProfileSurf(2)*convF(2); 
                outFig{6}(ii,jj) = stimProfileSurf(3)*convF(2); 
        end
    end

    Gass(1) = 0.1;

    
    ylab = {'AXL pan-pY \div Total (RU)', 'Total AXL (RU)', 'Surface AXL (RU)'};
    conv = [8, 3000, 0.25];
    
    for ii = 1:10%size(params,1)
        for xx = 1:2:5
            subplot(4,4,(xx+1)/2+13);
            semilogx(Gass,outFig{xx}(ii,:)' / conv((xx+1)/2),'b');
            hold on;
            semilogx(Gass,outFig{xx+1}(ii,:)' / conv((xx+1)/2),'k');
            axis([min(Gass) max(Gass)+0.01 0 1.5]);
            fixFig(gca);
            set(gca,'Xtick',[0.1 0.25 1 4 16 64]);
            set(gca,'XtickLabel',[0.0 0.25 1 4 16 64]);
            ylabel(ylab{(xx+1)/2});
            xlabel('Exogenous {\Delta}Gla Gas6 [nM]');
        end
    end
    
    
    subplot(4,4,14);
    Y = pY(:,1) ./ tot(:,1) * 3000 / conv(1);
    
    YE = sqrt((pYE(:,1) ./ pY(:,1)).^2 + (totE(:,1)./tot(:,1)).^2) .* Y;
    errorbar(GasConc, Y, YE,'b','LineWidth',2);
    
    Y = pY(:,2) ./ tot(:,2) * 3000 / conv(1);
    Y
    YE = sqrt((pYE(:,2) ./ pY(:,2)).^2 + (totE(:,2)./tot(:,2)).^2) .* Y;
    errorbar(GasConc, Y, YE,'k','LineWidth',2);
    
    
    
    subplot(4,4,15);
    errorbar(GasConc, tot(:,1)/conv(2), totE(:,1)/conv(2),'b','LineWidth',2);
    errorbar(GasConc, tot(:,2)/conv(2), totE(:,2)/conv(2),'k','LineWidth',2);
    
    
    subplot(4,4,16);
    errorbar(GasConc, surf(:,1)/conv(3), surfE(:,1)/conv(3),'b','LineWidth',2);
    errorbar(GasConc, surf(:,2)/conv(3), surfE(:,2)/conv(3),'k','LineWidth',2);
    
    
    
    set(gcf, 'Position', [100 100 1200 1000])
    
    %export_fig('fitParams','-pdf');
end

function plotData (names, minn, maxx, params, idx1, idx2)
    plot(params(:,idx1),params(:,idx2),'k.','MarkerSize',10);
    xlabel(names(idx1));
    ylabel(names(idx2));
    axis([minn(idx1) maxx(idx1) minn(idx2) maxx(idx2)]);
    fixFig(gca);
end

% This varies each parameter (except fElse) to examine the effect on
% phosphorylation
function varParams (params, names)
    xxx = logspace(-3.0,3.0,20);
    
    names{1} = 'k_{rb2}';
    names{2} = 'k_{f1}';
    names{3} = 'k_{r4}';

    stimProfile = @(x, pp) cLib_profile (30, x, 10);

    outter = zeros(size(params,1),9,length(xxx));

    IDXX = 1:3;

    for ii = 1:size(params,1)
        [B, BT] = stimProfile(params(ii,:));
        B = B / BT;

        for jj = 1:length(IDXX)
            parfor xx = 1:length(xxx)
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
            loglog(xxx,squeeze(outter(ii,jj,:))','k');
            hold on;
            loglog(xxx,ones(size(xxx)),'Color',[0.5 0.5 0.5]);
            title(names(IDXX(jj)));
            ylabel('Relative AXL pY \div Total');
            xlabel('Log_{10} Fold Change in Parameter');
            set(gca,'Xtick',10.^(-3:3:3))
            axis([min(xxx) max(xxx) 0.1 100]);
            fixFig(gca);
        end
        drawnow;
    end
    
    subplot(3,4,4);
    for ii = 1:size(params,1)
        Kd = (params(ii,1) ./ 0.06)/1000;
        Kd1 = 0.035/1000;
        x = logspace(-4, 1,1000);

        semilogx(x*1E3, Kd ./ (Kd + x),'k')
        hold on;    
    end
    
    hlin = semilogx(x*1E3, Kd1 ./ (Kd1 + x),'r');
    set(hlin,'linewidth',2);
    
    axis([min(x)*1E3 max(x)*1E3 0 1]);
    set(gca,'Xtick',10.^(-1:4));
    xlabel('Gas6 Concentration [nM]');
    fixFig(gca);
    ylabel('Fraction Sites Unbound');
    set(gcf, 'Position', [600 600 1200 700])
    
    
    doSpatialPred (params)
    
    export_fig('varParams','-pdf');
end

function fixFig (gcc)
    set(gcc,'FontName','Helvetica Neue');
    set(gcc,'FontSize',12);
    axis(gcc,'square');
end

% Real Gas6 local prediction

% This script generated Figure 2B

function doSpatialPred (paramsIn)

sshape = @(x,xx) cos(xx*pi/3).^x(1);
GasConc = 32;
GassF = @(x,xx) 2*sshape(x,xx)/max(sshape(x,xx))*(x(2));

shapeParam = 60;

D = zeros(1,14);
D(1) = 1;
D(7:end) = 1;

xx = linspace(0,1,60);
A = logspace(-1,3,116);

inC = 0;
aa = [1 0 0 0 0 0 inC 0];

ttt = [aa 0 0 0 0 0;...
       0 aa 0 0 0 0;...
       0 0 aa 0 0 0;...
       0 0 0 aa 0 0;...
       0 0 0 0 aa 0;...
       0 0 0 0 0 aa];

meanify = @(x,xx) 2*x'*xx'/length(xx);

for ii = 1:length(A)
    BB(:,ii) = GassF([A(ii) GasConc shapeParam],xx);
        
    meanConc(ii) = meanify(BB(:,ii),xx);
    peakConc(ii) = max(BB(:,ii));
end

subplot(3,4,5);
hold on;
for ii = 1:10:length(A)
    plot(xx,BB(:,ii),'Color',([1 1 1] / A(ii) / 10).^(0.1) - 0.2);
end

title('Gas6 Profile');
ylabel('Gas6 Concentration [nM]');
xlabel('Radius (Arbitrary Units)');
fixFig(gca);


subplot(3,4,6);
loglog(A,meanConc,'r');
hold on;
loglog(A,peakConc,'k');
axis([min(A) max(A) 0.1 100]);
set(gca,'Xtick',10.^(-1:3));
title('Gas6 Profile Attributes');
xlabel('Spatial Inhomogeneity Parameter');
ylabel('Gas6 Concentration [nM]');
fixFig(gca);

A = A(1:5:end);

for yyy = 1:size(paramsIn,1)

    params = paramsIn(yyy,:);

    ffSpat = @(x,y) cLib_diff_profile (10, params, GassF(x,xx), D*y);

    localPY = zeros(length(A),4);
    avgPY = zeros(length(A),4);

    parfor ii = 1:length(A)
        localPYz = zeros(1,4);
        avgPYz = zeros(1,4);
        
        [~, temp] = ffSpat([A(ii) GasConc],0);
        localPYz(1) = temp(1);
        avgPYz(1) = mean(temp.*xx);
        
        [~, temp] = ffSpat([A(ii) GasConc],1);
        localPYz(2) = temp(1);
        avgPYz(2) = mean(temp.*xx);

        [~, temp] = ffSpat([A(ii) GasConc],10);
        localPYz(3) = temp(1);
        avgPYz(3) = mean(temp.*xx);
        
        [spec{ii}, temp] = ffSpat([A(ii) GasConc],100);
        localPYz(4) = temp(1);
        avgPYz(4) = mean(temp.*xx);
        
        localPY(ii,:) = localPYz;
        avgPY(ii,:) = avgPYz;
        
        disp(ii);
    end
    
    subplot(3,4,7);
    if yyy == 1
        semilogx(A,localPY(:,1)/localPY(1,1),'k-');
    end
    hold on;
    semilogx(A,localPY(:,1)/localPY(1,1),'k-');
    semilogx(A,localPY(:,2)/localPY(1,2),'g-');
    semilogx(A,localPY(:,3)/localPY(1,3),'b-');
    semilogx(A,localPY(:,4)/localPY(1,4),'r-');
    axis([min(A) max(A) 0 2000]);
    set(gca,'Xtick',10.^(-1:3));
    ylabel('Center (r = 0) AXL pY');
    xlabel('Spatial Inhomogeneity Parameter');
    fixFig(gca);

    subplot(3,4,8);
    if yyy == 1
        semilogx(A,avgPY(:,1)/avgPY(1,1),'k-');
    end
    hold on;
    semilogx(A,avgPY(:,1)/avgPY(1,1),'k-');
    semilogx(A,avgPY(:,2)/avgPY(1,2),'g-');
    semilogx(A,avgPY(:,3)/avgPY(1,3),'b-');
    semilogx(A,avgPY(:,4)/avgPY(1,4),'r-');
    axis([min(A) max(A) 0 3]);
    set(gca,'Xtick',10.^(-1:3));
    ylabel('Average AXL pY');
    xlabel('Spatial Inhomogeneity Parameter');
    fixFig(gca);
    
    % Components
    
    subplot(3,4,9);

    %names = {'A','A1','A2','A12','D1','D2'};
    
    for ii = 1:length(A)
        temp = spec{ii};
        Bdetail(:,ii) = meanify(temp,xx); %#ok<AGROW>
        B(ii,:) = ttt*Bdetail(:,ii); %#ok<AGROW>
    end
    
    for ii = 1:size(B,2)
        B(:,ii) = B(:,ii) / B(1,ii);
    end

    if yyy == 1
        h = loglog(A,B);
    end
    
    hold on;
    h = loglog(A,B);
    set(h(1),'Color','k');
    set(h(2),'Color','b');
    set(h(3),'Color','g');
    set(h(4),'Color','r');
    set(h(5),'Color','c');
    set(h(6),'Color','m');
    axis([min(A) max(A) min(min(B)) max(max(B))]);
    ylabel({'Normalized Mean', 'Surface Abundance'});
    %legend(names);
    fixFig(gca);
    drawnow;
    
end

end

function [outter, pY, tot, surf] = cLib_diff_profile (tps, params, GasIn, Din)
    Nspecies = 13;

    if ~libisloaded('libOptimize')
        loadlibrary('libOptimize.dylib')
    end

    dataPtr = libpointer('doublePtr',1:(length(tps)*length(GasIn)*Nspecies));
    dataPtrpY = libpointer('doublePtr',1:(length(tps)*length(GasIn)));
    dataPtrTot = libpointer('doublePtr',1:(length(tps)*length(GasIn)));
    dataPtrSurf = libpointer('doublePtr',1:(length(tps)*length(GasIn)));

    x = calllib('libOptimize','diffCalc',dataPtr, dataPtrpY, dataPtrTot, ...
        dataPtrSurf, libpointer('doublePtr',GasIn), length(GasIn), ...
        libpointer('doublePtr',params), ...
        libpointer('doublePtr',tps), length(tps), ...
        libpointer('doublePtr',Din));

    if x == 0
        outter = dataPtr.Value;
        pY = dataPtrpY.Value;
        tot = dataPtrTot.Value;
        surf = dataPtrSurf.Value;

        outter = reshape(outter,length(GasIn),[],length(tps));
    else
        error('Error!');
    end
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