function paperFigures() 

[~, ~, ~, params] = loadParam;
clc;

params = params(end,:);

sshape = @(x,xx) cos(xx*pi/3).^x(1);
GasConc = -1.204;
GassF = @(x,xx) 2*sshape(x,xx)/mean(sshape(x,xx).*xx)*(10^x(2));

shapeParam = 60;

D = zeros(1,14);
D(1) = 1;
D(7:end) = 1;

DO_SPATIAL_PRED = 1;

% Real Gas6 local prediction

% This script generated Figure 2B

if DO_SPATIAL_PRED
    xx = linspace(0,1,40);
    A = logspace(-2,3,20);

    for ii = 1:length(A)
        BB(:,ii) = GassF([A(ii) GasConc shapeParam],xx);
    end

    subplot(2,2,1);
    plot(BB);
    title('Gas6 profile');
    
    % Next

    ffSpat = @(x,y) cLib_diff_profile (10, params, GassF(x,xx), D*y);

    localPY = zeros(length(A),4);
    avgPY = zeros(length(A),4);

    for ii = 1:length(A)
        localPYz = zeros(1,4);
        avgPYz = zeros(1,4);
        
        [~, temp] = ffSpat([A(ii) GasConc],0);
        localPYz(1) = temp(1);
        avgPYz(1) = mean(temp.*xx);
        
        [~, temp] = ffSpat([A(ii) GasConc],0.1);
        localPYz(2) = temp(1);
        avgPYz(2) = mean(temp.*xx);

        [~, temp] = ffSpat([A(ii) GasConc],1);
        localPYz(3) = temp(1);
        avgPYz(3) = mean(temp.*xx);
        
        [~, temp] = ffSpat([A(ii) GasConc],10);
        localPYz(4) = temp(1);
        avgPYz(4) = mean(temp.*xx);
        
        localPY(ii,:) = localPYz;
        avgPY(ii,:) = avgPYz;
        
        disp(ii);
    end
    
    subplot(2,2,2);
    hold off;
    semilogx(A,localPY(:,1)/localPY(1,1),'k-');
    hold on;
    semilogx(A,localPY(:,2)/localPY(1,2),'g-');
    semilogx(A,localPY(:,3)/localPY(1,3),'b-');
    semilogx(A,localPY(:,4)/localPY(1,4),'r-');
    axis([min(A) max(A) 0 max(max(localPY / localPY(1,1)))]);
    title('Local predictions');

    subplot(2,2,3);
    hold off;
    semilogx(A,avgPY(:,1)/avgPY(1,1),'k-');
    hold on;
    semilogx(A,avgPY(:,2)/avgPY(1,2),'g-');
    semilogx(A,avgPY(:,3)/avgPY(1,3),'b-');
    semilogx(A,avgPY(:,4)/avgPY(1,4),'r-');
    axis([min(A) max(A) 0 4]);
    title('Average predictions');
    
    
    % Components
    
    subplot(2,2,4);

    names = {'A','A1','A2','A12','D1','D2'};
    ttt = [1 0 0 0 0 0 0.5 0 0 0 0 0 0;...
           0 1 0 0 0 0 0 0.5 0 0 0 0 0;...
           0 0 1 0 0 0 0 0 0.5 0 0 0 0;...
           0 0 0 1 0 0 0 0 0 0.5 0 0 0;...
           0 0 0 0 1 0 0 0 0 0 0.5 0 0;...
           0 0 0 0 0 1 0 0 0 0 0 0.5 0];

    meanify = @(x) 2*x'*xx'/length(xx);
    
    for ii = 1:length(A)
        temp = ffSpat([A(ii) GasConc], 10);
        Bdetail(:,ii) = meanify(temp);
        B(ii,:) = ttt*Bdetail(:,ii);
        
        disp(ii);
    end
    
    %for ii = 1:size(B,2)
    %    B(:,ii) = B(:,ii) - B(1,ii); %#ok<SAGROW>
    %end

    semilogx(A,B);
    axis([min(A) max(A) min(min(B)) max(max(B))]);
    legend(names);
end

end

function [outter, pY, tot, surf] = cLib_diff_profile (tps, params, GasIn, Din)
    Nspecies = 13;

    if ~libisloaded('libOptimizeDiff')
        loadlibrary('libOptimizeDiff.dylib','BlasHeader.h')
    end

    dataPtr = libpointer('doublePtr',1:(length(tps)*length(GasIn)*Nspecies));
    dataPtrpY = libpointer('doublePtr',1:(length(tps)*length(GasIn)));
    dataPtrTot = libpointer('doublePtr',1:(length(tps)*length(GasIn)));
    dataPtrSurf = libpointer('doublePtr',1:(length(tps)*length(GasIn)));

    x = calllib('libOptimizeDiff','diffCalc',dataPtr, dataPtrpY, dataPtrTot, ...
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


