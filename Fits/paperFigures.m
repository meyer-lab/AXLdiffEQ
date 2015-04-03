clc; clear;

params = 10.^[2.3626, -3.8869, 3.5788, -0.8938, -1.3400, -2.9947, ...
    -0.1982, -1.9999,  4.1196, -0.9247,  0.7935];

figIDX = 1;

sshape = @(x,xx) cos(xx*pi/3).^x(1);
GasConc = -1.204;
GassF = @(x,xx) 2*sshape(x,xx)/mean(sshape(x,xx).*xx)*(10^x(2));

shapeParam = 60;

names = {'B1','B2','U1','U2','xFwd1','xRev1','xFwd3','xRev3','AXLint1','AXLint2',...
	'scaleA','kRec','kDeg','fElse','fD2','Gas1','AXL2'};

D = zeros(1,14);
D(1) = 1;

DO_SPATIAL_PRED = 1;

%% Real Gas6 local prediction

% This script generated Figure 2B

if DO_SPATIAL_PRED
    figure(figIDX);
    figIDX = figIDX + 1;
    
    xx = linspace(0,1,120);
    A = logspace(-2,3,16);

    for ii = 1:length(A)
        BB(:,ii) = GassF([A(ii) GasConc shapeParam],xx); %#ok<SAGROW>
    end

    subplot(2,2,1);
    plot(BB);
    title('Gas6 profile');
    
    % Next

    ffSpat = @(x,y) cLib_diff_profile_pY (30, params, GassF(x,xx), D*y, 0);

    localPY = zeros(length(A),4);
    avgPY = zeros(length(A),4);

    parfor ii = 1:length(A)
        localPYz = zeros(1,4);
        avgPYz = zeros(1,4);
        
        temp = ffSpat([A(ii) GasConc],0);
        localPYz(1) = temp(1);
        avgPYz(1) = mean(temp.*xx);
        
        temp = ffSpat([A(ii) GasConc],0.1);
        localPYz(2) = temp(1);
        avgPYz(2) = mean(temp.*xx);

        temp = ffSpat([A(ii) GasConc],1);
        localPYz(3) = temp(1);
        avgPYz(3) = mean(temp.*xx);
        
        temp = ffSpat([A(ii) GasConc],10);
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
    axis([min(A) max(A) 0 2.2]);
    title('Average predictions');
    
    
    % Components
    
    subplot(2,2,4);

    names = {'A','A1','A2','A12','D1','D2'};
    ttt = [1 0 0 0 0 0 1 0 0 0 0 0 0;...
           0 1 0 0 0 0 0 1 0 0 0 0 0;...
           0 0 1 0 0 0 0 0 1 0 0 0 0;...
           0 0 0 1 0 0 0 0 0 1 0 0 0;...
           0 0 0 0 1 0 0 0 0 0 1 0 0;...
           0 0 0 0 0 1 0 0 0 0 0 1 0];

    meanify = @(x) 2*x'*xx'/length(xx);

    ffSpat = @(x) cLib_diff_profile (30, params, GassF(x,xx), D*10);

    clc;
    
    parfor ii = 1:length(A)
        Bdetail(:,ii) = meanify(ffSpat([A(ii) GasConc]));
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