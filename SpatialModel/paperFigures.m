clc; clear;

params = 10.^[-2.3829 0.77811 -1.704 3.2353 -5 -1.1466 -0.5604 -0.024512 -2.2774 -3 ...
	-0.6862 -2.5648 -1 -1.2476 -3.923E-06 -2.1966 2.899];

figIDX = 1;

sshape = @(x,xx) cos(xx*pi/3).^x(1);
GasConc = -1.204;
GassF = @(x,xx) 2*sshape(x,xx)/mean(sshape(x,xx).*xx)*(10^x(2));

shapeParam = 60;

names = {'B1','B2','U1','U2','xFwd1','xRev1','xFwd3','xRev3','AXLint1','AXLint2',...
	'scaleA','kRec','kDeg','fElse','fD2','Gas1','AXL2'};

D = zeros(1,14);
D(2) = 1;

DO_ALL_OUTPUTS = 0;
DO_IMPAIRED_ENDO = 0;
DO_LONG_TIME_SCALE = 0;
SHOW_FIT = 0;
DO_SPATIAL_PRED = 0;
LOCAL_SENSITIVITY = 0;
SPATIAL_TIMECOURSE = 1;

%% Do spatial time course


if SPATIAL_TIMECOURSE
	tps = linspace(0,30,100);
	xx = linspace(0,1,80);
	
	GasConc = -1.204;
	
    subplot(1,2,1);
	
	pY = cLib_diff_profile_pYavg (tps, params, params(end), params(end-1), GassF([3 GasConc],xx), D*0.1, 1, 1, 0);
	plot(tps,pY / pY(1),'r');
	hold on;
	
	pY = cLib_diff_profile_pYavg (tps, params, params(end), params(end-1), GassF([0 -2.1966],xx), D*10, 1, 1, 0);
	plot(tps,pY / pY(1),'g');
	
	pY = cLib_profile (tps, params, params(end), params(end-1), 10^GasConc, 0);
	plot(tps,pY / pY(1),'k');
    
    
    
    subplot(1,2,2);
    
	pY = cLib_profile (tps, params, params(end), params(end-1), 10^GasConc, 2);
	plot(tps,pY / pY(1),'k');
    
end


%% Do local sensitivity analysis

if LOCAL_SENSITIVITY

	T = 30;

	figure(figIDX);
	figIDX = figIDX + 1;

	for ii = 1:15
		params2 = params;
		params2(ii) = params(ii)/10;
		outter(ii,1) = cLib_profile (T, params2, params2(end), params2(end-1), 1.25, 5);
		outterT(ii,1) = cLib_profile (T, params2, params2(end), params2(end-1), 1.25, 2);
	
	
		params2(ii) = params(ii);
		outter(ii,2) = cLib_profile (T, params2, params2(end), params2(end-1), 1.25, 5);
		outterT(ii,2) = cLib_profile (T, params2, params2(end), params2(end-1), 1.25, 2);
	

		params2(ii) = params(ii)*10;
		outter(ii,3) = cLib_profile (T, params2, params2(end), params2(end-1), 1.25, 5);
		outterT(ii,3) = cLib_profile (T, params2, params2(end), params2(end-1), 1.25, 2);
	end

	bar(outter);
	%figure;
	%bar(outterT);
end

%% Figure 2A - pY time course

if SHOW_FIT
    figure(figIDX);
    figIDX = figIDX + 1;
    
    tps = linspace(0,10,100);
    Gass = [1.25 0.25 0.05 0.01];

    for ii = 1:length(Gass)
        outter(ii,1:length(tps)) = cLib_profile (tps, params, params(end), params(end-1), Gass(ii), 1);
        outter(ii,:) = outter(ii,:) / outter(ii,1);
        outterS(ii,1:length(tps)) = cLib_profile (tps, params, params(end), params(end-1), Gass(ii), 4);
        outterS(ii,:) = outterS(ii,:) / outterS(ii,1);
    end

    subplot(2,2,1);
    plot(tps,outter');
    axis([min(tps) max(tps) 0 3]);

    subplot(2,2,2);
    plot(tps,outterS');
    axis([min(tps) max(tps) 0 1.2]);

    % Now concentration

    subplot(2,2,3);

    Gass = [0 logspace(log10(0.156/2),log10(2.5),50)];

    for ii = 1:length(Gass)
        pY(ii) = cLib_profile (240, params, params(end), params(end-1), Gass(ii), 0);
        ratioo(ii) = cLib_profile (240, params, params(end), params(end-1), Gass(ii), 1);
        totall(ii) = cLib_profile (240, params, params(end), params(end-1), Gass(ii), 2);
        surf(ii) = cLib_profile (240, params, params(end), params(end-1), Gass(ii), 4);
    end

    Gass(1) = Gass(2)/2;

    semilogx(Gass,pY / pY(1));
    hold on;
    semilogx(Gass,ratioo / ratioo(1));
    semilogx(Gass,totall / totall(1));
    axis([min(Gass),max(Gass),0,5]);

    subplot(2,2,4);

    semilogx(Gass,surf / surf(1));
    axis([min(Gass),max(Gass),0,1.2]);

    clear Gass pY ratioo totall surf tps outter outterS;
end

%% Figure Supp - All possible outputs

if DO_ALL_OUTPUTS

    figure(figIDX);
    figIDX = figIDX + 1;

    TPS = linspace(0,10,1000);

    detail = 5;

    Gass = logspace(-2,0.0969,detail);
    AXXL = logspace(-1,3,detail);
    autoC = logspace(-4,0,detail);

    for ii = 1:length(AXXL)
        for jj = 1:length(Gass)    
            for xx = 1:length(autoC)
                if ii == 1 && jj == 1 && xx == 1
                    YY(1,:) = cLib_profile (TPS, params, AXXL(ii), autoC(xx), Gass(jj), 0);
                    condition(1,:) = [AXXL(ii) autoC(xx) Gass(jj)];
                else
                    YY(end+1,:) = cLib_profile (TPS, params, AXXL(ii), autoC(xx), Gass(jj), 0);
                    condition(end+1,:) = [AXXL(ii) autoC(xx) Gass(jj)];
                end

                YY(end,:) = YY(end,:) / YY(end,1);
            end
        end
    end

    plot(TPS,YY')
    
    clear Gass AXXL autoC YY TPS detail;
end


%% Long time scale

if DO_LONG_TIME_SCALE
    figure(figIDX);
    figIDX = figIDX + 1;

    clear YY;

    TPS = 240;

    Gass = logspace(-2,1,50);
    AXXL = logspace(-1,3,4);

    color = 'kbrg';

    for ii = 1:length(AXXL)
        for jj = 1:length(Gass)
            YY(jj) = cLib_profile (240, params, AXXL(ii), params(end-1), Gass(jj), 1);
        end

        YY = YY / YY(1);

        semilogx(Gass, YY',color(ii));
        hold on;
    end
end

%% Impaired endocytosis prediction

if DO_IMPAIRED_ENDO
	
	h = waitbar(0,'Starting impaired endocytosis prediction...'); %#ok<*UNRCH>

    figure(figIDX);
    figIDX = figIDX + 1;

    xx = linspace(0,1,10);

    A = logspace(-2,0,100);

    for ii = 1:length(A)
        B(ii) = cLib_diff_profile_pYavg (30, params, params(end), params(end-1), GassF([-100 GasConc shapeParam],xx), zeros(1,14), A(ii), 1, 0); %#ok<SAGROW>
        
        waitbar(ii/length(A),h,'Endocytosis prediction...')
    end
    
    close(h)
    
    semilogx(A,B/B(end),'r-');
end

%% Real Gas6 local prediction

% This script generated Figure 2B

if DO_SPATIAL_PRED
    figure(figIDX);
    figIDX = figIDX + 1;
    
    xx = linspace(0,1,100);
    A = logspace(-2,3,10);

    for ii = 1:length(A)
        BB(:,ii) = GassF([A(ii) GasConc shapeParam],xx); %#ok<SAGROW>
    end

    subplot(2,2,1);
    plot(BB);
    title('Gas6 profile');
    
    % Next

    ffSpat = @(x,y) cLib_diff_profile_pY (30, params, params(end), params(end-1), GassF(x,xx), D*y, 1, 1, 0);

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
    ttt = [0 1 0 0 0 0 0 1 0 0 0 0 0 0;...
           0 0 1 0 0 0 0 0 1 0 0 0 0 0;...
           0 0 0 1 0 0 0 0 0 1 0 0 0 0;...
           0 0 0 0 1 0 0 0 0 0 1 0 0 0;...
           0 0 0 0 0 1 0 0 0 0 0 1 0 0;...
           0 0 0 0 0 0 1 0 0 0 0 0 1 0];

    meanify = @(x) 2*x'*xx'/length(xx);

    ffSpat = @(x,y) cLib_diff_profile (30, params, params(end), params(end-1), GassF(x,xx), D*10, 1, 1);

    clc;
    
    for ii = 1:length(A)
        B(ii,:) = ttt*meanify(ffSpat([A(ii) GasConc],1)); %#ok<SAGROW>
        disp(ii);
    end
    
    for ii = 1:size(B,2)
        B(:,ii) = B(:,ii) - B(1,ii); %#ok<SAGROW>
    end

    semilogx(A,B);
    axis([min(A) max(A) min(min(B)) max(max(B))]);
    legend(names);
end