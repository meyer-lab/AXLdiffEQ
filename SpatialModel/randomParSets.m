clc; clear;
rng('shuffle');

minn = log10([0.0006,0.0006,1E-5,1E-5,... % 'B1','B2','U1','U2'
    1E-5,1E-5,1E-5,1E-5,... % 'xFwd1','xRev1','xFwd3','xRev3'
    1E-6,1E-3,1E-24,... % 'AXLint1','AXLint2','scaleA'
    1E-3,1E-4,1E-5,1,1E-3,500]); % 'kRec','kDeg','fElse','fD2','Gas1','AXL2'

maxx = log10([6,6,1E5,1E5,... % 'B1','B2','U1','U2'
    1E3,1E5,1E3,1E5,... % 'xFwd1','xRev1','xFwd3','xRev3'
    1,10,1E-24,... % 'AXLint1','AXLint2','scaleA'
    0.1,0.1,1,1,1E-3,500]); % 'kRec','kDeg','fElse','fD2','Gas1','AXL2'

names = {'B1','B2','U1','U2','xFwd1','xRev1','xFwd3','xRev3',...
    'AXLint1','AXLint2','scaleA','kRec','kDeg','fElse','fD2',...
    'Gas1','AXL2','KD1','KD2','Y'};

KD1 = @(x) x(3)-x(1);
KD2 = @(x) x(4)-x(2);

IDX = 1;

N = 5000;
h = waitbar(0,'Starting.');
    
while IDX < N
	pp = minn + (rand(size(minn)) .* (maxx - minn));
	
	if KD1(pp) > 0
		continue;
	elseif KD2(pp) < 1
		continue;
	elseif KD2(pp) > 3
		continue;
	end
	
	try
		temp = cLib_profile ([0 10], 10.^pp, 10^pp(end), 10^pp(end-1), 1.25, 0);
		pY(IDX) = temp(2) / temp(1);
	catch
		continue;
	end
	
	outter(IDX,1:(length(pp)+4)) = [pp KD1(pp) KD2(pp) KD2(pp)-KD1(pp) log10(pY(IDX))];
	
	IDX = IDX + 1;
	
	waitbar(IDX/N,h,'Prediction...')
end

close(h)

inclIDX = [1:10 12:20];

IDXhi = outter(:,end) > 10;
namess = names(inclIDX);


plot(outter(:,inclIDX)','k')
hold on;
plot(outter(IDXhi,inclIDX)','y')

%mdl = fitglm(real(outter(:,5:(end-2))),real(outter(:,end)),'interactions','Distribution','normal','VarNames',names(5:end))