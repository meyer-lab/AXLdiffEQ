clc; clear;

xx = linspace(0,1,50);
Gass = @(x) 2*(cos((10^x(1))*xx + x(3)) + 1)/mean((cos((10^x(1))*xx + x(3)) + 1) .* xx)*(10^x(2));


D = [0; randi([0,1],[13,1])*0.1];

t = linspace(0,60,10);

names = {'Gas6','AXL','A1','A2','A12','D1','D2','AXLi','A1i','A2i','A12i','D1i','D2i','Gasi'};

outter = cLib_diff_profile (t, 359, 0.06, Gass([1.5 0 0]), D);


for ii = 1:14
    subplot(4,4,ii);
    imagesc(squeeze(outter(:,ii,:)));
end
