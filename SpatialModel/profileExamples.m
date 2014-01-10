clc; clear;

xx = linspace(0,1,10000);
sshape = @(x) cos((10^x(1))*xx).^x(3);

GasConc = -1.204;
Gass = @(x) 2*sshape(x)/mean(sshape(x) .* xx)*(10^x(2));

hold off;
plot(xx,Gass([1 GasConc 60]),'r');
hold on;
plot(xx,Gass([1 GasConc 20]),'g');
plot(xx,Gass([1 GasConc 2]),'b');
plot(xx,Gass([-1.5 GasConc 60]),'r-');

