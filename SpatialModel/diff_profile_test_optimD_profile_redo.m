clc; clear;

xx = linspace(0,1,1000);
sshape = @(x) cos(xx*pi/3).^x(1);

GasConc = -1.204;
Gass = @(x) 2*sshape(x)/mean(sshape(x).*xx)*(10^x(2));

A = logspace(-2,3,20);

for ii = 1:length(A)
    BB(:,ii) = Gass([A(ii) -1.204]);
end

plot(BB);