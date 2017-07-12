


for K = 5

    ff = @(xx) (K*xx + 1 - xx).^2;
    numm = @(xx) K^2 * xx + 1 - xx;



    fplot(@(x) numm(10.^x) / ff(10.^x) ,[-3 0])
    axis([-3 0 0 2]);
    hold on;
end