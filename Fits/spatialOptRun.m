clc; clear;

D = zeros(1,14);
D(2) = 100;
D = [0.04 1 D];


opts = gaoptimset('Display','iter', 'TimeLimit',60*5);

patternsearch(@(x) -spatialOpt(x),D,[],[],[],[],[0; zeros(15,1)],[1; 100; 0; 100*ones(13,1)],[],opts)