clc; clear;

%function Fitting(nameIn)
nameIn = 'new_cochran2_restricted';
% This script runs extensive optimization

pwClear;

ii = 1;

baseFolder = '/Users/aaron/Documents/MATLAB/github/AXLdiffEQ/Potterswheel/';

model = nameIn;
dataFolder = [baseFolder 'data/'];

global postTag cellName siEXP;

siEXP = 0;

% postTag = 'A172';
% cellName = 'A172';
% pwAddModel(model);
% pwSelect(ii);
% pwAddData([dataFolder 'A172.xls'],1,0);
% ii = ii + 1;
% 
postTag = 'A549';
cellName = 'A549';
pwAddModel(model);
pwSelect(ii);
pwAddData([dataFolder 'A549.xls'],1,0);
ii = ii + 1;
% 
% postTag = 'A172_longT';
% cellName = 'A172';
% pwAddModel(model);
% pwSelect(ii);
% pwAddData([dataFolder 'A172_longT.xls'],1,0);
% ii = ii + 1;

postTag = 'A549_longT';
cellName = 'A549';
pwAddModel(model);
pwSelect(ii);
pwAddData([dataFolder 'A549_longT.xls'],1,0);
ii = ii + 1;

postTag = 'A549_SI';
cellName = 'A549';
siEXP = 1;
pwAddModel(model);
pwSelect(ii);
pwAddData([dataFolder 'A549_SI.xls'],1,0);
ii = ii + 1;
% 
% postTag = 'U87_SI';
% cellName = 'U87';
% siEXP = 1;
% pwAddModel(model);
% pwSelect(ii);
% pwAddData([dataFolder 'U87_SI.xls'],1,0);
% ii = ii + 1;

postTag = 'BT549_SI';
cellName = 'BT549';
siEXP = 1;
pwAddModel(model);
pwSelect(ii);
pwAddData([dataFolder 'BT549_SI.xls'],1,0);
ii = ii + 1;

pwSelect('all');



rng shuffle;
randInd = randi(1e6);
pwSetQRNGIndex(randInd);


symbols = ['a':'z' 'A':'Z' '0':'9'];
randName = symbols (randi(numel(symbols),[1 8]));


mkdir([baseFolder 'temp/' nameIn '_' randName]);
cd([baseFolder 'temp/' nameIn '_' randName]);

pwSetIntegrationStartTime(-1E6, 1);
pwSetPlotStartTime(-10,1);
pwCombine;
pwCloseFigures;

pwShowFitting(1);


