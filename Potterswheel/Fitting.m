clc; clear;

%function Fitting(nameIn)
nameIn = 'new_v4c';
% This script runs extensive optimization

pwClear;

baseFolder = '/Users/aaron/Documents/MATLAB/PotterAxlModel/';

model = nameIn;
dataFolder = [baseFolder 'data/'];

global postTag cellName;

postTag = 'A172';
cellName = 'A172';
pwAddModel(model);
pwSelect(1);
pwAddData([dataFolder 'A172.xls'],1,0);

postTag = 'A549';
cellName = 'A549';
pwAddModel(model);
pwSelect(2);
pwAddData([dataFolder 'A549.xls'],1,0);

postTag = 'A172_longT';
cellName = 'A172';
pwAddModel(model);
pwSelect(3);
pwAddData([dataFolder 'A172_longT.xls'],1,0);

postTag = 'A549_longT';
cellName = 'A549';
pwAddModel(model);
pwSelect(4);
pwAddData([dataFolder 'A549_longT.xls'],1,0);

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

