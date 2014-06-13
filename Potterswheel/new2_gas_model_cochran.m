function m = new2_gas_model_cochran()
% PottersWheel model definition file

global postTag cellName siEXP;
gasName = ['Gas_' cellName];
axlExp = ['AXLexp_' cellName];
m = pwGetEmptyModel();
m.modelFormat = 3.0; m.ID = ['AXL_' postTag]; m.name = ['AXL_' postTag];
%m.amountBasedRates = true;

%% Dynamic variables
% m = pwAddX(m, ID, startValue, type, minValue, maxValue, unit, compartment, name, description, typeOfStartValue)
m = pwAddX(m,'AXL'      ,1.5e5,'fix',0, 1E20,'#/cell','cell');
m = pwAddX(m,'A<1|2|12>',0,'fix',0, 1E20,'#/cell','cell');
m = pwAddX(m,'D<1|2>'   ,0,'fix',0, 1E20,'#/cell','cell');

m = pwAddX(m,'AXLi'      ,1.5e5,'fix',  0, 1E20,'#/cell','endoMem');
m = pwAddX(m,'A<1|2|12>i',0,'fix',  0, 1E20,'#/cell','endoMem');
m = pwAddX(m,'D<1|2>i'   ,0,'fix',  0, 1E20,'#/cell','endoMem');
m = pwAddX(m,'Gasi'      ,0,'fix',  0, 1E20,'#/cell','endoMem');

%% Reactions2
% m = pwAddR(m, ID, reactants, products, modifiers,  type, options, rateSignature, parameters, description, name, fast)
if siEXP
    m = pwAddR(m,'R01',{'AXL'},{'A1'},{'Gas'},'C',[],'k1*r1*(m1+k3)-k2*p1',{'B1','U1',gasName});
    m = pwAddR(m,'R02',{'AXL'},{'A2'},{'Gas'},'C',[],'k1*r1*(m1+k3)-k2*p1',{'B2','U2',gasName});
    m = pwAddR(m,'R03',{'A1'},{'A12'},{'Gas'},'C',[],'k1*r1*(m1+k3)-k2*p1',{'B2','U2',gasName});
    m = pwAddR(m,'R04',{'A2'},{'A12'},{'Gas'},'C',[],'k1*r1*(m1+k3)-k2*p1',{'B1','U1',gasName});
else
    m = pwAddR(m,'R01',{'AXL'},{'A1'},{'Gas'},'C',[],'k1*r1*(m1*k3)-k2*p1',{'B1','U1',gasName});
    m = pwAddR(m,'R02',{'AXL'},{'A2'},{'Gas'},'C',[],'k1*r1*(m1*k3)-k2*p1',{'B2','U2',gasName});
    m = pwAddR(m,'R03',{'A1'},{'A12'},{'Gas'},'C',[],'k1*r1*(m1*k3)-k2*p1',{'B2','U2',gasName});
    m = pwAddR(m,'R04',{'A2'},{'A12'},{'Gas'},'C',[],'k1*r1*(m1*k3)-k2*p1',{'B1','U1',gasName});
end
m = pwAddR(m,'R05',{'AXL','A1'},{'D1'},[],'C',[],'k1*r1*r2-k2*p1',{'xFwd1','xRev1'});
m = pwAddR(m,'R07',{'AXL','A12'},{'D2'},[],'C',[],'k1*r1*r2-k2*p1',{'xFwd3','xRev3'});
m = pwAddR(m,'R06',{'AXL','A2'},{'D1'},[],'C',[],'k1*k2*r1*r2/k3-k4*p1*k5/k6',{'xFwd1','B1','B2','xRev1','U1','U2'}); 
% xRev2 is xRev1*U1/U2 _________ xFwd2 is xFwd1*B1/B2
m = pwAddR(m,'R08',{'A1','A1'}, {'D2'},[],'C',[],'k1*k2*r1*r2/k3-k4*p1*k5/k6',{'xFwd3','B2','B1','xRev3','U2','U1'}); 
% xRev4 is xRev3*U2/U1 _________ xFwd4 is xFwd3*B2/B1
m = pwAddR(m,'R09',{'A2','A2'}, {'D2'},[],'C',[],'k1*k2*r1*r2/k3-k4*k5*p1/k6',{'xFwd3','B1','B2','xRev3','U1','U2'});
% xRev5 is xRev3*U1/U2 _________ xFwd5', 'xFwd3*B1/B2
if siEXP
    m = pwAddR(m,'R11',{'D1'},{'D2'},{'Gas'},'C',[],'k1*r1*k2*(m1*k7)/k3-k4*k5*p1/k6',{'xFwd3','B2','xFwd1','xRev3','U2','xRev1',gasName});
else
    m = pwAddR(m,'R11',{'D1'},{'D2'},{'Gas'},'C',[],'k1*r1*k2*(m1+k7)/k3-k4*k5*p1/k6',{'xFwd3','B2','xFwd1','xRev3','U2','xRev1',gasName});
end
% xFwd6', 'xFwd3*B2/xFwd1'); _________ xRev6', 'xRev3*U2/xRev1');

% Internalization
m = pwAddR(m,'R12',{'AXL'},{'AXLi'},[],'C',[],'(k1 + k2*k3)*r1 - k4*(1-k5)*p1',{'AXLint1','AXLint2','scale_A','kRec','fElse'});
m = pwAddR(m,'R13',{'A1'}, {'A1i'}, [],'C',[],'(k1 + k2*k3)*r1 - k4*(1-k5)*p1',{'AXLint1','AXLint2','scale_A','kRec','fElse'});
m = pwAddR(m,'R14',{'A2'}, {'A2i'}, [],'C',[],'(k1 + k2*k3)*r1 - k4*(1-k5)*p1',{'AXLint1','AXLint2','scale_A','kRec','fElse'});
m = pwAddR(m,'R15',{'A12'},{'A12i'},[],'C',[],'(k1 + k2*k3)*r1 - k4*(1-k5)*p1',{'AXLint1','AXLint2','scale_A','kRec','fElse'});
m = pwAddR(m,'R16',{'D1'}, {'D1i'}, [],'C',[],'(k1 + k2*k3)*r1 - k4*(1-k5)*p1',{'AXLint1','AXLint2','scale_A','kRec','fElse'});
m = pwAddR(m,'R17',{'D2'}, {'D2i'}, [],'C',[],'(k1 + k2)*r1 - k3*(1-k4)*p1',{'AXLint1','AXLint2','kRec','fD2'});
m = pwAddR(m,'R19',[],{'AXL'},      [],'C',[],'k1',{axlExp});

m = pwAddR(m,'R26',{'AXLi'},[],[],'C',[],'k1*k2*r1',{'kDeg','fElse'});
m = pwAddR(m,'R27',{'A1i'}, [],[],'C',[],'k1*k2*r1',{'kDeg','fElse'});
m = pwAddR(m,'R28',{'A2i'}, [],[],'C',[],'k1*k2*r1',{'kDeg','fElse'});
m = pwAddR(m,'R29',{'A12i'},[],[],'C',[],'k1*k2*r1',{'kDeg','fElse'});
m = pwAddR(m,'R30',{'D1i'}, [],[],'C',[],'k1*k2*r1',{'kDeg','fElse'});
m = pwAddR(m,'R31',{'D2i'}, [],[],'C',[],'k1*k2*r1',{'kDeg','fD2'});

m = pwAddR(m,'R32',{'AXLi','Gasi'},{'A1i'},[],'C',[],'k1*r1*r2/623-k2*p1',{'B1','U1'});
m = pwAddR(m,'R33',{'AXLi','Gasi'},{'A2i'},[],'C',[],'k1*r1*r2/623-k2*p1',{'B2','U2'});
m = pwAddR(m,'R34',{'A1i','Gasi'},{'A12i'},[],'C',[],'k1*r1*r2/623-k2*p1',{'B2','U2'});
m = pwAddR(m,'R35',{'A2i','Gasi'},{'A12i'},[],'C',[],'k1*r1*r2/623-k2*p1',{'B1','U1'});
m = pwAddR(m,'R36',{'AXLi','A1i'}, {'D1i'},[],'C',[],'k1*r1*r2-k2*p1',{'xFwd1','xRev1'});
m = pwAddR(m,'R37',{'AXLi','A12i'},{'D2i'},[],'C',[],'k1*r1*r2-k2*p1',{'xFwd3','xRev3'});
m = pwAddR(m,'R38',{'AXLi','A2i'}, {'D1i'},[],'C',[],'k1*k2*r1*r2/k3-k4*p1*k5/k6',{'xFwd1','B1','B2','xRev1','U1','U2'}); 
m = pwAddR(m,'R39',{'A1i','A1i'},  {'D2i'},[],'C',[],'k1*k2*r1*r2/k3-k4*p1*k5/k6',{'xFwd3','B2','B1','xRev3','U2','U1'}); 
m = pwAddR(m,'R40',{'A2i','A2i'},  {'D2i'},[],'C',[],'k1*k2*r1*r2/k3-k4*k5*p1/k6',{'xFwd3','B1','B2','xRev3','U1','U2'});
m = pwAddR(m,'R41',{'D1i','Gasi'}, {'D2i'},[],'C',[],'k1*r1*k2*r2/k3/623-k4*k5*p1/k6',{'xFwd3','B2','xFwd1','xRev3','U2','xRev1'});
m = pwAddR(m,'R42',{'Gasi'},[],[],'C',[],'k1*r1',{'kDeg'});


%% Observables
% m = pwAddY(m, ID, rhs, errorModel, noiseType, unit, name, description)

EndoMemFrac = 0.5;

m = pwAddY(m,'GASamt',['scale_gas*(A1+A2+2*A12+D1+2*D2+Gasi + ' mat2str(EndoMemFrac) '*(A1i+A2i+2*A12i+D1i+2*D2i))']);

if ~isempty(strfind(postTag,'longT'))
    m = pwAddY(m,'dimAXL',['scale_obs1*(2*D2 + 2*' mat2str(EndoMemFrac) '*D2i +scale_A*(AXL+A1+A2+A12+2*D1 + ' mat2str(EndoMemFrac) '*(AXLi+A1i+A2i+A12i+2*D1i)))/(A1+A2+A12+2*D1+2*D2+AXL + ' mat2str(EndoMemFrac) '*(A1i+A2i+A12i+2*D1i+2*D2i+AXLi))']);
    
    m = pwAddY(m,'tAXL',['(A1+A2+A12+2*D1+2*D2+AXL + ' mat2str(EndoMemFrac) '*(A1i+A2i+A12i+2*D1i+2*D2i+AXLi))/135.2']); % Scaling for fg/mg
else
    m = pwAddY(m,'dimAXL',['scale_obs1*(2*D2 + 2*' mat2str(EndoMemFrac) '*D2i +scale_A*(AXL+A1+A2+A12+2*D1 + ' mat2str(EndoMemFrac) '*(AXLi+A1i+A2i+A12i+2*D1i)))']);

    if siEXP
        m = pwAddY(m,'tAXL',['scale_tot*(A1+A2+A12+2*D1+2*D2+AXL + ' mat2str(EndoMemFrac) '*(A1i+A2i+A12i+2*D1i+2*D2i+AXLi))/135.2']); % Scaling for fg/mg
    end
end

%% Observation parameters
m = pwAddS(m, 'scale_obs1', 2543,  'local',  1E-6, 1E6);
m = pwAddS(m, 'scale_gas', 2543,  'local',  1E-6, 1E6);

if ~isempty(strfind(postTag,'longT'))
    m = pwAddS(m, 'scale_tot',    0.014545, 'local', 1E-6, 1);
end
    
m = pwAddS(m, 'scale_A',    0.014545, 'global', 1E-6, 1);

%% Driving input, Compartments
% m = pwAddU(m, ID, uType, uTimes, uValues, compartment, name, description, u2Values, alternativeIDs)
% m = pwAddC(m, *ID, *size, outside, spatialDim, name, unit, constant, designerProps, classname, description)

m = pwAddU(m,'Gas','steps',[-1000 0],[0 1.25]);
m = pwAddC(m,'cell',1);
m = pwAddC(m,'endoMem',0.5);

%% K - Dynamic parameters
% m = pwAddK(m, *ID, *value, fitSetting, minValue, maxValue, unit, name, description)
m = pwAddK(m,'B1',1.2,'fix'); % Limits of 10^4 to 10^8
m = pwAddK(m,'B2',0.30673,'global',0.0006,0.6);
m = pwAddK(m,'U1',0.042,'fix');
m = pwAddK(m,'U2',3.3865,'global');
m = pwAddK(m,'xFwd1',0.0040926);
m = pwAddK(m,'xRev1',0.26397);
m = pwAddK(m,'xFwd3',0.00081641);
m = pwAddK(m,'xRev3',0.041195);

m = pwAddK(m,'AXLint1',0.0020244,'global',1E-6,1);
m = pwAddK(m,'AXLint2',0.99998,'global',1E-3,100);
m = pwAddK(m,'fElse',0.00089687,'global',1E-5,1);
m = pwAddK(m,'fD2',0.76196,'global',1E-1,1);
m = pwAddK(m,'kRec',0.0026825,'global',1E-3,1E-1);
m = pwAddK(m,'kDeg',0.055024,'global',1E-4,1E-1);

m = pwAddK(m,axlExp,39.571,'global',1,1E5);
m = pwAddK(m,gasName,0.015242,'global',1E-6,100);

%% Constraints
% m = pwAddCS(m, *ID, *lhs, *operator, *rhs, lambda)
% m = pwAddCS(m,'CS2','U2/B2','<','1000',1000);
% m = pwAddCS(m,'CS3','U2/B2','>','10'  ,1000);
