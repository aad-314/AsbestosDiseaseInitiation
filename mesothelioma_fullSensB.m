%Sensitivity Analysis Code used for AsbestosDisease 6-27-2025 
%Verified on Dec 27, 2025 at 1:52pm

clear all
close all
clf

%%%%%%%%%%%%%%%%%%%%%% PARAMETERS and INPUTS %%%%%%%%%%%%%%%%%%%%%%

kRM= .9; %recruitment rate from monocytes VALID
kF=.001; %MD fiber ingestion rate
qao=.01;
cRM= .0167; %clearance rate for RM VALID
kSS= 21; % in rate of polarization A (see auxiliary equations)
kS=  5.156; % in rate of polarization A (see auxiliary equations)
rM=  1.2e5; %normal recruitment for MO - VALID
k21= 1; % in increase recruitment rate from PIS (Voropaeva)
k22=1; % in increase recruitment rate from PIS (Voropaeva)sat
cM= 0.5; % clearance for MO VALID
cMF=1;  %MF removal rate
q=  0.1; % incomplete phagocytosis rate and also rate of sequestering fibers


k29=1; % death rate for M1 - 
k30=1; % death rate for M1 -under ILten 
k31=1; %sat ilten
k46= .512; % production rate of TGF by M2 (Voropaeva)max rate
k48= .693; % clearance rate of TGF  (Voropaeva)
cp= 1; % production rate of TNF by M1,M2, and RM
kTNF= 1; %in denominator production rate of TNF by M1,M2, and RM
k42=.1; %removal rate for TNF
k58= 10; % in term for production of IL1 by M1 (Voropaeva)max rate
k59=1; % in term for production of IL1 by M1 (Voropaeva)
k60=50; % now Tcell carrying capacity, adjust for immune response of T cells. If set to 200, default patient will not develop cancer.
k62=.001; % production rate of IL1 by M2 (Voropaeva)recycled as rate of MP


k81=1; % in term for production of ILten by M1 - Voropaeva?
k82=1; % in term for production of ILten by M1 - Voropaeva?
k83=1; % in term for production of ILten by M2 - Voropaeva?max rate
k84=1; % in term for production of ILten by M2 - Voropaeva?
k85=.1; % in term for production of ILten by M2 - Voropaeva?
c2= 1; %  death rate MNF, uptake of fiber
c1= 1; % M1 iron uptake rate
rFe= 1; %removal rate for Fe
cFen=1;
kFe= 1;
wdF= 1;
eF= 1;
q1= .0001; %rate of M to MM conversion
q3= .01;
s1= 100;
Cm= 1e10;
q2= 1;
h1= 0;
h2= 1; %unknown rate of dna damage
kTC= .0001;
nM=  0.0001; %cancer onset, dna damage rate
nMM= 0.0001; %unknown, cancer onset
pC=   .006; % proliferation rate of cancer cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   I= 0;%default (none)
% I= 100;
 I= 250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gROS= 1;
hROS= 1;
cROS= 1; %k85 Voropaeva
q4=1;
sTGF= 1;
sILten=1;
q5= .1;
cTC= 0.0001; % clearance rate of T cells
IFB=100;
pFB= 1;
cFB= 10; % return of MM to M
kTGFFB= 1; 
kMMFB=1;
wdFB=1;
dMFB= 1;
cECMFB= 1;
cECMMFB= 1;
c2ECMMFB= .001;
kFBTGF= .1;
wdECM= 1;
wdECMQ= 1;
cSC= 1;
th=.15;
cQM2= .001;
wdQ= 1;
%%%%%%%%%%%%%%%%%%%%%%%% EXPOSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EX=1095;% duration of exposure 1095=3yrs, 1460=4yrs, 2190=6 yrs, 3650=10 yrs, 4015=11 years, adjust as needed

parVec =[ kRM kF qao cRM  kSS kS rM k21 k22 cM cMF q k29 k30 k31 k46 k48 cp kTNF k42 k58 k59 k60 k62 k81 k82 k83 k84 k85 c2 c1 rFe cFen kFe wdF eF q1 q3 s1 Cm q2 h1 h2 kTC nM nMM pC I gROS hROS cROS q4 sTGF sILten q5 cTC IFB pFB cFB kTGFFB kMMFB wdFB dMFB cECMFB cECMMFB c2ECMMFB kFBTGF wdECM  wdECMQ  cSC th cQM2 wdQ EX ];
numPars = length(parVec);

%% simple tornado plot - vary each parameter up and down by 1%
% plot whole tumor volume and tnf alpha levels

sensitivity = 0.01;%0.1; % 10% sensitivity
lowcancer = zeros(1,numPars);
upcancer = zeros(1,numPars);

lowscar= zeros(1,numPars);
upscar = zeros(1,numPars);

for par = 1:numPars
    parVecHigh = parVec;
    parVecLow = parVec;
    
    parVecHigh(par) = (1+sensitivity)*parVec(par);
    parVecLow(par) = (1-sensitivity)*parVec(par);
    
    % Time
    Tf = 3000;
    tspan = 0:Tf;
    
    
    %%%%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RM0= 4.5527e+06; %resident macrophages (DMW 2017 paper)equilibrium when I=0
    MO0= 8.4478e+04; %Monocytes 
    MF0=0; %macrophage with ingested fiber
    MP0=0; %macrophage with incomplete phagocytosis
    M10=37.0296; %polarized macrophage M1
    M20=1.6942e+03; % polarized macrophage M2
    TGF0=1.2517e+03; %TGF-beta
    TNF0=0; %TNF-alpha
    IL10=411.1113; %Interlukin 1 beta (IL-1 beta)
    ILten0=124.1315; %Interlukin 10 (IL-10)
    Fe0=0; %extra-cellular iron
    FeM10=0;
    FeM20=0;
    M0=9.9977e+09;
    MM0=1.2514e+08;
    MNF0=0;
    MMNF0=0;
    MD0=0;
    MDNA0=0;
    MMDNA0=0;
    CM0=0;
    F0=0;
    ROS0=0;
    EXR0=1.0000e-03;
    TC0= 0.2021;
    FB0=1.0569;
    MFB0=1.0560;
    ECM0=6.2365e-04;
    SC0=0;
    Q0=1.6932;
    FS0=0;
    
    
    
    
    
    
    x0 = [RM0 MO0 MF0 MP0 M10 M20 TGF0 TNF0 IL10 ILten0 Fe0 FeM10 FeM20 M0 MM0 MNF0 MMNF0 MD0 MDNA0 MMDNA0 CM0 F0 ROS0 EXR0 TC0 FB0 MFB0 ECM0 SC0 Q0 FS0];
    
    [tl,xl] = ode15s(@mesothelioma_model,tspan,x0,[],parVecLow);
    [tc,xc] = ode15s(@mesothelioma_model,tspan,x0,[], parVec);
    [tu,xu] = ode15s(@mesothelioma_model,tspan,x0,[], parVecHigh);
    
    
    lRM= xl(:,1); %resident macrophages
    lMO= xl(:,2); %Monocytes   
    lMF= xl(:,3); %macrophage with ingested fiber
    lMP= xl(:,4); %macrophage with incomplete phagocytosis
    lM1= xl(:,5); %polarized macrophage M1
    lM2=xl(:,6);  % polarized macrophage M2
    lTGF= xl(:,7); %TGF-beta
    lTNF= xl(:,8); %TNF-alpha
    lIL1= xl(:,9); %Interlukin 1 beta (IL-1 beta)
    lILten= xl(:,10); %Interlukin 10 (IL-10)
    lFe= xl(:,11); %extra-cellular iron
    lFeM1=xl(:,12); %iron sequestered in M1 - 11b
    lFeM2=xl(:,13); %iron sequestered in M2 - 11c
    lM=xl(:,14); %unexposed mesothelial cells - 12
    lMM=xl(:,15);  %Mesothelial cells in mesenchymal form -13
    lMNF=xl(:,16); %Mesothelial cells near fibers -14
    lMMNF=xl(:,17); %mesenchymal form near  fibers -15
    lMD=xl(:,18); %dead mesothelial cells (after engulfing fibers) - 16
    lMDNA=xl(:,19);  %older cells from MNF with DNA damage - 17
    lMMDNA=xl(:,20); %older cells from MMNF with DNA damage - 18
    lCM=xl(:,21); %cancer cells - 19
    lF=xl(:,22);%free fibers - 20
    lROS=xl(:,23);%reactive oxygen species - 21
    lEXR=xl(:,24);%cumulative exposure at time t - 22
    lTC=xl(:,25);% T cells - 23
    lFB=xl(:,26);%fibroblasts - 24
    lMFB=xl(:,27);%myofibroblasts - 25
    lECM=xl(:,28);%extracellular matrix - 26
    lSC=xl(:,29);% scar tissue - 27
    lQ=xl(:,30);  % degrading enzyme- 28
    lFS=xl(:,31);%sequestered fibers not iron coated
    
    
    
    cRM= xc(:,1); %resident macrophages
    cMO= xc(:,2); %Monocytes   
    cMF= xc(:,3); %macrophage with ingested fiber
    cMP= xc(:,4); %macrophage with incomplete phagocytosis
    cM1= xc(:,5); %polarized macrophage M1
    cM2=xc(:,6);  % polarized macrophage M2
    cTGF= xc(:,7); %TGF-beta
    cTNF= xc(:,8); %TNF-alpha
    cIL1= xc(:,9); %Interlukin 1 beta (IL-1 beta)
    cILten= xc(:,10); %Interlukin 10 (IL-10)
    cFe= xc(:,11); %extra-cellular iron
    cFeM1=xc(:,12); %iron sequestered in M1 - 11b
    cFeM2=xc(:,13); %iron sequestered in M2 - 11c
    cM=xc(:,14); %unexposed mesothelial cells - 12
    cMM=xc(:,15);  %Mesothelial cells in mesenchymal form -13
    cMNF=xc(:,16); %Mesothelial cells near fibers -14
    cMMNF=xc(:,17); %mesenchymal form near  fibers -15
    cMD=xc(:,18); %dead mesothelial cells (after engulfing fibers) - 16
    cMDNA=xc(:,19);  %older cells from MNF with DNA damage - 17
    cMMDNA=xc(:,20); %older cells from MMNF with DNA damage - 18
    cCM=xc(:,21); %cancer cells - 19
    cF=xc(:,22);%free fibers - 20
    cROS=xc(:,23);%reactive oxygen species - 21
    cEXR=xc(:,24);%cumulative exposure at time t - 22
    cTC=xc(:,25);% T cells - 23
    cFB=xc(:,26);%fibroblasts - 24
    cMFB=xc(:,27);%myofibroblasts - 25
    cECM=xc(:,28);%extracellular matrix - 26
    cSC=xc(:,29);% scar tissue - 27
    cQ=xc(:,30);  % degrading enzyme- 28
    cFS=xc(:,31);%sequestered fibers not iron coated
    
    
    
    %% Upper
    
    
    uRM= xu(:,1); %resident macrophages
    uMO= xu(:,2); %Monocytes   
    uMF= xu(:,3); %macrophage with ingested fiber
    uMP= xu(:,4); %macrophage with incomplete phagocytosis
    uM1= xu(:,5); %polarized macrophage M1
    uM2=xu(:,6);  % polarized macrophage M2
    uTGF= xu(:,7); %TGF-beta
    uTNF= xu(:,8); %TNF-alpha
    uIL1= xu(:,9); %Interlukin 1 beta (IL-1 beta)
    uILten= xu(:,10); %Interlukin 10 (IL-10)
    uFe= xu(:,11); %extra-cellular iron
    uFeM1=xu(:,12); %iron sequestered in M1 - 11b
    uFeM2=xu(:,13); %iron sequestered in M2 - 11c
    uM=xu(:,14); %unexposed mesothelial cells - 12
    uMM=xu(:,15);  %Mesothelial cells in mesenchymal form -13
    uMNF=xu(:,16); %Mesothelial cells near fibers -14
    uMMNF=xu(:,17); %mesenchymal form near  fibers -15
    uMD=xu(:,18); %dead mesothelial cells (after engulfing fibers) - 16
    uMDNA=xu(:,19);  %older cells from MNF with DNA damage - 17
    uMMDNA=xu(:,20); %older cells from MMNF with DNA damage - 18
    uCM=xu(:,21); %cancer cells - 19
    uF=xu(:,22);%free fibers - 20
    uROS=xu(:,23);%reactive oxygen species - 21
    uEXR=xu(:,24);%cumulative exposure at time t - 22
    uTC=xu(:,25);% T cells - 23
    uFB=xu(:,26);%fibroblasts - 24
    uMFB=xu(:,27);%myofibroblasts - 25
    uECM=xu(:,28);%extracellular matrix - 26
    uSC=xu(:,29);% scar tissue - 27
    uQ=xu(:,30);  % degrading enzyme- 28
    uFS=xu(:,31);%sequestered fibers not iron coated
    
    %values for CM
    lcancer = lCM;
    ccancer = cCM;
    ucancer = uCM;
    
    lowcancer(par) = lcancer(Tf);
    cencancer  = ccancer(Tf);
    upcancer(par) = ucancer(Tf);
    
    %scar
    lscar = lSC;
    cscar = cSC;
    uscar = uSC;
    
    lowscar(par) = lscar(Tf);
    censcar  = cscar(Tf);
    upscar(par) = uscar(Tf);
end
parNames ={ 'kRM', 'kF', 'qao', 'cRM',  'kSS', 'kS', 'rM', 'k21', 'k22', 'cM', 'cMF', 'q', 'k29', 'k30', 'k31', 'k46', 'k48', 'cp', 'kTNF', 'k42', 'k58', 'k59', 'k60', 'k62', 'k81', 'k82', 'k83', 'k84', 'k85', 'c2', 'c1', 'rFe', 'cFen', 'kFe', 'wdF', 'eF', 'q1', 'q3', 's1', 'Cm', 'q2', 'h1', 'h2', 'kTC', 'nM', 'nMM', 'pC', 'I', 'gROS', 'hROS', 'cROS', 'q4', 'sTGF', 'sILten', 'q5', 'cTC', 'IFB', 'pFB', 'cFB', 'kTGFFB', 'kMMFB', 'wdFB', 'dMFB', 'cECMFB', 'cECMMFB', 'c2ECMMFB', 'kFBTGF', 'wdECM',  'wdECMQ',  'cSC', 'th', 'cQM2', 'wdQ', 'EX' };

[lowcancer, indcancer]=sort(lowcancer,'descend');
upcancer=upcancer(indcancer);
namescancer = parNames(indcancer);

[lowscar, indscar]=sort(lowscar,'descend');
upscar=upscar(indscar);
namesscar = parNames(indscar);

figure(1)
hcancer = barh(upcancer);
hold on

xmin=min([min(lowcancer),min(upcancer)]);
xmax=max([max(lowcancer),max(upcancer)]);
xlim([0.975*xmin 1.025*xmax])

barh(lowcancer)
bhcancer = get(hcancer,'BaseLine');
set(bhcancer,'BaseValue',cencancer);
title("Cancer cells, t=3000, I=" + I)
set(gca,'yticklabel',parNames)
set(gca,'Ytick',1:numPars,'YTickLabel',1:numPars)
set(gca,'yticklabel',namescancer)
xlabel('cancer  ')

figure(2)
hscar = barh(upscar);
hold on

xmin=min([min(lowscar),min(upscar)]);
xmax=max([max(lowscar),max(upscar)]);
xlim([0.975*xmin  1.025*xmax])
barh(lowscar)
bhscar = get(hscar,'BaseLine');
set(bhscar,'BaseValue',censcar);
title("Fibrosis or plaque, t=3000, I=" + I)
set(gca,'yticklabel',parNames)
set(gca,'Ytick',1:numPars,'YTickLabel',1:numPars)
set(gca,'yticklabel',namesscar)
xlabel('Fibrosis or plaque')
