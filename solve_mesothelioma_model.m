% Full model - asbestos and mesotheliomaMD
function [] = solve_mesothelioma_model
 
close all
clear all
  

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
cM=  0.5; % clearance for MO VALID
cMF=1;%MF removal rate
q=  0.1; % incomplete phagocytosis rate and also rate of sequestering fibers

k29=1; % death rate for M1
k30=1; % death rate for M1 -under ILten 
k31=1; %sat ilten
k46= .512; % production rate of TGF by M2 (Voropaeva)max rate
k48= .693; % clearance rate of TGF  (Voropaeva)
cp= 1; % production rate of TNF by M1,M2, and RM
kTNF= 1; %in denominator production rate of TNF by M1,M2, and RM
k42=0.1; %removal rate for TNF
k58= 10; % in term for production of IL1 by M1 (Voropaeva)max rate
k59=1; % in term for production of IL1 by M1 (Voropaeva)
k60=50; % Tcell carrying capacity, adjust for immune response of T cells.
k62=0.001; % production rate of IL1 by M2 (Voropaeva)recycled as rate of MP

k81=1; % in term for production of ILten by M1
k82=1; % in term for production of ILten by M1 - Voropaeva?
k83=1; % in term for production of ILten by M2 - Voropaeva?max rate
k84=1; % in term for production of ILten by M2 - Voropaeva?
k85=.1; % in term for production of ILten by M2 - Voropaeva?
c2= 1; %  death rate MNF, uptake of fiber
c1= 1; % M1 iron uptake rate
rFe= 1;%.72; %removal rate for Fe
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
kTC= 0.0001; 
nM=  0.0001; %cancer onset, dna damage rate
nMM= 0.0001; %unknown, cancer onset
pC=  0.006; % proliferation rate of cancer cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment each of these in turn and re-run to get all the figures
I=1000;
%I=100;
%I=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gROS= 1;
hROS= 1;
cROS= 1; %k85 Voropaeva
q4=1; 
sTGF= 1; 
sILten=1; 
q5= 0.1;
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
kFBTGF= 0.1; 
wdECM= 1; 
wdECMQ= 1; 
cSC= 1;
th=0.15;
cQM2= 0.001; 
wdQ= 1; 
%%%%%%%%%%%%%%%%%%%%%%%% EXPOSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EX=1095;% duration of exposure 1095=3yrs

parameter_vector =[ kRM kF qao cRM  kSS kS rM k21 k22 cM cMF q k29 k30 k31 k46 k48 cp kTNF k42 k58 k59 k60 k62 k81 k82 k83 k84 k85 c2 c1 rFe cFen kFe wdF eF q1 q3 s1 Cm q2 h1 h2 kTC nM nMM pC I gROS hROS cROS q4 sTGF sILten q5 cTC IFB pFB cFB kTGFFB kMMFB wdFB dMFB cECMFB cECMMFB c2ECMMFB kFBTGF wdECM  wdECMQ  cSC th cQM2 wdQ EX ];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  TIME RANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tf = 4000; %time range
tspan = [0:.001:9000] ;


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
 
%%%%%%%%%%%%%%%%%%%%%%%% SOLVER CHOICE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       [t,x] = ode15s(@mesothelioma_model,tspan,x0,[],parameter_vector);


%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RM= x(:,1); %resident macrophages
MO= x(:,2); %Monocytes   
MF= x(:,3); %macrophage with ingested fiber
MP= x(:,4); %macrophage with incomplete phagocytosis
M1= x(:,5); %polarized macrophage M1
M2= x(:,6);  % polarized macrophage M2
TGF= x(:,7); %TGF-beta
TNF= x(:,8); %TNF-alpha
IL1= x(:,9); %Interlukin 1 beta (IL-1 beta)
ILten= x(:,10); %Interlukin 10 (IL-10)
Fe= x(:,11); %extra-cellular iron
FeM1= x(:,12); %iron sequestered in M1 - 11b
FeM2= x(:,13); %iron sequestered in M2 - 11c
M= x(:,14); %unexposed mesothelial cells - 12
MM= x(:,15);  %Mesothelial cells in mesenchymal form -13
MNF= x(:,16); %Mesothelial cells near fibers -14
MMNF= x(:,17); %mesenchymal form near coated fibers -15
MD= x(:,18); %dead mesothelial cells (after engulfing fibers) - 16
MDNA= x(:,19);  %older cells from MNF with DNA damage - 17
MMDNA= x(:,20); %older cells from MMNF with DNA damage - 18
CM= x(:,21); %cancer cells - 19
F= x(:,22);%free fibers - 20
ROS= x(:,23);%reactive oxygen species - 21
EXR= x(:,24);%cumulative exposure at time t - 22
TC= x(:,25);% T cells - 23
FB= x(:,26);%fibroblasts - 24
MFB= x(:,27);%myofibroblasts - 25
ECM= x(:,28);%extracellular matrix - 26
SC= x(:,29);% scar tissue - 27
Q= x(:,30);  % degrading enzyme- 28
FS=x(:,31); %sequestered fibers

 MesDead = c2.*(MNF+MMNF).*F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 SC(9000)/CM(9000)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHOW VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%


% RM(8000) %resident macrophages
% MO(8000)  %Monocytes  
% MF(8000) %macrophage with ingested fiber
% MP(8000) %macrophage with incomplete phagocytosis
% M1(8000) %polarized macrophage M1
% M2(8000) % polarized macrophage M2
% TGF(8000) %TGF-beta
% TNF(8000) %TNF-alpha
% IL1(8000) %Interlukin 1 beta (IL-1 beta)
% ILten(8000) %Interlukin 10 (IL-10)
% Fe(8000) %extra-cellular iron
% FeM1(8000) %iron sequestered in M1 - 11b
% FeM2(8000) %iron sequestered in M2 - 11c
% M(8000) %unexposed mesothelial cells - 12
% MM(8000)  %Mesothelial cells in mesenchymal form -13
% MNF(8000) %Mesothelial cells near fibers -14
% MMNF(8000) %mesenchymal form near coated fibers -15
% MD(8000) %dead mesothelial cells (after engulfing fibers) - 16
% MDNA(8000)  %older cells from MNF with DNA damage - 17
% MMDNA(8000) %older cells from MMNF with DNA damage - 18
% CM(8000) %cancer cells - 19
% F(8000)%free fibers - 20
% ROS(8000)%reactive oxygen species - 21
% EXR(8000)%cumulative exposure at 0 t - 22
% TC(8000)% T cells - 23
% FB(8000)%fibroblasts - 24
% MFB(8000)%25
% ECM(8000)%extracellular matrix - 26
% SC(8000)% scar tissue - 27
% Q(8000) % degrading enzyme- 28
% FS(8000) %sequestered fibers
% 
% 



%%%%%%%%%%%%%%%%%%%%%%%  FIGURES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if I == 1000
    % Fig 2, bottom left
    figure(2) % M1, M2
        axes('FontSize', 20)
        plot(tspan, M1, '-.k','LineWidth', 1.5)  
        hold on
        plot(tspan, M2,'k','LineWidth', 1.5) 
        xlim([0 40])
        extraInputs = {'fontsize',14}; % name, value pairs
        xlabel('Time',extraInputs{:})
        ylabel('M1, M2',extraInputs{:})
        lgd=legend('M1','M2');
        fontsize(lgd, 14, "points")
        title('Polarized Macrophages (M1, M2), I=1000' )

    % Fig 3, top left
    figure(6) % IL1,TGF,ILten/100,TNF/100
        axes('FontSize', 20)
        plot(tspan, IL1,'k','LineWidth', 1.5) % 
        hold on
        plot(tspan, TGF,'-.b','LineWidth', 1.5) % 
        hold on
        plot(tspan,0.01.*ILten,'--r','LineWidth', 1.5) % 
        hold on
        plot(tspan,0.01.*TNF,'-.k','LineWidth', 1.5) 
        xlim([0 40])
        extraInputs = {'fontsize',14}; % name, value pairs
        xlabel('Time',extraInputs{:})
        ylabel(' IL1,  TGF, ILten, TNF',extraInputs{:})
        lgd=  legend(' IL1', 'TGF ', 'ILten/100', 'TNF/100');
        fontsize(lgd, 14, "points")
    
    % Fig 3, top right
    figure(7) % IL1,TGF, M1
        axes('FontSize', 20)
        plot(tspan, IL1,'k','LineWidth', 1.5) % 
        hold on
        plot(tspan, TGF,'-.b','LineWidth', 1.5) % 
        hold on
        plot(tspan,0.01.*ILten,'--r','LineWidth', 1.5) % 
        hold on
        plot(tspan,0.01.*TNF,'-.k','LineWidth', 1.5) 
        xlim([1050 1150])
        extraInputs = {'fontsize',14}; % name, value pairs
        xlabel('Time',extraInputs{:})
        ylabel(' IL1,  TGF, ILten, TNF',extraInputs{:})
        lgd=  legend(' IL1', 'TGF ', 'ILten/100', 'TNF/100');
        fontsize(lgd, 14, "points")
    
    % Fig 3, bottom left
    figure(8) % IL1,TGF, M1
        axes('FontSize', 20)
        plot(tspan, IL1,'k','LineWidth', 1.5) % 
        hold on
        plot(tspan, TGF,'-.b','LineWidth', 1.5) % 
        hold on
        plot(tspan,0.01.*ILten,'--r','LineWidth', 1.5) % 
        hold on
        plot(tspan,0.01.*TNF,'-.k','LineWidth', 1.5) 
        xlim([0 4000])
        extraInputs = {'fontsize',14}; % name, value pairs
        xlabel('Time',extraInputs{:})
        ylabel(' IL1,  TGF, ILten, TNF',extraInputs{:})
        lgd=  legend(' IL1', 'TGF ', 'ILten/100', 'TNF/100');
        fontsize(lgd, 14, "points")

    % Fig 5, top left
    figure(18) %FB, MFB
        axes('FontSize', 20)
        plot(tspan, FB,'-.k','LineWidth', 1.5)  
        hold on
        plot(tspan, MFB,'k','LineWidth', 1.5)  
        xlim([0 10])
        extraInputs = {'fontsize',14}; % name, value pairs
        xlabel('Time',extraInputs{:})
        ylabel('   FB, MFB  ',extraInputs{:})
        lgd=legend(  'FB', 'MFB ' );
        fontsize(lgd, 14, "points")
    
    % Fig 5, top right
    figure(19) % ECM, SC
        axes('FontSize', 20)
        plot(tspan, ECM,'-.k','LineWidth', 1.5)  
        hold on
        plot(tspan, SC,'k','LineWidth', 1.5)  
        xlim([0 10]) 
        extraInputs = {'fontsize',14}; % name, value pairs
        xlabel('Time',extraInputs{:})
        ylabel('   ECM,   SC  ',extraInputs{:})
        lgd =legend(  'ECM',   'SC' );
        fontsize(lgd, 14, "points")
elseif I == 100
    % Fig 2, top right
    figure(2) % M1, M2
        axes('FontSize', 20)
        plot(tspan, M1, '-.k','LineWidth', 1.5)  
        hold on
        plot(tspan, M2,'k','LineWidth', 1.5) 
        xlim([0 40])
        extraInputs = {'fontsize',14}; % name, value pairs
        xlabel('Time',extraInputs{:})
        ylabel('M1, M2',extraInputs{:})
        lgd=legend('M1','M2');
        fontsize(lgd, 14, "points")
        title('Polarized Macrophages (M1, M2), I=100' )
else
    % Fig 2, top left
    figure(2) % M1, M2
        axes('FontSize', 20)
        plot(tspan, M1, '-.k','LineWidth', 1.5)  
        hold on
        plot(tspan, M2,'k','LineWidth', 1.5) 
        xlim([0 40])
        extraInputs = {'fontsize',14}; % name, value pairs
        xlabel('Time',extraInputs{:})
        ylabel('M1, M2',extraInputs{:})
        lgd=legend('M1','M2');
        fontsize(lgd, 14, "points")
        title('Polarized Macrophages (M1, M2), Low input' )
end

% 
% figure(12) % CM
%     axes('FontSize', 20)
%       plot(tspan, CM,'k','LineWidth', 1.5)  
%        % xlim([0 40])
%        % xlim([0 2000])% to play with
%        % xlim([1050 1150])
%        xlim([0 9000])
%     extraInputs = {'fontsize',14}; % name, value pairs
%     xlabel('Time',extraInputs{:})
%      ylabel('    CM  ',extraInputs{:})
%        lgd=legend(   'CM' );
%        fontsize(lgd, 14, "points")
%  %    title('Cancerous cells (CM), raised kTC and k60, long term  dynamics')
% %title('Cancerous cells (CM), default, long term  dynamics')
% 
% figure(13) % SC
%     axes('FontSize', 20)
%       plot(tspan, SC,'k','LineWidth', 1.5)  
%          % xlim([0 40])
%         % xlim([1050 1150])
%         xlim([0 9000])
%     extraInputs = {'fontsize',14}; % name, value pairs
%     xlabel('Time',extraInputs{:})
%      ylabel('    SC  ',extraInputs{:})
%        lgd=legend(   'SC' );
%        fontsize(lgd, 14, "points")
%    % title('Plaque/Fibrosis) (SC), raised kTC and k60, long term  dynamics')
%    %title('Plaque/Fibrosis) (SC), default, long term  dynamics')
% % 
% figure(14) % TC - T cells
%     axes('FontSize', 20)
%       plot(tspan, TC,'k','LineWidth', 1.5)  
%         % xlim([0 40])
%        % xlim([1050 1150])
%        xlim([0 9000])
%     extraInputs = {'fontsize',14}; % name, value pairs
%     xlabel('Time',extraInputs{:})
%      ylabel('     TC  ',extraInputs{:})
%        lgd=legend(   ' TC' );
%        fontsize(lgd, 14, "points")
%     %  title('T cells (TC), default run, long term  dynamics')


 end
