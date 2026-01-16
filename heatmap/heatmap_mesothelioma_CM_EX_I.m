% Heatmap for CM and SC for various inputs of I and Ex
% I=j, and Ex=c where j and c are varied in code below
 
tic;
clear all


counter1 = 0;
sol1 = [];
sol4 = zeros(1,4);
for j = 0:100:1000
    counter1 = counter1 + 1;
    disp(j)
    counter2 = 0;
    for c = 0:200:2000
        counter2 = counter2 + 1;
        %parameters and inputs
        kRM= .9; %recruitment rate from monocytes VALID
        kF=.001; %MD fiber ingestion rate
        qao=.01;
        cRM= .0167; %clearance rate for RM VALID
        kSS= 21; % in rate of polarization A (see auxiliary equations)
        kS=  5.156; % in rate of polarization A (see auxiliary equations)
        rM=  1.2e5; %normal recruitment for MO - VALID
        k21= 1; % in increase recruitment rate from PIS (Voropaeva)
        k22=1; % in increase recruitment rate from PIS (Voropaeva)sat
        cM=   0.5; % clearance for MO VALID
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
        q1= 1;
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
        I= j;
        gROS= 1;
        hROS= 1;
        cROS= 1;
        q4=1; 
        sTGF= 1;
        sILten=1;
        q5= .1;
        cTC= 0.0001; % clearance rate of T cells
        IFB=100;
        pFB= 1; 
        cFB= 10; 
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
        wdQ= 1;%4.32; 
        EX=c;
        % Create vector of parameter values that we will pass to the ode solver:
        
        parameter_vector =[ kRM kF qao cRM  kSS kS rM k21 k22 cM cMF q k29 k30 k31 k46 k48 cp kTNF k42 k58 k59 k60 k62 k81 k82 k83 k84 k85 c2 c1 rFe cFen kFe wdF eF q1 q3 s1 Cm q2 h1 h2 kTC nM nMM pC I gROS hROS cROS q4 sTGF sILten q5 cTC IFB pFB cFB kTGFFB kMMFB wdFB dMFB cECMFB cECMMFB c2ECMMFB kFBTGF wdECM  wdECMQ  cSC th cQM2 wdQ EX ];
        
        %initial conditions
        RM0=  4.6129e6; %resident macrophages (DMW 2017 paper)equilibrium when I=0
        MO0=  8.5595e4; %Monocytes    
        MF0=0; %macrophage with ingested fiber
        MP0=0; %macrophage with incomplete phagocytosis
        M10=11.2692; %polarized macrophage M1
        M20=155.3275;  % polarized macrophage M2
        TGF0=114.7586; %TGF-beta
        TNF0=0; %TNF-alpha
        IL10=9.0059; %Interlukin 1 beta (IL-1 beta)
        ILten0=124.1315; %Interlukin 10 (IL-10)
        Fe0=0; %extra-cellular iron
        FeM10=0;
        FeM20=0;
        M0= 9.9996e9;
        MM0=1.1476e11;
        MNF0=0;
        MMNF0=0;
        MD0=0;
        MDNA0=0;
        MMDNA0=0;
        CM0=0;
        F0=0;
        ROS0=0.001;
        EXR0=0;
        TC0= 0.4492;
        FB0=3.2702;
        MFB0=3.2420;
        ECM0=0.0213;
        SC0=0;
        Q0= 0.1521;
        FS0=0;
        
        x0 = [RM0 MO0 MF0 MP0 M10 M20 TGF0 TNF0 IL10 ILten0 Fe0 FeM10 FeM20 M0 MM0 MNF0 MMNF0 MD0 MDNA0 MMDNA0 CM0 F0 ROS0 EXR0 TC0 FB0 MFB0 ECM0 SC0 Q0 FS0];
        
        Tf = 9000;
        tspan = 0:Tf;
        
        sol = ode15s(@mesothelioma_model,tspan,x0,[],parameter_vector);
        
        x1 = deval(sol,tspan);
        x = x1';
        
        
        sol1(counter2,1) = j;
        sol1(counter2,2) = c;
        
        
        sol1(counter2, 3) = (x(end,21)); %x(21) for CM
        
        sol1(counter2,4) = (x(end,29));  %x(29) for SC
        
        
        clear x t parameter_vector x0 k
    end
    sol2(counter1) = j;
    sol3(counter1) = c;
    TimeMatrix1(:, counter1) = sol1(:, 3); % X end
    TimeMatrix2(:, counter1) = sol1(:, 4); % Y end
     
    sol4 = [sol4; sol1];
    clear sol1 t parameter_vector x0
    
end
sol4(1,:) = [];

%----------------------------------------------
figure(1)
[X,Y] = meshgrid(0:100:1000,0:200:2000);%j range, c range
h = surf(X,Y,TimeMatrix1);
view(0,90)
colormap(hot)
colorbar
h1 = xlabel('I=j');
h2 = ylabel('Ex=c');
h3 = title('CM at time 9000');
% xlim([1000 201000])
set( findobj(gca,'type','line'), 'LineWidth', 5);
set(gca, 'FontSize', 10,'LineWidth',4)
box on;
