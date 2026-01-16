function dxdt = mesothelioma_model(t, x, par)
    dxdt = NaN(30,1);
    %variables
    RM= x(1); %resident macrophages
    MO= x(2); %Monocytes   
    MF= x(3); %macrophage with ingested fiber
    MP= x(4); %macrophage with incomplete phagocytosis
    M1= x(5); %polarized macrophage M1
    M2=x(6);  % polarized macrophage M2
    TGF= x(7); %TGF-beta
    TNF= x(8); %TNF-alpha
    IL1= x(9); %Interlukin 1 beta (IL-1 beta)
    ILten= x(10); %Interlukin 10 (IL-10)
    Fe= x(11); %extra-cellular iron
    FeM1=x(12); %iron sequestered in M1 - 11b
    FeM2=x(13); %iron sequestered in M2 - 11c
    M=x(14); %unexposed mesothelial cells - 12
    MM=x(15);  %Mesothelial cells in mesenchymal form -13
    MNF=x(16); %Mesothelial cells near fibers -14
    MMNF=x(17); %mesenchymal form near  fibers -15
    MD=x(18); %dead mesothelial cells (after engulfing fibers) - 16
    MDNA=x(19);  %older cells from MNF with DNA damage - 17
    MMDNA=x(20); %older cells from MMNF with DNA damage - 18
    CM=x(21); %cancer cells - 19
    F=x(22);%free fibers - 20
    ROS=x(23);%reactive oxygen species - 21
    EXR=x(24);%cumulative exposure at time t - 22
    TC=x(25);% T cells - 23
    FB=x(26);%fibroblasts - 24
    MFB=x(27);%myofibroblasts - 25
    ECM=x(28);%extracellular matrix - 26
    SC=x(29);% scar tissue - 27
    Q=x(30);  % degrading enzyme- 28
    FS=x(31);%sequestered fibers not iron coated
    
    
    
    
    
    
    %parameters
    kRM= par(1); %recruitment rate from monocytes 
    kF= par(2); %fiber ingestion rate 
    qao=par(3); % in denominator of fiber ingestion term
    cRM= par(4); %clearance rate for RM
    kSS= par(5);
    kS =  par(6);
    rM= par(7); %normal recruitment for MO
    k21= par(8);% in increase recruitment rate from PIS
    k22=par(9);% in increase recruitment rate from PIS
    cM=  par(10);% clearance for MO
    cMF= par(11);
    q=   par(12);
    k29= par(13); 
    k30= par(14); 
    k31= par(15); 
    k46= par(16); 
    k48= par(17);
    cp=  par(18);
    kTNF= par(19); 
    k42= par(20); 
    k58= par(21);
    k59= par(22);
    k60= par(23); %now used as TC carrying capacity
    k62= par(24);% now as MP rate of increase
    k81=par(25);
    k82=par(26);
    k83=par(27);
    k84=par(28);
    k85=par(29);
    c2= par(30); % death of MF, MMF
    c1= par(31); % fibroblast clearance
    rFe=par(32);
    cFen=par(33); %half sat for M2 iron release
    kFe= par(34); 
    wdF= par(35); 
    eF= par(36); %clearance of Q
    q1= par(37); 
    q3= par(38); 
    s1= par(39); 
    Cm= par(40); 
    q2= par(41); 
    h1= par(42); 
    h2= par(43); 
    kTC= par(44); 
    nM= par(45); 
    nMM= par(46); 
    pC= par(47); 
    I= par(48); 
    gROS= par(49); 
    hROS= par(50); 
    cROS= par(51); 
    q4= par(52); 
    sTGF= par(53); 
    sILten= par(54); 
    q5= par(55); 
    cTC= par(56); 
    IFB= par(57); 
    pFB= par(58); 
    cFB= par(59); %scaling for M MM transition
    kTGFFB= par(60); 
    kMMFB= par(61); 
    wdFB= par(62); 
    dMFB= par(63); 
    cECMFB= par(64); 
    cECMMFB= par(65); 
    c2ECMMFB= par(66); 
    kFBTGF= par(67); 
    wdECM= par(68); 
    wdECMQ= par(69); 
    cSC= par(70); 
    th= par(71); 
    cQM2= par(72); 
    wdQ= par(73); 
    EX= par(74);
    
    
    
    
    S= (TNF +IL1 +ILten + FS);% + .1);% sum of pro-inflammatory and anti-inflammatory signals 
    
    A=(.001.*kSS.*S./(kS+S)); % rate of polarization with max kSS
    p=((TNF+IL1 + FS)./(kS+S)); % probability of choosing polarization M1
    
    p2=((ILten)./(kS+S)) ;
    
    MesDead = c2.*(MNF+MMNF).*F;
    
    
    %differential equations
    dxdt(1) = kRM.*MO- kF.*MD.*RM./(qao+MD) - cRM.*RM;
    
    dxdt(2) = rM   - A.*p.*MO - A.*(p2).*MO -kRM.*MO  -cM.*MO ;
    
    dxdt(3) = kF.*(RM+M1).*MD./(qao+MD) - cMF.*MF;
    
    dxdt(4) =  k62.*cMF.*MF; %MP
    
    dxdt(5) = A.*p.*MO - A.*(p2).*M1 + A.*p.*M2 - kF.*M1.*MD./(qao+MD) - k29.*M1;%M1
    
    dxdt(6) = A.*(p2).*MO + A.*(p2).*M1-A.*p.*M2 -k29.*M2;%M2 + k30.*ILten.*M2./(k31+ILten)
    
    dxdt(7) = k46.*M2 - k48.*TGF;
    
    dxdt(8) = cp.*(M1+M2+RM).*MD./(kTNF+MD) - k42.*TNF;
    
    dxdt(9) = k58.*M1./((k59+ILten) )  - k85.*IL1;%IL1 after voropaeva2 eq 5
    
    dxdt(10) =  k81.*TGF.*M1./(k82+TGF)+  k83.*M2./(k84+ILten)-k85.*ILten;%ILten after voropaeva2 eq 6
    
    
    dxdt(11) = 0;
    
    dxdt(12) = 0;
    
    dxdt(13) = 0;
    
    dxdt(14) = - wdF.*max(min(F,M).*M./(M+MM),0)  -(q1.*M.*TGF-cFB.*MM)./(s1+TGF) +q3.*(MM+MMNF).*(1-M./Cm).*TGF./(s1+TGF);%M - 12
    
    dxdt(15) = (q1.*M.*TGF-cFB.*MM)./(s1+TGF)- wdF.*max(min(F,MM).*MM./(M+MM),0); %MM - 13
    
    dxdt(16) = - c2.*MNF + wdF.*max(min(F,M).*M./(M+MM),0)- h2.*EXR.*MNF; %MNF - 14
    
    dxdt(17) = - c2.*MMNF + wdF.*max(min(F,MM).*MM./(M+MM),0)- h2.*EXR.*MMNF; %MMNF - 15
    
    dxdt(18)=c2.*(MNF+MMNF).*F -kF.*(RM+M1)*MD./(qao+MD);
    
    dxdt(19) = h2.*EXR.*MNF- h1.*EXR.*MDNA-kTC.*MDNA.*TC - nM.*EXR.*MDNA; %MDNA - 17
    
    dxdt(20) = h2.*EXR.*MMNF- h1.*EXR.*MMDNA-kTC.*MMDNA.*TC - nMM.*EXR.*MMDNA;  %MMDNA - 18
    
    dxdt(21) =  nM.*EXR.*MDNA  + nMM.*EXR.*MMDNA  - kTC.*CM.*TC + pC.*CM; %CM - 19
    
    dxdt(22) = I.*(heaviside( EX -t))   -c2.*(MNF+MMNF).*F -q.*F ; 
    
    dxdt(23) = hROS.*MP -cROS.*ROS; %ROS - 21
    
    dxdt(24) = ROS; %EXR - 22
    
    dxdt(25) = (q4.*(TNF.*IL1)./((sTGF+TGF).*(sILten+ILten))+q5.*TC.*CM ).*(1-TC./k60)-cTC.*TC; %TC - 23
    
    dxdt(26) = IFB.*p-pFB.*(TGF./(TGF+kTGFFB)).*FB- c1.*FB; %FB - 24
    
    dxdt(27) = pFB.*(TGF./(TGF+kTGFFB)).*FB- c1.*MFB; %MFB - 25
    
    dxdt(28) =  (c2ECMMFB.*(TGF./(kFBTGF+TGF))).*MFB -wdECMQ.*Q.*ECM;
    
    dxdt(29) = cSC.*MFB.*(max(0,ECM-th));
    
    dxdt(30) =  cQM2.*M2 - wdECMQ.*Q.*ECM-eF.*Q; %Q - 28
    
    dxdt(31)= q.*F;

end
 


