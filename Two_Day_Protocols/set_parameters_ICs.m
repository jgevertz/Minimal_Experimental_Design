function [p, ICs] = set_parameters_ICs()
    %% Define variables and set initial conditions
    % 1 = Div0 (?), 2 = Dp0 (drug in plasma), 3 = DR0 (drug-target), 
    % 4 = Tr0 (target in central), 5 = x10 (growing tumor), 6 = x20 (stage 1 tumor death), 
    % 7 = x30 (stage 2 tumor death), 8 = x40 (stage 3 tumor death), 9 = Ttme0 (target in TME), 
    % 10 = Dtme0 (free drug in TME), 11 = DRtme0 (drug-target complex in TME), 
    % 12 = Dt0 (drug in peirpheral compartment)
    numVar = 12;
    ICs = zeros(1,numVar);
    p.Tr0 = 10; 
    ICs(4) = p.Tr0; 
    ICs(5) = 38; % *****fit
    ICs(9) = 4.3*p.Tr0; % Ttme0
    p.Ttme0 = ICs(9);
    
    %% For unit conversions
    p.MW_L = 150e3; % molecular weight for conversions
    p.MW   = p.MW_L/1000;% convert to ml
    
    %% PK parameters
    p.k01 = 0;
    p.V1  = 70; % (6/7)* mL/kg
    p.V2  = 33; % mL/kg
    p.Cl1 = 5; %0.25*20; % ***** ml/kg/h: CHANGED to get slower decay in plasma
    p.Cl2 = 22; % ml/kg/h
    p.k10 = p.Cl1/p.V1;
    p.k12 = p.Cl2/p.V1;
    p.k21 = p.Cl2/p.V2;
    p.KD   = .027; % pembro KD is 27 pM, or 0.027 nM    
    p.kon  = 0.005; % ***** 0.001; CHANGED made larger to get drug into TME much quicker 
    p.koff = p.KD*p.kon;
    p.k1t   = 0.3; % 1
    p.kt1   = 0.3; % 1
    p.delta = 0.0001; %very small
    p.scale_decay = 0.5; % ***** CHANGED FROM 0
        
    %% TME parameters
    p.prob = 2; % *****fit
    p.kint = 4.4; % *****fit
    p.scale_ksynt = 75; % *****fit
    p.ksynt=p.scale_ksynt*ICs(9)*p.kint; %target synthesis in TME
    p.kintP = p.kint; % assuming they are equal
    p.ksyn= ICs(4)*p.kintP; % target synthesis in plasma
    p.Kx = 500; %  ***** new parameter for saturating ksynt term
    p.K = 10000; 
    p.lam3 = 0.148542; % *****fit
    p.k1 = 100; % not used
    p.ic50 = 43;
    p.k2 = 0.38; % *****fit
end