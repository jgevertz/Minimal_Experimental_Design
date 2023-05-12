function dy=pembro(t,x,params_to_fit,p)
    konT = params_to_fit(1);
    ksynt_scaled = params_to_fit(2);
    % which changes some other parameters if we want to guarantee 
    % amount of target in plasma and tumor = initial amount in absence of
    % drug
    Tpp0 = 10; 
    Ttme0 = 4.3*Tpp0; 
    p.kint = (ksynt_scaled/p.scale_ksynt)/Ttme0;
%     fprintf('konT = %f, ksynt = %f which changes kint = %f compared to kintP = %f\n',...
%         konT,ksynt_scaled,p.kint,p.kintP);

    Div=x(1);  
    Dp=x(2);  
    DR=x(3);
    Tpp=x(4);
    x1=x(5);
    x2=x(6);
    x3=x(7);
    x4=x(8);
    Ttme=x(9); %target in TME
    Dtme=x(10); %free drug in TME
    DRtme=x(11); %drug-target complex in TME
    Dt=x(12);
    wt=x1+x2+x3+x4;
    TOtme=100*DRtme./(DRtme+Ttme);
    
    %%%equations
    dDiv  = -p.k01*Div;
    %central compartment
    dDp=p.k01*Div-p.k10*Dp-p.kon*Dp*Tpp+p.koff*DR... % drug in plasma
        -p.k12*Dp+p.k21*(p.V2/p.V1)*Dt...%turnover and on-off in plasma
        -p.k1t*Dp+p.kt1*(wt+p.delta)*Dtme/p.V1; %%%% NEW PIECE, with which we replaced "-TumDisp*wt/p.V1" 
    dDR    = p.kon*Dp*Tpp-p.koff*DR-p.kintP*DR; % drug-target complex in plasma
    dTpp   = p.ksyn-p.kintP*Tpp-p.kon*Dp*Tpp+p.koff*DR; %target in plasma
    dDt    = p.k12*(p.V1/p.V2)*Dp-p.k21*Dt - p.scale_decay*p.k10*Dt; %drug in peripheral: add decay term

    %TME
    dDtme= p.k1t*p.V1*Dp/(wt+p.delta)-p.kt1*Dtme...%%%% NEW PIECE, with which we replaced "-TumDisp*wt/p.V1" 
        -konT*Dtme*Ttme+p.koff*DRtme - p.scale_decay*p.k10*Dtme; %drug in TME
    %% synthesis proportional to tumor volume
    %dTtme  = p.ksynt*(wt/x10) -p.kint*Ttme-p.prob*p.kon*Dtme*Ttme+p.koff*DRtme; %target in TME
    %% synthesis increases with tumor volume, but plateaus: 
    %% still analyze ksynt (synthesis*max rate) and konT = prob*kon (apparent affinity)
    %% possible third parameter: pk2 (related to max kill rate)
    dTtme  = ksynt_scaled*wt/(wt+p.Kx) -p.kint*Ttme-konT*Dtme*Ttme+p.koff*DRtme; %target in TME
    dDRtme = konT*Dtme*Ttme-p.koff*DRtme-p.kint*DRtme; % drug-target complex in TME
    %tum_growth = p.lamda0*x1/((1.0+(p.lamda0/p.lamda1*wt)^p.psi)^(1.0/p.psi));
    tum_growth = p.lam3*x1*(1-wt/p.K);
    kill_term  = p.k2*TOtme*x1/(p.ic50+TOtme); % max kill rate = k2>lam3 so that, at dose  = 10, we can get some shrinkage
    dx1 = tum_growth - kill_term;
    dx2 = 0*kill_term - p.k1*x2;
    dx3 = p.k1*(x2 -x3);
    dx4 = p.k1*(x3 -x4);

    dy=[dDiv;dDp;dDR;dTpp;dx1;dx2;dx3;dx4;dTtme;dDtme;dDRtme;dDt];
end