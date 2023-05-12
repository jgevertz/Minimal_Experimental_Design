%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Fitting model to simulated percent TO in tme data                     % 
% Authors: Jana Gevertz and Irina Kareva. Updated: 2/27/23              %   
% - Pembro model accounts for distribution in plasma, peripheral and    %
%   TME. Parameters have been calibrated to PK and %TO in TME data from %
%   Lindauer et al (2017) paper,, and TGI data on Keytruda              %
%   (https://www.ema.europa.eu/en/documents/assessment-report/...       %
%   keytruda-epar-public-assessment-report_en.pdf)                      %
% - Reads in simulated percent target occupancy in tme data             %
%   (TO_sim_data.mat) that was generated in generate_TO_data.m.         %
%   While 62 days were created, truncates to 31 days, as this is all    %
%   we can reasonably expect to have from data.                         %
% - Goal: Find best-fit of two model parameters to simulated %TO in TME %
%   data:                                                               %
%    1) konT = rate that drug binds free target in TME                  %
%    2) ksynt = rate of free target synthesis in TME                    %
%   Given this is simulated data, we already know the mean(konT) and    %
%   mean(ksynt) used to generate the simulated data. The fitting        %
%   algorithm is initialized using these mean values, which should be   %
%   quite close to the best-fit values for the average of the           %
%   simulated samples                                                   %
% - Then create profile likelihood curves for konT and ksynt,           %
%   constraining konT (binding rate in TME) to be anywhere from 1/10    %
%   binding rate in plasma to 10x the binding rate in plasma.           %
% - Found that these parameters are practically identifiable given the  %
%   simulated %TO in TME data!                                          %
% - Local time when profile has 50 parameter sets: @21 minutes          %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc; tic; 
diary on;

%% Read in data and visualize
load Generate_TO_Data/TO_sim_data.mat % simulated data of % TO in tme
t_day_all = t_day; 
TO_tme_all = TO_tme;
clear t_day TO_tme;
Tmax = 30;
t_day = t_day_all(1:Tmax+1); % only use data up to 31 days
TO_tme = TO_tme_all(:,1:Tmax+1);
TO_tme_avg = mean(TO_tme);
TO_tme_std = std(TO_tme);
figure;
plot(t_day,TO_tme,'.-');
xlabel('Time (d)','FontSize',14);
ylabel('Percent target occupancy in TME','FontSize',14);
title('Simulated Data','FontSize',16);

%% Set baseline parmaeters, protocol and dosing
[p, ICs] = set_parameters_ICs(); 
dose_flat = 10; % in mg/kg
dose      = dose_flat*10^6/p.MW; % mg/kg, convert to ng/ml
intvl     = 3.5; % days
dosenum   = 5;
T = 61; % end time for solving ODEs
ivdose=dose/p.V1; %if needs conversion
[tDose,cDose,ti] = set_protocol(intvl,ivdose,dosenum,Tmax);

%% Parameters to fit
konT = params_mean(1);
ksynt = params_mean(2); 
params_to_fit = [konT ksynt];
fit_initial = objective(params_to_fit,tDose,cDose,ICs,Tmax,TO_tme_avg,...
    TO_tme_std,ti,1,p); % last input of 1 means I want to visualize solution
fprintf('Model fit to percent TO in tme simulated data, using mean of distributions that generated the data: %f\n',fit_initial);

%% Find best-fit parameters
minParams = zeros(1,length(params_to_fit)); % all parameters have lower bound of 0
options = optimoptions(@fmincon,'MaxFunctionEvaluations', 3*10000); 
fun = @(z)objective(z,tDose,cDose,ICs,Tmax,TO_tme_avg,TO_tme_std,ti,0,p); % 0 = don't visualize solution
[best_fit_params, best_fit] = fmincon(fun,params_to_fit,[],[],[],[],minParams,[],[],options);
fprintf('konT changed from initial value of %f to %f\n',konT,best_fit_params(1));
fprintf('ksynt changed from initial value of %f to %f\n',ksynt,best_fit_params(2));
fprintf('Fit function changed from initial value of %f to %f\n',fit_initial,best_fit);
fit_best = objective(best_fit_params,tDose,cDose,ICs,Tmax,TO_tme_avg,...
    TO_tme_std,ti,1,p); % last input of 1 means I want to visualize solution

%% Profile likelihood over konT_domain and ksynt_domain
num_pts = 50; % 12 was default but 50 gives very high resolution
min_konT = 0.1*p.kon; % min value of fit parameter, konT: 1/10 binding rate in TME as plasma
max_konT = 10*p.kon; % max value of fit parameter, konT: 10x binding rate in TME as plasma
konT_domain = linspace(0.005,0.02,num_pts);
ksynt_domain = linspace(12500,16000,num_pts); 
params = [konT_domain; ksynt_domain];
minParams_profiles = zeros(1,length(best_fit_params)-1); % default lower bound of 0
maxParams_profiles = zeros(1,length(best_fit_params)-1); 
profiles = zeros(length(best_fit_params),num_pts); 
param_fit = zeros(length(best_fit_params),num_pts); 
                    % k=1 fixed first param konT so the fit parameter is ksynt
                    % k=2 fixed second param ksynt so the fit parameter is konT
% Loop over every parameter that was fit
for k = 1:length(best_fit_params)
    % When k = 1 we are fixing value of first parameter, konT
    % When k = 2 we are fixing value of second parameter, ksynt
    if k == 1
        minParams_profiles = 0; 
        maxParams_profiles = inf; % ksynt can be anything from 0 to inf
    else % k = 2
        minParams_profiles = min_konT; % min value of fit parameter, konT
        maxParams_profiles = max_konT; % max value of fit parameter, konTf
    end
    % Initialize parameter to be fit at best-fit value for all parameters besides p(k)
    params_to_fit = best_fit_params(1:end ~= k); 
    for i=1:num_pts
        fprintf('Fitting model, fixing p(%d) = %f and initializing other p = %f\n',...
            k,params(k,i),params_to_fit); 
        fun_profile = @(z)objective_profile(z,k,params(k,i),tDose,...
            cDose,ICs,Tmax,TO_tme_avg,TO_tme_std,ti,p);
        [best_fit_params_profiles, best_fit_profile] = ...
            fmincon(fun_profile,params_to_fit,[],[],[],[],minParams_profiles,...
            maxParams_profiles,[],options);
        profiles(k,i) = best_fit_profile; 
        param_fit(k,i) = best_fit_params_profiles;
        for j = 1:length(best_fit_params_profiles)
            if(j<k)
                param_num = j; 
            else
                param_num = j+1;
            end
            fprintf('\tBest-fit p(%d) = %f with fit = %f\n',...
                param_num,best_fit_params_profiles(j),best_fit_profile); 
        end
    end
end

%% Now plot each profile likelihood curve along with confidence interval threshold
num_fit = length(best_fit_params)-1;
threshold_fval = chi2inv(0.95,num_fit)/2 + best_fit;  
for k = 1:length(best_fit_params)
    figure; 
    set(gca,'fontsize', 14) 
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.65]);
    subplot(1,2,1)
    plot(params(k,:),profiles(k,:),'-o','LineWidth',3); hold on; 
    threshold_plot = threshold_fval*ones(size(params(k,:))); 
    plot(params(k,:),threshold_plot,'--','LineWidth',2); hold off;
    ylabel('Cost function','FontSize',16);
    title('Profile Likelihood Curve','FontSize',16);      
    if k==1 
        xlabel('k_o_n_T','FontSize',16);
        fname_fig = 'konT_profile_bounds';
    elseif k==2
        xlabel('k_s_y_n_t','FontSize',16);
        fname_fig = 'ksynt_profile_bounds';
    end
    xlim([min(params(k,:)) max(params(k,:))]);
    profile_min = min(profiles,[],'all');
    profile_max = max(profiles,[],'all');
    ylim([profile_min profile_max]);
    
    subplot(1,2,2)
    % k=1 fixed first param konT so the fit parameter is ksynt
    % k=2 fixed second param ksynt so the fit parameter is konT
    if k==1 
        semilogy(params(k,:),param_fit(k,:),'o','LineWidth',2); 
        xlim([konT_domain(1), konT_domain(end)]); 
        xlabel('k_o_n_T (fixed)','FontSize',16);
        ylabel('k_s_y_n_t (fit)','FontSize',16);
    elseif k==2
        semilogx(params(k,:),param_fit(k,:),'o','LineWidth',2); 
        hold on;
        yline(min_konT,'--r','LineWidth',2); 
        h1 = yline(max_konT,'--r','LineWidth',2,'DisplayName','Bounds');
        hold off;
        xlim([ksynt_domain(1), ksynt_domain(end)]); 
        ylim([min_konT,max_konT]);
        xlabel('k_s_y_n_t (fixed)','FontSize',16);
        ylabel('k_o_n_T (fit)','FontSize',16);
        legend(h1,'FontSize',14);
    end
    hold off;
    title(['Optimal Value of Fit Parameter'],'FontSize',16);   
    
    saveas(gcf,[fname_fig,'.fig']);
    saveas(gcf,[fname_fig,'.png']);
end

toc;
diary off;
save fit_profiles_bounds.mat best_fit_params konT_domain ksynt_domain ...
    best_fit params profiles param_fit


%%%%%%%%% Functions%%%%%%%%%%%%%
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

function [tDose,cDose,ti] = set_protocol(intvl,ivdose,dosenum,T)
    tDose = zeros(1,dosenum); cDose = zeros(1,dosenum);
    for i=1:dosenum
        tDose(i)=intvl*(i-1);
        cDose(i)=ivdose;
    end
    ti = [tDose max(T)];
end

function sol_all = solve_model(tDose,cDose,ICs,params_to_fit,T,p) 
    ti = [tDose max(T)];
    options_DE=odeset('RelTol',1.0e-6);
    tPoints = []; cPoints = []; 
    Dp0 = ICs(2);
    for i = 1:length(tDose)
        %fprintf('\tSolving model for dose %d\n',i);
        if i == 1
            ICs(2) = ICs(2) + cDose(i);
        else
            Dp0 = Dp0 + cDose(i);
            ICs = [Div0 Dp0 DR0 Tr0 x10 x20 x30 x40 Ttme0 Dtme0 DRtme0 Dt0]; 
        end
        tSpan = [ti(i) ti(i+1)];
        sol_all{i} = ode23s(@(t,x) pembro(t,x,params_to_fit,p),tSpan,ICs,options_DE);
        Div0 = sol_all{i}.y(1,end); Dp0 = sol_all{i}.y(2,end); DR0 = sol_all{i}.y(3,end);...
            Tr0 = sol_all{i}.y(4,end); x10 = sol_all{i}.y(5,end); x20 = sol_all{i}.y(6,end);...
            x30 = sol_all{i}.y(7,end); x40 = sol_all{i}.y(8,end); Ttme0 = sol_all{i}.y(9,end);...
            Dtme0 = sol_all{i}.y(10,end); DRtme0 = sol_all{i}.y(11,end); Dt0 = sol_all{i}.y(12,end);
    end
end

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

%% Fit function (to be minimized)
function SSE_divideby_Var = objective(params_to_fit,tDose,cDose,ICs,T,...
    TO_tme_avg,TO_tme_std,ti,visualize,p)
    sol = solve_model(tDose,cDose,ICs,params_to_fit,T,p);
    if visualize == 1 % display model solution
        Time = []; Ttme = []; Dtme = []; DRtme = []; 
        for j = 1:size(sol,2)
            Time = [Time sol{j}.x];
            Ttme = [Ttme sol{j}.y(9,:)]; 
            Dtme = [Dtme sol{j}.y(10,:)]; 
            DRtme = [DRtme sol{j}.y(11,:)]; 
        end
        TO_tme = 100*DRtme./(DRtme+Ttme); % TO in the TME
        figure;
        errorbar(0:1:T,TO_tme_avg,TO_tme_std,'.-'); hold on;
        plot(Time,TO_tme,'--'); hold off;
        xlabel('Time (d)','FontSize',14)
        ylabel('Percent target occupancy in TME','FontSize',14)
        title('Average Trajectory with STD Error Bars','FontSize',14)
        xlim([0,T+1])
    end 
    %% Extract relevant model output for fitting at "experimental" time points
    t_halfday = []; Ttme = []; Dtme = []; DRtme = []; 
    for j = 1:length(ti)-1
        if j == 1
            t_halfday = [t_halfday ti(j):0.5:ti(j+1)];
            Ttme  = [Ttme deval(sol{j}, ti(j):0.5:ti(j+1), 9)]; 
            Dtme  = [Dtme deval(sol{j}, ti(j):0.5:ti(j+1), 10)]; 
            DRtme = [DRtme deval(sol{j}, ti(j):0.5:ti(j+1), 11)]; 
        else % don't double count repeated days
            t_halfday = [t_halfday ti(j)+0.5:0.5:ti(j+1)];
            Ttme  = [Ttme deval(sol{j}, ti(j)+0.5:0.5:ti(j+1), 9)]; 
            Dtme  = [Dtme deval(sol{j}, ti(j)+0.5:0.5:ti(j+1), 10)]; 
            DRtme = [DRtme deval(sol{j}, ti(j)+0.5:0.5:ti(j+1), 11)]; 
        end
    end
    TO_tme_halfday = 100*DRtme./(DRtme+Ttme); % TO in the TME ****
    t_day = t_halfday(1:2:end);
    TO_tme = TO_tme_halfday(:,1:2:end);
        
    % Goodness of fit: divide by 2 to match likelihood formula, 
    % remove t = 0 which has std of 0 and causes nan
    SSE_divideby_Var = 0.5*sum(((TO_tme_avg(2:end)-TO_tme(2:end)).^2)./(TO_tme_std(2:end).^2));
end

%% Fit function (to be minimized) for profile likelihood
function SSE_divideby_Var = objective_profile(params_to_fit,k,params_fixed,...
    tDose,cDose,ICs,T,TO_tme_avg,TO_tme_std,ti,p)
    
    if k == 1 % k = 1 means fixing value of first parameter, konT
        params = [params_fixed params_to_fit];
    else % k = 2 means fixing value of second parameter, ksynt
        params = [params_to_fit params_fixed];
    end

    sol = solve_model(tDose,cDose,ICs,params,T,p);
    %% Extract relevant model output to compare to (simulated) %TO in TME data
    t_halfday = []; Ttme = []; Dtme = []; DRtme = []; 
    for j = 1:length(ti)-1
        if j == 1
            t_halfday = [t_halfday ti(j):0.5:ti(j+1)];
            Ttme  = [Ttme deval(sol{j}, ti(j):0.5:ti(j+1), 9)]; 
            Dtme  = [Dtme deval(sol{j}, ti(j):0.5:ti(j+1), 10)]; 
            DRtme = [DRtme deval(sol{j}, ti(j):0.5:ti(j+1), 11)]; 
        else % don't double count repeated days
            t_halfday = [t_halfday ti(j)+0.5:0.5:ti(j+1)];
            Ttme  = [Ttme deval(sol{j}, ti(j)+0.5:0.5:ti(j+1), 9)]; 
            Dtme  = [Dtme deval(sol{j}, ti(j)+0.5:0.5:ti(j+1), 10)]; 
            DRtme = [DRtme deval(sol{j}, ti(j)+0.5:0.5:ti(j+1), 11)]; 
        end
    end
    TO_tme_halfday = 100*DRtme./(DRtme+Ttme); % TO in the TME ****
    % Extract data at each day (as compared to every half day)
    t_day = t_halfday(1:2:end);
    TO_tme = TO_tme_halfday(:,1:2:end);
    
    % Goodness of fit: divide by 2 to match likelihood formula, 
    % remove t = 0 which has std of 0 and causes nan
    SSE_divideby_Var = 0.5*sum(((TO_tme_avg(2:end)-TO_tme(2:end)).^2)./(TO_tme_std(2:end).^2));
end