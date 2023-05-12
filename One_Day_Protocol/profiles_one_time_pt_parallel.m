%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Profile likelihood for all possible 1-day experimental protocols      %  
% Authors: Jana Gevertz and Irina Kareva. Last revision: 2/27/2023      %
% - Reads in simulated percent target occupancy in tme data             %
%   TO_sim_data.mat) that was generated in generate_TO_data.m and       %
%   explore practical identifiability of two model parameters to        %
%   simulated percent TO in tme data:                                   %
%    1) konT = rate that drug binds target in TME                       %
%    2) ksynt = target synthesis in TME                                 %
% - In this code, practical identifiability is assessed via profile     %
%   likelihood method. Code outputs:                                    %
%    - Plot of average profile for each parameter, averaged over all    %
%      possible daily samplings of a single time point                  %
%    - Plot of all profiles (each possible daily sampling of a single   %
%      time point) for each parameter                                   %
%    - .mat file with profile data: parameter ranges, value of          %
%      objective at each value of profiled parameter, best-fit value of %
%      non-profiled parameter. Results can be further visualized by     %
%      running additional_plots.m                                       %
% - Addresses the following question: "For target occupancy in the      %
%   TME, is sampling only one time point sufficient?" We show the       %
%   answer is "no" - there is not a single time point we can sample so  %
%   that konT and ksynt are practically identifiable                    %
% - Ran locally, meaning only over 4 parallel pools, using 50 parameter %
%   values per profile likelihood curve. Run time for narrower konT     %
%   bounds is @4 hours and for broader bounds is @4.5 hours.            %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc; tic;
diary on

%% Read in and process simulated %TO in TME data
load ../Generate_TO_Data/TO_sim_data.mat % simulated data of % TO in tme
t_day_all = t_day; 
TO_tme_all = TO_tme;
clear t_day TO_tme;
Tmax = 30;
t_day = t_day_all(1:Tmax+1); % only use data up to 31 days
TO_tme = TO_tme_all(:,1:Tmax+1);
TO_tme_avg = mean(TO_tme);
TO_tme_std = std(TO_tme);

%% Read in profile likelihood output using data at all time points
load ../fit_profiles.mat
params_allData = params;

%% Set baseline parmaeters, protocol and dosing
[p, ICs] = set_parameters_ICs(); 
dose_flat = 10; % in mg/kg
dose      = dose_flat*10^6/p.MW; % mg/kg, convert to ng/ml
intvl     = 3.5; % days
dosenum   = 5;
ivdose=dose/p.V1; %if needs conversion
[tDose,cDose,ti] = set_protocol(intvl,ivdose,dosenum,Tmax);

%% Profile likelihood: define parameter range
num_pts = 50; 
min_konT = 0.1*p.kon; % min value of fit parameter, konT: half binding rate in TME as plasma
max_konT = 10*p.kon; % max value of fit parameter, konT: 5x binding rate in TME as plasma
clear konT_domain
konT_domain = linspace(min_konT,max_konT,num_pts);
min_ksynt = 0.8*min(ksynt_domain);
max_ksynt = 1.75*max(ksynt_domain);
clear ksynt_domain
ksynt_domain = linspace(min_ksynt,max_ksynt,num_pts);
params = [konT_domain; ksynt_domain];
minParams_profiles = zeros(1,length(best_fit_params)-1); % all parameters default to lower bound of 0
maxParams_profiles = zeros(1,length(best_fit_params)-1); % upper bound - set depending on k
options = optimoptions('fmincon','OptimalityTolerance',1e-10); %, 'Display','iter'
Nsamples = length(t_day); 
profiles_all = {};
param_fit_all = {}; % k=1 fixed first param konT so the fit parameter is ksynt
                    % k=2 fixed second param ksynt so the fit parameter is konT
for k = 1:length(best_fit_params)
    profiles_all{k} = zeros(Nsamples-1,num_pts); % remove t = 0 as sampling option
    param_fit_all{k} = zeros(Nsamples-1,num_pts); % remove t = 0 as sampling option
end

%% Open parallel pool 
npools = 4; % max 4 at home
parpool('local', npools); % Open distributed processing poolobj = gcp('nocreate'); % If no pool, don’t create 
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

%% Now consider every possible 1-day protocol for collecting %TO in TME
size_sample = 1; 
profiles_all_temp1 = zeros(length(best_fit_params),num_pts); 
param_fit_all_temp1 = zeros(length(best_fit_params),num_pts); 
profiles_all_temp2 = zeros(length(best_fit_params),num_pts); 
param_fit_all_temp2 = zeros(length(best_fit_params),num_pts); 
parfor w = 2:Nsamples
    profiles_temp1 = zeros(1,num_pts); 
    param_fit_temp1 = zeros(1,num_pts); 
    profiles_temp2 = zeros(1,num_pts); 
    param_fit_temp2 = zeros(1,num_pts); 

    % can stop solving ODE after all doses are given and we've reached sampled time point
    T = max(t_day(w),tDose(end)+1); 
    fprintf('Sample %d: uses data at t = %d (days)\n',w,t_day(w)); 
    % Loop over every parameter that was fit
    for k = 1:length(best_fit_params)
        % When k = 1 we are fixing value of first parameter, konT
        % When k = 2 we are fixing value of second parameter, ksynt
        if k == 1
            minParams_profiles = 0; 
            maxParams_profiles = inf; % ksynt can be anything from 0 to inf
        else % k = 2
            minParams_profiles = min_konT; % min value of fit parameter, konT
            maxParams_profiles = max_konT; % max value of fit parameter, konT
        end
        % Initialize parameter to be fit at best-fit value for all parameters besides p(k)
        params_to_fit = best_fit_params(1:end ~= k); 
        for i=1:num_pts
            if k == 1
                fprintf('Fitting model, fixing konT = %f and initializing ksynt = %f with lb = %f, ub = %f\n',...
                    params(k,i),params_to_fit,minParams_profiles,maxParams_profiles); 
            else  % k = 2
                fprintf('Fitting model, fixing ksynt = %f and initializing konT = %f with lb = %f, ub = %f\n',...
                    params(k,i),params_to_fit,minParams_profiles,maxParams_profiles); 
            end
            
            fun_profile = @(z)objective_profile(z,k,params(k,i),tDose,...
                cDose,ICs,Tmax,TO_tme_avg,TO_tme_std,ti,p,0,w); %t_day(w));            
            [best_fit_params_profiles, best_fit_profile,exitflag,output] = ...
                fmincon(fun_profile,params_to_fit,[],[],[],[],...
                minParams_profiles,maxParams_profiles,[],options);
            if k == 1
                profiles_temp1(i) = best_fit_profile;
                param_fit_temp1(i) = best_fit_params_profiles;
            else
                profiles_temp2(i) = best_fit_profile;
                param_fit_temp2(i) = best_fit_params_profiles;
            end
            %profiles_all{k}(w-1,i) = best_fit_profile; 
            %param_fit_all{k}(w-1,i) = best_fit_params_profiles;
            % Temporary to debugging
            fit_temp = objective_profile(best_fit_params_profiles,k,...
                params(k,i),tDose,cDose,ICs,Tmax,TO_tme_avg,TO_tme_std,...
                ti,p,1,w); %t_day(w)); 
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
    
    profiles_all_temp1(w-1,:) = profiles_temp1; 
    param_fit_all_temp1(w-1,:) = param_fit_temp1;
    profiles_all_temp2(w-1,:) = profiles_temp2; 
    param_fit_all_temp2(w-1,:) = param_fit_temp2;
    
end
profiles_all{1} = profiles_all_temp1; 
param_fit_all{1} = param_fit_all_temp1;
profiles_all{2} = profiles_all_temp2; 
param_fit_all{2} = param_fit_all_temp2;

delete(gcp('nocreate')); % close parallel pool

%% Plots
num_fit = length(best_fit_params)-1;
threshold_fval = chi2inv(0.95,num_fit)/2 + best_fit;  
avg_profiles = zeros(length(best_fit_params),num_pts); 
std_profiles = zeros(length(best_fit_params),num_pts); 
cmap = parula(Nsamples-1);
for k = 1:length(best_fit_params) 
    %% Profile likelihood from each random sampling of size_sample time 
    %% points, and compare to profile likelihood curve for all time points
	figure; hold on; 
    for w = 1:Nsamples-1 
        if k==1 
            h1 = plot(params(k,:),profiles_all{k}(w,:),'-','LineWidth',2,...
                'Color', cmap(w,:));
            h1.Color(4) = 0.5;     
            h2 = plot(params_allData(k,:),profiles(k,:),'-k','LineWidth',3,...
                'DisplayName','All time points');
            threshold_plot = threshold_fval*ones(size(params(k,:))); 
            h3 = plot(params(k,:),threshold_plot,'--k','LineWidth',2,...
                'DisplayName','Confidence Interval');
            %xlim([konT_domain(1), konT_domain(end)]); 
            xlabel('k_o_n_T','FontSize',16);
            fname_fig = 'profile_konT_samples_all';
            fname_fig2 = 'profile_konT_samples_avg';
        elseif k==2
            h1 = loglog(params(k,:),profiles_all{k}(w,:),'-','LineWidth',2,...
                'Color', cmap(w,:));    
            h1.Color(4) = 0.5;
            h2 = loglog(params_allData(k,:),profiles(k,:),'-k','LineWidth',3,...
                'DisplayName','All time points');
            threshold_plot = threshold_fval*ones(size(params(k,:))); 
            h3 = loglog(params(k,:),threshold_plot,'--k','LineWidth',2,...
                'DisplayName','Confidence Interval');
            %xlim([ksynt_domain(1), ksynt_domain(end)]); 
            xlabel('k_s_y_n_t','FontSize',16);
            fname_fig = 'profile_ksynt_samples_all';
            fname_fig2 = 'profile_ksynt_samples_avg';
        end
    end
    hold off;
    legend([h2 h3],'FontSize',16,'Position',[0.45 0.65 0.1 0.2]); 
    colorbar;
    set(gca, 'CLim', [1, Nsamples])
    subtitle('Color Indicates Day Sampled','FontSize',14,'FontWeight','Bold')
    xlim([min(params(k,:)) max(params(k,:))]);
    ylabel('Cost function','FontSize',16);
    title(['Profile Likelihood Curves (Color) Using ' num2str(size_sample) ...
        '/' num2str(Nsamples-1) ' Days'],'FontSize',16);   
    saveas(gcf,[fname_fig,'.fig']);
    saveas(gcf,[fname_fig,'.png']);
    
    %% Average of profile likelihood (with error bars) from random samplings of
    %% size_sample time points, and compare to profile likelihood for all time points
    avg_profiles(k,:) = mean(profiles_all{k});
    std_profiles(k,:) = std(profiles_all{k});
    figure; hold on;
    errorbar(params(k,:), avg_profiles(k,:),std_profiles(k,:),'-o',...
        'LineWidth',3)
    if k==1 
        xlabel('k_o_n_T','FontSize',16);
        fname_fig3 = 'profile_konT_samples_avg';
    elseif k==2
            xlabel('k_s_y_n_t','FontSize',16);
        fname_fig3 = 'profile_ksynt_samples_avg';
    end
    plot(params_allData(k,:),profiles(k,:),'-ok','LineWidth',3);
    plot(params(k,:),threshold_plot,'--k','LineWidth',2);
    xlim([min(params(k,:)) max(params(k,:))]);
    ylim([0 inf])
    ylabel('Cost function','FontSize',16);
    hold off;
    title(['Profile Likelihood Curves Using ' num2str(size_sample) ...
    '/' num2str(Nsamples-1) ' Time Points'],'FontSize',16);   
    legend('Sampling average','All time points','Confidence Interval',...
        'FontSize',16,'Position',[0.45 0.65 0.1 0.2]); 
    saveas(gcf,[fname_fig3,'.fig']);
    saveas(gcf,[fname_fig3,'.png']);
end

toc;
diary off
save profiles_sample_one_pt.mat t_day params profiles_all param_fit_all ...
    Nsamples num_pts min_konT max_konT %min_ksynt max_ksynt 

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

%% Fit function (to be minimized) for profile likelihood
function SSE_divideby_Var = objective_profile(params_to_fit,k,params_fixed,...
    tDose,cDose,ICs,T,TO_tme_avg,TO_tme_std,ti,p,flag,t_sample)
    
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
    %idx = find(t_day == Tprof);

    %% Now pull out data at only time points in t_sample
    % Goodness of fit: divide by 2 to match likelihood formula
    %t_day_prof = t_sample - 1; 
    SSE_divideby_Var = 0.5*sum(((TO_tme_avg(t_sample)-TO_tme(t_sample)).^2)./(TO_tme_std(t_sample).^2));
    %SSE_divideby_Var = 0.5*sum(((TO_tme_avg(idx)-TO_tme(idx)).^2)./(TO_tme_std(idx).^2));
    if flag == 1
        fprintf('\tAnd an error of %f\n', SSE_divideby_Var); 
    end
end 