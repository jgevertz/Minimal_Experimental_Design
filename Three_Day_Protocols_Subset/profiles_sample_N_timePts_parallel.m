%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Profile likelihood for a random third of all possible N-day           %
% experimental protocols. User must set N (size_sample on line 75)      %
% Authors: Jana Gevertz and Irina Kareva. Last revision: 2/27/2023      %
% - Reads in simulated percent target occupancy in tme data             %
%   (TO_sim_data.mat) that was generated in generate_TO_data.m and      %
%   explore practical identifiability of two model parameters to        %
%   simulated percent TO in tme data:                                   %
%    1) konT = rate that drug binds target in TME                       %
%    2) ksynt = target synthesis in TME                                 %
% - In this code, practical identifiability is assessed via profile     %
%   likelihood method. Code outputs:                                    %
%    - Plot of average profile for each parameter, averaged over        %
%      Nsamples randomly sampled sets of size_sample time points        %
%    - Plot of all profiles (each random sampling of size_sample        %
%      time points) for each parameter                                  %
%    - .mat file with profile data: parameter ranges, value of          %
%      objective at each value of profiled parameter, best-fit value of %
%      non-profiled parameter. Results can be further visualized by     %
%      running analyze_samples.m                                        %
% - Run on cluster with 63 parallel pools using 50 parameter values per %
%   profile likelihood curve.                                           %
%   Run time for size_sample = 2: @1.5h using tighter bounds and        %
%                                 @2h w/weaker constraints on konT      %
%   Run time for size_sample = 3: @13.75h using tighter bounds and      %
%                                 @16.75h w/weaker constraints on konT  %
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
tf = Tmax+1;
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
num_pts = 50; %4 for quick local test 
min_konT = 0.1*p.kon; % 0.5* min value of fit parameter, konT: half binding rate in TME as plasma
max_konT = 10*p.kon; % 5* max value of fit parameter, konT: 5x binding rate in TME as plasma
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

%% Declare settings to get Nsamples random samplings of size_sample time points
size_sample = 3; 
Nsamples_all = ceil(0.33*nchoosek(Tmax,size_sample));  % 33% of all possible samplings
%Nsamples_all = 8; % 8 for quick local test
t_sample_all = zeros(Nsamples_all,size_sample);
for w = 1:Nsamples_all
    % randperm(tf,size_sample) selects size_sample integers between 1 and Tmax = 30
    % I want to skip sampling first index, so sampling from 1 to 30 and then adding 1
    % So this gives me random integers between 2 and 31. Index = N
    % corresponds to day = N - 1
    t_sample_all(w,:) = 1+sort(randperm(Tmax,size_sample));
end
% Ensures each sampling is unique
t_sample = unique(t_sample_all,'rows');
Nsamples = size(t_sample,1);
profiles_all = {};
param_fit_all = {}; % k=1 fixed first param konT so the fit parameter is ksynt
                    % k=2 fixed second param ksynt so the fit parameter is konT
for k = 1:length(best_fit_params)
    profiles_all{k} = zeros(Nsamples,num_pts); % remove t = 0 as sampling option
    param_fit_all{k} = zeros(Nsamples,num_pts); % remove t = 0 as sampling option
end

%% Open parallel pool and declare variables
npools = 63; % max 4 at home, 63 on amd cluster
parpool('local', npools); % Open distributed processing poolobj = gcp('nocreate'); % If no pool, don’t create 
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end
profiles_all_temp1 = zeros(length(best_fit_params),num_pts); 
param_fit_all_temp1 = zeros(length(best_fit_params),num_pts); 
profiles_all_temp2 = zeros(length(best_fit_params),num_pts); 
param_fit_all_temp2 = zeros(length(best_fit_params),num_pts); 

parfor w = 1:Nsamples
%     TO_tme_avg_sample = TO_tme_avg(t_sample(w,:));
%     TO_tme_std_sample = TO_tme_std(t_sample(w,:));
    profiles_temp1 = zeros(1,num_pts); 
    param_fit_temp1 = zeros(1,num_pts); 
    profiles_temp2 = zeros(1,num_pts); 
    param_fit_temp2 = zeros(1,num_pts); 

    % can stop solving ODE after all doses are given and we've reached sampled time point
    T = max(t_day(t_sample(w,end)),tDose(end)+1); 
    fprintf('Sample %d: uses data at:\n',w); 
    for i = 1:size_sample
        fprintf('\tindex = %d which gives t = %d (days)\n', t_sample(w,i),...
            t_day(t_sample(w,i)));
    end
    
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
            if k == 1
                fprintf('Fitting model, fixing konT = %f and initializing ksynt = %f with lb = %f, ub = %f\n',...
                    params(k,i),params_to_fit,minParams_profiles,maxParams_profiles); 
            else  % k = 2
                fprintf('Fitting model, fixing ksynt = %f and initializing konT = %f with lb = %f, ub = %f\n',...
                    params(k,i),params_to_fit,minParams_profiles,maxParams_profiles); 
            end
            fun_profile = @(z)objective_profile(z,k,params(k,i),tDose,...
                cDose,ICs,Tmax,TO_tme_avg,TO_tme_std,ti,p,0,t_sample(w,:));
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
                ti,p,1,t_sample(w,:));
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
    
    profiles_all_temp1(w,:) = profiles_temp1; 
    param_fit_all_temp1(w,:) = param_fit_temp1;
    profiles_all_temp2(w,:) = profiles_temp2; 
    param_fit_all_temp2(w,:) = param_fit_temp2;
    
end
profiles_all{1} = profiles_all_temp1; 
param_fit_all{1} = param_fit_all_temp1;
profiles_all{2} = profiles_all_temp2; 
param_fit_all{2} = param_fit_all_temp2;

delete(gcp('nocreate')); % close parallel pool

%% Assign color in plots by spacing between time points - only works when size_sample = 2
if size_sample==2
    % Sort by spacing between two sampled points
    t_sample_spacing = t_sample(:,2)-t_sample(:,1); 
    [t_spacing_sort, t_spacing_idx] = sort(t_sample_spacing); 
    % If spacings are equal, sort by starting time point
    t_spacing_min = t_spacing_sort(1);
    t_spacing_max = t_spacing_sort(end);
    for t_spacing = t_spacing_min:1:t_spacing_max
        fprintf('Looking for t_spacing = %f\n',t_spacing); 
        tfind = find(t_spacing_sort == t_spacing);
        if(isempty(tfind))
            fprintf('\tEmpty\n');
        else
            if(length(tfind)>1)
                %tfind
                tidx = t_spacing_idx(tfind);
                tstart = t_sample(tidx,1);
                [tstart_srt tstart_srt_idx] = sort(tstart);
                % Now all indices with same spacing need to be resorted by tsort
                for j = 1:length(tfind)
                    t_spacing_idx(tfind(j)) = tidx(tstart_srt_idx(j));
                end
            else
                fprintf('\tOnly one sampling with this spacing\n'); 
            end
        end
    end
else
    t_spacing_idx = 1:1:Nsamples; % no sorting of samplings for plots
end

%% Plots
num_fit = length(best_fit_params)-1;
% scale_fit = tf/size_sample; % to compensate for a fit done to tf data points versus only size_sample
% best_fit = best_fit/scale_fit
threshold_fval = chi2inv(0.95,num_fit)/2 + best_fit;  
avg_profiles = zeros(length(best_fit_params),num_pts); 
std_profiles = zeros(length(best_fit_params),num_pts); 
cmap = parula(Nsamples);
Ncol = ceil((Nsamples+2)/12); % at most 12 per column
for k = 1:length(best_fit_params) 
    %% Profile likelihood from each random sampling of size_sample time 
    %% points, and compare to profile likelihood curve for all time points
    figure; 
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.75, 0.75]);
    hold on; 
    if k == 1
        for w = 1:Nsamples
            w_srt = t_spacing_idx(w); 
            h1 = semilogy(params(k,:),profiles_all{k}(w_srt,:),'-','LineWidth',2,...
                'Color', cmap(w,:),'DisplayName',...
                sprintf('%d&%dd',t_day(t_sample(w_srt,1)),t_day(t_sample(w_srt,2))));
            h1.Color(4) = 0.5;  
        end
        h2 = semilogy(params_allData(k,:),profiles(k,:),'-k','LineWidth',3,...
            'DisplayName','All time points');
        threshold_plot = threshold_fval*ones(size(params(k,:))); 
        h3 = semilogy(params(k,:),threshold_plot,'--k','LineWidth',2,...
            'DisplayName','Confidence Interval');
        xlabel('k_o_n_T','FontSize',16);
        fname_fig = 'profile_konT_samples_all';
        fname_fig2 = 'profile_konT_samples_avg';
    else
        for w = 1:Nsamples
            w_srt = t_spacing_idx(w); 
            h1 = loglog(params(k,:),profiles_all{k}(w_srt,:),'-','LineWidth',2,...
                'Color', cmap(w,:),'DisplayName',...
                sprintf('%d&%dd',t_day(t_sample(w_srt,1)),t_day(t_sample(w_srt,2))));
            h1.Color(4) = 0.5;
        end
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
    
    hold off;
    legend('FontSize',16,'Location','EastOutside','NumColumns',Ncol);
    subtitle('Sorted by spacing between doses','FontSize',14,'FontWeight','Bold')
    xlim([min(params(k,:)) max(params(k,:))]);
    ylabel('Cost function','FontSize',16);
    title(['Profile Likelihood Curves (Color) Randomly Sampling ' num2str(size_sample) ...
        '/' num2str(tf-1) ' Days'],'FontSize',16);   
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
    title(['Profile Likelihood Curves Randomly Sampling ' num2str(size_sample) ...
    '/' num2str(tf-1) ' Time Points'],'FontSize',16);   
    legend('Sampling average','All time points','Confidence Interval',...
        'FontSize',16,'Location','NorthEast'); 
    saveas(gcf,[fname_fig3,'.fig']);
    saveas(gcf,[fname_fig3,'.png']);
end

toc;
diary off
save profiles_sample_Npts.mat t_day size_sample t_sample params profiles_all ...
    param_fit_all Nsamples num_pts min_konT max_konT