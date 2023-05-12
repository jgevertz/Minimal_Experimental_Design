%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Compare model output (TGI) across parameter sets within 95%           %
% confidence interval for only the 2-day experimental protocols that    %
% correspond to both parameters being practically identifiable          %
% Author: Jana Gevertz. Last revision: 5/12/2023                        %
% - Reads in output and solves model at each parameterization           %
% - Compares predicted % TO in TME as a function of time across         %
%   parameterizations                                                   %
% - And compares predicted tumor trajectory as a function of time       %
%   across parameterizations                                            %
% - For each time point, the first set of plots shows all possible      %
%   trajectories, whereas the second set of plots shows the plausible   %
%   range extracted from all the trajectories                           %
% - Run time: 115 seconds                                               %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close all; tic;

%% Read in data, set ICs and set baseline parameters
[time_d, conc_ug_ml_1mpk, time_d_10mpk, conc_ug_ml_10mpk, ...
    raw_data_plasma_plasma_mg_L, raw_TME_RO, time_days_TGI, ...
    vol_mm3_control, vol_mm3_2mpk, vol_mm3_10mpk] = read_data(); 
[p, ICs] = set_parameters_ICs(); 
dose_flat = 10; % in mg/kg
dose      = dose_flat*10^6/p.MW; % mg/kg, convert to ng/ml
intvl     = 3.5; % days
dosenum   = 5;
ivdose=dose/p.V1; %if needs conversion
Tmax = 30;
[tDose,cDose,ti] = set_protocol(intvl,ivdose,dosenum,Tmax);

%% read in data from profile likelihood curves
load fit_profiles.mat
num_params = length(best_fit_params);
num_fit = num_params-1;
threshold_fval = chi2inv(0.95,num_fit)/2 + best_fit;
path = 'Randomly_Sample_Npts/All_2Pts/';
fname = [path 'profiles_sample_Npts.mat'];
fname2 = [path 'analyze_sampling.mat'];
load(fname);
load(fname2);

%% Read in data for % TO in TME
load Generate_TO_Data/TO_sim_data.mat % simulated data of % TO in tme
t_day_all = t_day; 
TO_tme_all = TO_tme;
clear t_day TO_tme;
t_day = t_day_all(1:Tmax+1); % only use data up to 31 days
TO_tme = TO_tme_all(:,1:Tmax+1);
TO_tme_avg = mean(TO_tme);
TO_tme_std = std(TO_tme);

%% Loop over desired time samplings
count_time = 1; 
params_to_test = {}; 
time_tested_params = {}; 
tumor = {}; 
TO_tme = {}; 
time_tested_params2 = {}; 
tumor2 = {}; 
TO_tme2 = {}; 
t_sample = t_sample - 1; % adjust to match day, which is what's in tpts_2id
for t = 1:1:size(t_sample,1) % loop through every protocol
    % Check if protocol index is in tpts_2id (both identifiable)
    check_identifiable = sum(ismember(tpts_2id,t_sample(t,:),'rows')); 
    if check_identifiable > 0 % both identifiable
        count = 0;
        fprintf('%d. Protocol with t1 = %d and t2 = %d is identifiable\n',...
            t,t_sample(t,1),t_sample(t,2)); 
        %% Loop over the number of profiled parameters
        for k = 1:size(params,1)

            %% Loop over the number of values tested for each profiled parameter
            for i = 1:size(params,2) 
                if k == 1
                    konT = params(1,i);
                    ksynt = param_fit_all{1}(t,i);
                    profile_value = profiles_all{1}(t,i);
                else % k = 2
                    ksynt = params(2,i);
                    konT = param_fit_all{2}(t,i);
                    profile_value = profiles_all{2}(t,i);
                end
                fprintf('\tk = %d: konT = %f, ksynt = %f, profile_value = %f\n',...
                    k,konT,ksynt,profile_value);

                %% Check that parameterization is within 95% CI
                if(profile_value<threshold_fval)
                    count = count+1; 
                    params_to_test{count_time}(count,1) = konT;
                    params_to_test{count_time}(count,2) = ksynt;
                    sol = solve_model(tDose,cDose,ICs,params_to_test{count_time}(count,:),Tmax,p);
                    Time = []; x1 = []; Ttme = []; Dtme = []; DRtme = []; 
                    for j = 1:size(sol,2)
                        Time = [Time sol{j}.x];
                        x1 = [x1 sol{j}.y(5,:)]; 
                        Ttme = [Ttme sol{j}.y(9,:)]; 
                        Dtme = [Dtme sol{j}.y(10,:)]; 
                        DRtme = [DRtme sol{j}.y(11,:)]; 
                    end
                    TOtme = 100*DRtme./(DRtme+Ttme); % TO in the TME
                    time_tested_params{count_time}{count} = Time; 
                    tumor{count_time}{count} = x1; 
                    TO_tme{count_time}{count} = TOtme; 

                    Time2 = []; x12 = []; Ttme2 = []; Dtme2 = []; DRtme2 = []; 
                    for j = 1:length(ti)-1
                        if j == 1
                            Time2 = [Time2 ti(j):0.1:ti(j+1)];
                            x12 = [x12 deval(sol{j}, ti(j):0.1:ti(j+1), 5)]; 
                            Ttme2  = [Ttme2 deval(sol{j}, ti(j):0.1:ti(j+1), 9)]; 
                            Dtme2  = [Dtme2 deval(sol{j}, ti(j):0.1:ti(j+1), 10)]; 
                            DRtme2 = [DRtme2 deval(sol{j}, ti(j):0.1:ti(j+1), 11)]; 
                        else % don't double count repeated days
                            Time2 = [Time2 ti(j)+0.1:0.1:ti(j+1)];
                            x12 = [x12 deval(sol{j}, ti(j)+0.1:0.1:ti(j+1), 5)]; 
                            Ttme2  = [Ttme2 deval(sol{j}, ti(j)+0.1:0.1:ti(j+1), 9)]; 
                            Dtme2  = [Dtme2 deval(sol{j}, ti(j)+0.1:0.1:ti(j+1), 10)]; 
                            DRtme2 = [DRtme2 deval(sol{j}, ti(j)+0.1:0.1:ti(j+1), 11)]; 
                        end
                    end
                    TOtme2 = 100*DRtme2./(DRtme2+Ttme2); % TO in the TME **** 
                    time_tested_params2{count_time}{count} = Time2; 
                    tumor2{count_time}{count} = x12; 
                    TO_tme2{count_time}{count} = TOtme2; 

                end
            end
        end
    
        %% Plot tumor trajectories over 95% CI
        cmap = parula(size(tumor{count_time},2));
        figure; 
        set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.65]);
        sgtitle(['Parameterizations within 95% CI when sampling at day t_1 = ' ...
            num2str(t_sample(t,1)) ', t_2 = ' num2str(t_sample(t,2))], ...
            'FontSize',16,'FontWeight','bold');   
        subplot(1,2,1)
        hold on;
        for i = 1:size(tumor{count_time},2) 
            plot(time_tested_params{count_time}{i},tumor{count_time}{i},'Color',cmap(i,:));
        end
        plot(time_days_TGI,vol_mm3_10mpk,'*k','MarkerSize',10);
        hold off; 
        xlabel('Time (days)','FontSize',14)
        ylabel('Tumor size','FontSize',14)
        title('Predicted tumor trajectories','FontSize',14);

        %% Plot % TO in TME trajectories over 95% CI
        subplot(1,2,2)
        hold on;
        for i = 1:size(tumor{count_time},2) 
            plot(time_tested_params{count_time}{i},TO_tme{count_time}{i},'Color',cmap(i,:));
        end
        errorbar(t_day,TO_tme_avg,TO_tme_std,'ok'); 
        hold off; 
        xlabel('Time (days)','FontSize',14)
        ylabel('Percent target occupancy in TME','FontSize',14)
        title('Predicted %TO in TME trajectories','FontSize',14);
        fname_fig2 = [path 'predictionsAll_withinCI_days' ...
            num2str(t_sample(t,1)) '_' num2str(t_sample(t,2))];
        saveas(gcf,[fname_fig2,'.fig']);
        saveas(gcf,[fname_fig2,'.png']);

        %% Find indices corresponding extrema for tumor and %TO in TME
        [find_min_tumor, find_max_tumor] = find_min_max(tumor2,count_time)
        t_bw = [time_tested_params2{count_time}{find_min_tumor}, fliplr(time_tested_params2{count_time}{find_min_tumor})];
        y_bw = [tumor2{count_time}{find_max_tumor}, fliplr(tumor2{count_time}{find_min_tumor})];

        [find_min_TOtme, find_max_TOtme] = find_min_max(tumor2,count_time)
        t_bw2 = [time_tested_params2{count_time}{find_min_TOtme}, fliplr(time_tested_params2{count_time}{find_min_TOtme})];
        y_bw2 = [TO_tme2{count_time}{find_max_TOtme}, fliplr(TO_tme2{count_time}{find_min_TOtme})];
        grayColor = [.7 .7 .7];

        %% Plot plausible ranges (as compared to all trajectories)
        figure; 
        set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.65]);
        sgtitle(['Parameterizations within 95% CI when sampling at day t_1 = ' ...
            num2str(t_sample(t,1)) ', t_2 = ' num2str(t_sample(t,2))], ...
            'FontSize',16,'FontWeight','bold');   
        subplot(1,2,1)
        fill(t_bw, y_bw,'b','FaceAlpha',0.35) %, grayColor); 
        hold on;
        plot(time_days_TGI,vol_mm3_10mpk,'*k','MarkerSize',10);
        hold off;
        xlabel('Time (days)','FontSize',14)
        ylabel('Tumor size','FontSize',14)
        title('Predicted range for tumor trajectory','FontSize',14);

        subplot(1,2,2)
        fill(t_bw2, y_bw2,'b','FaceAlpha',0.35); %, grayColor); 
        hold on;
        errorbar(t_day,TO_tme_avg,TO_tme_std,'ok'); 
        hold off; 
        xlabel('Time (days)','FontSize',14)
        ylabel('Percent target occupancy in TME','FontSize',14)
        title('Predicted range for %TO in TME trajectory','FontSize',14);
        fname_fig3 = [path 'predictionRange_withinCI_days' ...
            num2str(t_sample(t,1)) '_' num2str(t_sample(t,2))];
        saveas(gcf,[fname_fig3,'.fig']);
        saveas(gcf,[fname_fig3,'.png']);
        count_time = count_time+1; 
    end
end

fname_out = [path 'predictions_across_params_identifiable_2days.mat'];
save(fname_out,'params_to_test','time_tested_params','tumor','TO_tme');
toc

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

function [time_d, conc_ug_ml_1mpk, time_d_10mpk, conc_ug_ml_10mpk, ...
    raw_data_plasma_plasma_mg_L, raw_TME_RO, time_days_TGI, ...
    vol_mm3_control, vol_mm3_2mpk, vol_mm3_10mpk] = read_data()
    %% PK data: just to get a visual
    time_d            =  [0.000802773	0.013147279	0.023592631	0.025	0.299919666	0.69494388	1.727134555	3.403138733	4.56447194	6.496862026	6.633601177	6.389559776	6.408551325	6.680130472	7.339137214	7.482523407	9.024637166	9.162325895	9.169922514	11.35584978	11.23335429	12.91030804	13.05274466	13.33951704	13.35850859	13.62154154	15.81506543	16.09329162	16.85865103	17.75505213	18.03042959	19.3095104];
    conc_ug_ml_1mpk   =  [15.72690153	13.80115848	12.1716836	10.54030956	9.064666224	7.440888804	6.419143481	4.962491691	3.794511442	2.341657962	21.01035039	19.08080904	16.11812743	13.75178046	10.94672871	8.578482575	8.008736112	6.529294464	5.344221821	4.339568892	3.448865255	1.844079385	19.62396733	14.88747507	11.92479347	10.89165321	8.701927642	5.298642104	5.902573355	6.064001519	3.105118222	3.568512012];

    time_d_10mpk      =  [	0	0.12738835	0.12738835	0.12738835	0.38216542	0.891719561	1.019108096	3.057324656	3.057324656	4.585987076	6.751592172	6.878980707	6.751592172	6.751592172	6.751592172	7.006369242	7.515923382	9.171974338	9.299362873	9.426751408	11.2101909	11.2101909	12.86624185	13.37579599	13.24840746	13.37579599	13.63057306	13.50318453	13.7579616	14.01273867	14.01273867	14.01273867	15.92356669	15.92356669	15.92356669	17.70700618	17.96178325	17.96178325	19.23566861 19.61783421	19.74522275	];
    conc_ug_ml_10mpk  =  [135.9955872	121.3313351	114.7953874	91.91957038	83.83289622	72.56150868	57.89725657	53.66137963	43.85745806	52.52695558	30.35885267	190.5312019	159.44382	141.4699638	128.3980684	97.4355772	82.89621581	70.36551351	62.23720911	57.37687856	41.61983265	35.08388493	45.42899963	262.915782	230.1944132	202.4582657	176.3977353	160.0162358	153.5635486	130.770992	116.0651097	98.09125349	64.4019816	51.33008617	44.79413846	51.91290954	40.55826152	35.65630074	34.43861621 60.70729778	52.57899338	];

    raw_data_plasma_plasma_mg_L = [0.120858857	0.128113512	0.139820104	0.168984979	0.171465684	0.173982805	0.179128445	0.184426271	0.192668238	0.558258627	0.609270466	0.636498623	0.636498623	0.694659776	0.736357319	0.758135504	0.780557792	0.889953035	1.029575567	1.156887528	1.482129242	1.739828053	1.871336544	2.012785376	2.397434968	2.733435691	2.773562614	2.855592302	3.027001654	4.232356855	4.35753121	4.486407667	4.486407667	4.552268276	5.832070763	6.004557564	6.092704661	7.048574036	7.152047329	7.58135504	8.643882621	10.60025849	10.91376715	10.91376715	10.91376715	10.91376715	13.58036343	13.98201038	14.82129242	15.48365257	17.14656836	17.65368784	17.65368784	21.64925909	25.41343037	25.78650034	26.16504699	28.14278602	30.71438135	39.34927263	39.92692121	39.92692121	39.92692121	47.55708096	51.15178099	55.01819384	62.72899858	63.64986229	63.64986229];
    raw_TME_RO                  = [33.88083296	8.870995257	14.48292333	3.713080169	48.78137546	11.90896994	12.77332333	14.06894356	27.22531191	5.728388322	27.94408998	25.14358502	30.10315375	47.35928668	32.69985329	38.73940883	17.1778862	29.90844678	34.01458028	29.27792373	60.99150432	61.00151262	61.00606185	44.40683749	51.5336586	51.32621379	33.42954951	28.47180047	36.45387651	41.65000512	26.34185176	46.18194639	52.21968224	39.92948696	78.97460393	55.04111365	36.28191568	69.28292788	64.324269	54.19313749	18.40617785	48.17632809	72.11345775	64.56628794	45.80618013	40.63097797	69.32387094	82.91059629	43.45331923	49.06251777	50.57832065	62.2243452	54.46154196	41.10500756	62.67835818	58.15096614	60.52384364	71.52569745	76.27509184	48.47384762	51.06235855	35.53675208	25.61761461	60.34551389	32.10208468	45.26027273	59.06900042	63.38257873	45.26937118];

    
    %% TGI: H1299 strain
    time_days_TGI = [0	3	7	10	14	17	21	25];
    vol_mm3_control = [57.17217621	62.84941329	170.8501862	259.8335477	355.0138932	694.0222425	1077.648878	1833.067907	];	
    vol_mm3_2mpk=[	57.18550306	69.27295383	81.12052135	93.18131843	175.5412366	309.3960939	596.8828503	1416.404017	];
    vol_mm3_10mpk=[	57.18550306	43.64542589	49.04279917	41.9129357	47.33696268	104.2559288	128.9239237	243.3215835	];   
end

function [find_min, find_max] = find_min_max(Z,count_time)
    find_min = -1; min_size = 10^10; 
    for j = 1:size(Z{count_time},2) % number of trajectories within 95% CI
        size_compare = Z{count_time}{j}(end);
        if size_compare<min_size
            min_size = size_compare;
            find_min = j; 
        end
    end
    fprintf('At count = %d\n',count_time); 
    fprintf('\tMin size of %f at index %d\n',min_size,find_min); 
    %% Find index corresponding to largest tumor
    find_max = -1; max_size = -1; 
    for j = 1:size(Z{count_time},2) % number of trajectories within 95% CI
        size_compare = Z{count_time}{j}(end);
        if size_compare>max_size
            max_size = size_compare;
            find_max = j; 
        end
    end
    fprintf('\tMax size of %f at index %d\n',max_size,find_max); 
end
