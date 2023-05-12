%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% More profile likelihood plots for all 1-day experimental protocols    %  
% Jana Gevertz, Last revision: 2/25/2023                                %
% - Reads in best-fit parameters and profile likelihood curves when     %
%   percent target occupancy in TME is sampled daily                    %
%   (from fit_profiles.mat)                                             %
% - Also reads in profile likelihood curves when we only have a percent %
%   target occupancy in TME at a single day                             %
%   (from profiles_sample_one_pt.mat)                                   %
% - If user selects that they want to visualize profile likelihood      %
%   curve when only a single day is sampled, and they can enter the     %
%   desired day                                                         %
% - If user selects that they want to visualize profile likelihood      %
%   curves over a range of single-day samples, they can enter the       %
%   desired day range                                                   %
% - Figure 1: konT profiles and best-fit parameter values of ksynt      %
% - Figure 2: ksynt profiles and best-fit parameter values of konT      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; %close all;

%% Load needed data
load ../fit_profiles.mat
params_allData = params;
clear params; 
load profiles_sample_one_pt.mat
konT_domain = params(1,:);
ksynt_domain = params(2,:);
Nsamples_all = Nsamples;

%% Specify what you want to plot
prompt = "Enter 1 if you want to look at one profile, 2 if you want to look at multiple profiles: ";
answ = input(prompt);
fprintf('For following questions, enter integer between %d and %d\n',1,Nsamples-1); 
if answ == 1
    prompt2 = "Enter the index of the profile you want (time = index+1 days): ";
    index = input(prompt2);
    if((index<1)||(index>Nsamples))
        fprintf('Entered index = %d at prompt2, but can only enter values between %d and %d\n',...
            index,1,Nsamples); 
        stop
    end
elseif answ == 2
    prompt2 = "Enter the the smallest index of the profile you want (time = index+1 days): ";
    index_min = input(prompt2);
    prompt3 = "Enter the the largest index of the profile you want (time = index+1 days): ";
    index_max = input(prompt3);
    if((index_min<1)||(index_min>Nsamples)||(index_max<1)||(index_max>Nsamples)||(index_min>=index_max))
        fprintf('Entered min_index = %d at prompt2, max_index = %d at prompt3, ',...
            index_min,index_max);
        fprintf('but can only enter values between %d and %d, and need min<max\n',...
            1,Nsamples); 
        stop
    end
else
    fprintf('Entered ans = %d at prompt, but can only enter 1 or 2\n',answ); 
    stop
end

%% Visualize profile when we sample at precisely one specified time point
if answ == 1
    sampling = index; 
    %day = t_day(sampling+1)/24;
    num_fit = length(best_fit_params)-1;
    threshold_fval = chi2inv(0.95,num_fit)/2 + best_fit;   
    cmap = parula(Nsamples-1);
    for k = 1:length(best_fit_params) 
        figure; 
        set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.65]);
        subplot(1,2,1)
        if k==1 
            semilogy(params(k,:),profiles_all{k}(sampling,:),'-o','LineWidth',2); 
            hold on;
            semilogy(params_allData(k,:),profiles(k,:),'-o','LineWidth',2);
            threshold_plot = threshold_fval*ones(size(params(k,:))); 
            semilogy(params(k,:),threshold_plot,'--k','LineWidth',2);
            xlim([konT_domain(1), konT_domain(end)]); 
            xlabel('k_o_n_T','FontSize',16);
            legend(['Sample at ' num2str(t_day(sampling+1)) ' Days'],'All time points',...
                'Confidence Interval','FontSize',16,'Location','SouthEast');
            fname_fig4 = ['Output/profile_konT_sample_' num2str(t_day(sampling+1)) 'days'];
        elseif k==2
            loglog(params(k,:),profiles_all{k}(sampling,:),'-o','LineWidth',2); 
            hold on;
            loglog(params_allData(k,:),profiles(k,:),'-o','LineWidth',2);
            threshold_plot = threshold_fval*ones(size(params(k,:))); 
            loglog(params(k,:),threshold_plot,'--k','LineWidth',2);
            xlim([ksynt_domain(1), ksynt_domain(end)]); 
            xlabel('k_s_y_n_t','FontSize',16);
            legend(['Sample at ' num2str(t_day(sampling+1)) ' days'],'All time points',...
                'Confidence Interval','FontSize',16,'Location','SouthEast');
            fname_fig4 = ['Output/profile_ksynt_sample_' num2str(t_day(sampling+1)) 'days'];
        end
        hold off;
        ylabel('Cost function','FontSize',16);
        title(['Profile Likelihood Curve: Sample at ' num2str(t_day(sampling+1)) ' Days'],'FontSize',16);      

        subplot(1,2,2)
        % k=1 fixed first param konT so the fit parameter is ksynt
        % k=2 fixed second param ksynt so the fit parameter is konT
        if k==1 
            semilogy(params(k,:),param_fit_all{k}(sampling,:),'o','LineWidth',2); 
            xlim([konT_domain(1), konT_domain(end)]); 
            xlabel('k_o_n_T (fixed)','FontSize',16);
            ylabel('k_s_y_n_t (fit)','FontSize',16);
        elseif k==2
            semilogx(params(k,:),param_fit_all{k}(sampling,:),'o','LineWidth',2); 
            hold on;
            yline(min_konT,'--r','LineWidth',2); 
            h1 = yline(max_konT,'--r','LineWidth',2,'DisplayName','Bounds');
            hold off;
            xlim([ksynt_domain(1), ksynt_domain(end)]); 
            ylim([min_konT,max_konT]);
            xlabel('k_s_y_n_t (fixed)','FontSize',16);
            ylabel('k_o_n_T (fit)','FontSize',16);
            legend(h1,'FontSize',16,'Location','SouthEast');
        end
        hold off;
        title(['Optimal Value of Fit Parameter: Sample at ' num2str(t_day(sampling+1)) ' Days'],...
            'FontSize',16);   
        saveas(gcf,[fname_fig4,'.fig']);
        saveas(gcf,[fname_fig4,'.png']);
    end
else
    Nsamples = index_max-index_min+1; 
    num_fit = length(best_fit_params)-1;
    threshold_fval = chi2inv(0.95,num_fit)/2 + best_fit;  
    num_pts = size(params,2);
    avg_profiles = zeros(length(best_fit_params),num_pts); 
    std_profiles = zeros(length(best_fit_params),num_pts); 
    cmap = parula(Nsamples);
    for k = 1:length(best_fit_params) 
        %% Profile likelihood from selected random sampling of size_sample time 
        %% points, and compare to profile likelihood curve for all time points
        figure; 
        set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.65]);
        subplot(1,2,1)
        for w = index_min:1:index_max
            index_cmap = w-index_min+1; 
            if k==1 
                h1 = semilogy(params(k,:),profiles_all{k}(w,:),...
                    '-','LineWidth',2,'Color', cmap(index_cmap,:));
                hold on;
                h1.Color(4) = 0.5;     
                h2 = semilogy(params_allData(k,:),profiles(k,:),'-k','LineWidth',3,...
                    'DisplayName','All time points');
                threshold_plot = threshold_fval*ones(size(params(k,:))); 
                h3 = semilogy(params(k,:),threshold_plot,'--k',...
                    'LineWidth',2,'DisplayName','Confidence Interval');
                xlabel('k_o_n_T','FontSize',16);
                legend([h2 h3],'FontSize',16,'Location','SouthEast')
                fname_fig = ['Output/profile_konT_samples_' num2str(index_min) '_to_' num2str(index_max)];
            elseif k==2
                h1 = loglog(params(k,:),profiles_all{k}(w,:),...
                    '-','LineWidth',2,'Color', cmap(index_cmap,:));
                hold on;
                h1.Color(4) = 0.5;
                h2 = loglog(params_allData(k,:),profiles(k,:),'-k',...
                    'LineWidth',3,'DisplayName','All time points');
                threshold_plot = threshold_fval*ones(size(params(k,:))); 
                h3 = loglog(params(k,:),threshold_plot,'--k','LineWidth',2,...
                    'DisplayName','Confidence Interval');
                xlabel('k_s_y_n_t','FontSize',16);
                legend([h2 h3],'FontSize',16,'Location','SouthEast');
                fname_fig = ['Output/profile_ksynt_samples_' num2str(index_min) '_to_' num2str(index_max)];
            end
        end
        hold off;
        h = colorbar;
        set(h, 'XTick', index_min:1:index_max)
        set(gca, 'CLim', [index_min, index_max])
        subtitle('Color Indicates Day Sampled','FontSize',14,'FontWeight','Bold')
        xlim([min(params(k,:)) max(params(k,:))]);
        ylabel('Cost function','FontSize',16);
        title(['Profile Likelihood Curves (Color) Using 1/' ...
            num2str(Nsamples_all-1) ' Days'],'FontSize',16);   
        
        if k==1 
            subplot(1,2,2)
            for w = index_min:1:index_max
                index_cmap = w-index_min+1; 
                semilogy(params(k,:),param_fit_all{k}(w,:),'o','LineWidth',2,...
                    'Color', cmap(index_cmap,:)); 
                hold on;
            end
            hold off;
            xlim([konT_domain(1), konT_domain(end)]); 
            xlabel('k_o_n_T (fixed)','FontSize',16);
            ylabel('k_s_y_n_t (fit)','FontSize',16); 
        else  
            subplot(1,2,2)
            for w = index_min:1:index_max
                index_cmap = w-index_min+1; 
                semilogx(params(k,:),param_fit_all{k}(w,:),'o','LineWidth',2,...
                    'Color', cmap(index_cmap,:)); 
                hold on;
            end
            yline(min_konT,'--r','LineWidth',2); 
            ylim([min_konT,max_konT]);
            h1 = yline(max_konT,'--r','LineWidth',2,'DisplayName','Bounds');
            legend(h1,'FontSize',16,'Location','SouthEast');
            xlim([ksynt_domain(1), ksynt_domain(end)]); 
            xlabel('k_s_y_n_t (fixed)','FontSize',16);
            ylabel('k_o_n_T (fit)','FontSize',16);
        end
        hold off;
        h = colorbar;
        set(h, 'XTick', index_min:1:index_max)
        set(gca, 'CLim', [index_min, index_max])
        title(['Optimal Value of Fit Parameter: Sample from ' num2str(index_min) ...
            ' to ' num2str(index_max) ' Days'],'FontSize',16);  
        subtitle('Color Indicates Day Sampled','FontSize',14,'FontWeight','Bold')
        saveas(gcf,[fname_fig,'.fig']);
        saveas(gcf,[fname_fig,'.png']);

    end
end
