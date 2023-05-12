%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% 12/20/22: Statistics of profile likelihood on randomly sampled data   %
% - Uses the output from profile_sample_N_timePts_parallel.m, which     %
%   gives the profile likelihood curve when the data is sampled at      %
%   size_sample number of time points Nsamples times.                   %
% - Identifies just the sampled time points for which both parameters   %
%   are practically identifiable, makes three plots:                    %
%   Figure 1: the profile likelihood of konT in all identifiable cases  %
%   Figure 2: the profile likelihood of ksynt in all identifiable cases %
%   Figure 3: scatter plot of sampled time points that result in both   %
%             parameters being practically identifiable
% - Note: when >2 points are sampled, the code uses gplotmatrix to      %
%   visualize the higher dimensional samplings. I can't seem to control %
%   all features of this plot (thickness of stairs, xlim, ylim, where   %
%   legend goes), so this plot needs to be prettied up manually and     %
%   re-saved after doing so                                             %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc;
%% Comment if not running directly after running OV_profile_likelihood_randomSampling.m
load ../Generate_TO_Data/TO_sim_data.mat
TO_tme_avg = mean(TO_tme);
TO_tme_std = std(TO_tme);
load ../fit_profiles.mat
params_allData = params;
clear params; 
path = 'All_2Pts/' ;
%path = 'Sample_3Pts/' ;
%path = 'All_Domain_3Pts/';
fname = [path 'profiles_sample_Npts.mat'];
load(fname);
tf = length(t_day)
Nsamples = size(t_sample,1);
size_sample = size(t_sample,2);
num_params = length(best_fit_params);
num_fit = num_params-1;
threshold_fval = chi2inv(0.95,num_fit)/2 + best_fit;
practically_identifiable = zeros(Nsamples,num_params); 
tpts_2id = [];
tpts_1id = [];
tpts_0id = [];

count_2id = 0; count_1id =0; count_0id = 0; 
for j = 1:Nsamples
    if ((profiles_all{1}(j,1)>=threshold_fval)&&(profiles_all{1}(j,end)>=threshold_fval)) % p1 identifiable
        if ((profiles_all{2}(j,1)>=threshold_fval)&&(profiles_all{2}(j,end)>=threshold_fval)) %p2 identifiable
            count_2id = count_2id+1;
            practically_identifiable(j,1) = 1; % practically identifiable  
            practically_identifiable(j,2) = 1; % practically identifiable 
            for i = 1:size_sample
                % evaluate t_day at appropriate index, and divide by 24 to convert from hrs to days
                tpts_2id(count_2id,i) = t_day(t_sample(j,i)); 
            end
        else % p2 not identifiable
            count_1id = count_1id+1;
            practically_identifiable(j,1) = 1; % practically identifiable  
            practically_identifiable(j,2) = -1; % not practically identifiable 
            for i = 1:size_sample
                tpts_1id(count_1id,i) = t_day(t_sample(j,i)); 
            end
        end
    else %p1 not identifiable
        if ((profiles_all{2}(j,1)>=threshold_fval)&&(profiles_all{2}(j,end)>=threshold_fval)) %p2 identifiable
            count_1id = count_1id+1;
            practically_identifiable(j,1) = -1; % not practically identifiable  
            practically_identifiable(j,2) = 1; % practically identifiable 
            for i = 1:size_sample
                tpts_1id(count_1id,i) = t_day(t_sample(j,i)); 
            end
        else % p2 not identifiable
            count_0id = count_0id+1;
            practically_identifiable(j,1) = -1; % not practically identifiable  
            practically_identifiable(j,2) = -1; % not practically identifiable 
            for i = 1:size_sample
                tpts_0id(count_0id,i) = t_day(t_sample(j,i)); 
            end
        end
    end
end
num_identifiable = size(tpts_2id,1);
fprintf('Number of cases where parameters remain practically identifiable: %d\n',num_identifiable);

%% Plot all cases where parameters are practically identifiable
if(size(tpts_2id,1)>0) 
    cmap = parula(num_identifiable);
    Ncol = ceil(num_identifiable/14); % at most 12 per column
    for i = 1:num_params
        count = 1;
        figure; 
        set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.45, 0.55]);
        hold on;
        for j = 1:Nsamples
            if (practically_identifiable(j,1) == 1) && (practically_identifiable(j,2) == 1)
                if(size_sample==2)
                    h1 = plot(params(i,:),profiles_all{i}(j,:),'-','LineWidth',2,...
                        'Color', cmap(count,:),'DisplayName',...
                        sprintf('%d&%dd',tpts_2id(count,1),tpts_2id(count,2)));
                elseif(size_sample==3)
                    h1 = plot(params(i,:),profiles_all{i}(j,:),'-','LineWidth',2,...
                        'Color', cmap(count,:),'DisplayName',...
                        sprintf('%d,%d,%dd',tpts_2id(count,1),tpts_2id(count,2),tpts_2id(count,3)));
                else
                    h1 = plot(params(i,:),profiles_all{i}(j,:),'-','LineWidth',2,...
                        'Color', cmap(count,:)); 
                end
                h1.Color(4) = 0.5;
                count = count+1; 
            end
        end
        h2 = plot(params_allData(i,:),profiles(i,:),'-k','LineWidth',3,...
            'DisplayName','All days');
        threshold_plot = threshold_fval*ones(size(params(i,:))); 
        plot(params(i,:),threshold_plot,'--k','LineWidth',2,'DisplayName','CI');
        hold off;
        ylabel('Cost function','FontSize',16);
        legend('FontSize',16,'Location','EastOutside','NumColumns',Ncol);
        %legend('FontSize',16,'Location','NorthEast','NumColumns',Ncol);
        title(['Practically Identifiable Profiles Using ' num2str(size_sample) ...
                '/' num2str(tf-1) ' Time Pts'],'FontSize',16);      
        if i==1 
            xlabel('k_o_n_T','FontSize',16);
            fname_fig = [path 'profile_konT_identifiable'];
        elseif i==2
            xlabel('k_s_y_n_t','FontSize',16);
            fname_fig = [path 'profile_ksynt_identifiable'];
        end
        saveas(gcf,[fname_fig,'.fig']);
        saveas(gcf,[fname_fig,'.png']);
    end
end

%% Scatter plot of times sampled for which parameters are practically identifiable
if(size_sample==2)
    figure; hold on;
    if(size(tpts_2id,1)>0) 
        scatter(tpts_2id(:,1),tpts_2id(:,2),60,"filled",'DisplayName','Both identifiable');
    end
    if(size(tpts_1id,1)>0) 
        scatter(tpts_1id(:,1),tpts_1id(:,2),60,'d',"filled",'DisplayName','One identifiable');
    end
    if(size(tpts_0id,1)>0) 
        scatter(tpts_0id(:,1),tpts_0id(:,2),60,'s',"filled",'DisplayName','None identifiable'); 
    end
    hold off;
    xlabel('First sampled time point','FontSize',16);
    ylabel('Second sampled time point','FontSize',16);
    xlim([0 max(t_sample,[],'all')]); 
    ylim([0 max(t_sample,[],'all')]); 
    title([num2str(num_identifiable) ' unique practically identifiable samplings of ' ...
        num2str(Nsamples)],'FontSize',16);
    legend('FontSize',16,'Location','SouthEast'); 
    fname_fig = [path 'identifiable_time_samples'];
    saveas(gcf,[fname_fig,'.fig']);
    saveas(gcf,[fname_fig,'.png']);
else
    tpts_0idA = [tpts_0id zeros(size(tpts_0id,1),1)];
    tpts_1idA = [tpts_1id ones(size(tpts_1id,1),1)];
    tpts_2idA = [tpts_2id 2*ones(size(tpts_2id,1),1)];
    tpts_all = [tpts_0idA; tpts_1idA; tpts_2idA];
    Y = [];
    for i = 1:size_sample
        Y(:,i) = tpts_all(:,i);
    end
    labeledIdent = categorical(tpts_all(:,size_sample+1),[2 1 0],{'Two identifiable',...
        'One identifiable','None identifiable'}); 
    group = {labeledIdent};
    colors = lines(3);
    ynames = {};
    for i = 1:size_sample
        temp_label = ['t_' num2str(i)]; 
        ynames{i}=temp_label; 
    end
    fig1 = figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]); % fit to screen parameters
    [h,ax,bigAx] = gplotmatrix(Y,[],group,colors,[],20,[],'stairs',ynames);
    set(ax,'FontSize',14);
    title([num2str(num_identifiable) ' unique practically identifiable samplings of ' num2str(Nsamples)],...
        'FontSize',16,'FontWeight','bold');
    h1 = findobj('Tag','legend');
    set(h1, 'Location','NorthWest')
    fname_fig = [path 'identifiable_time_samples'];
    saveas(gcf,[fname_fig,'.fig']);
    saveas(gcf,[fname_fig,'.png']);
    
    % If sample_size = 3, look at spacing between points as well
    if size_sample==3
        if num_identifiable ~= Nsamples % have some non-identifiable samples
            spacing_2id(:,1) = tpts_2id(:,2)-tpts_2id(:,1); 
            spacing_2id(:,2) = tpts_2id(:,3)-tpts_2id(:,2); 
            spacing_1id(:,1) = tpts_1id(:,2)-tpts_1id(:,1); 
            spacing_1id(:,2) = tpts_1id(:,3)-tpts_1id(:,2); 
            spacing_0id(:,1) = tpts_0id(:,2)-tpts_0id(:,1); 
            spacing_0id(:,2) = tpts_0id(:,3)-tpts_0id(:,2); 
            figure; hold on;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); 
            if(size(spacing_2id,1)>0) 
                scatter(spacing_2id(:,1),spacing_2id(:,2),60,"filled",'DisplayName','Both identifiable');
            end
            if(size(spacing_1id,1)>0) 
                scatter(spacing_1id(:,1),spacing_1id(:,2),60,'d',"filled",'DisplayName','One identifiable');
            end
            if(size(spacing_0id,1)>0) 
                scatter(spacing_0id(:,1),spacing_0id(:,2),60,'s',"filled",'DisplayName','None identifiable'); 
            end
            hold off;
            xlabel('Spacing between first two sampled time points','FontSize',16);
            ylabel('Spacing between last two sampled time point','FontSize',16);
            xlim([0 max(t_sample,[],'all')]); 
            ylim([0 max(t_sample,[],'all')]); 
            title([num2str(num_identifiable) ' unique practically identifiable samplings of ' ...
                num2str(Nsamples)],'FontSize',16);
            legend('FontSize',16,'Location','NorthEast'); 
            fname_fig = [path 'identifiable_time_samples_spacing'];
            saveas(gcf,[fname_fig,'.fig']);
            saveas(gcf,[fname_fig,'.png']);
        end
    end
end
if(size_sample==3)
    num_ident = size(tpts_2id(:,1)',2); % tpts_2id because 2 parameters fit here  
                                    % and want case when all fit params are 
                                    % identifiable
    tbuild = [];
    if num_identifiable ~= Nsamples % have some non-identifiable samples
        for i = 1:size_sample
            tbuild(i,:) = [tpts_0id(:,i)' tpts_1id(:,i)' tpts_2id(:,i)'];
        end
    else % all identifiable
        for i = 1:size_sample
            tbuild(i,:) = [tpts_2id(:,i)'];
        end
    end
    t = [];
    for i = 1:size_sample
        t = [t tbuild(i,:)'];
    end
    num_not_ident = size(t,1) - num_ident;
    % classify as 0 if both not identifiable. 1 if both identifiable
    classif = [zeros(1,num_not_ident) ones(1,num_ident)]; 
    
    % Separate data for plotting
    indNotIdent = find(classif==0);
    tNotIdent = t(indNotIdent,:);
    classifNotIdent = classif(indNotIdent); 
    sNotIdent = 50*ones(size(classifNotIdent));

    indIdent = find(classif==1); 
    tIdent = t(indIdent,:);
    classifIdent = classif(indIdent); 
    sIdent = 50*ones(size(classifIdent));

    
    figure;
    scatter3(tNotIdent(:,1),tNotIdent(:,2),tNotIdent(:,3),sNotIdent,classifNotIdent,'r','filled'); % non-identifiable training data
    hold on
    scatter3(tIdent(:,1),tIdent(:,2),tIdent(:,3),sIdent,classifIdent,'b','filled'); % identifiable training data
    legend('Not all identifiable','All identifiable','Location',...
        'northeast','FontSize',16); 
    hold off
    xlabel('t_1','FontSize',16);
    ylabel('t_2','FontSize',16);
    zlabel('t_3','FontSize',16);
    %zlim([0 max(y,[],'all')]); 
    title([num2str(num_ident) ' unique practically identifiable samplings of ' ...
        num2str(size(t,1))],'FontSize',16);
    fname_fig = [path '/identifiable_time_samples_3D'];
    saveas(gcf,[fname_fig,'.fig']);
    saveas(gcf,[fname_fig,'.png']);
    
    %% Now narrow in on the region where most practically identifiable 
    %% samplings appear: t1 = 4 days
    find_t1eq4 = find(tpts_2id(:,1)==4);
    fprintf('Of %d samplings that are practically identifiable, %d of them have t_1 = 4\n',...
        size(tpts_2id,1),size(find_t1eq4,1)); 
    fprintf('The second sampling always occurs between days %d and %d\n',...
        min(tpts_2id(:,2)),max(tpts_2id(:,2))); 
    fprintf('The third sampling always occurs between days %d and %d\n',...
        min(tpts_2id(:,3)),max(tpts_2id(:,3))); 

    % Count number of samplings where both parameters are practically 
    % identifiability within region: t1 = 4, 5<=t2<=10, t3>=16
    count = 0;
    for i = 1:size(tpts_2id(:,1),1)
       if tpts_2id(i,1) == 4
           if (tpts_2id(i,2)>=5)&&(tpts_2id(i,2)<=10)
               if tpts_2id(i,3)>=16
                   count = count+1; 
               end
           end
       end
    end
    fprintf('Of %d samplings that are practically identifiable, %d of them fall in defined region\n',...
        size(tpts_2id,1),count);  
end

fname_out = [path 'analyze_sampling.mat'];
save(fname_out,'tpts_2id','tpts_1id','tpts_0id');