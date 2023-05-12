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
    
    %% Now pull out data at only time points in t_sample
    % Goodness of fit: divide by 2 to match likelihood formula
    SSE_divideby_Var = 0.5*sum(((TO_tme_avg(t_sample)-TO_tme(t_sample)).^2)./(TO_tme_std(t_sample).^2));
    if flag == 1
        fprintf('\tAnd an error of %f\n', SSE_divideby_Var); 
    end
end