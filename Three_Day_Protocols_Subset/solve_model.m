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
