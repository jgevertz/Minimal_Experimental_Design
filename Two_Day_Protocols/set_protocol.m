function [tDose,cDose,ti] = set_protocol(intvl,ivdose,dosenum,T)
    tDose = zeros(1,dosenum); cDose = zeros(1,dosenum);
    for i=1:dosenum
        tDose(i)=intvl*(i-1);
        cDose(i)=ivdose;
    end
    ti = [tDose max(T)];
end