% Sigmoid function derided from mean F-V curves 
% From the allen brain database (https://celltypes.brain-map.org/data)

function [sigmE,sigmPV,sigmSOM,sigmVIP] = get_sigmoid(speedOn,plotOn)

    if speedOn
        %-----------E---------------
        sigmE.a=0.1278;
        sigmE.b=32.11;
        sigmE.c=31.4;
        %-----------PV--------------
        sigmPV.a=0.1422;
        sigmPV.b=40.03;
        sigmPV.c=166.8;
        %-----------SOM-------------
        sigmSOM.a=0.07937;
        sigmSOM.b=42.01;
        sigmSOM.c=56.95;
        %-----------VIP-------------
        sigmVIP.a=0.07041;
        sigmVIP.b=37.86;
        sigmVIP.c=38.53;
        return    
    end

    %-----read file-----
    T = readtable(fullfile('data','cell_types_specimen_details.csv'));
    id=T.specimen__id;
    line=T.line_name;
    location=T.structure__name;
    locationShort=T.structure__acronym;
    threshold=T.ef__threshold_i_long_square;
    slope=T.ef__f_i_curve_slope;  % Hz/pA
    avgRate=T.ef__avg_firing_rate;
    resistance=T.ef__ri;
    

    %-----------E---------------
    idx=find((contains(line,'Cux2-CreERT2') | ...   %L23
              contains(line,'Rbp4-Cre_KL100') | ... %L5
              contains(line,'Rorb-IRES2-Cre') | ... %L4
              contains(line,'Scnn1a-Tg2-Cre') | ... %L4, L5
              contains(line,'Scnn1a-Tg3-Cre') | ... %L4, L5
              contains(line,'Tlx3-Cre_PL56') | ...  %L5 IT
              contains(line,'Sim1-Cre_KJ18') | ...  %L5 PT
              contains(line,'Ntsr1-Cre_GN220')) ... %L6
              & contains(locationShort,'VIS'));
    E.id=id(idx);
    E.line=line(idx);
    E.location=location(idx);
    E.locationShort=locationShort(idx);
    E.threshold=threshold(idx);
    E.slope=slope(idx);
    E.avgRate=avgRate(idx);
    E.resistance=resistance(idx);    
    [meanFVcurves_e,V,fit_e,sigmE]=get_curve(E,150,plotOn);  %maxf: 100-200Hz

    %------------PV--------------
    idx=find(contains(line,'Pvalb-IRES-Cre') & contains(locationShort,'VIS'));
    PV.id=id(idx);
    PV.line=line(idx);
    PV.location=location(idx);
    PV.locationShort=locationShort(idx);
    PV.threshold=threshold(idx);
    PV.slope=slope(idx);
    PV.avgRate=avgRate(idx);
    PV.resistance=resistance(idx);    
    [meanFVcurves_pv,V,fit_pv,sigmPV]=get_curve(PV,250,plotOn); %maxf: 200-300Hz

    %------------SST--------------
    idx=find(contains(line,'Sst-IRES-Cre') & contains(locationShort,'VIS'));
    SST.id=id(idx);
    SST.line=line(idx);
    SST.location=location(idx);
    SST.locationShort=locationShort(idx);
    SST.threshold=threshold(idx);
    SST.slope=slope(idx);
    SST.avgRate=avgRate(idx);
    SST.resistance=resistance(idx);    
    [meanFVcurves_sst,V,fit_sst,sigmSOM]=get_curve(SST,100,plotOn); %maxf: 50-150Hz
    
    %-----------VIP---------------
    idx=find(contains(line,'Vip-IRES-Cre') & contains(locationShort,'VIS'));
    VIP.id=id(idx);
    VIP.line=line(idx);
    VIP.location=location(idx);
    VIP.locationShort=locationShort(idx);
    VIP.threshold=threshold(idx);
    VIP.slope=slope(idx);
    VIP.avgRate=avgRate(idx);
    VIP.resistance=resistance(idx);    
    [meanFVcurves_vip,V,fit_vip,sigmVIP]=get_curve(VIP,150,plotOn); %maxf: 100-200Hz

    if plotOn
        %----------------
        figure;
        plot(V,[meanFVcurves_e,meanFVcurves_pv,meanFVcurves_sst,meanFVcurves_vip],'linewidth',2); hold on
        plot(V,[fit_e,fit_pv,fit_sst,fit_vip],'k--','linewidth',1);
        xlabel('mV'); ylabel('Hz')
        legend('E','PV','SST','VIP','fit')
        set(gcf,'position',[0 0 300 300])
    
        %----------------
        figure;
        subplot(211);
        plot(V,[meanFVcurves_e,meanFVcurves_pv,meanFVcurves_sst,meanFVcurves_vip],'linewidth',2); xlim([0 35]);ylim([0 35]);
        legend('E','PV','SST','VIP')
        subplot(212);
        plot(V,[fit_e,fit_pv,fit_sst,fit_vip],'--','linewidth',2); xlim([0 35]);ylim([0 35]);
        legend('E-fit','PV-fit','SST-fit','VIP-fit')
        
        set(gcf,'position',[0 0 300 600])
    
        %----------------
        figure;
        scatter(E.threshold.*E.resistance/1e3,E.slope./E.resistance*1e3,20,'filled');hold on;
        scatter(PV.threshold.*PV.resistance/1e3,PV.slope./PV.resistance*1e3,20,'filled');
        scatter(SST.threshold.*SST.resistance/1e3,SST.slope./SST.resistance*1e3,20,'filled');
        scatter(VIP.threshold.*VIP.resistance/1e3,VIP.slope./VIP.resistance*1e3,20,'filled');
        xlabel('threshold (mV)');
        ylabel('slope (Hz/mV)');
        legend('E','PV','SST','VIP')
        set(gcf,'position',[0 0 300 300])
    end

end
%==========================================================================
function [meanFVcurves,V,fitcurve,myfit]=get_curve(celltype,maxf,plotOn)

    avg2max=2;%2
    V=(0:1:200)'; % PSP (mV)
    FVcurves=zeros(length(V),length(celltype.id));
    for i=1:length(celltype.id)
       I=V/celltype.resistance(i)*1e3; % I(pA)
       tmp=celltype.slope(i)*(I-celltype.threshold(i));
       tmp=max(tmp,0);
       tmp=min(tmp,celltype.avgRate(i)*avg2max);  % first cut
       %tmp=min(tmp,maxf);                         % second cut
       FVcurves(:,i)=tmp;
    end
    FVcurves(:,isnan(celltype.avgRate))=[]; % remove, because no avgRate
    meanFVcurves=mean(FVcurves,2);
    %----- fit------
    myfittype = fittype("c./(1+exp(-a*(x-b)))",...
                   dependent="y",independent="x",...
                   coefficients=["a" "b" "c"]);
    [myfit, gof] = fit(V,meanFVcurves,myfittype,'StartPoint',[1 50 100],'Lower',[0 0 0],'upper',[3, 100 300]);
    a=myfit.a;
    b=myfit.b;
    c=myfit.c;
    fitcurve=c./(1+exp(-a*(V-b)));
    %----------------
    if plotOn
        figure;
        subplot(421);histogram(celltype.threshold); title('threshold');xlabel('pA');
        subplot(423);histogram(celltype.slope); title('slope');xlabel('Hz/pA');
        subplot(425);histogram(celltype.avgRate); title('avgRate');xlabel('Hz');
        subplot(427);histogram(celltype.resistance); title('resistance R');xlabel('Mohm');
        subplot(122);plot(V,FVcurves,'color',[1 1 1]*0.7); hold on;
                     h1=plot(V,meanFVcurves,'k','linewidth',2);
                     h2=plot(V,fitcurve,'r--','linewidth',1);
                     legend([h1;h2],{'mean','fit'})

                    xlabel('[V = I*R] mV');ylabel('Hz')
                    title({sprintf('F-I curve (%d cells) ',length(celltype.id)),...
                           sprintf('fit(th=%g;r=%g;max=%g)',b,a,c)})
    end
    

end