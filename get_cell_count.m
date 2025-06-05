% Cell counts (PV,SOM,VIP)
% [ref] Gonchar, Y., Wang, Q., & Burkhalter, A. H. (2008). 
%       Multiple distinct subtypes of GABAergic neurons in mouse visual 
%       cortex identified by triple immunostaining. Frontiers in neuroanatomy.
% Cell counts (E)
% [ref] Fang, R., Xia, C., Close, J. L., Zhang, M., He, J., Huang, 
%       Z., ... & Zhuang, X. (2022). Conservation and divergence of cortical 
%       cell organization in human and mouse revealed by MERFISH. Science.
%
% Ce: Cell count of E   (L2/3, L4, L5, L6)
% Cp: Cell count of PV  (L2/3, L4, L5, L6)
% Cs: Cell count of SST (L2/3, L4, L5, L6)
% Cv: Cell count of VIP (L2/3, L4, L5, L6)
%
function [Ce,Cp,Cs,Cv]=get_cell_count(count,plotOn)

    if nargin<1
        count=1000;   % set total cell count
    end

    if nargin<2
        plotOn=0;
    end

    % [1] Fang, R., Xia, C., Close, J. L., Zhang, M., He, J., Huang, Z., ... & Zhuang, X. (2022). 
    %     Conservation and divergence of cortical cell organization in human 
    %     and mouse revealed by MERFISH. Science, 377(6601), 56-62.
    %
    proportion = [10,7,4,3]/24; % PV, SOM, VIP, LAMP5 (mouse visual) (Sum = 100%) 
    EIratio=[6.5, 7.5, 4, 7];   % L2/3, L4, L5, L6 (mouse visual)


    
    % [2] Gonchar, Y., Wang, Q., and Burkhalter, A. (2007). Multiple distinct 
    %     subtypes of GABAergic neurons in mouse visual cortex identified by 
    %     triple immunostaining. Front. Neuroanat. 1:3. doi: 10.3389/neuro.05.003.2007
    %
    PV_distribition=[.32,.18,.40,.22]; % L2/3, L4, L5, L6 (mouse visual) (Sum > 1 !!!)
    PV_distribition=PV_distribition/sum(PV_distribition); % re-size (Sum = 100%)
    SOM_distribition=[.31,.14,.30,.20];% L2/3, L4, L5, L6 (mouse visual) (Sum = 95%)
    VIP_distribition=[.57,.13,.11,.12];% L2/3, L4, L5, L6 (mouse visual) (Sum = 93%)
    LAMP5_distribition=[.25,.25,.25,.25];% L2/3, L4, L5, L6 LAMP5 (!!!assuming uniform because lack of ref)
    
    
    
    %-----------------------------------------------
    I_count_all=100;   % assume I={PV,SOM,VIP,LAMP5}

    PV_count_all=I_count_all*proportion(1);
    SOM_count_all=I_count_all*proportion(2);
    VIP_count_all=I_count_all*proportion(3);
    LAMP5_count_all=I_count_all*proportion(4);
    
    PV_count_layers = PV_count_all*PV_distribition;          % L2/3, L4, L5, L6
    SOM_count_layers = SOM_count_all*SOM_distribition;       % L2/3, L4, L5, L6
    VIP_count_layers = VIP_count_all*VIP_distribition;       % L2/3, L4, L5, L6
    LAMP5_count_layers = LAMP5_count_all*LAMP5_distribition; % L2/3, L4, L5, L6
    E_count_layers = (PV_count_layers+SOM_count_layers+...
                      VIP_count_layers+LAMP5_count_layers).*EIratio;
    
    E_count_all=sum(E_count_layers);
    %---------adjust cell count-------------
    adj=count/(E_count_all+I_count_all);

%     E_count_all=E_count_all*adj;
%     I_count_all=I_count_all*adj;
%     PV_count_all=PV_count_all*adj;
%     SOM_count_all=SOM_count_all*adj;
%     VIP_count_all=VIP_count_all*adj;
%     LAMP5_count_all=LAMP5_count_all*adj;

    E_count_layers=E_count_layers*adj;
    PV_count_layers=PV_count_layers*adj;
    SOM_count_layers=SOM_count_layers*adj;
    VIP_count_layers=VIP_count_layers*adj;
    LAMP5_count_layers=LAMP5_count_layers*adj;

    %----------------------------------    
%     PV_count_pops = [PV_count_layers(1)+PV_count_layers(2), PV_count_layers(3)+PV_count_layers(4)]; %PV1(L2/3/4), PV2(L5/6)
%     SOM_count_pops = [SOM_count_layers(1)+SOM_count_layers(2), SOM_count_layers(3)+SOM_count_layers(4)]; %SOM1(L2/3/4), SOM2(L5/6)
%     VIP_count_pops = [VIP_count_layers(1)+VIP_count_layers(2), VIP_count_layers(3)+VIP_count_layers(4)]; %VIP1(L2/3/4), VIP2(L5/6)
%     E_count_pops = [E_count_layers(1), E_count_layers(2), E_count_layers(3)+E_count_layers(4)]; %E1(L2/3), E3(L4), E2(L5/6)
    
    %-------------output-------------- 
%     Ce=E_count_pops;    % E1(L2/3), E3(L4), E2(L5/6);
%     Cp=PV_count_pops;   % PV1(L2/3/4), PV2(L5/6);
%     Cs=SOM_count_pops;  % SOM1(L2/3/4), SOM2(L5/6);
%     Cv=VIP_count_pops;  % VIP1(L2/3/4), VIP2(L5/6);

    Ce=E_count_layers;  % L2/3, L4, L5, L6
    Cp=PV_count_layers; % L2/3, L4, L5, L6
    Cs=SOM_count_layers; % L2/3, L4, L5, L6
    Cv=VIP_count_layers; % L2/3, L4, L5, L6

    if plotOn
        %-----print-----
        fprintf('Ce:[%g, %g, %g, %g]\n',Ce)
        fprintf('Cp:[%g, %g, %g, %g]\n',Cp)
        fprintf('Cs:[%g, %g, %g, %g]\n',Cs)
        fprintf('Cv:[%g, %g, %g, %g]\n',Cv)
        
        %-----plot-----
        figure;
        subplot(221); pie(proportion); title('Inh. cell proportion [1]');legend({'PV', 'SOM', 'VIP', 'LAMP5'},'location','northeastoutside')
        subplot(223); bar(EIratio); xlim([0.5,4.5]); set(gca,'xtick',1:4,'xticklabel',{'L2/3','L4','L5','L6'});title('E/I ratio [1]');
        %subplot(222);bar([PV_count_layers',SOM_count_layers',VIP_count_layers',LAMP5_count_layers']);
        subplot(222);bar([PV_distribition',SOM_distribition',VIP_distribition',LAMP5_distribition']*100);
                      legend({'PV', 'SOM', 'VIP', 'LAMP5(?)'},'location','eastoutside');set(gca,'xtick',1:4,'xticklabel',{'L2/3','L4','L5','L6'});title('Inh. cell dist. [2]'); ylabel('[%]')
        subplot(224);bar([Cp',Cs',Cv',LAMP5_count_layers',Ce']);
                      legend({'PV', 'SOM', 'VIP', 'LAMP5','E'},'location','eastoutside');set(gca,'xtick',1:4,'xticklabel',{'L2/3','L4','L5','L6'});
                      title(sprintf('Assume total cell count = %d cells', count));
    end
end