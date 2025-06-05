% Output:
% [to x from]
% Wee {E1,E2,E3a,E3b,E4}    x {E1,E2,E3a,E3b,E4}
% Wpe {PV1,PV2,PV3,PV4}     x {E1,E2,E3a,E3b,E4}
% Wse {SOM1,SOM2,SOM3,SOM4} x {E1,E2,E3a,E3b,E4}
% Wve {VIP1,VIP2,VIP3,VIP4} x {E1,E2,E3a,E3b,E4}
% Wep 
% Wpp 
% Wsp 
% Wvp 
% Wes 
% Wps 
% Wss
% Wvs 
% Wev 
% Wpv 
% Wsv 
% Wvv
% , where
% E1: E in L2/3
% E2: E in L4
% E3a: E in L5 ET neurons
% E3b: E in L5 IT neurons
% E4: E in L6
%
%
% [1] Billeh, Yazan N., et al. "Systematic integration of structural and
%     functional data into multi-scale models of mouse primary visual cortex."
%     Neuron 106.3 (2020): 388-403.
% [2] Campagnola, Luke, et al. "Local connectivity and synaptic dynamics 
%     in mouse and human neocortex." Science 375.6585 (2022): eabj5861.
function [W2,C2,labels,nP,nE,nPV,nSOM,nVIP]=get_W17(plotOn)

    if nargin<1
        plotOn=0;
    end
    labels={'E2/3','E4','E5ET','E5IT','E6',...
            'PV2/3','PV4','PV5','PV6',...
            'SOM2/3','SOM4','SOM5','SOM6',...
            'VIP2/3','VIP4','VIP5','VIP6'};

    % Connection probabilities (to x from) [2]
    %  1    2     3      4     5   6    7   8    9    10  11   12  13   14 15   16   17
    % E2/3,PV2/3,SST2/3,VIP2/3,E4,PV4,SST4,VIP4,E5ET,E5IT,PV5,SST5,VIP5,E6,PV6,SST6,VIP6
    [P,~]=compagnola_conn(0); % 17 x 17   
    % scale P by lateral spread (sigma)
    sigmaEE=127;    % um  (E->E)
    sigmaEI=99.84;  % um  (E->I)
    sigmaIE=96.6;   % um  (I->E)
    sigmaII=126.77; % um  (I->I)
    iE=[1,5,9,10,14];
    iI=[2,3,4,6,7,8,11,12,13,15,16,17];
    sigma=zeros(16);
    sigma(iE,iE)=sigmaEE; % um  (E->E)
    sigma(iI,iE)=sigmaEI; % um  (E->I)
    sigma(iE,iI)=sigmaIE; % um  (I->E)
    sigma(iI,iI)=sigmaII; % um  (I->I)
    sigma=sigma/max(sigma(:));
    P=P.*sigma;

    % Synaptic strengths (to x from) [1] duplicate E5 for E5ET and E5IT
    %  1    2     3      4     5   6    7   8    9    10  11   12  13   14 15   16   17
    % E2/3,PV2/3,SST2/3,VIP2/3,E4,PV4,SST4,VIP4,E5ET,E5IT,PV5,SST5,VIP5,E6,PV6,SST6,VIP6
    S     =[0.36,	1.49,	0.86,	1.31,	0.34,	1.39,	0.69,	0.91,	0.74,0.74,	1.32,	0.53,	0,	0,	0,	0,	0;...
            0.48,	0.68,	0.42,	0.41,	0.56,	0.68,	0.42,	0.41,	0.2,0.2,	0.79,	0,	0,	0,	0,	0,	0;...
            0.31,	0.5,	0.15,	0.52,	0.3,	0.5,	0.15,	0.52,	0.22,0.22,	0,	0,	0,	0,	0,	0,	0;...
            0.28,	0.18,	0.32,	0.37,	0.29,	0.18,	0.32,	0.37,	0,0,	    0,	0,	0,	0,	0,	0,	0;...
            0.78,	1.39,	0.69,	0.91,	0.83,	1.29,	0.51,	0.51,	0.63,0.63,	1.25,	0.52,	0.91,	0.96,	0,	0,	0;...
            0.56,	0.68,	0.42,	0.41,	0.64,	0.68,	0.42,	0.41,	0.73,0.73,	0.94,	0.42,	0.41,	0,	0,	0,	0;...
            0.3,	0.5,	0.15,	0.52,	0.29,	0.5,	0.15,	0.52,	0.28,0.28,	0.45,	0.28,	0.52,	0,	0,	0,	0;...
            0.29,	0.18,	0.32,	0.37,	0.29,	0.18,	0.32,	0.37,	0,0,	    0.18,	0.33,	0.37,	0,	0,	0,	0;...
            0.47,	1.25,	0.52,	0.91,	0.38,	1.25,	0.52,	0.91,	0.75,0.75,	1.2,	0.52,	1.31,	0.4,	2.5,	0.52,	1.31;...
            0.47,	1.25,	0.52,	0.91,	0.38,	1.25,	0.52,	0.91,	0.75,0.75,	1.2,	0.52,	1.31,	0.4,	2.5,	0.52,	1.31;...        
            0,	0.51,	0,	0,	0,	0.94,	0.42,	0.41,	0.81,	1.19,	0.41,0.41,	0.41,	0.81,	1.19,	0.41,	0.41;...
            0.25,	0,	0.39,	0,	0.28,	0.45,	0.28,	0.52,	0.27,	0.4,0.4,	0.4,	0.52,	0.27,	0.4,	0.4,	0.52;...
            0,	0,	0,	0,	0.29,	0.18,	0.33,	0.37,	0.28,	0.18,	0.33,0.33,	0.37,	0.28,	0.18,	0.33,	0.37;...
            0,	0,	0,	0,	0,	0,	0,	0,	0.23,	2.5,	0.52,	1.31,	0.94,0.94,	3.8,	0.52,	1.31;...
            0.81,	0,	0,	0,	0.81,	0,	0,	0,	0.81,	1.19,	0.41,	0.41,0.41,	0.81,	1.19,	0.41,	0.41;...
            0,	0,	0,	0,	0,	0,	0,	0,	0.27,	0.4,	0.4,	0.52,	0.27,0.27,	0.4,	0.4,	0.52;...
            0,	0,	0,	0,	0,	0,	0,	0,	0.28,	0.18,	0.33,	0.37,	0.28,0.28,	0.18,	0.33,	0.37]'; % to x from

    % Cell count 
    % Assume cell count of E5ET and E5IT are equally half of E5
    [Ce,Cp,Cs,Cv]=get_cell_count(1,0); 
    C = [Ce(1),Cp(1),Cs(1),Cv(1), ...
         Ce(2),Cp(2),Cs(2),Cv(2), ...
         Ce(3)/2,Ce(3)/2,Cp(3),Cs(3),Cv(3),...
         Ce(4),Cp(4),Cs(4),Cv(4)];
    
    
    PS=P.*S;                                                      
    W = PS.*repmat(C,[17 1]);
    
   
    nP   = 17;
    nE   = 5;
    nPV  = 4;
    nSOM = 4;
    nVIP = 4;
       
    iE = [1,5,9,10,14]; % E2/3,E4,E5ET,E5IT,E6
    iP = [2,6,11,15];   % PV2/3,PV4,PV5,PV6
    iS = [3,7,12,16];   % SOM2/3,SOM4,SOM5,SOM6
    iV = [4,8,13,17];   % VIP2/3,VIP4,VIP5,VIP6
    
    Wee = W(iE,iE);
    Wpe = W(iP,iE);
    Wse = W(iS,iE);
    Wve = W(iV,iE);
    Wep = W(iE,iP);
    Wpp = W(iP,iP);
    Wsp = W(iS,iP);
    Wvp = W(iV,iP);
    Wes = W(iE,iS);
    Wps = W(iP,iS);
    Wss = W(iS,iS);
    Wvs = W(iV,iS);
    Wev = W(iE,iV);
    Wpv = W(iP,iV);
    Wsv = W(iS,iV);
    Wvv = W(iV,iV);
    
    W2 = [Wee, -Wep, -Wes, -Wev;...
          Wpe, -Wpp, -Wps, -Wpv;...
          Wse, -Wsp, -Wss, -Wsv;...
          Wve, -Wvp, -Wvs, -Wvv];
    
    C2 = C([iE,iP,iS,iV]); 


    % plot connectivity
    if plotOn
        plot_conn(Wee,Wpe,Wse,Wve,Wep,Wpp,Wsp,Wvp,Wes,Wps,Wss,Wvs,Wev,Wpv,Wsv,Wvv); 
        set(gcf,'Name','W','NumberTitle','off')
    end

    % reduce L4E self-excitation (from mouse to human)
    W2(2,2) = W2(2,2)*0; 
    

end
%==========================================================================
function plot_conn(Wee,Wpe,Wse,Wve,Wep,Wpp,Wsp,Wvp,Wes,Wps,Wss,Wvs,Wev,Wpv,Wsv,Wvv)

  xtext={'E2/3','E4','E5ET','E5IT','E6','PV2/3','PV4','PV5','PV6','SOM2/3','SOM4','SOM5','SOM6','VIP2/3','VIP4','VIP5','VIP6'};
    figure
    % E
    ymax=max([Wee(:);Wpe(:);Wse(:);Wve(:)])*1.1;
    ttext={'E2/3','E4','E5ET','E5IT','E6'};
    h=[2,7,11,12,17];
    for i=1:5 
       subplot(4,5,h(i))
       bar([Wee(:,i);Wpe(:,i);Wse(:,i);Wve(:,i)]); ylim([0 ymax]);xlim([0.5,17.5])
       hold on;plot([5.5;9.5;13.5]*[1,1],ylim,'k')
       set(gca,'xtick',1:length(xtext),'xticklabel',xtext)
       title(sprintf('from %s',ttext{i}))
       grid on
    end
    % PV
    ymax=max([Wep(:);Wpp(:);Wsp(:);Wvp(:)])*1.1;
    ttext={'PV2/3','PV4','PV5','PV6'};
    h=[3,8,13,18];
    for i=1:4 
       subplot(4,5,h(i))
       bar([Wep(:,i);Wpp(:,i);Wsp(:,i);Wvp(:,i)]); ylim([0 ymax]);xlim([0.5,17.5])
       hold on;plot([5.5;9.5;13.5]*[1,1],ylim,'k')
       set(gca,'xtick',1:length(xtext),'xticklabel',xtext)
       title(sprintf('from %s',ttext{i}))
       grid on
    end
    % SOM
    ymax=max([Wes(:);Wps(:);Wss(:);Wvs(:)])*1.1;
    ttext={'SOM2/3','SOM4','SOM5','SOM6'};
    h=[4,9,14,19];
    for i=1:4 
       subplot(4,5,h(i))
       bar([Wes(:,i);Wps(:,i);Wss(:,i);Wvs(:,i)]); ylim([0 ymax]);xlim([0.5,17.5])
       hold on;plot([5.5;9.5;13.5]*[1,1],ylim,'k')
       set(gca,'xtick',1:length(xtext),'xticklabel',xtext)
       title(sprintf('from %s',ttext{i}))
       grid on
    end
    % VIP
    ymax=max([Wev(:);Wpv(:);Wsv(:);Wvv(:)])*1.1;
    ttext={'VIP2/3','VIP4','VIP5','VIP6'};
    h=[5,10,15,20];
    for i=1:4 
       subplot(4,5,h(i))
       bar([Wev(:,i);Wpv(:,i);Wsv(:,i);Wvv(:,i)]); ylim([0 ymax]);xlim([0.5,17.5])
       hold on;plot([5.5;9.5;13.5]*[1,1],ylim,'k')
       set(gca,'xtick',1:length(xtext),'xticklabel',xtext)
       title(sprintf('from %s',ttext{i}))
       grid on
    end
end