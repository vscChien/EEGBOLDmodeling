% simEEG [nE x Time]
% simBOLD [1 x Time]
function [simEEG, simBOLD, x] = sim_eeg_bold_opticmode(PSP,time,nE,nPV,nSOM,nVIP,nP,C,dt,OpticMode,dipConf)

    if nargin<10
        OpticMode=0; % Optic stimulation (i.e., Iext is not included in the BOLD simulation)
    end
    % ----- simulate EEG ------
    if nargin<11
        dipEtoE   = -1;    % E -> E
        dipPVtoE  = -0.1;  % PV -> E
        dipSOMtoE =  1;    % SOM -> E
        dipVIPtoE =  1;    % VIP -> E  
        dipTHtoE  =  0.1;  % Th -> E
    else
        dipEtoE   = dipConf(1);  % E -> E
        dipPVtoE  = dipConf(2);  % PV -> E
        dipSOMtoE = dipConf(3);  % SOM -> E
        dipVIPtoE = dipConf(4);  % VIP -> E  
        dipTHtoE  = dipConf(5);  % Th -> E        
    end
    

    % dipole configuration [nE x (nP+1)]
    dipoles=[ones(nE,nE)*dipEtoE, ...
             ones(nE,nPV)*dipPVtoE, ...
             ones(nE,nSOM)*dipSOMtoE,...
             ones(nE,nVIP)*dipVIPtoE,...
             ones(nE,1)*dipTHtoE];   
    
    % weighted by cell count of E pops
    cellcount=repmat(C(1:nE)',[1,nP+1]);
    dipoles=dipoles.*cellcount; 
    
    % simEEG [nE x Time]
    simEEG=zeros(nE,length(time));
    for i=1:nE
        simEEG(i,:)=dipoles(i,:)*abs(squeeze(PSP(i,:,:))); %[1x(nP+1)]*[(nP+1)xTime]
    end

    % ----- simulate BOLD ------
    % consider all EPSPs
    if OpticMode
        EPSP = squeeze(sum(PSP(:,1:nE,:),2)); % [nP x Time] EPSP from E and Ext
        %disp('OpticaMode: Iext is excluded from the BOLD simulation.')
    else
        EPSP = squeeze(sum(PSP(:,[1:nE,nP+1],:),2)); % [nP x Time] EPSP from E and Ext
    end
    x = C*EPSP; % weighted by cell count
    simBOLD=sim_bold(x,0,dt);  
 
end