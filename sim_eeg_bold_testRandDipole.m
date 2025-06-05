% simMEG [nE x Time]
% simBOLD [1 x Time]
function [simEEG, simBOLD,x] = sim_eeg_bold_testRandDipole(PSP,time,nE,nPV,nSOM,nVIP,nP,C,dt)

    % ----- simulate EEG ------
    dipEtoE   = rand*2-1;%-1;    % E -> E
    dipPVtoE  = rand*2-1;%-0.1;  % PV -> E
    dipSOMtoE = rand*2-1;% 1;    % SOM -> E
    dipVIPtoE = rand*2-1;% 1;    % VIP -> E 
    dipTHtoE  = rand*2-1;% 0.1;  % Th -> E
    
    % dipoles [nE x (nP+1)]
    dipoles=[ones(nE,nE)*dipEtoE, ...
             ones(nE,nPV)*dipPVtoE, ...
             ones(nE,nSOM)*dipSOMtoE,...
             ones(nE,nVIP)*dipVIPtoE,...
             ones(nE,1)*dipTHtoE];   % [nE x (nP+1)]
    
    % weighted by cell count
    cellcount=repmat(C(1:nE)',[1,nP+1]);
    dipoles=dipoles.*cellcount; % dipole weighted by cell count of E pops
    
    % simEEG [nE x Time]
    simEEG=zeros(nE,length(time));
    for i=1:nE
        simEEG(i,:)=dipoles(i,:)*abs(squeeze(PSP(i,:,:))); %[1x18]*[18xTime]
    end
   
    % ----- simulate BOLD ------
    % consider all EPSPs
    EPSP = squeeze(sum(PSP(:,[1:nE,nP+1],:),2)); % [nP x Time] EPSP from E and Ext
    x = C*EPSP; % weighted by cell count
    simBOLD=sim_bold(x,0,dt);  

end