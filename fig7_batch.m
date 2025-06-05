% >> fig7_batch(); % run simulations and save result.
% >> fig7();       % plot result.
%
function fig7_batch()

    % input combinations
    rangeA=[5,6,7,8,9]/2;
    rangeB=(0:50)/2;
    rangeC=(0:10)/2;
    [As,Bs,Cs] =ndgrid(rangeA,rangeB,rangeC);
    all_spect  =zeros(numel(As),2501); 
    all_freq   =zeros(numel(As),2501); 
    all_pop_spect = zeros(numel(As),17,2501);
    all_psp_spect = zeros(numel(As),17,2501);
    all_rate_mean = zeros(numel(As),2501);
    all_rate_std = zeros(numel(As),2501);

    baseline=simulation(0,0,0);

    %parfor i=1:numel(As)
    for i=1:numel(As)  

        fprintf('%d/%d...',i,numel(As))
        r=simulation(As(i),Bs(i),Cs(i));

        % -----record-----
        all_spect(i,:)=r.spect'; 
        all_freq(i,:)=r.freq';
        all_pop_spect(i,:,:)=r.spect_pop'; % spect of 17 rates
        all_psp_spect(i,:,:)=r.spect_psp'; % spect of 17->E PSPs
        all_rate_mean(i,:)=r.rate_mean;    % mean of 17 rates
        all_rate_std(i,:)=r.rate_std;      % std of 17 rates
  
        fprintf('done.\n')
    end
    freq=all_freq(1,:);
    
    filename=fullfile('data',sprintf('fig7_batch_result_%s.mat',datestr(now,'yyyymmdd-HHMM')));
    save(filename,'As','Bs','Cs','rangeA','rangeB','rangeC',...
         'freq','all_spect',...
         'all_pop_spect','all_psp_spect',...
         'all_rate_mean','all_rate_std','baseline');    

    fprintf('%s saved.\n',filename)

end

%==========================================================================
function r=simulation(a,b,c)

    duration = 600;   % sec per run
    dt       = 1e-3;  % simulation time step (sec)
    time_discard=[30 duration];  % sec

    % Input weights (to 17 populations)
    Ilat = [1; %to E2/3
            0; %to E4
            0; %to E5ET
            0; %to E5IT
            0; %to E6
            0; %to PV2/3 
            0; %to PV4
            0; %to PV5
            0; %to PV6
            1; %to SOM2/3
            0; %to SOM4
            0; %to SOM5
            0; %to SOM6
            0; %to VIP2/3 
            0; %to VIP4    
            0; %to VIP5 
            0];%to VIP6 
    
    Imod = [0; %to E2/3
            0; %to E4
            0; %to E5ET
            0; %to E5IT
            0; %to E6
            0; %to PV2/3 
            0; %to PV4
            0; %to PV5
            0; %to PV6
            0; %to SOM2/3
            0; %to SOM4
            0; %to SOM5
            0; %to SOM6
            1; %to VIP2/3 
            0; %to VIP4    
            0; %to VIP5 
            0];%to VIP6 
    
    Itha = [0; %to E2/3
            1; %to E4
            0; %to E5ET
            0; %to E5IT
            0; %to E6
            0; %to PV2/3 
            1; %to PV4
            0; %to PV5
            0; %to PV6
            0; %to SOM2/3
            0; %to SOM4
            0; %to SOM5
            0; %to SOM6
            0; %to VIP2/3 
            0; %to VIP4    
            0; %to VIP5 
            0];%to VIP6 


    noiseLevel=1;    
    Iext  = (Ilat*a + Imod*b + Itha*c)*ones(1,duration/dt);  

    % Connectivity W
    [W,C,~,nP,nE,nPV,nSOM,nVIP]=get_W17(0);
    g  = 20*1e3; 
    W = [W*g, ones(nP,1)*1e3];
    
    
    % Synaptic time constants
    tau = get_tau(nP,nE,nPV,nSOM,nVIP);
    
    % Sigmoid functions
    sigmParam = get_sigmParam(nE,nPV,nSOM,nVIP);
    
    % simulation
    [rate, time, PSP] = model_jr(W, tau, sigmParam, Iext, noiseLevel, dt);

    % simEEG and simBOLD
    [simEEG, simBOLD, x] = sim_eeg_bold(PSP,time,nE,nPV,nSOM,nVIP,nP,C,dt);
     simEEG = sum(simEEG,1);
  
    %-----EEG power spectrum-----
    tidx=find(time>time_discard(1) & time<=time_discard(2)); 
    windowLength = 5000; 
    data = simEEG(tidx)-mean(simEEG(tidx));
    [spect,freq] = compute_spectrum(data,1/dt,windowLength); 
    %----- rate mean & std
    rate_mean=mean(rate(:,tidx),2);  % mean of 17 rates
    rate_std=std(rate(:,tidx),[],2); % std of 17 rates
    %-----firing rate power spectrum-----
    data = rate(:,tidx)-repmat(rate_mean,[1 length(tidx)]);
    [spect_pop,~]=compute_spectrum(data',1/dt,windowLength); 
    %-----PSP power spectrum-----
    tmpPSP=squeeze(mean(PSP(1:5,1:17,:),1));
    data = tmpPSP(:,tidx)-repmat(mean(tmpPSP(:,tidx),2),[1 length(tidx)]);
    [spect_psp,~]=compute_spectrum(data',1/dt,windowLength); 

    % -----record-----
    r.freq  = freq;
    r.spect = spect; % EEG
    r.spect_pop = spect_pop; % rate
    r.spect_psp = spect_psp; % PSP
    r.rate_mean = rate_mean;
    r.rate_std  = rate_std;
end