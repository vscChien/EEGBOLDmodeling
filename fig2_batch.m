% >> fig2_batch(); % run simulations and save result. (will take a couple hours)
% >> fig2();       % plot fig.2
%
function fig2_batch()

    duration = 600;   % sec per run
    dt       = 1e-3;  % simulation time step (sec)
    dt2      = 1;     % downsampled time step (sec)
    time_discard=30;  % only take stable bold signals (t>30 sec)

    % input combinations
    rangeA=[5,6,7,8,9]/2;
    rangeB=(0:50)/2;
    rangeC=(0:10)/2;
    [As,Bs,Cs] =ndgrid(rangeA,rangeB,rangeC);
    all_spect  =zeros(numel(As),2501); 
    all_freq   =zeros(numel(As),2501); 
    all_time_ds=zeros(numel(As),(duration-time_discard)/dt2);
    all_alpha_env_hrf_ds=zeros(numel(As),(duration-time_discard)/dt2);
    all_alpha_env_hrf_bp=zeros(numel(As),(duration-time_discard)/dt2);
    all_gamma_env_hrf_ds=zeros(numel(As),(duration-time_discard)/dt2);
    all_gamma_env_hrf_bp=zeros(numel(As),(duration-time_discard)/dt2);

    all_simBOLD_ds=zeros(numel(As),(duration-time_discard)/dt2);
    all_simBOLD_bp=zeros(numel(As),(duration-time_discard)/dt2);

    r=simulation(0,0,0,duration,time_discard,dt,dt2);
    spect_baseline=r.spect';

    %parfor i=1:numel(As)
    for i=1:numel(As)  

        fprintf('%d/%d...',i,numel(As))
        r=simulation(As(i),Bs(i),Cs(i),duration,time_discard,dt,dt2);

        % -----record-----
        all_spect(i,:)=r.spect'; 
        all_freq(i,:)=r.freq';
        all_time_ds(i,:)=r.time_ds;
        all_alpha_env_hrf_ds(i,:)=r.alpha_env_hrf_ds; 
        all_alpha_env_hrf_bp(i,:)=r.alpha_env_hrf_bp;
        all_gamma_env_hrf_ds(i,:)=r.gamma_env_hrf_ds; 
        all_gamma_env_hrf_bp(i,:)=r.gamma_env_hrf_bp;
        all_simBOLD_ds(i,:)=r.simBOLD_ds; 
        all_simBOLD_bp(i,:)=r.simBOLD_bp;
    
        fprintf('done.\n')
    end
    freq=all_freq(1,:);
    time_ds=all_time_ds(1,:);
    
    filename=sprintf('fig2_batch_result_%s.mat',datestr(now,'yyyymmdd-HHMM'));
    save(fullfile('data',filename),...
         'As','Bs','Cs','rangeA','rangeB','rangeC',...
         'freq','all_spect','spect_baseline',...
         'time_ds',...
         'all_alpha_env_hrf_ds','all_alpha_env_hrf_bp',...
         'all_gamma_env_hrf_ds','all_gamma_env_hrf_bp',...
         'all_simBOLD_ds','all_simBOLD_bp');    

    fprintf('%s saved.\n',filename)

end

%==========================================================================
function r=simulation(a,b,c, duration, time_discard,dt,dt2)

    time_discard=[time_discard duration]; 
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
  
    %-----filtering-----
    alpha = filter_bp(simEEG,[8 12],1/dt); % alpha
    alpha_env = abs(hilbert(alpha));    
    alpha_env_hrf = sim_bold(alpha_env);

    gamma = filter_bp(simEEG,[35 45],1/dt); % gamma
    gamma_env = abs(hilbert(gamma));    
    gamma_env_hrf = sim_bold(gamma_env);

    %-----power spectrum-----
    tidx=find(time>time_discard(1) & time<=time_discard(2)); 
    windowLength = 5000; % timepoints
    data = simEEG(tidx)-mean(simEEG(tidx));
    [spect,freq] = compute_spectrum(data,1/dt,windowLength); 

    %-----downsample----- 
    time_ds=0:dt2:duration;
    simBOLD_ds=interp1(time,simBOLD,time_ds);
    alpha_env_hrf_ds=interp1(time,alpha_env_hrf,time_ds);
    gamma_env_hrf_ds=interp1(time,gamma_env_hrf,time_ds);
    time_ds(isnan(simBOLD_ds))=[];
    simBOLD_ds(isnan(simBOLD_ds))=[];
    alpha_env_hrf_ds(isnan(alpha_env_hrf_ds))=[];
    gamma_env_hrf_ds(isnan(gamma_env_hrf_ds))=[];

    %-----filtering-----
    BPorder=2; 
    simBOLD_bp = filter_bp(simBOLD_ds,[0.008 0.09],1/dt2,BPorder);             
    alpha_env_hrf_bp = filter_bp(alpha_env_hrf_ds,[0.008 0.09],1/dt2,BPorder); 
    gamma_env_hrf_bp = filter_bp(gamma_env_hrf_ds,[0.008 0.09],1/dt2,BPorder);

    % -----record-----
    tidx2=find(time_ds>time_discard(1) & time_ds<=time_discard(2));
    r.spect = spect;
    r.freq  = freq;
    r.time_ds          = time_ds(tidx2);
    r.simBOLD_ds       = simBOLD_ds(tidx2);
    r.alpha_env_hrf_ds = alpha_env_hrf_ds(tidx2);
    r.gamma_env_hrf_ds = gamma_env_hrf_ds(tidx2);
    r.simBOLD_bp       = simBOLD_bp(tidx2);
    r.alpha_env_hrf_bp = alpha_env_hrf_bp(tidx2);
    r.gamma_env_hrf_bp = gamma_env_hrf_bp(tidx2);

end