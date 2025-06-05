% >> for i=1:10,fig4_batch(10);end  % run 100 simulations and save result. (a couple days)
% >> fig4();       % plot fig4
% >> fig4_S1(1);   % plot fig4_S1A
% >> fig4_S1(2);   % plot Fig4_S1B

function fig4_batch(nRuns)

    duration = 600;  % sec per run
    time_discard=[200 550]; 

    % input combinations
    rangeA=[5,6,7,8,9]/2;
    rangeB=(0:50)/2;
    rangeC=(0:10)/2;
    [As,Bs,Cs]=ndgrid(rangeA,rangeB,rangeC);

    all_r_alphaEnvHrfds_simBOLDds = zeros(numel(As),nRuns);
    all_r_alphaEnvHrfbp_simBOLDbp = zeros(numel(As),nRuns);
    all_r_gammaEnvHrfds_simBOLDds = zeros(numel(As),nRuns);
    all_r_gammaEnvHrfbp_simBOLDbp = zeros(numel(As),nRuns);
    all_r_alpha_x = zeros(numel(As),nRuns);
    all_r_gamma_x = zeros(numel(As),nRuns);
    all_r_alphaEnv_x = zeros(numel(As),nRuns);
    all_r_gammaEnv_x = zeros(numel(As),nRuns);
    all_alphaPower= zeros(numel(As),nRuns);
    all_gammaPower= zeros(numel(As),nRuns);


    for run =1:nRuns 
        %parfor i=1:numel(As)
        for i=1:numel(As)  
            fprintf('%d/%d...',i,numel(As))
            r=simulation(As(i),Bs(i),Cs(i),duration,time_discard);

            % -----record-----  
            all_r_alpha_x(i,run) = r.all_r_alpha_x;
            all_r_gamma_x(i,run) = r.all_r_gamma_x;
            all_r_alphaEnv_x(i,run) = r.all_r_alphaEnv_x; 
            all_r_gammaEnv_x(i,run) = r.all_r_gammaEnv_x; 
            all_r_alphaEnvHrfds_simBOLDds(i,run) = r.all_r_alphaEnvHrfds_simBOLDds;
            all_r_alphaEnvHrfbp_simBOLDbp(i,run) = r.all_r_alphaEnvHrfbp_simBOLDbp;
            all_r_gammaEnvHrfds_simBOLDds(i,run) = r.all_r_gammaEnvHrfds_simBOLDds;
            all_r_gammaEnvHrfbp_simBOLDbp(i,run) = r.all_r_gammaEnvHrfbp_simBOLDbp;        
            fprintf('done.\n')
        end
    end % nRuns

    
    filename=sprintf('fig4_batch_result_%s.mat',datestr(now,'yyyymmdd-HHMM'));
    save(fullfile('data',filename),...
        'As','Bs','Cs','rangeA','rangeB','rangeC',...
        'all_r_alphaEnvHrfds_simBOLDds',...
        'all_r_alphaEnvHrfbp_simBOLDbp',...
        'all_r_gammaEnvHrfds_simBOLDds',...
        'all_r_gammaEnvHrfbp_simBOLDbp',...
        'all_r_alpha_x',...
        'all_r_gamma_x',...
        'all_r_alphaEnv_x',... 
        'all_r_gammaEnv_x',...
        'all_alphaPower',...
        'all_gammaPower');    
    fprintf('%s saved.\n',filename) 

end

%==========================================================================
function r=simulation(a,b,c, duration, time_discard)

    dt       = 1e-3; % sec
    dt2      = 1;    % sec

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
    %
    time_ds(isnan(simBOLD_ds))=[];
    simBOLD_ds(isnan(simBOLD_ds))=[];
    alpha_env_hrf_ds(isnan(alpha_env_hrf_ds))=[];
    gamma_env_hrf_ds(isnan(gamma_env_hrf_ds))=[];
    %-----bandpass-----
    BPorder=2;
    simBOLD_bp = filter_bp(simBOLD_ds,[0.008 0.09],1/dt2,BPorder); 
    alpha_env_hrf_bp = filter_bp(alpha_env_hrf_ds,[0.008 0.09],1/dt2,BPorder); 
    gamma_env_hrf_bp = filter_bp(gamma_env_hrf_ds,[0.008 0.09],1/dt2,BPorder); 
    %-----record-----
    tidx2=find(time_ds>time_discard(1) & time_ds<=time_discard(2));
    r.all_r_alpha_x = corr(x(tidx)',alpha(tidx)');
    r.all_r_gamma_x = corr(x(tidx)',gamma(tidx)');
    r.all_r_alphaEnv_x = corr(x(tidx)',alpha_env(tidx)'); 
    r.all_r_gammaEnv_x = corr(x(tidx)',gamma_env(tidx)'); 
    r.all_r_alphaEnvHrfds_simBOLDds = corr(simBOLD_ds(tidx2)',alpha_env_hrf_ds(tidx2)');
    r.all_r_alphaEnvHrfbp_simBOLDbp = corr(simBOLD_bp(tidx2)',alpha_env_hrf_bp(tidx2)');
    r.all_r_gammaEnvHrfds_simBOLDds = corr(simBOLD_ds(tidx2)',gamma_env_hrf_ds(tidx2)');
    r.all_r_gammaEnvHrfbp_simBOLDbp = corr(simBOLD_bp(tidx2)',gamma_env_hrf_bp(tidx2)');

end