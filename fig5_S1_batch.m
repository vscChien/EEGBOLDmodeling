
function fig5_S1_batch(ex)

    switch ex
        case 1
            Ilat=4.5;
            Imod=1; 
            Itha=1; 
        case 2
            Ilat=4.5;
            Imod=17; 
            Itha=1;  
    end

    steps=0:0.05:1;
    % input conditions
    for i=1:length(steps)
        cond(i).Ilat_range=[Ilat, Ilat];
        cond(i).Imod_range=[Imod-steps(i),Imod+steps(i)]; 
        cond(i).Itha_range=[Itha-steps(i),Itha+steps(i)]; 
    end  

    %------ plot single trial ------
    % load random-walk input
    [input,duration,dt]=load_input();  

    % check saved result
    noiseLevel=1;
    resultfilename=fullfile('data',...
        sprintf('fig5_S1_[%g,%g,%g]_10min_noise=%g.mat',Ilat,Imod,Itha,noiseLevel));  
    if exist(resultfilename,'file')
        fprintf('Loading %s...',resultfilename);
        load(resultfilename);
        fprintf('done.\n');
    else
        nConds= length(cond); 
        nRuns = 100; 
        %----------record%---------------
        all_r_alphaEnvHrfds_simBOLDds = zeros(nConds,nRuns);
        all_r_alphaEnvHrfbp_simBOLDbp = zeros(nConds,nRuns);
        all_r_gammaEnvHrfds_simBOLDds = zeros(nConds,nRuns);
        all_r_gammaEnvHrfbp_simBOLDbp = zeros(nConds,nRuns);
        all_r_alphaEnv_x = zeros(nConds,nRuns);
        all_r_gammaEnv_x = zeros(nConds,nRuns);
        all_r_rate=zeros(nConds,19,19,nRuns); % correlation of rates
        %---------------------------------
        for c=1:nConds
            for run =1:nRuns 
                fprintf('run%d/%d...cond%d/%d...\n',run,nRuns,c,nConds)
                r = simulation(input(:,:,run),duration,dt,noiseLevel,...
                               cond(c).Ilat_range,cond(c).Imod_range,cond(c).Itha_range,0);
                all_r_alphaEnv_x(c,run) = r.r_alphaEnv_x; 
                all_r_gammaEnv_x(c,run) = r.r_gammaEnv_x; 
                all_r_alphaEnvHrfds_simBOLDds(c,run) = r.r_alphaEnvHrfds_simBOLDds;
                all_r_alphaEnvHrfbp_simBOLDbp(c,run) = r.r_alphaEnvHrfbp_simBOLDbp;
                all_r_gammaEnvHrfds_simBOLDds(c,run) = r.r_gammaEnvHrfds_simBOLDds;
                all_r_gammaEnvHrfbp_simBOLDbp(c,run) = r.r_gammaEnvHrfbp_simBOLDbp;
                all_r_rate(c,:,:,run)=r.r_rate;
            
            end % nRuns
        end % nConds
        
        save(resultfilename,...
            'all_r_alphaEnv_x','all_r_gammaEnv_x',...
            'all_r_alphaEnvHrfds_simBOLDds','all_r_alphaEnvHrfbp_simBOLDbp',...
            'all_r_gammaEnvHrfds_simBOLDds','all_r_gammaEnvHrfbp_simBOLDbp',...
            'all_r_rate','cond','noiseLevel');
    end

end

%==========================================================================
function r=simulation(input,duration,dt,noiseLevel,Ilat_range,Imod_range,Itha_range,plotOn,rmbase,filename)

    dt2      = 1;     % downsampled time step (sec)
    time_discard=[200 duration-50]; 
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


    Iext  = Ilat*rescale(input(:,1),Ilat_range(1),Ilat_range(2))' + ...
            Imod*rescale(input(:,2),Imod_range(1),Imod_range(2))' + ...
            Itha*rescale(input(:,3),Itha_range(1),Itha_range(2))'; 

    % Connectivity W
    [W,C,~,nP,nE,nPV,nSOM,nVIP]=get_W17(0);
    g  = 20*1e3; 
    W = [W*g, ones(nP,1)*1e3];
    
    
    % Synaptic time constants
    tau = get_tau(nP,nE,nPV,nSOM,nVIP);
    
    % Sigmoid functions
    sigmParam = get_sigmParam(nE,nPV,nSOM,nVIP);
    
    [rate, time, PSP] = model_jr(W, tau, sigmParam, Iext, noiseLevel, dt);
    [simEEG, simBOLD, x] = sim_eeg_bold(PSP,time,nE,nPV,nSOM,nVIP,nP,C,dt);
     simEEG=sum(simEEG,1);


    %------ filtering ------
    alpha = filter_bp(simEEG,[8 12],1/dt); % alpha
    alpha_env = abs(hilbert(alpha));       % amplitude of alpha
    alpha_env_hrf = sim_bold(alpha_env);   % generate BOLD

    gamma = filter_bp(simEEG,[35 45],1/dt); % gamma
    gamma_env = abs(hilbert(gamma));        % amplitude of gamma
    gamma_env_hrf = sim_bold(gamma_env);    % generate BOLD

    %------ correlation (rate) ------
    %-----power spectrum-----
    tidx=find(time>time_discard(1) & time<=time_discard(2)); 
    tmp=[rate;simEEG;x];
    r_rate=corr(tmp(:,tidx)');

    %-----downsample----- 
    time_ds=0:dt2:duration;
    simBOLD_ds=interp1(time(tidx),simBOLD(tidx),time_ds);
    alpha_env_hrf_ds=interp1(time(tidx),alpha_env_hrf(tidx),time_ds);
    gamma_env_hrf_ds=interp1(time(tidx),gamma_env_hrf(tidx),time_ds);
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
    % -----record-----
    tidx2=find(time_ds>time_discard(1) & time_ds<=time_discard(2));
    r.r_alphaEnv_x = corr(x(tidx)',alpha_env(tidx)'); 
    r.r_gammaEnv_x = corr(x(tidx)',gamma_env(tidx)'); 
    r.r_alphaEnvHrfds_simBOLDds = corr(simBOLD_ds(tidx2)',alpha_env_hrf_ds(tidx2)');
    r.r_alphaEnvHrfbp_simBOLDbp = corr(simBOLD_bp(tidx2)',alpha_env_hrf_bp(tidx2)');
    r.r_gammaEnvHrfds_simBOLDds = corr(simBOLD_ds(tidx2)',gamma_env_hrf_ds(tidx2)');
    r.r_gammaEnvHrfbp_simBOLDbp = corr(simBOLD_bp(tidx2)',gamma_env_hrf_bp(tidx2)');
    r.r_rate=r_rate;

end

%==========================================================================
% load random-walk input
function [input,duration,dt]=load_input()

    speed    = 5; 
    mvwindow=0*1e3; 
    inputfilename=fullfile('data',sprintf('randWalk_input_speed%d_mv%d.mat',speed,mvwindow/1e3)); 
    if exist(inputfilename,'file')
        fprintf('Loading %s...',inputfilename);
        load(inputfilename);
        fprintf('done.\n');
    else
        duration = 600; % sec
        dt       = 0.001;
        nTrials=100;        
        input=zeros(duration/dt,3,nTrials); % [time x {I_lat,I_mod,I_tha} x trials]
        for i=1:3
            for j=1:nTrials
                input(:,i,j)=gen_randwalk_input(duration,dt,speed)';                                     
            end
        end
        if mvwindow>0 % smooth input
            for i=1:3
                for j=1:nTrials
                    input(:,i,j)=movmean(input(:,i,j),mvwindow);
                end
            end
        end
        fprintf('Saving %s...',inputfilename);
        save(inputfilename,'input','duration','dt','speed','nTrials','mvwindow');
        fprintf('done.\n');
    end
end
