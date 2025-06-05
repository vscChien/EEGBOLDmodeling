% >> for i=1:10,fig4_batch(10);end  % run 100 simulations and save result. (a couple days)
% >> fig4();       % plot fig4
% >> fig4_S1(1);   % plot fig4_S1A
% >> fig4_S1(2);   % plot Fig4_S1B

function fig4_S1(ex)

    if nargin<1
        ex=1;
    end
    
    switch ex
        case 1 % example 1: alpha-BOLD correlation < 0 
            Ilat=4.5;
            Imod=17;
            Itha=0.5;
        case 2 % example 2: alpha-BOLD correlation > 0 
            Ilat=4.5;
            Imod=14;
            Itha=4;
    end
    
    resultfilename=fullfile('data',sprintf('fig4_S1_result_[%g,%g,%g].mat',Ilat,Imod,Itha));
   
    
    if exist(resultfilename,'file')
        fprintf('Loading %s...',resultfilename);
        load(resultfilename);
        fprintf('done.\n');
    else
        nRuns_simulation=100;
        nRuns_randdipole=100;
        %-------initialization------
        all_r_simEEG_x=zeros(nRuns_simulation,nRuns_randdipole);
        all_r_alpha_x=zeros(nRuns_simulation,nRuns_randdipole);
        all_r_gamma_x=zeros(nRuns_simulation,nRuns_randdipole);
        all_r_alphaEnv_x=zeros(nRuns_simulation,nRuns_randdipole);
        all_r_gammaEnv_x=zeros(nRuns_simulation,nRuns_randdipole);
        all_r_alphaEnvHrfds_simBOLDds=zeros(nRuns_simulation,nRuns_randdipole);
        all_r_alphaEnvHrfbp_simBOLDbp=zeros(nRuns_simulation,nRuns_randdipole);
        all_r_gammaEnvHrfds_simBOLDds=zeros(nRuns_simulation,nRuns_randdipole);
        all_r_gammaEnvHrfbp_simBOLDbp=zeros(nRuns_simulation,nRuns_randdipole);
        
        for i=1:nRuns_simulation

            fprintf('-----Simulation %d/%d-----\n',i,nRuns_simulation)
            [PSP,time,info]=simulation(Ilat,Imod,Itha);
            stat=random_dipole_analysis(PSP,time,info,nRuns_randdipole);
            all_r_simEEG_x(i,:)=stat.all_r_simEEG_x;
            all_r_alpha_x(i,:)=stat.all_r_alpha_x;
            all_r_gamma_x(i,:)=stat.all_r_gamma_x;
            all_r_alphaEnv_x(i,:)=stat.all_r_alphaEnv_x;
            all_r_gammaEnv_x(i,:)=stat.all_r_gammaEnv_x;
            all_r_alphaEnvHrfds_simBOLDds(i,:)=stat.all_r_alphaEnvHrfds_simBOLDds;
            all_r_alphaEnvHrfbp_simBOLDbp(i,:)=stat.all_r_alphaEnvHrfbp_simBOLDbp;
            all_r_gammaEnvHrfds_simBOLDds(i,:)=stat.all_r_gammaEnvHrfds_simBOLDds;
            all_r_gammaEnvHrfbp_simBOLDbp(i,:)=stat.all_r_gammaEnvHrfbp_simBOLDbp;
   
        end
    
        fprintf('Saving %s...',resultfilename);
        save(resultfilename,'Ilat','Imod','Itha','info','nRuns_simulation','nRuns_randdipole',...
            'all_r_simEEG_x','all_r_alpha_x','all_r_gamma_x','all_r_alphaEnv_x',...
            'all_r_gammaEnv_x','all_r_alphaEnvHrfds_simBOLDds','all_r_alphaEnvHrfbp_simBOLDbp',...
            'all_r_gammaEnvHrfds_simBOLDds','all_r_gammaEnvHrfbp_simBOLDbp');   
    
        fprintf('done.\n');
    end
      
    figure;
    t=tiledlayout(5,2,'TileSpacing','compact','Padding','loose');
    nexttile();
    plot_histogram(all_r_simEEG_x,'r(EEG,z)',[0 1]);
    nexttile(3);
    plot_histogram(all_r_alpha_x,'r(alpha,z)',[0 1]);
    nexttile(4);
    plot_histogram(all_r_gamma_x,'r(gamma,z)',[0 1]);
    nexttile(5);
    plot_histogram(all_r_alphaEnv_x,'r(alpha-env,z)',[0 .5]);
    nexttile(6);
    plot_histogram(all_r_gammaEnv_x,'r(gamma-env,z)',[0 .5]);
    nexttile(7);
    plot_histogram(all_r_alphaEnvHrfds_simBOLDds,'r(alpha-env-hrf,BOLD)',[0 .3]);
    nexttile(8);
    plot_histogram(all_r_gammaEnvHrfds_simBOLDds,'r(gamma-env-hrf,BOLD)',[0 .3]);
    nexttile(9);
    plot_histogram(all_r_alphaEnvHrfbp_simBOLDbp,'r(alpha-env-hrf-bp,BOLD-bp)',[0 .2]);
    nexttile(10);
    plot_histogram(all_r_gammaEnvHrfbp_simBOLDbp,'r(gamma-env-hrf-bp,BOLD-bp)',[0 .2]);
    width=10; 
    height=15;
    set(gcf,'units','centimeters','position',[2 2 width height])
    xlabel(t,'Correlation');
    ylabel(t,'Probability');
    figfilename=fullfile('figures',sprintf('fig4_S1[%g,%g,%g].png',Ilat,Imod,Itha));
    saveas(gcf,figfilename);
end

%==========================================================================
function plot_histogram(data,ttext,ylimits)

    nbins=20;
    h=histogram(data(:),nbins,'Normalization','probability');hold on; 
    h=histogram(data(:,1),h.BinEdges,'Normalization','probability'); 
    xlim([-1 1]*min(1,max(abs(data(:)))*1.1));
    if ~isempty(ylimits),ylim(ylimits);end
    title(ttext);
    r1=mean(data(:,1));r2=mean(data(:));
    plot([1;1]*r1,ylim,'r','linewidth',1);
    plot([1;1]*r2,ylim,'b','linewidth',1);
    xlimits=xlim;
    ylimits=ylim;
    text(xlimits(1)*0.95,ylimits(2)*0.9,sprintf('$r_1$= %g',round(r1,2)),'color','r','fontsize',8,'Interpreter','latex');
    text(xlimits(1)*0.95,ylimits(2)*0.7,sprintf('$r_2$= %g',round(r2,2)),'color','b','fontsize',8,'Interpreter','latex');
    plot([0;0],ylim,'k');ylim(ylimits);
    set(gca,'fontsize',8,'FontName','calibri');
end
%==========================================================================
function [PSP,time,info]=simulation(a,b,c)

    duration = 600;   % 600 sec per run
    dt       = 1e-3;  % simulation time step (sec)
    dt2      = 1;     % downsampled time step (sec)
    time_discard=[200 550]; 
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
    [~, time, PSP] = model_jr(W, tau, sigmParam, Iext, noiseLevel, dt);


    info.nE=nE;
    info.nPV=nPV;
    info.nSOM=nSOM;
    info.nVIP=nVIP;
    info.nP=nP;
    info.C=C;
    info.dt=dt;
    info.time_discard=time_discard;
    info.dt2=dt2;
    info.duration=duration;
end
%==========================================================================
function stat=random_dipole_analysis(PSP,time,info,nRun)

    nE=info.nE;
    nPV=info.nPV;
    nSOM=info.nSOM;
    nVIP=info.nVIP;
    nP=info.nP;
    C=info.C;
    duration=info.duration;
    dt=info.dt;
    dt2=info.dt2;
    time_discard=info.time_discard;


    fprintf('simEEG: %d/%d\n',1,nRun);
    simEEGs=zeros(nRun,length(time));
    [simEEG, simBOLD, x] = sim_eeg_bold(PSP,time,nE,nPV,nSOM,nVIP,nP,C,dt);
    simEEGs(1,:)=sum(simEEG,1);
    for i=2:nRun
        fprintf('simEEG: %d/%d\n',i,nRun);
        [simEEG, simBOLD, x] = sim_eeg_bold_testRandDipole(PSP,time,nE,nPV,nSOM,nVIP,nP,C,dt);
        simEEGs(i,:)=sum(simEEG,1);
    end

    stat.all_r_simMEG_x = zeros(1,nRun);
    stat.all_r_alpha_x = zeros(1,nRun);
    stat.all_r_gamma_x = zeros(1,nRun);
    stat.all_r_alphaEnv_x = zeros(1,nRun);
    stat.all_r_gammaEnv_x = zeros(1,nRun);
    stat.all_r_alphaEnvHrfds_simBOLDds = zeros(1,nRun);
    stat.all_r_alphaEnvHrfbp_simBOLDbp = zeros(1,nRun);
    stat.all_r_gammaEnvHrfds_simBOLDds = zeros(1,nRun);
    stat.all_r_gammaEnvHrfbp_simBOLDbp = zeros(1,nRun);

    for i=1:nRun
       fprintf('Correlation: %d/%d\n',i,nRun);
       alpha = filter_bp(simEEGs(i,:),[8 12],1/dt); % alpha
       alpha_env = abs(hilbert(alpha));    % amplitude of alpha
       alpha_env_hrf = sim_bold(alpha_env);% use alpha to generate BOLD
       gamma = filter_bp(simEEGs(i,:),[35 45],1/dt); % gamma
       gamma_env = abs(hilbert(gamma));    % amplitude of gamma
       gamma_env_hrf = sim_bold(gamma_env);% use alpha to generate BOLD

       % -----downsample-----
       tidx=find(time>time_discard(1) & time<=time_discard(2));   
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
       stat.all_r_simEEG_x(i) = corr(x(tidx)',simEEGs(i,tidx)');
       stat.all_r_alpha_x(i) = corr(x(tidx)',alpha(tidx)');
       stat.all_r_gamma_x(i) = corr(x(tidx)',gamma(tidx)');
       stat.all_r_alphaEnv_x(i) = corr(x(tidx)',alpha_env(tidx)');
       stat.all_r_gammaEnv_x(i) = corr(x(tidx)',gamma_env(tidx)');
       stat.all_r_alphaEnvHrfds_simBOLDds(i) = corr(simBOLD_ds(tidx2)',alpha_env_hrf_ds(tidx2)');
       stat.all_r_alphaEnvHrfbp_simBOLDbp(i) = corr(simBOLD_bp(tidx2)',alpha_env_hrf_bp(tidx2)');
       stat.all_r_gammaEnvHrfds_simBOLDds(i) = corr(simBOLD_ds(tidx2)',gamma_env_hrf_ds(tidx2)');
       stat.all_r_gammaEnvHrfbp_simBOLDbp(i) = corr(simBOLD_bp(tidx2)',gamma_env_hrf_bp(tidx2)');
    end
    
end