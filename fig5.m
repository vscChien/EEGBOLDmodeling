
function fig5()

    % input conditions
    cond(1).Ilat_range=[4.5,4.5];
    cond(1).Imod_range=[0,2];
    cond(1).Itha_range=[0,2];

    cond(2).Ilat_range=[4.5,4.5];
    cond(2).Imod_range=[12.5,14.5];
    cond(2).Itha_range=[0,2];    

    cond(3).Ilat_range=[4.5,4.5];
    cond(3).Imod_range=[16,18];
    cond(3).Itha_range=[0,2];  


    %------ plot input conditions ------
    plot_fig5a(cond);

    %------ plot single trial ------
    % load random-walk input
    [input,duration,dt]=load_input();  

    % check saved result
    noiseLevel=1;
    resultfilename=fullfile('data',sprintf('fig5_result_10min_noise=%g.mat',noiseLevel));    
    if exist(resultfilename,'file')
        fprintf('Loading %s...',resultfilename);
        load(resultfilename);
        fprintf('done.\n');
    else
        nConds= length(cond); 
        nRuns = 100; 
        %----------record---------------
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
        
        save(resultfilename,'all_r_alphaEnv_x','all_r_gammaEnv_x',...
                               'all_r_alphaEnvHrfds_simBOLDds','all_r_alphaEnvHrfbp_simBOLDbp',...
                               'all_r_gammaEnvHrfds_simBOLDds','all_r_gammaEnvHrfbp_simBOLDbp',...
                               'all_r_rate','cond','noiseLevel');
    end
    
    plot_histogram(cond, all_r_alphaEnvHrfbp_simBOLDbp,all_r_gammaEnvHrfbp_simBOLDbp,'fig5C.png');

    plotOn=1; rmbase=1;
    for i=1:3
        simulation(input(:,:,1),duration,dt,noiseLevel,cond(i).Ilat_range,cond(i).Imod_range,cond(i).Itha_range,plotOn,rmbase,sprintf('fig5B_%d.png',i));
    end

end
%==========================================================================
function plot_histogram(cond, all_r_alphaEnvHrfbp_simBOLDbp,all_r_gammaEnvHrfbp_simBOLDbp,filename)
        
        figure;
        t=tiledlayout(1,3,'TileSpacing','compact','Padding','loose');
        for i=1:length(cond)
            nexttile();
            histogram(all_r_alphaEnvHrfbp_simBOLDbp(i,:),-1:0.05:1,'Normalization','probability');hold on; 
            histogram(all_r_gammaEnvHrfbp_simBOLDbp(i,:),-1:0.05:1,'Normalization','probability');           
            xlim([-1 1]);ylim([0 0.5]);
            plot([0;0],ylim,'k');
             

            mean_r_alpha=mean(all_r_alphaEnvHrfbp_simBOLDbp(i,:));
            mean_r_gamma=mean(all_r_gammaEnvHrfbp_simBOLDbp(i,:));
            plot([1;1]*mean_r_alpha,ylim,'b');
            plot([1;1]*mean_r_gamma,ylim,'r');
            title(sprintf('Case %d:',i))
            text(round(mean_r_alpha,2)-0.22,0.3,sprintf('$r_{\\alpha}$=%g',round(mean_r_alpha,2)),...
                'color','b','fontsize',8,'FontName','calibri','Interpreter','latex');
            text(min(0.25,round(mean_r_gamma,2)-0.22),0.4,sprintf('$r_{\\gamma}$=%g',round(mean_r_gamma,2)),...
                'color','r','fontsize',8,'FontName','calibri','Interpreter','latex');



            if i~=1,set(gca,'yticklabel',[]);end
            title(sprintf('Case %d',i));
            set(gca,'fontsize',8,'FontName','calibri')
        end      
        xlabel(t,'Correlation','fontsize',10,'FontName', 'calibri'); 
        ylabel(t,'Probability','fontsize',10,'FontName', 'calibri');
        legend({'Alpha-BOLD','Gamma-BOLD'},'NumColumns', 1,'Location','southoutside'...
                ,'fontsize',8,'FontName', 'calibri')

        width=15; 
        height=6;
        set(gcf,'units','centimeters','position',[2 2 width height])
        saveas(gcf,fullfile('figures',filename));
end

%==========================================================================
function plot_single_trial(simMEG,simBOLD_bp,alpha_env_hrf_bp,gamma_env_hrf_bp,...
                           time_discard,time_ds,time,dt,Iext,rmbase,filename)

    figure;
    t=tiledlayout(3,1,'TileSpacing','compact','Padding','loose');
    nexttile();
    plot(time,Iext(1,:),'linewidth',1);hold on; grid on;
    plot(time,Iext(14,:),'linewidth',1); 
    plot(time,Iext(2,:),'linewidth',1); 
    xlim(time_discard); ylim([0 20])
    legend('I_{lat}','I_{mod}','I_{tha}','location','northeastoutside');
    ylabel('Rate (Hz)');title('Input')
    set(gca,'fontsize',8,'FontName','calibri','xticklabel',[])
    
    
    nexttile();
    [S,F,T] = spectrogram(simMEG-mean(simMEG),2000,1000,2000,1/dt);
    tmp=10*log(abs(S));
    if rmbase % remove baseline
        load(fullfile('data','fig5_baseline.mat'),'baseline');
        tmp=tmp-repmat(baseline,[1, size(tmp,2)]);
    else % remove mean
        tmp=tmp-repmat(mean(tmp,2),[1, size(tmp,2)]);   
    end
    mvstep=10; 
    tmp_m = movmean(tmp,mvstep,2);
    imagesc(T,F,tmp_m);set(gca,'ydir','normal');
    ylim([0.5 50]);xlim(time_discard);
    clim([-1 1]*25);
    colormap(jet);
    pos=get(gca,'position');
    colorbar('position',[0.92 pos(2)+0.02 0.01 0.15]); 
    hold on;plot(xlim,[1;1]*[10,40],'k:'); 
    ylabel('Freq (Hz)');title('Spectrogram')
    set(gca,'fontsize',8,'FontName','calibri','xticklabel',[],...
            'ytick',0:10:50,'yticklabel',0:10:50)

    nexttile();
    plot(time_ds,zscore(alpha_env_hrf_bp),'linewidth',1); hold on;grid on;
    plot(time_ds,zscore(gamma_env_hrf_bp),'linewidth',1); 
    plot(time_ds,zscore(simBOLD_bp),'linewidth',1); 
    xlim(time_discard);
    legend('Alpha','Gamma','BOLD','location','northeastoutside');

    ylabel('z-score');
    title(sprintf('$r_{\\alpha} = %g;\\; r_{\\gamma} = %g$',...
                            round(corr(simBOLD_bp',alpha_env_hrf_bp'),2),...
                            round(corr(simBOLD_bp',gamma_env_hrf_bp'),2)),'Interpreter','latex');
    xlabel('Time (sec)');ylim([-3 3])
    set(gca,'fontsize',8,'FontName','calibri')

    width=10; 
    height=9;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures',filename));   
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

    if plotOn
        plot_single_trial(simEEG,simBOLD_bp,alpha_env_hrf_bp,gamma_env_hrf_bp,...
                          time_discard,time_ds,time,dt,Iext,rmbase,filename);
    end

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
        fprintf('Generating rand inputs...\n');
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

%==========================================================================
function plot_fig5a(cond)

    figure;
    t=tiledlayout(1,1,'TileSpacing','compact','Padding','loose');    
    load(fullfile('data','fig2_kmeans8.mat'),'k','tag','rangeA','rangeB','rangeC','As'); 
    tag=reshape(tag,size(As));
    nexttile();
    imagesc(rangeB,rangeC,squeeze(tag(ismember(rangeA,4.5),:,:))');clim([0 k]+0.5);
    set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri', ...
        'xtick',0:5:25,'xticklabel',0:5:25)
    xlabel('I_{mod}');ylabel('I_{tha}')
    title(sprintf('I_{lat} = %g',rangeA(ismember(rangeA,4.5))));
    colormap(brewermap(k,'Paired'))
    hold on;
    for i=1:length(cond)
            rectangle('Position',[cond(i).Imod_range(1),...
                                  cond(i).Itha_range(1),...
                                  diff(cond(i).Imod_range),...
                                  diff(cond(i).Itha_range)],'LineWidth',1,'EdgeColor','k')
    end

    colorbar('ytick',1:k,'yticklabel',1:k);
    width=7; 
    height=5.5;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig5A.png'));
end