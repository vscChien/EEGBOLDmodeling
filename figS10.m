% Show how BOLD signal changes (%) in response to optical stimulation on (1)PV2/3, (2)SOM2/3, and (3)VIP2/3 
% Input condition: I_lat=4.5, I_mod=5, I_tha=0
% Optical stimulations: 1s, 5s, 10s, 20s. (Each for 10 times.)
% >> figS10(0);   % plot null response
% >> figS10(1);   % plot figS10A
% >> figS10(2);   % plot figS10B
% >> figS10(3);   % plot figS10C
%
%
% Show how BOLD signal changes (%) in response to optical stimulation on (1)PV2/3, (2)SOM2/3, and (3)VIP2/3 
% Optical stimulation: 20s. 
% Scan different input conditions: I_lat=[2.5-4.5], I_mod=[0-25], I_tha=[0-5]
% >> figS10_scan(0);   % plot null response
% >> figS10_scan(1);   % plot figS10D
% >> figS10_scan(2);   % plot figS10E
% >> figS10_scan(3);   % plot figS10F


function figS10(cond)

    if nargin<1
        cond=1;
    end

    noiseLevel=1;
    duration = 90;   % sec per run
    dt       = 1e-3;  % simulation time step (sec)  


    % Connectivity W
    [W,C,~,nP,nE,nPV,nSOM,nVIP]=get_W17(0);
    g  = 20*1e3; 
    W = [W*g, ones(nP,1)*1e3];


    % external input
    a=4.5;  % level of lateral 
    b=5;  % level of modulatory
    c=0;  % level of thalamic
    Iext0 = gen_Iext(a,b,c,duration,dt);

    Ts=[1,5,10,20];
    nRuns=10;
    simBOLDs=zeros(duration/dt,length(Ts),nRuns);
    rate=zeros(17,duration/dt,nRuns);

    switch cond
        case 0 % no stim
            label='null';
            EPSV=[0,0,0,0]; % optical stim strength
            figname1='figS10_null1.png'; 
            figname2='figS10_null1.svg'; 
        case 1 % stim on PV
            label='PV2/3';
            EPSV=[0,2,0,0]; % optical stim strength
            figname1='figS10A.png';
            figname2='figS10A.svg';
        case 2 % stim on SOM
            label='SOM2/3';
            EPSV=[0,0,0.3,0]; % optical stim strength
            figname1='figS10B.png';
            figname2='figS10B.svg';
        case 3 % stim on VIP
            label='VIP2/3';
            EPSV=[0,0,0,7]; % optical stim strength
            figname1='figS10C.png';
            figname2='figS10C.svg';
    end

    for i=1:length(Ts)
        for n=1:nRuns
            fprintf('%d, run %d\n',Ts(i),n)
            Tstart = 30; %sec
            Tend   = Tstart+Ts(i); %sec
            Istm = gen_Istm(duration,dt,Tstart,Tend,EPSV);
            Iext = Iext0+Istm;
        
        
                
            % Synaptic time constants
            tau = get_tau(nP,nE,nPV,nSOM,nVIP);
            
            % Sigmoid functions
            sigmParam = get_sigmParam(nE,nPV,nSOM,nVIP);
            
            % simulation
            [rate(:,:,n), time, PSP] = model_jr(W, tau, sigmParam, Iext, noiseLevel, dt);
        
            % simEEG and simBOLD
            OpticMode=1; % OpticaMode: Iext is excluded from the BOLD simulation.
            [simEEG, simBOLDs(:,i,n), x] = sim_eeg_bold_opticmode(PSP,time,nE,nPV,nSOM,nVIP,nP,C,dt,OpticMode);
             %simEEGs(:,i) = sum(simEEG,1);
        end
    end

    %---------plot------------
    %plot1(simBOLDs,simEEGs,rate,time,a,b,c);
    plot3(simBOLDs,rate,time,a,b,c,label);
    saveas(gcf,fullfile('figures',figname1));
    saveas(gcf,fullfile('figures',figname2));
end
%==========================================================================
function plot1(simBOLDs,simEEGs,rate,time,a,b,c)
    tidx=time>25;
    figure;
    subplot(422);plot(time(tidx)-30,rate(1:5,tidx)');   legend({'E2/3','E4','E5ET','E5IT','E6'},'Location','eastoutside');title('E')
    subplot(424);plot(time(tidx)-30,rate(6:9,tidx)');   legend({'P2/3','P4','P5','P6'},'Location','eastoutside');title('PV')
    subplot(426);plot(time(tidx)-30,rate(10:13,tidx)'); legend({'S2/3','S4','S5','S6'},'Location','eastoutside');title('SOM')
    subplot(428);plot(time(tidx)-30,rate(14:17,tidx)'); legend({'V2/3','V4','V5','V6'},'Location','eastoutside');title('VIP')

    subplot(221);plot(time(tidx)-30,simEEGs(tidx,:));
    title(sprintf('baseline[%g, %g, %g] + SOM2/3 \n simEEG',a,b,c));
    legend({'1s','5s','10s','20s'},'Location','eastoutside')

    subplot(223);plot(time(tidx)-30,simBOLDs(tidx,:));title('simBOLD');legend({'1s','5s','10s','20s'},'Location','eastoutside')

end
%==========================================================================
% simBOLDs [90000 timepoints x 4 conditions]
% rate [17 populations x 90000 timepoints] % last condition
% time [1 x 90000]
function plot2(simBOLDs,rate,time,a,b,c,label)

    meanRate= [mean(rate(:,time>25 & time < 30),2),...% baseline
               mean(rate(:,time>40 & time < 45),2)];  % during stimulus

    tidx=time>25;
    figure;%set(gcf,'position',[0 0 800 200])
    t = tiledlayout(1,6,'TileSpacing','compact','Padding','compact');
    width=20; 
    height=6;
    set(gcf,'units','centimeters','position',[2 2 width height])

    nexttile([1 2])
    meanBOLDs=mean(simBOLDs(time>25 & time<30,:),1);
    tmp=(simBOLDs-repmat(meanBOLDs,[length(time),1]))./repmat(meanBOLDs,[length(time),1])*100;
    plot(time(tidx)-30,tmp(tidx,:),'LineWidth',1.5);
    title(sprintf('I_{ext} = [%g, %g, %g]\n stim. on %s',a,b,c,label));
    lgd=legend({'1s','5s','10s','20s'},'Location','southeast');lgd.ItemTokenSize = [8,1];
    xlabel('Time (sec)');ylabel('BOLD (%)');box off;xlim([-5 60])

    nexttile();
    bar(meanRate(1:5,:));
    set(gca,'XTickLabel',{'E2/3','E4','E5ET','E5IT','E6'});
    lgd=legend({'stim Off','stim On'},'Location','northoutside');lgd.ItemTokenSize = [8,1];
    ylabel('Mean firing rate (Hz)'); box off;
    nexttile();
    bar(meanRate(6:9,:));
    set(gca,'XTickLabel',{'PV2/3','PV4','PV5','PV6'});box off;
    nexttile();
    bar(meanRate(10:13,:));
    set(gca,'XTickLabel',{'SOM2/3','SOM4','SOM5','SOM6'});box off;
    nexttile();
    bar(meanRate(14:17,:));
    set(gca,'XTickLabel',{'VIP2/3','VIP4','VIP5','VIP6'});box off;


   
    

end
%==========================================================================
function plot3(simBOLDs, rate, time, a, b, c, label)

    % ================================
    % Inputs:
    % simBOLDs: [T x 4 x nRuns]
    % rate:     [17 x T x nRuns]
    % time:     [1 x T]
    % ================================
    T = length(time);
    nRuns = size(simBOLDs,3);

    % ---- Time windows ----
    stim_onset    = 30; % stim: [30~50]s
    stim_offset   = 50; % 
    baseline_idx = time > (stim_onset-10) & time < stim_onset;  % 20~30  
    stimOn_idx     = time > (stim_offset-10) & time < (stim_offset-5); % 40~45
    stimOff_idx    = time > (time(end)-10) & time < (time(end)-5); % 80~85
    plot_idx     = time > (stim_onset-10); % >20
    % ==============================================================
    % ==== Compute mean/stdev firing rate for baseline vs stimulus
    % ==============================================================

    Rate_stimOn  = reshape(rate(:,stimOn_idx,:),17,[]); % [17 x t*nRuns]
    Rate_stimOff = reshape(rate(:,stimOff_idx,:),17,[]);     % [17 x t*nRuns]

    meanRate = [ mean(Rate_stimOff,2) , mean(Rate_stimOn,2) ];   % [17 x 2]
    stdRate  = [ std(Rate_stimOff,[],2) , std(Rate_stimOn,[],2) ]; % [17 x 2]

    % ==============================================================
    % ==== Prepare figure and layout
    % ==============================================================
    figure;
    t = tiledlayout(1,6,'TileSpacing','compact','Padding','compact');
    set(gcf,'Units','centimeters','Position',[2 2 20 7]);

    % ==============================================================
    % ==== BOLD time series: mean + shaded STD
    % ==============================================================

    nexttile([1 2]);

    % Convert BOLD to percent signal change for each run
    meanBOLDs = squeeze(mean(simBOLDs(baseline_idx,:,:),1)); % [4 x nRuns]
    tmp = zeros(sum(plot_idx),4,nRuns);

    for r = 1:nRuns
        baseline = meanBOLDs(:,r)';
        tmp(:,:,r) = ((simBOLDs(plot_idx,:,r) - repmat(baseline,[sum(plot_idx),1])) ./ repmat(baseline,[sum(plot_idx),1])) * 100;
    end

    % Mean & STD across runs
    meanBOLD = mean(tmp,3);       % [Tplot x 4]
    stdBOLD  = std(tmp,0,3);      % [Tplot x 4]

    hold on;
    colors = lines(4);
    h=[];
    for cond = 1:4
        tvec = time(plot_idx)' - stim_onset;

        % Shaded area (STD)
        upperbound=meanBOLD(:,cond) + stdBOLD(:,cond);
        lowerbound=meanBOLD(:,cond) - stdBOLD(:,cond);
        fill([tvec(1:500:end); flipud(tvec(1:500:end))], ...
             [upperbound(1:500:end); flipud(lowerbound(1:500:end))], ...
             colors(cond,:), 'FaceAlpha',0.2, 'EdgeColor','none');

        % Mean line
        h(cond)=plot(tvec, meanBOLD(:,cond),'Color',colors(cond,:),'LineWidth',1.5);        
    end
    title(sprintf('I_{ext} = [%g, %g, %g]\n Stim. on %s',a,b,c,label));
    set(gca,'fontsize',10,'FontName', 'calibri')
    lgd=legend(h,{'1s','5s','10s','20s'},'fontsize',8,'FontName', 'calibri','Location','southeast');lgd.ItemTokenSize = [8,1];
    
    xlabel('Time (sec)','fontsize',10,'FontName', 'calibri');
    ylabel('BOLD (%)','fontsize',10,'FontName', 'calibri');

    box off; xlim([tvec(1) 60]);
    %ax = gca;
    %ax.FontSize = 10;   % axis tick label ��2
    % ==============================================================
    % ===== Bar plots for rate changes: mean + STD whiskers
    % ==============================================================
    groups = {1:5, 6:9, 10:13, 14:17};
    labels = {
        {'E2/3','E4','E5ET','E5IT','E6'}
        {'PV2/3','PV4','PV5','PV6'}
        {'SOM2/3','SOM4','SOM5','SOM6'}
        {'VIP2/3','VIP4','VIP5','VIP6'}
    };

    for k = 1:4
        nexttile();
        idx = groups{k};

        bb = bar(meanRate(idx,:));hold on;
        hold on;

        % Add error bars
        [ngroups, nbars] = size(meanRate(idx,:));
        x = nan(nbars, ngroups);
        for i = 1:nbars
            x(i,:) = bb(i).XEndPoints;
        end

        errorbar(x', meanRate(idx,:), stdRate(idx,:), 'k','LineStyle','none','LineWidth',1,'CapSize',0);

        set(gca,'XTickLabel', labels{k},'fontsize',8,'FontName', 'calibri');
        box off;

        if k == 1
            lgd=legend({'stim Off','stim On'},'fontsize',8,'FontName', 'calibri','Location','northoutside');
            lgd.ItemTokenSize = [8,1];
            ylabel('Mean firing rate (Hz)','fontsize',10,'FontName', 'calibri');
        end
    end
    
end

%==========================================================================
% a: level of lateral 
% b: level of modulatory
% c: level of thalamic
function Iext = gen_Iext(a,b,c,duration,dt)

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
   
    Iext  = (Ilat*a + Imod*b + Itha*c)*ones(1,duration/dt);  
end

%==========================================================================
% a: level of lateral 
% b: level of modulatory
% c: level of thalamic
function Istm = gen_Istm(duration,dt,Tstart,Tend,EPSV)

    % Input weights (to 17 populations)
    w    = [EPSV(1); %to E2/3
            0; %to E4
            0; %to E5ET
            0; %to E5IT
            0; %to E6
            EPSV(2); %to PV2/3 
            0; %to PV4
            0; %to PV5
            0; %to PV6
            EPSV(3); %to SOM2/3
            0; %to SOM4
            0; %to SOM5
            0; %to SOM6
            EPSV(4); %to VIP2/3 
            0; %to VIP4    
            0; %to VIP5 
            0];%to VIP6 
      
    input = zeros(1,duration/dt);
    time  = dt:dt:duration;
    input(time>Tstart & time<=Tend)=1;
    Istm  = w*input;  
end