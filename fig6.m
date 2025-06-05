% >> fig6(1); % plot Fig.6
% >> fig6(2); % plot Fig.S4 example 1
% >> fig6(3); % plot Fig.S5 example 2
% >> fig6(4); % plot Fig.S6 example 3

function fig6(fig)

    if nargin<1
        fig=1;
    end

    switch fig   
        case 1 
            Ilat_range=[2.5,3];
            Imod_range=[0,25];
            Itha_range=[0,5];
            mvwindow  = 0*1e3; 
        case 2 
            Ilat_range=[4.5,4.5];
            Imod_range=[12.5,14.5];
            Itha_range=[0,2];
            mvwindow  = 0*1e3; 
        case 3 
            Ilat_range=[2.5,3];
            Imod_range=[0,25];
            Itha_range=[0,5];
            mvwindow  = 5*1e3; 
        case 4 
            Ilat_range=[3.5,4.5];
            Imod_range=[0,25];
            Itha_range=[0,5];
            mvwindow  = 0*1e3; 
    end
    % load random-walk input
    [input,duration,dt,speed,nTrials,inputfilename]=load_input(mvwindow); 
    inputTag = input2tag(input,Ilat_range,Imod_range,Itha_range);

    % check saved result
    resultfilename=fullfile('data',...
        sprintf('fig6_result_randspeed%d_mv%d_lat%g%g_mod%g%g_tha%g%g.mat',speed,mvwindow/1e3,Ilat_range,Imod_range,Itha_range));
    if exist(resultfilename,'file')
        fprintf('Loading %s...',resultfilename);
        load(resultfilename);
        fprintf('done.\n');
    else % run simulations and save result
        nRuns = 100; 
        %----------record---------------
        tiles=0:10:100; 
        r_overall = zeros(200,nRuns); % 200 fbins of spectrogram (0.5:0.5:100 Hz)
        r_pertile = zeros(200,length(tiles)-1,nRuns); 
        stat(nRuns).alphaz=[]; %zscored alpha;
        stat(nRuns).gammaz=[]; %zscored gamma;
        stat(nRuns).LFPpowerz=[];
        stat(nRuns).simBOLDbpz=[];
        stat(nRuns).PLabel=[];
        stat(nRuns).inputTag=[];
        %---------------------------------
        for run =1:nRuns 
            fprintf('run%d/%d...\n',run,nRuns)
            [r,F,P,stat(run)] = simulation(input(:,:,run),duration,dt,...
                                 Ilat_range,Imod_range,Itha_range,run==1,inputTag(:,run),[]);
            r_overall(:,run)=r.r_overall;
            r_pertile(:,:,run)=r.r_pertile;  
        end % nRuns
        
        save(resultfilename,'r_overall','r_pertile','tiles','F','P',...
                      'Ilat_range','Imod_range','Itha_range','speed',...
                      'inputfilename','resultfilename','stat');
    end
    
    plot_magri(mean(r_overall,2),mean(r_pertile,3),F,P,tiles,'fig6C.png');
    trialidx=1;plotOn=1;
    simulation(input(:,:,trialidx),duration,dt,Ilat_range,Imod_range,Itha_range,plotOn,inputTag(:,trialidx),'fig6A.png');
    plot_spect_stat(stat);

end
%==========================================================================
% input0: [time x 3 x nTrials]
function inputTag = input2tag(input0,Ilat_range,Imod_range,Itha_range)
    %-----K-means result------    
    load(fullfile('data','fig2_kmeans8.mat'));
    tag=reshape(tag,size(As));
    %-------rescale input------
    nTimepoints=size(input0,1);
    nTrials=size(input0,3);
    input=zeros(size(input0));
    for i=1:nTrials
        input(:,1,i)=rescale(input0(:,1,i),Ilat_range(1),Ilat_range(2))*2-5; 
        input(:,2,i)=rescale(input0(:,2,i),Imod_range(1),Imod_range(2))*2;
        input(:,3,i)=rescale(input0(:,3,i),Itha_range(1),Itha_range(2))*2;
    end
    %-------find closet input index------
    inputidx=round(input)+1;
    %-------return input tag------
    inputTag=zeros(nTimepoints,nTrials);
    for t=1:nTimepoints
        for i=1:nTrials
            inputTag(t,i)=tag(inputidx(t,1,i),inputidx(t,2,i),inputidx(t,3,i));
        end
    end
end
%==========================================================================
function plot_spect_stat(stat) 
    alphaz=[];
    gammaz=[];
    LFPpowerz=[];
    PLabel=[];
    simBOLDbpz=[];
    inputTag=[];
    for i=1:length(stat)
        alphaz=[alphaz;stat(i).alphaz];
        gammaz=[gammaz;stat(i).gammaz];
        LFPpowerz=[LFPpowerz;stat(i).LFPpowerz];
        PLabel=[PLabel;stat(i).PLabel];
        simBOLDbpz=[simBOLDbpz;stat(i).simBOLDbpz];
        inputTag=[inputTag;stat(i).inputTag];
    end

    % ----- plot hist2, ellipse, and line-----
    figure;
    t=tiledlayout(3,5,'TileSpacing','compact','Padding','loose');
    a=simBOLDbpz;
    b=alphaz;

    nexttile();
    xbins=-4:0.4:4;
    ybins=-4:0.4:4;
    histogram2(a,b,xbins,ybins,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
    %scatter(a,b,'k.')
    hold on;grid on;
    h= plot_ellipse(a,b);h.Color='w';
    climits=get(gca,'clim');
    set(gca,'yticklabel',[],'xticklabel',[],'fontsize',8,'FontName','calibri');
    Fit = polyfit(a,b,1);  
    plot(xlim,polyval(Fit,xlim),'w')
    text(-3.5,-3,sprintf('r= %g',round(corr(a,b),3)),'color','w','fontsize',8,'FontName','calibri');
    title('Overall')
    pos=get(gca,'position');
    colorbar('position',[0.28 pos(2)+0.02 0.01 0.15]);
    daspect([1 1 1]);

    tiles=0:10:100;
    for i=1:10
        nexttile(5+i);
        histogram2(a(PLabel==i),b(PLabel==i),...
                  xbins,ybins,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
        clim(climits/5)
        hold on;grid on;
        h= plot_ellipse(a,b);h.Color='w';
        h= plot_ellipse(a(PLabel==i),b(PLabel==i));h.Color='r';

        Fit = polyfit(a(PLabel==i),b(PLabel==i),1);  
        plot(xlim,polyval(Fit,xlim),'m')
        text(-3.5,-3,sprintf('r= %g',round(corr(a(PLabel==i),b(PLabel==i)),3)),'color','m','fontsize',8,'FontName','calibri')
        if ~ismember(i,[6])  
            set(gca,'yticklabel',[],'xticklabel',[]);
        else
            xlabel('BOLD');ylabel('Alpha');
        end
        set(gca,'fontsize',8,'FontName','calibri');
        title(sprintf('%d-%d%%',tiles([i,i+1])));
        daspect([1 1 1]);
    end
    pos=get(gca,'position');
    colorbar('position',[0.92 pos(2)+0.02 0.01 0.15]);
    width=12; 
    height=10;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig6D.png'));

    % -----plot cluster labels----- 
    figure;
    t=tiledlayout(3,5,'TileSpacing','compact','Padding','loose');
    nexttile();
    a=inputTag;
    h=histogram(a,0.5:8.5);ymax=sum(h.Values);
    hb=bar(1:8,h.Values/ymax);hb.FaceColor='flat';hb.CData=brewermap(8,'Paired');
    xlim([0 9])
    ylabel('Portion');ylimits=ylim;
    set(gca,'xticklabel',[],'fontsize',8,'FontName','calibri');
    title('Overall');axis('square')
    for i=1:10
        nexttile(5+i);
        h=histogram(a(PLabel==i),0.5:8.5);
        hb=bar(1:8,h.Values/ymax);hb.FaceColor='flat';hb.CData=brewermap(8,'Paired');  
        if ~ismember(i,[6])  
            set(gca,'yticklabel',[],'xticklabel',[]);
        else
            xlabel('Cluster ID');ylabel('Portion');
        end
        ylim([0 ylimits(2)/6]);
        set(gca,'fontsize',8,'FontName','calibri');
        title(sprintf('%d-%d%%',tiles([i,i+1])));
        axis('square')
        xlim([0 9])
    end
    width=12; 
    height=10;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig6E.png'));

    % -----plot hist2 of {BOLD, overal}x{alpha,gamma,overal}----- 
    figure;  
    t=tiledlayout(2,3,'TileSpacing','compact','Padding','loose');
    nexttile();
    cl=hist2_line(LFPpowerz,alphaz,'LFP','Alpha',[]);set(gca,'xticklabel',[],'yticklabel',[]);   
    nexttile();
    hist2_line(LFPpowerz,gammaz,'LFP','Gamma',cl); set(gca,'xticklabel',[],'yticklabel',[]);   
    nexttile();
    hist2_line(alphaz,gammaz,'Alpha','Gamma',cl); set(gca,'xticklabel',[],'yticklabel',[]);    
    nexttile();
    hist2_line(simBOLDbpz,LFPpowerz,'BOLD','LFP',cl);
    nexttile();
    hist2_line(simBOLDbpz,gammaz,'BOLD','Gamma',cl);set(gca,'xticklabel',[],'yticklabel',[]); 
    nexttile();
    hist2_line(simBOLDbpz,alphaz,'BOLD','Alpha',cl);set(gca,'xticklabel',[],'yticklabel',[]); 
    h= plot_ellipse(simBOLDbpz,alphaz);h.Color='w';
    pos=get(gca,'position');
    colorbar('position',[0.92 pos(2)+0.02 0.01 0.15]);

    width=8; 
    height=7;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig6B.png'));
    
end
%==========================================================================
function [r,F,P,stat]=simulation(input,duration,dt,Ilat_range,Imod_range,Itha_range,plotOn,inputTag,filename)

    tiles=0:10:100;
    dt2      = 1; 
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


    noiseLevel=1;
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
  
    %----- bandpass -----
    time_ds=0:dt2:duration;
    simBOLD_ds=interp1(time,simBOLD,time_ds);
    time_ds(isnan(simBOLD_ds))=[];
    simBOLD_ds(isnan(simBOLD_ds))=[];
    BPorder=2;
    simBOLD_bp = filter_bp(simBOLD_ds,[0.008 0.09],1/dt2,BPorder);  
    simBOLD_bp = interp1(time_ds,simBOLD_bp,time);

    %------ collect result ------
    [r.r_overall,r.r_pertile,F,P,stat]=plot_spectrogram(simEEG,simBOLD,simBOLD_bp,time,dt,time_discard,tiles,Iext,plotOn,inputTag,filename);

end
%==========================================================================
function [r_overall,r_pertile,F,P,stat]=plot_spectrogram(simEEG,simBOLD,simBOLDbp,time,dt,time_discard,tiles,Iext,plotOn,inputTag,filename)

    rmBase=true; % 1: remove baseline; 0: remove mean
    WINDOW=2000;
    mvstep=10;  
    delay=3.5; % EEG-BOLD latency (sec)
    frange=[0.5 100]; % Hz
    [S,F,T] = spectrogram(simEEG,WINDOW,WINDOW/2,WINDOW,1/dt);
    fidx=F>=frange(1) & F<=frange(2);
    tidx=find(T>time_discard(1) & T<time_discard(2));
    tmp=10*log(abs(S)); 
    
    if rmBase % remove baseline
        load(fullfile('data','fig5_baseline.mat'),'baseline');
        tmp=tmp-repmat(baseline,[1, size(tmp,2)]);    
    else  % remove mean
        tmp=tmp-repmat(mean(tmp,2),[1, size(tmp,2)]); 
    end
    Sm = movmean(tmp,mvstep,2);
    simBOLD0=interp1(time,simBOLDbp,T); 
    simBOLD2=interp1(time-delay,simBOLD,T);   
    simBOLDbp2=interp1(time-delay,simBOLDbp,T);
    T=T(tidx);
    F=F(fidx);
    Sm=Sm(fidx,tidx);
    LFPpowerz=zscore(sum(Sm)); 
    simBOLDz=zscore(simBOLD2(tidx));
    simBOLDbpz=zscore(simBOLDbp2(tidx));

    % ------ calcualte correlations ------
    r_overall=corr(Sm',simBOLDbpz');
    P=prctile(LFPpowerz,tiles);
    PLabel=zeros(size(LFPpowerz));
    for i=1:length(P)-1
        PLabel(LFPpowerz>=P(i) & LFPpowerz<P(i+1))=i;
    end
    r_pertile=zeros(size(Sm,1),length(P)-1);
    for i=1:length(P)-1
        r_pertile(:,i)=corr(Sm(:,PLabel==i)',simBOLDbpz(PLabel==i)');
    end
    
    %-----collect spectral statistics-----   
    stat.alphaz=zscore(Sm(F==10,:))'; %zscored alpha;
    stat.gammaz=zscore(Sm(F==40,:))'; %zscored gamma;
    stat.LFPpowerz=LFPpowerz';
    stat.simBOLDbpz=simBOLDbpz'; % left shifted 3.5 msec
    stat.PLabel=PLabel';
    stat.inputTag=interp1(time',inputTag,T','nearest');

    if plotOn
        figure;
        t=tiledlayout(3,1,'TileSpacing','compact','Padding','loose');
        nexttile();%subplot(311)
        plot(time,Iext(1,:),'linewidth',1);hold on; grid on;
        plot(time,Iext(14,:),'linewidth',1); 
        plot(time,Iext(2,:),'linewidth',1); legend('I_{lat}','I_{mod}','I_{tha}');
        plot(xlim,[0;0],'k','linewidth',2)
        xlim(time_discard);ylim([-5 25]);ylabel('Rate (Hz)')
        k=8;
        imagesc(time,-4.5:-.5,repmat(inputTag',5,1));clim([0,k]+.5);
        colormap(gca,brewermap(k,'Paired'));
        set(gca,'fontsize',8,'FontName','calibri',...
            'xticklabel',[],'ytick',0:10:25,'yticklabel',0:10:25)
        title('Input')
        legend('I_{lat}','I_{mod}','I_{tha}','location','northeastoutside');

        nexttile();
        imagesc(T,F,Sm);set(gca,'ydir','normal','xticklabel',[]);
        if rmBase  % remove baseline
            clim([-1 1]*25);
        else % remove mean
            clim([-1 1]*max(abs(Sm(:)))*0.7);
        end
        xlim(time_discard); ylim([0.5 50]);
        ylabel('Freq (Hz)');
        hold on; plot(xlim,[1;1]*[10,40],'k:')
        set(gca,'fontsize',8,'FontName','calibri',...
                'ytick',0:10:50,'yticklabel',0:10:50)
        title('Spectrogram'); 
        colormap(gca,jet);
        pos=get(gca,'position');
        colorbar('position',[0.92 pos(2)+0.02 0.01 0.15]); 

        nexttile();
        plot(T,LFPpowerz,'linewidth',1); hold on;grid on
        plot(T,zscore(simBOLD0(tidx)),'linewidth',1);
        set(gca,'fontsize',8,'FontName','calibri');
        title('LFP power & BOLD')
        xlim(time_discard);
        legend({'LFP pw.','BOLD'},'location','northeastoutside');
        ylabel('z-score');
        xlabel('Time (sec)');
        width=13; 
        height=10;
        set(gcf,'units','centimeters','position',[2 2 width height])
        if ~isempty(filename)
            saveas(gcf,fullfile('figures',filename));
        end
    end
end

%==========================================================================
% load random-walk input
function [input,duration,dt,speed,nTrials,inputfilename]=load_input(mvwindow)

    speed    = 5; 
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

%==========================================================================
function plot_magri(r_overall,r_pertile,F,P,tiles,filename)
  
    figure;
    width=8; 
    height=10;
    set(gcf,'units','centimeters','position',[2 2 width height])

    flimits=[F(2) F(end)]; 
    subplot(2,10,1);
    imagesc(1,F,r_overall);ylim(flimits);
    set(gca,'fontsize',8,'FontName','calibri','xtick',[]);
    ylabel('Frequency (Hz)');
    title('Overall');   
    clim([-1 1]*max(abs(r_pertile(:))));

    subplot(2,10,2:8);
    imagesc(1:length(P)-1,F,r_pertile); ylim(flimits); colormap(coolwarm);
    clim([-1 1]*max(abs(r_pertile(:))));
    set(gca,'xtick',(1:length(tiles))-0.5,...
            'xticklabel',sprintfc('%g%%',tiles),'yticklabel',[],...
            'fontsize',8,'FontName','calibri')
    xlabel('LFP power - percentile');
    title('Fixed LFP power');
    colorbar('SouthOutside','position',[0.08 0.5 0.2 0.02]);

    subplot(2,10,[9 10]); 
    plot(F,r_overall,'r','linewidth',1);hold on;
    plot(F,mean(r_pertile,2),'k','linewidth',2);   
    plot(xlim,[1;1]*0,'k--'); title('Average')
    plot([1;1]*[10,40],ylim,'k:'); title('Average')
    camroll(-90);
    xlim(flimits);
    ylimit=ylim;
    set(gca,'xticklabel',[],'YAxisLocation','right',...
        'ytick',[0 round(ylimit(2),2)/2],'yticklabel',[0 round(ylimit(2),2)/2],...
        'fontsize',8,'FontName','calibri');

    if ~isempty(filename)
        saveas(gcf,fullfile('figures',filename));
    end

end

%==========================================================================
function climits=hist2_line(a,b,labelA,labelB,climits)
    xbins=-4:0.4:4;
    ybins=-4:0.4:4;
    histogram2(a,b,xbins,ybins,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');  
    hold on;
    Fit = polyfit(a,b,1);plot(xlim,polyval(Fit,xlim),'w');
    text(-3.5,-3,sprintf('r= %g',round(corr(a,b),3)),'color','w','fontsize',8,'FontName','calibri')
    xlabel(labelA);ylabel(labelB);
    set(gca,'fontsize',8,'FontName','calibri'); 
    if isempty(climits)
        climits=clim;
    else
        clim(climits);
    end
    daspect([1 1 1]);
end
