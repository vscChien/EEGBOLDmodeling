% >> fig7_batch(); % run simulations and save result.
% >> fig7fig8();   % plot result.
%
function fig7fig8()

    %---------------------------------------
    % load power spectrum (rate and psp)
    % all_pop_spect 2805x17x2501
    % all_psp_spect 2805x17x2501
    % all_rate_mean 2805x17
    % all_rate_std 2805x17
    % all_spect 2805x2501
    % freq 1x2501
    % rangeA [5,6,7,8,9]/2
    % rangeB (0:50)/2;
    % rangeC (0:10)/2;
    % baseline
    resultname=fullfile('data','fig7_batch_result.mat');
    fprintf('Loading %s...',resultname);
    load(resultname);
    disp('done.')
    bl=baseline;
    %---------------------------------------
    % load EEG-BOLD correlation (constant input condition) 
    % all_r_gammaEnvHrfbp_simBOLDbp 2805 x 100
    % all_r_alphaEnvHrfbp_simBOLDbp 2805 x 100
    load(fullfile('data','fig4_batch_result.mat'));
    r_alphaBOLD=mean(all_r_alphaEnvHrfbp_simBOLDbp,2); % 2805 x 1
    r_gammaBOLD=mean(all_r_gammaEnvHrfbp_simBOLDbp,2); % 2805 x 1
    %---------------------------------------
    load(fullfile('data','fig2_kmeans8.mat'),'tag','k'); 
    
    %--------calculate spatial correlation---------
    r_EEG_PSP = spatial_correlation(all_spect,all_psp_spect,bl,freq,tag,k);
    % --------plot-------------
    plot_r(r_EEG_PSP,all_spect,bl.spect,freq,tag,k,'fig7.png'); 
    plot_rr(r_alphaBOLD,all_pop_spect,freq,all_rate_mean,all_rate_std,[-1 1]*0.6,'fig8_alpha-BOLD.png');
    plot_rr(r_gammaBOLD,all_pop_spect,freq,all_rate_mean,all_rate_std,[-1 1]*0.7,'fig8_gamma-BOLD.png');
    plot_1x17_map(all_rate_std,rangeA,rangeB,rangeC,'fig8S1_std.png'); 
    plot_1x17_map(all_rate_mean,rangeA,rangeB,rangeC,'fig8S1_mean.png'); 
    plot_5x17_map(all_psp_spect,bl.spect_psp,freq,rangeA,rangeB,rangeC,'fig8S1_psp.png'); 
    plot_5x17_map(all_pop_spect,bl.spect_pop,freq,rangeA,rangeB,rangeC,'fig8S1_rate.png'); 

end



%==========================================================================
% spect [N x freq]; N conditions in tagi
function [spect,freq]=get_spect(all_spect,baseline_spect,freq,tag,tagi)
    spect=all_spect-repmat(baseline_spect',[size(all_spect,1),1]); % remove baseline
    spect=spect(:,freq>=0 & freq<=50);  
    spect=spect(tag==tagi,:);
    freq=freq(freq>=0 & freq<=50);
end
%==========================================================================
% spect [N x 17 x freq]
function [spect,freq]=get_psp_spect(all_psp_spect,baseline_spect0,freq,tag,tagi)
    baseline_spect=permute(baseline_spect0,[3 2 1]);
    spect=all_psp_spect-repmat(baseline_spect,[size(all_psp_spect,1),1]); % remove baseline
    spect=spect(:,:,freq>=0 & freq<=50);
    spect=spect(tag==tagi,:,:);
    freq=freq(freq>=0 & freq<=50);  
end

%==========================================================================
function plot_r(r,all_spect,baseline_spect,freq,tag,k,filename)

    yticklabeltext={'E2/3','E4','E5ET','E5IT','E6','P2/3','P4','P5','P6','S2/3','S4','S5','S6','V2/3','V4','V5','V6'};

    idx=[1,6,10,14,...
         2,7,11,15,...
         3,4,8,12,16,...
         5,9,13,17];
    figure;
    t=tiledlayout(6,4,'TileSpacing','compact','Padding','loose');
    for tagi=1:k
        %--------spectrum----
        if tagi<=4, nexttile(tagi);else,nexttile(tagi+8);end
        [target,fs] = get_spect(all_spect,baseline_spect,freq,tag,tagi); % [N x freq]
        plot(fs,target,'color',[1,1,1]*0.7); hold on; 
        plot(fs,mean(target),'k');
        title(sprintf('Cluster %d',tagi));xlim([0 50])
        
        plot([1;1]*[4,8,12,30],ylim,'k'); 
        plot(xlim,[0;0],'k');
        if ismember(tagi,[1,5]),ylabel('Power (dB)');end  
        set(gca,'xtick',[2.5,5.5,10,20,40],'xticklabel',...
                 {'$\delta$','$\theta$','$\alpha$','$\beta$','$\gamma$'},'TickLabelInterpreter', 'latex',...
                'fontsize',8,'FontName', 'calibri')
        %--------correlation----
        if tagi<=4, nexttile(tagi+4,[2 1]);else,nexttile(tagi+12,[2 1]);end
        imagesc(fs,1:17,r(idx,:,tagi));clim([-1 1]);
        set(gca,'ytick',1:17,'yticklabel',yticklabeltext(idx),...
             'xtick',0:10:50,'xticklabel',0:10:50,'fontsize',8,'FontName', 'calibri');
        xlabel('Frequency (Hz)');
        if ismember(tagi,[1,5]),ylabel('Populations');end
        hold on; plot([1;1]*[4,8,12,30],ylim,'k'); 
        plot(xlim,[1;1]*[4.5,8.5,13.5],'k'); 

    end
    colormap(coolwarm);
    colorbar('position',[0.92 0.12 0.01 0.1]);
    width=20; 
    height=20;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures',filename));
end

%==========================================================================
function plot_rr(r_alphaBOLD,all_pop_spect,freq,all_rate_mean,all_rate_std,ylimits,filename)

    xticklabeltext={'E2/3','E4','E5ET','E5IT','E6','P2/3','P4','P5','P6','S2/3','S4','S5','S6','V2/3','V4','V5','V6'};
    idx=[1,6,10,14,...
         2,7,11,15,...
         3,4,8,12,16,...
         5,9,13,17];

    r1 = corr(r_alphaBOLD,all_rate_mean(:,idx));
    r2 = corr(r_alphaBOLD,all_rate_std(:,idx));
    r3 = zeros(17,2501);
    for i=1:17
        r3(i,:) = corr(r_alphaBOLD,squeeze(all_pop_spect(:,idx(i),:)));
    end
    figure;
    t=tiledlayout(3,1,'TileSpacing','compact','Padding','loose');
    nexttile
    bar([r1;r2]');xlim([0.5 17.5]);ylim(ylimits)
    ylabel('Correlation');
    hold on;plot([1;1]*[4.5,8.5,13.5],ylim,'k'); 
    legend({'r_{mean}','r_{std}'},'location','northoutside')
    set(gca,'xtick',1:17,'xticklabel',xticklabeltext(idx),'fontsize',8,'FontName', 'calibri');
 
    nexttile([2 1])
    imagesc(1:17,freq,r3');ylim([0 50]);clim([-1 1]*max(abs(r3(:))));
    set(gca,'ytick',[2.5,5.5,10,20,40],'yticklabel',...
                 {'$\delta$','$\theta$','$\alpha$','$\beta$','$\gamma$'},'TickLabelInterpreter', 'latex',...
                 'xtick',1:17,'xticklabel',xticklabeltext(idx),'fontsize',8,'FontName', 'calibri');
    hold on; plot(xlim,[1;1]*[4,8,12,30],'k'); 
    plot([1;1]*[4.5,8.5,13.5],ylim,'k'); 
    colormap(coolwarm);
    colorbar('position',[0.92 0.12 0.02 0.2]);
    width=10; 
    height=10;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures',filename));

end
%==========================================================================
function plot_5x17_map(data0,baseline0,freq,rangeA,rangeB,rangeC,filename)
  
    ttext={'E2/3','E4','E5ET','E5IT','E6','PV2/3','PV4','PV5','PV6','SOM2/3','SOM4','SOM5','SOM6','VIP2/3','VIP4','VIP5','VIP6'};   
    idx=[1,6,10,14,...
         2,7,11,15,...
         3,4,8,12,16,...
         5,9,13,17];

    frange=[30,50;...  % gamma
            13,30;...  % beta
            8,12; ...  % alpha
            4,7;...    % theta
            0.5,4];    % delta
    baseline=permute(baseline0,[3 2 1]); % 2805 x 17 x 2501
    data=data0-repmat(baseline,[size(data0,1),1]); % remove baseline
    figure;
    t = tiledlayout(size(frange,1),17,'TileSpacing','compact','Padding','loose');
    width=24; 
    height=8;
    set(gcf,'units','centimeters','position',[2 2 width height])
    for f=1:size(frange,1) 
        fband=mean(data(:,:,freq>=frange(f,1) & freq<=frange(f,2)),3);
        fband=reshape(fband,[length(rangeA),length(rangeB),length(rangeC),17]);
        tmp=squeeze(fband(end,:,:,:)); %[51x11x17]
        for pop=1:17
            nexttile()
            tmp=squeeze(fband(end,:,:,idx(pop))); %[51x11] 
            imagesc(rangeB,rangeC,tmp'); 
            text(1,2,sprintf('%g',round(max(tmp(:)),2)),...
                'color','w','fontsize',7,'FontName','calibri');
            text(1,0.5,sprintf('%g',round(min(tmp(:)),2)),...
                'color','w','fontsize',7,'FontName','calibri');

            set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri')
            axis square
            if f==1, title(ttext{idx(pop)});end
            if f~=size(frange,1) || pop~=1
                set(gca,'xticklabel',[],'yticklabel',[]);
            else
                xlabel('I_{mod}','fontsize',8,'FontName', 'calibri')
                ylabel('I_{tha}','fontsize',8,'FontName', 'calibri')
            end
        end 
        
    end
    pos=get(gca,'position');
    hc=colorbar('position',[pos(1)+pos(3)+0.005 pos(2) 0.005 0.08],'ytick',[]); title(hc,'dB')
    saveas(gcf,fullfile('figures',filename));
end

%==========================================================================
function plot_1x17_map(data0,rangeA,rangeB,rangeC,filename)

    data=reshape(data0,length(rangeA),length(rangeB),length(rangeC),[]); 
    data=squeeze(data(rangeA==4.5,:,:,:)); 
    ttext={'E2/3','E4','E5ET','E5IT','E6','PV2/3','PV4','PV5','PV6','SOM2/3','SOM4','SOM5','SOM6','VIP2/3','VIP4','VIP5','VIP6'};
    idx=[1,6,10,14,...
         2,7,11,15,...
         3,4,8,12,16,...
         5,9,13,17];


    figure;
    t=tiledlayout(1,17,'TileSpacing','compact','Padding','loose');
    width=24; 
    height=8;
    set(gcf,'units','centimeters','position',[2 2 width height])
    
    for pop=1:17

        nexttile(); 
        tmp=data(:,:,idx(pop))';
        imagesc(rangeB,rangeC,tmp);hold on; 
        text(1,2,sprintf('%g',round(max(tmp(:)),2)),...
                'color','w','fontsize',7,'FontName','calibri');
        text(1,0.5,sprintf('%g',round(min(tmp(:)),2)),...
                'color','w','fontsize',7,'FontName','calibri');
        set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri')
        if pop~=1
            set(gca,'yticklabel',[],'xticklabel',[]);
        else
            xlabel('I_{mod}','fontsize',8,'FontName', 'calibri')
            ylabel('I_{tha}','fontsize',8,'FontName', 'calibri')
        end           
        title(ttext{idx(pop)})
        axis('square')

        pos=get(gca,'position');
    end
    pos=get(gca,'position');
    hc=colorbar('position',[pos(1)+pos(3)+0.005 pos(2) 0.005 0.08],'ytick',[]); title(hc,'Hz')
    saveas(gcf,fullfile('figures',filename));
end

%==========================================================================
function r = spatial_correlation(all_spect,all_psp_spect,bl,freq,tag,k)

    %----- spatial correlation (spect_EEG,spect_PSP)-----
    r = zeros(17,251,8); %[pops x freqs x clusters]
    for tagi=1:k         
        [target,fs] = get_spect(all_spect,bl.spect,freq,tag,tagi);        % spect_EEG  
        regressors = get_psp_spect(all_psp_spect,bl.spect_psp,freq,tag,tagi); % spect_PSP
        for p=1:17
            for f=1:length(fs)
                y=target(:,f);
                x=regressors(:,p,f);
                r(p,f,tagi)=corr(y(:),x(:));
            end
        end  
    end
end