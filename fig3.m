
function fig3()

    load(fullfile('data','fig2_batch_result.mat'))
    % rangeA [1x5]: I_lat (lateral input strength)
    % rangeB [1x51]: I_mod (modulatory input strength)
    % rangeC [1x11]: I_tha (thalamic input strength)
    % As [5x51x11]: I_lat of 2805 InputConditions
    % Bs [5x51x11]: I_mod of 2805 InputConditions
    % Cs [5x51x11]: I_tha of 2805 InputConditions
    % all_spect [2805x2501]: power spectrum (in dB)
    % spect_baseline [1x2501]: baseline power spectrum (in dB)
    % freq [1x2501]
    % all_alpha_env_hrf_bp: [2805x570] time series of one run (t=31-600 sec, sampling 1Hz)
    % all_alpha_env_hrf_ds: [2805x570] time series of one run (t=31-600 sec, sampling 1Hz)
    % all_gamma_env_hrf_bp: [2805x570] time series of one run (t=31-600 sec, sampling 1Hz)
    % all_gamma_env_hrf_ds: [2805x570] time series of one run (t=31-600 sec, sampling 1Hz)
    % all_simBOLD_bp: [2805x570] time series of one run (t=31-600 sec, sampling 1Hz)
    % all_simBOLD_ds: [2805x570] time series of one run (t=31-600 sec, sampling 1Hz)    
    % all_r_bp [2805x100]: [nInputConditions x nRuns] alpha-bold correlation (bandpass)
    % all_r_ds [2805x100]: [nInputConditions x nRuns] alpha-bold correlation (without bandpass)

    % ------ crop ------ 
    selectA=(5:9)/2;
    all_spect=crop(all_spect,rangeA,rangeB,rangeC,selectA);
    all_BOLD=crop(all_simBOLD_ds,rangeA,rangeB,rangeC,selectA);
    As=As(ismember(rangeA,selectA),:,:);
    Bs=Bs(ismember(rangeA,selectA),:,:);
    Cs=Cs(ismember(rangeA,selectA),:,:);
    rangeA=rangeA(ismember(rangeA,selectA));

    % ------ plot ------     
    frange=[30,50;...  % gamma
            13,30;...  % beta
            8,12; ...  % alpha
            4,7;...    % theta
            0.5,4];    % delta

    load(fullfile('data','fig2_kmeans8.mat'),'k','tag'); 
    plot_fig4_S1ac(all_BOLD,rangeA,rangeB,rangeC,As,k,tag);
    plot_fig4_S1b(all_spect,freq,all_BOLD,rangeA,As,k,tag,frange);
    plot_fig4_S1d(all_spect,freq,all_BOLD,rangeA,As,k,tag,frange);

end
%==========================================================================
function data=crop(data,rangeA,rangeB,rangeC,selectA)
    data=reshape(data,length(rangeA),length(rangeB),length(rangeC),[]);
    data=data(ismember(rangeA,selectA),:,:,:);  
    rangeA=rangeA(ismember(rangeA,selectA));
    data=reshape(data,length(rangeA)*length(rangeB)*length(rangeC),[]);  
end

%==========================================================================
function plot_fig4_S1ac(all_BOLD,rangeA,rangeB,rangeC,As,k,tag)

    mean_BOLD=mean(all_BOLD,2);
    tmp0=cell(2,1);
    tmp0{1}=reshape(mean_BOLD,[length(rangeA),length(rangeB),length(rangeC)]);
    std_BOLD=std(all_BOLD,[],2);
    tmp0{2}=reshape(std_BOLD,[length(rangeA),length(rangeB),length(rangeC)]);
    
    figure;
    t = tiledlayout(3,length(rangeA),'TileSpacing','compact','Padding','loose');
    %-----kmeans-----
    idx=reshape(tag,size(As)); % sorted clusters
    for i=1:length(rangeA)
        nexttile
        imagesc(rangeB,rangeC,squeeze(idx(i,:,:))');clim([0 k]+0.5);
        set(gca,'ydir','normal','fontsize',8,'FontName','calibri')
        if i~=1,set(gca,'yticklabel',[]);end
        title(sprintf('I_{lat} = %g',rangeA(i)));
        colormap(gca,brewermap(k,'Paired'));
    end
    pos=get(gca,'position');
    colorbar('position',[0.92 pos(2)+0.04 0.02 0.2],'ytick',1:k,'yticklabel',1:k);
    
    %-----BOLD-----
    for f=1:2
        tmp=tmp0{f}; 
        tmp=tmp/max(tmp(:));
        for i=1:length(rangeA)
            nexttile;
            imagesc(rangeB,rangeC,squeeze(tmp(i,:,:))');clim([0 1]);
            set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri')
            colormap(gca,'parula');
            if i~=1,set(gca,'yticklabel',[]);end
            if i==length(rangeA)
                pos=get(gca,'position');
                h=colorbar('position',[0.92 pos(2)+0.04 0.02 0.15]);
            end 
        end
    end
    xlabel(t,'I_{mod}','fontsize',10,'FontName', 'calibri')
    ylabel(t,'I_{tha}','fontsize',10,'FontName', 'calibri')
    width=10; 
    height=8;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig3_S1AC.png')); 
end

%==========================================================================
function plot_fig4_S1b(all_spect,freq,all_BOLD,rangeA,As,k,tag,frange)

    feature=zeros(numel(As),5);
    for f=1:size(frange,1)
        feature(:,f)=mean(all_spect(:,freq>=frange(f,1) & freq<=frange(f,2)),2);
    end

    figure;
    t=tiledlayout(5,5,'TileSpacing','compact','Padding','loose');
    ftext={'Gamma','Beta','Alpha','Theta','Delta'};
    tmp=mean(all_BOLD,2); % mean BOLD
    tmp = rescale(tmp,0,1);
    cmap=brewermap(k,'Paired');
    for f=1:5 % freq
        tmp2=rescale(feature(:,f),0,1);
        for a=1:length(rangeA) % A
           nexttile;
           for i=1:k
              scatter(tmp(tag==i & As(:)==rangeA(a)),tmp2(tag==i & As(:)==rangeA(a)),2,cmap(i,:),'filled');hold on;  
           end
           xlim([0 1]);
           ylim([0 1]);
           grid on;box on;
           set(gca,'xticklabel',[],'yticklabel',[],'fontsize',8,'FontName', 'calibri');
           if a==1,ylabel(ftext{f});end
           if f==1,title(sprintf('I_{lat} = %g',rangeA(a)));end
        end
        
    end
    xlabel(t,'BOLD amplitude (mean)','fontsize',10,'FontName', 'calibri')
    width=10; 
    height=10;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig3_S1B.png')); 
end

%==========================================================================
function plot_fig4_S1d(all_spect,freq,all_BOLD,rangeA,As,k,tag,frange)

    feature=zeros(numel(As),5);
    for f=1:size(frange,1)
        feature(:,f)=mean(all_spect(:,freq>=frange(f,1) & freq<=frange(f,2)),2);
    end
     
    figure;
    t=tiledlayout(5,5,'TileSpacing','compact','Padding','loose');
    ftext={'Gamma','Beta','Alpha','Theta','Delta'};
    tmp=std(all_BOLD,[],2); 
    tmp = rescale(tmp,0,1);
    cmap=brewermap(k,'Paired');
    for f=1:5 % freq
        tmp2=rescale(feature(:,f),0,1);
        for a=1:length(rangeA) % A
           nexttile;
           for i=1:k
              scatter(tmp(tag==i & As(:)==rangeA(a)),tmp2(tag==i & As(:)==rangeA(a)),2,cmap(i,:),'filled');hold on;  
           end
           xlim([0 1]);
           ylim([0 1]);
           grid on;box on;
           set(gca,'xticklabel',[],'yticklabel',[],'fontsize',8,'FontName', 'calibri');
           if a==1,ylabel(ftext{f});end
           if f==1,title(sprintf('I_{lat} = %g',rangeA(a)));end
        end
     
    end
    xlabel(t,'BOLD amplitude (std)','fontsize',10,'FontName', 'calibri')
    width=10; 
    height=10;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig3_S1D.png')); 
end


