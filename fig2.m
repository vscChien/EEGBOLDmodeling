% >> fig2_batch(); % run simulations and save result. (will take a couple hours)
% >> fig2();       % plot fig.2
%
function fig2()

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
    
    plot_fig2a(all_spect,spect_baseline,freq,rangeA,rangeB,rangeC,frange);   
    plot_fig2bcd(all_spect,spect_baseline,freq,rangeA,rangeB,rangeC,As,Bs,Cs,frange);
    plot_fig2_S1(all_spect,freq,As,frange);

    
end
%==========================================================================
function data=crop(data,rangeA,rangeB,rangeC,selectA)
    data=reshape(data,length(rangeA),length(rangeB),length(rangeC),[]);
    data=data(ismember(rangeA,selectA),:,:,:);  
    rangeA=rangeA(ismember(rangeA,selectA));
    data=reshape(data,length(rangeA)*length(rangeB)*length(rangeC),[]);  
end
%==========================================================================
function plot_fig2a(all_spect,spect_baseline,freq,rangeA,rangeB,rangeC,frange)

    all_spect_rb=all_spect-repmat(spect_baseline,[size(all_spect,1),1]); % remove baseline
    
    figure;
    t = tiledlayout(size(frange,1),length(rangeA),'TileSpacing','compact','Padding','loose');
    for f=1:size(frange,1)
    
        fband=mean(all_spect_rb(:,freq>=frange(f,1) & freq<=frange(f,2)),2);
        fband=reshape(fband,[length(rangeA),length(rangeB),length(rangeC)]);
        
        climit=[min(fband(:)),max(fband(:))];
        for i=1:length(rangeA)
            nexttile()
            imagesc(rangeB,rangeC,squeeze(fband(i,:,:))');clim(climit);
            set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri')

            if f==1, title(sprintf('I_{lat} = %g',rangeA(i)));end
            if i~=1
                set(gca,'yticklabel',[]);
            end
            if i==length(rangeA)
                pos=get(gca,'position');
                h=colorbar('position',[0.92 pos(2)+0.01 0.02 0.12]);
                h.Ticks = [max(0,ceil(climit(1))),floor(climit(2))] ; 
                h.TickLabels = [max(0,ceil(climit(1))),floor(climit(2))];
            end 
        end    
    end
    xlabel(t,'I_{mod}','fontsize',10,'FontName', 'calibri')
    ylabel(t,'I_{tha}','fontsize',10,'FontName', 'calibri')

    width=10; 
    height=12;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig2A.png')); 
 
end
%==========================================================================
function plot_fig2bcd(all_spect,spect_baseline,freq,rangeA,rangeB,rangeC,As,Bs,Cs,frange)

    feature=zeros(numel(As),5);
    for f=1:size(frange,1)
        feature(:,f)=mean(all_spect(:,freq>=frange(f,1) & freq<=frange(f,2)),2);
    end
    % --------Fig2 B and C---------
    k=8; 
    tmp0=zscore(all_spect(:,freq>=0.5 & freq<=50),1); % for display
    F=freq(freq>=0.5 & freq<=50);
    tmp=zscore(feature); % for k-means
    rng(2); 
    [tag0,~,~]= kmeans(tmp0,k);
    kfeatures=zeros(k,5);
    d=vecnorm([As(:),Bs(:),Cs(:)/2],2,2);
    dfeatures=zeros(k,1); % a cluster i's distance to zero in input space. 
    for i=1:k
      kfeatures(i,:)=mean(tmp(tag0==i,:),1);
      dfeatures(i,:)=min(d(tag0==i,:));
    end
    %-----reorder labels-----
    tag=tag0;
    if k==8
        kidx=[5,2,6,8,7,4,1,3]; % manual sort for k==8, rng(2);
    else
        [kfeatures2,kidx]=sortrows(dfeatures); % sort by input conditions
    end
    for i=1:k
        tag(tag0==kidx(i))=i;
    end
    if k==8 && ~exist(fullfile('data','fig2_kmeans8.mat'),'file')
        save(fullfile('data','fig2_kmeans8.mat'),...
            'k','tag','rangeA','rangeB','rangeC','As'); 
    end 

    % --------Fig2CD---------
    figure
    tiledlayout(3,ceil(k/2),'TileSpacing','compact','Padding','loose');
    [tag_sorted,i]=sort(tag);
    nexttile([1 ceil(k/2)])
    T=1:numel(As);
    imagesc(T,F,tmp0(i,:)');set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri');
    hold on;plot([1;1]*(find(diff(tag_sorted))'+0.5),ylim,'k')
    tmp=find(diff(tag_sorted))';
    set(gca,'xtick',[0,tmp]+diff([0,tmp,T(end)])/2,'xticklabel',1:k);
    ylim([F(1) F(end)]);clim([-1 1]*3);colormap(jet);
    hold on;plot(xlim,[1;1]*[10,40],'k:');
    title('Clusters of power spectra');ylabel('Freq (Hz)');
    set(gca,'ytick',[10,40],'yticklabel',[10,40]);
    pos=get(gca,'position');
    colorbar('position',[0.92 pos(2)+0.045 0.02 0.1]);

    for i=1:k
        nexttile;
        plot(freq,spect_baseline,'k','linewidth',1.5); hold on;
        plot(freq,mean(all_spect(tag==i,:),1),'r','linewidth',1.5);
        title(sprintf('Cluster %d',i));
        set(gca,'xtick',[10,40],'xticklabel',[10,40]);
        xlim([0.5 50]);grid on;
        ylim([-55 -25])
        if i==1 || i==floor(k/2)+1
            ylabel('Power (dB)');
        else
            set(gca,'yticklabel',[]);
        end
        if i>floor(k/2)
            xlabel('Freq (Hz)');
        end
        set(gca,'fontsize',8,'FontName', 'calibri');
    end
    width=10; 
    height=9;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig2CD.png')); 
    
    %--------Fig2B---------
    figure;
    t = tiledlayout(1,length(rangeA),'TileSpacing','compact','Padding','loose');
    idx=reshape(tag,size(As)); 
    for i=1:length(rangeA)
        nexttile
        imagesc(rangeB,rangeC,squeeze(idx(i,:,:))');clim([0 k]+0.5);
        set(gca,'ydir','normal','fontsize',8,'FontName','calibri')
        if i~=1
            set(gca,'yticklabel',[]);
        end
        title(sprintf('I_{lat} = %g',rangeA(i)));
    end
    xlabel(t,'I_{mod}','fontsize',10,'FontName', 'calibri');
    ylabel(t,'I_{tha}','fontsize',10,'FontName', 'calibri');
    
    pos=get(gca,'position');
    colorbar('position',[0.92 0.4 0.02 0.4],'ytick',1:k,'yticklabel',1:k);
    colormap(brewermap(k,'Paired'));
    width=10; 
    height=3.8;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig2B.png')); 
    
end
%==========================================================================
function plot_fig2_S1(all_spect,freq,As,frange)

    feature=zeros(numel(As),5);
    for f=1:size(frange,1)
        feature(:,f)=mean(all_spect(:,freq>=frange(f,1) & freq<=frange(f,2)),2);
    end
    Ks=20;
    all_sumd=zeros(1,Ks);
    uniqueSignedFeatures=zeros(1,Ks);
    tmp0=zscore(all_spect(:,freq>=0.5 & freq<=50)); % for k-means
    tmp=zscore(feature); 
    for k=1:Ks
        rng(2)
        [idx,~,sumd]= kmeans(tmp0,k);
        all_sumd(k)=sum(sumd);
        %-----calculate uniqueSignedFeatures
        kfeatures=zeros(k,5);
        for i=1:k
            kfeatures(i,:)=sign(mean(tmp(idx==i,:),1));
        end
        uniqueSignedFeatures(k)=size(unique(kfeatures,'rows'),1);
    end
    figure;plot(all_sumd,'-o');grid on;
    xlabel('Number of clusters'); ylabel('Sums of point-to-centroid distances'); 
    title('K-means')
    yyaxis right;
    plot(uniqueSignedFeatures,'-o');ylabel('Number of unique signed features');
    hold on;plot([1;15],[1;15],':')
    set(gca,'fontsize',10,'FontName', 'calibri')
    width=10; 
    height=8;
    set(gcf,'units','centimeters','position',[2 2 width height])
    saveas(gcf,fullfile('figures','fig2_S1.png')); 

end