% >> for i=1:10,fig4_batch(10);end  % run 100 simulations and save result. (take a couple days)
% >> fig4();       % plot fig4
% >> fig4_S1(1);   % plot fig4_S1A
% >> fig4_S1(2);   % plot Fig4_S1B

resultname=fullfile('data','fig4_batch_result.mat');
if ~exist(resultname,'file')
    filename=fullfile('data','fig4_batch_result_*'); 
    files=dir(filename);
    all_r_alphaEnvHrfds_simBOLDds=[];
    all_r_alphaEnvHrfbp_simBOLDbp=[];
    all_r_gammaEnvHrfds_simBOLDds=[];
    all_r_gammaEnvHrfbp_simBOLDbp=[];
    all_r_alpha_x=[];
    all_r_gamma_x=[];
    all_r_alphaEnv_x=[];
    all_r_gammaEnv_x=[];
    all_alphaPower=[];
    all_gammaPower=[];
    for i=1:length(files)
        tmp=load(fullfile('data',files(i).name));
        all_r_alphaEnvHrfds_simBOLDds=[all_r_alphaEnvHrfds_simBOLDds,tmp.all_r_alphaEnvHrfds_simBOLDds];
        all_r_alphaEnvHrfbp_simBOLDbp=[all_r_alphaEnvHrfbp_simBOLDbp,tmp.all_r_alphaEnvHrfbp_simBOLDbp];
        all_r_gammaEnvHrfds_simBOLDds=[all_r_gammaEnvHrfds_simBOLDds,tmp.all_r_gammaEnvHrfds_simBOLDds];
        all_r_gammaEnvHrfbp_simBOLDbp=[all_r_gammaEnvHrfbp_simBOLDbp,tmp.all_r_gammaEnvHrfbp_simBOLDbp];
        all_r_alpha_x=[all_r_alpha_x,tmp.all_r_alpha_x];
        all_r_gamma_x=[all_r_gamma_x,tmp.all_r_gamma_x];
        all_r_alphaEnv_x=[all_r_alphaEnv_x,tmp.all_r_alphaEnv_x];
        all_r_gammaEnv_x=[all_r_gammaEnv_x,tmp.all_r_gammaEnv_x];
        all_alphaPower=[all_alphaPower,tmp.all_alphaPower];
        all_gammaPower=[all_gammaPower,tmp.all_gammaPower];
    end
    rangeA=tmp.rangeA;
    rangeB=tmp.rangeB;
    rangeC=tmp.rangeC;
    As=tmp.As;
    Bs=tmp.Bs;
    Cs=tmp.Cs;
    dipConf=[-1,-0.1,1,1,0.1];
    noiseLevel=1;
    disp('-------')
    fprintf('input dimension: [%d, %d, %d]\n',size(As))
    fprintf('dipole setting: [%g, %g, %g, %g, %g]\n',dipConf)
    fprintf('noise level: %g\n',noiseLevel)
    disp('-------')
    save(resultname,'all_r_alphaEnvHrfds_simBOLDds','all_r_alphaEnvHrfbp_simBOLDbp',...
        'all_r_gammaEnvHrfds_simBOLDds','all_r_gammaEnvHrfbp_simBOLDbp',...
        'all_r_alpha_x','all_r_gamma_x','all_r_alphaEnv_x','all_r_gammaEnv_x',...
        'all_alphaPower','all_gammaPower','rangeA','rangeB','rangeC',...
        'As','Bs','Cs','noiseLevel','dipConf')
else
    fprintf('Loading %s...',resultname);
    load(resultname);
    fprintf('done.\n')
end
plot_correlation_map(all_r_alphaEnvHrfbp_simBOLDbp,As,rangeA,rangeB,rangeC,[-1 1]*.06,'fig4A.png')
plot_correlation_map(all_r_gammaEnvHrfbp_simBOLDbp,As,rangeA,rangeB,rangeC,[-1 1]*.2,'fig4B.png')
plot_histogram(all_r_alphaEnvHrfbp_simBOLDbp,all_r_gammaEnvHrfbp_simBOLDbp,...
               As,rangeA,rangeB,rangeC,'fig4C.png'); 

%==========================================================================
function plot_correlation_map(data,As,rangeA,rangeB,rangeC,climits,filename)
    r_mean=reshape(mean(data,2),size(As));
    r_std=reshape(std(data,[],2),size(As));

    figure;
    t=tiledlayout(2,5,'TileSpacing','compact','Padding','loose');
    width=10; 
    height=6;
    set(gcf,'units','centimeters','position',[2 2 width height])
    
    if isempty(climits)
        climits=[-1 1]*max(abs(r_mean(:))); 
    end
    for i=1:5 
        nexttile; 
        h=imagesc(rangeB,rangeC,squeeze(r_mean(i,:,:))');clim(climits); hold on;
        set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri')
        if i~=1,set(gca,'yticklabel',[]);end
        title(sprintf('I_{lat}=%g',rangeA(i)))
        if i==5, scatter(15,0.5,"g^",'filled');end 
        if i==5, scatter(14,4,"y^",'filled');end
        axis('square')
    end
    pos=get(gca,'position');
    colormap(coolwarm);
    colorbar('position',[0.92 pos(2)+0.02 0.01 0.15]);
    
    % mean/std
    tmp=r_mean./r_std;
    climits=[-1 1]*max(abs(tmp(:)));
    for i=1:5 
        nexttile; 
        h=imagesc(rangeB,rangeC,squeeze(tmp(i,:,:))');clim(climits);
        if i~=1,set(gca,'yticklabel',[]);end
        set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri')
        axis('square')
    end
    pos=get(gca,'position');
    colormap(coolwarm);
    colorbar('position',[0.92 pos(2)+0.02 0.01 0.15]);
    xlabel(t,'I_{mod}','fontsize',10,'FontName', 'calibri')
    ylabel(t,'I_{tha}','fontsize',10,'FontName', 'calibri')
    saveas(gcf,fullfile('figures',filename));
end

%==========================================================================
function plot_histogram(all_r_alphaEnvHrfbp_simBOLDbp,all_r_gammaEnvHrfbp_simBOLDbp,As,rangeA,rangeB,rangeC,filename)
        examples=[4.5,17,0.5;...
                  4.5,14,4];
        ex(1)=sub2ind(size(As),find(rangeA==examples(1,1)),find(rangeB==examples(1,2)),find(rangeC==examples(1,3)));
        ex(2)=sub2ind(size(As),find(rangeA==examples(2,1)),find(rangeB==examples(2,2)),find(rangeC==examples(2,3)));
        figure;
        t=tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        for i=1:2
            nexttile();
            histogram(all_r_alphaEnvHrfbp_simBOLDbp(ex(i),:),-0.5:0.05:0.5,'Normalization','probability');hold on;
            histogram(all_r_gammaEnvHrfbp_simBOLDbp(ex(i),:),-0.5:0.05:0.5,'Normalization','probability');
            xlim([-1 1]*0.5);
            ylim([0 0.35]);
            plot([0;0],ylim,'k');
            hold on;
            mean_r_alpha=mean(all_r_alphaEnvHrfbp_simBOLDbp(ex(i),:));
            mean_r_gamma=mean(all_r_gammaEnvHrfbp_simBOLDbp(ex(i),:));
            plot([1;1]*mean_r_alpha,ylim,'b','linewidth',1.5);
            plot([1;1]*mean_r_gamma,ylim,'r','linewidth',1.5);
            title(sprintf('Case %d: [I_{lat}, I_{mod}, I_{tha}] = [%g, %g, %g]',i,examples(i,:)))
            text(round(mean_r_alpha,2)-0.22,0.32,sprintf('$r_{\\alpha}$=%g',round(mean_r_alpha,2)),...
                'color','b','fontsize',8,'FontName','calibri','Interpreter','latex');
            text(round(mean_r_gamma,2)+.02,0.32,sprintf('$r_{\\gamma}$=%g',round(mean_r_gamma,2)),...
                'color','r','fontsize',8,'FontName','calibri','Interpreter','latex');
            set(gca,'fontsize',8,'FontName','calibri')
        end      
        xlabel(t,'Correlation','fontsize',10,'FontName', 'calibri'); 
        ylabel(t,'Probability','fontsize',10,'FontName', 'calibri');
        legend({'Alpha-BOLD','Gamma-BOLD'},'NumColumns', 1,'Location','southoutside'...
                ,'fontsize',8,'FontName', 'calibri')
        width=10; 
        height=12;
        set(gcf,'units','centimeters','position',[2 2 width height])
        saveas(gcf,fullfile('figures',filename));
end
