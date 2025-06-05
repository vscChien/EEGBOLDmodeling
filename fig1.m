
load(fullfile('data','fig1_EmpDataVis.mat')); 
[~,idx]=sort(median(EmpDataVis.individual));
individual=EmpDataVis.individual(:,idx);


figure;
nbins=10;
histogram(EmpDataVis.meanMapMinVoxSubs,nbins);hold on;
histogram(EmpDataVis.meanMapMaxVoxSubs,nbins);

xlim([-1 1]*0.5);
ylim([0 18]);
plot([0;0],ylim,'k');

mean1=mean(EmpDataVis.meanMapMinVoxSubs);
mean2=mean(EmpDataVis.meanMapMaxVoxSubs);
plot([1;1]*mean1,ylim,'b','linewidth',1.5);
plot([1;1]*mean2,ylim,'r','linewidth',1.5);
text(round(mean1,2)-0.29,16,sprintf('r_1=%g',round(mean1,2)),...
                'color','b','fontsize',8,'FontName','calibri');
text(round(mean2,2)+.02,16,sprintf('r_2=%g',round(mean2,2)),...
                'color','r','fontsize',8,'FontName','calibri');
set(gca,'fontsize',8,'FontName','calibri')
xlabel('Alpha-BOLD correlation','fontsize',10,'FontName', 'calibri'); 
ylabel('Count','fontsize',10,'FontName', 'calibri');
title('Histogram')
width=5; 
height=5;
set(gcf,'units','centimeters','position',[2 2 width height])
saveas(gcf,fullfile('figures','fig1B.png')); 


%--------------------------------------------------------------------------
figure;
t=tiledlayout(1,15,'TileSpacing','compact','Padding','loose');
nexttile()
plot([0.5;1.5],[0;0],'k'); hold on;
m=median(EmpDataVis.meanMap);
high = quantile(EmpDataVis.meanMap,0.95);
low = quantile(EmpDataVis.meanMap,0.05);
high2=max(EmpDataVis.meanMap);
low2=min(EmpDataVis.meanMap);
for i=1
  fill([-.5,-.5,.5,.5]+i,[low(i),high(i),high(i),low(i)],[1 1 1]*0.8,'linestyle','none');
  plot([1;1]*i,[low2(i);high2(i)],'k-');
end
scatter(1,m,[],'r','_','linewidth',1.5);scatter(1,low2,[],'k','_');scatter(1,high2,[],'k','_');
xlim([0 2]);ylim([-1 1]*max(abs(individual(:)))*1.1);
plot(xlim,[0;0],'k');
title('Group');
set(gca,'xtick',[],'fontsize',8,'FontName','calibri')
ylabel('Alpha-BOLD correlation','fontsize',10,'FontName', 'calibri')

m=median(individual);
high = quantile(individual,0.95);
low = quantile(individual,0.05);
high2=max(individual);
low2=min(individual);
nexttile(2,[1 14])
hold on;
for i=1:72
  fill([-.5,-.5,.5,.5]+i,[low(i),high(i),high(i),low(i)],[1 1 1]*0.8,'linestyle','none');
  plot([1;1]*i,[low2(i);high2(i)],'k-');
end
xlim([0.5 72.5]);ylim([-1 1]*max(abs(individual(:)))*1.1);
scatter(1:72,m,[],'r','_','linewidth',1.5);scatter(1:72,low2,[],'k','_');scatter(1:72,high2,[],'k','_');
xlim([0 73])
plot(xlim,[0;0],'k');
set(gca,'ytick',[],'xtick',[1,67,70],'xticklabel',[1,67,70],'fontsize',8,'FontName','calibri')
title('Individual')
xlabel('Subject ID','fontsize',10,'FontName', 'calibri')
box on;
width=15; 
height=5;
set(gcf,'units','centimeters','position',[2 2 width height])
saveas(gcf,fullfile('figures','fig1C.png')); 


