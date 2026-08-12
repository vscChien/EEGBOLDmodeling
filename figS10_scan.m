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


function figS10_scan(cond)

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

    % Synaptic time constants
    tau = get_tau(nP,nE,nPV,nSOM,nVIP);
            
    % Sigmoid functions
    sigmParam = get_sigmParam(nE,nPV,nSOM,nVIP);

    Ts=20;%[1,5,10,20];

    switch cond
        case 0 % no stim
            label='null';
            EPSV=[0,0,0,0]; % optical stim strength
            filename='figS10_scan_null.mat'; 
            figname1='figS10_scan_null.png'; 
            figname2='figS10_scan_null.svg';
        case 1 % stim on PV
            label='PV2/3';
            EPSV=[0,2,0,0]; % optical stim strength
            filename='figS10_scan_PV.mat';
            figname1='figS10D.png';
            figname2='figS10D.svg';              
        case 2 % stim on SOM
            label='SOM2/3';
            EPSV=[0,0,0.3,0]; % optical stim strength
            filename='figS10_scan_SOM.mat';
            figname1='figS10E.png';
            figname2='figS10E.svg';                
        case 3 % stim on VIP
            label='VIP2/3';
            EPSV=[0,0,0,7]; % optical stim strength
            filename='figS10_scan_VIP.mat';
            figname1='figS10F.png';
            figname2='figS10F.svg';               
    end
    filename=fullfile('data',filename);

    if exist(filename,'file')
        fprintf('Loading result....%s.\n',filename)
        load(filename)
    else
        % external input
        nRuns=1;
        rangeA=[5,6,7,8,9]/2;
        rangeB=(0:50)/2;
        rangeC=(0:10)/2;
        [As,Bs,Cs] =ndgrid(rangeA,rangeB,rangeC);
        all_BOLDchanges  =zeros(numel(As),nRuns);
        for i=1:numel(As)
    
            fprintf('%d/%d...',i,numel(As))
            a=As(i);%4.5;% level of lateral 
            b=Bs(i);%5;  % level of modulatory
            c=Cs(i);%0;  % level of thalamic
            Iext0 = gen_Iext(a,b,c,duration,dt);
        
     
        
            for n=1:nRuns
                fprintf('%d, run %d\n',Ts,n)
                Tstart = 30; %sec
                Tend   = Tstart+Ts; %sec
                Istm = gen_Istm(duration,dt,Tstart,Tend,EPSV);
                Iext = Iext0+Istm;
    
                
                % simulation
                [~, time, PSP] = model_jr(W, tau, sigmParam, Iext, noiseLevel, dt);
            
                % simEEG and simBOLD
                OpticMode=1; % OpticaMode: Iext is excluded from the BOLD simulation.
                [~, simBOLD, ~] = sim_eeg_bold_opticmode(PSP,time,nE,nPV,nSOM,nVIP,nP,C,dt,OpticMode);
                    
                all_BOLDchanges(i,n) = cal_BOLD_change(simBOLD,time);
            end
        end % A
        save(filename,'all_BOLDchanges','As','Bs','Cs','rangeA','rangeB','rangeC');  
        fprintf('%s saved.\n',filename);
    end
    plot_result(all_BOLDchanges,rangeA,rangeB,rangeC);
    saveas(gcf,fullfile('figures',figname1));
    saveas(gcf,fullfile('figures',figname2));
end
%==========================================================================
function plot_result(all_BOLDchanges,rangeA,rangeB,rangeC)

    tmp=reshape(mean(all_BOLDchanges,2),[length(rangeA),length(rangeB),length(rangeC)]);
    climit=[-1,1]*0.3;%[-1,1]*max(abs(tmp(:)));%
    figure;
    t = tiledlayout(1,length(rangeA),'TileSpacing','compact','Padding','loose');
    width=10; 
    height=3.8;
    set(gcf,'units','centimeters','position',[2 2 width height])
    for i=1:length(rangeA)
            nexttile()
            imagesc(rangeB,rangeC,squeeze(tmp(i,:,:))');clim(climit);
            set(gca,'ydir','normal','fontsize',8,'FontName', 'calibri')


            title(sprintf('I_{lat} = %g',rangeA(i)));
            if i~=1
                set(gca,'yticklabel',[]);
            end
            if i==length(rangeA)
                pos=get(gca,'position');
                h=colorbar('position',[0.91 pos(2)+0.20 0.01 0.35],...
                            'xtick',[-0.3,0,0.3],'xticklabel',{'-30%','0%','30%'},...
                            'fontsize',7);
            end 
            colormap(coolwarm())
            %-----------label Iext=[4.5, 5, 0]-------------
            if rangeA(i)==4.5
                hold on;scatter(5,0,10,'k^')
            end
    end
    xlabel(t,'I_{mod}','fontsize',10,'FontName', 'calibri')
    ylabel(t,'I_{tha}','fontsize',10,'FontName', 'calibri')

end
%==========================================================================
function BOLDchange = cal_BOLD_change(simBOLD, time)

    % ================================
    % Inputs:
    % simBOLD: [1 x T]
    % time:     [1 x T]
    % ================================
    stim_onset    = 30; % stim: [30~50]s
    stim_offset   = 50; % 
    baseline_idx = time > (stim_onset-10) & time < stim_onset;  % 20~30  
    stimOn_idx     = time > (stim_offset-10) & time < (stim_offset-5); % 40~45

    % Convert BOLD to percent signal change for each run
    meanBOLD = mean(simBOLD(baseline_idx)); 
    BOLDchange =(simBOLD-meanBOLD)/meanBOLD;
    BOLDchange = mean(BOLDchange(stimOn_idx));
   
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