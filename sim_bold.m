% input
%    z: [timepoints x 1] synaptic activity
% output
%    y: [timepoints x 1] BOLD signal
%
% Also check sim_bold_ode.m for a same implementation
%
% [ref]  
% https://doi.org/10.1006/nimg.2000.0630
% https://doi.org/10.1016/S1053-8119(03)00202-7
%    
function y=sim_bold(z, plotOn,dt)

    if nargin==1
        plotOn=0;
    end

    if nargin<1
        plotOn=1;
        % -----input----- 
        z=zeros(30*1000,1); 
        z((1:100)+1000)=1;  

    end
    if nargin<3
        dt      = 1e-3; % s
    end
    
    
    z_t=(1:length(z))*dt;


    switch 2
        case 1 % Friston2000
            % -----parameters----- 
            % rCBF component (linear)
            epsi    = 0.5;
            taus    = 0.8;  % s
            tauf    = 0.4;  % s
            % Ballon component (nonlinear)
            E0      = 0.8;  % resting oxygen extraction fraction
            tau0    = 1;    % mean transit time
            alpha0  = 0.2;  % stiffness exponent
            V0      = 0.02; % resting blood volume fraction
        case 2 % https://github.com/neurolib-dev/neurolib/blob/master/neurolib/models/bold/timeIntegration.py
            % -----parameters----- 
            % rCBF component (linear)
            epsi    = 1e-3;%0.5;
            taus    = 1/0.65;  %1.5385 sec
            tauf    = 1/0.41;  %2.4390 sec
            % Ballon component (nonlinear)
            E0      = 0.34;    % resting oxygen extraction fraction
            tau0    = 0.98;    % mean transit time
            alpha0  = 0.32;    % stiffness exponent
            V0      = 0.02;    % resting blood volume fraction
    end
    u=z;              % synaptic activity
    s=zeros(size(u)); % vasodilatory signal
    f=ones(size(u));  % inflow (rCBF)
    v=ones(size(u));  % blood volume
    q=ones(size(u));  % deoxyhemoglobin content
        

    for t=1:length(z)-1

        ds = epsi*u(t) - s(t)/taus - (f(t)-1)/tauf;
        df = s(t);

        dv = ( f(t) - v(t)^(1/alpha0) )/tau0;
        E  = 1-(1-E0)^(1/f(t));
        dq = ( f(t)*E/E0 - v(t)^(1/alpha0)*q(t)/v(t) )/tau0;


        s(t+1) = s(t)+ds*dt;
        f(t+1) = f(t)+df*dt;
        v(t+1) = v(t)+dv*dt;
        q(t+1) = q(t)+dq*dt;

    end

    % -----BOLD signal-----
    k1 = 7*E0;
    k2 = 2;
    k3 = 2*E0 - 0.2;
    y = V0*(k1*(1-q) + k2*(1-q./v) + k3*(1-v));


    % -----plot-----
    if plotOn
        figure;
        subplot(231);plot(z_t,z); title('input'); xlabel('Time (sec)');grid on;
        subplot(232);plot(z_t,s); title('s'); xlabel('Time (sec)');grid on;
        subplot(233);plot(z_t,[v',q']); xlabel('Time (sec)');  title('v & q'); legend('v','q');grid on;
        subplot(235);plot(z_t,f); xlabel('Time (sec)'); title('f');grid on;
        subplot(236);plot(z_t,y); xlabel('Time (sec)'); title('BOLD');grid on;

        subplot(234);
        fs=1/dt;
        data=y(20000:end); 
        window_length=length(data);
        [spects,freqs]=compute_spectrum(data,fs,window_length); 
        plot(freqs,spects,'linewidth',1); xlim([0 0.3]); xlabel('Hz');ylabel('power');title('BOLD Spectrum');
    end


end
