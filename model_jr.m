function [rate, time, PSP] = model_jr(W, tau, sigmParam, Iext, noiseLevel, dt)
 
    nP=size(W,1); 
    T=size(Iext,2);
    time = (1:T)*dt;

         
    rate = zeros(nP,T);
    PSP  = zeros(nP,nP+1,T);   
    
    u=zeros(nP,nP+1); 
    v=zeros(nP,nP+1);   

    for t = 1:T    

        % PSP to rate
        r = sigm(sum(u,2),sigmParam);      

        % save rate and PSP     
        rate(:,t) = r;  
        PSP(:,:,t)= u;

        % update potential
        m=[repmat(r',[nP 1]),Iext(:,t)];
        du = v;
        dv = (W.*m + noiseLevel*randn(size(W))/sqrt(dt))./tau - 2*v./tau - u./tau.^2;
        v = v + dv*dt;
        u = u + du*dt;

    end 
    
end

%==========================================================================
function r=sigm(x,param)
    a=param(:,1);
    b=param(:,2);
    c=param(:,3);
    r=c./(1+exp(-a.*(x-b))); 
end