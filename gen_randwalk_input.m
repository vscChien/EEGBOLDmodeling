function input = gen_randwalk_input(duration,dt,speed)
    input= ones(1,duration/dt);
    walk = randn(1,duration/dt)*speed*dt; 
    input(1)=rand; % initial [0,1]
    for t=1:length(input)-1
       tmp=input(t)+walk(t);
       tmp(tmp>1)=1;tmp(tmp<0)=0;
       input(t+1)=tmp;
    end
end