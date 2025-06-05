function sigmParam=get_sigmParam(nE,nPV,nSOM,nVIP)
    speedOn=1;
    plotOn=0;
    [sigmE,sigmPV,sigmSOM,sigmVIP] = get_sigmoid(speedOn,plotOn); 
    sigmParam=[repmat([sigmE.a,sigmE.b,sigmE.c],      [nE,1]);    % E
               repmat([sigmPV.a,sigmPV.b,sigmPV.c],   [nPV,1]);   % PV
               repmat([sigmSOM.a,sigmSOM.b,sigmSOM.c],[nSOM,1]);  % SOM
               repmat([sigmVIP.a,sigmVIP.b,sigmVIP.c],[nVIP,1])]; % VIP
end