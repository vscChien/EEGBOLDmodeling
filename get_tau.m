
function tau = get_tau(nP,nE,nPV,nSOM,nVIP)
    tauE=6*1e-3;    % sec
    tauPV=3*1e-3;   % sec
    tauSOM=20*1e-3; % sec
    tauVIP=15*1e-3; % sec
    tauExt=6*1e-3;  % sec
    tau   = repmat([ones(1,nE)*tauE,ones(1,nPV)*tauPV,ones(1,nSOM)*tauSOM,ones(1,nVIP)*tauVIP,tauExt],[nP,1]);
end