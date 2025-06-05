
function out=filter_bp(data,cut,fs,order)
    if nargin<4
        order=4; % Filter order (adjust as needed)
    end

    lowcut = cut(1);  % Lower cutoff frequency in Hz
    highcut = cut(2); % Upper cutoff frequency in Hz
    
    % Normalize the cutoff frequencies with respect to the Nyquist frequency
    nyquist = fs / 2;
    lowcut = lowcut / nyquist;
    highcut = highcut / nyquist;
    
    % Butterworth bandpass filter
    %order = 4;  % Filter order (adjust as needed)
    [b, a] = butter(order, [lowcut, highcut], 'bandpass');
    
    % Apply the bandpass filter to the BOLD signal
    out = filtfilt(b, a, data);
    
end