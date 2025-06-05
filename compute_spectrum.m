%  Input
%    data: [timepoints x ch]
%  Output
%    spects: [F x ch], where F= [window_length/2]+1 
%    freqs: [F x 1]

function [dB,freqs,power]=compute_spectrum(data,fs,window_length)

    [power, freqs] = pwelch(data, window_length, 0, window_length, fs);
    dB=10*log10(power);
    
end