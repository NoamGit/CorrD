function [ outProcess ] = LowPassAcorr( obj, inProcess, cutoff_freq )
%[ outProcess ] = LowPassAcorr( obj, inProcess, cutoff_freq ) this function 
%applies low pass filtering with specified filter properties to
%Autocorrelation function. The filter is embedded in the function

% define FIR Window lowpass filter for a preproceesing step of the estimation
    fc = cutoff_freq;
    Wn = (2/obj.Fs)*fc; 
    if 2*obj.maxlags <= 300 % data must be more than 3 times the filter order
        ord = floor(obj.maxlags*2/3)-1;
    else
             ord = 100; % default max filter order
    end
    
    filt_coeff = fir1(ord,Wn,'low',kaiser(ord+1,1));
    outProcess = filtfilt(filt_coeff, 1, inProcess); %  filter data with provided FIR filter
end

