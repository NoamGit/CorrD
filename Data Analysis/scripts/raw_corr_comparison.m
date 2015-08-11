%% Here we compare raw data - STA ac and GCP ac
% we compare preprocessed data to see if visulally, there is any chance for
% fittnes of the signals

% Load data
D = LoadDataIntoProcess( 1 ); % loads Segev's data
process = D.Data;

% find CGP ac only left-side without peak
kernNAME = 'gaussSine';           
% maxlags = ceil((2 * sqrt(2*log(10)) * 1)/process.dt); % maxlags are set according to the width of envelope 
maxlags = ceil((0.5 * sqrt(2*log(10)) * 1)/process.dt); % we divide by by 4 because of the high sampling rate 
process.maxlags = maxlags;
process = process.PreProcessCP( );
CGP_Corr_RAW = cell(process.numChannels,1);
for k = 1:process.numChannels
    CGP_Corr_RAW{k} = normax( process.acorr.CGP_Corr_RAW{k}(1:maxlags) );
end

% find STA ac only left-side without peak
dt_STIM = process.stimulus.dt;
Fs_STIM = 1/dt_STIM;
time_STIM = linspace(0,  length(process.stimulus.Yuncorr)*dt_STIM,  size(process.stimulus.Yuncorr,1));
STA_LENGTH = 35;
process = process.CalcSTA( STA_LENGTH, 5);
STA_Corr_RAW = cell(process.numChannels,1);
for k = 1:process.numChannels
    realSTA_rs = resample( process.STA(k).realSTA, process.Fs, Fs_STIM ); % realSTA up sampled to GCP AC resolution
    sta_corr = xcov( realSTA_rs, maxlags ,'unbiased' );
    STA_Corr_RAW{k} = normax( sta_corr(1:maxlags) );
end

% find  lambda ac only left-side without peak
L_Corr_RAW = cell(process.numChannels,1);
for k = 1:process.numChannels
    L_Corr_RAW{k} = normax( process.acorr.lambda_Corr_RAW{k}(1:maxlags) );
end

% compare results
timeaxis = fliplr(process.dt .* (1:maxlags));
num2disp = 6; % number of plots to display in a figure
numfigure = ceil(process.numChannels/6);
for k = 1:numfigure
    figure(k);
    for n = 1:num2disp
        dataNum = n+(k-1)*num2disp;
        try
            s(n) = subplot(3,2,n);
            plot(timeaxis, CGP_Corr_RAW{dataNum}, timeaxis, STA_Corr_RAW{dataNum});
            hold on ; plot(timeaxis, L_Corr_RAW{dataNum}, '--g'); hold off
        catch err
            break;
        end
        axis([0 timeaxis(1) -inf inf]);
        title(s(n),[num2str(dataNum),' cell']); xlabel('[sec]');ylabel('norm corr');
    end
    legend('CGP corr','STA corr','Rate corr');
end

%%  Here we try to find a better function f* s.t. CGP_Corr_RAW = f*(L_Corr_RAW)
% try only with good data

cellIndx = [3 6 11 12 13 17 18 19 20];
process = process.PreProcessCP('other', 8);
CGP_Corr_RAW_2 = cell(process.numChannels,1);
for k = 1:process.numChannels
    CGP_Corr_RAW_2{k} = normax( process.acorr.CGP_Corr_RAW{k}(1:maxlags) );
end

% compare results
timeaxis = fliplr(process.dt .* (1:maxlags));
num2disp = 4; % number of plots to display in a figure
numfigure = ceil(process.numChannels/6);
for k = 1:2
    figure(k);
    for n = 1:num2disp
        dataNum = cellIndx(n+(k-1)*num2disp);
        try
            s(n) = subplot(3,2,n);
            plot(timeaxis, CGP_Corr_RAW_2{dataNum},'--g', timeaxis, STA_Corr_RAW{dataNum},'r');
            hold on ; plot(timeaxis, CGP_Corr_RAW{dataNum}); hold off
        catch err
            break;
        end
        axis([0 timeaxis(1) -inf inf]);
        title(s(n),[num2str(dataNum),' cell']); xlabel('[sec]');ylabel('norm corr');
    end
    legend('CGP corr','STA corr');
end