function [ sta ] = compute_sta( stimulus, spiketimes, numSAMPLES, sta_lwe , dt_STIM)
%COMPUTE_STA Calculates the spike-triggered average for a neuron that
%           is driven by a stimulus defined in stim. The spike-
%           triggered average is computed over num_timesteps timesteps.
%           inputs: 
%                   sta_lwe - STA left window edge (5 in Segev's  implementation)
%  Description:     To do this, compute the average of all of the vectors
%           starting 300 ms (exclusive) before a spike and ending at the time of
%           the event (inclusive). Each of these vectors defines a list of
%           samples that is contained within a window of 300 ms before the each
%           spike. The average of these vectors should be completed in an
%           element-wise manner.
% inspired by - Computational NS Coursera class
    
    % This command finds the indices of all of the spikes that occur
    % after numSAMPLES samples into the recording and before the end of the stimulus
    
    spike_SAMPLES = spiketimes/dt_STIM;
    spike_SAMPLES_rel = round(   spike_SAMPLES(  spike_SAMPLES > numSAMPLES &...         % all spikes after a STA size window (windowSize)
        ( spike_SAMPLES + numSAMPLES < length(stimulus) )  )   );                    % all spikes after the last sample of stimulus - windowSize
    
    % simple spike-triggered average computation of the spikes found using the find command.
    sta_x = arrayfun(@(x) linspace( x-(numSAMPLES-sta_lwe-1), x + sta_lwe, numSAMPLES ), spike_SAMPLES_rel,'UniformOutput',false);
    spikesArrays = cellfun(@(x) stimulus(x),sta_x,'UniformOutput',false);        
%     sta = mean(cell2mat(spikesArrays'),2);
    sta = 1 * mean(cell2mat(spikesArrays'),2); % the linear kernel should have a negative sign (?ASK SHY?)
%     disp('finished one sta esimation');
end

