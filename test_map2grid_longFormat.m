close all
clear
clc

figure_num = 0;

%% constant finding
%% [6]      %% [3.7]
%1 % 0:8    
            % 1,2 % 5:13
%2 % 9:17                    % 1,2,3 % 5:22
            % 2,3 % 14:22
%3 % 18:26
            % 3,4 % 23:31
%4 % 27:35
            % 4,5 % 32:40
%5 % 36:44
            % 5,6 % 41:49
%6 % 45:53
            % 6,7 % 50:59
%7 % 54:62
            % 7,8 % 60:63
%% [10]
%8 % 63 (some error)
%%      

% update peak for ea
UE_preIdxs = [60];                                        
array_timeDomainSignal = [];
snr = 0;

tic
for preIdx_index = 1:length(UE_preIdxs)

    % prachConfig.PrachConfigurationIndex = randi([67 250], 1, 1);
    % prachConfig.RootSequenceIndex = randi([0 137], 1, 1);

    prachConfig.PrachConfigurationIndex = 27;
    prachConfig.RootSequenceIndex = 39;

    prachConfig.PreambleIndex = UE_preIdxs(preIdx_index);

    prachConfig.SubcarrierSpacing = 1.25;  % delta_f_RA
    prachConfig.Set = 'Unrestricted';
    % prachConfig.ZeroCorrelationZoneConfig = randi([0 15], 1, 1);
    prachConfig.ZeroCorrelationZoneConfig = 8;
    prachConfig.FrequencyRange = 'FR1';
    prachConfig.SpectrumType = 'Unpaired';

    prachConfig.PrachFreqStart = 0;
    fprintf('PrachConfigurationIndex %d\n', prachConfig.PrachConfigurationIndex);
    fprintf('RootSequenceIndex %d\n', prachConfig.RootSequenceIndex);
    fprintf('ZeroCorrelationZoneConfig %d\n', prachConfig.ZeroCorrelationZoneConfig);
    

    %% Carrier Config
    % N_grid_size_mui, clause 5.3.2, TS 38.211
    carrierConfig.n_UL_RB = 273;
    
    carrierConfig.SubcarrierSpacing = 30;
    carrierConfig.numElementPerResourceBlock = 12;
    carrierConfig.numFrame = 2;
    
    %%
    PrachConfigFR1UnpairedSpectrum = get_Table6332x(prachConfig, carrierConfig);
    
    %% get N_CS
    N_CS = get_N_CS(prachConfig, PrachConfigFR1UnpairedSpectrum);

    %% Physical root sequence number (u)
    [u, u_arr] = get_u(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

    %% Cyclic shift
    [C_v, C_v_arr] = get_C_v(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

    %% zadoffChu sequence
    x_u = zadoffChu(u, PrachConfigFR1UnpairedSpectrum.L_RA);
    x_uv = circshift(x_u, [0 -C_v]);
    y_uv = fft(x_uv);

    %% Modulation
    grid = map2grid_longFormat(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
    figure(1)
    plotResourceGrid(grid);
    % timeDomain_signal = PRACH_modulation_main(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
    % 
    % %% tdl config, add noise
    % signalOut_TDL_whithoutNoise = timeDomain_signal;
    % % sum_timeDomainSignal = timeDomain_signal;
    % % 
    % % tdl = nrTDLChannel;
    % % 
    % % tdl.SampleRate = prachConfig.SubcarrierSpacing * 1000 * 4096;
    % % tdl.TransmissionDirection = 'Uplink';
    % % tdl.MaximumDopplerShift = 100;  % 100 Hz
    % % tdl.DelayProfile = 'TDL-C';
    % % tdl.DelaySpread = 300e-9;
    % % tdl.NumTransmitAntennas = 1;
    % % tdl.NumReceiveAntennas = 4;
    % % 
    % % signalOut_TDL = tdl(sum_timeDomainSignal);
    % % signalOut_TDL = awgn(signalOut_TDL, snr, 'measured');
    % % signalOut_TDL = gpuArray(signalOut_TDL);
    % 
    % %% Demodulation
    % isPlot = 0;
    % isPrintResult = 1;
    % [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr] = PRACH_demodulation_main(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, signalOut_TDL_whithoutNoise, isPlot, isPrintResult);
    % figure(1)
    % plotResourceGrid(grid);
    % % figure(2)
    % % plot(abs(timeDomain_signal));
    % % % plot(timeDomain_signal);
end