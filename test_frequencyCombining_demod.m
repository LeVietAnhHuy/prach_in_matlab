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
UE_preIdxs = [0:63];
array_timeDomainSignal = [];
constant = 27;
snr_range = flip(-40:-29);
gain = 8;

% Table 8.4.1.1-1 in TS 38.141-1v15.03, page 166
timeErrorTolerance = 0.26; % micro second, for AWGNs
% timeErrorTolerance = 1.77; % micro second, for TDLC300-100

tic
for snr = snr_range
    numDetectedPreambleIndex = 0;
    for preIdx_index = 1:length(UE_preIdxs)

        % prachConfig.PrachConfigurationIndex = randi([67 250], 1, 1);
        prachConfig.PrachConfigurationIndex = 158;

        prachConfig.RootSequenceIndex = randi([0 137], 1, 1);
        % prachConfig.RootSequenceIndex = 39;

        prachConfig.PreambleIndex = UE_preIdxs(preIdx_index);

        prachConfig.SubcarrierSpacing = 30;  % delta_f_RA
        prachConfig.Set = 'Unrestricted';
        prachConfig.ZeroCorrelationZoneConfig = randi([0 15], 1, 1);
        % prachConfig.ZeroCorrelationZoneConfig = 15;
        prachConfig.FrequencyRange = 'FR1';
        prachConfig.SpectrumType = 'Unpaired';

        prachConfig.PrachFreqStart = 0;
        % fprintf('PrachConfigurationIndex %d\n', prachConfig.PrachConfigurationIndex);
        % fprintf('RootSequenceIndex %d\n', prachConfig.RootSequenceIndex);
        % fprintf('ZeroCorrelationZoneConfig %d\n', prachConfig.ZeroCorrelationZoneConfig);


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
        % grid = map2grid_longFormat(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
        grid = mapPrachSymbol2ResourceGrid_main(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);

        timeDomain_signal = PRACH_modulation_main(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
        timeDomain_signal = signalScale(timeDomain_signal, gain);
        %% Testing procedure in Clause 8.4.1.4  in TS 38.141-1v15.03
        %% Part 1: Delay profiles
        %% Step 1: AWGN
        % Table 8.4.1.4.2-1 in TS 38.141-1v15.03
        % snr_dbm = 100;
        % snr_dbW = snr_dbm - 30;
        snr_dbW = snr;
        timeDomain_signal_AWGN = awgn(timeDomain_signal, snr_dbW, 'measured');

        %% Step 3: Multipath Fading
        TDLC300_100 = nrTDLChannel("DelayProfile", "TDLC300", ...
            "MaximumDopplerShift", 100, ...
            "SampleRate", prachConfig.SubcarrierSpacing * 1000 * 4096, ...
            "TransmissionDirection", "Uplink", ...
            "NumTransmitAntennas", 1, ...
            "NumReceiveAntennas", 1);
        %% Received Signal
        timeDomain_signal_TDL = TDLC300_100(timeDomain_signal);
        timeDomain_signal_AWGN_TDL = TDLC300_100(timeDomain_signal_AWGN);

        %% Plot
        % figure(1)
        % plotResourceGrid(grid);
        % figure(2)
        % plot(abs(timeDomain_signal));
        % figure(3)
        % plot(abs(timeDomain_signal_AWGN));
        % figure(4)
        % plot(abs(timeDomain_signal_TDL));
        % figure(5)
        % plot(abs(timeDomain_signal_AWGN_TDL));

        %% Demodulation
        isPlot = 0;
        isPrintResult = 0;
        %% Frequency Contiguous Combining
        %% without noise
        % [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal, isPlot, isPrintResult);
        %% with AWGN 
        [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal_AWGN, isPlot, isPrintResult);
        %% with Multipath Fading 
        % [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal_TDL, isPlot, isPrintResult);
        %% with AWGN + Multipath Fading 
        % [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal_AWGN_TDL, isPlot, isPrintResult);
        

        %% Frequency Random Combining
        %% without noise
        % [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyRandomCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal, isPlot, isPrintResult);
        %% with AWGN 
        % [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyRandomCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal_AWGN, isPlot, isPrintResult);
        %% with Multipath Fading 
        % [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyRandomCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal_TDL, isPlot, isPrintResult);
        %% with AWGN + Multipath Fading 
        % [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyRandomCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal_AWGN_TDL, isPlot, isPrintResult);

        if ~isempty(detected_preambleIndex_arr)
            detectedPreambleIndex = mode(detected_preambleIndex_arr);
            timeErrorTolerance = mode(nTA_microSecond_arr);

            if detectedPreambleIndex == UE_preIdxs(preIdx_index)
                numDetectedPreambleIndex = numDetectedPreambleIndex + 1;
            else
                disp(detectedPreambleIndex);
                disp(UE_preIdxs(preIdx_index))
            end
        end
    end
    acc = (numDetectedPreambleIndex / 64) * 100;
    fprintf('\n%ddB: %d%%', snr, acc);
end
