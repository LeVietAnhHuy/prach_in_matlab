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

snr_arr = [-50:0];
performance_arr = [];
preIdx_arr = 0:63;
epoch = 1;
tic

for snr = snr_arr
    performance = 0;
    prog = 0;
    fprintf('Case %d/%d (%2ddB): %3d%%\n', epoch, length(snr_arr), snr, prog);
    for preIdx = preIdx_arr
        prachConfig.PrachConfigurationIndex = 158;
        prachConfig.RootSequenceIndex = 39;
        prachConfig.L_RA = 139;

        prachConfig.PreambleIndex = preIdx;

        prachConfig.SubcarrierSpacing = 30;  % delta_f_RA
        prachConfig.Set = 'Unrestricted';
        prachConfig.ZeroCorrelationZoneConfig = 8;
        prachConfig.FrequencyRange = 'FR1';
        prachConfig.SpectrumType = 'Unpaired';

        prachConfig.PrachFreqStart = 0;

        %% Carrier Config
        % N_grid_size_mui, clause 5.3.2, TS 38.211
        carrierConfig.n_UL_RB = 273;

        carrierConfig.SubcarrierSpacing = 30;
        carrierConfig.numElementPerResourceBlock = 12;
        carrierConfig.numFrame = 1;

        %% get N_CS
        N_CS = get_N_CS(prachConfig);

        %% Physical root sequence number (u)
        [u, u_arr] = get_u(prachConfig, N_CS);

        %% Cyclic shift
        [C_v, C_v_arr] = get_C_v(prachConfig, N_CS);

        %% zadoffChu sequence
        x_u = zadoffChu(u, prachConfig.L_RA);
        x_uv = circshift(x_u, [0 -C_v]);
        y_uv = fft(x_uv);

        %%
        PrachConfigFR1UnpairedSpectrum = get_Table6332x(prachConfig, carrierConfig);

        %% Modulation
        grid = mapPrachSymbol2ResourceGrid(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
        timeDomain_signal = PRACH_modulation(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

        %% Scale signal
        sum_timeDomainSignal = signalScale(timeDomain_signal);

        %% tdl config, add noise
        signalOut_TDL_whithoutNoise = timeDomain_signal;

        tdl = nrTDLChannel;

        tdl.SampleRate = prachConfig.SubcarrierSpacing * 1000 * 4096;
        tdl.TransmissionDirection = 'Uplink';
        tdl.MaximumDopplerShift = 100;  % 100 Hz
        tdl.DelayProfile = 'TDL-C';
        tdl.DelaySpread = 300e-9;
        tdl.NumTransmitAntennas = 1;
        tdl.NumReceiveAntennas = 4;

        signalOut_TDL = tdl(sum_timeDomainSignal);
        signalOut_TDL = awgn(signalOut_TDL, snr, 'measured');
        signalOut_TDL = gpuArray(signalOut_TDL);

        %% Demodulation
        [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, output_preambleIndex_arr] = PRACH_demodulation_8(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, signalOut_TDL);
        if length(output_preambleIndex_arr) == 0;
            output_preambleIndex_arr = -100;
        end
        % disp(output_preambleIndex_arr);
        numCorrect = 0;
        for preIdx_idx = 1:length(output_preambleIndex_arr)
            if(preIdx == output_preambleIndex_arr(preIdx_idx))
                numCorrect = numCorrect + 1;
            end
        end
        performance = performance + numCorrect / length(output_preambleIndex_arr);
        prog = ((preIdx + 1) / length(preIdx_arr)) * 100;
        % disp(prog);
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b%3.0f%%', prog);
        
    end
    performance_arr = [performance_arr; performance / 3];
    epoch = epoch + 1;
    fprintf('\n');

end

disp("Writing results to csv file...")
output_path = "evaluation_results.xlsx";
evaluation_results = table(snr_arr.', performance_arr);
writetable(evaluation_results, ...
           output_path, ...
           'Range','D1');

disp("Done!");
toc

