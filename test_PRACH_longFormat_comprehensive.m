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
prachConfigIndex_arr = [0:66];

zeroCorrConfig_arr = [0:15];


array_timeDomainSignal = [];
snr = 0;

threshholdConstant_arr = [];
ifftCircShift_arr = [];
ConfigurationIndex_arr = [];
RootSequenceIndex_arr = [];
truePreambleIndex_arr = [];
detectedPreambleIndex_arr = [];
frame_arr = [];
subframe_arr = [];
slot_arr = [];
startingSymBol_arr = [];
TA_arr = [];

filename = strcat('test_longFormat_prachConfigIndex', num2str(prachConfigIndex_arr(1)), num2str(prachConfigIndex_arr(end)), ...
    '_zeroCorrConfig_', num2str(zeroCorrConfig_arr(1)), num2str(zeroCorrConfig_arr(end)), ...
    '.csv');
path_filename = strcat('test_results\', filename);
col_names = {'threshholdConstant_arr' ...
             'ifftCircShift_arr' ...
             'ConfigurationIndex_arr' ...
             'RootSequenceIndex_arr' ...
             'truePreambleIndex_arr' ...
             'detectedPreambleIndex_arr' ...
             'frame_arr' ...
             'slot_arr' ...
             'startingSymBol_arr' ...
             'TA_arr'};
writecell(col_names, path_filename);
% dlmwrite(path_filename, col_names);
% type(path_filename);
tic
prog = 0;
fprintf('Computation Progress: %3d%%\n', prog);
for preIdx_index = 1:length(UE_preIdxs)

    prog = (100*(preIdx_index/length(UE_preIdxs)));
    for prachConfigIndex = prachConfigIndex_arr
        for zeroCorrConfig = zeroCorrConfig_arr

            % prachConfig.PrachConfigurationIndex = randi([40 66], 1, 1);
            % prachConfig.RootSequenceIndex = randi([0 838], 1, 1);

            % prachConfig.PrachConfigurationIndex = 158;
            % prachConfig.RootSequenceIndex = 39;

            prachConfig.PrachConfigurationIndex = prachConfigIndex;
            prachConfig.RootSequenceIndex = randi([0 837], 1, 1);

            prachConfig.PreambleIndex = UE_preIdxs(preIdx_index);

            prachConfig.SubcarrierSpacing = 1.25;  % delta_f_RA
            if prachConfigIndex > 39
                prachConfig.SubcarrierSpacing = 5;  % delta_f_RA
            end

            prachConfig.Set = 'Unrestricted';
            % prachConfig.ZeroCorrelationZoneConfig = randi([0 15], 1, 1);
            % prachConfig.ZeroCorrelationZoneConfig = 15;
            prachConfig.ZeroCorrelationZoneConfig = zeroCorrConfig;
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
            grid = mapPrachSymbol2ResourceGrid(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
            timeDomain_signal = PRACH_modulation(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

            %% tdl config, add noise
            signalOut_TDL_whithoutNoise = timeDomain_signal;
            % sum_timeDomainSignal = timeDomain_signal;
            % 
            % tdl = nrTDLChannel;
            % 
            % tdl.SampleRate = prachConfig.SubcarrierSpacing * 1000 * 4096;
            % tdl.TransmissionDirection = 'Uplink';
            % tdl.MaximumDopplerShift = 100;  % 100 Hz
            % tdl.DelayProfile = 'TDL-C';
            % tdl.DelaySpread = 300e-9;
            % tdl.NumTransmitAntennas = 1;
            % tdl.NumReceiveAntennas = 4;
            % 
            % signalOut_TDL = tdl(sum_timeDomainSignal);
            % signalOut_TDL = awgn(signalOut_TDL, snr, 'measured');
            % signalOut_TDL = gpuArray(signalOut_TDL);

            %% Demodulation
            isPlot = 0;
            isPrintResult = 0;
            [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr] = PRACH_demodulation_4(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, signalOut_TDL_whithoutNoise, isPlot, isPrintResult);
            [num_row_output, ~] = size(output);

            for row_output = 1:num_row_output
                detectedPreambleIndex = output{row_output, 10};
                if detectedPreambleIndex ~= preIdx_index
                    % threshholdConstant_arr = [threshholdConstant_arr output{row_output, 2}];
                    % ifftCircShift_arr = [ifftCircShift_arr output{row_output, 4}];
                    % ConfigurationIndex_arr = [ConfigurationIndex_arr prachConfigIndex];
                    % RootSequenceIndex_arr = [RootSequenceIndex_arr prachConfig.RootSequenceIndex];
                    % truePreambleIndex_arr = [truePreambleIndex_arr preIdx_index];
                    % detectedPreambleIndex_arr = [detectedPreambleIndex_arr detectedPreambleIndex];
                    % frame_arr = [frame_arr output{row_output, 6}];
                    % subframe_arr = [subframe_arr output{row_output, 8}];
                    % % slot_arr = [slot_arr, ];
                    % % startingSymBol_arr = [];
                    % TA_arr = [TA_arr, output{row_output, 12}];

                    result = [output{row_output, 2}, ... 
                              output{row_output, 4}, ...
                              prachConfigIndex, ...
                              prachConfig.RootSequenceIndex, ...
                              preIdx_index, ...
                              detectedPreambleIndex, ...
                              output{row_output, 6}, ...
                              output{row_output, 8}, ...
                              output{row_output, 10}, ...
                              output{row_output, 12}];
                    dlmwrite(path_filename, result, '-append');
                    % type(path_filename);
                end
            end
            % figure(1)
            % plotResourceGrid(grid);
            % figure(2)
            % plot(abs(timeDomain_signal));
            % % plot(timeDomain_signal);
        end
    end
    fprintf('\b\b\b\b%3.0f%%',prog);
end

testing_result = table(threshholdConstant_arr.', ifftCircShift_arr.', ConfigurationIndex_arr.', RootSequenceIndex_arr.', truePreambleIndex_arr.', detectedPreambleIndex_arr.', frame_arr.', subframe_arr.', TA_arr.');
% filename = strcat('test_longFormat_prachConfigIndex', num2str(prachConfigIndex_arr(1)), num2str(prachConfigIndex_arr(end)), ...
%     '_zeroCorrConfig_', num2str(zeroCorrConfig_arr(1)), num2str(zeroCorrConfig_arr(end)), ...
%     string(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')), ...
%     '.xls');
% path_filename = strcat('test_results\', filename);
% fprintf('\nSaving to  %s', path_filename);
% writetable (testing_result, filename,'Sheet', 1,'Range','D1');
fprintf('  DONE!\n');
% toc