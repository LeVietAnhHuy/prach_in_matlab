close all
clear
clc

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
preIdx = 2;
array_timeDomainSignal = [];
constant = 1.5;
start_snr = -45;
end_snr = -39;

% snr_range = flip(start_snr:end_snr);
snr_range = start_snr:end_snr;
gain = 1;

% Table 8.4.1.1-1 in TS 38.141-1v15.03, page 166
% timeErrorTolerance = 0.26; % micro second, for AWGNs
% timeErrorTolerance = 1.77; % micro second, for TDLC300-100

% prachConfig.PrachConfigurationIndex = randi([67 250], 1, 1);
prachConfig.PrachConfigurationIndex = 158;

% prachConfig.RootSequenceIndex = randi([0 137], 1, 1);
prachConfig.RootSequenceIndex = 39;

prachConfig.PreambleIndex = preIdx;

prachConfig.SubcarrierSpacing = 30;  % delta_f_RA
prachConfig.Set = 'Unrestricted';
% prachConfig.ZeroCorrelationZoneConfig = randi([0 15], 1, 1);
prachConfig.ZeroCorrelationZoneConfig = 8;
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
carrierConfig.numFrame = 1;

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

timeDomain_signal = PRACH_modulation_main_plot1(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
timeDomain_signal = signalScale(timeDomain_signal, gain);

%% Step 3: Multipath Fading
TDLC300_100 = nrTDLChannel("DelayProfile", "TDLC300", ...
    "MaximumDopplerShift", 100, ...
    "SampleRate", prachConfig.SubcarrierSpacing * 1000 * 4096, ...
    "TransmissionDirection", "Uplink", ...
    "NumTransmitAntennas", 1, ...
    "NumReceiveAntennas", 1);
%% Testing procedure in Clause 8.4.1.4  in TS 38.141-1v15.03
%% Part 1: Delay profiles
%% Step 1: AWGN
% Table 8.4.1.4.2-1 in TS 38.141-1v15.03
% snr_dbm = 100;
% snr_dbW = snr_dbm - 30;
% timeErrorTolerance = 0.26; % micro second, for AWGNs
timeErrorTolerance = 1.77; % micro second, for TDLC300-100
% timeErrorTolerance = 77; % micro second, for TDLC300-100

% Get current date and time
currentTime = datetime('now');

% Convert to desired format
formattedText = datestr(currentTime, 'ddmmyy_HHMMSS');
% Convert components to text
dateText = datestr(currentTime, 'dd-mm-yyyy HH:MM:SS');

% Define the file name
fileDir = 'test_results/';
fileName = strcat(fileDir, 'log_testTDLC_', formattedText, '.txt');

% Open the file for writing
fileID = fopen(fileName, 'w');

% Check if the file opened successfully
if fileID == -1
    error('Cannot open file for writing.');
end


notDetected_test_id_arr = [];
num_test = 1000;

fprintf(fileID, 'Test time: %s\n', dateText);
fprintf(fileID, 'Number of test: %d\n', num_test);
fprintf(fileID, 'Signal gain: %d\n', gain);
fprintf(fileID, 'Threshold constant: %d\n', constant);

fprintf(fileID, '\nNoise:\n');
fprintf(fileID, '\tMulti-path fading: TDLC300-100\n');
fprintf(fileID, '\tAWGN: %ddB -> %ddB\n', start_snr, end_snr);

fprintf(fileID, '\nConfiguration:\n');
fprintf(fileID, '\tPreambleIndex: %d\n', prachConfig.PreambleIndex);
fprintf(fileID, '\tPrachConfigurationIndex: %d\n', prachConfig.PrachConfigurationIndex);
fprintf(fileID, '\tRootSequenceIndex: %d\n', prachConfig.RootSequenceIndex);
fprintf(fileID, '\tZeroCorrelationZoneConfig: %d\n', prachConfig.ZeroCorrelationZoneConfig);
fprintf(fileID, '\tFrequencyRange: %s\n', prachConfig.FrequencyRange);
fprintf(fileID, '\tSpectrumType: %s\n', prachConfig.SpectrumType);
fprintf(fileID, '\tPrachFreqStart: %d\n', prachConfig.PrachFreqStart);

numDetectedPreambleIndex = 0;
% timeErrorTolerance = 0.26; % micro second, for AWGNs
timeErrorTolerance = 1.77; % micro second, for TDLC300-100

fprintf(fileID, '\tCase 3: AWGN + TDLC300-100\n');
fprintf('Case 6: AWGN + TDLC300-100\n');
for snr = snr_range
    numDetectedPreambleIndex = 0;
    for test_id = 1:num_test
        fprintf('TEST %d:\n', test_id);
        fprintf(fileID, 'TEST %d:\n', test_id);

        figure_num = 1;

        snr_dbW = snr;

        %% Received Signal
        timeDomain_signal_AWGN = awgn(timeDomain_signal, snr_dbW, 'measured');
        timeDomain_signal_AWGN_TDL = TDLC300_100(timeDomain_signal_AWGN);

        %% Demodulation
        isPlot = 0;
        isPrintResult = 0;
        %% with AWGN + Multipath Fading
        [zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyRandomCombining_PRACH_demodulation_plot1(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal_AWGN_TDL, isPlot, isPrintResult, u, C_v, figure_num, fileID);
        if ~isempty(detected_preambleIndex_arr)
            detectedPreambleIndex = mode(detected_preambleIndex_arr);
            detected_timeErrorTolerance = mode(nTA_microSecond_arr);

            fprintf(fileID, '\tdetectedPreambleIndex = %d; detected_timeErrorTolerance = %.3f\n', detectedPreambleIndex, detected_timeErrorTolerance);

            if detectedPreambleIndex == preIdx
                numDetectedPreambleIndex = numDetectedPreambleIndex + 1;
            end
        else
            notDetected_test_id_arr = [notDetected_test_id_arr, test_id];
        end

        figure_num = figure_num + 1;

    end

    acc = (numDetectedPreambleIndex / num_test) * 100;

    fprintf('%ddB: %.3f%%\n', snr, acc);
    fprintf(fileID, '\t%ddB: %.3f%%\n', snr, acc);

    fprintf(fileID, 'Not detected Test: \n');
    for notDetectedText_idx = notDetected_test_id_arr
        fprintf(fileID, '%d\n', notDetectedText_idx);
    end
    disp(notDetected_test_id_arr);

    % if acc == 0
    %     break;
    % end
end
fclose('all');

