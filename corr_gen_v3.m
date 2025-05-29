close all
clear
clc

x_u_fft_arr = [];
time_domain_sig_arr = [];
us = [];

pre_idxs = 0:63;

for pre_idx_th = 1:length(pre_idxs)

    prachConfig.PrachConfigurationIndex = 168;
    prachConfig.RootSequenceIndex = 39;

    prachConfig.PreambleIndex = pre_idxs(pre_idx_th);

    prachConfig.SubcarrierSpacing = 30;  % delta_f_RA
    prachConfig.Set = 'Unrestricted';
    % prachConfig.ZeroCorrelationZoneConfig = randi([0 15], 1, 1);
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

    constant = 27;

    % Table 8.4.1.1-1 in TS 38.141-1v15.03, page 166
    timeErrorTolerance = 0.26; % micro second, for AWGNs

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

    x_u_fft = fft(x_u);
    x_u_fft_arr = [x_u_fft_arr, x_u_fft'];

    %% Modulation
    % grid = map2grid_longFormat(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
    grid = mapPrachSymbol2ResourceGrid_main(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);

    [timeDomain_signal, ux, start_idx_preamble_arr, end_idx_preamble_arr] = PRACH_modulation_main(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
    % plot(abs(timeDomain_signal));
    us = [us, ux];
    time_domain_sig_arr = [time_domain_sig_arr, timeDomain_signal];
end

PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

non_zero_indices = find(timeDomain_signal ~= 0);

% Get first and last non-zero index
first_index = non_zero_indices(1);
last_index = non_zero_indices(end);

snr_arr = -50:5:0;

% Quỳnh Anh 
% snr_arr = -50:5:-35;

% Hoàng
% snr_arr = -30:5:0;

folderName = 'mat_pair_fft_data_MIMO';

if ~exist(folderName, 'dir')
    mkdir(folderName);
end

num_slot_contain_prach = length(PrachConfigFR1UnpairedSpectrum.subframeNumber);
num_samp1snr = 120000;
position = 2;
nfft = 4096;
nifft = 1024;
ifftCircShift = -10;

tdl = nrTDLChannel;
tdl.SampleRate = prachConfig.SubcarrierSpacing * 1000 * 4096;
tdl.TransmissionDirection = 'Uplink';
tdl.MaximumDopplerShift = 100;  % 100 Hz
tdl.DelayProfile = 'TDL-C';
tdl.DelaySpread = 300e-9;
tdl.NumTransmitAntennas = 1;
tdl.NumReceiveAntennas = 8;

for snr = snr_arr
    
    matfile_names = {};
    mat_data_objects = {};
    
    for matfile_idx = 1:tdl.NumReceiveAntennas
        matfile_names{end + 1} = strcat(folderName, '/pair_fft_data_', num2str(snr), 'dB_rx', num2str(matfile_idx), '.mat');
        mat_data_objects{end + 1} = matfile(matfile_names{matfile_idx}, 'Writable', true);  % Append matfile object
    end
    
    next_rows = ones(tdl.NumReceiveAntennas, 1);
    
    % nextRow = 1;

    while next_rows(1) <= num_samp1snr
        pre_idx_add_noises = ones(tdl.NumReceiveAntennas, 1);

        for pre_idx = 1:length(pre_idxs)
            signalOut_TDL = tdl(time_domain_sig_arr(:, pre_idx));   % 8x1228800
            signalOut_TDL = awgn(signalOut_TDL, snr, 'measured');

            for rx_idx = 1:tdl.NumReceiveAntennas
                signalOut_TDL_1rx = signalOut_TDL(:, rx_idx);  % 1x1228800
                signalOut_TDL_1rx = reshape(signalOut_TDL_1rx, length(signalOut_TDL_1rx) / num_slot_contain_prach, num_slot_contain_prach); % 122880x10
                signalOut_TDL_1rx = signalOut_TDL_1rx(start_idx_preamble_arr(1):end_idx_preamble_arr(1), :); % 51088x10
 
                % signalOut_TDL = time_domain_sig_arr(:, pre_idx);
                % signalOut_TDL = repmat(signalOut_TDL, 1, tdl.NumReceiveAntennas);
                % signalOut_TDL = signalOut_TDL(first_index:last_index, :);
                signalOut_TDL_1rx = signalOut_TDL_1rx((PrachOFDMInfo.CyclicPrefix + 1):end, :);
                signalOut_TDL_1rx = reshape(signalOut_TDL_1rx, [nfft, PrachConfigFR1UnpairedSpectrum.PrachDuration, num_slot_contain_prach]);
                signalOut_TDL_1rx = fft(signalOut_TDL_1rx, nfft, 1);
                signalOut_TDL_1rx = fftshift(signalOut_TDL_1rx, 1);
                signalOut_TDL_1rx = signalOut_TDL_1rx(position + (0:(PrachConfigFR1UnpairedSpectrum.L_RA - 1)), :, :);
                signalOut_TDL_1rx = signalOut_TDL_1rx * sqrt(PrachConfigFR1UnpairedSpectrum.L_RA);
                
                signalOut_TDL_1rx = permute(signalOut_TDL_1rx, [2 1 3]);
 
                %% No peak
                if(pre_idx_add_noises(rx_idx) > length(C_v_arr))
                    
                    temp_signalOut_TDL_1rx = signalOut_TDL_1rx;

                    x_u = zadoffChu(us(mod(pre_idx + 20, length(pre_idxs))), PrachConfigFR1UnpairedSpectrum.L_RA);

                    reduced_pre_idx_noise = length(C_v_arr);

                    reduced_pre_idx_mat = zeros(size(signalOut_TDL_1rx, 1), 1, size(signalOut_TDL_1rx, 3)) + reduced_pre_idx_noise;

                    x_u_fft = fft(x_u);
                    x_u_fft = x_u_fft.';

                    x_u_fft_mat = repmat(x_u_fft, 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
                    x_u_fft_mat = repmat(reshape(x_u_fft_mat, [PrachConfigFR1UnpairedSpectrum.L_RA, PrachConfigFR1UnpairedSpectrum.PrachDuration, 1]), 1, 1, num_slot_contain_prach);
                    
                    x_u_fft_mat = permute(x_u_fft_mat, [2 1 3]);

                    signalOut_TDL_1rx = cat(2, signalOut_TDL_1rx, x_u_fft_mat, reduced_pre_idx_mat);

                    [r, c, ~] = size(signalOut_TDL_1rx);

                    for k = 1:size(signalOut_TDL_1rx, 3)
                        slice = signalOut_TDL_1rx(:, :, k);
                        rows  = next_rows(rx_idx) : next_rows(rx_idx) + r - 1;

                        mat_data_objects{rx_idx}.data(rows, 1:c) = slice;
                        next_rows(rx_idx) = next_rows(rx_idx) + r;
                    end

                    fprintf('%d/%d, reduced_pre_idx = noise, true_pre_idx = noise, snr = %d(dB), rx_idx = %d\n', next_rows(rx_idx), num_samp1snr, snr, rx_idx);
                end

                %% Contain peak
                x_u = zadoffChu(us(pre_idx), PrachConfigFR1UnpairedSpectrum.L_RA);
                reduced_pre_idx = mod(pre_idx - 1, length(C_v_arr));
                reduced_pre_idx_mat = zeros(size(signalOut_TDL_1rx, 1), 1, size(signalOut_TDL_1rx, 3)) + reduced_pre_idx;

                x_u_fft = fft(x_u);
                x_u_fft = x_u_fft.';
                % x_u_fft_mat = repmat(x_u_fft_arr(:, pre_idx), 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);

                x_u_fft_mat = repmat(x_u_fft, 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
                x_u_fft_mat = repmat(reshape(x_u_fft_mat, [PrachConfigFR1UnpairedSpectrum.L_RA, PrachConfigFR1UnpairedSpectrum.PrachDuration, 1]), 1, 1, num_slot_contain_prach);

                x_u_fft_mat = permute(x_u_fft_mat, [2 1 3]);
                
                if(pre_idx_add_noises(rx_idx) > length(C_v_arr))
                    signalOut_TDL_1rx = cat(2, temp_signalOut_TDL_1rx, x_u_fft_mat, reduced_pre_idx_mat);
                    pre_idx_add_noises(rx_idx) = 1;

                    clear temp_signalOut_TDL_1rx;
                else
                    signalOut_TDL_1rx = cat(2, signalOut_TDL_1rx, x_u_fft_mat, reduced_pre_idx_mat);
                end

                [r, c, ~] = size(signalOut_TDL_1rx);

                for k = 1:size(signalOut_TDL_1rx, 3)
                    slice = signalOut_TDL_1rx(:, :, k);
                    rows  = next_rows(rx_idx) : next_rows(rx_idx) + r - 1;

                    mat_data_objects{rx_idx}.data(rows, 1:c) = slice;
                    next_rows(rx_idx) = next_rows(rx_idx) + r;
                end

                fprintf('%d/%d, reduced_pre_idx = %d, true_pre_idx = %d, snr = %d(dB), rx_idx = %d\n', next_rows(rx_idx), num_samp1snr, reduced_pre_idx, pre_idxs(pre_idx), snr, rx_idx);

                pre_idx_add_noises(rx_idx) = pre_idx_add_noises(rx_idx) + 1;
                % break;

                % figure(1)
                % plot(ifft_corr_fft(1, :, 1));
            end
            % break;
        end
    end
    % break;
end