close all
clear
clc

x_u_fft_arr = [];
time_domain_sig_arr = [];
us = [];

pre_idxs = 0:63;

for pre_idx_th = 1:length(pre_idxs)

    prachConfig.PrachConfigurationIndex = 158;
    prachConfig.RootSequenceIndex = 39;
    
    prachConfig.PreambleIndex = pre_idxs(pre_idx_th);
    
    prachConfig.SubcarrierSpacing = 30;  % delta_f_RA
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
    
    [timeDomain_signal, ux] = PRACH_modulation_main(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
    us = [us, ux];
    time_domain_sig_arr = [time_domain_sig_arr, timeDomain_signal];
end

PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

non_zero_indices = find(timeDomain_signal ~= 0);

% Get first and last non-zero index
first_index = non_zero_indices(1);
last_index = non_zero_indices(end);

% snr_arr = -30:5:0;
snr_arr = 15:5:30;
num_samp1snr = 150000;
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
tdl.NumReceiveAntennas = 2;

for snr = 1:length(snr_arr)
    num_sample = 0;
    while num_sample <= num_samp1snr
        for pre_idx = 1:length(pre_idxs)
        
        % 
        signalOut_TDL = tdl(time_domain_sig_arr(:, pre_idx));
        signalOut_TDL = awgn(signalOut_TDL, snr, 'measured');
        % signalOut_TDL = time_domain_sig_arr(:, pre_idx);
        % signalOut_TDL = repmat(signalOut_TDL, 1, tdl.NumReceiveAntennas);
        signalOut_TDL = signalOut_TDL(first_index:last_index, :);
        signalOut_TDL = signalOut_TDL((PrachOFDMInfo.CyclicPrefix + 1):end, :);
        signalOut_TDL = reshape(signalOut_TDL, [nfft, PrachConfigFR1UnpairedSpectrum.PrachDuration, tdl.NumReceiveAntennas]);
        signalOut_TDL = fft(signalOut_TDL, nfft, 1);
        signalOut_TDL = fftshift(signalOut_TDL, 1);
        signalOut_TDL = signalOut_TDL(position + (0:(PrachConfigFR1UnpairedSpectrum.L_RA - 1)), :, :);
        signalOut_TDL = signalOut_TDL * sqrt(PrachConfigFR1UnpairedSpectrum.L_RA);

        x_u = zadoffChu(us(pre_idx), PrachConfigFR1UnpairedSpectrum.L_RA); 
        x_u_fft = fft(x_u);
        x_u_fft = x_u_fft.';
        % x_u_fft_mat = repmat(x_u_fft_arr(:, pre_idx), 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
        
        x_u_fft_mat = repmat(x_u_fft, 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
        x_u_fft_mat = repmat(reshape(x_u_fft_mat, [PrachConfigFR1UnpairedSpectrum.L_RA, PrachConfigFR1UnpairedSpectrum.PrachDuration, 1]), 1, 1, tdl.NumReceiveAntennas);
        % signalOut_TDL = fft(signalOut_TDL);
        signalOut_TDL = conj(signalOut_TDL);
        % x1 = signalOut_TDL(:, 1).';
        % x2 = x_u_fft_mat(:, 1).';
        % corr_fft1 = signalOut_TDL(:, 1).' .* x_u_fft_mat(:, 1).';
        corr_fft = signalOut_TDL .* x_u_fft_mat;

        ifft_corr_fft = zeros(nifft, PrachConfigFR1UnpairedSpectrum.PrachDuration, tdl.NumReceiveAntennas);
        ifft_corr_fft(1:PrachConfigFR1UnpairedSpectrum.L_RA, :, :) = corr_fft;
        ifft_corr_fft = ifft(ifft_corr_fft, nifft, 1);
        ifft_corr_fft = abs(ifft_corr_fft);
        ifft_corr_fft = circshift(ifft_corr_fft, ifftCircShift, 1);
        
        figure(1)
        plot(ifft_corr_fft(:, 1, 1));
        figure(2)
        plot(ifft_corr_fft(:, 1, 2));

        end

    end

end