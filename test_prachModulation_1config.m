prachConfig.PrachConfigurationIndex = 158;
prachConfig.RootSequenceIndex = 39;

prachConfig.PreambleIndex = 60;

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

%% Modulation
% grid = map2grid_longFormat(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
grid = mapPrachSymbol2ResourceGrid_main(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);

timeDomain_signal = PRACH_modulation_main(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
isPlot = 0;
isPrintResult = 1;
[zadoffChuSeq_slot, preambleIndex_timeAdvance_arr, array_ifft_corr_fft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal, isPlot, isPrintResult)
snr_db = -45;
timeDomain_signal_awgn = awgn(timeDomain_signal, snr_db);
figure(1)
plotResourceGrid(grid);
figure(2)
plot(real(timeDomain_signal));
figure(3)
plot(abs(timeDomain_signal_awgn))