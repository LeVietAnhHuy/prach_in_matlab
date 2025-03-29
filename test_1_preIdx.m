preIdx = 1;

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
[zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, C_v_arr, preambleIndex_arr] = PRACH_demodulation_8(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, signalOut_TDL);
disp(preambleIndex_arr);