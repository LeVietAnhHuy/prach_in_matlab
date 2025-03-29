close all
clear
clc

% %% Prach Config
% % Logical root sequence number (i)
% prachConfig.PrachConfigurationIndex = 158;
% prachConfig.RootSequenceIndex = 39;
% prachConfig.L_RA = 139;
% % prachConfig.PreambleIndex = 0;
% prachConfig.SubcarrierSpacing = 30;  % delta_f_RA
% prachConfig.Set = 'Unrestricted';
% prachConfig.ZeroCorrelationZoneConfig = 8;
% prachConfig.FrequencyRange = 'FR1';
% prachConfig.SpectrumType = 'Unpaired';
% 
% % n_RA_start or nPrachFreqStart
% prachConfig.PrachFreqStart = 0;
% 
% %% Carrier Config
% % N_grid_size_mui, clause 5.3.2, TS 38.211
% carrierConfig.n_UL_RB = 273;
% 
% carrierConfig.SubcarrierSpacing = 30;
% carrierConfig.numElementPerResourceBlock = 12;
% carrierConfig.numFrame = 2;
% %% Code Gen
% PrachConfigFR1UnpairedSpectrum = get_Table6332x(prachConfig, carrierConfig);
% % get N_CS
% N_CS = get_N_CS(prachConfig);
% 
% % Physical root sequence number (u)
% [u, u_arr] = get_u(prachConfig, N_CS);
% 
% % Cyclic shift
% [C_v, C_v_arr] = get_C_v(prachConfig, N_CS);
% 
% x_u = zadoffChu(u, prachConfig.L_RA);
% x_uv = circshift(x_u, [0 -C_v]);
% y_uv = fft(x_uv);
% 
% % Normalization
% y_uv1 = y_uv / sqrt(length(y_uv));
% 
% % Table 6.3.3.2-1, TS 38.211
% % k_bar = 2;
% 
% % Clause 5.3.2, TS 38.211
% k1 = prachConfig.PrachFreqStart * 12 - carrierConfig.n_UL_RB * 6;
% first = k1 + PrachConfigFR1UnpairedSpectrum.k_bar + 1;
% 
% nifft = 4096;
% ifftin = zeros(nifft, 1);
% 
% % Zero padding
% position = nifft / 2 + 0 + first;
% ifftin(nifft / 2 + (0:(prachConfig.L_RA - 1)) + first) = y_uv1;
% 
% shiftfft = fftshift(ifftin);
% ifftout = ifft(fftshift(ifftin), nifft);
% 
% 
% 
% % 12 Prach symbols in time domain.
% SEQ = repmat(ifftout, PrachConfigFR1UnpairedSpectrum.PrachDuration, 1);
% 
% PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
% 
% % Get Cyclic Prefix of Prach in time domain
% CP = ifftout((end - PrachOFDMInfo.CyclicPrefix + 1):end, 1);
% 
% % Get Guard Period of Prach in time domain
% GP = zeros(PrachOFDMInfo.GuardPeriod, 1);
% 
% % Get null remain sample of 1 slot
% twoOFDM = zeros(nifft*2, 1);
% remain = zeros(PrachOFDMInfo.PathProfile, 1);
% 
% preamble = [CP; SEQ; GP; zeros(61440 - (PrachOFDMInfo.Sequence + PrachOFDMInfo.CyclicPrefix + PrachOFDMInfo.GuardPeriod), 1)];

% % Define a carrier configuration object
% carrier = nrCarrierConfig;
% carrier.SubcarrierSpacing = 30;
% carrier.NSizeGrid = 273;
% 
% prach = nrPRACHConfig;
% prach.ConfigurationIndex = 158;
% prach.DuplexMode = 'TDD';
% prach.SubcarrierSpacing = 30;
% prach.ActivePRACHSlot = 1;
% prach.NPRACHSlot = 15;
% prach.ZeroCorrelationZone = 8;
% prach.PreambleIndex = 9;
% 
% [preambleIdx, offset] = nrPRACHDetect(carrier,prach, preamble);

% figure(1)
% plot(real(preamble), 'r');
% title('Preamble in Time domain');
% xlabel('Number of samples'); ylabel('Amplitude');
% xlim([0 61440]);

% grid = mapPrachSymbol2ResourceGrid(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
% cmap = parula(64);
% figure(1);
% axRG = axes;
% image(100*abs(grid));
% axis(axRG,'xy');
% colormap(axRG,cmap);
% title(axRG,sprintf('PRACH Resource Grid (Size in RE [%s])',strjoin(string(size(grid)),' ')));
% xlabel(axRG,'Symbols'); ylabel(axRG,'Subcarriers in RE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cmap = parula(64);
% figure('Name','Resource Grid');
% axRG = axes;
% figure(3);
% grid_in_RB = ResourceGrid_in_ResourceElement2ResourceBlock(grid, 0:1);
% image(abs(grid_in_RB));
% axis(axRG,'xy');
% colormap(axRG,cmap);
% title(axRG,sprintf('PRACH Resource Grid (Size in RB [273 560])'));
% xlabel(axRG,'Symbols'); ylabel(axRG,'Subcarriers in RB');

% figure(4)

% 
% figure(2)
% timeDomain_signal = PRACH_modulation(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
% plot(real(timeDomain_signal));

% output = PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, timeDomain_signal);
figure_num = 0;

%1 % 0:8
            % 1,2 % 5:13
%2 % 9:17                    % 1,2,3 % 5:22
            % 2,3 % 14:22
%3 % 18:26
            % 3,4 % 23:31
%4 % 27:35

%5 % 36:44

%6 % 45:53

%7 % 54:62

%8 % 63

% update peak for ea
UE_preIdxs = [0:9, 18:26];                                        
array_timeDomainSignal = [];

tic
for preIdx_index = 1:length(UE_preIdxs)

    prachConfig.PrachConfigurationIndex = 158;
    prachConfig.RootSequenceIndex = 39;
    prachConfig.L_RA = 139;

    prachConfig.PreambleIndex = UE_preIdxs(preIdx_index);

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
    
    %% put time domain signal for each preambleIndex into a matrix, row by row
    array_timeDomainSignal = [array_timeDomainSignal, timeDomain_signal];
end

%% sum all time domain signal
sum_timeDomainSignal = sum(array_timeDomainSignal, 2);

%% tdl config, add noise

signalOut_TDL_whithoutNoise = sum_timeDomainSignal;

tdl = nrTDLChannel;

tdl.SampleRate = prachConfig.SubcarrierSpacing * 1000 * 4096;
tdl.TransmissionDirection = 'Uplink';
tdl.MaximumDopplerShift = 100;  % 100 Hz
tdl.DelayProfile = 'TDL-C';
tdl.DelaySpread = 300e-9;
tdl.NumTransmitAntennas = 1;
tdl.NumReceiveAntennas = 4;

signalOut_TDL = tdl(sum_timeDomainSignal);
signalOut_TDL = awgn(signalOut_TDL, 10, 'measured');
signalOut_TDL = gpuArray(signalOut_TDL);

% figure_num = figure_num + 1;
% figure(figure_num)
% plot(real(signalOut_TDL));

% disp(preIdx);

%% Demodulation
[zadoffChuSeq_slot, output, corr_iffft, ifftSeq_array, mean_u_ifft, pos, C_v_arr] = PRACH_demodulation_4(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, signalOut_TDL);

% %% get window sample
% posDiff = [1];
% for posIdx = 2:2:length(pos)
%     if posIdx >= length(pos)
%         break;
%     end
%     posDiff = [posDiff (pos(posIdx) + pos(posIdx + 1)) / 2];
% end
% posDiff = [posDiff pos(length(pos))];
% 
% %% put C_v index into plot
% C_V_index_text_pos = [];
% for textIdx = 1:(length(posDiff) - 1)
%     C_V_index_text_pos = [C_V_index_text_pos (posDiff(textIdx) + posDiff(textIdx + 1)) / 2];
% end
% xticks(posDiff);

% [numSignal, numElement] = size(corr_iffft);
% 
% 
% for signalIndex = 1:numSignal  
%     figure_num = figure_num + 1;
%     figure(figure_num)
%     findpeaks(corr_iffft(signalIndex, :), MinPeakHeight=(max(corr_iffft(signalIndex, :)) / 1.1));
% end
% 
% xticks(pos);
% axis tight
% 
% [numElement, numSignal] = size(ifftSeq_array);


% for signalIndex = 1:numSignal 
%     figure_num = figure_num + 1;
%     figure(figure_num)
%     findpeaks(ifftSeq_array(:, signalIndex), MinPeakHeight=(max(ifftSeq_array(:, signalIndex)) / 1.1));
% end
% xticks(pos);
% axis tight

toc
disp('----')

% tic
% for preIdx = 0:63
% 
%     prachConfig.PrachConfigurationIndex = 168;
%     prachConfig.RootSequenceIndex = 39;
%     prachConfig.L_RA = 139;
% 
%     prachConfig.PreambleIndex = preIdx;
% 
%     % prachConfig.PreambleIndex = 0;
%     prachConfig.SubcarrierSpacing = 30;  % delta_f_RA
%     prachConfig.Set = 'Unrestricted';
%     prachConfig.ZeroCorrelationZoneConfig = 2;
%     prachConfig.FrequencyRange = 'FR1';
%     prachConfig.SpectrumType = 'Unpaired';
% 
%     % n_RA_start or nPrachFreqStart
%     prachConfig.PrachFreqStart = 0;
% 
%     %% Carrier Config
%     % N_grid_size_mui, clause 5.3.2, TS 38.211
%     carrierConfig.n_UL_RB = 273;
% 
%     carrierConfig.SubcarrierSpacing = 30;
%     carrierConfig.numElementPerResourceBlock = 12;
%     carrierConfig.numFrame = 2;
% 
%     % get N_CS
%     N_CS = get_N_CS(prachConfig);
% 
%     % Physical root sequence number (u)
%     [u, u_arr] = get_u(prachConfig, N_CS);
% 
%     % Cyclic shift
%     [C_v, C_v_arr] = get_C_v(prachConfig, N_CS);
% 
%     x_u = zadoffChu(u, prachConfig.L_RA);
%     x_uv = circshift(x_u, [0 -C_v]);
%     y_uv = fft(x_uv);
%     PrachConfigFR1UnpairedSpectrum = get_Table6332x(prachConfig, carrierConfig);
%     grid = mapPrachSymbol2ResourceGrid(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
%     timeDomain_signal = PRACH_modulation(grid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
%     % timeDomain_signal = gpuArray(timeDomain_signal);
% 
%     % % save("timeDomain_signal.mat", "timeDomain_signal");
%     % x = load("timeDomain_signal.mat");
%     % timeDomain_signal = timeDomain_signal + x.timeDomain_signal;
% 
%     % figure_num = figure_num + 1;
%     % figure(figure_num)
%     % plot(real(timeDomain_signal));
% 
%     % disp(preIdx);
%     % output = PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, timeDomain_signal);
%     % disp('----')
% 
%     % cmap = parula(64);
%     % figure('Name','Resource Grid');
%     % axRG = axes;
% 
%     % figure(3);
%     % grid_in_RB = ResourceGrid_in_ResourceElement2ResourceBlock(grid, 0:1);
%     % image(abs(grid_in_RB));
%     % axis(axRG,'xy');
%     % colormap(axRG,cmap);
%     % title(axRG,sprintf('PRACH Resource Grid (Size in RB [273 560])'));
%     % xlabel(axRG,'Symbols'); ylabel(axRG,'Subcarriers in RB');
% 
%     grid = mapPrachSymbol2ResourceGrid(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv);
%     % cmap = parula(64);
%     % figure_num = figure_num + 1;
%     % figure(figure_num);
%     % axRG = axes;
%     % image(100*abs(grid));
%     % axis(axRG,'xy');
%     % colormap(axRG,cmap);
%     % title(axRG,sprintf('PRACH Resource Grid (Size in RE [%s])',strjoin(string(size(grid)),' ')));
%     % xlabel(axRG,'Symbols'); ylabel(axRG,'Subcarriers in RE');
% 
%     signalOut_TDL = timeDomain_signal;
% 
%     tdl = nrTDLChannel;
% 
%     %% tdl config
%     tdl.SampleRate = prachConfig.SubcarrierSpacing * 1000 * 4096;
%     tdl.TransmissionDirection = 'Uplink';
%     tdl.MaximumDopplerShift = 100;  % 100 Hz
%     tdl.DelayProfile = 'TDL-C';
%     tdl.DelaySpread = 300e-9;
%     tdl.NumTransmitAntennas = 1;
%     tdl.NumReceiveAntennas = 4;
% 
%     signalOut_TDL = tdl(timeDomain_signal);
%     signalOut_TDL = awgn(signalOut_TDL, 10, 'measured');
%     signalOut_TDL = gpuArray(signalOut_TDL);
% 
%     % figure_num = figure_num + 1;
%     % figure(figure_num)
%     % plot(real(signalOut_TDL));
% 
%     disp(preIdx);
%     [zadoffChuSeq_slot, output, corr_iffft, pos, C_v_arr] = PRACH_demodulation_4(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, signalOut_TDL);
% 
%     % %% get window sample
%     % posDiff = [1];
%     % for posIdx = 2:2:length(pos)
%     %     if posIdx >= length(pos)
%     %         break;
%     %     end
%     %     posDiff = [posDiff (pos(posIdx) + pos(posIdx + 1)) / 2];
%     % end
%     % posDiff = [posDiff pos(length(pos))];
%     % 
%     % %% put C_v index into plot
%     % C_V_index_text_pos = [];
%     % for textIdx = 1:(length(posDiff) - 1)
%     %     C_V_index_text_pos = [C_V_index_text_pos (posDiff(textIdx) + posDiff(textIdx + 1)) / 2];
%     % end
%     % plot(corr_iffft);
%     % xticks(posDiff);
%     toc
%     disp('----')
% end



