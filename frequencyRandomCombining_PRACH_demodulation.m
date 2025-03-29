function [zadoffChuSeq_slot, preambleIndex_timeAdvance_arr, array_ifft_corr_fft, ifftSeq_array, mean_u_ifft, C_v_arr, detected_preambleIndex_arr, nTA_microSecond_arr] = frequencyRandomCombining_PRACH_demodulation(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, constant, timeErrorTolerance, timeDomain_signal, isPlot, isPrintResult)
% 27 December 2024

default_numSample_perSlot = 30720;
BW = 98.28e3;
numResourceElement_freqDomain = BW / prachConfig.SubcarrierSpacing;

sampleTime_microSecond = (10 / (length(timeDomain_signal) / carrierConfig.numFrame)) * 10^3;

%N_RA_RB Number of resource blocks occupied which is given by the
% parameter allocation expressed in number of RBs for PUSCH.
% TS 38.211 Table 6.3.3.2-1.
N_RA_RB = PrachConfigFR1UnpairedSpectrum.N_RA_RB;

%% Get frequency locations of PRACH on resource grid
% Based on Clause 5.4.2 in TS 38.211 version 18.2.0

% Find k0_mui

N_start_mui_grid = 0; % Common resource block indicated by higher-signalling
N_size_mui_grid = carrierConfig.n_UL_RB;
N_RB_sc = 12;

% mui_0 is the largest  value among
% the subcarrier spacing configurations
% by the higher-layer parameter scs-SpecificCarrierList;

% At this stage, we assume
mui = log2(carrierConfig.SubcarrierSpacing / 15);
mui_0 = mui;
N_start_mui_0_grid = N_start_mui_grid;
N_size_mui_0_grid = N_size_mui_grid;

k_0_mui = (N_start_mui_grid + N_size_mui_grid / 2)*N_RB_sc - ...
    (N_start_mui_0_grid + N_size_mui_0_grid / 2)*N_RB_sc*2^(mui_0 - mui);

% Find value of final operand in curly bracket in k_1 formula
if any(PrachConfigFR1UnpairedSpectrum.L_RA == [571, 1151]) && prachConfig.FrequencyRange == "FR1"
    %RBSetOffset Starting RB index of the uplink RB set for this PRACH transmission occasion (0...274)
    %   Specify the RB index of the uplink RB set for this PRACH
    %   transmission occasion. It is defined as the difference between
    %   the start CRB of uplink RB sets related to "n_RA^start + n_RA"
    %   and "n_RA^start", respectively. This property is used in the
    %   computation of the PRACH indices, as discussed in TS 38.211
    %   Section 5.3.2, in the case of FR1 and LRA 571 or 1151. It must
    %   be an integer scalar in the range 0...274.
    %   The default is 0.
    % (specify in prachConfig)
    RBSetOffset = 0;
    rbsetOffset = RBSetOffset;
else % (LRA == 139,839) || (LRA == 571, 1151 && FR2)
    %FrequencyIndex Index of the PRACH transmission occasion in frequency domain (0...7)
    %   Specify the frequency index of the PRACH transmission occasion.
    %   It is defined as an integer in the range {0:M-1}, in which M is
    %   the higher layer parameter "msg1-FDM" defined in TS 38.331
    %   Section 6.3.2. M can assume values in {1, 2, 4, 8}.
    %   FrequencyIndex is referred to as "n_RA" in TS 38.211 Sections
    %   5.3.2 and 6.3.3.2. This property is inactive for an FR1 carrier
    %   with LRA = {571, 1151}.
    %   The default is 0.
    % (specify in prachConfig)
    FrequencyIndex = 0;
    rbsetOffset = FrequencyIndex*N_RB_sc;
end

% find k_1

% find RBOffset - first element in second operand in k_1 formula
%RBOffset Starting RB index of the initial uplink BWP (0...274)
%   Specify the RB index where the initial uplink BWP starts
%   relative to the carrier resource grid. It must be an integer
%   scalar in the range 0...274.
%   The default is 0.
% (specify in prachConfig)
RBOffset = 0;

% find FrequencyStart - "n_RA^start" - first element in fourth operand in k_1 formula
%FrequencyStart Frequency offset of lowest PRACH transmission occasion in frequency domain with respect to PRB 0 (0...274)
%   Specify the frequency offset, as defined by the higher layer
%   parameter "msg1-FrequencyStart" and referred to as "n_RA^start"
%   in TS 38.211 Section 5.3.2.
%   The default is 0.
FrequencyStart = 0;

k_1 = k_0_mui + RBOffset*N_RB_sc - N_size_mui_grid*N_RB_sc / 2 + FrequencyStart*N_RB_sc + rbsetOffset;

% Find K
K = carrierConfig.SubcarrierSpacing / prachConfig.SubcarrierSpacing;

% PRACH-specific frequency shift
first = K*k_1 + PrachConfigFR1UnpairedSpectrum.k_bar;
% Set the origin of the frequency domain to the middle of the grid
zeroFreq = double(numResourceElement_freqDomain/2);
PrachStartingResourceElementIndex_freqDomain = ceil(first + zeroFreq);

PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

frame_mod_x = mod(0:(carrierConfig.numFrame - 1), PrachConfigFR1UnpairedSpectrum.x);
frame_contain_prach = find(frame_mod_x == PrachConfigFR1UnpairedSpectrum.y) - 1;

%% Compute Average Antenna
num_RX_antenna = size(timeDomain_signal, 2);
avg_timeDomainSignal = sum(timeDomain_signal, 2) / num_RX_antenna;
timeDomain_signal = avg_timeDomainSignal;

%% Initialize necessary arrays
preambleIndex_timeAdvance_arr = [];
zadoffChuSeq_slot = [];
array_ifft_corr_fft = [];
ifftSeq_array = [];
mean_u_ifft = [];
preamble_only_3dim_arr = [];
preamble_info_arr = [];

allframe = 0:(carrierConfig.numFrame - 1);

%% Find nfft
switch PrachConfigFR1UnpairedSpectrum.preambleFormat
    case {'0'}
        nfft = PrachOFDMInfo.Sequence;
    case {'1', '2', '3'}
        nfft = PrachOFDMInfo.Sequence / PrachConfigFR1UnpairedSpectrum.PrachDuration;
    otherwise
        grid_column = carrierConfig.n_UL_RB * carrierConfig.numElementPerResourceBlock;
        min_nfft = 128;

        while min_nfft < grid_column
            min_nfft = min_nfft * 2;
        end

        nfft = min_nfft;
end

ifftCircShift = -10;

%% Extract preamble from received signal
switch PrachConfigFR1UnpairedSpectrum.preambleFormat
    case {'0', '1', '2', '3'}
        % totalSample = length(timeDomain_signal);
        % [~, numSubframe] = size(resourceGrid);
        PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
        % Number of sample of 1 slot for prachConfig.SubcarrierSpacing
        % numSample_perSlot = default_numSample_perSlot * (prachConfig.SubcarrierSpacing / 15);
        numSample_persubframe = default_numSample_perSlot * (carrierConfig.SubcarrierSpacing / 15) * 2;
        numSample_perFrame = numSample_persubframe * 10;
        % numSlot_perFrame = 10 * (prachConfig.SubcarrierSpacing / 15);
        % totalSlot = carrierConfig.numFrame * numSlot_perFrame;

        frame_mod_x = mod(0:(carrierConfig.numFrame - 1), PrachConfigFR1UnpairedSpectrum.x);
        frame_contain_prach = find(frame_mod_x == PrachConfigFR1UnpairedSpectrum.y) - 1;

        for frame = allframe
            if ismember(frame, frame_contain_prach)
                zadoffChuSeq_symbol = [];
                for mappingSubframe = PrachConfigFR1UnpairedSpectrum.subframeNumber
                    startingSampleSymbol = frame * numSample_perFrame + numSample_persubframe * mappingSubframe + PrachConfigFR1UnpairedSpectrum.startingSymbol * 2048 * (carrierConfig.SubcarrierSpacing / 15) * 2 + 1;
                    preamble_only = timeDomain_signal((startingSampleSymbol + PrachOFDMInfo.CyclicPrefix):(startingSampleSymbol + PrachOFDMInfo.CyclicPrefix + PrachConfigFR1UnpairedSpectrum.PrachDuration * nfft - 1));
                    preamble_only = reshape(preamble_only, [nfft, PrachConfigFR1UnpairedSpectrum.PrachDuration]);

                    preamble_only_3dim_arr = cat(3, preamble_only_3dim_arr, preamble_only);

                    preamble_info = [frame mappingSubframe];
                    preamble_info_arr = [preamble_info_arr; preamble_info];
                end
            end
        end

    otherwise
        % Check for n_RA_slot, TS 38.211, section 5.3.2
        if (prachConfig.SubcarrierSpacing == 30 || prachConfig.SubcarrierSpacing == 120) && PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe == 1
            n_RA_slot = 1;
        else
            n_RA_slot = [0, 1];
        end

        % Number of sample of 1 slot for prachConfig.SubcarrierSpacing
        numSample_perSlot = default_numSample_perSlot * (prachConfig.SubcarrierSpacing / 15);
        numSlot_perFrame = 10 * (prachConfig.SubcarrierSpacing / 15);
        totalSlot = carrierConfig.numFrame * numSlot_perFrame;

        slot_contain_prach = [];
        for frameIndex = 1:length(frame_contain_prach)
            for subFrameNumber_index = 1:length(PrachConfigFR1UnpairedSpectrum.subframeNumber)
                for n_RA_slot_index = 1:length(n_RA_slot)
                    subframeNumber_arr = PrachConfigFR1UnpairedSpectrum.subframeNumber;
                    slot_contain_prach = [slot_contain_prach ((frame_contain_prach(frameIndex) * 10 + subframeNumber_arr(subFrameNumber_index)) * 2 + n_RA_slot(n_RA_slot_index))];

                end
            end
        end

        start_symbol_contain_prach = PrachConfigFR1UnpairedSpectrum.startingSymbol + (0:PrachConfigFR1UnpairedSpectrum.numTimeDomainPrachOccasionsWithinAPrachSlot - 1) * PrachConfigFR1UnpairedSpectrum.PrachDuration;

        for slot_index = 1:length(slot_contain_prach)
            startingSample = slot_contain_prach(slot_index) * numSample_perSlot + 1;

            zadoffChuSeq_symbol = [];
            for symbol_index = 1:length(start_symbol_contain_prach)
                startingSampleSymbol = startingSample + start_symbol_contain_prach(symbol_index) * nfft;

                preamble_only = timeDomain_signal((startingSampleSymbol + PrachOFDMInfo.CyclicPrefix):(startingSampleSymbol + PrachOFDMInfo.CyclicPrefix + PrachConfigFR1UnpairedSpectrum.PrachDuration * nfft - 1));
                preamble_only = reshape(preamble_only, [nfft, PrachConfigFR1UnpairedSpectrum.PrachDuration]);
                preamble_only_3dim_arr = cat(3, preamble_only_3dim_arr, preamble_only);

                preamble_info = [slot_contain_prach(slot_index) start_symbol_contain_prach(symbol_index)];
                preamble_info_arr = [preamble_info_arr; preamble_info];
            end
        end
end

[~, ~, num_preamble_only] = size(preamble_only_3dim_arr);

%% Find Multiples of PrachDuration
prachDuration = PrachConfigFR1UnpairedSpectrum.PrachDuration;
multiples = [];
for multiple = 1:prachDuration
    if mod(prachDuration, multiple) == 0
        multiples = [multiples, multiple];
    end
end

%% Method 1: Average preambles in time domain
%% Modulation

detected_preambleIndex_arr = [];
nTA_microSecond_arr = [];

for multipleIndex = 1:length(multiples)
    for preamble_only_idx = 1:num_preamble_only
        preamble_only = preamble_only_3dim_arr(:, :, preamble_only_idx);

        multiple = multiples(multipleIndex);

        blockSize = prachDuration / multiple;
        
        numPrachDuration = 1:prachDuration;
        suffle_numPrachDuration = numPrachDuration(randperm(length(numPrachDuration)));

        prachBlockIndex_arr = reshape(suffle_numPrachDuration, [blockSize, multiple]);
        prachBlockSumFreq_arr = [];
        for prachBlockIndex = 1:multiple
            prachBlockSumFreq = [];
            for prachIndex = 1:blockSize
                prachBlock_fft = fft(preamble_only(:, prachBlockIndex_arr(prachIndex, prachBlockIndex)), nfft);
                prachBlockSumFreq = [prachBlockSumFreq, prachBlock_fft];
            end
            % prachBlockSumFreq = sum(prachBlockSumFreq, 2);
            prachBlockSumFreq = sum(prachBlockSumFreq, 2) / blockSize;
            prachBlockSumFreq_arr = [prachBlockSumFreq_arr, prachBlockSumFreq];
        end

        for prachSumFreqIndex = 1:multiple
            prachSumFreq_fftshift = fftshift(prachBlockSumFreq_arr(:, prachSumFreqIndex));
            position = PrachStartingResourceElementIndex_freqDomain;
            y_uv = prachSumFreq_fftshift(position:(position + PrachConfigFR1UnpairedSpectrum.L_RA - 1));
            y_uv = y_uv * sqrt(length(y_uv));
            y_uv = y_uv.';
            x_uv = ifft(y_uv);

            x_uv_avg = x_uv;

            % add noise later!

            N_CS = get_N_CS(prachConfig, PrachConfigFR1UnpairedSpectrum);
            [~, u_arr] = get_u(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);
            [~, C_v_arr] = get_C_v(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

            numPreablePerSeq = zeros(length(u_arr));
            for u_index = 1:length(u_arr)
                if u_index < length(u_arr)
                    numPreablePerSeq(u_index) = length(C_v_arr);
                else
                    numPreablePerSeq(u_index) = 64 - (u_index-1)*length(C_v_arr);
                end
            end

            for u_index = 1:length(u_arr)

                switch PrachConfigFR1UnpairedSpectrum.preambleFormat
                    case {'0', '1', '2', '3'}
                        ifftLength = 2048;
                    otherwise
                        ifftLength = 4096;
                end
                x_u = zadoffChu(u_arr(u_index), PrachConfigFR1UnpairedSpectrum.L_RA);
                x_u_fft = fft(x_u);
                x_uv_fft = fft(x_uv_avg);
                corr_fft = conj(x_uv_fft) .* x_u_fft;
                ifft_corr_fft = zeros(1, ifftLength);
                ifft_corr_fft(1:PrachConfigFR1UnpairedSpectrum.L_RA) = corr_fft;
                ifft_corr_fft = ifft(ifft_corr_fft);
                ifft_corr_fft = abs(ifft_corr_fft);
                ifft_corr_fft = circshift(ifft_corr_fft, ifftCircShift);

                %% Plot ifft_corr_fft
                if isPlot ~= 0
                    figure(u_arr(u_index))
                    plot(ifft_corr_fft);
                    title(['Sequence ' num2str(u_index)]);
                    hold on
                end
                %%
                
                %% Plot initial threshold
                maxPeak_threshold = mean(ifft_corr_fft) * constant;
                if isPlot ~= 0
                    yline(maxPeak_threshold, "LineWidth", 2, "Color", 'red');
                end
                %%
                switch PrachConfigFR1UnpairedSpectrum.preambleFormat
                    case {'0', '1', '2', '3'}
                        nifft = 2048;
                    otherwise
                        nifft = 4096;
                end
                % if u_index == 4
                %     a = 0;
                % end
                % nPreambleSeq = length(C_v_arr); % Number of windows
                nPreambleSeq = numPreablePerSeq(u_index);
                overSampling = nifft / PrachConfigFR1UnpairedSpectrum.L_RA;
                startI = ceil(overSampling * N_CS * (nPreambleSeq - 1)) + 1;
                endI = nifft;

                peakWindowSum = 0;
                numPeakWindow = 0;

                for ncv = 0:nPreambleSeq - 1

                    [maxPeakWindow, ~] = max(ifft_corr_fft(startI:endI));

                    if maxPeakWindow > maxPeak_threshold
                        peakWindowSum = peakWindowSum + maxPeakWindow;
                        numPeakWindow = numPeakWindow + 1;
                    end

                    startI = mod((endI + 1), nifft);
                    endI = ceil(overSampling * (ncv + 1) * N_CS);
                end

                %% Plot new threshold
                new_maxPeak_threshold = ((sum(ifft_corr_fft) - peakWindowSum) / (nifft - numPeakWindow)) * constant;
                % Check for new threshold, it must be larger than initial one
                if new_maxPeak_threshold <= maxPeak_threshold
                    new_maxPeak_threshold = maxPeak_threshold;
                end

                if isPlot ~= 0
                    yline(new_maxPeak_threshold, "LineWidth", 2, "Color", 'green');
                end
                %%
                windowBin_arr = [];

                startI = ceil(overSampling * N_CS * (nPreambleSeq - 1)) + 1;
                endI = nifft;

                for ncv = 0:nPreambleSeq - 1

                    % save window position
                    windowBin_arr = [windowBin_arr, startI, endI];

                    [maxPeakWindow, maxPeakWindow_idx] = max(ifft_corr_fft(startI:endI));

                    if maxPeakWindow > new_maxPeak_threshold

                        nIndex = ncv;
                        endWin = endI;
                        nTa = floor((endWin - (maxPeakWindow_idx + startI - 1)) / 8); % Timing advance

                        nTA_time_microSecond = sampleTime_microSecond * nTa;

                        start_preambleIndex = length(C_v_arr) * (u_index - 1);
                        preambleIndex = start_preambleIndex + nIndex;
                        
                        if nTA_time_microSecond <= timeErrorTolerance
                            detected_preambleIndex_arr = [detected_preambleIndex_arr preambleIndex];
                            nTA_microSecond_arr = [nTA_microSecond_arr nTA_time_microSecond];
                        end

                        switch PrachConfigFR1UnpairedSpectrum.preambleFormat
                            case {'0', '1', '2', '3'}
                                result = {'constant', constant, 'ifftCircShit', ifftCircShift, 'frame', preamble_info_arr(preamble_only_idx, 1), 'subframe', preamble_info_arr(preamble_only_idx, 2), 'preIdx', preambleIndex, 'TA', nTa};
                            otherwise
                                result = {'constant', constant, 'ifftCircShit', ifftCircShift, 'slot', preamble_info_arr(preamble_only_idx, 1), 'startSym', preamble_info_arr(preamble_only_idx, 2), 'preIdx', preambleIndex, 'TA', nTa};
                        end


                        if isPrintResult ~= 0
                            disp(result);
                        end

                        preambleIndex_timeAdvance_arr = [preambleIndex_timeAdvance_arr; result];
                    end

                    startI = mod((endI + 1), nifft);
                    endI = ceil(overSampling * (ncv + 1) * N_CS);
                end

                % Get unique values of window bin
                windowBin_arr = unique(windowBin_arr);
                windowBin_arr = windowBin_arr(1:2:end);
                % Sort increment
                windowBin_arr = sort(windowBin_arr);

                %% Plot window bins
                if isPlot ~= 0
                    for windowBin_idx = 1:length(windowBin_arr)
                        xline(windowBin_arr(windowBin_idx));
                    end
                    axis tight

                    xticks(windowBin_arr);
                    legend('Signal', 'Initial threshold', 'New threshold');
                end
            end
        end
    end
       
end

 





