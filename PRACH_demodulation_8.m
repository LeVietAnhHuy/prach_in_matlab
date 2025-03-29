function [zadoffChuSeq_slot, preambleIndex_timeAdvance_arr, array_ifft_corr_fft, ifftSeq_array, mean_u_ifft, C_v_arr, preambleIndex_arr] = PRACH_demodulation_8(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, timeDomain_signal)
%PRACH_DEMODULATION Summary of this function goes here
%   Detailed explanation goes here

    % Number of sample of 1 slot (= 1ms) for 15kHz
    default_numSample_perSlot = 30720;
    % Clause 5.3.2, TS 38.211
    k1 = prachConfig.PrachFreqStart * 12 - carrierConfig.n_UL_RB * 6;
    first = k1 + PrachConfigFR1UnpairedSpectrum.k_bar + 1;
    
    PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

    grid_column = carrierConfig.n_UL_RB * carrierConfig.numElementPerResourceBlock;
    min_nfft = 128;

    while min_nfft < grid_column
        min_nfft = min_nfft * 2;
    end

    nfft = min_nfft;

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

    frame_mod_x = mod(0:(carrierConfig.numFrame - 1), PrachConfigFR1UnpairedSpectrum.x);
    frame_contain_prach = find(frame_mod_x == PrachConfigFR1UnpairedSpectrum.y) - 1;
    
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
    
    %% Compute Average Antenna
    num_RX_antenna = size(timeDomain_signal, 2);
    avg_timeDomainSignal = sum(timeDomain_signal, 2) / num_RX_antenna;
    timeDomain_signal = avg_timeDomainSignal; 

    %%
    preambleIndex_timeAdvance_arr = [];

    preambleIndex_arr = [];

    zadoffChuSeq_slot = [];

    array_ifft_corr_fft = [];
    ifftSeq_array = [];
    
    mean_u_ifft = [];

    numFigure = 3;
    for slot_index = 1:length(slot_contain_prach)
        startingSample = slot_contain_prach(slot_index) * numSample_perSlot + 1;
        
        zadoffChuSeq_symbol = [];
        for symbol_index = 1:length(start_symbol_contain_prach)
            startingSampleSymbol = startingSample + start_symbol_contain_prach(symbol_index) * nfft;
            
            preamble_only = timeDomain_signal((startingSampleSymbol + PrachOFDMInfo.CyclicPrefix):(startingSampleSymbol + PrachOFDMInfo.CyclicPrefix + PrachConfigFR1UnpairedSpectrum.PrachDuration * nfft - 1));
            preamble_only = reshape(preamble_only, [nfft, PrachConfigFR1UnpairedSpectrum.PrachDuration]);
            
            x_uv_arr = [];
            for prach_index = 1:PrachConfigFR1UnpairedSpectrum.PrachDuration

                preamble_only_fft = fft(preamble_only(:, prach_index), nfft);
                preamble_only_fft_fftshift = fftshift(preamble_only_fft);
                position = nfft / 2 + prachConfig.PrachFreqStart + first;
                y_uv = preamble_only_fft_fftshift(position:(position + prachConfig.L_RA - 1));
                y_uv = y_uv * sqrt(length(y_uv));
                y_uv = y_uv.';
                x_uv = ifft(y_uv);

                x_uv_arr = [x_uv_arr; x_uv];
            end
            x_uv_avg = sum(x_uv_arr) / PrachConfigFR1UnpairedSpectrum.PrachDuration;
            zadoffChuSeq_symbol = [zadoffChuSeq_symbol; x_uv_avg];

            % add noise later!
            average_preamble_only = sum(preamble_only, 2) / PrachConfigFR1UnpairedSpectrum.PrachDuration;

            N_CS = get_N_CS(prachConfig);
            [~, u_arr] = get_u(prachConfig, N_CS);
            [~, C_v_arr] = get_C_v(prachConfig, N_CS);

            for u_index = 1:length(u_arr)
       
                ifftLength = 1024;
                x_u = zadoffChu(u_arr(u_index), prachConfig.L_RA);
                x_u_fft = fft(x_u);
                x_uv_fft = fft(x_uv_avg);
                corr_fft = conj(x_uv_fft) .* x_u_fft;
                ifft_corr_fft = zeros(1, ifftLength);
                ifft_corr_fft(1:prachConfig.L_RA) = corr_fft;
                ifft_corr_fft = ifft(ifft_corr_fft);
                ifft_corr_fft = abs(ifft_corr_fft);
                ifft_corr_fft = circshift(ifft_corr_fft, -6);
                
                % %% Plot ifft_corr_fft
                % figure(u_arr(u_index))
                % plot(ifft_corr_fft);
                % hold on
                % %%
                
                constant = 10;
                
                %% Plot initial threshold
                maxPeak_threshold = mean(ifft_corr_fft) * constant;
                % yline(maxPeak_threshold, "LineWidth", 2, "Color", 'red');
                %% 

                nifft = 1024;
                nPreambleSeq = length(C_v_arr); % Number of windows
                overSampling = nifft / prachConfig.L_RA;
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
                new_maxPeak_threshold = ((sum(ifft_corr_fft) - peakWindowSum) / (1024 - numPeakWindow)) * constant;
                % yline(new_maxPeak_threshold, "LineWidth", 2, "Color", 'green');
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

                        start_preambleIndex = length(C_v_arr) * (u_index - 1);
                        preambleIndex = start_preambleIndex + nIndex;
                        

                        result = {'slot', slot_contain_prach(slot_index), 'startSym', start_symbol_contain_prach(symbol_index), 'preIdx', preambleIndex, 'TA', nTa};
                        preambleIndex_arr = [preambleIndex_arr, preambleIndex];
                        % disp(result);
                        preambleIndex_timeAdvance_arr = [preambleIndex_timeAdvance_arr; result];

                    end

                    startI = mod((endI + 1), nifft);
                    endI = ceil(overSampling * (ncv + 1) * N_CS);
                end  
                
                % get unique values of window bin
                windowBin_arr = unique(windowBin_arr);
                windowBin_arr = windowBin_arr(1:2:end);
                % sort increment
                windowBin_arr = sort(windowBin_arr);

                % %% Plot window bins
                % for windowBin_idx = 1:length(windowBin_arr)
                %     xline(windowBin_arr(windowBin_idx));
                % end
                % axis tight
                % %%

                % xticks(windowBin_arr);
                % legend('Signal', 'Initial threshold', 'New threshold');
            end
        end
        zadoffChuSeq_slot = cat(3, zadoffChuSeq_slot, zadoffChuSeq_symbol.');
    end
end

