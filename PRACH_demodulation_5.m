function [zadoffChuSeq_slot, preambleIndex_timeAdvance_arr, array_ifft_corr_fft, ifftSeq_array, mean_u_ifft, pos, C_v_arr] = PRACH_demodulation_5(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, timeDomain_signal)
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
    zadoffChuSeq_slot = [];

    array_ifft_corr_fft = [];
    ifftSeq_array = [];
    
    mean_u_ifft = [];

    numFigure = 3;
    for slot_index = 1:length(slot_contain_prach)
        startingSample = slot_contain_prach(slot_index) * numSample_perSlot + 1;
        
        zadoffChuSeq_symbol = [];
        
        papr_u_ifftSeq_arr = [];

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

            maxPA = 0;
            u_index_used = 1;
                
            for u_index = 1:length(u_arr)
                x_u = zadoffChu(u_arr(u_index), prachConfig.L_RA);
                
                corr_avg = zeros(1, prachConfig.L_RA);
                for x_uv_index = 1:PrachConfigFR1UnpairedSpectrum.PrachDuration
                    corr_avg = corr_avg + x_uv_arr(x_uv_index, :) .* conj(x_u);
                end
                % corr = x_uv .* conj(x_u);
                corr_avg = corr_avg ./ PrachConfigFR1UnpairedSpectrum.PrachDuration;
                ifftSeq = zeros(1024, 1);
                ifftSeq(1:139) = corr_avg.';
                ifftSeq = abs(ifft(ifftSeq));
                %%
                ifftSeq = circshift(ifftSeq, -4);
                % numFigure = numFigure + 1;
                % figure(numFigure)
                figure(u_arr(u_index))

                mean_u_ifft = [mean_u_ifft, max(ifftSeq) / mean(ifftSeq)];
                
                papr_u_ifftSeq_arr = [papr_u_ifftSeq_arr, ifftSeq];


                plot(ifftSeq);
                xticks(ceil(1:(1024 / length(u_arr)):1024));
                axis tight

            end

            [~, correct_u_pos]  = max(mean_u_ifft);
            
            ifftSeq = papr_u_ifftSeq_arr(correct_u_pos, :);

                %%
            maxPA =0;
            % currentPA = max(ifftSeq) / mean(ifftSeq);
            currentPA = max(ifftSeq) / 1.1;

            if(currentPA > maxPA)
                maxPA = currentPA;
                u_index_used = correct_u_pos;
                % maxPA_ifftSeq = fftshift(ifftSeq);

                [~, maxPeakPos_u] = findpeaks(ifftSeq, MinPeakHeight=currentPA);
                numUE = length(maxPeakPos_u);

                ifftLength = 1024;
                u_correct = u_index_used;
                x_u = zadoffChu(u_arr(u_correct), prachConfig.L_RA);
                x_u_fft = fft(x_u);
                x_uv_fft = fft(x_uv_avg);
                corr_fft = conj(x_uv_fft) .* x_u_fft;
                ifft_corr_fft = zeros(1, ifftLength);
                ifft_corr_fft(1:prachConfig.L_RA) = corr_fft;
                ifft_corr_fft = ifft(ifft_corr_fft);
                ifft_corr_fft = abs(ifft_corr_fft);
                ifft_corr_fft = circshift(ifft_corr_fft, -3);

                %%
                array_ifft_corr_fft = [array_ifft_corr_fft; ifft_corr_fft];
                %%

                % windowLength = ceil(ifftLength / length(C_v_arr));
                maxPeak_threshold = max(ifft_corr_fft) / 2;
                [~, maxPeakPos_C_v] = findpeaks(ifft_corr_fft, MinPeakHeight=maxPeak_threshold);
                % [maxPeak, maxPeakPosition] = max(ifft_corr_fft);
                
                for maxPeakPos_C_v_Index = 1:length(maxPeakPos_C_v)
                    pos = [];
                    nifft = 1024;
                    nPreambleSeq = length(u_arr); % Number of windows
                    overSampling = nifft / prachConfig.L_RA;
                    startI = ceil(overSampling * N_CS * (nPreambleSeq - 1)) + 1;
                    endI = nifft;
                    nIndex = 0;
                    endWin = 0;
                    for ncv = 0:nPreambleSeq-1
                        if ((maxPeakPos_C_v(maxPeakPos_C_v_Index) >= startI) && (maxPeakPos_C_v(maxPeakPos_C_v_Index) <= endI))
                            nIndex = ncv;
                            endWin = endI;
                        end
                        startI = mod((endI + 1), nifft);
                        endI = ceil(overSampling * (ncv + 1) * N_CS);
                    end
                    nTa = floor((endWin - maxPeakPos_C_v(maxPeakPos_C_v_Index)) / 8); % Timing advance

                    pos = ceil(1:(ifftLength / length(u_arr)):ifftLength);
                    if ~ismember(ifftLength, pos)
                        pos = [pos ifftLength];
                    end

                    start_preambleIndex = length(u_arr) * (u_index_used - 1);
                    preambleIndex = start_preambleIndex + nIndex;


                    result = {'slot', slot_contain_prach(slot_index), 'startSym', start_symbol_contain_prach(symbol_index), 'preIdx', preambleIndex, 'TA', nTa};
                    disp(result);
                    preambleIndex_timeAdvance_arr = [preambleIndex_timeAdvance_arr; result];
                
                end
                ifftSeq_array = [ifftSeq_array, ifftSeq];
            end
            

            % ifftLength = 1024;
            % u_correct = u_index_used;
            % x_u = zadoffChu(u_arr(u_correct), prachConfig.L_RA);
            % x_u_fft = fft(x_u);
            % x_uv_fft = fft(x_uv_avg);
            % corr_fft = conj(x_uv_fft) .* x_u_fft;
            % ifft_corr_fft = zeros(1, ifftLength);
            % ifft_corr_fft(1:prachConfig.L_RA) = corr_fft;
            % ifft_corr_fft = ifft(ifft_corr_fft);
            % ifft_corr_fft = abs(ifft_corr_fft);
            % ifft_corr_fft = circshift(ifft_corr_fft, -3);
            % 
            % %%
            % array_ifft_corr_fft = [array_ifft_corr_fft; ifft_corr_fft];
            % %%
            % 
            % % windowLength = ceil(ifftLength / length(C_v_arr));
            % [maxPeak, maxPeakPosition] = max(ifft_corr_fft);
            % 
            % pos = [];
            % nifft = 1024;
            % nPreambleSeq = length(u_arr); % Number of windows
            % overSampling = nifft / prachConfig.L_RA;
            % startI = ceil(overSampling * N_CS * (nPreambleSeq-1)) + 1;
            % endI = nifft;
            % nIndex = 0;
            % endWin = 0;
            % for ncv = 0:nPreambleSeq-1
            %     if ((maxPeakPosition >= startI) && (maxPeakPosition <= endI))
            %         nIndex = ncv;
            %         endWin = endI;
            %     end
            %     startI = mod((endI + 1), nifft);
            %     endI = ceil(overSampling * (ncv + 1) * N_CS);
            % end
            % nTa = floor((endWin - maxPeakPosition) / 8); % Timing advance
            % 
            % pos = ceil(1:(ifftLength / length(u_arr)):ifftLength);
            % if ~ismember(ifftLength, pos) 
            %     pos = [pos ifftLength];
            % end
            % 
            % start_preambleIndex = length(u_arr) * (u_index_used - 1);
            % preambleIndex = start_preambleIndex + nIndex;
            % 
            % 
            % result = {'slot', slot_contain_prach(slot_index), 'startSym', start_symbol_contain_prach(symbol_index), 'preIdx', preambleIndex, 'TA', nTa};
            % disp(result);
            % preambleIndex_timeAdvance_arr = [preambleIndex_timeAdvance_arr; result];

        end
        
        zadoffChuSeq_slot = cat(3, zadoffChuSeq_slot, zadoffChuSeq_symbol.');
    end
end


