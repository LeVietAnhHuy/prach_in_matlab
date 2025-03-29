function [zadoffChuSeq_slot, preambleIndex_timeAdvance_arr] = PRACH_demodulation_2(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, timeDomain_signal)
%PRACH_DEMODULATION Summary of this function goes here
%   Detailed explanation goes here
    
    timeDomain_signal = circshift(timeDomain_signal, [-1000, 0]);
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
    preambleIndex_timeAdvance_arr = [];
    zadoffChuSeq_slot = [];
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
            
            % % Define a carrier configuration object
            % carrier = nrCarrierConfig;
            % carrier.SubcarrierSpacing = 30;
            % carrier.NSizeGrid = 273;
            % 
            % % PRACH configuration
            % prach = nrPRACHConfig;
            % prach.FrequencyRange = 'FR1';   % Frequency range ('FR1', 'FR2')
            % prach.DuplexMode = 'TDD';       % Duplex mode ('FDD', 'TDD', 'SUL')
            % prach.ConfigurationIndex = 86; % Configuration index (0...255)
            % prach.SubcarrierSpacing = 30;   % Subcarrier spacing (1.25, 5, 15, 30, 60, 120)
            % prach.FrequencyIndex = 1;       % Index of the PRACH transmission occasions in frequency domain (0...7)
            % prach.TimeIndex = 0;            % Index of the PRACH transmission occasions in time domain (0...6)
            % prach.ActivePRACHSlot = 1;      % Active PRACH slot number within a subframe or a 60 kHz slot (0, 1)
            % prach.NPRACHSlot = 15;
            % prach.SequenceIndex = 39;
            % prach.ZeroCorrelationZone = 8;
            % prach.PreambleIndex = 60;
            % 
            % % Store the PRACH configuration and additional parameters in the
            % % waveconfig structure
            % waveconfig.PRACH.Config = prach;
            % waveconfig.PRACH.AllocatedPreambles = 'all'; % Index of the allocated PRACH preambles
            % waveconfig.PRACH.Power = 1;  
            % [preambleIdx, offset] = nrPRACHDetect(carrier,prach, average_preamble_only);
            % result = {'slot', slot_contain_prach(slot_index), 'starting symbol', start_symbol_contain_prach(symbol_index), 'preamble index', preambleIdx, 'time advance', offset};
            % preambleIndex_timeAdvance_arr = [preambleIndex_timeAdvance_arr; result];
            % one_prachSequence = timeDomain_signal((startingSampleSymbol + PrachOFDMInfo.CyclicPrefix):(startingSampleSymbol + PrachOFDMInfo.CyclicPrefix + nfft - 1));
            
            % preamble_only_fft = fft(average_preamble_only, nfft);
            % preamble_only_fft_fftshift = fftshift(preamble_only_fft);
            % position = nfft / 2 + prachConfig.PrachFreqStart + first;
            % y_uv = preamble_only_fft_fftshift(position:(position + prachConfig.L_RA - 1));
            % y_uv = y_uv * sqrt(length(y_uv));
            % y_uv = y_uv.';
            % x_uv = ifft(y_uv);

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
                

                currentPA = max(ifftSeq) / mean(ifftSeq);
                if(currentPA > maxPA)
                    maxPA = currentPA;
                    u_index_used = u_index;
                    maxPA_ifftSeq = fftshift(ifftSeq);
                end
            end
            u_test = u_index_used;
            u_used = u_arr(u_index_used);
            C_v_index_used = 0;

            x_u = zadoffChu(u_used, prachConfig.L_RA);
            x_uv_fft = fft(x_u);
            x_uv_detecting_fft = fft(x_uv);
            corr_fft = conj(x_uv_current_fft) .* x_uv_fft;
            window_length = ceil(length(x_uv) / length(C_v_index_used));

            maxPeak = max(corr_fft);
            maxPeakPosition = find(corr_fft == maxPeak);
            while maxPeakPosition < 1024
                maxPeakPosition = maxPeakPosition + window_length;
                
            end
            timingAdvance = (maxPeakPosition - length(C_v_index_used)) / 8;
            
             % peakThreshold = max(ifftSeq) / 2;
             %    peakThreshold = mean(ifftSeq) * 4;
             %    windowLength = floor(1024 / length(C_v_arr));
             %    timingAdvance = 0;
             %    if (max(ifftSeq(1: windowLength)) > peakThreshold) && (max(ifftSeq((end - windowLength): end)) > peakThreshold)
             %        C_v_index_used = C_v_index;
             % 
             %        ifftSeq_ifftShift = ifftshift(ifftSeq);
             %        peakPosition = find(ifftSeq_ifftShift == max(ifftSeq_ifftShift));
             % 
             %        while peakPosition < 1024
             %            peakPosition = peakPosition + windowLength;
             %        end
             % 
             %        timingAdvance = (peakPosition - 1024) / 8;
             % 
             %        break    
             %    end
             % 
             %    if(currentPA > maxPA)
             %        maxPA = currentPA;
             %        C_v_index_used = C_v_index;
             %        maxPA_ifftSeq = fftshift(ifftSeq);
             %    end

            % for C_v_index = 1:length(C_v_arr)
            % 
            %     x_u = zadoffChu(u_used, prachConfig.L_RA);
            %     x_uv_current = circshift(x_u, [0 -C_v_arr(C_v_index)]);
            % 
            %     % x_uv_current_fft = zeros(1024, 1);
            %     % x_uv_current_fft(1:prachConfig.L_RA) = x_uv_current;
            % 
            %     % x_uv_current_fft = fft(x_uv);
            %     x_uv_current_fft = fft(x_uv_current);
            % 
            %     % x_uv_fft = zeros(1024, 1);
            %     % x_uv_fft(1:prachConfig.L_RA) = x_uv;
            % 
            %     x_uv_fft = fft(x_u);
            %     % x_uv_fft = x_u;
            % 
            %     corr_fft = conj(x_uv_current_fft) .* x_uv_fft;
            %     figure(10)
            %     plot(real(ifft(corr_fft)));
            %     hold on;
            %     % plot(ifft(corr_fft));
            % 
            %     % corr_avg = zeros(1, prachConfig.L_RA);
            %     % for x_uv_index = 1:PrachConfigFR1UnpairedSpectrum.PrachDuration
            %     %     corr_avg = corr_avg + x_uv_arr(x_uv_index, :) .* conj(x_uv_current);
            %     % end
            %     % % corr = x_uv .* conj(x_u);
            %     % corr_avg = corr_avg ./ PrachConfigFR1UnpairedSpectrum.PrachDuration;
            %     % corr = x_uv .* conj(x_uv_current);
            %     % ifftSeq = zeros(1024, 1);
            %     % ifftSeq(1:139) = corr_avg;
            %     % ifftSeq = abs(ifft(ifftSeq));
            %     % ifftSeq_ifftShift = ifftshift(ifftSeq);
            %     % 
            %     % % figure(C_v_index)
            %     % % % plot(ifftSeq_ifftShift);
            %     % % % hold on 
            %     % % plot(ifftSeq);
            %     % % 
            %     % % title(C_v_arr(C_v_index));
            %     % % legend('yes','no')
            %     % 
            %     % peakThreshold = max(ifftSeq) / 2;
            %     % peakThreshold = mean(ifftSeq) * 4;
            %     % windowLength = floor(1024 / length(C_v_arr));
            %     % timingAdvance = 0;
            %     % if (max(ifftSeq(1: windowLength)) > peakThreshold) && (max(ifftSeq((end - windowLength): end)) > peakThreshold)
            %     %     C_v_index_used = C_v_index;
            %     % 
            %     %     ifftSeq_ifftShift = ifftshift(ifftSeq);
            %     %     peakPosition = find(ifftSeq_ifftShift == max(ifftSeq_ifftShift));
            %     % 
            %     %     while peakPosition < 1024
            %     %         peakPosition = peakPosition + windowLength;
            %     %     end
            %     % 
            %     %     timingAdvance = (peakPosition - 1024) / 8;
            %     % 
            %     %     break    
            %     % end
            % 
            %     % if(currentPA > maxPA)
            %     %     maxPA = currentPA;
            %     %     C_v_index_used = C_v_index;
            %     %     maxPA_ifftSeq = fftshift(ifftSeq);
            %     % end
            % end

            % figure(prachConfig.PreambleIndex + 1)
            % plot(real(maxPA_ifftSeq));

            start_preambleIndex = length(u_arr) * (u_index_used - 1);
            preambleIndex = start_preambleIndex + C_v_index_used - 1;


            result = {'slot', slot_contain_prach(slot_index), 'starting symbol', start_symbol_contain_prach(symbol_index), 'preamble index', preambleIndex, 'TA', timingAdvance};
            disp(result);
            preambleIndex_timeAdvance_arr = [preambleIndex_timeAdvance_arr; result];


            % % disp(C_V_arr(C_v_index_used));
            % % figure(3);
            % % plot(maxPA_ifftSeq);
            % threshold = mean(maxPA_ifftSeq) * 30;
            % max_C_v_index = 0;
            % pos = [];
            % timingAdvance = 0;
            % for C_v_index = 0:(length(u_arr) - 1)
            % 
            %     if C_v_index == 0
            %         startPosWindow = 1;
            %     else
            %         startPosWindow = endPosWindow + 1;
            %     end
            % 
            %     endPosWindow = startPosWindow + ceil(1024 / length(u_arr)) - 1;
            % 
            %     if C_v_index == length(u_arr) - 1
            %         endPosWindow = 1024;
            %     end
            % 
            %     max_ifftSeq_inWindow = max(maxPA_ifftSeq(startPosWindow:endPosWindow));
            %     if max_ifftSeq_inWindow > threshold
            %         max_C_v_index = C_v_index;
            %         timingAdvance = (endPosWindow - find(maxPA_ifftSeq == max_ifftSeq_inWindow)) / 8;
            %     end
            % 
            %     pos = [pos startPosWindow endPosWindow];
            % end
            % 
            % start_preambleIndex = length(u_arr) * (u_index_used - 1);
            % preambleIndex = start_preambleIndex + max_C_v_index + 1;
            % 
            % result = {'slot', slot_contain_prach(slot_index), 'starting symbol', start_symbol_contain_prach(symbol_index), 'preamble index', preambleIndex}; % 'time advance', timingAdvance};
            % disp(result);
            % figure(prachConfig.PreambleIndex + 1)
            % plot(maxPA_ifftSeq);
            % preambleIndex_timeAdvance_arr = [preambleIndex_timeAdvance_arr; result];
        end
        
        zadoffChuSeq_slot = cat(3, zadoffChuSeq_slot, zadoffChuSeq_symbol.');
    end
    
    % PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
    % preambleIndexConfig = get_Table6332x(prachConfig, carrierConfig);
    % 
    % preamble_only = preamble((PrachOFDMInfo.CyclicPrefix + 1):(PrachOFDMInfo.CyclicPrefix + preambleIndexConfig.PrachDuration * nfft));
    % preamble_only = reshape(preamble_only, [nfft, preambleIndexConfig.PrachDuration]);
    % 
    % % add noise later!
    % average_preamble_only = sum(preamble_only, 2) / preambleIndexConfig.PrachDuration;
    % 
    % preamble_only_fft = fft(average_preamble_only, nfft);
    % preamble_only_fft_fftshift = fftshift(preamble_only_fft);
    % 
    % k_bar = 2;
    % 
    % % Clause 5.3.2, TS 38.211
    % k1 = prachConfig.PrachFreqStart * 12 - carrierConfig.n_UL_RB * 6;
    % first = k1 + k_bar + 1;
    % position = nfft / 2 + 0 + first;
    % y_uv = preamble_only_fft_fftshift(position:(position + prachConfig.L_RA - 1));
    % y_uv = y_uv * sqrt(length(y_uv));
    % y_uv = y_uv.';
    % x_uv = ifft(y_uv);
    % 
    % N_CS = get_N_CS(prachConfig);
    % [~, u_arr] = get_u(prachConfig, N_CS);
    % 
    % maxPA = 0;
    % u_index_used = 0;
    % 
    % for u_index = 1:length(u_arr)
    %     x_u = zadoffChu(u_arr(u_index), prachConfig.L_RA);
    %     corr = x_uv .* conj(x_u);
    %     ifftSeq = zeros(1024, 1);
    %     ifftSeq(1:139) = corr;
    %     ifftSeq = abs(ifft(ifftSeq));
    %     currentPA = max(ifftSeq) / mean(ifftSeq);
    %     if(currentPA > maxPA)
    %         maxPA = currentPA;
    %         u_index_used = u_index;
    %         maxPA_ifftSeq = fftshift(ifftSeq);
    %     end
    % end
    % 
    % % plot(maxPA_ifftSeq);
    % threshold = mean(maxPA_ifftSeq) * 30;
    % max_C_v_index = 0;
    % pos = [];
    % for C_v_index = 0:(length(u_arr) - 1)
    % 
    %     if C_v_index == 0
    %         startPosWindow = 1;
    %     else
    %         startPosWindow = endPosWindow + 1;
    %     end
    % 
    %     endPosWindow = startPosWindow + ceil(1024 / length(u_arr)) - 1;
    % 
    %     if C_v_index == length(u_arr) - 1
    %         endPosWindow = 1024;
    %     end
    % 
    %     max_ifftSeq_inWindow = max(maxPA_ifftSeq(startPosWindow:endPosWindow));
    %     if max_ifftSeq_inWindow > threshold
    %         max_C_v_index = C_v_index;
    %         timingAdvance = (endPosWindow - find(maxPA_ifftSeq == max_ifftSeq_inWindow)) / 8;
    %     end
    % 
    %     pos = [pos startPosWindow endPosWindow];
    % end
    % 
    % start_preambleIndex = length(u_arr) * (u_index_used - 1);
    % preambleIndex = start_preambleIndex + max_C_v_index + 1;
end

