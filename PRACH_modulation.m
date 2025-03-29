function [timeDomain_signal] = PRACH_modulation(resourceGrid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Number of sample of 1 slot (= 1ms) for 15kHz
default_numSample_perSlot = 30720;

% Clause 5.3.2, TS 38.211
k1 = prachConfig.PrachFreqStart * 12 - carrierConfig.n_UL_RB * 6;
first = k1 + PrachConfigFR1UnpairedSpectrum.k_bar + 1;

[grid_column, ~] = size(resourceGrid);
min_nfft = 128;

while min_nfft < grid_column
    min_nfft = min_nfft * 2;
end

nfft = min_nfft;

%     nfft_norm = 2048 * (15 / prachConfig.SubcarrierSpacing);
%     nfft = nfft_norm *

% Check for n_RA_slot, TS 38.211, section 5.3.2
if (prachConfig.SubcarrierSpacing == 30 || prachConfig.SubcarrierSpacing == 120) && PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe == 1
    n_RA_slot = 1;
else
    n_RA_slot = [0, 1];
end


switch PrachConfigFR1UnpairedSpectrum.preambleFormat
    case {'0', '1', '2'}
        [~, numSubframe] = size(resourceGrid);
        PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
        % Number of sample of 1 slot for prachConfig.SubcarrierSpacing
        % numSample_perSlot = default_numSample_perSlot * (prachConfig.SubcarrierSpacing / 15);
        numSample_persubframe = default_numSample_perSlot * (carrierConfig.SubcarrierSpacing / 15) * 2;
        % numSlot_perFrame = 10 * (prachConfig.SubcarrierSpacing / 15);
        % totalSlot = carrierConfig.numFrame * numSlot_perFrame;

        frame_mod_x = mod(0:(carrierConfig.numFrame - 1), PrachConfigFR1UnpairedSpectrum.x);
        frame_contain_prach = find(frame_mod_x == PrachConfigFR1UnpairedSpectrum.y) - 1;
        

        % get N_CS
        N_CS = get_N_CS(prachConfig, PrachConfigFR1UnpairedSpectrum);

        % Physical root sequence number (u)
        [u, ~] = get_u(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

        % Cyclic shift
        [C_v, ~] = get_C_v(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

        x_u = zadoffChu(u, PrachConfigFR1UnpairedSpectrum.L_RA);
        x_uv = circshift(x_u, [0 -C_v]);
        y_uv = fft(x_uv);

        % Simulation
        x_u_ref = fft(x_u);
        corr = conj(y_uv).*x_u_ref;
        ifftIn =  zeros(1, 1024);
        ifftIn(1:length(corr)) = corr;
        
        nfft = PrachOFDMInfo.Sequence;
        switch PrachConfigFR1UnpairedSpectrum.preambleFormat
            case {'0'}
                nfft = nfft;
            otherwise
                nfft = nfft / PrachConfigFR1UnpairedSpectrum.PrachDuration;
        end

        % Normalization
        y_uv = y_uv / sqrt(length(y_uv));

        ifftin = zeros(nfft, 1);

        % Zero padding
        ifftin(nfft / 2 + (0:(PrachConfigFR1UnpairedSpectrum.L_RA - 1)) + first) = y_uv;

        ifftout = ifft(fftshift(ifftin), nfft);

        % n Prach symbols in time domain.
        SEQ = repmat(ifftout, PrachConfigFR1UnpairedSpectrum.PrachDuration, 1);

        % Get Cyclic Prefix of Prach in time domain
        CP = SEQ((end - PrachOFDMInfo.CyclicPrefix + 1):end, 1);

        % % Get Guard Period of Prach in time domain
        % GP = zeros(PrachOFDMInfo.GuardPeriod, 1);

        % prachSequence = [CP; SEQ; GP];
        prachSequence = [CP; SEQ];
        timeDomain_signal = [];
        allframe = 0:(numSubframe / 10  - 1);
        for frame = allframe
            frame_sample = zeros(10*numSample_persubframe, 1);
            if ismember(frame, frame_contain_prach)
                for mappingSubframe = PrachConfigFR1UnpairedSpectrum.subframeNumber
                    mappingSubframe_index = numSample_persubframe * mappingSubframe + PrachConfigFR1UnpairedSpectrum.startingSymbol * 2048 * (carrierConfig.SubcarrierSpacing / 15) * 2 + 1;
                    frame_sample(mappingSubframe_index:(mappingSubframe_index + length(prachSequence) - 1)) = prachSequence;
                    % % frame = zeros(numSample_perSlot, 1);
                    % for subframe = 10*frame + (1:10)
                    %     mapping_frame = zeros(numSample_persubframe, 1);
                    %     if ismember(subframe, 10*frame + PrachConfigFR1UnpairedSpectrum.subframeNumber + 1)
                    %         mapping_frame(1:length(prachSequence)) = prachSequence;
                    %     end
                    %     timeDomain_signal = [timeDomain_signal; mapping_frame];
                    % end      
                end
            end

             timeDomain_signal = [timeDomain_signal; frame_sample];
        end
    case '3'
        [~, numSubframe] = size(resourceGrid);
        PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);
        % Number of sample of 1 slot for prachConfig.SubcarrierSpacing
        % numSample_perSlot = default_numSample_perSlot * (prachConfig.SubcarrierSpacing / 15);
        numSample_persubframe = default_numSample_perSlot * (carrierConfig.SubcarrierSpacing / 15) * 2;
        % numSlot_perFrame = 10 * (prachConfig.SubcarrierSpacing / 15);
        % totalSlot = carrierConfig.numFrame * numSlot_perFrame;

        frame_mod_x = mod(0:(carrierConfig.numFrame - 1), PrachConfigFR1UnpairedSpectrum.x);
        frame_contain_prach = find(frame_mod_x == PrachConfigFR1UnpairedSpectrum.y) - 1;
        

        % get N_CS
        N_CS = get_N_CS(prachConfig, PrachConfigFR1UnpairedSpectrum);

        % Physical root sequence number (u)
        [u, ~] = get_u(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

        % Cyclic shift
        [C_v, ~] = get_C_v(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

        x_u = zadoffChu(u, PrachConfigFR1UnpairedSpectrum.L_RA);
        x_uv = circshift(x_u, [0 -C_v]);
        y_uv = fft(x_uv);

        x_u_ref = fft(x_u);
        corr = conj(y_uv).*x_u_ref;
        ifftIn =  zeros(1, 1024);
        ifftIn(1:length(corr)) = corr;
        
        nfft = PrachOFDMInfo.Sequence;
        switch PrachConfigFR1UnpairedSpectrum.preambleFormat
            case {'0'}
                nfft = nfft;
            otherwise
                nfft = nfft / PrachConfigFR1UnpairedSpectrum.PrachDuration;
        end

        % Normalization
        y_uv = y_uv / sqrt(length(y_uv));

        ifftin = zeros(nfft, 1);

        % Zero padding
        ifftin(nfft / 2 + (0:(PrachConfigFR1UnpairedSpectrum.L_RA - 1)) + first) = y_uv;

        ifftout = ifft(fftshift(ifftin), nfft);

        % n Prach symbols in time domain.
        SEQ = repmat(ifftout, PrachConfigFR1UnpairedSpectrum.PrachDuration, 1);

        % Get Cyclic Prefix of Prach in time domain
        CP = SEQ((end - PrachOFDMInfo.CyclicPrefix + 1):end, 1);

        % % Get Guard Period of Prach in time domain
        % GP = zeros(PrachOFDMInfo.GuardPeriod, 1);

        % prachSequence = [CP; SEQ; GP];
        prachSequence = [CP; SEQ];
        timeDomain_signal = [];
        allframe = 0:(numSubframe / 40  - 1);
        for frame = allframe
            frame_sample = zeros(10*numSample_persubframe, 1);
            if ismember(frame, frame_contain_prach)
                for mappingSubframe = PrachConfigFR1UnpairedSpectrum.subframeNumber
                    mappingSubframe_index = numSample_persubframe * mappingSubframe + 1;
                    frame_sample(mappingSubframe_index:(mappingSubframe_index + length(prachSequence) - 1)) = prachSequence;
                    % % frame = zeros(numSample_perSlot, 1);
                    % for subframe = 10*frame + (1:10)
                    %     mapping_frame = zeros(numSample_persubframe, 1);
                    %     if ismember(subframe, 10*frame + PrachConfigFR1UnpairedSpectrum.subframeNumber + 1)
                    %         mapping_frame(1:length(prachSequence)) = prachSequence;
                    %     end
                    %     timeDomain_signal = [timeDomain_signal; mapping_frame];
                    % end  
                end
            end

             timeDomain_signal = [timeDomain_signal; frame_sample];
        end

    otherwise
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

        % get N_CS
        N_CS = get_N_CS(prachConfig, PrachConfigFR1UnpairedSpectrum);

        % Physical root sequence number (u)
        [u, ~] = get_u(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

        % Cyclic shift
        [C_v, ~] = get_C_v(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS);

        x_u = zadoffChu(u, PrachConfigFR1UnpairedSpectrum.L_RA);
        x_uv = circshift(x_u, [0 -C_v]);
        y_uv = fft(x_uv);
        
        x_u_ref = fft(x_u);
        corr = conj(y_uv).*x_u_ref;
        ifftIn =  zeros(1, 1024);
        ifftIn(1:length(corr)) = corr;

        % Normalization
        y_uv = y_uv / sqrt(length(y_uv));

        ifftin = zeros(nfft, 1);

        % Zero padding
        ifftin(nfft / 2 + (0:(PrachConfigFR1UnpairedSpectrum.L_RA - 1)) + first) = y_uv;

        ifftout = ifft(fftshift(ifftin), nfft);

        % n Prach symbols in time domain.
        SEQ = repmat(ifftout, PrachConfigFR1UnpairedSpectrum.PrachDuration, 1);

        PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

        % Get Cyclic Prefix of Prach in time domain
        CP = SEQ((end - PrachOFDMInfo.CyclicPrefix + 1):end, 1);

        % Get Guard Period of Prach in time domain
        GP = zeros(PrachOFDMInfo.GuardPeriod, 1);

        % Get null remain sample of 1 slot
        twoOFDM = zeros(nfft * 2, 1);
        remain = zeros(PrachOFDMInfo.PathProfile, 1);

        prachSequence = [CP; SEQ; GP];
        timeDomain_signal = [];
        for slot_index = 0:(totalSlot - 1)
            % if slot_index ~= 35
            %     continue
            % end

            slot = zeros(numSample_perSlot, 1);
            if ismember(slot_index, slot_contain_prach)
                for symbol_index = 0:13
                    if ismember(symbol_index, start_symbol_contain_prach)
                        startingSample_symbol = symbol_index  * nfft + 1;
                        slot(startingSample_symbol:(startingSample_symbol + length(prachSequence) - 1)) = prachSequence;
                    end
                end
            end
            timeDomain_signal = [timeDomain_signal; slot];
        end
end
end

