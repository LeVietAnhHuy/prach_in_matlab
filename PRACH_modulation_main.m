function [timeDomain_signal, u] = PRACH_modulation_main(resourceGrid, prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum)

% % Number of sample of 1 slot (= 1ms) for 15kHz
% default_numSample_perSlot = 30720;
% 
% % Clause 5.3.2, TS 38.211
% k1 = prachConfig.PrachFreqStart * 12 - carrierConfig.n_UL_RB * 6;
% first = k1 + PrachConfigFR1UnpairedSpectrum.k_bar + 1;

default_numSample_perSlot = 30720;
BW = 98.28e3;
numResourceElement_freqDomain = BW / prachConfig.SubcarrierSpacing;


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

first = PrachStartingResourceElementIndex_freqDomain;

PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

%% Find size of fft based on Preamble Format
switch PrachConfigFR1UnpairedSpectrum.preambleFormat
    case '0'
        nfft = PrachOFDMInfo.Sequence;
    case {'1', '2', '3'}
        nfft = PrachOFDMInfo.Sequence / PrachConfigFR1UnpairedSpectrum.PrachDuration;
    otherwise
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
end
%%


%% Find number of subframe per frame of Long Format
switch PrachConfigFR1UnpairedSpectrum.preambleFormat
    case {'0', '2'}
        subframe_factor = 10;
    case '1'
        subframe_factor = 6;
    case '3'
        subframe_factor = 40;
end
%%


%% Calculate position where Pramble appears based on Preamble Format
switch PrachConfigFR1UnpairedSpectrum.preambleFormat
    case {'0', '1', '2', '3'}
        [~, numSubframe] = size(resourceGrid);
        % Number of sample of 1 slot for prachConfig.SubcarrierSpacing
        % numSample_perSlot = default_numSample_perSlot * (prachConfig.SubcarrierSpacing / 15);
        numSample_persubframe = default_numSample_perSlot * (carrierConfig.SubcarrierSpacing / 15) * 2;
        % numSlot_perFrame = 10 * (prachConfig.SubcarrierSpacing / 15);
        % totalSlot = carrierConfig.numFrame * numSlot_perFrame;

        frame_mod_x = mod(0:(carrierConfig.numFrame - 1), PrachConfigFR1UnpairedSpectrum.x);
        frame_contain_prach = find(frame_mod_x == PrachConfigFR1UnpairedSpectrum.y) - 1;

        allframe = 0:(numSubframe / subframe_factor  - 1);
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
end
%%


%% Generate Preamble Sequences in time domain
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
ifftin((0:(PrachConfigFR1UnpairedSpectrum.L_RA - 1)) + first) = y_uv;
% ifftin = fftshift(ifftin);
% ifftin = fftshift(ifftin);
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
%%

%% Map Preamble Sequences to uplink frames
switch PrachConfigFR1UnpairedSpectrum.preambleFormat
    case {'0', '1', '2', '3'}
        for frame = allframe
            frame_sample = zeros(10*numSample_persubframe, 1);
            if ismember(frame, frame_contain_prach)
                for mappingSubframe = PrachConfigFR1UnpairedSpectrum.subframeNumber
                    mappingSubframe_index = numSample_persubframe * mappingSubframe + 1;
                    frame_sample(mappingSubframe_index:(mappingSubframe_index + length(prachSequence) - 1)) = prachSequence; 
                end
            end

             timeDomain_signal = [timeDomain_signal; frame_sample];
        end
    otherwise
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
%%

end