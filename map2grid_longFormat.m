function grid = map2grid_longFormat(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, y_uv)
%MAP2GRID_LONGFORMAT Summary of this function goes here
%   Detailed explanation goes here
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

    %%
    % k1 = prachConfig.PrachFreqStart * 12 - carrierConfig.n_UL_RB * 6;
    % first = k1 + PrachConfigFR1UnpairedSpectrum.k_bar + 1;
    % PrachStartingResourceElementIndex_freqDomain = abs(first);

    zadoffChuSequence = y_uv.';
    PrachMatrix = repmat(zadoffChuSequence, 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
    grid_block_empty = zeros(numResourceElement_freqDomain, PrachConfigFR1UnpairedSpectrum.PrachDuration);
    grid_block_prach = grid_block_empty;
    grid_block_prach(PrachStartingResourceElementIndex_freqDomain + (1:839), 1:PrachConfigFR1UnpairedSpectrum.PrachDuration) = PrachMatrix;

    samplingRate_ratio = (carrierConfig.SubcarrierSpacing / 15) * 2;

    PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum);

    
    % grid = createResourceGrid(carrierConfig, prachConfig, PrachConfigFR1UnpairedSpectrum);
    blockSize = PrachConfigFR1UnpairedSpectrum.PrachDuration;

    samplesPerSubframe = 30720; % in units of T_s
    SubframesPerPRACHSlot = ceil(((PrachOFDMInfo.Sequence / samplingRate_ratio) + (PrachOFDMInfo.CyclicPrefix / samplingRate_ratio)) / samplesPerSubframe);

    % SubframesPerPRACHSlot = 2;
    PRACHSlotsPerPeriod = lcm(PrachConfigFR1UnpairedSpectrum.x * 10, round(max([1 SubframesPerPRACHSlot]))) / SubframesPerPRACHSlot;
    % PRACHSlotsPerPeriod = lcm(PrachConfigFR1UnpairedSpectrum.x * 10, round(max([1 carrierConfig.numFrame]))) / SubframesPerPRACHSlot;
    % Determine the starting subframes of the nominal PRACH slots    
    % nominalSFs = [0 cumsum(repmat(SubframesPerPRACHSlot, 1, PRACHSlotsPerPeriod))];
    % nominalSFs = [0 cumsum(repmat(SubframesPerPRACHSlot, 1, PrachConfigFR1UnpairedSpectrum.subframeNumber * 10))];
    nominalSFs = [0 cumsum(repmat(SubframesPerPRACHSlot, 1, carrierConfig.numFrame * 10))];

    y = PrachConfigFR1UnpairedSpectrum.y;
    x = PrachConfigFR1UnpairedSpectrum.x;
    startSFs = PrachConfigFR1UnpairedSpectrum.subframeNumber;
    startSFs = startSFs + (y * 10);
    % startSFs = startSFs(:) + (x * 10 * (0:(SubframesPerPRACHSlot - 1)));
    startSFs = startSFs(:) + (x * 10 * (0:(carrierConfig.numFrame - 1)));
    startSFs = startSFs(:).';

    % nominalSFs = [0 cumsum(repmat(SubframesPerPRACHSlot, 1, length(startSFs)))];

    % Establish which nominal PRACH slots are active, and their 
    % corresponding PRACH slot indices
    activeSlotFn = @(n)any(startSFs >= n & startSFs <= (n + SubframesPerPRACHSlot - 1));
    activePRACHSlot = arrayfun(activeSlotFn, nominalSFs);
    activeIndex = find(activePRACHSlot);

    % activeSlotInGrid = blockSize * activeIndex - 1;
        
    % for activeSlot = activeSlotInGrid
    %     grid(PrachStartingResourceElementIndex_freqDomain:(PrachStartingResourceElementIndex_freqDomain + PrachConfigFR1UnpairedSpectrum.L_RA - 1), activeSlot + (0:PrachConfigFR1UnpairedSpectrum.PrachDuration - 1)) = PrachMatrix;
    % end

    numPRACHSlots = floor(carrierConfig.numFrame * 10 / SubframesPerPRACHSlot);
    grid = [];

    for PRACHSlot = 1:numPRACHSlots
        if ismember(PRACHSlot, activeIndex)
            grid = [grid, grid_block_prach];
        else
            grid = [grid, grid_block_empty];
        end
    end   
end

