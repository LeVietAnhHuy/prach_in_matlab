function resourceGrid = createResourceGrid(carrierConfig, prachConfig, PrachConfigFR1UnpairedSpectrum)
%CREATERESOURCEGRID Summary of this function goes here
%   Detailed explanation goes here
    BW = 98.28e3;

    numSubFrame = carrierConfig.numFrame * 10;
    if PrachConfigFR1UnpairedSpectrum.preambleFormat == '1'
        % numSubFrame = carrierConfig.numFrame * 6;
        numSubFrame = round(carrierConfig.numFrame * 10 * (2 / 3));
        
    end

    
    % numResourceElement_freqDomain = carrierConfig.n_UL_RB * carrierConfig.numElementPerResourceBlock * (30 / prachConfig.SubcarrierSpacing) * (carrierConfig.SubcarrierSpacing / 15);
    numResourceElement_freqDomain = BW / prachConfig.SubcarrierSpacing;
    switch PrachConfigFR1UnpairedSpectrum.preambleFormat
        case {'0', '1', '2'}
            resourceGrid = zeros(numResourceElement_freqDomain, numSubFrame);
        case '3'
            resourceGrid = zeros(numResourceElement_freqDomain, numSubFrame * 4);
        otherwise
            % For Normal Cylic Prefix
            numResourceElement_per_Slot = 14;
            mui = log2(carrierConfig.SubcarrierSpacing / 15);

            numSlot = numSubFrame * 2^mui;
            numResourceElement_timeDomain = numResourceElement_per_Slot * numSlot;

            resourceGrid = zeros(numResourceElement_freqDomain, numResourceElement_timeDomain);
    end
end

