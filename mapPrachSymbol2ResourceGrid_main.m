function ResourceGridWithPrachSymbol = mapPrachSymbol2ResourceGrid_main(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, zadoffChuSequence)
%MAPPRACHSYMBOL2RESOURCEGRID Summary of this function goes here
%   Detailed explanation goes here
    
    numFrame = carrierConfig.numFrame;
    frame = 0:numFrame - 1;
    ResourceGridWithPrachSymbol = createResourceGrid(carrierConfig, prachConfig, PrachConfigFR1UnpairedSpectrum);
    zadoffChuSequence = zadoffChuSequence.';
    k1 = prachConfig.PrachFreqStart * 12 - carrierConfig.n_UL_RB * 6;
    first = k1 + PrachConfigFR1UnpairedSpectrum.k_bar + 1;
    
    %% Find number of symbol per subframe
    switch PrachConfigFR1UnpairedSpectrum.preambleFormat
        case {'0', '1', '2'}
            numSybolperSubframe = 1;
        case '3'
            numSybolperSubframe = 4;
    end
    %%
    
    %% Map Preamble Sequence onto resource grid
    switch PrachConfigFR1UnpairedSpectrum.preambleFormat
        case {'0', '3'}
            for frameIndex = frame
                if mod(frameIndex, PrachConfigFR1UnpairedSpectrum.x) ~= PrachConfigFR1UnpairedSpectrum.y
                    continue;
                end

        %         if PrachConfigFR1UnpairedSpectrum.preambleFormat
                for subFrameNumber_index = 1:length(PrachConfigFR1UnpairedSpectrum.subframeNumber)
                    subFrameIndex_in_currentFrameIndex = PrachConfigFR1UnpairedSpectrum.subframeNumber(subFrameNumber_index) + 10 * frameIndex;
                    totalResourceElement_before_subframeNumber = numSybolperSubframe * subFrameIndex_in_currentFrameIndex;
                    if isnan(PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe)
                        numPrachSlotsWithinASubframe = 0;
                    else
                        numPrachSlotsWithinASubframe = PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe;
                    end
                    % PrachStartingResourceElementIndex_timeDomain = totalResourceElement_before_subframeNumber + 14 * numPrachSlotsWithinASubframe + PrachConfigFR1UnpairedSpectrum.startingSymbol + 1;
                    PrachStartingResourceElementIndex_timeDomain = totalResourceElement_before_subframeNumber + 14 * numPrachSlotsWithinASubframe;
                    PrachStartingResourceElementIndex_freqDomain = abs(first);

                    % 12 identical columns of Zadoff-Chu Sequences
                    PrachMatrix = repmat(zadoffChuSequence, 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
                    ResourceGridWithPrachSymbol(PrachStartingResourceElementIndex_freqDomain:(PrachStartingResourceElementIndex_freqDomain + PrachConfigFR1UnpairedSpectrum.L_RA - 1), PrachStartingResourceElementIndex_timeDomain + (1:PrachConfigFR1UnpairedSpectrum.PrachDuration)) = PrachMatrix;
                end
            end
        case '1'
            for frameIndex = frame
                if mod(frameIndex, PrachConfigFR1UnpairedSpectrum.x) ~= PrachConfigFR1UnpairedSpectrum.y
                    continue;
                end

        %         if PrachConfigFR1UnpairedSpectrum.preambleFormat
                for subFrameNumber_index = 1:length(PrachConfigFR1UnpairedSpectrum.subframeNumber)
                    mappingSubframe_otherFormats = 10 * frameIndex + (PrachConfigFR1UnpairedSpectrum.subframeNumber(subFrameNumber_index) + 1);
                    mappingSubframe_format1 = mappingSubframe_otherFormats * 2 / 3;
                    if mappingSubframe_format1 == fix(mappingSubframe_format1)
                        mappingSubframe_format1 = mappingSubframe_format1 - 1;
                    else
                        mappingSubframe_format1 = round(mappingSubframe_format1);
                    end

                  
                    totalResourceElement_before_subframeNumber = numSybolperSubframe * mappingSubframe_format1;
                    if isnan(PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe)
                        numPrachSlotsWithinASubframe = 0;
                    else
                        numPrachSlotsWithinASubframe = PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe;
                    end
                    % PrachStartingResourceElementIndex_timeDomain = totalResourceElement_before_subframeNumber + 14 * numPrachSlotsWithinASubframe + PrachConfigFR1UnpairedSpectrum.startingSymbol + 1;
                    PrachStartingResourceElementIndex_timeDomain = totalResourceElement_before_subframeNumber + 14 * numPrachSlotsWithinASubframe;
                    PrachStartingResourceElementIndex_freqDomain = abs(first);

                    % 12 identical columns of Zadoff-Chu Sequences
                    PrachMatrix = repmat(zadoffChuSequence, 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
                    ResourceGridWithPrachSymbol(PrachStartingResourceElementIndex_freqDomain:(PrachStartingResourceElementIndex_freqDomain + PrachConfigFR1UnpairedSpectrum.L_RA - 1), PrachStartingResourceElementIndex_timeDomain + (0:(PrachConfigFR1UnpairedSpectrum.PrachDuration - 1))) = PrachMatrix;
                end
            end
        case '2'
            for frameIndex = frame
                if mod(frameIndex, PrachConfigFR1UnpairedSpectrum.x) ~= PrachConfigFR1UnpairedSpectrum.y
                    continue;
                end

        %         if PrachConfigFR1UnpairedSpectrum.preambleFormat
                for subFrameNumber_index = 1:length(PrachConfigFR1UnpairedSpectrum.subframeNumber)
                    mappingSubframe = 10 * frameIndex + (PrachConfigFR1UnpairedSpectrum.subframeNumber(subFrameNumber_index) + 1);
                    if mod(mappingSubframe, 4) == 0
                        upperBound = mappingSubframe;
                        lowerBound = upperBound - 3;
                    elseif mod(mappingSubframe, 4) == 1
                        lowerBound = mappingSubframe;
                        upperBound = lowerBound + 3;
                    else
                        remainder = mod(mappingSubframe, 4);
                        lowerBound = mappingSubframe - (remainder - 1);
                        upperBound = lowerBound + 3;
                    end
                    mappingSubframe = lowerBound;

                    totalResourceElement_before_subframeNumber = numSybolperSubframe * mappingSubframe;
                    if isnan(PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe)
                        numPrachSlotsWithinASubframe = 0;
                    else
                        numPrachSlotsWithinASubframe = PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe;
                    end
                    % PrachStartingResourceElementIndex_timeDomain = totalResourceElement_before_subframeNumber + 14 * numPrachSlotsWithinASubframe + PrachConfigFR1UnpairedSpectrum.startingSymbol + 1;
                    PrachStartingResourceElementIndex_timeDomain = totalResourceElement_before_subframeNumber + 14 * numPrachSlotsWithinASubframe;
                    PrachStartingResourceElementIndex_freqDomain = abs(first);

                    % 12 identical columns of Zadoff-Chu Sequences
                    PrachMatrix = repmat(zadoffChuSequence, 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
                    ResourceGridWithPrachSymbol(PrachStartingResourceElementIndex_freqDomain:(PrachStartingResourceElementIndex_freqDomain + PrachConfigFR1UnpairedSpectrum.L_RA - 1), PrachStartingResourceElementIndex_timeDomain + (0:(PrachConfigFR1UnpairedSpectrum.PrachDuration - 1))) = PrachMatrix;
                end
            end
        otherwise
            for frameIndex = frame
                if mod(frameIndex, PrachConfigFR1UnpairedSpectrum.x) ~= PrachConfigFR1UnpairedSpectrum.y
                    continue;
                end

                for subFrameNumber_index = 1:length(PrachConfigFR1UnpairedSpectrum.subframeNumber)
                    subFrameIndex_in_currentFrameIndex = PrachConfigFR1UnpairedSpectrum.subframeNumber(subFrameNumber_index) + 10 * frameIndex;
                    totalResourceElement_before_subframeNumber = 14 * 2 * subFrameIndex_in_currentFrameIndex;
                    if isnan(PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe)
                        numPrachSlotsWithinASubframe = 0;
                    else
                        numPrachSlotsWithinASubframe = PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe;
                    end
                    PrachStartingResourceElementIndex_timeDomain = totalResourceElement_before_subframeNumber + 14 * numPrachSlotsWithinASubframe + PrachConfigFR1UnpairedSpectrum.startingSymbol + 1;
                    PrachStartingResourceElementIndex_freqDomain = abs(first);

                    % 12 identical columns of Zadoff-Chu Sequences
                    PrachMatrix = repmat(zadoffChuSequence, 1, PrachConfigFR1UnpairedSpectrum.PrachDuration);
                    ResourceGridWithPrachSymbol(PrachStartingResourceElementIndex_freqDomain:(PrachStartingResourceElementIndex_freqDomain + PrachConfigFR1UnpairedSpectrum.L_RA - 1), PrachStartingResourceElementIndex_timeDomain:(PrachStartingResourceElementIndex_timeDomain + PrachConfigFR1UnpairedSpectrum.PrachDuration - 1)) = PrachMatrix;
                end
            end
    end
    %%

    ResourceGridWithPrachSymbol = ResourceGridWithPrachSymbol;
end