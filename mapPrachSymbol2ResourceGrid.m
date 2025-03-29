function ResourceGridWithPrachSymbol = mapPrachSymbol2ResourceGrid(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum, zadoffChuSequence)
%MAPPRACHSYMBOL2RESOURCEGRID Summary of this function goes here
%   Detailed explanation goes here
    
    numFrame = carrierConfig.numFrame;
    frame = 0:numFrame - 1;
    ResourceGridWithPrachSymbol = createResourceGrid(carrierConfig, prachConfig, PrachConfigFR1UnpairedSpectrum);
    zadoffChuSequence = zadoffChuSequence.';
    k1 = prachConfig.PrachFreqStart * 12 - carrierConfig.n_UL_RB * 6;
    first = k1 + PrachConfigFR1UnpairedSpectrum.k_bar + 1;
    switch PrachConfigFR1UnpairedSpectrum.preambleFormat
        case '0'
            for frameIndex = frame
                if mod(frameIndex, PrachConfigFR1UnpairedSpectrum.x) ~= PrachConfigFR1UnpairedSpectrum.y
                    continue;
                end

        %         if PrachConfigFR1UnpairedSpectrum.preambleFormat
                for subFrameNumber_index = 1:length(PrachConfigFR1UnpairedSpectrum.subframeNumber)
                    subFrameIndex_in_currentFrameIndex = PrachConfigFR1UnpairedSpectrum.subframeNumber(subFrameNumber_index) + 10 * frameIndex;
                    totalResourceElement_before_subframeNumber = subFrameIndex_in_currentFrameIndex;
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
        case {'1', '2'}
            for frameIndex = frame
                if mod(frameIndex, PrachConfigFR1UnpairedSpectrum.x) ~= PrachConfigFR1UnpairedSpectrum.y
                    continue;
                end
                % disp(frameIndex);

        %         if PrachConfigFR1UnpairedSpectrum.preambleFormat
                for subFrameNumber_index = 1:length(PrachConfigFR1UnpairedSpectrum.subframeNumber)
                    subFrameIndex_in_currentFrameIndex = PrachConfigFR1UnpairedSpectrum.subframeNumber(subFrameNumber_index) + 10 * frameIndex;
                    totalResourceElement_before_subframeNumber = subFrameIndex_in_currentFrameIndex;
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


        case '3'
            for frameIndex = frame
                if mod(frameIndex, PrachConfigFR1UnpairedSpectrum.x) ~= PrachConfigFR1UnpairedSpectrum.y
                    continue;
                end

        %         if PrachConfigFR1UnpairedSpectrum.preambleFormat
                for subFrameNumber_index = 1:length(PrachConfigFR1UnpairedSpectrum.subframeNumber)
                    subFrameIndex_in_currentFrameIndex = PrachConfigFR1UnpairedSpectrum.subframeNumber(subFrameNumber_index) + 10 * frameIndex;
                    totalResourceElement_before_subframeNumber = 4 * subFrameIndex_in_currentFrameIndex;
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
            
        otherwise
            for frameIndex = frame
                if mod(frameIndex, PrachConfigFR1UnpairedSpectrum.x) ~= PrachConfigFR1UnpairedSpectrum.y
                    continue;
                end

        %         if PrachConfigFR1UnpairedSpectrum.preambleFormat
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
    
    ResourceGridWithPrachSymbol = ResourceGridWithPrachSymbol;
end

