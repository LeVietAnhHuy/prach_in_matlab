function ResourceGrid_in_ResourceBlock = ResourceGrid_in_ResourceElement2ResourceBlock(ResourceGrid_in_ResourceELement, frameDisplay)
%PLOTRESOURCEGRID_IN_RESOURCEBLOCL Summary of this function goes here
%   Detailed explanation goes here
    
    [numRE_inFreq, numRE_inTime] = size(ResourceGrid_in_ResourceELement);

    ResourceGrid_in_ResourceBlock = zeros(numRE_inFreq / 12, numRE_inTime);
    
    for RE_index_inFreq = 1:12:numRE_inFreq
        RB_index_inFreq = fix(RE_index_inFreq / 12);
  
        for frame_index_inTime = 0:(numRE_inTime / 280 - 1)
            if not(ismember(frame_index_inTime, frameDisplay))
                continue
            end
            
            for displayResourceELement_inTime = (frame_index_inTime * 280 + 1):(frame_index_inTime * 280 + 280)
                if ResourceGrid_in_ResourceELement(RE_index_inFreq + 1, displayResourceELement_inTime) ~= 0
                    ResourceGrid_in_ResourceBlock(RB_index_inFreq + 1, displayResourceELement_inTime) = 200;
                end
            end
        end
    end
end

