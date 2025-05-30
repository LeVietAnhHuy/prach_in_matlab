function PrachConfigFR1UnpairedSpectrum = get_Table6332x(prachConfig, carrierConfig)
%GETTABLE6332_3 Summary of this function goes here
%   Detailed explanation goes here
    
    switch prachConfig.FrequencyRange
        case 'FR1'
            switch prachConfig.SpectrumType
                case 'Unpaired'
                    preambleFormat_column = [string(zeros(1, 28)), ...
                                             string(ones(1, 6)), ...
                                             string(ones(1, 6) .* 2), ...
                                             string(ones(1, 27) .* 3), ...
                                             "A" + string(ones(1, 20)), ...
                                             "A" + string(ones(1, 23) .* 2), ...
                                             "A" + string(ones(1, 23) .* 3), ...
                                             "B" + string(ones(1, 12)), ...
                                             "B" + string(ones(1, 24) .* 4), ...
                                             "C" + string(zeros(1, 20)), ...
                                             "C" + string(ones(1, 22) .* 2), ...
                                             "A1/B" + string(ones(1, 15)), ...
                                             "A2/B" + string(ones(1, 15) .* 2), ...
                                             "A3/B" + string(ones(1, 15) .* 3), ...
                                             string(zeros(1, 7))];
                    
                    x_column = [16, 8, 4, ones(1, 4) .* 2, ones(1, 21), ...
                                16, 8, 4, ones(1, 2) .* 2, 1, 16, 8, 4, ...
                                ones(1, 2) .* 2, 1, 16, 8, 4, ones(1, 4) .* 2, ...
                                ones(1, 20), 16, 8, 4, ones(1, 7) .* 2, ...
                                ones(1, 10), 16, 8, 4, ones(1, 6) .* 2, ...
                                ones(1, 2), 2, ones(1, 11), 16, 8, 4, ...
                                ones(1, 6) .* 2, ones(1, 2), 2, ones(1, 11), ...
                                4, ones(1, 4) .* 2, ones(1, 7), 16, 8, 4, ... 
                                ones(1, 7) .* 2, ones(1, 14), 16, 8, 4, ...
                                ones(1, 7) .* 2, ones(1, 10), 16, 8, 4, ...
                                ones(1, 7) .* 2, 8, 4, ones(1, 10), ...
                                ones(1, 6) .* 2, ones(1, 9), ones(1, 5) .* 2, ...
                                ones(1, 10), ones(1, 6) .* 2, ones(1, 9), ...
                                16, 8, 4, ones(1, 4) .* 2];

                    y_column = [ones(1, 3), 0, 1, 0, 1, zeros(1, 21), ...
                                ones(1, 3), 0, 1, 0, ones(1, 3), 0, 1, 0, ...
                                ones(1, 3), 0, 1, 0, 1, zeros(1, 20), ...
                                ones(1, 10), zeros(1, 10), ones(1, 9), ...
                                zeros(1, 2), 1, zeros(1, 11), ones(1, 9), ...
                                zeros(1, 2), 1, zeros(1, 11), ones(1, 5), ...
                                zeros(1, 7), ones(1, 10), zeros(1, 14), ...
                                ones(1, 10), zeros(1, 10), ones(1, 12), ...
                                zeros(1, 10), ones(1, 6), zeros(1, 9), ...
                                ones(1, 5), zeros(1, 10), ones(1, 6), ...
                                zeros(1, 9), ones(1, 3), 0, 1, 0, 1];

                    subframeNumber_column = [{9}, 9, 9, 9, 9, 4, 4, 9, 8, ...
                                             7, 6, 5, 4, 3, 2, [1, 6], ...
                                             [1, 6], [4, 9], [3, 8], [2, 7], ...
                                             [8, 9], [4, 8, 9], [3, 4, 9], ...
                                             [7, 8, 9], [3, 4, 8, 9], 6:9, ...
                                             [1, 4, 6, 9], 1:2:9, 7, 7, 7, ...
                                             7, 7, 7, 6, 6, 6, 6, 6, 6, 9, ...
                                             9, 9, 9, 9, 4, 4, 9, 8, 7, 6, ...
                                             5, 4, 3, 2, [1, 6], [1, 6], ...
                                             [4, 9], [3, 8], [2, 7], [8, 9], ...
                                             [4, 8, 9], [3, 4, 9], [7, 8, 9], ...
                                             [3, 4, 8, 9], [1, 4, 6, 9], 1:2:9, ...
                                             9, 9, 9, 9, [4, 9], [7, 9], ...
                                             [7, 9], [8, 9], [4, 9], [2:4, 7:9], ...
                                             9, 9, 9, [8, 9], [4, 9], [7, 9], ...
                                             [3, 4, 8, 9], [3, 4, 8, 9], ...
                                             1:2:9, 0:9, 9, 9, 9, [7, 9], [8, 9], ...
                                             [7, 9], [4, 9], [4, 9], [2:4, 7:9], ...
                                             2, 7, 9, 9, 9, 9, [2, 7], [8, 9], ...
                                             [4, 9], [7, 9], [3, 4, 8, 9], ...
                                             [3, 4, 8, 9], 1:2:9, 0:9, 9, 9, 9, ...
                                             [4, 9], [7, 9], [7, 9], [4, 9], ...
                                             [8, 9], [2:4, 7:9], 2, 7, 9, 9, 9, 9, ...
                                             [2, 7], [8, 9], [4, 9], [7, 9], ...
                                             [3, 4, 8, 9], [3, 4, 8, 9], ...
                                             1:2:9, 0:9, 9, 9, [7, 9], [4, 9], ...
                                             [4, 9], 9, 9, 9, [8, 9], [4, 9], ...
                                             [7, 9], 1:2:9, 9, 9, 9, 9, 9, ...
                                             [7, 9], [4, 9], [4, 9], [8, 9], ...
                                             [2:4, 7:9], 1, 2, 4, 7, 9, 9, 9, ...
                                             [4, 9], [7, 9], [8, 9], [3, 4, 8, 9], ...
                                             1:2:9, 0:9, 0:9, 9, 9, 9, 9, ...
                                             [8, 9], [7, 9], [7, 9], [4, 9], ...
                                             [4, 9], [2:4, 7:9], 9, 9, 9, ...
                                             [8, 9], [4, 9], [7, 9], [3, 4, 8, 9], ...
                                             [3, 4, 8, 9], 1:2:9, 0:9, 9, 9, 9, 9, ...
                                             [8, 9], [7, 9], [7, 9], [4, 9], ...
                                             [4, 9], [2:4, 7:9], 9, 9, 9, 9, ...
                                             9, [8, 9], [4, 9], [7, 9], [3, 4, 8, 9], ...
                                             [3, 4, 8, 9], [2:4, 7:9], 0:9, ...
                                             9, [4, 9], [7, 9], [7, 9], [4, 9], ...
                                             [8, 9], 9, 9, 9, [8, 9], [4, 9], [7, 9], ...
                                             [3, 4, 8, 9], 1:2:9, 0:9, 9, ...
                                             [4, 9], [7, 9], [4, 9], [8, 9], ...
                                             9, 9, 9, [8, 9], [4, 9], [7, 9], ...
                                             [3, 4, 8, 9], [3, 4, 8, 9], 1:2:9, ...
                                             0:9, 9, [4, 9], [7, 9], [7, 9], [4, 9], ...
                                             [8, 9], 9, 9, 9, [8, 9], [4, 9], [7, 9], ...
                                             [3, 4, 8, 9], 1:2:9, 0:9, 7, 7, ...
                                             7, 7, 7, 2, 2];

                    startingSymbol_column = [zeros(1, 16), 7, zeros(1, 20), ...
                                             7, 7, 7, zeros(1, 16), 7, zeros(1, 14), ...
                                             7, 7, zeros(1, 5), 7, zeros(1, 3), ...
                                             7, zeros(1, 3), 7, zeros(1, 5), ...
                                             9, 9, zeros(1, 6), 9, zeros(1, 4), ...
                                             9, zeros(1, 3), 9, zeros(1, 3), ...
                                             7, 7, zeros(1, 8), 7, zeros(1, 4), ...
                                             7, zeros(1, 3), 7, 2, 2, 2, 8, ...
                                             2, 2, 8, 2, 2, 2, 8, 2, 0, 0, 2, 0, ...
                                             2, 2, 2, zeros(1, 8), 2, 0, 2, 2, ...
                                             0, 2, 2, 0, ones(1, 7) .* 2, 8, 8, ...
                                             2, 2, 2, 8, 2, 2, 2, 8, 2, 2, 2, ...
                                             8, ones(1, 6) .* 2, 8, 8, 2, 2, ...
                                             8, 8, 2, 8, 2, 2, 2, 8, 2, 2, 2, ...
                                             8, 2, 8, 8, 2, 2, 2, 2, 8, 2, 2, 2, 8, ...
                                             2, 2, 8, 0, 6, 6, 0, 0, 0, 6, 0, 0, 0, ...
                                             6, 0, 0, 0, 6, 0, 2, 0, 2, 0, 0, 0, ...
                                             2, 0, 0, 0, 2, 0, 0, 2, zeros(1, 7)];

                    numberPrachSlotWithinSubframe_column = [NaN(1, 67), 2, 2, ones(1, 5), 2, 2, ...
                                                            1, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, ...
                                                            2, 1, 1, 2, 1, 1, 2, ones(1, 4), 2, ...
                                                            1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2, ...
                                                            1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, ...
                                                            1, 1, 2, 1, 1, 1, 2, ones(1, 6), 2, ...
                                                            2, 1, 1, 2, 1, 1, 1, 2, 2, ones(1, 5), ...
                                                            2, 2, ones(1, 7), 2, 1, 1, 2, 1, 1, 2, ...
                                                            1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, ...
                                                            1, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, ...
                                                            1, 1, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 1, ...
                                                            1, 2, ones(1, 6), 2, 2, 2, 1, 1, 2, 1, ...
                                                            1, 2, ones(1, 5), 2, 2, 2, 1, 1, 2, 1, ...
                                                            1, 1, 2, ones(1, 6), 2, 2, 2, 1, 1, 2, ...
                                                            1, 1, 2, 1, 1, NaN(1, 7)];

                    numberTimeDomainPrachOccasionsWithinPrachSlot_column = [NaN(1, 67), 6, 6, 6, 6, 3, 3, ones(1, 5) .* 6, ...
                                                                            3, 6, 6, 6, 3, 6, 6, 6, ones(1, 6) .* 3, 1, 1, ...
                                                                            ones(1, 6) .* 3, 1, 3, 3, 3, 3, 1, 3, 3, 3, 1, ...
                                                                            2, 2, 2, 1, 1, ones(1, 8) .* 2, 1, 2, 2, 2, 2, ...
                                                                            1, 2, 2, 2, 1, 6, 6, 6, 3, 6, 6, 3, 6, 6, 6, 3, ...
                                                                            6, ones(1, 24), ones(1, 6) .* 6, 3, 3, 6, 6, 6, ...
                                                                            3, 6, 6, 6, 3, 6, 6, 6, 3, ones(1, 6) .* 2, 1, ...
                                                                            1, 2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, ...
                                                                            6, 3, 3, 6, 6, 6, 6, 3, 6, 6, 6, 3, 6, 6, 3, ...
                                                                            3, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, ...
                                                                            ones(1, 16) .* 2, NaN(1, 7)];

                    prachDuration_column = [zeros(1, 67), ones(1, 20) .* 2,  ones(1, 23) .* 4, ...
                                            ones(1, 23) .* 6, ones(1, 12) .* 2, ones(1, 24) .* 12, ...
                                            ones(1, 20) .* 2, ones(1, 22) .* 6, ones(1, 15) .* 2, ...
                                            ones(1, 15) .* 4, ones(1, 15) .* 6, zeros(1, 7)];
                
                case {'Pair', 'Supplementary'}

                otherwise 
                    warning('"prachConfig.spectrumType" must be one of {"Unpair", "Pair", "Supplementary"}}');
                    return;
            end
        case 'FR2'

        otherwise
           warning('"prachConfig.frequencyRange" must be one of {"FR1", "FR2"}}'); 
           return;
    end
    
    %% Table 6.3.3.2-1: Supported combinations of  Δf^RA and Δf , and the
    % corresponding value of k bar
    L_RA_column  = [ones(1, 6) .* 839, ...
                    ones(1, 10) .* 139, ...
                    ones(1, 3) .* 571, ...
                    ones(1, 3) .* 1151];

    delta_f_RA_forPrach_column = [ones(1, 3) .* 1.25, ones(1, 3) .* 5, ...
                                  ones(1, 3) .* 15, ones(1, 3) .* 30, ones(1, 2) .* 60, ones(1, 2) .* 120, ...
                                  ones(1, 3) .* 30, ...
                                  ones(1, 3) .* 15];

    delta_f_forPUSCH_column = [15, 30, 60, 15, 30, 60, ...
                               15, 30, 60, 15, 30, 60, 60, 120, 60, 120, ...
                               15, 30, 60, ...
                               15, 30, 60];

    % allocation expressed in number of RBs for PUSCH 
    N_RA_RB_column = [6, 3, 2, 24, 12, 6, ...
                      12, 6, 3, 24, 12, 6, 12, 6, 24, 12, ...
                      96, 48, 24, ...
                      96, 48, 24];

    k_bar_column = [7, 1, 133, 12, 10, 7, ...
                    ones(1, 10) .* 2, ...
                    2, 2, 2, ...
                    1, 2, 1];
    
    PrachConfigFR1UnpairedSpectrum.preambleFormat = preambleFormat_column(prachConfig.PrachConfigurationIndex + 1);
    
    switch PrachConfigFR1UnpairedSpectrum.preambleFormat
        case '0'
            PrachConfigFR1UnpairedSpectrum.L_RA = 839;
            PrachConfigFR1UnpairedSpectrum.PrachDuration = 1;
        case '1'
            PrachConfigFR1UnpairedSpectrum.L_RA = 839;
            PrachConfigFR1UnpairedSpectrum.PrachDuration = 2;

        case {'2', '3'}
            PrachConfigFR1UnpairedSpectrum.L_RA = 839;
            PrachConfigFR1UnpairedSpectrum.PrachDuration = 4;
        otherwise
            PrachConfigFR1UnpairedSpectrum.L_RA = 139;
            PrachConfigFR1UnpairedSpectrum.PrachDuration = prachDuration_column(prachConfig.PrachConfigurationIndex + 1);
    end

    L_RA_index = find(L_RA_column == PrachConfigFR1UnpairedSpectrum.L_RA);
    prach_delta_f_index = find(delta_f_RA_forPrach_column == prachConfig.SubcarrierSpacing);
    PUSCH_delta_f_index = find(delta_f_forPUSCH_column == carrierConfig.SubcarrierSpacing);
    k_bar_index = intersect(intersect(L_RA_index, prach_delta_f_index), PUSCH_delta_f_index);
    
    if isempty(k_bar_index)
        warning('No k_bar and N_RA_RB for PRACH_SCS = %d and PUSCH_SCS = %d', prach_delta_f_index, PUSCH_delta_f_index);
    end


    %% Outputs
    PrachConfigFR1UnpairedSpectrum.PrachConfigurationIndex = prachConfig.PrachConfigurationIndex;
    
    PrachConfigFR1UnpairedSpectrum.x = x_column(prachConfig.PrachConfigurationIndex + 1);
    PrachConfigFR1UnpairedSpectrum.y = y_column(prachConfig.PrachConfigurationIndex + 1);
    PrachConfigFR1UnpairedSpectrum.subframeNumber = subframeNumber_column{1, prachConfig.PrachConfigurationIndex + 1};
    PrachConfigFR1UnpairedSpectrum.startingSymbol = startingSymbol_column(prachConfig.PrachConfigurationIndex + 1);
    PrachConfigFR1UnpairedSpectrum.numPrachSlotsWithinASubframe = numberPrachSlotWithinSubframe_column(prachConfig.PrachConfigurationIndex + 1);
    PrachConfigFR1UnpairedSpectrum.numTimeDomainPrachOccasionsWithinAPrachSlot = numberTimeDomainPrachOccasionsWithinPrachSlot_column(prachConfig.PrachConfigurationIndex + 1);
    
    PrachConfigFR1UnpairedSpectrum.N_RA_RB = N_RA_RB_column(k_bar_index);
    PrachConfigFR1UnpairedSpectrum.k_bar = k_bar_column(k_bar_index);
end

