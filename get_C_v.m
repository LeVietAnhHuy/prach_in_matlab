function [C_v, C_v_arr] = get_C_v(prachConfig, PrachConfigFR1UnpairedSpectrum, N_CS)
%GETC_V Summary of this function goes here
%   Detailed explanation goes here
    
    switch prachConfig.Set
        case 'Unrestricted'
            if N_CS == 0
                C_v = 0;
                C_v_arr = NaN;
                return;
            else
                v = 0:(floor(PrachConfigFR1UnpairedSpectrum.L_RA / N_CS) - 1);
                C_v_arr = v * N_CS;
            end

        case 'RestrictedTypeA'
            disp('do later')
        otherwise
            warning('set must be "Unrestricted" or "RestrictedTypeA"');
            return
    end

    C_v_preamble_arr = zeros(1, prachConfig.PreambleIndex + 1);
    
    C_v_arr_index = 1;
    for i = 1:length(C_v_preamble_arr)
        if C_v_arr_index > length(C_v_arr)
            C_v_arr_index = 1;
        end
        C_v_preamble_arr(i) = C_v_arr(C_v_arr_index);
        C_v_arr_index = C_v_arr_index + 1;
    end
    
    C_v = C_v_preamble_arr(prachConfig.PreambleIndex + 1);
end

