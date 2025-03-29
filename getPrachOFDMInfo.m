function PrachOFDMInfo = getPrachOFDMInfo(prachConfig, carrierConfig, PrachConfigFR1UnpairedSpectrum)
%GETPRACHOFDMINFO Summary of this function goes here
%   Detailed explanation goes here
    
    nfft = 4096;
    prachNumerology = log2(prachConfig.SubcarrierSpacing / 15);
    kappa = 64;
    T_C = 0.509 * 10^(-6); % in millisecond
    default_samplingRate = 15000 * 2048;
    samplingRate = (prachConfig.SubcarrierSpacing * 1000)  * nfft;
    samplingRate_ratio = samplingRate / default_samplingRate;

    switch PrachConfigFR1UnpairedSpectrum.L_RA
        case 839
            samplingRate_ratio = (carrierConfig.SubcarrierSpacing / 15) * 2;
            format_column = ["0", "1", "2", "3"];
            N_u_column = [24576, 2 * 24576, 4 * 24576, 4 * 6144];
            N_RA_CP_column = [3168, 21024, 4688, 3168];
            
            % Table 2: Chakrapani, Arvind. (2020). On the design details of
            %          SS/PBCH, Signal Generation and PRACH in 5G-NR. 
            %          10.36227/techrxiv.12465743.v2. 
            N_RA_GP_column = [2976, 21904, 4528, 2976];
            pathProfile_column = [192, 512, 192, 192];

        case {139, 1151, 571}
            format_column = ["A1", "A2", "A3", "B1", "B2", "B3", "B4", "C0", "C2"];
            N_u_column = [2 * 2048 * 2^(-prachNumerology),  ...
                          4 * 2048 * 2^(-prachNumerology),  ...
                          6 * 2048 * 2^(-prachNumerology),  ...
                          2 * 2048 * 2^(-prachNumerology),  ...
                          4 * 2048 * 2^(-prachNumerology),  ...
                          6 * 2048 * 2^(-prachNumerology),  ...
                          12 * 2048 * 2^(-prachNumerology), ...
                          2048 * 2^(-prachNumerology),      ...
                          4 * 2048 * 2^(-prachNumerology)];
            N_RA_CP_column = [288 * 2^(-prachNumerology),  ...
                              576 * 2^(-prachNumerology),  ...
                              864 * 2^(-prachNumerology),  ...
                              216 * 2^(-prachNumerology),  ...
                              360 * 2^(-prachNumerology),  ...
                              504 * 2^(-prachNumerology),  ...
                              936 * 2^(-prachNumerology),  ...
                              1240 * 2^(-prachNumerology), ...
                              2048 * 2^(-prachNumerology)];

            % Table 2: Chakrapani, Arvind. (2020). On the design details of
            %          SS/PBCH, Signal Generation and PRACH in 5G-NR. 
            %          10.36227/techrxiv.12465743.v2. 
            N_RA_GP_column = [0, ...
                              0, ...
                              0, ...
                              72 * 2^(-prachNumerology),   ...
                              216 * 2^(-prachNumerology),  ...
                              360 * 2^(-prachNumerology),  ...
                              792 * 2^(-prachNumerology),  ...
                              1096 * 2^(-prachNumerology), ...
                              2916 * 2^(-prachNumerology)];
            pathProfile_column = [96 * 2^(-prachNumerology),  ...
                                  114 * 2^(-prachNumerology),  ...
                                  114 * 2^(-prachNumerology),  ...
                                  96 * 2^(-prachNumerology),  ...
                                  114 * 2^(-prachNumerology),  ...
                                  114 * 2^(-prachNumerology),  ...
                                  114 * 2^(-prachNumerology),  ...
                                  114 * 2^(-prachNumerology), ...
                                  114 * 2^(-prachNumerology)];
        otherwise
            warning('"prachConfig.L_RA" must be one of {839, 139, 1151, 571}}');
            return;
    end
    % Use split if format is something like "A1/B1", etc
    preamble_format = split(PrachConfigFR1UnpairedSpectrum.preambleFormat, '/');
    preamble_format = preamble_format(1);
    
    format_index = find(format_column == preamble_format);
    N_u = N_u_column(format_index) * samplingRate_ratio;
    N_RA_CP = N_RA_CP_column(format_index) * samplingRate_ratio;
    N_RA_GP = N_RA_GP_column(format_index) * samplingRate_ratio;
    pathProfile = pathProfile_column(format_index) * samplingRate_ratio;
    
    % Since delta_f_RA = 30 and format B4 cross 0.5ms => n = 1
    % Clause 5.3.2, in TS 38.211
    switch prachConfig.SubcarrierSpacing 
        case {15, 30, 60, 120}
            time_interval = (N_u + N_RA_CP) * kappa * T_C;

            if time_interval > 0.5
                n = 1;
            else 
                n = 0;
            end
        otherwise
            n = 0;
    end

    %
    N_RA_CP_l = N_RA_CP + n * 16 * samplingRate_ratio;

    % OFDMInfo = getOFDMInfo(carrierConfig.n_UL_RB, prachConfig.SubcarrierSpacing);
    % % Take 15kHz as reference subcarrier frequency
    % delta_f_ref = 15e3;
    % N_f_ref = 2048;
    % Ts = 1 / (delta_f_ref * N_f_ref);
    % 
    % % Reference Sample Rate
    % sampleRate_ref = 1 / Ts;
    % 
    % % for L_RA = 139, table 6.3.3.1-2, in TS 38.211
    % numSamplerPerPreamble = 2048 * 2^-OFDMInfo.mui;
    % 
    % % OFDM symbol sample, clause 5.3.1, in 38.211
    % N_mui_u = 2048 * 2^-OFDMInfo.mui;
    % 
    % % Get number of sample for Prach Cyclic Prefix and Prach Symbof, format B4, 
    % % table 6.3.3.1-2, in TS 38.211 (write function)
    % N_RA_CP = 936 * 2^-OFDMInfo.mui;
    % N_u = 12 * 2048 * 2^-OFDMInfo.mui;
    % 
    % 
    % N_RA_CP_l = N_RA_CP + n * 16;
    % 
    % % Get number of sample for Prach Guard Period, format B4, in < PRACH
    % % Cell Dimensioning >,
    % % in 'https://www.sharetechnote.com/html/5G/5G_RACH.html'
    % N_RA_GP = 792 * 2^-OFDMInfo.mui;
    % 
    % %
    % remainSample = 144;
    % 
    % % Adjust Cyclic Prefix, Guard Period, and remain from 
    % % reference size_ifft = 2048, reference delta_f_RA = 15
    % % for new size_ifft = 4096 and delta_f_ref = 30
    % N_RA_CP_l = N_RA_CP_l * OFDMInfo.OFDMSampleRate / sampleRate_ref;
    % N_RA_GP = N_RA_GP * OFDMInfo.OFDMSampleRate / sampleRate_ref;
    % remainSample = remainSample * OFDMInfo.OFDMSampleRate / sampleRate_ref;
    
    PrachOFDMInfo.Sequence = N_u;
    PrachOFDMInfo.CyclicPrefix = N_RA_CP_l;
    PrachOFDMInfo.GuardPeriod = N_RA_GP;
    PrachOFDMInfo.PathProfile = pathProfile;

end

