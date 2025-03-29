function N_CS = get_N_CS(prachConfig, PrachConfigFR1UnpairedSpectrum)
%GETN_CS Summary of this function goes here
%   Detailed explanation goes here

    if (prachConfig.ZeroCorrelationZoneConfig > 15) || (prachConfig.ZeroCorrelationZoneConfig < 0)
        warning('zeroCorrelationZoneConfig must be in range [0, 15]');
        return;
    end
    
    switch prachConfig.SubcarrierSpacing
        case 1.25
            switch prachConfig.Set
                case 'Unrestricted'
                    Unrestricted = [0; 13; 15; 18; 22; 26; 32; 38; 46; ...
                                    59; 76; 93; 119; 167; 279; 419];
                    N_CS = Unrestricted(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                case 'RestrictedTypeA'
                    Restricted_typeA = [15; 118; 22; 26; 32; 38; 46; 55; ...
                                        68; 82; 100; 128; 158; 202; 237; NaN];
                    N_CS = Restricted_typeA(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                case 'RestrictedTypeB'
                    Restricted_typeB = [15; 18; 22; 26; 32; 38; 46; 55; ...
                                        68; 82; 100; 118; 137; NaN; NaN; NaN];
                    N_CS = Restricted_typeB(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                otherwise
                    warning(['Type of set must be "Unrestricted", ' ...
                             '"RestrictedTypeA", or "RestrictedTypeB"']);
                    return;
            end

        case 5
            switch prachConfig.Set
                case 'Unrestricted'
                    Unrestricted = [0; 13; 26; 33; 38; 41; 49; 55; 64; ...
                                    76; 93; 119; 139; 209; 279; 419];
                    N_CS = Unrestricted(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                case 'RestrictedTypeA'
                    Restricted_typeA = [36; 57; 72; 81; 89; 94; 103; 112; ...
                                        121; 132; 137; 152; 173; 195; 216; 237];
                    N_CS = Restricted_typeA(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                case 'RestrictedTypeB'
                    Restricted_typeB = [36; 57; 60; 63; 65; 68; 71; 77; ...
                                        81; 85; 97; 109; 122; 137; NaN; NaN];
                    N_CS = Restricted_typeB(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                otherwise
                    warning(['Type of set must be "Unrestricted", ' ...
                             '"RestrictedTypeA", or "RestrictedTypeB"']);
                    return;
            end

        case {15, 30, 60, 120}
            switch PrachConfigFR1UnpairedSpectrum.L_RA
                case 139
                    Unrestricted = [0; 2; 4; 6; 8; 10; 12; 13; 15; ...
                                    17; 19; 23; 27; 34; 46; 69];
                    N_CS = Unrestricted(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                case 571
                    Restricted_typeA = [0; 8; 10; 12; 15; 17; 21; 25; ...
                                        31; 40; 51; 63; 81; 114; 190; 285];
                    N_CS = Restricted_typeA(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                case 1151
                    Restricted_typeB = [0; 17; 21; 25; 30; 35; 44; 52; ...
                                        63; 82; 104; 127; 164; 230; 383; 575];
                    N_CS = Restricted_typeB(prachConfig.ZeroCorrelationZoneConfig + 1);
                    return;

                otherwise
                    warning(['Type of set must be "Unrestricted", ' ...
                             '"RestrictedTypeA", or "RestrictedTypeB"']);
                    return;
            end

        otherwise
            warning('delta_f_RA must be one of {1.25, 5, 15, 30, 60, 120}');
            return;
    end
end

