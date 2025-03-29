function OFDMInfo = getOFDMInfo(numResourceBlock_per_ResourceGrid, delta_f_RA)
%GETOFDMINFO Summary of this function goes here
%   Detailed explanation goes here
    
    mui = log2(delta_f_RA / 15);
    OFDMSymbolPerSlot = 14; % normal Cyclic Prefix
    numResourceElement_per_ResourceGrid_freq = numResourceBlock_per_ResourceGrid * 12;
    % Choose IDFT size based on three rules:
    %   * power of 2 size
    %   * maximum occupancy (K / nfft) of 85%
    %   * minimum of 128 (so that cyclic prefix lengths are
    %     always integer)
    % size_idft = max(power(2,ceil(log2(numResourceElement_per_ResourceGrid_freq / 0.85))),128);
    
    size_ifft = 4096;
    OFDMSampleRate = delta_f_RA * 1e3 * size_ifft;
    slotPerSubframe = 2^mui;
    
    % Normal Cyclic Prefix length for l ~= 0 and l ~= 7 * 2^mui, in clause
    % 5.3.1
    % N_mui_CP_l = 144 * 2^-mui;
    
    % Normal Cyclic Prefix length for l = 0 and l = 7 * 2^mui, in clause
    % 5.3.1
    % N_mui_CP_l1 = 144 * 2^-mui + 16;

    % Adjust for new IDFT size = 4096 since standard IDFT size = 2048
    % N_mui_CP_l = N_mui_CP_l * 2048 / size_ifft;
    
    % Table F.5.3-2, in TS 38.101-1 V17.5.0
    N_mui_CP_l = 288;
    
    % Normal Cyclic Prefix length for l = 0 and l = 7 * 2^mui, in clause
    % H.2.4, table H.2.4-1, in TP to TS 38.176-1: Annex G and H: In-channel TX test
    N_mui_CP_l1 = N_mui_CP_l + size_ifft / 64;
    
    OFDMInfo.mui = mui;
    OFDMInfo.OFDMSymbolPerSlot = OFDMSymbolPerSlot;
    OFDMInfo.numResourceElement_per_ResourceGrid_freq = numResourceElement_per_ResourceGrid_freq;
    OFDMInfo.size_ifft = size_ifft;
    OFDMInfo.OFDMSampleRate = OFDMSampleRate;
    OFDMInfo.slotPerSubframe = slotPerSubframe;
    OFDMInfo.firstCyclicPrefix = N_mui_CP_l;
    OFDMInfo.subsequentCylicPrefix = N_mui_CP_l1;
end

