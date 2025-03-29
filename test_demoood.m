close all; clear; clc;
prachConfig.PrachConfigurationIndex = 158;
prachConfig.RootSequenceIndex = 39;
prachConfig.L_RA = 139;

prachConfig.PreambleIndex = 60;

prachConfig.SubcarrierSpacing = 30;  % delta_f_RA
prachConfig.Set = 'Unrestricted';
prachConfig.ZeroCorrelationZoneConfig = 8;
prachConfig.FrequencyRange = 'FR1';
prachConfig.SpectrumType = 'Unpaired';

% n_RA_start or nPrachFreqStart
prachConfig.PrachFreqStart = 0;

%% Carrier Config
% N_grid_size_mui, clause 5.3.2, TS 38.211
carrierConfig.n_UL_RB = 273;

carrierConfig.SubcarrierSpacing = 30;
carrierConfig.numElementPerResourceBlock = 12;
carrierConfig.numFrame = 2;

% get N_CS
N_CS = get_N_CS(prachConfig);

% Physical root sequence number (u)
[u, u_arr] = get_u(prachConfig, N_CS);

% Cyclic shift
[C_v, C_v_arr] = get_C_v(prachConfig, N_CS);

x_u = zadoffChu(u, prachConfig.L_RA);
x_uv = circshift(x_u, [0 -C_v]);
y_uv = fft(x_uv);

C_v_index_used = 0;
maxPA = 0;
% figure(1)
for C_v_index =  1:length(C_v_arr)
    x_uv_current = circshift(x_u, [0 -C_v_arr(C_v_index)]);
    corr = x_uv .* conj(x_uv_current);
    ifftSeq = zeros(1024, 1);
    ifftSeq(1:139) = corr;
    ifftSeq = abs(ifft(ifftSeq));
    ifftSeq_ifftShift = fftshift(ifftSeq);
    
    figure(C_v_index)
    plot(ifftSeq_ifftShift);
    hold on 
    plot(ifftSeq);

    title(C_v_arr(C_v_index));
    legend('yes','no')
    % hold on

    currentPA = max(ifftSeq) / mean(ifftSeq);

    peakThreshold = max(ifftSeq) / 2;
    windowLength = ceil(1024 / length(C_v_arr));
    timingAdvance = 0;
    if (max(ifftSeq(1: windowLength)) > peakThreshold) && (max(ifftSeq((end - windowLength): end)) > peakThreshold)
        C_v_index_used = C_v_index;

        ifftSeq_ifftShift = ifftshift(ifftSeq);
        peakPosition = find(ifftSeq_ifftShift == max(ifftSeq_ifftShift));
        

        while peakPosition < 1024
            peakPosition = peakPosition + windowLength;
        end

        timingAdvance = (peakPosition - 1024) / 8;

        break    
    end
    % if(currentPA > maxPA)
    %     maxPA = currentPA;
    %     C_v_index_used = C_v_index;
    %     maxPA_ifftSeq = fftshift(ifftSeq);
    % end
end

% figure(2)
% for C_v_index = 1:length(C_v_arr)
%     x_uv_current = circshift(x_u, [0 -C_v_arr(C_v_index)]);
%     corr = x_uv .* conj(x_uv_current);
%     ifftSeq = zeros(1024, 1);
%     ifftSeq(1:139) = corr;
%     ifftSeq = abs(ifft(ifftSeq));
%     ifftSeq_ifftShift = ifftSeq;
% 
%     plot(ifftSeq_ifftShift);
%     hold on
%     currentPA = max(ifftSeq) / mean(ifftSeq);
%     if(currentPA > maxPA)
%         maxPA = currentPA;
%         C_v_index_used = C_v_index;
%         maxPA_ifftSeq = fftshift(ifftSeq);
%     end
% end