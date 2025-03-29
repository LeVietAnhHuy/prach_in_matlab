aveconfig = [];

waveconfig.NumSubframes = 10*2; % Number of 1 ms subframes in generated waveform

waveconfig.DisplayGrids = 1;          % Display the resource grids and the generated waveforms

 

% Define a carrier configuration object

carrier = nrCarrierConfig;

carrier.SubcarrierSpacing = 120;

carrier.NSizeGrid = 66;

 

% Store the carrier into the waveconfig structure

waveconfig.Carriers = carrier;

 

% PRACH configuration

prach = nrPRACHConfig;

prach.FrequencyRange = 'FR2';   % Frequency range ('FR1', 'FR2')

prach.DuplexMode = 'TDD';         % Duplex mode ('FDD', 'TDD', 'SUL')

prach.ConfigurationIndex = 70;    % Configuration index (0...255)

prach.SubcarrierSpacing = 120;   % Subcarrier spacing (1.25, 5, 15, 30, 60, 120)

prach.FrequencyIndex = 1;         % Index of the PRACH transmission occasions in frequency domain (0...7)

prach.TimeIndex = 0;                 % Index of the PRACH transmission occasions in time domain (0...6)

prach.ActivePRACHSlot = 1;        % N^RA_slot. 38.211 - 5.3.2 (0, 1)

 

% Store the PRACH configuration and additional parameters in the

% waveconfig structure

waveconfig.PRACH.Config = prach;

waveconfig.PRACH.AllocatedPreambles = 'all';   % Index of the allocated PRACH preambles

waveconfig.PRACH.Power = 0;                       % PRACH power scaling in dB

 

% NOTE : hNRPRACHWaveformGeneratorFR2 is modified version of the original helper function :

%            hNRPRACHWaveformGenerator

[waveform,gridset,winfo] = hNRPRACHWaveformGenerator(waveconfig);

 

hFig=figure();

set(hFig, 'Position', [100 500 800 250]);

set(gcf,'color','w')

    

subplot(2,6,[1 2 7 8]);

plot(real(winfo.WaveformResources.PRACH.Resources(1).PRACHSymbols), ...
      imag(winfo.WaveformResources.PRACH.Resources(1).PRACHSymbols),'ro')

daspect([1 1 1]);

axis([-15 15 -15 15]);

str = sprintf("FR2(120Khz),ConfigIndex = %d",prach.ConfigurationIndex);

title(str);

 

sLength = length(winfo.WaveformResources.PRACH.Resources(1).PRACHSymbols);

 

subplot(2,6,[3 4 5 6]);

plot(real(winfo.WaveformResources.PRACH.Resources(1).PRACHSymbols),'r-')

ylim([-15 15]);

xlim([0 sLength]);

str = sprintf("Real[PRACH] : sequence length = %d",sLength);

title(str);

 

subplot(2,6,[9 10 11 12]);

plot(imag(winfo.WaveformResources.PRACH.Resources(1).PRACHSymbols),'r-')

ylim([-15 15]);

xlim([0 sLength]);

str = sprintf("Real[PRACH] : sequence length = %d",sLength);

title(str);

 

 

winfo.WaveformResources.PRACH.Resources(1).PRACHIndicesInfo

 

sCount = length(winfo.WaveformResources.PRACH.Resources);

 

sList = [];

for i = 1:sCount

   sList = [sList winfo.WaveformResources.PRACH.Resources(i).NPRACHSlot];  

end    

 

sList