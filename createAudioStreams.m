function createAudioStreams(configuration,source,decorrON)

% HALL:                 Pollack = 1, Tanna = 2
% SOURCE:               High hat = 1, Faure Cello = 2
% LOUDSPEAKER SYSTEM:   See getLoudspeakerConfiguration()
% CONFIGURATION:        See getParameterSet()
% DFM ON/OFF decorrON - turn decorrelation off if you want it to sound bad.

% 0th order is an excellent bad case - sounds awful.
zeroOrder = 0;

% Get the configuration used to design the filters:
config = getParameterSet(configuration);
[ virMicDir, numLdspkr ] = getLoudspeakerConfiguration(config.ldspkrArray);

if(config.numChannels ~= numLdspkr)
    disp('Something is wrong with the configuration specification.');
    return;
end;

% Synthesize filename - used in all files - serves as versioning sanity check:
filename = strcat('dfm',num2str(configuration));

cd('RoomImpulses');

% Load impulses with direct sound clipped off.
if(config.hall == 1)
    approximateDelay    = 2048; % Chop out the propagation time.
    gainCompDB          = 11;   % dB comp of impulses.
    timeNoiseFloor      = 2.6;  % Time in seconds that IR goes to noise
    [ bformatW, Fs ]    = audioread('BfrmtW_NearBrightCenter_pollack.wav');
    [ bformatX, ~ ]     = audioread('BfrmtX_NearBrightCenter_pollack.wav');
    [ bformatY, ~ ]     = audioread('BfrmtY_NearBrightCenter_pollack.wav');
    [ bformatZ, ~ ]     = audioread('BfrmtZ_NearBrightCenter_pollack.wav');
elseif(config.hall == 2)
    approximateDelay    = 874;  % Chop out the propagation time.
    gainCompDB          = 20;   % dB comp of impulses.
    timeNoiseFloor      = 1.6;  % Time in seconds that IR goes to noise
    [ bformatW, Fs ]    = audioread('BfrmtW_NearCenter_tanna.wav');
    [ bformatX, ~ ]     = audioread('BfrmtX_NearCenter_tanna.wav');
    [ bformatY, ~ ]     = audioread('BfrmtY_NearCenter_tanna.wav');
    [ bformatZ, ~ ]     = audioread('BfrmtZ_NearCenter_tanna.wav');
end;
cd ..

gain = 10^(gainCompDB/20);
samplesNoiseFloor = Fs * timeNoiseFloor;

if(0)
    % This may be of some use to some users, but doesn't apply to DFM3:
    indexDiffuse = computeEchoDensity(bformatW,approximateDelay);
else
    indexDiffuse = approximateDelay + (config.numTaps/2);
end;

% Clip the RIR at indexDiffuse as the symmetric filters will provide ramp:
timeCompW = gain .* bformatW(indexDiffuse:(samplesNoiseFloor));
timeCompX = gain .* bformatX(indexDiffuse:(samplesNoiseFloor));
timeCompY = gain .* bformatY(indexDiffuse:(samplesNoiseFloor));
timeCompZ = gain .* bformatZ(indexDiffuse:(samplesNoiseFloor));

figure();
signal = gain .* bformatW(approximateDelay:samplesNoiseFloor); 
plot(20*log10(abs(signal)),'r');
hold on; grid on;
maxDB = 20*log10(abs(max(signal)));
plot([ indexDiffuse indexDiffuse ],[maxDB -120],'o-');
legend({'Original W','Clipped Start Point'})


%% Load source material
cd('SourceMaterial');
if(source == 1)
    [ sourceAudio, FsSrc ] = audioread('looseZildjanHighHat.wav');
elseif(source == 2)
    [ sourceAudio, FsSrc ] = audioread('faureCello.wav');
elseif(source == 3)
    [ sourceAudio, FsSrc ] = audioread('unitImpulse.wav');
end;

if(Fs ~= FsSrc)
    disp('Sample Rate of IR and SRC do not match. USING IR SAMPLE RATE');
end;
cd ..

% Pad the source for fftfilt buffer voodoo.
impulseTail     = zeros(length(timeCompW),1);
srcPadded       = [ sourceAudio; impulseTail ];

% Set up buffers for DFM convolution:
impulseLength   = length(srcPadded);
channel         = zeros(numLdspkr,impulseLength);

% Convolve RIR with source material.
streamW = fftfilt((1/1).*timeCompW,srcPadded);
streamX = fftfilt((1/1).*timeCompX,srcPadded);
streamY = fftfilt((1/1).*timeCompY,srcPadded);
streamZ = fftfilt((1/1).*timeCompZ,srcPadded);

% Load specifed range of virtual array decorrelation filters.
cd ImpulsesAndFilters
if(decorrON)
    load(filename);
end;
cd ..

% DFM uses linear combinations of omni and dipoles to form cardioids microphones
bFormatGain     = zeros(numLdspkr,4);

if(zeroOrder)
    omniGain = 10^(3/20);   % By ear gain bump
    for ii = 1:numLdspkr    % Zap X,Y,Z for 0-order.
        bFormatGain(ii,1)   = omniGain;
        bFormatGain(ii,2)   = 0.00;
        bFormatGain(ii,3)   = 0.00;
        bFormatGain(ii,4)   = 0.00;
    end;
else
    omniGain = 1.0;    
    for ii = 1:numLdspkr      % compute gains for W,X,Y,Z:
        bFormatGain(ii,1)   = omniGain;
        bFormatGain(ii,2)   = cos(virMicDir(ii,1)) * cos(virMicDir(ii,2));
        bFormatGain(ii,3)   = sin(virMicDir(ii,1)) * cos(virMicDir(ii,2));
        bFormatGain(ii,4)   = sin(virMicDir(ii,2));
    end;
end;

for ii = 1:numLdspkr      % apply for W,X,Y,Z:
    channel(ii,:)       = streamW' * bFormatGain(ii,1);
    channel(ii,:)       = channel(ii,:) + streamX' * bFormatGain(ii,2);
    channel(ii,:)       = channel(ii,:) + streamY' * bFormatGain(ii,3);
    channel(ii,:)       = channel(ii,:) + streamZ' * bFormatGain(ii,4);
end;

if(decorrON)
    for ii = 1:numLdspkr
        channel(ii,:)       = fftfilt(decorrFilt(ii,:),channel(ii,:));
    end;
else
    filename = strcat(filename,'OFF');
end;

cd('OutputAudioWAV');
filename2 = strcat(filename,'Hall',int2str(config.hall),'Src',int2str(source),'.wav');
audiowrite(filename2,[srcPadded channel(1:numLdspkr,:)'],Fs);
cd ..



