function synthesizeBinaural(configuration,source,decorrON)

debugRIR = 0;
debugNoggin = 0;

% Get the configuration used to design the filters:
config = getParameterSet(configuration);
[ ~, numLdspkr ] = getLoudspeakerConfiguration(config.ldspkrArray);

if(config.numChannels ~= numLdspkr)
    disp('Something is wrong with the configuration specification.');
    return;
end;

% Synthesize filename - used in all files - serves as versioning sanity check:
filename = strcat('dfm',num2str(configuration));
if(~decorrON)
    filename = strcat(filename,'OFF');
end;

% Copy over the frequently used stuff:
Fs              = config.Fs;
speedOfSound    = config.sos;
numChannels     = config.numChannels;
nogginRadius    = config.arrayDistance;

cd('OutputAudioWAV');
filename1 = strcat(filename,'Hall',int2str(config.hall),'Src',int2str(source),'.wav');
audio           = audioread(filename1);
directSound     = audio(:,1);
diffuseSound    = audio(:,(2:numChannels+1));
cd ..

audioLength     = length(directSound);
% Add a small amount for ITD and a huge amount for reflection propagation.
audioLengthPadded = audioLength+16384;

% All audio files need to be the same length or MATLAB whines
diffuseOutput   = zeros(2,audioLengthPadded);

if(config.isBinaural == 1)
    diffuseOutput(:,1:length(diffuseSound))  = diffuseSound';
    % Simplistic radiation power model: -8 dB DC, pole at 0.9700
    b = [ 0.9897 -0.9803 ]; a = [ 1 -0.9700 ];
    diffuseOutput = filter(b,a,diffuseOutput);
else
    disp('Diffuse Channels are specified as loudspeaker channels.');
    return;
end;

% Synthesize Direct:
directAzimuth = 0;
% Pad for ITD calculations
binauralDirect  = zeros(2,audioLengthPadded);

[ lFilt, LEFTtimeDelay ]   = computeSHM(directAzimuth+110,Fs,speedOfSound,nogginRadius);
[ rFilt, RIGHTtimeDelay ] = computeSHM(directAzimuth-110,Fs,speedOfSound,nogginRadius);

LEFTIndexes     = LEFTtimeDelay:(LEFTtimeDelay+audioLength-1);
RIGHTIndexes    = RIGHTtimeDelay:(RIGHTtimeDelay+audioLength-1);

binauralDirect(1,LEFTIndexes) = filter([lFilt.b0 lFilt.b1],[lFilt.a0 lFilt.a1],directSound(:,1));
binauralDirect(2,RIGHTIndexes) = filter([rFilt.b0 rFilt.b1],[rFilt.a0 rFilt.a1],directSound(:,1));

if(debugNoggin)
    [H,W] = freqz([lFilt.b0 lFilt.b1],[lFilt.a0 lFilt.a1],4096);
    
    freqVec = (Fs / (2*pi)) .* W;
    subplot(2,1,1);
    semilogx(freqVec,20*log10(abs(H)));
    grid on;
    subplot(2,1,2);
    semilogx(freqVec,angle(H));
    grid on;
end;

% Simplistic radiation power model: -8 dB DC, pole at 0.9700
b = [ 0.9897 -0.9803 ]; a = [ 1 -0.9700 ];
directSoundHPF = filter(b,a,directSound(:,1));
% Also a lowpass filter for the reflections
b = [ 0.5511 -0.3511 ]; a = [ 1 -0.8000 ];
directSoundBPF = filter(b,a,directSoundHPF);

% Synthesize Reflections - there is no rigorous physical meaning to these:
bounceAzimuth  = [ 32, -37, 21, -27,  11, -13, 71, -83, 12, -20, 67, -71, 17, -23, 11, -37 ];
roomScale = config.roomScale;
bounceDelaySmpls    = roomScale .* [ 907, 1051, 1171, 1279, 1423, 1493, 1607, 1741, 1933, 2087, 2473, 2797, 2971, 3257, 3617, 4241];
bounceDelaySmpls    = floor(bounceDelaySmpls);
bounceDistanceMeters = (343/Fs) .* bounceDelaySmpls;
bounceGain = 1 ./ bounceDistanceMeters;
reflectionOutput    = zeros(2,audioLengthPadded);
binauralReflection  = zeros(2,audioLengthPadded);

for i = 1:16
    [ lFilt, LeftITD ]    = computeSHM(bounceAzimuth(i)+110,Fs,speedOfSound,nogginRadius);
    [ rFilt, RightITD ]   = computeSHM(bounceAzimuth(i)-110,Fs,speedOfSound,nogginRadius);
    
    LEFTtimeDelay   = bounceDelaySmpls(i) + LeftITD;
    RIGHTtimeDelay  = bounceDelaySmpls(i) + RightITD;
    
    LEFTIndexes     = LEFTtimeDelay:(LEFTtimeDelay+audioLength-1);
    RIGHTIndexes    = RIGHTtimeDelay:(RIGHTtimeDelay+audioLength-1);
    
    binauralReflection(1,LEFTIndexes) = bounceGain(i) .* filter([lFilt.b0 lFilt.b1],[lFilt.a0 lFilt.a1],directSoundBPF(:,1));
    binauralReflection(2,RIGHTIndexes) = bounceGain(i) .* filter([rFilt.b0 rFilt.b1],[rFilt.a0 rFilt.a1],directSoundBPF(:,1));
    
    reflectionOutput  = reflectionOutput + binauralReflection;
end;

filename2 = strcat(filename,'Hall',int2str(config.hall),'Src',int2str(source));

cd('BinauralAudio');
filename3 = strcat(filename2,'Spat','DIRECT.wav');
audiowrite(filename3,binauralDirect',Fs);
filename3 = strcat(filename2,'Spat','REFLECT.wav');
audiowrite(filename3,reflectionOutput',Fs);
filename3 = strcat(filename2,'Spat','DIFFUSE.wav');
audiowrite(filename3,diffuseOutput',Fs);
cd ..

if(debugRIR)
    figure();
    time = (1:length(diffuseOutput))./48000;
    plot(time,20*log10(abs(diffuseOutput(1,:))))
    hold on; grid on;
    directLimited = max(20*log10(abs(binauralDirect(1,:))),-90);
    plot(time,directLimited,'Black')
    reflectLimited = max(20*log10(abs(reflectionOutput(1,:))),-90);
    plot(time,reflectLimited,'Green')
end;

return;

%% Spherical Head Model based on Brown and Duda's approximation
function [filt, sampleOffset] = computeSHM(earAngleDegrees,Fs,speedOfSound,nogginRadius)

earAngle = 2*pi*earAngleDegrees/360;
% Pre-warp frequency before the Bilinear transform.
% This value follows Oppenheim and Schaffer as EQ Cookbook had a poor fit.
warp = 2/(1/Fs);
% Brown + Duda 98:
alphaMin = 0.1;
thetaConst = 6/5;
tau =  2*nogginRadius/speedOfSound;

% Compute alpha
alpha   = (1+alphaMin/2) + (1-alphaMin/2)*cos(earAngle*thetaConst);

% Create the spherical head shadowing filter:
filt.b0 = (1+alpha*tau*warp)/(1+tau*warp);
filt.b1 = (1-alpha*tau*warp)/(1+tau*warp);
filt.a1 = (1-tau*warp)/(1+tau*warp);
filt.a0 = 1;

% Compute time offset:
if(earAngle < pi/2)
    timeOffset = -nogginRadius/speedOfSound * cos(abs(earAngle));
else
    timeOffset = nogginRadius/speedOfSound * (abs(earAngle) - (pi/2));
end;
% Add a fixed 1ms to account for negative ITD.
timeOffset = timeOffset + 0.001;
sampleOffset = round(timeOffset * Fs); % seconds * (samples/second) = samples
return;



