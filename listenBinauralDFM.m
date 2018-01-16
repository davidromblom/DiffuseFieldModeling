function listenBinauralDFM(configuration,source,decorrON,directGainDB,reflectionGainDB,roomGainDB)

% Constants
LEFT    = 1;
RIGHT   = 2;

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

% Copy frequently used stuff into local variables for readibility.
Fs              = config.Fs;
speedofSound    = config.sos;
numTaps         = config.numTaps;
numChannels     = config.numChannels;

cd('BinauralAudio');

filename2 = strcat(filename,'Hall',int2str(config.hall),'Src',int2str(source));

filename3 = strcat(filename2,'Spat','DIRECT.wav');
audio       = audioread(filename3);
directLeft  = audio(:,LEFT);
directRight = audio(:,RIGHT);

filename3 = strcat(filename2,'Spat','REFLECT.wav');
audio       = audioread(filename3);
reflectLeft    = audio(:,LEFT);
reflectRight   = audio(:,RIGHT);

filename3 = strcat(filename2,'Spat','DIFFUSE.wav');
audio       = audioread(filename3);
roomLeft    = audio(:,LEFT);
roomRight   = audio(:,RIGHT);

cd ..

directGainLinear = 10^(directGainDB/20);
reflectGainLinear = 10^(reflectionGainDB/20);
roomGainLinear = 10^(roomGainDB/20);

output = zeros(length(audio),2);

output(:,1) = (roomGainLinear .* roomLeft) + (reflectGainLinear .* reflectLeft) + (directGainLinear .* directLeft);
output(:,2) = (roomGainLinear .* roomRight) + (reflectGainLinear .* reflectRight) + (directGainLinear .* directRight);

p = audioplayer(output, config.Fs); 
timeToPause = length(roomLeft) / Fs;
play(p); pause(timeToPause);

return;

