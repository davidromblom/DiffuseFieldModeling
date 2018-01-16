function config = getParameterSet(configuration)

% DFM parameters:
config.numChannels   = 0;    % number of equally spaced channels.
config.arrayDistance = 0;    % distance of the loudspeaker (and virtual cardioid) array.
config.filterRNG     = 1;    % RNG seed for filters.
config.fieldRNG      = 1;    % RNG seed for field simulation.


% Parameters that do not change:
config.Fs            = 48000;
config.sos           = 343;
config.radialSamples = 64;

switch configuration
    case 01, % Faux binaural
        config.controlFreq      = [0   10  20  50  100 200 500 1000 2000 5000 10000 20000 24000];
        config.controlRT60      = [1.8 1.8 1.8 1.8 1.8 1.8 1.8 1.8  1.7  1.6  0.8   0.4   0.4 ];
        config.hall             = 1; % Pollack
        config.numChannels      = 2;
        config.ldspkrArray      = 4;
        config.arrayDistance    = 0.10;
        config.numTaps          = 32768;
        config.truncationGain   = 0;
        config.filterRNG        = 7;
        config.isBinaural       = 1;
        config.roomScale        = 2.0;
    case 02, % Faux binaural
        config.controlFreq      = [0   10  20  50  100 200 500 1000 2000 5000 10000 20000 24000];
        config.controlRT60      = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0  1.0  0.9  0.8   0.4   0.4 ];
        config.hall             = 2; % Tanna
        config.numChannels      = 2;
        config.ldspkrArray      = 4;
        config.arrayDistance    = 0.10;
        config.numTaps          = 16384;
        config.truncationGain   = 0;
        config.filterRNG        = 7;
        config.isBinaural       = 1;
        config.roomScale        = 1.6;
    case 21, % 16 channel horizontal - Pollack
        config.controlFreq      = [0   10  20  50  100 200 500 1000 2000 5000 10000 20000 24000];
        config.controlRT60      = [1.8 1.8 1.8 1.8 1.8 1.8 1.8 1.8  1.7  1.6  0.8   0.4   0.4 ];
        config.hall             = 2; % Pollack
        config.numChannels      = 16;
        config.ldspkrArray      = 2;
        config.arrayDistance    = 1.5;
        config.numTaps          = 32768;
        config.truncationGain   = 0;
        config.filterRNG        = 7;
        config.isBinaural       = 0;
        config.roomScale        = 2.0;
    case 22, % 16 channel horizontal - Tanna
        config.controlFreq      = [0   10  20  50  100 200 500 1000 2000 5000 10000 20000 24000];
        config.controlRT60      = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0  1.0  0.9  0.8   0.4   0.4 ];
        config.hall             = 2; % Tanna
        config.numChannels      = 16;
        config.ldspkrArray      = 2;
        config.arrayDistance    = 1.5;
        config.numTaps          = 32768;
        config.truncationGain   = 0;
        config.filterRNG        = 7;
        config.isBinaural       = 0;
        config.roomScale        = 1.6;
    otherwise
        disp('Nope Smarty Pants')
        
end;

