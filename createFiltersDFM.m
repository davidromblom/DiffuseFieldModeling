function createFiltersDFM(configuration,externalRNG)

% Setup and constants
global  jjj;    jjj     = sqrt(-1);
global  c;      c       = 343;
global rho;     rho     = 1.2;

% Helps with reading the code:
X = 1; Y = 2; Z = 3; numDim = 3;
% makes zeros and ones declarations easier to read.
flat = 1;

config = getParameterSet(configuration);
[ virMicDir, numLdspkr ] = getLoudspeakerConfiguration(config.ldspkrArray);

if(config.numChannels ~= numLdspkr)
    disp('Something is wrong with the configuration specification.');
    return;
end;

% Copy frequently used stuff into local variables for readibility.
Fs              = config.Fs;
%speedofSound    = config.sos;
numTaps         = config.numTaps;
numChannels     = config.numChannels;

nyquistIndex    = ((numTaps/2)+1);
%bins2plot       = nyquistIndex/1;
nDFT            = numTaps;
freqRes         = Fs/nDFT;                          % (1/T) = cycles/sec per bin.
freqVec         = freqRes * (0:nyquistIndex-1);     % 0.0 = DC (0Hz) to Nyquist.

decorrFilt     = zeros(numChannels,numTaps);

filename = strcat('dfm',num2str(configuration));

validateSpacing = 0;
timeAlign = 1;

if(externalRNG == 0)
    rng(config.filterRNG);
else
    rng(externalRNG);
end;

%% Frequency spacing
controlFrequencies  = config.controlFreq;
controlRT6O         = config.controlRT60;

% Assume that freq spacing will not be less than 1 bin.
maxNumDesignFreqs = numTaps;
designFreqL = zeros(maxNumDesignFreqs,flat);
timeEnvL = zeros(maxNumDesignFreqs,numTaps);
tempBuffer = zeros(1,numTaps);
designFreqCurrentHz = 16; % Hz
longestWindow = 0;
earliestStartIndex = numTaps;
maxDesignFreq = 8000;

meanDiffDB = 10^(2.5/10);   % For exponential distribution in p^2, the dB mean drops 2.5 dB.

% Define an autocorrelation level approximately the width of the main lobe.
edgeDBlevelFAC = -20;
edgeLevelFAC = 10^(edgeDBlevelFAC/20);
filterCenter = numTaps / 2;
numDesignFreq = 0;
while(designFreqCurrentHz < maxDesignFreq);
    numDesignFreq = numDesignFreq + 1;
    RT60 = interp1(controlFrequencies,controlRT6O,designFreqCurrentHz);
    % Relationship is 13.82 for squared quantities, 6.9 for linear quantities.
    % The FAC equation below is from AP who uses 13.82 in his derivations.
    tau = RT60/13.82;
    
    % Estimate the bandwidth of the FAC function. (in radians / seconds.)
    bandwidthFACrad = (1./tau) .* sqrt((1/edgeLevelFAC) - 1);
    % Bandwidth is two sided and defines Hamming window.
    % Normalized freq (Oppenheim + Schafer) [ 0 pi ] = [ DC Nyquist ]
    bdwtHammDscrtRad = (2 * bandwidthFACrad) / Fs;
    
    % Main lobe width (norm freq) of a Hamm window is 8pi / M (edge to edge)
    lengthHamm = (8 * pi) / bdwtHammDscrtRad;
    % Type I FIR filter - symmetric around an integer. M = even = L + 1.  window = hamming(L)
    lengthHamm = 1 + 2*floor(lengthHamm/2);
    longestWindow = max(longestWindow,lengthHamm);
     
    if(timeAlign)
        % Symmetry point of window:
        alignmentPoint = (lengthHamm + 1) / 2;
        timeEnv = hamming(lengthHamm);
        %timeEnv = hann(lengthHamm);
        shiftStart = filterCenter - alignmentPoint + 1; % + 1 for symmetry point at numTaps/2
        shiftEnd = shiftStart + length(timeEnv) - 1;    % - 1 for inclusive index
        tempBuffer(flat,shiftStart:shiftEnd) = timeEnv';
        earliestStartIndex = min(shiftStart,earliestStartIndex);
    else
        timeEnv = hamming(lengthHamm);
        % timeEnv = hann(lengthHamm);
        tempBuffer(flat,1:length(timeEnv)) = timeEnv';
    end;

    % Overlap of filters determined empirically:
    %deltaFreqCTrad = (16/16) * (Fs * bdwtHammDscrtRad / 2);
    deltaFreqCTrad = (14/16) * (Fs * bdwtHammDscrtRad / 2);
    deltaFreqHz = deltaFreqCTrad / (2*pi);

    % Store for filter design below
    designFreqL(numDesignFreq,flat) = designFreqCurrentHz;
    timeEnvL(numDesignFreq,1:length(tempBuffer)) = tempBuffer(flat,:);
    
    designFreqCurrentHz = designFreqCurrentHz + deltaFreqHz;
end;

expMag     = zeros(numChannels,numDesignFreq);
expPhase   = zeros(numChannels,numDesignFreq);
% Per bin decay prototype:
decayBandpass  = zeros(numDesignFreq,numTaps);
superimposedBandpass    = zeros(numChannels,numTaps);

coherentGain = 0.54;
%hannCoherentGain = 0.5;
norm = (1/longestWindow) * (1/coherentGain);

% Create the highpass filter to take over after the decorrelation filters.
% 0.995 factor to nudge the corner frequency down for a nice transition.
highpass = fir1(longestWindow+1,(0.995*designFreqCurrentHz/(Fs/2)),'high');

%% Spatial Structure:
% Diffuse Field Simulation - Q random plane waves.
numPlaneWave    = 1024;
planeWave       = (2 .* rand(numDim,numPlaneWave)) -1;
phaseRandom     = 2 * pi * rand(flat,numDesignFreq,numPlaneWave);
%phaseRandom     = 2 * pi * rand(flat,flat,numPlaneWave);
planeMagnitude  = ones(flat,numPlaneWave);
pwNormalization = 1/sqrt(numPlaneWave);

% Normalize the randomized k vectors:
normPW          = zeros(flat,numPlaneWave);
normPW(flat,:)  = sqrt(planeWave(X,:).^2 + planeWave(Y,:).^2 + planeWave(Z,:).^2);  % length.
planeWave(X,:)  = planeWave(X,:) ./ normPW(flat,:);                % normalize x.
planeWave(Y,:)  = planeWave(Y,:) ./ normPW(flat,:);                % normalize y.
planeWave(Z,:)  = planeWave(Z,:) ./ normPW(flat,:);                % normalize z.
clear normPW;


%% Virtual Array

% Cylindrical coordinates stored in rectangular matrices.
% theta (incRadianArray * ii) proceeds counter clockwise from the positive
% x axis (convention of this code - right horizontal of your screen.)
% virtArrayPressure has three rings, the assumed KH circle and an inner
% and outer ring for computation of the radial pressure gradient.

AZ = 1; EL = 2;
insideRing = 1; ring = 2; outsideRing = 3; numRings = 3;
virtArrayCoord  = zeros(numDim,numChannels,numRings);
virtArrayPressure  = zeros(flat,numChannels,numRings);

inOrOutDeltaR   = 0.0001;  % Distance for radial (normal) dipole spatial derivative.

if(config.isBinaural)
    % Binaural - the diffuse field interference is different b/c there is a rigid noggin in the way.
    % This is a work around to maintain the current paradigm and should be improved in the future.
    arrayDistance = 2 * config.arrayDistance; 
else
    arrayDistance = config.arrayDistance;
end;

insideDeltaR    = arrayDistance - inOrOutDeltaR;
outsideDeltaR   = arrayDistance + inOrOutDeltaR;

% Create the cartesian coordinates for the inside, KH, and outside surfaces.
for ii = 1:numChannels
    % Inner ring
    virtArrayCoord(X,ii,insideRing)    = cos(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * insideDeltaR;
    virtArrayCoord(Y,ii,insideRing)    = sin(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * insideDeltaR;
    virtArrayCoord(Z,ii,insideRing)    = sin(virMicDir(ii,EL)) * insideDeltaR;
    % Array
    virtArrayCoord(X,ii,ring)           = cos(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * arrayDistance;
    virtArrayCoord(Y,ii,ring)           = sin(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * arrayDistance;
    virtArrayCoord(Z,ii,ring)           = sin(virMicDir(ii,EL)) * arrayDistance;
    % Outer ring
    virtArrayCoord(X,ii,outsideRing)    = cos(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * outsideDeltaR;
    virtArrayCoord(Y,ii,outsideRing)    = sin(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * outsideDeltaR;
    virtArrayCoord(Z,ii,outsideRing)    = sin(virMicDir(ii,EL)) * outsideDeltaR;
end;

% Values associated with the virtual array:
radialDeltaPinside      = zeros(flat,numChannels,flat);
radialDeltaPoutside     = zeros(flat,numChannels,flat);
radialDeltaP            = zeros(flat,numChannels,flat);
pressureGradnKH         = zeros(flat,numChannels,flat);
complexTransfer         = zeros(flat,numChannels,flat);

for ll = 1:numDesignFreq
    % Calculate spatial frequency for each "design" bin
    omega   = 2 * pi * designFreqL(ll,flat);
    k       = omega / c;
    
    % virtArrayPressure(flat,numChannels,numRings);
    for qq = 1:numPlaneWave
        virtArrayPressure(flat,:,:) = virtArrayPressure(flat,:,:) + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * ...
            (phaseRandom(flat,ll,qq) + k .* (planeWave(X,qq) * virtArrayCoord(X,:,:) ...
            + planeWave(Y,qq) * virtArrayCoord(Y,:,:) + planeWave(Z,qq) * virtArrayCoord(Z,:,:))));
    end;
    
    radialDeltaPinside(flat,:,flat)      =  virtArrayPressure(flat,:,ring) - virtArrayPressure(flat,:,insideRing);
    radialDeltaPoutside(flat,:,flat)     =  virtArrayPressure(flat,:,outsideRing) - virtArrayPressure(flat,:,ring);
    radialDeltaP(flat,:,flat)            =  0.5 * (radialDeltaPinside(flat,:,flat) + radialDeltaPoutside(flat,:,flat));
    radialDeltaP(flat,:,flat)            = (1/inOrOutDeltaR) * radialDeltaP(flat,:,flat);
    
    pressureGradnKH(flat,:,flat)        = radialDeltaP(flat,:,flat);
    cardioidKH = 0.5 .* (virtArrayPressure(flat,:,ring) - (pressureGradnKH(flat,:,flat)/(jjj*k)));
    
    sincKR0 = sinc((k*arrayDistance)/pi);
    % Assume the origin cardioid contributes a fully correlated value scaled by sinc(kr0).
    % We don't want to touch the values from the RIR, so this component is sinc(kr0) ( * gain = 1.0)
    % Assume the rest of the field to come from the simulated edge values scaled by (1-sinc(kr0)).    
    complexTransfer(flat,:,flat) = (meanDiffDB * (1-sincKR0) .* cardioidKH(flat,:,flat)) + (sincKR0);
    
    % Turn field values into mag and phase:
    expMag(:,ll)      = abs(complexTransfer(flat,:,flat));
    expPhase(:,ll)    = angle(complexTransfer(flat,:,flat));

    % Bound for audio quality.
    expMag(:,ll)    = min(expMag(:,ll),10^(3/20));
    expMag(:,ll)    = max(expMag(:,ll),10^(-20/20));
    
    % Must be zero'd out for every frequency.
    virtArrayPressure  = zeros(flat,numChannels,numRings);
end;

if(validateSpacing)
    % TEST FOR CHANNEL SPACING:
    expMag     = ones(numChannels,numDesignFreq);
    expPhase   = zeros(numChannels,numDesignFreq);
end;

% Create a bank of prototype bandpass filters - the temporal decay depends upon frequency.
% 1) Defined in the time domain as a window modulated by center frequency.
% 2) Apply random (but structured) amplitude and phase values in time.
% 3) Superimpose the prototype bandpass filters.
for ii = 1:numChannels
    for ll = 1:numDesignFreq
        % Modulator scaled by random amplitude (exponential distribution) and random phase  (uniform distribution.)
        accumulator     = ((2*pi*designFreqL(ll,flat)/Fs)*(1:numTaps)) + expPhase(ii,ll);
        % Multiplication by 2 accounts for missing negative frequencies.
        modulator       = expMag(ii,ll) * 2 .* cos(accumulator);
        
        % Create the temporal response of the bandpass.
        decayBandpass(ll,(1:numTaps))   =  timeEnvL(ll,1:numTaps) .* modulator(1:numTaps);
        
        superimposedBandpass(ii,1:numTaps) ...
            = superimposedBandpass(ii,1:numTaps) + decayBandpass(ll,(1:numTaps));
    end;

    if(timeAlign)
        superimposedBandpass(ii,:) = norm .* superimposedBandpass(ii,:);
        
        % Get rid of the unnecessary delay and align with the highpass.
        startIndex = earliestStartIndex;
        endIndex = startIndex + longestWindow - 1;
        shiftIndexes = startIndex:1:endIndex;
        duration = length(shiftIndexes);
        % Shift the bank of bandpasses into output buffer.
        decorrFilt(ii,1:duration) = superimposedBandpass(ii,shiftIndexes);
        % Add in the highpass
        decorrFilt(ii,1:length(highpass)) = decorrFilt(ii,1:length(highpass)) + highpass;
    else
        % The gain 'norm' is defined in the frequency spacing code.
        superimposedBandpass(ii,:) = norm .* superimposedBandpass(ii,:);
        superimposedBandpass(ii,1:length(highpass)) ...
            = superimposedBandpass(ii,1:length(highpass)) + highpass;
        decorrFilt(ii,1:numTaps) = superimposedBandpass(ii,1:numTaps);
    end;

end;


if(1)
    cd ImpulsesAndFilters
    save(filename, 'decorrFilt','numDesignFreq','maxDesignFreq','designFreqL','expMag','expPhase');
    cd ..
end;

if(1)
    cd PrettyPictures; cd Data
    % Save the filters.
    localName = strcat(filename,'_Filters');
    save(localName,'decorrFilt');    
    cd ..; cd ..
end;

viewFilters = 1;
if(viewFilters)
    numTaps = numTaps * 4;
    nyquistIndex    = ((numTaps/2)+1);
    bins2plot       = nyquistIndex/1;
    nDFT            = numTaps;
    freqRes         = Fs/nDFT;                          % (1/T) = cycles/sec per bin.
    freqVec         = freqRes * (0:nyquistIndex-1);     % 0.0 = DC (0Hz) to Nyquist.
    
    viewRange = 1:numChannels;
    for ii = viewRange;
        scnsize = get(0,'ScreenSize');
        figure('Position',scnsize-[-30 -40 60 120]);
        
        subplot(2,1,1);
        evalFreq1 = fft(decorrFilt(ii,:),numTaps);
        semilogx(freqVec(1:bins2plot),20*log10(abs(evalFreq1(1:bins2plot))),'magenta');
        hold on;
        semilogx(designFreqL(1:numDesignFreq,flat),20*log10(expMag(ii,:)),'o-');
        grid on;
        axis([20 2*maxDesignFreq -40 20]);
        subplot(2,1,2);
        semilogx(freqVec(1:bins2plot),angle(evalFreq1(1:bins2plot)),'magenta');
        hold on;
        semilogx(designFreqL(1:numDesignFreq,flat),(expPhase(ii,:)),'o-');
        grid on;
        axis([20 2*maxDesignFreq -4 4]);
        %         pause(16);
        %         close;
    end;
end;


