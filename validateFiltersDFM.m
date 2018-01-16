function validateFiltersDFM(configuration,FACbaseFreq,FAC1orSAC2orFILT3)


%% Setup and constants
jjj     = sqrt(-1);
c       = 343;
rho     = 1.2;

% readability utilities:
flat = 1;

config = getParameterSet3(configuration);
[ virMicDir, numLdspkr ] = getLoudspeakerConfiguration(config.ldspkrArray);

if(config.numChannels ~= numLdspkr)
    disp('Something is wrong with the configuration specification.');
    return;
end;

% Copy frequently used stuff into local variables for readibility.
Fs              = config.Fs;
speedofSound    = config.sos;
numTaps         = config.numTaps;
numChannels     = config.numChannels;

filename = strcat('dfm',num2str(configuration));
cd ImpulsesAndFilters
load(strcat(filename,'.mat'));
cd ..

cd RoomImpulses; cd PrunedRIR
if(config.hall == 1)
    load('pollackBrightWXYZclipped.mat');
elseif(config.hall == 2)
    load('tannaWXYZclipped.mat');
else
    disp('Did not find requested RIR.');
    return;
end;
cd ../..


nyquistIndex    = ((numTaps/2)+1);
bins2plot       = nyquistIndex;
nDFT            = numTaps;
freqRes         = Fs/nDFT;                          % (1/T) = cycles/sec per bin.
freqVec         = freqRes * (0:nyquistIndex-1);     % 0.0 = DC (0Hz) to Nyquist.
kVec            = freqVec .* (2*pi/c);

%% View frequency autocorrelation (FAC)
if(FAC1orSAC2orFILT3 == 1) 
    scnsize = get(0,'ScreenSize');
    plotAll = figure('Position',scnsize-[-30 -40 500 500]);
    
    controlFrequencies  = config.controlFreq;
    controlRT6O         = config.controlRT60;
    
    upsampleSpectrum = 16;
    numTapsUP = numTaps * upsampleSpectrum;
    nyquistIndexUP  = ((numTapsUP/2)+1);
    bins2plot       = nyquistIndexUP/1;
    nDFTup          = numTapsUP;
    freqResUP       = Fs/nDFTup;                        % (1/T) = cycles/sec per bin.
    freqVecUP       = freqResUP * (0:nyquistIndexUP-1); % 0.0 = DC (0Hz) to Nyquist.
    
    binStart = floor(FACbaseFreq/freqResUP);
    numBins = 1024;
    FACarray = zeros(numChannels,numBins+1);
    sumSpec = zeros(1,numBins+1);
    for ii = 1:numChannels
        tempSpec = fft(decorrFilt(ii,1:numTaps),numTapsUP);
        tempFreqSeq = abs(tempSpec(binStart:(binStart+numBins)));
        [ FACarray(ii,:), ~, ~ ] = autocorr(tempFreqSeq,numBins);
        sumSpec(1,:) = sumSpec(1,:) + tempFreqSeq;
    end;
    averageFAC = mean(FACarray(:,:),1);
    meanSpec = sumSpec(1,:) ./ numChannels;
    
        RT60 = interp1(controlFrequencies,controlRT6O,FACbaseFreq);
    % Relationship is 13.82 for squared quantities, 6.9 for linear quantities.
    tau = RT60 / 13.82;
    
    deltaOmega = 2*pi*freqVecUP(binStart:(binStart+numBins));
    deltaOmega = deltaOmega - 2*pi*freqVecUP(binStart);

    physicalFAC = 1./(1+(deltaOmega*tau).^2);
    plot(freqVecUP(binStart:(binStart+numBins)),physicalFAC,'k');
    grid on; hold on;
    axis([FACbaseFreq (FACbaseFreq+10) 0 1 ]);

    % DFM curve second for the paper's consistency
    plot(freqVecUP(binStart:(binStart+numBins)),abs(averageFAC),'k--');   

    % Compare to physical FAC from RIR of target room
    card(1,:) = timeCompW + timeCompX;
    card(2,:) = timeCompW - timeCompX;
    card(3,:) = timeCompW + timeCompY;
    card(4,:) = timeCompW - timeCompY;
    card(5,:) = timeCompW + timeCompZ;
    card(6,:) = timeCompW - timeCompZ;
    num2av = 6;

    physSpec = fft(card(1,:),numTapsUP);
    physFreqSeq = abs(physSpec(binStart:(binStart+numBins)));
    physFAC = autocorr(physFreqSeq,numBins);
    for ff = 2:num2av
        physSpec = fft(card(ff,:),numTapsUP);
        physFreqSeq = abs(physSpec(binStart:(binStart+numBins)));
        physFAC = physFAC + autocorr(physFreqSeq,numBins);
    end;

    avPhysFAC = physFAC ./ num2av;

    plot(freqVecUP(binStart:(binStart+numBins)),abs(avPhysFAC),'k:');
    legend1 = legend('Theoretical','DFM','Measured');
    set(legend1,...
    'Position',[0.691025641025641 0.643333333333333 0.201616023137019 0.250000000000002],...
    'FontSize',12);
    xlabel('Frequency Hz','FontSize',18);
    ylabel('Frequency Autocorrelation','FontSize',18);



%% View Spatial Autocorrelation (SAC)
elseif (FAC1orSAC2orFILT3 == 2)   
    channels = 1:numChannels;
    pairwise = nchoosek(channels,2);
    numPairs = nchoosek(numChannels,2);
    
    ccByPair = zeros(numPairs,nyquistIndex);
    krVec    = zeros(numPairs,nyquistIndex);
    spacing  = zeros(numPairs,flat);
    
    for pair = 1:numPairs
        temp1           = fft(decorrFilt(pairwise(pair,1),1:numTaps),numTaps);
        firstChannel    = temp1(1:nyquistIndex);
        temp2           = fft(decorrFilt(pairwise(pair,2),1:numTaps),numTaps);
        secondChannel   = temp2(1:nyquistIndex);
        crossCorrNumer  = firstChannel .* conj(secondChannel);
        crossCorrDenom  = sqrt(abs(firstChannel).^2 .* abs(secondChannel).^2);
        ccByPair(pair,:)= crossCorrNumer ./ crossCorrDenom;
        % kr spectrum depends on distance between the pair:
        virtArrayCoordX1    = cos(virMicDir(pairwise(pair,1),1)) * cos(virMicDir(pairwise(pair,1),2)) * config.arrayDistance;
        virtArrayCoordX2    = cos(virMicDir(pairwise(pair,2),1)) * cos(virMicDir(pairwise(pair,2),2)) * config.arrayDistance;
        deltaX = virtArrayCoordX1 - virtArrayCoordX2;
        virtArrayCoordY1    = sin(virMicDir(pairwise(pair,1),1)) * cos(virMicDir(pairwise(pair,1),2)) * config.arrayDistance;
        virtArrayCoordY2    = sin(virMicDir(pairwise(pair,2),1)) * cos(virMicDir(pairwise(pair,2),2)) * config.arrayDistance;
        deltaY = virtArrayCoordY1 - virtArrayCoordY2;
        virtArrayCoordZ1    = sin(virMicDir(pairwise(pair,1),2)) * config.arrayDistance;
        virtArrayCoordZ2    = sin(virMicDir(pairwise(pair,2),2)) * config.arrayDistance;
        deltaZ = virtArrayCoordZ1 - virtArrayCoordZ2;
        distance = sqrt(deltaX^2 + deltaY^2 + deltaZ^2);
        krVec(pair,:)       = kVec .* distance;
        spacing(pair,flat)  = distance;
    end;
    
    maxDesignFreq = 3400; % add to config?
    maxDesignK = 2*pi*maxDesignFreq / c;
    numBinsAv = 8;
    numAvBands = floor(nyquistIndex/numBinsAv);
    ccCoalated = zeros(flat,numAvBands);
    krVecCoalated = zeros(flat,numAvBands);
    counter = 0;
    for pair = 1:numPairs
        if(spacing(pair,flat) < 1.5)
            counter = counter + 1;
            for bb = 1:numAvBands
                binStart = (bb-1) * numBinsAv + 1;
                binEnd = binStart + numBinsAv - 1;
                indexes = find((krVec(pair,:) > kVec(binStart)) & (krVec(pair,:) < kVec(binEnd)) & (krVec(pair,:) < maxDesignK));
                bandAv = mean(ccByPair(pair,indexes));
                ccCoalated(flat,bb) = ccCoalated(flat,bb) + bandAv;
                krVecCoalated(flat,bb) = (kVec(binStart) + kVec(binEnd))/2;
            end;
        end;
    end;
    ccCoalated(flat,:) = ccCoalated(flat,:)./counter;

    scnsize = get(0,'ScreenSize');
    plotAll = figure('Position',scnsize-[-30 -40 500 500]);

    sincKR = sinc((1.0/pi) * krVecCoalated(flat,:));
    semilogx(krVecCoalated(flat,:),abs(sincKR),'k','LineWidth',2);
    grid on; hold on;
    axis([.1 100 0 1]);

    semilogx(krVecCoalated(flat,1:floor(numAvBands/4)),abs(ccCoalated(flat,1:floor(numAvBands/4))),'k--');
    
    legend1 = legend('Theoretical','DFM');
    set(legend1,...
    'Position',[0.691025641025641 0.643333333333333 0.201616023137019 0.250000000000002],...
    'FontSize',12);
    xlabel('Spatial Frequency times Spacing - kr','FontSize',18);
    ylabel('Spatial Autocorrelation','FontSize',18);

elseif (FAC1orSAC2orFILT3 == 3)
   
    numTaps = numTaps * 4;
    nyquistIndex    = ((numTaps/2)+1);
    bins2plot       = nyquistIndex/1;
    nDFT            = numTaps;
    freqRes         = Fs/nDFT;                          % (1/T) = cycles/sec per bin.
    freqVec         = freqRes * (0:nyquistIndex-1);     % 0.0 = DC (0Hz) to Nyquist.
    
    viewRange = 1:numChannels;
    for ii = viewRange;
        scnsize = get(0,'ScreenSize');
        figure('Position',scnsize-[-30 -40 60 360]);
        
        evalFreq1 = fft(decorrFilt(ii,:),numTaps);
        semilogx(freqVec(1:bins2plot),20*log10(abs(evalFreq1(1:bins2plot))),'k','LineWidth',1.5);
        hold on;
        semilogx(designFreqL(1:numDesignFreq,flat),20*log10(expMag(ii,:)),'ko');
        grid on;
        % axis([100 2000 -25 5]);
        axis([20 1000 -25 5]);
        xlabel('Frequency Hz','FontSize',18); ylabel('dB Magnitude','FontSize',18);
        legend1 = legend('Filter Response','Target Response');
        set(legend1,'Location','Southwest','FontSize',24);
    end;

end;

% Plot the spectrum of each cardioid:
if(0)
    figure();
    physSpec = fft(card(1,:),numTapsUP);
    semilogx(freqVecUP(1:nyquistIndexUP),20*log10(abs(physSpec(1:nyquistIndexUP))));
    hold on; grid on;
    avSpec = physSpec;
    for ff = 2:6
        physSpec = fft(card(ff,:),numTapsUP);
        semilogx(freqVecUP(1:nyquistIndexUP),20*log10(abs(physSpec(1:nyquistIndexUP))));
        hold on; grid on;
        avSpec = avSpec + physSpec;
    end;
    avSpec = avSpec ./ 6;
    semilogx(freqVecUP(1:nyquistIndexUP),20*log10(abs(avSpec(1:nyquistIndexUP))),'black');
    axis([20 1000 -30 5]);
end;



















% WTF - this is impossible with binaural.......
if(0)
    %for baseFreq = numBinsToAverage * (0:nyquistIndex-1)
    for bin = 1:nyquist
        if(numChannels == 2)
            SACarrayTemp(1,bin) = tempSpec(1,bin) .* tempSpec(2,bin);
            SACarrayTemp(2,bin) = tempSpec(2,bin) .* tempSpec(1,bin);
        else
            [ SACarrayTemp(:,bin), ~, ~ ] = autocorr(tempSpec(:,bin),numChannels-1);
        end;
    end;
    % baseFreqDOWN = floor(baseFreq / upsampleSpectrum) + 1;
    averageSAC(flat,baseFreqDOWN) = mean(mean(SACarrayTemp(:,:),2),1);
    %end;
    figure();
    plot(freqVec,abs(averageSAC(flat,:)));
    grid on; hold on;
end;






