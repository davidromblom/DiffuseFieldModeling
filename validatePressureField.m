function validatePressureField(frequencyRequested,configuration,externalFieldRNG)


%% Setup and constants
global  jjj;    jjj     = sqrt(-1);
global  omega;                          % value computed below. 
global  c;      c       = 343;
global  k;      k       = omega / c; 
global rho;     rho     = 1.2;

verbosePlots    = 0;
paperPlots      = 1;
saveGraphingData= 1;

% Helps with reading the code:
X = 1; Y = 2; Z = 3; numDim3 = 3; numDim2 = 2;
flat = 1; % makes zeros and ones declarations easier to read.

config = getParameterSet3(configuration);
[ virMicDir, numLdspkr ] = getLoudspeakerConfiguration(config.ldspkrArray);

if(config.numChannels ~= numLdspkr)
    disp('Something is wrong with the configuration specification.');
    return;
end;
filename = strcat('dfm',num2str(configuration));

%% All sorts of parameters:

Fs              = config.Fs;
numTaps         = config.numTaps;
numChannels     = config.numChannels;

% Set deterministic values for the origin measurements.
overrideMeasurements = 0;

radiationKH         = 1; % double layer = 1
radiationOmni       = 2; % omni radiation = 2 
radiationCardioid   = 3; % cardioid radiation = 3

% Set to 1 as test case, will always be +x (hard coded below.)
numPlaneWave    = 1024; % = 1;
nonIdealGain    = 1;  % Gain for +x to demonstrate non-ideal energy distribution. 

% Determine if plane wave idealization should include random amplitude.
unityField      = 1;

% Cartesian grid for graphing.  
edgeDistance    = config.arrayDistance;    % meters - fit array within the square grid
                                        
global spaceSmplCart;   spaceSmplCart   = 128;  % smpls in square grid.  

% Cylindrical coordinates for spatial autocorrelation.
global radianSmplsFIELD;     radianSmplsFIELD     = 64;  % samples of 2pi
global radialSmplsFIELD;     radialSmplsFIELD     = 64;   % effectively the number of rings around the origin.  

% Kirchhoff-Helmholtz in cylindrical coordinates:
% Note createSpaceFilters makes 256 channels, currently the limit (arbitrary.)
global numPortholes;   numPortholes   = config.numChannels;  % samples of 2pi - number of portholes.  Can be different from CC sampling. 
global radialDistance;  radialDistance  = config.arrayDistance;  % KH center distance.  

inOrOutDeltaR   = 0.0001;  % Distance for radial (normal) dipole spatial derivative.
spatialGradient = inOrOutDeltaR;


%% Define Cartesian and Cylindrical Space:
global      spaceCoordCART;      spaceCoordCART   = zeros(numDim2,spaceSmplCart,spaceSmplCart);
increment                   = (2*edgeDistance)/spaceSmplCart;

% Cylindrical coordinates need to be stored in square or rectangular matrices.  
% theta (incRadian) proceeds counter clockwise from the positive x axis.  
global      spaceCoordCYL;      spaceCoordCYL        = zeros(numDim2,radianSmplsFIELD,radialSmplsFIELD);
incRadial                   = edgeDistance/radialSmplsFIELD;
incRadian                   = 2 * pi / radianSmplsFIELD;


%% Create filters and adjust frequency to nearest bin.  

cd ImpulsesAndFilters
load(strcat(filename,'.mat'));
cd ..
[numFilters tapsEqualsNFFT ] = size(decorrFilt);

freqSpacing = Fs / tapsEqualsNFFT;
binIndex = round(frequencyRequested / freqSpacing)+1; % +1 to account for bin 1 = DC
frequency = (binIndex-1) * freqSpacing;
omega   = 2 * pi * frequency;


%% Diffuse Field Simulation 
%   The plane waves are specified in cartesian coordinates.
%   The rest of the algorithm works in cylindrical coordinates.
%   The plane wave field is idealized.  Q random incident plane waves, infinite
%   constant phase surfaces normal to each propagation vector.  

if(externalFieldRNG == 0)
    rng(config.fieldRNG);
else
    rng(externalFieldRNG);
end;

planeWave       = (2 .* rand(numDim3,numPlaneWave)) -1;       
phaseRandom     = 2 * pi * rand(flat,numPlaneWave);

if(unityField)
    planeMagnitude  = ones(flat,numPlaneWave);
else
    planeMagnitude  = (2 .* rand(flat,numPlaneWave)) -1;
end;

% Normalize the randomized k vectors:
normPW          = zeros(flat,numPlaneWave);
normPW(1,:)     = sqrt(planeWave(X,:).^2 + planeWave(Y,:).^2 + planeWave(Z,:).^2);  % length.
planeWave(X,:)  = planeWave(X,:) ./ normPW(flat,:);                % normalize x.
planeWave(Y,:)  = planeWave(Y,:) ./ normPW(flat,:);                % normalize y.
planeWave(Z,:)  = planeWave(Z,:) ./ normPW(flat,:);                % normalize y.
clear normPW;

% x-directed pw is a very good means of understanding the algorithm.
if(numPlaneWave == 1)
    planeWave(X,:)  = 1;
    planeWave(Y,:)  = 0;
    planeWave(Z,:)  = 0;
    phaseRandom = 0;
end;

% Apply spatial frequency to all random plane wave directions.  
planeWave       = k .* planeWave;
pwNormalization = 1/sqrt(numPlaneWave);

% Make all +x plane waves larger in amplitude
for ii = 1:numPlaneWave
    if(planeWave(1,ii) > 0.0)
        planeMagnitude(X,ii) = planeMagnitude(X,ii) * nonIdealGain;
    end;
end;

%% Cartesian Calculations of Idealized Diffuse Field

for i = 1:spaceSmplCart
    spaceCoordCART(X,i,1:spaceSmplCart)    = -edgeDistance + increment * (1:spaceSmplCart);
    spaceCoordCART(Y,1:spaceSmplCart,i)    = -edgeDistance + increment * (1:spaceSmplCart);
end;

idealFieldCART     = zeros(flat,spaceSmplCart,spaceSmplCart); 

for ii = 1:numPlaneWave
    idealFieldCART(flat,:,:) = ...
        idealFieldCART(flat,:,:) + (planeMagnitude(flat,ii)*pwNormalization) .* exp(jjj * ... 
        (phaseRandom(flat,ii) + (planeWave(X,ii) * spaceCoordCART(X,:,:) + planeWave(Y,ii) * spaceCoordCART(Y,:,:))));
end;
% The field is computed on the horizontal plane, so the Z component is flattened on the xy plane. 


%% Multipole projection at origin

spatialDeltaMeters = inOrOutDeltaR;

pressureOrigin  = 0.0;
pressurePlusX   = 0.0;
pressureMinusX  = 0.0;
pressurePlusY   = 0.0;
pressureMinusY  = 0.0;
pressurePlusZ   = 0.0;
pressureMinusZ  = 0.0;

for qq = 1:numPlaneWave
    pressureOrigin  = pressureOrigin    ...
        + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * (phaseRandom(flat,qq)));
    pressurePlusX   = pressurePlusX     ...
        + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * (phaseRandom(flat,qq) + (planeWave(X,qq) * spatialDeltaMeters)));
    pressureMinusX  = pressureMinusX    ...
        + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * (phaseRandom(flat,qq) + (planeWave(X,qq) * -spatialDeltaMeters)));
    pressurePlusY   = pressurePlusY     ...
        + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * (phaseRandom(flat,qq) + (planeWave(Y,qq) * spatialDeltaMeters)));
    pressureMinusY  = pressureMinusY    ...
        + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * (phaseRandom(flat,qq) + (planeWave(Y,qq) * -spatialDeltaMeters)));    
     pressurePlusZ   = pressurePlusZ     ...
         + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * (phaseRandom(flat,qq) + (planeWave(Z,qq) * spatialDeltaMeters)));
     pressureMinusZ  = pressureMinusZ    ...
         + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * (phaseRandom(flat,qq) + (planeWave(Z,qq) * -spatialDeltaMeters)));
end;

dpdx =  (1/2)/spatialGradient * ((pressurePlusX - pressureOrigin) + (pressureOrigin - pressureMinusX));
dpdy =  (1/2)/spatialGradient * ((pressurePlusY - pressureOrigin) + (pressureOrigin - pressureMinusY));
dpdz =  (1/2)/spatialDeltaMeters * ((pressurePlusZ - pressureOrigin) + (pressureOrigin - pressureMinusZ));

if(overrideMeasurements)
    dpdx = 0.0;
    dpdy = 0.0;
    pressureOrigin = 1.0;
end;

%% Cylindrical Calculations of Idealized Diffuse Field

for ii = 1:radianSmplsFIELD
    for jj = 1:radialSmplsFIELD
        % x axis projection
        spaceCoordCYL(X,ii,jj)    = cos(incRadian * ii) * (incRadial * jj);
        % y axis projection
        spaceCoordCYL(Y,ii,jj)    = sin(incRadian * ii) * (incRadial * jj);
    end;
end;

idealFieldCYL  = zeros(flat,radianSmplsFIELD,radialSmplsFIELD); 

for ii = 1:numPlaneWave
    idealFieldCYL(flat,:,:) = idealFieldCYL(flat,:,:) + (planeMagnitude(flat,ii)*pwNormalization) .* exp(jjj * ...
        (phaseRandom(flat,ii) + (planeWave(X,ii) * spaceCoordCYL(X,:,:) + planeWave(Y,ii) * spaceCoordCYL(Y,:,:))));
end;


%% K/H Projection

% Cylindrical coordinates stored in rectangular matrices.  
% theta (incRadianKH * ii) proceeds counter clockwise from the positive x axis.  
% virtArrayPressure has three rings, the assumed KH circle and an inner 
% and outer ring for computation of the radial pressure gradient.

AZ = 1; EL = 2;
insideRing = 1; ring = 2; outsideRing = 3; numRings = 3;
virtArrayCoord  = zeros(numDim3,numPortholes,numRings);

insideDeltaR    = radialDistance - inOrOutDeltaR;
outsideDeltaR   = radialDistance + inOrOutDeltaR;
incRadianKH     = 2 * pi / numPortholes;

% Create the cartesian coordinates for the inside, KH, and outside surfaces.
for ii = 1:numPortholes
    % Inner ring
    virtArrayCoord(X,ii,insideRing)    = cos(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * insideDeltaR;
    virtArrayCoord(Y,ii,insideRing)    = sin(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * insideDeltaR;
    virtArrayCoord(Z,ii,insideRing)    = sin(virMicDir(ii,2)) * insideDeltaR;
    % Array
    virtArrayCoord(X,ii,ring)           = cos(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * config.arrayDistance;
    virtArrayCoord(Y,ii,ring)           = sin(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * config.arrayDistance;
    virtArrayCoord(Z,ii,ring)           = sin(virMicDir(ii,2)) * config.arrayDistance;
    % Outer ring
    virtArrayCoord(X,ii,outsideRing)    = cos(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * outsideDeltaR;
    virtArrayCoord(Y,ii,outsideRing)    = sin(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL)) * outsideDeltaR;
    virtArrayCoord(Z,ii,outsideRing)    = sin(virMicDir(ii,EL)) * outsideDeltaR;
end;

virtArrayPressure  = zeros(flat,numPortholes,numRings); 
for qq = 1:numPlaneWave
    virtArrayPressure(flat,:,:) = virtArrayPressure(flat,:,:) + (planeMagnitude(flat,qq)*pwNormalization) .* exp(jjj * ...
        (phaseRandom(flat,qq) + (planeWave(X,qq) * virtArrayCoord(X,:,:) ...
        + planeWave(Y,qq) * virtArrayCoord(Y,:,:) + planeWave(Z,qq) * virtArrayCoord(Z,:,:))));
end;

radialDeltaPinside      = zeros(flat,numPortholes,flat); 
radialDeltaPoutside     = zeros(flat,numPortholes,flat); 
radialDeltaP            = zeros(flat,numPortholes,flat); 
pressureGradnKH         = zeros(flat,numPortholes,flat); 

radialDeltaPinside(flat,:,flat)      =  virtArrayPressure(flat,:,ring) - virtArrayPressure(flat,:,insideRing);
radialDeltaPoutside(flat,:,flat)     =  virtArrayPressure(flat,:,outsideRing) - virtArrayPressure(flat,:,ring);
radialDeltaP(flat,:,flat)            =  0.5 * (radialDeltaPinside(flat,:,flat) + radialDeltaPoutside(flat,:,flat));
radialDeltaP(flat,:,flat)            = (1/inOrOutDeltaR) * radialDeltaP(flat,:,flat);

% The gradient calculations are "looking outward" and we want "looking inward" -
pressureGradnKH(flat,:,flat)        = -radialDeltaP(flat,:,flat);

%% K/H re-radiation monopoles and dipoles.

[c2cRefCART] = kirchhoffHelmholtzCART(virtArrayPressure(flat,:,ring),...
    pressureGradnKH(flat,:,flat),virtArrayCoord(:,:,ring),radiationCardioid);

[c2cRefCYL] = kirchhoffHelmholtzCYL(virtArrayPressure(flat,:,ring),pressureGradnKH(flat,:,flat),virtArrayCoord(:,:,ring),radiationCardioid);

[fullKHCART] = kirchhoffHelmholtzCART(virtArrayPressure(flat,:,ring),...
    pressureGradnKH(flat,:,flat),virtArrayCoord(:,:,ring),radiationKH);


%% Driving Functions for not decorrelated case:

pressureCombined          = zeros(flat,numPortholes,flat); 
gradientCombined          = zeros(flat,numPortholes,flat); 

    % This needs to be consistent with the K/H array, that is, inward looking.  
    % dpdx is done positive x minus negative x, so this needs to be reversed.
    derivativeSign     = -1; 
    
    for ii = 1:numPortholes
        pressureCombined(ii)    = pressureOrigin * exp(jjj* k * radialDistance); 
        weightedX = dpdx * cos(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL));
        weightedY = dpdy * sin(virMicDir(ii,AZ)) * cos(virMicDir(ii,EL));
        weightedZ = dpdz * sin(virMicDir(ii,EL));
        gradientCombined(ii)    = weightedX + weightedY + weightedZ; % DR - SANITY CHECK.
        gradientCombined(ii)    = derivativeSign * gradientCombined(ii)  * exp(jjj* k * radialDistance);
    end;

    [outwardCardioidRebuildCART] = kirchhoffHelmholtzCART(pressureCombined,gradientCombined,virtArrayCoord(:,:,ring),radiationCardioid);
    [outwardCardioidRebuildCYL]  = kirchhoffHelmholtzCYL(pressureCombined,gradientCombined,virtArrayCoord(:,:,ring),radiationCardioid);


%% Apply Filters to Create Portholes

pressurePorthole          = zeros(flat,numPortholes,flat); 
presGradPorthole          = zeros(flat,numPortholes,flat); 

for ii = 1:numPortholes
    decorrFreq              = fft(decorrFilt(ii,:),tapsEqualsNFFT,2);   
    pressurePorthole(ii)    = pressureCombined(ii) * decorrFreq(binIndex);
    presGradPorthole(ii)    = gradientCombined(ii) * decorrFreq(binIndex);
end;

[portholeRebuildCART] = kirchhoffHelmholtzCART(pressurePorthole,presGradPorthole,virtArrayCoord(:,:,ring),radiationCardioid);
[portholeRebuildCYL]  = kirchhoffHelmholtzCYL(pressurePorthole,presGradPorthole,virtArrayCoord(:,:,ring),radiationCardioid);


%% Spatial Autocorrelation Calculations
% abs(cmplxNumber).^2 = cmplxNumber .* conj(cmplxNumber) 

meanIdeal       = squeeze(mean(mean(idealFieldCYL(1,:,:),2)));
deviation       = squeeze(idealFieldCYL(1,:,:)) - meanIdeal;
idealCC         = mean(pressureOrigin .* conj(deviation(:,1:radialSmplsFIELD)),1) ...
    ./ sqrt(abs(pressureOrigin).^2 .* mean(abs(deviation(:,1:radialSmplsFIELD)).^2,1));

% K/H reconstruction SAC.
pressureOrigin2  = c2cRefCART(1,(spaceSmplCart/2),(spaceSmplCart/2));
meanFullKH       = squeeze(mean(mean(c2cRefCYL(1,:,:),2)));
deviation        = squeeze(c2cRefCYL(1,:,:)) - meanFullKH;
fullKHCC         = mean(pressureOrigin2 .* conj(deviation(:,1:radialSmplsFIELD)),1) ...
    ./ sqrt(abs(pressureOrigin2).^2 .* mean(abs(deviation(:,1:radialSmplsFIELD)).^2,1));

% outwardCardioid SAC:
% Need origin pressure - the point at spaceSmplCart/2 is 0,0 for even numbers.
pressureOrigin2  = outwardCardioidRebuildCART(1,(spaceSmplCart/2),(spaceSmplCart/2));
meanAmbKH       = squeeze(mean(mean(outwardCardioidRebuildCYL(1,:,:),2)));
deviation        = squeeze(outwardCardioidRebuildCYL(1,:,:)) - meanAmbKH;
ambKHCC         = mean(pressureOrigin2 .* conj(deviation(:,1:radialSmplsFIELD)),1) ...
    ./ sqrt(abs(pressureOrigin2).^2 .* mean(abs(deviation(:,1:radialSmplsFIELD)).^2,1));

% Porthole SAC
% Need origin pressure - the point at spaceSmplCart/2 is 0,0 for even numbers.
pressureOrigin2= portholeRebuildCART(1,(spaceSmplCart/2),(spaceSmplCart/2));
meanPort       = squeeze(mean(mean(portholeRebuildCYL(1,:,:),2)));
deviation      = squeeze(portholeRebuildCYL(1,:,:)) - meanPort;
portKHCC         = mean(pressureOrigin2 .* conj(deviation(:,1:radialSmplsFIELD)),1) ...
    ./ sqrt(abs(pressureOrigin2).^2 .* mean(abs(deviation(:,1:radialSmplsFIELD)).^2,1));


%% End of calculations - below is plotting or saving data.  Radiation functions at the very bottom.  

if(verbosePlots)
    kr = (1:radialSmplsFIELD) * incRadial * k;
    figure();

    subplot(2,2,1); semilogx(kr,squeeze(abs(idealCC))); hold on;
    semilogx(kr,squeeze(real(idealCC)),'black');
    semilogx(kr,squeeze(imag(idealCC)),'red');
    title('Reference Field'); xlabel('kr -- cycles'); ylabel('SACC');
    axis([kr(1) kr(radialSmplsFIELD) -1 1 ]); grid on; hold off;

    subplot(2,2,2); semilogx(kr,squeeze(abs(fullKHCC))); hold on;
    semilogx(kr,squeeze(real(fullKHCC)),'black');
    semilogx(kr,squeeze(imag(fullKHCC)),'red');
    title('K/H Reconstruction'); xlabel('kr -- cycles'); ylabel('SACC');
    axis([kr(1) kr(radialSmplsFIELD) -1 1 ]); grid on; hold off;

    subplot(2,2,3); semilogx(kr,squeeze(abs(portKHCC))); hold on;
    semilogx(kr,squeeze(real(portKHCC)),'black');
    semilogx(kr,squeeze(imag(portKHCC)),'red');
    title('Porthole'); xlabel('kr -- cycles'); ylabel('SACC');
    axis([kr(1) kr(radialSmplsFIELD) -1 1 ]); grid on; hold off;

    subplot(2,2,4); semilogx(kr,squeeze(abs(ambKHCC))); hold on;
    semilogx(kr,squeeze(real(ambKHCC)),'black');
    semilogx(kr,squeeze(imag(ambKHCC)),'red');
    title('outwardCardioid'); xlabel('kr -- cycles'); ylabel('SACC');
    axis([kr(1) kr(radialSmplsFIELD) -1 1 ]); grid on; hold off;
end;



%% Plot Reconstructed Fields
if(verbosePlots)
    colorRange = [-1 1];
    
    figure();
    subplot(2,2,1)
    pcolor(squeeze(spaceCoordCART(1,:,:)),squeeze(spaceCoordCART(2,:,:)),squeeze(real(idealFieldCART(1,:,:))));
    shading(gca,'interp');
    titleString = strcat('Reference Field.  Actual Frequency Used =',num2str(frequency));
    title(titleString);
    caxis(colorRange); grid on; set(gca,'layer','top'); colorbar;
    
    subplot(2,2,2)
    pcolor(squeeze(spaceCoordCART(1,:,:)),squeeze(spaceCoordCART(2,:,:)),squeeze(real(c2cRefCART(1,:,:))));
    shading(gca,'interp');
    title('Kirchhoff Helmholtz Ideal Array Reconstruction');
    caxis(colorRange); grid on; set(gca,'layer','top'); colorbar;
    
    subplot(2,2,3)
    pcolor(squeeze(spaceCoordCART(1,:,:)),squeeze(spaceCoordCART(2,:,:)), ...
        squeeze(real(portholeRebuildCART(1,:,:))));
    shading(gca,'interp');
    title('Diffuse Field Synthesis Reconstruction');
    caxis(colorRange); grid on; set(gca,'layer','top'); colorbar;
    
    subplot(2,2,4)
    pcolor(squeeze(spaceCoordCART(1,:,:)),squeeze(spaceCoordCART(2,:,:)),squeeze(real(outwardCardioidRebuildCART(1,:,:))));
    shading(gca,'interp');
    title('Transposed Microphone Without Decorrelation');
    caxis(colorRange); grid on; set(gca,'layer','top'); colorbar;
end;


if(saveGraphingData)
    cd PrettyPictures; cd Data
        % Save the spatial autocorrelation data.
        localName = strcat(filename,'_SAC',num2str(frequencyRequested));
        save(localName,'idealCC','fullKHCC','ambKHCC','portKHCC');
        % Save the barnacles on the fields - 
        localName = strcat(filename,'_FIELDS',num2str(frequencyRequested));
        save(localName,'idealFieldCART','fullKHCART','c2cRefCART','outwardCardioidRebuildCART','portholeRebuildCART','spaceCoordCART');

        edgeMeasuredPres = squeeze(abs(virtArrayPressure(1,:,2)));
        edgeMeasuredGrad = squeeze(abs(pressureGradnKH(1,:,1)));

        % Save the edge values
        localName = strcat(filename,'_EDGE',num2str(frequencyRequested));
        save(localName,'edgeMeasuredPres','edgeMeasuredGrad'); 

    cd .. ; cd ..
end;

if(paperPlots)
    colorRange = [-1 1];
    scnsize = get(0,'ScreenSize');
    plotAll = figure('Position',scnsize-[-30 -40 200 500]);

    plotHandel = subplot(1,3,1);
    pcolor(squeeze(spaceCoordCART(1,:,:)),squeeze(spaceCoordCART(2,:,:)),squeeze(real(c2cRefCART(1,:,:))));
    shading(gca,'interp');
    title('Measured');  axis square;
    xlabel('meters');
    caxis(colorRange); grid on; set(gca,'layer','top'); 
    
    subplot(1,3,2)    
    pcolor(squeeze(spaceCoordCART(1,:,:)),squeeze(spaceCoordCART(2,:,:)),squeeze(real(portholeRebuildCART(1,:,:))));
    shading(gca,'interp');
    title('DFM'); axis square;
    xlabel('meters');
    caxis(colorRange); grid on; set(gca,'layer','top');
    
     % Create colorbar
    h2 = colorbar('peer',plotHandel,'Position',...
        [0.0601223241590216 0.133333333333334 0.0206880733944954 0.766666666666668]);
    LatexReP = ['$Re[\hat{\mathbf{p}}]$'];
    set(get(h2,'title'),'string',LatexReP,'Interpreter','Latex','fontSize',12);

    subplot(1,3,3)    
    pcolor(squeeze(spaceCoordCART(1,:,:)),squeeze(spaceCoordCART(2,:,:)),squeeze(real(outwardCardioidRebuildCART(1,:,:))));
    shading(gca,'interp');
    title('None'); axis square;
    xlabel('meters');
    caxis(colorRange); grid on; set(gca,'layer','top');
    
end;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN FUNCTION.  RADIATION FUNCTIONS BELOW.




function [fieldCART] ...
    = kirchhoffHelmholtzCART(radialPressure,pressureGradient,portholeCoords,radiationCondition)

global spaceSmplCart;
global numPortholes;
global spaceCoordCART;
global jjj; global omega; global c; global k; global rho;
global radialDistance;

X = 1; Y = 2; Z = 3; numDim3 = 3; numDim2 = 2;
flat = 1; % makes zeros and ones declarations easier to read.

monopoleKHCART      = zeros(flat,spaceSmplCart,spaceSmplCart); 
dipoleKHCART        = zeros(flat,spaceSmplCart,spaceSmplCart); 
singleKHCART        = zeros(flat,spaceSmplCart,spaceSmplCart); 
distanceKHCART      = zeros(flat,spaceSmplCart,spaceSmplCart);
dotKHCART           = zeros(flat,spaceSmplCart,spaceSmplCart);

RadiationOrder0      = zeros(flat,spaceSmplCart,spaceSmplCart); 
RadiationOrder1      = zeros(flat,spaceSmplCart,spaceSmplCart); 

for ii = 1:numPortholes    
    % Distance from secondary source location to arbitrary grid location:
    distanceKHCART(flat,:,:)   = sqrt((portholeCoords(X,ii) - spaceCoordCART(X,:,:)).^2 ...
                                    + (portholeCoords(Y,ii) - spaceCoordCART(Y,:,:)).^2 ...;
                                    + (portholeCoords(Z,ii) - 0.0).^2);

    % Set a minimum distace to avoid pole singularity.  
    distanceKHCART(flat,:,:)   = max(distanceKHCART(flat,:,:),.001);
    
    % Green's function (definition from G/D, not AP)
    RadiationOrder0(flat,:,:)   = (1./(4 .* pi .* (distanceKHCART(flat,:,:)))) .* exp(j .* (k .* distanceKHCART(flat,:,:)));
    % Approximation of 1st derivative => jk * G(r)
    RadiationOrder1(flat,:,:)   = (j*k)./(4 .* pi .* (distanceKHCART(flat,:,:))) .* exp(j .* (k .* distanceKHCART(flat,:,:)));

    % Dot product of inward surface normal and vector from secondary source location to arbitrary grid location.
    % AP uses vector from observation to source.  I use source to observation.  
    dotKHCART(1,:,:)      = ... 
        (-portholeCoords(X,ii)/radialDistance) .* ((spaceCoordCART(X,:,:) - portholeCoords(Z,ii)) ./ distanceKHCART(flat,:,:)) ...
        + (-portholeCoords(Y,ii)/radialDistance) .* ((spaceCoordCART(Y,:,:) - portholeCoords(Y,ii)) ./ distanceKHCART(flat,:,:)) ...
        + (-portholeCoords(Z,ii)/radialDistance) .* ((0.0 - portholeCoords(Z,ii)) ./ distanceKHCART(flat,:,:));
    
    % Explicit double layer K/H
    if(radiationCondition == 1) 
        monopoleKHCART(flat,:,:)   = (1/(j*k)) .* pressureGradient(1,ii) .* RadiationOrder0(flat,:,:);
        
        % Monopole driving dipole sources:
        dipoleKHCART(flat,:,:)     = (radialPressure(1,ii) .* (1/(j*k)) .* dotKHCART(1,:,:)) .* RadiationOrder1(flat,:,:);
    
        % Addition is correct here given change of vector AP to DR.
        singleKHCART(flat,:,:)   = singleKHCART(flat,:,:) + (monopoleKHCART(flat,:,:) + dipoleKHCART(flat,:,:));
        

    % Cardioid to Omni
    elseif(radiationCondition == 2)
        drivingFunction = (radialPressure(1,ii) + (1/(j*k)) .* pressureGradient(1,ii));
        % Omni Radiation:
        singleKHCART(flat,:,:)   = singleKHCART(flat,:,:) ...
            +  drivingFunction .* (RadiationOrder0(flat,:,:));

    % Cardioid to Cardioid
    elseif(radiationCondition == 3)
        drivingFunction = (radialPressure(1,ii) + (1/(j*k)) .* pressureGradient(1,ii));

        % Cardioid Radiation:
        singleKHCART(flat,:,:)   = singleKHCART(flat,:,:) ...
            +  (1/2) * drivingFunction .* (RadiationOrder0(flat,:,:) + (1/(j*k)).* (dotKHCART(1,:,:) .* RadiationOrder1(flat,:,:)));
 
    end;
end;

% dS = 2pi * radialDistance/ numPortholes.  EH uses r d(theta) as dS, AP uses dS.  
% Dividing by the the number of secondary sources yields correct results.
% The compensation by j*k is due to my division of 1/jk in the C2C formulation.  

% The stationary phase Equation from RN is sqrt((2*pi)/(j*k)).
stationaryPhase     = sqrt(-(2*pi)/(j*k)); % This works for C2C.  
integralGain        = j*k * stationaryPhase * (2*pi*radialDistance/(numPortholes)) * (-1);
fieldCART(1,:,:)    = integralGain .* singleKHCART(flat,:,:);

return;
 

function [fieldCYL] = kirchhoffHelmholtzCYL(radialPressure,pressureGradient,portholeCoords,radiationCondition)

global numPortholes 
global radialSmplsFIELD;
global radianSmplsFIELD;
global spaceCoordCYL;
global jjj; global omega; global c; global k; global rho;
global radialDistance;

% Helps with reading the code:
X = 1; Y = 2; Z = 3; numDim = 3; flat = 1;

numRadian = radianSmplsFIELD;
numRadial = radialSmplsFIELD;

monopoleKHCYL       = zeros(flat,radianSmplsFIELD,radialSmplsFIELD); 
dipoleKHCYL         = zeros(flat,radianSmplsFIELD,radialSmplsFIELD); 
singleKHCYL         = zeros(flat,radianSmplsFIELD,radialSmplsFIELD); 
fieldCYL            = zeros(flat,radianSmplsFIELD,radialSmplsFIELD); 
distanceKHCYL       = zeros(flat,radianSmplsFIELD,radialSmplsFIELD);
dotKHCYL            = zeros(flat,radianSmplsFIELD,radialSmplsFIELD);
 
RadiationOrder0      = zeros(flat,radianSmplsFIELD,radialSmplsFIELD); 
RadiationOrder1      = zeros(flat,radianSmplsFIELD,radialSmplsFIELD); 

for ii = 1:numPortholes    
    % Distance from secondary source location to arbitrary grid location:
    distanceKHCYL(flat,:,:)   = sqrt((portholeCoords(X,ii) - spaceCoordCYL(X,:,:)).^2 ...
                                    + (portholeCoords(Y,ii) - spaceCoordCYL(Y,:,:)).^2 ...;
                                    + (portholeCoords(Z,ii) - 0.0).^2);
    
    % Set a minimum distace to avoid pole singularity.  
    distanceKHCYL(flat,:,:)   = max(distanceKHCYL(1,:,:),.001);
    
    % Green's function (definition from G/D, not AP)
    RadiationOrder0(flat,:,:)   = (1./(4 .* pi .* (distanceKHCYL(flat,:,:)))) .* exp(j .* (k .* distanceKHCYL(flat,:,:)));
    % Approximation of 1st derivative => jk * G(r)
    RadiationOrder1(flat,:,:)   = (j*k)./(4 .* pi .* (distanceKHCYL(flat,:,:))) .* exp(j .* (k .* distanceKHCYL(flat,:,:)));

    % Dot product of inward surface normal and vector from secondary source location to arbitrary grid location.
    % AP uses vector from observation to source.  I use source to observation.  
    dotKHCYL(flat,:,:)      = ... 
        (-portholeCoords(X,ii)/radialDistance) .* ((spaceCoordCYL(X,:,:) - portholeCoords(X,ii)) ./ distanceKHCYL(flat,:,:)) ...
        + (-portholeCoords(Y,ii)/radialDistance) .* ((spaceCoordCYL(Y,:,:) - portholeCoords(Y,ii)) ./ distanceKHCYL(flat,:,:)) ...
        + (-portholeCoords(Z,ii)/radialDistance) .* ((0.0 - portholeCoords(Z,ii)) ./ distanceKHCYL(flat,:,:));
    
    % Explicit double layer K/H
    if(radiationCondition == 1) 
        monopoleKHCYL(flat,:,:)   = monopoleKHCYL(flat,:,:) ...
            + (1/(j*k)) .* pressureGradient(flat,ii) .* RadiationOrder0(flat,:,:);
        
        % Monopole driving dipole sources:
        dipoleKHCYL(flat,:,:)     = dipoleKHCYL(flat,:,:) + ...
            (radialPressure(1,ii) .* (1/(j*k)) .* dotKHCYL(flat,:,:)) .* RadiationOrder1(flat,:,:);
    
        singleKHCYL(flat,:,:)   = singleKHCYL(flat,:,:) + (monopoleKHCYL(flat,:,:) + dipoleKHCYL(flat,:,:));
        

    % Cardioid to Omni
    elseif(radiationCondition == 2)
        drivingFunction = (radialPressure(1,ii) + (1/(j*k)) .* pressureGradient(1,ii));
        % Omni Radiation:
        singleKHCYL(flat,:,:)   = singleKHCYL(flat,:,:) ...
            +  drivingFunction .* (RadiationOrder0(flat,:,:));

    % Cardioid to Cardioid
    elseif(radiationCondition == 3)
        drivingFunction = (radialPressure(1,ii) + (1/(j*k)) .* pressureGradient(1,ii));

        % Cardioid Radiation:
        singleKHCYL(flat,:,:)   = singleKHCYL(flat,:,:) ...
            +  (1/2) * drivingFunction .* (RadiationOrder0(flat,:,:) + (1/(j*k)).* (dotKHCYL(flat,:,:) .* RadiationOrder1(flat,:,:))); 
    end;
end;

% dS = 2pi * radialDistance/ numPortholes.  EH uses r d(theta) as dS, AP uses dS.  
% Dividing by the the number of secondary sources yields correct results.
% The compensation by k is due to my division of 1/jk in the C2C formulation.  

% The stationary phase Equation from RN is sqrt((2*pi)/(j*k)).
stationaryPhase     = sqrt(-(2*pi)/(j*k)); % This works for C2C.  
integralGain        = j*k * stationaryPhase * (2*pi*radialDistance/(numPortholes)) * (-1);
fieldCYL(1,:,:)     = integralGain .* singleKHCYL(1,:,:);

return;


