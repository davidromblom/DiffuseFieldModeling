function indexDiffuse = computeEchoDensity2(bformatW,approximateDelay)

lengthNED   = 1025;
computeNED  = length(bformatW)/4;
nED         = zeros(1,computeNED);
windowNED   = hann(lengthNED);
windowNED   = windowNED ./(sum(windowNED));

for i = approximateDelay:1:(length(bformatW)/8)
    
    % MATLAB doc says that t corresponds to the center of the window.
    hopCenter   = i;
    % Time aperature of the window:
    leftEdge    = hopCenter - floor(lengthNED/2);
    rightEdge   = hopCenter + floor(lengthNED/2);
    
    % MATLAB was whining about these as (integer-valued) doubles.
    leftEdge    = cast(leftEdge,'int32');
    rightEdge    = cast(rightEdge,'int32');
    
    % Compute standard deviation of the window
    windowedSig     = (bformatW(leftEdge:rightEdge).^2) .* windowNED;
    sumOfSquares    = sum(windowedSig);
    sigma           = sqrt(sumOfSquares);
        
    for j = 1:1:lengthNED
        if(abs(bformatW(leftEdge + j - 1,1)) > sigma)
            nED(i) = nED(i) + windowNED(j);
        end;
    end;
    % Coefficient from Abel and Huang.
    nED(i)  = (1/0.3173) * nED(i);
        
end;

figure();
viewIndices = approximateDelay:1:(length(bformatW)/8);
plot(viewIndices,nED(viewIndices));
title('Normalized Echo Density vs. Hop Index');
grid on; hold on;

% Long filter - control rate signal coded at audio rate.
order = 512;
preroll = order;
a = zeros(1,order);
a(1) = 1.0;
b = ones(1,order) / order;

nEDfiltered = filter(b,a,nED);
plot(viewIndices,nEDfiltered(viewIndices),'red');
    

indexDiffuse = find(nEDfiltered(1:computeNED) > 0.95,1,'first');
