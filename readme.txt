A few notes:
- BSD 3-clause license.

- Most of the functions resynthesize file names as a sort of versioning sanity check. Manually renaming things will likely cause things to break.

- This toolkit corresponds to Part 3 of the DFM papers, URL: http://www.aes.org/e-lib/browse.cfm?elib=19362

- Citation: David Romblom. Diffuse Field Modeling using Physically-Inspired Decorrelation Filters and B-Format Microphones: Improvements to the Filter Design Method. The Journal of the Audio Engineering Society, 65(11):943–953, November 2017.

- You’ll need to set up these subfolders:
—- ImpulsesAndFilters
—- PrettyPictures/Data
—- BinauralAudio
—- OutputAudioWAV

- Start by creating filters for configuration 1. Then createAudioStreams, then synthesizeBinaural, then listenBinaural.  After running through it once, it should be easy to create loudspeaker filters (no audition tool for that).


createFiltersDFM(configuration,rngSeed) - this function creates the ensemble of filters for the configuration number specified in getParameterSet(configuration). You will occasionally get fickle filter sets (boomy mid-low frequencies), trying a different random number seed usually fixes this.

getParameterSet(configuration) - contains a number of parameters for the filter generation and listening.  The number of channels and array distance are specified here while the angles are specified in getLoudspeakerConfiguration(). The RT60 estimates are specified as a series of breakpoints (freq,RT60est) and MATLAB’s interpolation function does the rest. “roomScale” is used in the listening functions to stretch the direct and ER. If you’re running a formal experiment, you should be more precise in your definitions. The “isBinaural” flag will multiply the array distance by 2 - spaced omni CC is similar to binaural IACC when the distance is doubled.

getLoudspeakerConfiguration(ldspkrSystem) - specifies the angles of the loudspeakers.  The order of your loudspeakers defines the order of the filters. The loudspeaker configuration is a field in the parameter set (field called ldspkrArray).

createAudioStreams(configuration,source,decorrON) - this convolves source material with the individual channel filter, which is itself the convolution of the B-Format RIR * DFM_i. If you’d like, you can turn the DFM filters on and off with the argument “decorrON”. For most arrays, leaving the filters off will be very unpleasant. The switch to drop the 1st-order harmonics is called “zeroOrder”.  This sounds surprisingly awful even with decorrelation filters on and deserves some more research.

synthesizeBinaural(configuration,source,decorrON) - this function synthesizes the direct and early reflections for binaural listening.  Note that there is not a corresponding loudspeaker function in this toolkit.  The diffuse output of “createAudioStreams” has a simple *HPF* source model to account for the fact that most instruments can’t push a lot of a air at low frequencies. The ER have both the HPF and an additional LPF.

listenBinauralDFM(configuration,source,decorrON,directGainDB,reflectionGainDB,roomGainDB) - this function plays the binaural D/ER/LF stems created in “synthesizeBinaural”. The level of the RIR largely determines these gains - for the Tanna and Pollack RIR I’ve used ~24B of reduction on the D/ER based on the RIR dB plot of the W signal.

validateFiltersDFM(configuration,FACbaseFreq,FAC1orSAC2orFILT3) - this function performs the validation of the filters used in Part 3, JAES 2017. The field FAC1orSAC2orFILT3 switches between frequency autocorrelation = 1, pairwise spatial autocorrelation = 2, or viewing the magnitude response = 3.  FAC requires a base frequency “FACbaseFreq”, leave a dummy variable here for other calculations. I used 200Hz in the paper. While this was surprisingly well behaved, I do recall some occasionally fickle results which I attributed to poor RT60 estimates.

validatePressureField(frequencyRequested,configuration,externalFieldRNG) - this function plots reconstructed pure tone pressure fields and also performs the SAC field calculation used in DFM1 and DFM3. There is a 2D implementation of the double layer Kirchhoff Helmholtz integral that is hopefully useful to some researchers.






