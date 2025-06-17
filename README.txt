README

Code for dyanmic spatiotemporal normalization model (D-STAN)
Angus Chapman [angusc@bu.edu]
January 2025

Reference:
"A dynamic spatiotemporal normalization model for continuous vision"
Angus F. Chapman & Rachel N. Denison

Instructions:
The model code is a library of MATLAB functions. Put the directory on your computer and run the code from this directory in MATLAB.

Main code used in the analyses reported in the manuscript are in the home directory, while core model code is located in the 'model' subdirectory.

Output data for the more computationally intensive analyses are available for download on our OSF page [https://osf.io/qy9pa/], and can be analyzed in demoResults.m after being placed in the 'output' subdirectory. (We exclude them from this repository due to their filesizes).

- demoResults.m: produces example analyses for each result, consistent with several figures shown in the manuscript
- testRandomSeq.m: performs reverse correlation analyses with random stimulus sequences (computationally intensive)
- fitRandFunctions.m: summarizes reverse correlation analyses and fits difference of Gammas to sensory responses (computationally intensive)
- testSurroundSuppression.m: performs surround suppression simulation varying contrast of center and surround stimuli
- testSubadditivity.m: analyses subadditive responses (computationally intensive)
- testResponseAdaptation_iden.m: analyses response adaptation for identical stimulus sequences (computationally intensive)
- testResponseAdaptation_noniden.m: analysis response adaptation for non-identical stimulus sequences (varying similiarity and suppressive tuning)
- testBackwardMasking.m: analyses backward masking (computationally intensive)
- testContrastDependentSuppression.m: analyses interactions between T1 and T2 contrast (computationally intensive)
