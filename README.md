# G.729 in Matlab

This project started as an experiment to learn about a basic CELP vocoder, by implementing the ITU standard in Matlab.  Standardized codecs have fixed and/or floating point reference code, but it is not always easy to learn from.

## Goals
1. Implement floating point model in Matlab, WITHOUT using ITU reference C-code (except for VQ tables)
2. Implement encoder, decoder, and basic bitstream syntax for G.729 only (no extesions, e.g. G.729A, G.729AB)

## Notes
1. Implementation (and Equation references) directly based on specification "ITU-T Rec. G.729 (06/2012) Coding of speech at 8 kbit/s using conjugate-structure algebraic-code-excited linear prediction (CS-ACELP)"  
2. Additional background material on CELP codecs and G.729 obtained from [Wai C Chu's book](http://www.amazon.com/Speech-Coding-Algorithms-Foundation-Standardized/dp/0471373125) 

## Matlab code notes
1. Many loops purposely written inefficiently in Matlab (C-style for loops, to make conversion easier) 
2. Some code using ZI notation for buffers that contain current frame data and past history.  ZI references the Matlab buffer index corresponding to time = 0.  Also personal style preference to aid in converting to C code (with corresponding C code pointer offsets)
3. Some details may be incorrect!  It was a learning experiment after all.  :)  
