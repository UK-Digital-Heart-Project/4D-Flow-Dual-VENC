# 4D Flow Dual VENC

Analysis pipeline for dual-Venc 4D-Flow cardiac MRI.

## Prerequisites

- MATLAB, with at least the Image Processing Toolbox installed.
- Windows, Linux or MacOS - there should not be any platform dependencies.  

## Introduction

Merged dual-Venc flow maps are created from separately acquired low-Venc and High-Venc phase-contrast MR scans.

## Input data

A. A high and a low-Venc modulus cine-stack.

B. A high and a low-Venc PC cine-stack.

Modulus and phase images acquired at the same Venc are simultaneous, but not necessarily images acquired at the low and high Vencs -
i.e., the acquisitions at different Vencs may be interleaved or sequential.

A complete analysis involves merging low and high-Venc data sets acquired with 3 velocity-encoding directions along orthogonal axes
(typically referred to the patient orientation rather than the LAB frame).

Velocities (in units of cm/s) are calculated from the input PC cine-stacks according to the equation:

Velocity = Intercept + Slope.Grayscale, 

where Intercept = - LowVenc and Slope = LowVenc/(2^11).

## Method

1. At each cardiac phase (epoch), the low-Venc modulus volume is co-registered to the high-Venc volume.
   First, a manual shift is applied in all 3 directions, to allow for gross subject motion, and the amount of shift recorded.   
   Then, a non-rigid co-registration is performed, and the displacement field saved for re-use.
   
2. The manual shift correction and the displacement field saved in step (1) are now used to co-register the corresponding epoch of the PC 
   low-Venc image-stack to its matching epoch in the high-Venc acquisition. 
   
   All volumes used in phases (1, 2) are interpolated to cubical voxels in the slice direction before the manual and automated 
   co-registration steps, and down-sampled back to the original z-locations immediately afterwards.
   
3. A *wrapping correction* is now calculated from the original high-Venc and the motion-corrected low-Venc cine-stacks.

4. This is rounded to the nearet multiple of 2 x LowVenc, and applied to the low-Venc cine-stack to create a *merged* result.
   Voxels with a wrapping correction in excess of 8 x LowVenc are regarded as background noise and left uncorrected.
   
5. The merged output is stored in a series of DICOM files at 16-bits precision (whereas the inputs were 12-bits).
   This has the *precision* (velocity resolution) of the low-Venc acquisition, with the *dynamic range* of the high-Venc acquisition,
   with any wrapping corrected *using a measurement rather than an estimate.* The velocity-to-noise ratio is also better than that of the 
   original high-Venc acquisition.

The final PC cine-stacks are generated according to the equation:

Grayscale = uint16((2^15).(Velocity/Venc + 1.0)), 

where Venc is the velocity limit appropriate to the cine-stack.

The Venc is written (explicitly or implicitly) to the output  header in at least 3 places:

- Header.csa.FlowVenc.
- Header.RescaleSlope and Header.RescaleIntercept.
- Header.SequenceName, e.g., "FL200" for Venc = 200 cm/s.

NOTE THAT, FOR _OUTPUT_ IMAGES, THE _FlowVenc_ MAY NOT BE WRITTEN RELIABLY, AND YOU SHOULD NOT DEPEND ON IT;
USE THE _RescaleSlope_ and _RescaleIntercept_ INSTEAD.
YOU CAN ALSO PARSE THE _SequenceName_ BY STRIPPING ANY NON-NUMERIC CHARACTERS AND CONVERTING THE REMAINING STRING TO A NUMBER.

NOTE, ALSO, THAT UNPREOCESSED _INPUT_ IMAGES NEED TO BE TREATED DIFFERENTLY,
SINCE THE _RescaleSlope_ and _RescaleIntercept_ ARE SET TO 2 AND - 4096, RESPECTIVELY, AND THESE ARE MEANINGLESS.
INSTEAD OF THESE, YOU SHOULD FETCH THE _FlowVenc_ USING THE FUNCTION _pft_FetchVencFromHeader_, AND USE THE CONVERSION:

Velocity = Intercept + Slope.Grayscale, 

where Intercept = - Venc,

Slope = 2.0*Venc/double(2^BS),

and BS is read from the BitsStored field (either 12 or 16) of the Dicom header.

## Outputs
   
The workflow creates an audit trail of several intermediate results, as well as the final merged PC cine-stack and a short text-mode summary.

a. The motion-corrected low-Venc modulus cine-stack.

b. The motion-corrected low-Venc PC cine-stack.

c. The manual shift corrections (as text files) and the displacement fields (as MAT files).

d. Screenshots of each epoch of the modulus and PC acquisition, before and after co-registration, in the format of a volume montage. 

e. The merged PC cine-stack.

f. The difference between the original high-Venc PC cine-stack and the final merged PC cine-stack - i.e., a "discrepancy" measure.

g. A simple text-mode summary.

Note that the input cine-stacks are sorted according to slice location and trigger time, and that the output cine-stacks are written
in the order [ Row, Column, Epoch, Slice ], with the last index varying most slowly according to the MATLAB convention.
This may mean that, whereas the input data show entire volumes at a given epoch, then step to the next epoch, the outputs show complete cines
of each slice, before advancing to the next slice.

## Installation

Clone this repo to a folder in your MATLAB workspace, then add all the directories to the path:

```addpath(genpath(pwd)); savepath;```

## Usage

Run the single main script at the outer level.

Follow the prompts to nominate source and target directories.

During the co-registration step, a GUI with sliders will prompt for manual shift corrections in 3 directions to align a low-Venc modolus epoch
to its high-Venc counterpart. A push-button applies the change. Corrections are carried over from one epoch to the next, so it may be sufficient
to make a careful correction for the first epoch, then hit "Apply" for the rest. In any case, the manually corrected low-Venc stack is merely the
starting point for the automated non-rigid co-registration.

## Test data

Anonymized DICOM data are available on request from the author:

ptokarcz@ic.ac.uk

## Citation

TBC.
