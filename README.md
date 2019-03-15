# 4D Flow Dual VENC
Analysis pipeline for dual VENC 4D Flow CMR

Merged dual-Venc flow maps are created from separately acquired low-Venc and High-Venc phase-contrast MR scans.

Inputs are:

A. A high and a low-Venc modulus cine-stack.

B. A high and a low-Venc PC cine-stack.

Modulus and phase images acquired at the same Venc are simultaneous, but not necessarily images acquired at the low and high Vencs -
i.e., the acquisitions at different Vencs may be interleaved or sequential.

The analysis workflow is as follows:

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
   
The workflow creates an audit trail of several intermediate results, as well as the final merged PC cine-stack and a short text-mode summary.

a. The motion-corrected low-Venc modulus cine-stack.

b. The motion-corrected low-Venc PC cine-stack.

c. The manual shift corrections (as text files) and the displacement fields (as MAT files).

d. Screenshots of each epoch of the modulus and PC acquisition, before and after co-registration, in the format of a volume montage. 

e. The merged PC cine-stack.

f. The difference between the original high-Venc PC image stack and the final merged PC cine-stack - i.e., a "discrepancy" measure.

g. A simple text-mode summary.
