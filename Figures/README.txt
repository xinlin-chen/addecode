Figures in '* recon performance' folders contain figures displaying reconstruction quality metrics (QM) as different parameters are evaluated

QM used:
RLNE (RLNE)
Peak signal-to-noise ratio (PSNR)
Structural similarity index (SSIM)

Parameters varied:
Acceleration (a): sampling acceleration rate
Calibration region (cr): size of calibration region in sampling pattern
Regularization step-size (rs): step-size during ESPiRIT and SENSE reconstruction
Kernel size (eks): size of ESPiRIT kernel
Number of maps (enm): number of sensitivity maps to generate (for ESPiRIT)

Plot descriptions:
Reconstruction quality for each slice, as well as mean reconstruction quality, were compared. Zero-filled and ESPiRIT reconstructions were compared in all figures. In some figures, SENSE reconstruction was also compared (denoted by 'sense' in the filename).

File naming conventions:
QM values for each slice and evaluated parameter values are saved in recon_qm_* as fig and png images
Mean QM values for each tested value of the parameter are saved in mean_recon_qm* as fig and png images
If SENSE reconstruction was used, 'sense' will then appear in the filename
After the headers listed above, the axis sliced is indicated next (e.g. 'X')
Then, the values of fixed parameters are listed, attached to the parameters' acronyms (encased in parentheses above). For example, if acceleration is varied, the values of cr, rs, eks and enm will be included in the filename.
If more than one range of values was evaluated for a parameter (e.g. step-size, if the first run did not properly demonstrate the effect of lower/higher step-sze), a _<int> is appended (e.g. *_1)

Notes:
Regularization method: in general, l1 was used. Figures demonstrating the effect of step-size for both l1 and l2 were generated. Any plots using l2 reg start with l2_*

True acceleration values were calculated by finding the proportion of non-zero values in the sampling pattern. Different levels of acceleration were obtained by changing the 'nominal' acceleration, i.e. the inputs to -y and -z in the 'poisson' function used by BART. These inputs affected acceleration in the y and z directions (min. 1). True acceleration values were impacted by the 'nominal' acceleration input as well as calibration region size.

Acceleration value when not evaluated: in general, a 'true' acceleration value of 4.09 (a 1, cr 50) was used for 127 (4 coil), and 7.81 (a 4, cr 40) for B04027 (2 coil)

