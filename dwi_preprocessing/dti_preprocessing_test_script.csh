#!/bin/csh

# directory
cd /mnt/d/MRI/Analysis/example_dti_brooke_2

## STEP(1): Denoising
#  see http://mrtrix.readthedocs.io/en/latest/dwi_preprocessing/denoising.html)

# Step(1a) The first step of the image-processing pipeline is image denoising:
echo 'Performing dwi_pa denoising using Mrtrix dwidenoise'
# dwidenoise dwi_pa.nii.gz -extent 3,3,3 dwi_pa_noise_corrected.nii.gz -noise dwi_pa_full_noise_map.nii.gz -info -force
dwidenoise dwi_pa.nii.gz -extent 3,3,3 dwi_pa_noise_corrected.nii.gz -noise dwi_pa_full_noise_map.nii.gz -info -force
# view
fslview_deprecated dwi_pa.nii.gz dwi_pa_noise_corrected.nii.gz dwi_pa_full_noise_map.nii

# Steo(1b) Check: What is the difference/are the residuals once you have removed the noise? (i.e. out-in, shows noise distribution)
mrcalc dwi_pa.nii.gz dwi_pa_noise_corrected.nii.gz -subtract residual_full_PA.nii -info -force
# mrview residual_full_PA.nii)
# view
fslview_deprecated residual_full_PA.nii.gz

# Step:Computation of Rician Noise
# This step isn't working at the moment
# matlab -nosplash -nodesktop -r dwi_extractor_pa <====
dwidenoise dwi_b2000_pa.nii.gz -extent 3,3,3 dwi_b2000_pa_nc.nii.gz -noise dwi_b2000_pa_noise_map.nii.gz -info -force
mrcalc dwi_pa_noise_corrected.nii.gz 2 -pow dwi_b2000_pa_noise_map.nii.gz 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if dwi_pa_noise_corrected.nii.gz -force

# STEP(2): Gibbs ringing
echo 'Gibbs Correction'
# DOUBLE CHECK THE AP / PA IS A TYPO OR OUTPUT OF dwi_extractor_pa
 mrdegibbs dwi_pa_noise_corrected.nii.gz dwi_pa_noise_gibbs_corrected.nii.gz -info -force

# STEP(3): Eddy current correction
# The following correction only needs to be performed when an Eddy current induced artifact is visible in the data. To detect this artifact, compare the DW images to the b0 images (no eddy currents in b0 images) and thus look for stretching or motion-like artifacts.
# Will see Eddy currents in some artery DWI data sets– Salman performed this correction number of times.

# When you perform this correction => visually check the result! If a geometric transformation occurs (stretching of the images) due to overcompensation by this correction (so in case there are in fact no eddy current induced artifacts), don’t perform this correction.
# The output includes std maps, which will tell you for which slice the correction was performed.
echo 'Eddy Current Correction'

# Generate index file
# Assumes all volumes acquired with same phase encoding direction
# This makes a whole bunch of assumptions at the moment and needs to be customised
# See https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#A--index
make_index_for_eddy.csh dwi_pa_noise_gibbs_corrected.nii.gz

# Generate acqparams.txt file
# Again this makes a whole bunch of assumptions at the moment and needs to be customised
make_acqp_for_eddy.csh dwi_pa_noise_gibbs_corrected.nii.gz


# eddy_openmp --imain=dwi_pa_noise_gibbs_corrected \
#             --mask=dwi_pa-mask \
#             --acqp=acqparams.txt \
#             --index=index_PA.txt \
#             --bvecs=bvecs_PA.bvec \
#             --bvals=bval_PA.bval \
#             --topup=topup_PA_AP_b0 \ # <= We haven't run topup do we need to do this?
#             --fwhm=0 \
#             --flm=quadratic \
#             --repol \
#             --ol_type=sw \
#             --out=eddy_corrected_data \
#             --data_is_shelled \ # <= this is a flag to turn off internal checks ... should make sure this is appropriate
#             --verbose

eddy_openmp --imain=dwi_pa_noise_gibbs_corrected \
            --mask=dwi_pa_mask \
            --acqp=acqparams.txt \
            --index=index.txt \
            --bvecs=$2 \
            --bvals=$3 \
            --fwhm=0 \
            --flm=quadratic \
            --repol \
            --ol_type=sw \
            --out=eddy_corrected_data \
            --data_is_shelled \
            --verbose

# # STEP(4) DW bias correction
# # “Bias field signal is a low-frequency and very smooth signal that corrupts MRI images.”
# # Inhomogeneity correction done in FSL.
# echo 'DW bias correction'
# dwibiascorrect -fsl eddy_corrected_data.nii.gz -fslgrad eddy_corrected_data.eddy_rotated_bvecs bval_PA.bval dwi_b1_bc.nii
