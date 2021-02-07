[![View mni2atlas on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/87047-mni2atlas)
# mni2atlas 
FSL Anatomical Labels for MNI Vector/ROI. Mni2atlas takes an ROI or a vector of coordinates (both in the MNI
space) and returns labels from different FSL atlases.

In **VECTOR** modality labels are returned in probability values (same results of FSL atlas tool).

In **ROI** modality the probability value reported for a label represents the frequency of that label in the roi for a given threshold of the FSL atlas
probability map (0, 25 or 50; default = 25).

## HOW TO USE
   mni2atlas(VECTOR/ROI) the first input can be a MNI vector or an ROI in the MNI space. Depending on the input the script switches between two
   different work modalities. With no other input the script will seek labels among all available fsl atalses.

   mni2atlas(VECTOR/ROI,ATLAS_SELECTOR) allows to choose among the following atlases:
   1. Juelich Histological Atlas
   2. Harvard-Oxford Cortical Structural Atlas
   3. Harvard-Oxford Subcortical Structural Atlas
   4. JHU ICBM-DTI-81 White Matter labels
   5. JHU White Matter tractography Atlas
   6. Oxford Thalamic Connectivity Atlas
   7. Cerebellar Atlas in MNI152 after FLIRT
   8. Cerebellar Atlas in MNI152 after FNIRT
   9. MNI Structural Atlas
   
   ATLAS_SELECTOR must be a row vector (i.e. [1,3,6]). Default value is [1:1:9]. You can also leave it as an empty vector (i.e. (VECTOR/ROI,[])).

   [ATLAS]=MNI2ATLAS(VECTOR/ROI,...) the script returns the structure ATLAS whit the following fields: .name (of the atlas), .labels (a cell vector). No stdout will be print.

   mni2atlas(VECTOR) prints on screen labels found for the MNI VECTOR position.

   mni2atlas(ROI) prints on screen labels found for the input ROI. ROI can be a preloaded (with load_nii) volume or the path of a nifti volume. 

### ADVANCED OPTIONS
   mni2atlas(ROI,ATLAS_SELECTOR,THR) THR allows to choose among 3 threshold levels: 0, 25, 50 (i.e,. 0, 25, 50). Default value is 25. Option available only under ROI modality.

   mni2atlas(VECTOR,ATLAS_SELECTOR,RESOLUTION) RESOLUTION allows to choose between ‘1mm’ or ‘2mm’ atlases. 1mm atlases performs better region identification but requires more loading time. Default value is ‘1mm’. Option available only under VECTOR modality.

## SYSTEM REQUIREMENTS
  [NifTI and ANALYZE tool](https://it.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) (version > 2012-10-12) 

## ACKNOWLEDGEMENTS
  This function uses some of the available [FSL atlases](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases)
