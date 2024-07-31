This is the code repository for manuscript "Efficient 3D FISP-MRF at 0.55 T using long spiral readouts and concomitant field effect mitigation", Zhibo Zhu, Nam G. Lee, Krishna S. Nayak

The programe was originally excuted on a Hyperplane 2U system with 4 Nvidia A100 GPUs. Excution on a Linux system and GPU(s) is therefore necessary.

# External dependencies:
* *siemens_to_ismrmrd*: See https://github.com/ismrmrd/siemens_to_ismrmrd.
* *ismrmrd*: Included.
* *pulseq*: Included.

# Installation (Linux system required):
* Download this repository.
* Install *siemens_to_ismrmrd* (Linux installation).
* Download data (See Data Availability, total size ~700 GB).

# Folder organization:
* *figure1* to *figure8*, *movieM3*: Figure generation codes.
* *src*: Image reconstruction source codes.\
  Example:
  * Presence of Siemens raw data (```.dat``` format), raw data in ISMRMRD data format (```.h5``` format, see *ismrmrd* section below) , Pulseq file (```.seq``` format, see *pulseq* section below), spiral trajectory file (```.h5``` format, see *pulseq* section below) and dictionaries with B1+ correction (```.mat``` format) are required.
  * Set B0 and B1 mapping DICOM pathes to *../data/volunteer/pulseq_mrf_20240329*.
  * Convert B0 and B1 mapping DICOMs to numerical maps.
  * Perform MRF image recontruction on slice 25 and 36 sequentially with 3 approaches: Plain gridding, subspace modelling + gridding and subsapce modeling + MaxGIRF encoding. All 6 readouts will be reconstructed.
  * Perform MRF dictionary matching for slice 25 and 36.
  ```
  cd ./src/pulseq_mrf_recon;
  b0_path = '../data/volunteer/pulseq_mrf_20240329/b0mapping';
  b1_path = '../data/volunteer/pulseq_mrf_20240329/b1mapping';
  b0map = convert_dicom_to_b0(b0_path, false);
  [b1map, indx] = conver_dicom_to_b1(b1_path, (0.5:0.05:1.5).', false);
  run pulseq_mrf_recon_3d_gpu.m;
  run pulseq_mrf_patternmatch.m;
  ```
  To modify paths:
  * Modify ```b0_path``` and ```b1_path```.
  * Modify ```data_dir``` in ```pulseq_mrf_recon_3d_gpu.m``` and ```pulseq_mrf_patternmatch.m```.

  To perform partial reconstructed on dataset with available readout durations:
  * Set ```opt_indx``` between 1 to 14 (phantom datasets) or within [1 3 5 8 11 14] (volunteer datasets) in ```pulseq_mrf_recon_3d_gpu.m``` and ```pulseq_mrf_patternmatch.m```.
  * Selections on dictionaries and sequence settings are automatic.

  To reconstruct different slices:
  * Change ```slice_loc``` in ```pulseq_mrf_recon_3d_gpu.m``` and ```pulseq_mrf_patternmatch.m```.

* *ismrmrd*: External libraries for ISMRMRD data format.\
  Example: Convert raw Siemens data in *data/volunteer/pulseq_mrf_20240401* to ISMRMRD format.
  ```
  cd ./src/pulseq_mrf_recon;
  % Linux system with siemens_to_ismrmrd installed is required.
  run convert_siemens_to_ismrmrd.m;
  ```
  To modify path:
  * Modify ```data_dir``` in ```convert_siemens_to_ismrmrd.m```.
  
* *pulseq*: Sequence source codes and external Pulseq libraries.\
  Example: Create all Pulseq sequences in ```.seq``` format and save spiral trajectories in ```.h5``` format.
  ```
  cd ./pulseq/lowfield_mrf;
  run make_mrf_sequence.m;
  ```
* *data*: Data storage folder. Suffix: &nbsp;*_iso* or no suffix: Data were acquired at isocenter. &nbsp;*_h75mm*: Data were acquired at Head 75 mm location.\
  *data/dictionaries/disc/pulsel_readout_experiments/*: All dictionary files.\
  *data/phantom/pulseq_mrf_20231113*: The ACR MRI phantom data.\
  *data/phantom/pulseq_mrf_20240306*: The ISMRM/NIST system phantom data.\
  *data/volunteer/*: Two healthy volunteer data.\
  *$subsfolder/b0mapping/*: B0 mapping dicom images.\
  *$subsfolder/b1mapping/*: B1 mapping dicom images.\
  *$subsfolder/subspace_3d/*: Subspace modelling + gridding reconstruction results.\
  *$subsfolder/maxgirf_3d/*: Subspace modelling + MaxGIRF encoding reconstruction results.\
  
  
# Data availability
https://drive.google.com/drive/folders/1UfJkXGt2ioBZ8tzAUWCpVHZqJ27lrU03?usp=sharing

# Contact
Correspondence to: zhibozhu@usc.edu, Zhibo Zhu
