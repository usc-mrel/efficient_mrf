This is the code repository for manuscript "Efficient 3D FISP-MRF at 0.55 T using long spiral readouts and concomitant field effect mitigation", Zhibo Zhu, Nam G. Lee, Krishna S. Nayak

The programe was originally excuted on a Hyperplane 2U system with 4 Nvidia A100 GPUs. Excution on a Linux system and GPU(s) is therefore necessary.

# Folder organization:
* *figure1* to *figure8*, *moviceM3*: Figure generation codes.
* *src*: Image reconstruction source codes.
  ```
  %%% Modify path accordingly.
  cd ./src/pulseq_mrf_recon;
  b0map = convert_dicom_to_b0(b0_path, false);
  [b1map, indx] = conver_dicom_to_b1(b1_path, (0.5:0.05:1.5).', false);
  run pulseq_mrf_recon_3d_gpu.m; % Modify path accordingly inside the script.
  run pulseq_mrf_patternmatch.m; % Modify path accordingly inside the script.
  ```
* *ismrmrd*: External libraries for ISMRMRD data format.
  ```
  cd ./src/pulseq_mrf_recon;
  % Linux system with siemens_to_ismrmrd installed is required.
  run convert_siemens_to_ismrmrd.m; % Modify path accordinglu inside the script.
  ```
* *pulseq*: Sequence source codes and external Pulseq libraries.
  ```
  cd ./pulseq/lowfield_mrf;
  run make_mrf_sequence.m;
  ```
* *data*: Data storage folder.
  
# Data availability
https://drive.google.com/drive/folders/1UfJkXGt2ioBZ8tzAUWCpVHZqJ27lrU03?usp=sharing

# Contact
Correspondence to: zhibozhu@usc.edu, Zhibo Zhu
