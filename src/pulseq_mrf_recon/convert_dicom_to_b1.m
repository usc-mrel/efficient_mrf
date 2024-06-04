function [b1map, b1_indx] = convert_dicom_to_b1(b1_path, b1_range, debug)

dicom_list = dir((fullfile(b1_path, '*.IMA')));
for ii = 1:length(dicom_list)
    dicom_name = fullfile(dicom_list(ii).folder, dicom_list(ii).name);
    img = dicomread(dicom_name);
    info = dicominfo(dicom_name);
    
    if ii == 1
        [Nx, Ny] = size(img);
        b1 = zeros(Nx, Ny, length(dicom_list));
        slice_loc = zeros(length(dicom_list), 1);
    end
    
    b1(:, :, ii) = double(img) / 800;
    slice_loc(ii) = info.SliceLocation;
end

[slice_loc, indx] = sort(slice_loc);
b1 = b1(:, :, indx);
b1map = imresize3d(b1, [256 256 48]) * 0.92;

if debug
    figure();
    imagesc3D(b1map);
end

%% Round to the nearest 0.05.
% b1_range = (0.5:0.05:1.5).';
b1_template = repmat(permute(b1_range, [2 3 4 1]), [256 256 48 1]);
b1_temp = abs(repmat(b1map, [1 1 1 length(b1_range)]) - b1_template);
[~, b1_indx] = min(b1_temp, [], 4);

if debug
    figure();
    imagesc3D(b1_indx);
end

save(fullfile(b1_path, 'b1map.mat'), 'b1map', 'b1_indx');