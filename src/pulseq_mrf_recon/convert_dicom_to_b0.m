function b0map = convert_dicom_to_b0(b0_path, debug)

dicom_list = dir((fullfile(b0_path, '*.IMA')));
for ii = 1:length(dicom_list)
    dicom_name = fullfile(dicom_list(ii).folder, dicom_list(ii).name);
    img = dicomread(dicom_name);
    info = dicominfo(dicom_name);
    
    if ii == 1
        [Nx, Ny] = size(img);
        b0map = zeros(Nx, Ny, length(dicom_list));
        slice_loc = zeros(length(dicom_list), 1);
    end
    
    b0map(:, :, ii) = double(img);
    slice_loc(ii) = info.SliceLocation;
end

[slice_loc, indx] = sort(slice_loc);
b0map = b0map(:, :, indx);

b0map = (b0map - 2048) * 2;
b0map = b0map / 2482; % [ppm]
b0map = b0map * 42.57; % [Hz]
b0map = permute(b0map, [2 1 3]);

if debug
    figure();
    imagesc3D(b0map);
end

save(fullfile(b0_path, 'b0map.mat'), 'b0map');