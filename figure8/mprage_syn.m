clear; clc; close all;

ni = [48 22 14 10 8 6];
nt = [1036 709 535 421 352 278];
ii = 1;
mat_name = sprintf('./maxgirf_3d/pulseq_mrf_maxgirf_slice36_B0__fa75_cycle2_ni%d_nt%d.mat', ni(ii), nt(ii));
load(mat_name, 'm0map', 't1map', 't2map');
load('./mask_slice36.mat');
mask = double(mask);

% Set up relaxation operators.
E1 = @(t, T1) exp(-t ./ T1);
E2 = @(t, T2) exp(-t ./ T2);
E1i = @(t, PD, T1) -PD .* (E1(t, T1) - (1 - E1(t, T1)));
fa = 20;
tr = 6.1;
te = 3.05;

% Remove csf.
info1 = dicominfo('./mprage_brain.IMA');
ti = 3000;
mprage = @(PD, T1, T2) E1i(ti, PD, T1) * sind(fa) .* (1 - E1(tr, T1)) ./ (1 - E1(tr, T1) * cosd(fa)); %  .* E2(te, T2);
brain = abs(mprage(m0map, t1map, t2map));
brain = mask .* brain;

brain = brain * 255 / max(brain(:));
brain = uint8(brain);

figure;
imagesc(abs(brain).');
axis image;
colormap gray;

info1.LargestImagePixelValue = 255;
info1.WindowCenter = 127;
info1.WindowWidth = 256;

dicomwrite(abs(brain).', 'mprage_brain.IMA', info1);

% Keep GM.
% info2 = dicominfo('./mprage/VOL836.MR.MRF_PULSEQ.0006.0048.2024.03.29.09.07.38.976261.175965267.IMA');
% ti = 350;
% mprage = @(PD, T1, T2) E1i(ti, PD, T1) * sind(fa) .* (1 - E1(tr, T1)) ./ (1 - E1(tr, T1) * cosd(fa)); %  .* E2(te, T2);
% gm = mprage(m0map, t1map, t2map);
% gm = mask .* gm;
% 
% figure;
% imagesc(abs(gm).');
% axis image;
% colormap gray;
% 
% gm = gm * 255 / max(gm(:));
% gm = uint8(gm);
% 
% info2.LargestImagePixelValue = 255;
% info2.WindowCenter = 127;
% info2.WindowWidth = 256;
% 
% dicomwrite(abs(gm).', 'mprage_gm.IMA', info2);

