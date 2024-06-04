function mrf_params(nr_interleaves, ratio)

fa = fisp_fa_schedule(1 / ratio);

fname = sprintf('mrf_params_ni%d_nt%d.csv', nr_interleaves, size(fa, 1));
headers = {'flip_angle_fraction_of_nominal', 'rf_phase_(deg)', 'extension_of_base_TR_ms', 'extension_of_base_TE_ms'};

fa = fa / 75;
phase = (mod((1:size(fa, 1)) - 1, 2) * 180)';
tr = zeros(size(fa));
te = zeros(size(fa));

mrf_params = cat(2, fa, phase, tr, te);
mrf_params = num2cell(mrf_params);
writecell(cat(1, headers, mrf_params), fname);