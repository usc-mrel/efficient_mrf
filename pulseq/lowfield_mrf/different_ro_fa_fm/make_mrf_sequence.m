% demo_pulseq_fisp_mrf_2d_spiral.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 03/17/2022, Last modified: 12/18/2022
% Spiral FISP modification, Zhibo.
% Add averages, Zhibo, 02/28/2023.
% Add adjustable number of interleaves, oversampling and averages, Zhibo,
% 03/02/2023.
% Add optional rewinder design, Zhibo, 03/06/ 2023.
% Add MRF set up, including under-sampling (single repetition),
% fully-sampling (all repetitions) and averages, Zhibo, 03/06/2023.
% Converted into an invokable function, Zhibo, 04/06/2023.
% Re-converted into a stand-alone script, 06/29/2023.
% Major recoding, started: 06/29/2023.

function make_mrf_sequence(nr_interleaves, fn_mrf_seq_path, grad_mode)
%% Define MRI scanner specifications
B0 = 0.55; % main field strength [T]
% grad_mode = 'fast';

switch grad_mode
    case 'fast'
        max_grad = 24;      % Max gradient strength [mT/m]
        max_slew = 180.18;  % Maximum slew rate [mT/m/ms]
        opt_data = [];
    case 'normal'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 100;     % Maximum slew rate [mT/m/ms]
    case 'whisper'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 50;      % Maximum slew rate [mT/m/ms]
    case 'freemax'          % A rough estimation of the derated Freemax, e.g., Gmax is divived by sqrt(3).
        max_grad = 14;
        max_slew = 40.54;
        opt_data = 'freemax';
end

%% Set system limits
sys = mr.opts('MaxGrad'       , max_grad, 'GradUnit', 'mT/m' , ...
    'MaxSlew'       , max_slew, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6 , ...
    'rfDeadTime'    , 100e-6, ...
    'adcDeadTime'   , 10e-6 , ...
    'B0'            , B0);

%% Define imaging parameters
%--------------------------------------------------------------------------
% Parameters for a spiral trajectory
%--------------------------------------------------------------------------
flip_angle      = 75;        % Flip angle [deg] (reference flip angle)
fov_read        = 300e-3;    % field of view [m]
base_resolution = 256;       % Base resolution
slice_thickness = 5e-3;      % slice thickness [m]

%--------------------------------------------------------------------------
% Parameters for a windowed sinc pulse (excitation)
%--------------------------------------------------------------------------
rf_length      = 2e-3;       % RF length [sec]
rf_apodization = 0.42;       % RF apodization
rf_tbw         = 8;          % RF time bandwidth product [sec*Hz]
rf_phase       = 180;        % RF phase [deg]

%--------------------------------------------------------------------------
% Parameters for an adiabatic inversion pulse
%--------------------------------------------------------------------------
T    = 20e-3;                % RF duration [sec]
A0   = 15;                   % maximum RF amplitude [uT]
beta = 800;                  % modulation angular frequency [rad/sec]
mu   = 4.9;                  % dimensionless parameter

%--------------------------------------------------------------------------
% Misc. parameters
%--------------------------------------------------------------------------
discard_pre       = 20;      % number of samples to be discarded at the beginning of acquisition
discard_post      = 20;      % number of samples to be discarded at the end of acquisition
% bandwidth         = 390;     % Readout bandwidth [Hz/px]
real_dwell_time   = 2.5e-6;  % sampling interval for k-space acquisition [sec]
dephasing         = 2;       % dephasing in cycles across the slice thickness [cycle]
fov_phase         = 100;     % FoV phase [%]

%% Calculate the real dwell time [sec]
%--------------------------------------------------------------------------
% real dwell time [sec]
% IDEA p219: dRealDwellTime denotes the dwell time with oversampling
%--------------------------------------------------------------------------
% round-down dwell time to 100 ns (sys.adcRasterTime  = 100 ns)
% real_dwell_time = round((1 / bandwidth) / (readout_os_factor * base_resolution) * 1e7) * 1e-7;

%% Define the relative path of a .csv file
%--------------------------------------------------------------------------
% A list of MRF setting, 1036x1.
% Column 1: FA ratio, actual_fa / 75.
% Column 2: RF phase, [0 180].
% Column 3: TR variation.
% Column 4: TE variation.
% fn_mrf_seq_path = 'mrf_params_ni48_nt1036.csv';

%% Define the number of time frames
nr_averages = 1;
nr_repetitions = 1;
% nr_interleaves = 48;
osf = 1;
% nr_time_frames = nr_interleaves * osf;

%% Read a .csv file
mrf_struct = importdata(fn_mrf_seq_path);

%--------------------------------------------------------------------------
% FA series [unitless]
%--------------------------------------------------------------------------
FA_series = mrf_struct.data(:, 1); % [unitless]
FA_series(FA_series == 0) = 0.001; % pure zero gives an error

%--------------------------------------------------------------------------
% TR series [sec]
%--------------------------------------------------------------------------
TR_base = 10e-3; % [sec]
% Set a constant TR for simplicity, Zhibo.
% TR_base is not in use due to a minimum TR set up.
TR_series = mrf_struct.data(:, 3) + TR_base;
TR_series = round(TR_series ./ sys.gradRasterTime) * sys.gradRasterTime;

%--------------------------------------------------------------------------
% TE series [msec]
%--------------------------------------------------------------------------
TE_series = mrf_struct.data(:, 4);

nr_time_frames = length(FA_series);

%% Define a series of rotation angles [rad]
phi_series = 360 / (nr_interleaves * osf) * pi / 180 * (0:nr_time_frames-1).'; % [rad]

%% Create a sequence object
seq = mr.Sequence(sys);
start_time = tic;

%% FISP-MRF pulse sequence
%--------------------------------------------------------------------------------------------------------
%                                                         TR
%    |         |<-------------------------------------------------------------------------------------->|
%    |         |            TE_series                                   TR_series                       |
%    |         |<------------------------------>|<----------------------------------------------------->|
%    |         |                           |    |                                  |    |               |
%    | |       |rf_ex  |                   |    |                                  |    |               |
%    | |      _|_      |                   |    |                                  |    |               |
% RF | |     / | \     |                   |    |                                  |    |               |
% ___|_|    /  |  \    |___________________|____|__________________________________|____|_______________|
%    | |\__/   |   \__/|                   |    |                                  |    |               |
%    | |       |       |                   |    |                                  |    |               |
%    | |               | |                 |    |                                  |    |   gz_spoil    |
%    | |               | |                 |    |                                  |    |   ______      |
%    | |_______|_______| |                 |    |                                  |    |  /      \     |
% Gz | /       | gz_ex \ |                 |    |                                  |    | /        \    |
% ___|/        |        \|          |______|____|__________________________________|____|/          \___|
%    |         |         \         /|      |    |                                  |    |               |
%    |         |         |\_______/ |      |    |                                  |    |               |
%    |         |         | gz_reph  |      |    |                                  |    |               |
%    |<----------------->|<-------->|      |    |                                  |    |               |
%    |  mr.calcDuration  |          |      |    |                                  |    |               |
%    |      (gz_ex)      |<--------------->|<-->| delay_pre             delay_post |<-->|<------------->|
%    |                   |     delayTE     |    |                       _          |    |    delayTR    |
%    |                   |                 | |  |          _           / \         |  | |               |
% Gx |                   |                 | |  | _       / \         /   \        |  | |               |
% ___|___________________|_________________|_|__|/ \     /   \       /     \       |__|_|_______________|
%    |                   |                 | |  |   \   /     \     /       \     /|  | |               |
%    |                   |                 | |  |    \_/       \   /         \___/ |  | |               |
%    |                   |                 | |  |               \_/                |  | |               |
%    |                   |                 | |  |                                  |  | |               |
%    |                   |                 | |  |                _                 |  | |               |
%    |                   |                 | |  |     _         / \           ___  |  | |               |
% Gy |                   |                 | |  |    / \       /   \         /   \ |  | |               |
% ___|___________________|_________________|_|__|   /   \     /     \       /     \|__|_|_______________|
%    |                   |                 | |  |\_/     \   /       \     /       |  | |               |
%    |                   |                 | |  |         \_/         \   /        |  | |               |
%    |                   |                 | |  |                      \_/         |  | |               |
%    |                   |   adcDeadTime ->| |<-|<-------------------------------->|->| |<- adcDeadTime |
%    |                   |    = 10 usec    | |  |           grad_duration          |  | |               |
%    |                   |                 | |xx|xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx| |               |
%    |                   |     discard_pre ->|  |<-              ADC             ->|  |<- discard_post  |
%    |                   |                 |                                            |               |
%    |<----------------->|<--------------->|<------------------------------------------>|<------------->|
%    |      block 1      |     block 2     |                  block 3                   |    block 4    |
%----o-------------------+-----------------+--------------------------------------------|-------------> t
%--------------------------------------------------------------------------------------------------------
% Assume TE is a multiple of 10 usec.
% delay_pre : delay_pre  = round((sys.adcDeadTime + adc_duration_pre)  / sys.gradRasterTime) * sys.gradRasterTime
% delay_post: delay_post = round((sys.adcDeadTime + adc_duration_post) / sys.gradRasterTime) * sys.gradRasterTime
% minTE     : minTE = gz_ex.flatTime / 2 + gz_ex.fallTime + mr.calcDuration(gz_reph)
% delayTE   : TE_series = gz_ex.flatTime / 2 + gz_ex.fallTime + delayTE + delay_pre
% delayTR   : delayTR = TR_series - (grad_duration + delay_post)

%% Calculate an adiabatic inversion pulse
%--------------------------------------------------------------------------
% Calculate the bandwidth of a hyperbolic secant pulse
%--------------------------------------------------------------------------
% [dimensionless] * [rad/sec] / [2pi rad/cycle] => [Hz]
BW = 2 * mu * beta / (2 * pi); % [Hz]

%--------------------------------------------------------------------------
% Calculate a time vector [sec]
%--------------------------------------------------------------------------
Nt = floor(T / sys.rfRasterTime);
t = ((-floor(Nt/2):ceil(Nt/2)-1).' + 0.5) * sys.rfRasterTime; % [sec]

%--------------------------------------------------------------------------
% Calculate the amplitude modulation A(t)
%--------------------------------------------------------------------------
am_shape = A0 * sech(beta * t); % [uT]

%--------------------------------------------------------------------------
% Calculate the phase modulation phi(t)
%--------------------------------------------------------------------------
pm_shape = mu * log(sech(beta * t)) + mu * log(A0);

%--------------------------------------------------------------------------
% Calculate a complex-valued hyperbolic secant pulse
%--------------------------------------------------------------------------
rf_shape = am_shape .* exp(1j * pm_shape); % [uT]

%% Create an adiabatic inversion pulse [Hz] and corresponding a slice-selection gradient [Hz/m]
% [Hz/T] * [T/1e6uT] * [uT] * [sec] * [2*pi rad/cycle] => [rad]
nr_slices = 48;
fov_z = nr_slices * slice_thickness;
flip = (sys.gamma * 1e-6 * 2 * pi) * abs(sum(rf_shape)) * sys.rfRasterTime; % [rad]
[rf_inv, gz_inv] = mr.makeArbitraryRf(rf_shape, flip, 'bandwidth', BW, 'sliceThickness', fov_z, 'system', sys);

%% Create a spoiler in the slice direction
% The parameters are obtained from the POET simulation of a (CWRU's 3D) FISP-MRF sequence.
sf2 = 9/20;
amplitude = (sys.gamma * 1e-3) * 8 * sf2; % [Hz/T] * [T/1e3mT] * [mT/m] => *1e-3 [Hz/m]
rise_time = 800 * 1e-6; % [sec]
flat_time = (9000 / sf2 - 800) * 1e-6; % [sec]
fall_time = 600 * 1e-6; % [sec]
gz_inv_spoil = mr.makeTrapezoid('z', 'amplitude', amplitude, 'riseTime', rise_time, 'flatTime', flat_time, 'fallTime', fall_time, 'system', sys);

%% Create a 90 degree slice selection pulse [Hz] and corresponding gradients [Hz/m]
[rf_ex, gz_ex, gz_reph] = mr.makeSincPulse(flip_angle * pi / 180, 'Duration', rf_length, 'SliceThickness', fov_z, 'apodization', rf_apodization, 'timeBwProduct', rf_tbw, 'system', sys);
% This a hardcodes scaling factor controlling Gz rephaser duration. Its
% value is 2.4 to match the default triangle Gz encoding waveform
% durations.
switch grad_mode
    case 'fast'
        sf = 2.4;
    case 'freemax'
        sf = 2.5;
end
% sf = 2.4;
gz_reph.riseTime = gz_reph.riseTime * sf;
gz_reph.fallTime = gz_reph.fallTime * sf;
gz_reph.amplitude = gz_reph.amplitude / sf;
ref_signal = rf_ex.signal(:);

% This the baseline Gz encoding waveform, functioning as a triangle
% waveform to achieve -kz_max encoding.
% The following other Gz encoding are scaled based on it.
deltak_phase = 1 / fov_z;
pe_steps = (-floor(nr_slices / 2):ceil(nr_slices / 2) - 1).' * deltak_phase;
gz_phase_encoding = mr.makeTrapezoid('z', 'Area', pe_steps(1), 'System', sys);

assert(abs(mr.calcDuration(gz_reph) - mr.calcDuration(gz_phase_encoding)) < eps);

%% Calculate a base variable-density spiral trajectory: k in [cycle/cm], g in [G/cm]
% nr_interleaves = 48;                         % number of interleaves
resolution     = fov_read / base_resolution; % resolution [m]
raster_time    = sys.gradRasterTime / 10;    % "raster time" between samples in sec, i.e., 0.000004 gives 4 usec gradient update

%--------------------------------------------------------------------------
% The next 4 variables are for variable density spirals
% They create a transition in the radial spacing as the kspace radius goes
% from 0 to 1, i.e.,
%     0 < kr < us_0, spacing = nyquist distance
%  us_0 < kr < us_1, spacing increases to us_r (affected by ustype)
%  us_1 < kr < 1   , spacing = us_r
% SET us_0 = 1 or us_r = 1 to avoid variable sampling
%
% rate of change in undersampling (see code below)
% ustype: 0 = linearly changing undersampling
%         1 = quadratically changing undersampling
%         2 = Hanning-type chnage in undersampling
%--------------------------------------------------------------------------
us_0 = 0.5;
us_1 = 1;
us_r = 2;
ustype = 0;

%--------------------------------------------------------------------------
% Some more variables determining the type of waveform
% For gtype: 0 = calculate through readout
%            1 = include grad rampdown
%            2 = include rewinder to end at k=0
%--------------------------------------------------------------------------
gtype  = 0;
sptype = 0; % 0 = Archimedean, 1 = Fermat

%--------------------------------------------------------------------------
% Set spparams array
%--------------------------------------------------------------------------
spparams = [sys.gamma * 1e-4; % [Hz/T] * [T/1e4G] => *1e-4 [Hz/G]
    fov_read * 1e2;   % [m] * [1e2cm/m] => *1e2 [cm]
    resolution * 1e2; % [m] * [1e2cm/m] => *1e2 [cm]
    max_slew * 1e2;  % [T/m/sec] * [1e4G/T] * [m/1e2cm] => *1e2 [G/cm/sec]
    max_grad * 1e-1; % [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    max_grad * 1e-1; % [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    raster_time;
    nr_interleaves;
    us_0;
    us_1;
    us_r;
    ustype;
    gtype;
    sptype];

%--------------------------------------------------------------------------
% This function takes parameters passed in spparams array and
% returns a single spiral arm calculated numerically.
% The k-space trajectory for the arm is in karray
% The corresponding gradient waveform is in garray
% The parameter garrlen are the lengths of karray and garray
% The parameter karrlen indicates when sampling stops
%
% FINAL ARRAYS are {kx(0), ky(0), kx(1), ky(1), ...} and
%                  {gx(0), gy(0), gx(1), gy(1), ...}
%--------------------------------------------------------------------------
[karray, garray, karrlen, garrlen] = spiralgen_jgp_11apr_mex(spparams);
if gtype
    % Gradient coordinate system: PRS = [PE,RO,SL] = [gu,gv,gw]
    gv = garray(1:2:garrlen*2) * 1e1; % RO [G/cm] * [T/1e4G] * [1e3mT/T] * [1e2cm/m] => *1e1 [mT/m]
    gu = garray(2:2:garrlen*2) * 1e1; % PE [G/cm] * [T/1e4G] * [1e3mT/T] * [1e2cm/m] => *1e1 [mT/m]
    
    %--------------------------------------------------------------------------
    % Subsample the gradient waveform
    %--------------------------------------------------------------------------
    gv_base_waveform = cat(1, 0, gv(1:10:end), 0);
    gu_base_waveform = cat(1, 0, gu(1:10:end), 0);
else
    % Rewinder design
    gv = garray(1:2:garrlen*2); % RO [G/cm]
    gu = garray(2:2:garrlen*2); % PE [G/cm]
    
    %--------------------------------------------------------------------------
    % Subsample the gradient waveform
    %--------------------------------------------------------------------------
    gv_base_waveform = cat(1, 0, gv(1:10:end));
    gu_base_waveform = cat(1, 0, gu(1:10:end));
    g = cat(2, gv_base_waveform, gu_base_waveform);
    N = length(g);
    
    [k, ~, ~, ~] = calcgradinfo(g, sys.gradRasterTime); % g: N x 2, m1: N x 2
    
    %--------------------------------------------------------------------------
    % starting estimate for the number of points in a possible solution
    % if 2x1, Nest(1) is upper, Nest(2) is lower
    %--------------------------------------------------------------------------
    Nest = 60;
    
    %--------------------------------------------------------------------------
    % 1 x D array; gradient value at start [G/cm]
    %--------------------------------------------------------------------------
    g0 = g(N,:);
    
    %--------------------------------------------------------------------------
    % 1 x D array; gradient value at end [G/cm]
    %--------------------------------------------------------------------------
    gf = [0 0];
    
    %--------------------------------------------------------------------------
    % Q x D array; moment of gradient (/cm, s/cm, s^2/cm...)
    %--------------------------------------------------------------------------
    flag_M1 = 0;
    if flag_M1
        moment = [-k(N,:); -m1(N,:)]; % zeroth (M0) and first (M1) moments
    else
        moment = -k(N,:); % zeroth moment (M0) only
    end
    
    %--------------------------------------------------------------------------
    % 1 x 1; Sampling time [sec]
    %--------------------------------------------------------------------------
    T = sys.gradRasterTime;
    
    %--------------------------------------------------------------------------
    % 1 x 1; Maximum gradient amplitude [G/cm]
    % [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    %--------------------------------------------------------------------------
    gmax = max_grad * 1e-1; % [G/cm]
    
    %--------------------------------------------------------------------------
    % 1 x 1;  Maximum slew rate [G/cm/s]
    % [mT/m/ms] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] * [1e3ms/s] => *1e2 [G/cm/s]
    %--------------------------------------------------------------------------
    smax = max_slew * 1e2;
    
    %--------------------------------------------------------------------------
    % 1 x 1; time at start of gradient [sec]
    %--------------------------------------------------------------------------
    t0 = N * T; % [sec]
    
    %--------------------------------------------------------------------------
    % Call mintinegrad()
    %--------------------------------------------------------------------------
    g_rewinder = mintimegrad(Nest, g0, gf, moment, T, gmax, smax, t0, 3);
    
    % Combine a spiral readout and a spiral rewinder
    g_spiral_waveform = cat(1, g, g_rewinder);
    plotgradinfo(g_spiral_waveform, T);
    
    gv_base_waveform = g_spiral_waveform(:, 1) * 1e1; % RO [G/cm] * [T/1e4G] * [1e3mT/T] * [1e2cm/m] => *1e1 [mT/m]
    gu_base_waveform = g_spiral_waveform(:, 2) * 1e1; % PE [G/cm] * [T/1e4G] * [1e3mT/T] * [1e2cm/m] => *1e1 [mT/m]
end

set(gcf, 'Position', [200 300 792 420]);
% saveas(gcf, sprintf('grad_ni%d_nt%d.png', nr_interleaves, nr_time_frames));

%% Calculate spiral gradient waveforms in the GCS [Hz/m] ([PE,RO,SL] = [y,x,z] in Pulseq)
grad_samples = length(gv_base_waveform);
g_base_waveform = zeros(grad_samples, 3, 'double');
% [mT/m] * [T/1e3mT] * [Hz/T] => * 1e-3 [Hz/m]
g_base_waveform(:,1) = gu_base_waveform * sys.gamma * 1e-3; % PE (gu)
g_base_waveform(:,2) = gv_base_waveform * sys.gamma * 1e-3; % RO (gv)

%% Calculate rotated spiral gradient waveforms in the GCS [Hz/m] [PE,RO,SL]
% Note that [y,x,z] in Pulseq corresponds to [PE,RO,SL]!
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
%         ^ RO(2)                   ^ RO(2)
%         |                         |  /
%         |                         | /
%         |                         |/  phi
%         +---------> PE(1) =>      +---------> PE(1)
%        /                         /
%       /                         /
%      v SL(3)                   v SL(3)
%
% All rotation matrices are right-handed.
% 1. Apply rotation about the SL direction by phi
%
% R_PE(a) = [1     0     0   ], R_RO(a) = [ cos(a) 0 sin(a)], R_SL(a) = [cos(a) -sin(a) 0]
%           [0 cos(a) -sin(a)]            [   0    1   0   ]            [sin(a)  cos(a) 0]
%           [0 sin(a)  cos(a)]            [-sin(a) 0 cos(a)]            [  0       0    1]
%--------------------------------------------------------------------------
g_gcs = zeros(grad_samples, 3, nr_time_frames, 'double');

for i = 1:nr_time_frames
    cos_phi = cos(phi_series(i));
    sin_phi = sin(phi_series(i));
    
    %----------------------------------------------------------------------
    % Calculate the rotation matrix about the SL direction by phi
    %----------------------------------------------------------------------
    R_SL = [cos_phi    -sin_phi    0 ;
        sin_phi     cos_phi    0 ;
        0           0       1];
    
    %----------------------------------------------------------------------
    % Calculate rotated spiral gradient waveforms
    %----------------------------------------------------------------------
    g_gcs(:,:,i) = (R_SL * g_base_waveform.').';
end

%% Flip the sign of gradients in the PE and SL directions [PE,RO,SL]
% This step is necessary to make use of coordinate transformations in Siemens
% for datasets acquired with Pulseq.
g_gcs(:,1,:) = -g_gcs(:,1,:); % PE
g_gcs(:,3,:) = -g_gcs(:,3,:); % SL

%% Create an ADC event
adc_duration_pre = real_dwell_time * discard_pre;
adc_duration_post = real_dwell_time * discard_post;

delay_pre  = round((sys.adcDeadTime + adc_duration_pre)  / sys.gradRasterTime) * sys.gradRasterTime;
delay_post = round((sys.adcDeadTime + adc_duration_post) / sys.gradRasterTime) * sys.gradRasterTime;

grad_duration = grad_samples * sys.gradRasterTime;
adc_samples = floor((adc_duration_pre + grad_duration + adc_duration_post) / real_dwell_time);
ro_duration = size(g, 1) * sys.gradRasterTime;
ro_samples = floor(ro_duration / real_dwell_time);

adc_delay = delay_pre - adc_duration_pre;
adc = mr.makeAdc(adc_samples, 'Dwell', real_dwell_time, 'delay', adc_delay, 'system', sys);

%% Calculate timing (need to decide on the block structure already)
minTE = gz_ex.flatTime / 2 + gz_ex.fallTime + mr.calcDuration(gz_reph) + 0.0001;
% minTE = 0.00252;
delayTE = TE_series + minTE + delay_pre - (gz_ex.flatTime / 2 + gz_ex.fallTime + delay_pre);
delayTR = 0.001 * ones(size(delayTE)); % 1 ms window for Gz dephasing gradients.
minTR = TE_series + minTE + delay_pre + grad_duration + delay_post + gz_ex.flatTime / 2 + gz_ex.riseTime + delayTR;
% delayTR = TR_series - TE_series - minTE - delay_pre - (grad_duration + delay_post) - (gz_ex.flatTime / 2 + gz_ex.riseTime);

for i = 1:nr_time_frames
    assert(delayTE(i) >= mr.calcDuration(gz_reph));
%     assert(delayTR(i) >= 0);
end

%% Create spiral (readout + rewinder) gradient events ([PE,RO,SL] = [y,x,z] in Pulseq)
%--------------------------------------------------------------------------
% Create a spiral gradient event on the PE direction (PRS)
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
gy_spiral = mr.makeArbitraryGrad('y', g_gcs(:,1,1) * 0, 'Delay', delay_pre, 'system', sys); % PE => 'y'

%--------------------------------------------------------------------------
% Create a spiral gradient event on the RO direction (PRS)
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
gx_spiral = mr.makeArbitraryGrad('x', g_gcs(:,2,1) * 0, 'Delay', delay_pre, 'system', sys); % RO => 'x'


%% Define sequence blocks for FISP data acquisitions
rf_phase_rad = rf_phase * pi / 180;

for k = 1:nr_averages
    for j = 1:nr_repetitions
        for p = 1:nr_slices
            count = 1;
            
            gz_phase = mr.scaleGrad(gz_phase_encoding, pe_steps(p)/pe_steps(1));
            gz_phase.amplitude = gz_phase.amplitude + gz_reph.amplitude;
            gz_phase.area = gz_phase.area + gz_reph.area;
            
            %% Create a spoiler event in the slice direction (rewinder + spoiler)
            %--------------------------------------------------------------------------
            % gradient amplitude in [Hz/m] = [Hz/T] * [mT/m]
            % gradient area in [Hz/m] * [sec] => [Hz/m*sec] = [cycle/m]
            % area = 4 / slice_thickness: 4 cycle dephasing across the slice thickness
            %--------------------------------------------------------------------------
            % gz_spoil = mr.makeTrapezoid('z', 'Area', gz_phase.area + dephasing / slice_thickness, 'system', sys);
            gz_spoil = mr.makeTrapezoid('z', 'Area', nr_slices * dephasing / fov_z - gz_ex.area - gz_phase.area, 'system', sys);
            
            %% Define sequence blocks for an adiabatic inversion pulse
            tstart = tic; fprintf('Defining blocks for an adiabatic inversion pulse... ');
            seq.addBlock(rf_inv, gz_inv);
            seq.addBlock(gz_inv_spoil);
            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
            
            for i = 1:nr_time_frames
                nr_interleave = i + j - 1;
                if nr_interleave > nr_time_frames
                    nr_interleave = mod(nr_interleave, nr_time_frames);
                end
                
                tstart = tic; fprintf('Defining blocks for FISP data acquisitions (%3d/%3d)(%3d/%3d)(%3d/%3d)(%5d/%5d)(%5d/%5d)... ', k, nr_averages, j, nr_repetitions, p, nr_slices, i, nr_time_frames, nr_interleave, nr_time_frames);
                
                %----------------------------------------------------------------------
                % Set the phase of an RF event and an ADC event
                %----------------------------------------------------------------------
                rf_ex.phaseOffset = rf_phase_rad * mod(count,2);
                adc.phaseOffset = rf_phase_rad * mod(count,2);
                
                %----------------------------------------------------------------------
                % Set the flip angle [degree]
                %----------------------------------------------------------------------
                rf_ex.signal = ref_signal * FA_series(i);
                
                %----------------------------------------------------------------------
                % Create a spiral gradient event on the PE direction (PRS)
                % [PE,RO,SL] = [y,x,z] in Pulseq
                %----------------------------------------------------------------------
                gy = gy_spiral;
                gy.waveform = g_gcs(:,1,nr_interleave); % PE => 'y'
                
                %----------------------------------------------------------------------
                % Create a spiral gradient event on the RO direction (PRS)
                % [PE,RO,SL] = [y,x,z] in Pulseq
                %----------------------------------------------------------------------
                gx = gx_spiral;
                gx.waveform = g_gcs(:,2,nr_interleave); % RO => 'x'
                
                %----------------------------------------------------------------------
                % Add a new block to the sequence
                %----------------------------------------------------------------------
                seq.addBlock(rf_ex, gz_ex);
                seq.addBlock(gz_phase, mr.makeDelay(delayTE(i)));
                seq.addBlock(gx, gy, adc, mr.makeDelay(delay_pre + grad_duration + delay_post));
                seq.addBlock(mr.align('right', mr.makeDelay(delayTR(i)), gz_spoil));
                count = count + 1;
                fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
            end
            
            %% Define a sequence block of constant delay
            tstart = tic; fprintf('Defining blocks for a constant delay... ');
            delay_const = 2; % [sec]
            seq.addBlock(mr.makeDelay(delay_const));
            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        end
    end
end

%% check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Prepare sequence export
seq.setDefinition('BaseResolution', base_resolution);
seq.setDefinition('Dephasing', dephasing);
seq.setDefinition('DiscardPre', discard_pre);
seq.setDefinition('DiscardPost', discard_post);
seq.setDefinition('FOV', [fov_read fov_read slice_thickness]);
seq.setDefinition('Interleaves', nr_interleaves);
seq.setDefinition('Repetitions', nr_repetitions);
seq.setDefinition('Averages', nr_averages);
seq.setDefinition('Frames', nr_time_frames);
seq.setDefinition('MaxGrad', max_grad);
seq.setDefinition('MaxSlew', max_slew);
seq.setDefinition('Name', 'FISP 3D spiral');
seq.setDefinition('RealDwellTime', real_dwell_time);
seq.setDefinition('Resolution', resolution);
seq.setDefinition('TE', minTE);
seq.setDefinition('TR', minTR);
seq.setDefinition('RasterTime', raster_time);
seq.setDefinition('gtype', gtype);
seq.setDefinition('sptype', sptype);
seq.setDefinition('us_0', us_0);
seq.setDefinition('us_1', us_1);
seq.setDefinition('us_r', us_r);
seq.setDefinition('ustype', ustype);

seq_filename = sprintf('fisp_mrf_3d_fa%d_cycle%d_ni%d_nt%d_%s.seq', flip_angle, dephasing, nr_interleaves, nr_time_frames, opt_data);
seq_path = fullfile(pwd, seq_filename);
seq.write(seq_path); % Write to a pulseq file

%% Save important design parameters into a separate .mat file.
params_fname = sprintf('fisp_mrf_3d_seq_params_ni%d_nt%d_%s.mat', nr_interleaves, nr_time_frames, opt_data);
save(params_fname, 'nr_interleaves', 'nr_time_frames', 'ro_samples', 'ro_duration', 'minTR', 'minTE', 'fn_mrf_seq_path');

%% Plot sequence and k-space diagrams
if 0
    seq.plot('timeRange', [0 10] * TR_base);
end
%seq.plot();

if 0
    % k-space trajectory calculation
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
    
    % plot k-spaces
    figure; plot(ktraj(1,:), ktraj(2,:), 'b'); % a 2D k-space plot
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    hold; plot(ktraj_adc(1,:), ktraj_adc(2,:), 'r.'); % plot the sampling points
    title('full k-space trajectory (k_x x k_y)');
end

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits
% rep = seq.testReport;
% fprintf([rep{:}]);

%% Write kspace trajectories.
traj_name = sprintf('traj_vdspiral_ni%d_nt%d_%s.h5', nr_interleaves, nr_time_frames, opt_data);
if ~isfile(traj_name)
    k_gcs = cumsum(g_gcs * (2 * pi) * sys.gradRasterTime, 1);
    k_gcs = permute(k_gcs, [2 1 3]);
    
    %% Write an .MRD format for k-space trajectories
    ismrmrd_traj_filename = traj_name;
    ismrmrd_traj_path = fullfile(pwd, ismrmrd_traj_filename);
    dset = ismrmrd.Dataset(ismrmrd_traj_path);
    
    %% Add a block of acquisitions at a time
    %--------------------------------------------------------------------------
    % ISMRMRD header from ismrmrd.h (Header for each MR acquisition)
    %--------------------------------------------------------------------------
    % uint16_t version;                                    /**< First unsigned int indicates the version */
    % uint64_t flags;                                      /**< bit field with flags */
    % uint32_t measurement_uid;                            /**< Unique ID for the measurement */
    % uint32_t scan_counter;                               /**< Current acquisition number in the measurement */
    % uint32_t acquisition_time_stamp;                     /**< Acquisition clock */
    % uint32_t physiology_time_stamp[ISMRMRD_PHYS_STAMPS]; /**< Physiology time stamps, e.g. ecg, breating, etc. */
    % uint16_t number_of_samples;                          /**< Number of samples acquired */
    % uint16_t available_channels;                         /**< Available coils */
    % uint16_t active_channels;                            /**< Active coils on current acquisiton */
    % uint64_t channel_mask[ISMRMRD_CHANNEL_MASKS];        /**< Mask to indicate which channels are active. Support for 1024 channels */
    % uint16_t discard_pre;                                /**< Samples to be discarded at the beginning of  acquisition */
    % uint16_t discard_post;                               /**< Samples to be discarded at the end of acquisition */
    % uint16_t center_sample;                              /**< Sample at the center of k-space */
    % uint16_t encoding_space_ref;                         /**< Reference to an encoding space, typically only one per acquisition */
    % uint16_t trajectory_dimensions;                      /**< Indicates the dimensionality of the trajectory vector (0 means no trajectory) */
    % float sample_time_us;                                /**< Time between samples in micro seconds, sampling BW */
    % float position[3];                                   /**< Three-dimensional spatial offsets from isocenter */
    % float read_dir[3];                                   /**< Directional cosines of the readout/frequency encoding */
    % float phase_dir[3];                                  /**< Directional cosines of the phase */
    % float slice_dir[3];                                  /**< Directional cosines of the slice direction */
    % float patient_table_position[3];                     /**< Patient table off-center */
    % ISMRMRD_EncodingCounters idx;                        /**< Encoding loop counters, see above */
    % int32_t user_int[ISMRMRD_USER_INTS];                 /**< Free user parameters */
    % float user_float[ISMRMRD_USER_FLOATS];               /**< Free user parameters */
    %--------------------------------------------------------------------------
    % Where EncodingCounters are defined as:
    % uint16_t kspace_encode_step_1;    /**< e.g. phase encoding line number */
    % uint16_t kspace_encode_step_2;    /**< e.g. partition encoding number */
    % uint16_t average;                 /**< e.g. signal average number */
    % uint16_t slice;                   /**< e.g. imaging slice number */
    % uint16_t contrast;                /**< e.g. echo number in multi-echo */
    % uint16_t phase;                   /**< e.g. cardiac phase number */
    % uint16_t repetition;              /**< e.g. dynamic number for dynamic scanning */
    % uint16_t set;                     /**< e.g. flow encoding set */
    % uint16_t segment;                 /**< e.g. segment number for segmented acquisition */
    % uint16_t user[ISMRMRD_USER_INTS]; /**< Free user parameters */
    %--------------------------------------------------------------------------
    % g_shot_gcs: 3 x grad_samples x nr_inner_shots x nr_shots
    acqblock = ismrmrd.Acquisition(nr_time_frames);
    
    %--------------------------------------------------------------------------
    % Set the header elements that don't change
    %--------------------------------------------------------------------------
    acqblock.head.version(:)               = 1;
    acqblock.head.number_of_samples(:)     = grad_samples;
    acqblock.head.available_channels(:)    = 1;
    acqblock.head.active_channels(:)       = 1;
    acqblock.head.discard_pre(:)           = 0;
    acqblock.head.discard_post(:)          = 0;
    acqblock.head.trajectory_dimensions(:) = 3;
    acqblock.head.sample_time_us(:)        = sys.gradRasterTime * 1e6;
    
    %% Loop over the acquisitions, set the header, set the data and append
    %--------------------------------------------------------------------------
    % Interleaves
    %--------------------------------------------------------------------------
    for interleave = 1:nr_time_frames
        scan_counter = interleave;
        tstart = tic; fprintf('Filling the MRD trajectory (%3d/%3d)... ', scan_counter, nr_time_frames);
        acqblock.head.scan_counter(scan_counter) = scan_counter;
        %         acqblock.head.idx.kspace_encode_step_1(scan_counter) = shot - 1;
        %         acqblock.head.idx.average(scan_counter) = avg - 1;
        %         acqblock.head.idx.contrast(scan_counter) = echo - 1;
        acqblock.traj{scan_counter} = single(k_gcs(:,:,interleave));
        acqblock.data{scan_counter} = zeros(grad_samples, 1, 'single');
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
    
    %--------------------------------------------------------------------------
    % Append the acquisition block
    %--------------------------------------------------------------------------
    dset.appendAcquisition(acqblock);
    
    %% Fill the xml header
    %--------------------------------------------------------------------------
    % We create a matlab struct and then serialize it to xml.
    % Look at the xml schema to see what the field names should be
    %--------------------------------------------------------------------------
    header = [];
    
    %--------------------------------------------------------------------------
    % Experimental Conditions (Required)
    %--------------------------------------------------------------------------
    header.experimentalConditions.H1resonanceFrequency_Hz = 0.55 * sys.gamma; % [T] * [Hz/T] = > [T]
    
    %--------------------------------------------------------------------------
    % Acquisition System Information (Optional)
    %--------------------------------------------------------------------------
    header.acquisitionSystemInformation.systemVendor = 'Pulseq';
    header.acquisitionSystemInformation.systemModel = 'v1.4.0';
    
    %--------------------------------------------------------------------------
    % The Encoding (Required)
    %--------------------------------------------------------------------------
    header.encoding.trajectory = 'spiral';
    header.encoding.encodedSpace.fieldOfView_mm.x = fov_read * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]
    header.encoding.encodedSpace.fieldOfView_mm.y = fov_read * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]
    header.encoding.encodedSpace.fieldOfView_mm.z = slice_thickness * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]
    header.encoding.encodedSpace.matrixSize.x = base_resolution;
    header.encoding.encodedSpace.matrixSize.y = base_resolution;
    header.encoding.encodedSpace.matrixSize.z = 1;
    
    % Recon Space (in this case same as encoding space)
    header.encoding.reconSpace = header.encoding.encodedSpace;
    
    % Encoding Limits
    header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_0.maximum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_0.center = 0;
    header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_1.maximum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_1.center = 0;
    header.encoding.encodingLimits.repetition.minimum = 0;
    header.encoding.encodingLimits.repetition.maximum = 0;
    header.encoding.encodingLimits.repetition.center = 0;
    
    %--------------------------------------------------------------------------
    % Trajectory Description (Optional)
    %--------------------------------------------------------------------------
    header.encoding.trajectoryDescription.identifier = 'spiralgen_jgp_11apr_mex';
    header.encoding.trajectoryDescription.userParameterLong = [];
    header.encoding.trajectoryDescription.userParameterDouble(1).value = spparams(1);
    header.encoding.trajectoryDescription.userParameterDouble(1).name = 'gamma_Hz_per_G';
    header.encoding.trajectoryDescription.userParameterDouble(2).value = spparams(2);
    header.encoding.trajectoryDescription.userParameterDouble(2).name = 'fov_cm';
    header.encoding.trajectoryDescription.userParameterDouble(3).value = spparams(3);
    header.encoding.trajectoryDescription.userParameterDouble(3).name = 'resolution_cm';
    header.encoding.trajectoryDescription.userParameterDouble(4).value = spparams(4);
    header.encoding.trajectoryDescription.userParameterDouble(4).name = 'MaxSlewRate_G_per_cm_per_s';
    header.encoding.trajectoryDescription.userParameterDouble(5).value = spparams(5);
    header.encoding.trajectoryDescription.userParameterDouble(5).name = 'MaxGradientAmplitude_G_per_cm';
    header.encoding.trajectoryDescription.userParameterDouble(6).value = spparams(6);
    header.encoding.trajectoryDescription.userParameterDouble(6).name = 'TotalMaxGradientAmplitude_G_per_cm';
    header.encoding.trajectoryDescription.userParameterDouble(7).value = spparams(7);
    header.encoding.trajectoryDescription.userParameterDouble(7).name = 'RasterTime_s';
    header.encoding.trajectoryDescription.userParameterDouble(8).value = spparams(8);
    header.encoding.trajectoryDescription.userParameterDouble(8).name = 'interleaves';
    header.encoding.trajectoryDescription.userParameterDouble(9).value = spparams(9);
    header.encoding.trajectoryDescription.userParameterDouble(9).name = 'us_0';
    header.encoding.trajectoryDescription.userParameterDouble(10).value = spparams(10);
    header.encoding.trajectoryDescription.userParameterDouble(10).name = 'us_1';
    header.encoding.trajectoryDescription.userParameterDouble(11).value = spparams(11);
    header.encoding.trajectoryDescription.userParameterDouble(11).name = 'us_r';
    header.encoding.trajectoryDescription.userParameterDouble(12).value = spparams(12);
    header.encoding.trajectoryDescription.userParameterDouble(12).name = 'ustype';
    header.encoding.trajectoryDescription.userParameterDouble(13).value = spparams(13);
    header.encoding.trajectoryDescription.userParameterDouble(13).name = 'gtype';
    header.encoding.trajectoryDescription.userParameterDouble(14).value = spparams(14);
    header.encoding.trajectoryDescription.userParameterDouble(14).name = 'sptype';
    header.encoding.trajectoryDescription.comment = 'Using spiral design by James Pipe (https://www.ismrm.org/mri_unbound/sequence.htm)';
    
    %% Serialize and write to the data set
    xmlstring = ismrmrd.xml.serialize(header);
    dset.writexml(xmlstring);
    
    %% Write the dataset
    dset.close();
end