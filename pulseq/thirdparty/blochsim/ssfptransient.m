%	== Example of transient SSFP response calculation using bloch.m ==
%	

%	--- Setup parameters ---
TR = .005;      % Sec.
Trf = 0.0001;   % 100 us "hard" RF pulse.
alpha = 60;     % Degrees.
gamma = 4258;   % Hz/G.
T1 = 1;         % Sec.
T2 = .2;        % Sec.
freq = [-200:200]; % Hz
N = 100;
Tpad = (TR-Trf)/2; % Sec.
plot = 1;

% 	--- Setup B1, just delta-function RF ---
%
%	-- Time intervals are [Tpad, Trf, Tpad]
t = [Tpad Trf Tpad];

%	-- b1 is non-zero during Trf.
b1 = [0 pi/180*alpha/Trf/gamma/2/pi 0];	% Gauss.



% === Calculate steady state, for comparison.
[mxss,myss,mzss] = bloch(b1,0*b1,t,T1,T2,freq,0,1);
mss = mxss+i*myss;


% === Start with alpha/2, from equilibrium. ===
[mx,my,mz] = bloch(+i*max(b1)/2,0,Trf,T1,T2,freq,0,0);


for n=1:N
	[mx,my,mz] = bloch(b1,0*b1,t,T1,T2,freq,0,0,mx,my,mz);

	if ((plot==1) | (n==N))

		sig = mx+i*my;		% In-plane signal
		subplot(2,1,1);
		magphase(freq(:),[sig(:) mss(:)]);
	
		drawnow;
		pause;
	end;
end;



