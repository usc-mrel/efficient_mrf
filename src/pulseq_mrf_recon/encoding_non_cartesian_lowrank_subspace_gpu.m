function outp = encoding_non_cartesian_lowrank_subspace_gpu(inp, csm, Omega, Phi, nufft_st, w, transpose_indicator)
% Written by Nam Gyun Lee
% Email: nmgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 02/09/2022, Last modified: 02/13/2022
% Optimized with reduced loops, Zhibo.
% Add GPU option, Zhibo.

%% Declare a persistent variable
persistent cg_iter;
if isempty(cg_iter)
    cg_iter = 0;
end

%% Define dimensions
[N1, N2, Nc] = size(csm);
N = N1 * N2; % total number of voxels
[Lt, M] = size(Phi);
NLt = N * Lt;
Np = numel(w);

%% Determine the operator type
if strcmp(transpose_indicator, 'transp')
    operator_type = 'adjoint';
elseif strcmp(transpose_indicator, 'notransp')
    operator_type = 'forward';
    cg_iter = cg_iter + 1;
end

outp = complex(zeros(NLt, 1, 'double', 'gpuArray'));

fprintf('Calculating the %s operator...\n', operator_type);

%--------------------------------------------------------------------------
% Calculate U
% U = vec_{N,Lt}^H(u): NLt x 1 -> N x Lt
%--------------------------------------------------------------------------
U = reshape(inp, [N1 N2 Lt]); % NLt x 1 -> N1 x N2 x Lt

start_time = tic;
for c = 1:Nc
    tstart = tic; fprintf('(CG=%2d): Calculating the %s operator (c=%2d/%2d)... ', cg_iter, operator_type, c, Nc);
    %----------------------------------------------------------------------
    % Calculate fully-sampled non-Cartesian k-space data (Np x M)
    % F * Sc * U: (Np x N) * (N x N) * (N x Lt) -> Np x Lt
    %----------------------------------------------------------------------
%     FScU = complex(zeros(Np, Lt, 'double', 'gpuArray')); % Np x Lt
    scale_factor = 1 / sqrt(prod(nufft_st{1, 1}.Nd));
%     for ell = 1:Lt
%         FScU(:,ell) = nufft(csm(:,:,c) .* U(:,:,ell), nufft_st{1, 1}) * scale_factor; % Np x 1
%     end

    % Verified optimization, ~7% computation time reduction, Zhibo.
    FScU = nufft_gpu(repmat(csm(:, :, c), [1 1 Lt]) .* U, nufft_st{1, 1}) * scale_factor;

    %----------------------------------------------------------------------
    % Multiply Phi
    % (F * Sc * U) * Pih: (Np x Lt) * (Lt x M) -> Np x M
    %----------------------------------------------------------------------
    FScUPhi = FScU * Phi; % Np x M

    %----------------------------------------------------------------------
    % Calculate a vector of k-space measurements (P x 1)
    % Ac(U) = Omega(F * Sc * U * Phi): Np x M -> P x 1
    %----------------------------------------------------------------------
    dc = FScUPhi(Omega); % P x 1
    %----------------------------------------------------------------------
    % Calculate undersampled non-Cartesian k-space data (Np x M)
    % Omega^H(b): P x 1 -> Np x M
    %----------------------------------------------------------------------
    FScUPhi = complex(zeros(Np, M, 'double', 'gpuArray')); % Np x M
    FScUPhi(Omega) = dc; % Np x M

    %----------------------------------------------------------------------
    % Multiply Phi^H
    % Omega^H(b) * Phi^H: (Np x M) * (M x Lt) => Np x Lt
    %----------------------------------------------------------------------
    FscUPhiPhiH = FScUPhi * Phi'; % Np x Lt

    %----------------------------------------------------------------------
    % Perform NUFFT reconstruction per time frame
    % F^H * (Omega^H(b) * Phi^H): (N x Np) * (Np x Lt) => N x Lt
    %----------------------------------------------------------------------
%     AcHb = complex(zeros(N1, N2, Lt, 'double', 'gpuArray'));
%     for ell = 1:Lt
%         AcHb(:,:,ell) = nufft_adj(FscUPhiPhiH(:,ell) .* w(:), nufft_st{1, 1}) * scale_factor;
%     end

    % Verified optimization, ~7% computation time reduction, Zhibo.
    AcHb = nufft_adj_gpu(FscUPhiPhiH .* repmat(w(:), [1 Lt]), nufft_st{1, 1}) * scale_factor;
    %----------------------------------------------------------------------
    % Sc^H * (F^H * (Omega^H(b) * Phi^H)): N x N * N x Lt => N x Lt
    %----------------------------------------------------------------------
    AcHb = bsxfun(@times, conj(csm(:,:,c)), AcHb); % N1 x N2 x Lt
    outp = outp + AcHb(:); % N x Lt => NLt x 1
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end
end