function b = calculate_b_subspace_gpu(Phi, d, Omega, csm, nufft_st, w)
% Optimized with reduced loops, Zhibo.
% Add GPU option, Zhibo.

%% Calculate b: A(u) = b
%--------------------------------------------------------------------------
% b = sum_c=1^Nc vec_{N,Lt}(Sc^H * F^H * Omega^H(dc) * Phi^H) in C^(NLt x 1)
% dc in C^(P x 1) : all the measured data from the cth coil
% F in C^(Np x N) : a fully-sampled Fourier encoding matrix
% Sc in C^(N x N) : a diagonal matrix containing the coil sensitivities for the cth coil
% Phi in C^(Lt x M) : the rank-Lt temporal subspace of a dictionary D
% Omega^H() : C^(P x 1) -> C^(Np x M) : a linear operator that undersamples the k-space and then concatenates the measured data into a vector
% vec_{N,Lt}(): C^(N x Lt) -> C^(NLt x 1): a linear operator that converts an N x Lt matrix to an NL x 1 vector
%--------------------------------------------------------------------------
% Other input:
% nufft_st: NUFFT structure by Jeff Fessler's toolbox.
% w: DCF.

d = gpuArray(permute(d, [1 2 4 3]));
[Ns, Narms, M, Nc] = size(d);
P = sum(Omega(:, 1));
Np = size(Omega, 1); % This is true only when fully sampling.
[N1, N2, ~] = size(csm);
Lt = size(Phi, 1);
d = reshape(d, [], Nc);

b = complex(zeros(N1 * N2 * Lt, 1, 'double', 'gpuArray')); % NL x 1
start_time = tic;
for c = 1 : Nc
    tstart = tic;
    fprintf('Calculating b in A(u) = b (%d/%d)... ', c, Nc);
    %----------------------------------------------------------------------
    % Calculate undersampled non-Cartesian k-space data (Np x M)
    % Omega^H(b): P x 1 -> Np x M
    %----------------------------------------------------------------------
    FScUPhi = complex(zeros(Np, M, 'double', 'gpuArray')); % Np x M
    FScUPhi(Omega) = d(:, c); % Np x M

    %----------------------------------------------------------------------
    % Multiply Phi^H
    % Omega^H(b) * Phi^H: (Np x M) * (M x Lt) => Np x Lt
    %----------------------------------------------------------------------
    FscUPhiPhiH = FScUPhi * Phi'; % Np x Lt
    %----------------------------------------------------------------------
    % Perform NUFFT reconstruction per time frame
    % F^H * (Omega^H(b) * Phi^H): N x Np * Np x Lt => N x Lt
    %----------------------------------------------------------------------
%     AcHb = complex(zeros(N1, N2, Lt, 'double', 'gpuArray'));
    scale_factor = 1 / sqrt(prod(nufft_st{1, 1}.Nd));
%     for ell = 1 : Lt
%         AcHb(:, :, ell) = nufft_adj(FscUPhiPhiH(:, ell) .* w(:), nufft_st{1, 1}) * scale_factor; % N1 x N2 x Lt
%     end
    
    % Verified optimization, Zhibo.
    AcHb = nufft_adj_gpu(FscUPhiPhiH .* repmat(w(:), [1 Lt]), nufft_st{1, 1}) * scale_factor;

    %----------------------------------------------------------------------
    % Sc^H * (F^H * (Omega^H(b) * Phi^H)): (N x N) * (N x Lt) => N x Lt
    %---------------------------------------------------------------------- 
    AcHb = bsxfun(@times, conj(csm(:, :, c)), AcHb); % N1 x N2 x Lt
    b = b + AcHb(:); % N x Lt => NLt x 1
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end
