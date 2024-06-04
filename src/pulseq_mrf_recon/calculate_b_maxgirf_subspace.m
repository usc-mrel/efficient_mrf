function b = calculate_b_maxgirf_subspace(Phi, Ue, Ve, d, Omega, csm, st, w)%% Calculate b: A(u) = b
%--------------------------------------------------------------------------
% b = sum_c=1^Nc vec_{N,Lt}(Sc^H * F^H * Omega^H(dc) * Phi^H) in C^(NLt x 1)
% dc in C^(P x 1) : all the measured data from the cth coil
% F in C^(Np x N) : a fully-sampled Fourier encoding matrix
% Sc in C^(N x N) : a diagonal matrix containing the coil sensitivities for the cth coil
% Phi in C^(Lt x M) : the rank-Lt temporal subspace of a dictionary D
% U in C^(Np x Ls): the rank-Ls spatial subspace of higher order encoding.
% V in C^(N x Ls): the rank-Ls temporal subspace of higer order encoding.
% Omega^H() : C^(P x 1) -> C^(Np x M) : a linear operator that undersamples the k-space and then concatenates the measured data into a vector
% vec_{N,Lt}(): C^(N x Lt) -> C^(NLt x 1): a linear operator that converts an N x Lt matrix to an NL x 1 vector
%--------------------------------------------------------------------------
% Other input:
% st in {Narms x 1}: NUFFT structure per interleave by Jeff Fessler's toolbox.
% w: DCF.

% Modified to combine with MaxGIRF encoding operator, Zhibo.
% Optimization to reduce loop, Zhibo.

d = permute(d, [1 4 2 3]); % Nk x M x Narms x Nc
[Ns, M, Narms, Nc] = size(d);
P = Ns; % Fully sampling.
Np = P; % This is true only when fully sampling.
[N1, N2, ~] = size(csm);
Lt = size(Phi, 1);
Ls = size(Ue, 2);
d = reshape(d, [], Narms, Nc);

b = complex(zeros(N1 * N2 * Lt, 1, 'double')); % NL x 1
start_time = tic;
for i = 1 : Narms
    for c = 1 : Nc
        tstart = tic;
        fprintf('Calculating b in A(u) = b (Interleave: %d/%d, coil: %d/%d)... ', i, Narms, c, Nc);
        %----------------------------------------------------------------------
        % Calculate undersampled non-Cartesian k-space data (Np x M)
        % Omega^H(b): P x 1 -> Np x M
        %----------------------------------------------------------------------
        FScUPhi = complex(zeros(Np, M, 'double')); % Np x M
        dic = d(:, i, c);
        FScUPhi(Omega(:, :, i)) = dic(abs(dic) ~= 0); % Np x M

        %----------------------------------------------------------------------
        % Multiply Phi^H
        % Omega^H(b) * Phi^H: (Np x M) * (M x Lt) => Np x Lt
        %----------------------------------------------------------------------
        FscUPhiPhiH = FScUPhi * Phi'; % Np x Lt
        %----------------------------------------------------------------------
        % Perform NUFFT reconstruction per time frame
        % F^H * (Omega^H(b) * Phi^H): N x Np * Np x Lt => N x Lt
        %----------------------------------------------------------------------
%         AcHb = complex(zeros(N1, N2, Lt, 'double'));
        scale_factor = 1 / sqrt(prod(st{i, 1}.Nd));
%         for ell = 1 : Lt
%             AHd = complex(zeros(N1 * N2, 1, 'double'));
%             for ell2 = 1 : Ls
%                 FHDuHd = nufft_adj((conj(Ue(:, ell2, i)) .* FscUPhiPhiH(:, ell)) .* w(:, i), st{i, 1}) * scale_factor; % N1 x N2 x Lt
%                 AHd = AHd + Ve(:,ell2,i) .* reshape(FHDuHd, [N1 * N2 1]);
%             end
%             AcHb(:, :, ell) = reshape(AHd, [N1 N2]);
%         end

        % Verified optimization, ~24% computation reduction in a 2D test
        % case.
        AHd = complex(zeros(N1 * N2, Lt, 'double'));
        for ell = 1 : Ls
            FHDuHd = nufft_adj((repmat(conj(Ue(:, ell, i)), [1 Lt]) .* FscUPhiPhiH) .* repmat(w(:, i), [1 Lt]), st{i, 1}) * scale_factor;
            AHd = AHd + Ve(:, ell, i) .* reshape(FHDuHd, [N1 * N2 Lt]);
        end
        AcHb = reshape(AHd, [N1 N2 Lt]);

        %----------------------------------------------------------------------
        % Sc^H * (F^H * (Omega^H(b) * Phi^H)): (N x N) * (N x Lt) => N x Lt
        %---------------------------------------------------------------------- 
        AcHb = bsxfun(@times, conj(csm(:, :, c)), AcHb); % N1 x N2 x Lt
        b = b + AcHb(:); % N x Lt => NLt x 1
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end