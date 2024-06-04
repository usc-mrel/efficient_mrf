function outp = encoding_lowrank_maxgirf_subspace(inp, csm, Omega, Phi, Ue, Ve, w, st, transpose_indicator)
% Modified by Zhibo.
% Combined with MaxGIRF encoding and optimized with reduced loops.

%% Declare a persistent variable
persistent cg_iter;
if isempty(cg_iter)
    cg_iter = 0;
end

%% Define dimensions
[N1, N2, Nc] = size(csm);
N = N1 * N2; % total number of voxels
[Lt, M] = size(Phi);
[Ns, Ls, Narms] = size(Ue);
NLt = N * Lt;
Np = size(w, 1);

%% Determine the operator type
if strcmp(transpose_indicator, 'transp')
    operator_type = 'adjoint';
elseif strcmp(transpose_indicator, 'notransp')
    operator_type = 'forward';
    cg_iter = cg_iter + 1;
end

outp = complex(zeros(NLt, 1, 'double'));

fprintf('Calculating the %s operator...\n', operator_type);

%--------------------------------------------------------------------------
% Calculate U
% U = vec_{N,Lt}^H(u): NLt x 1 -> N x Lt
%--------------------------------------------------------------------------
U = reshape(inp, [N1 N2 Lt]); % NLt x 1 -> N1 x N2 x Lt

start_time = tic;
for i = 1 : Narms
    for c = 1 : Nc
        tstart = tic; fprintf('(CG=%2d): Calculating the %s operator (i=%2d/%2d, c=%2d/%2d)... ', cg_iter, operator_type, i, Narms, c, Nc);
        %----------------------------------------------------------------------
        % Calculate fully-sampled non-Cartesian k-space data (Np x M)
        % F * Sc * U: (Np x N) * (N x N) * (N x Lt) -> Np x Lt
        %----------------------------------------------------------------------
%         EScU = complex(zeros(Np, Lt, 'double')); % Np x Lt
        scale_factor = 1 / sqrt(prod(st{i, 1}.Nd));
%         for ell = 1:Lt
%             EiScU = complex(zeros(Ns, 1, 'double'));
%             for ell2 = 1 : Ls
%                 ScU = csm(:,:,c) .* U(:,:,ell);
%                 ScU = reshape(ScU, [N 1]);
%                 EiScU = EiScU + Ue(:,ell2,i) .* nufft(reshape(conj(Ve(:,ell2,i)) .* ScU, [N1 N2]), st{i, 1}) * scale_factor; % Np x 1
%             end
%             EScU(:,ell) = EiScU;
%         end

        EScU = complex(zeros(Ns, Lt, 'double'));
        for ell = 1 : Ls
            ScU = repmat(csm(:, :, c), [1 1 Lt]) .* U;
%             ScU = reshape(ScU, [N Lt]);
            EScU = EScU + repmat(Ue(:, ell, i), [1 Lt]) .* nufft(repmat(reshape(conj(Ve(:,ell,i)), [N1 N2]), [1 1 Lt]) .* ScU, st{i, 1}) * scale_factor;
        end

        %----------------------------------------------------------------------
        % Multiply Phi
        % (F * Sc * U) * Pih: (Np x Lt) * (Lt x M) -> Np x M
        %----------------------------------------------------------------------
        EScUPhi = EScU * Phi; % Np x M

        %----------------------------------------------------------------------
        % Calculate a vector of k-space measurements (P x 1)
        % Ac(U) = Omega(F * Sc * U * Phi): Np x M -> P x 1
        %----------------------------------------------------------------------
        dic = EScUPhi(Omega(:, :, i)); % P x 1
        %----------------------------------------------------------------------
        % Calculate undersampled non-Cartesian k-space data (Np x M)
        % Omega^H(b): P x 1 -> Np x M
        %----------------------------------------------------------------------
        EScUPhi = complex(zeros(Np, M, 'double')); % Np x M
        EScUPhi(Omega(:, :, i)) = dic; % Np x M

        %----------------------------------------------------------------------
        % Multiply Phi^H
        % Omega^H(b) * Phi^H: (Np x M) * (M x Lt) => Np x Lt
        %----------------------------------------------------------------------
        EscUPhiPhiH = EScUPhi * Phi'; % Np x Lt

        %----------------------------------------------------------------------
        % Perform NUFFT reconstruction per time frame
        % F^H * (Omega^H(b) * Phi^H): (N x Np) * (Np x Lt) => N x Lt
        %----------------------------------------------------------------------
%         AcHb = complex(zeros(N1, N2, Lt, 'double'));
        %     for ell = 1:Lt
        %         AcHb(:,:,ell) = nufft_adj(EscUPhiPhiH(:,ell) .* w(:), st{i, 1}) * scale_factor;
        %     end
%         for ell = 1 : Lt
%             AHd = complex(zeros(N1 * N2, 1, 'double'));
%             for ell2 = 1 : Ls
%                 FHDuHd = nufft_adj((conj(Ue(:, ell2, i)) .* EscUPhiPhiH(:, ell)) .* w(:, i), st{i, 1}) * scale_factor; % N1 x N2 x Lt
%                 AHd = AHd + Ve(:,ell2,i) .* reshape(FHDuHd, [N1 * N2 1]);
%             end
%             AcHb(:, :, ell) = reshape(AHd, [N1 N2]);
%         end

        AHd = complex(zeros(N1 * N2, Lt, 'double'));
        for ell = 1 : Ls
            FHDuHd = nufft_adj((repmat(conj(Ue(:, ell, i)), [1 Lt]) .* EscUPhiPhiH) .* repmat(w(:, i), [1 Lt]), st{i, 1}) * scale_factor;
            AHd = AHd + Ve(:, ell, i) .* reshape(FHDuHd, [N1 * N2 Lt]);
        end
        AcHb = reshape(AHd, [N1 N2 Lt]);

        %----------------------------------------------------------------------
        % Sc^H * (F^H * (Omega^H(b) * Phi^H)): N x N * N x Lt => N x Lt
        %----------------------------------------------------------------------
        AcHb = bsxfun(@times, conj(csm(:,:,c)), AcHb); % N1 x N2 x Lt
        outp = outp + AcHb(:); % N x Lt => NLt x 1
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end
end