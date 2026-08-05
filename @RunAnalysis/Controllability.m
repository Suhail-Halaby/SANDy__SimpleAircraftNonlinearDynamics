function C = Controllability(obj, A, B, opts)
arguments
    obj %#ok<INUSA>
    A (:,:) double
    B (:,:) double
    opts.tol   (1,1) double = -1      % -1 => auto tolerance
    opts.Print (1,1) logical = true
    opts.Label (1,1) string  = ""
end

n = size(A,1);

% --- Kalman rank via SVD (robust to bad scaling) ---
Co  = ctrb(A, B);
s   = svd(Co);
if opts.tol < 0
    tol = max(size(Co)) * eps(max(s));   % same rule rank() uses
else
    tol = opts.tol;
end
r          = sum(s > tol);
fullRank   = (r == n);

% --- PBH test: [A - lambda*I, B] must be full row rank at every eigenvalue ---
ev        = eig(A);
pbhMinSig = zeros(numel(ev),1);
for i = 1:numel(ev)
    M            = [A - ev(i)*eye(n), B];
    pbhMinSig(i) = min(svd(M));         % ~0 => that mode is uncontrollable
end
[worstSig, iw] = min(pbhMinSig);

C = struct('rank',r, 'n',n, 'isControllable',fullRank, ...
    'svals',s, 'tol',tol, ...
    'eig',ev, 'pbhMinSig',pbhMinSig, ...
    'worstMode',ev(iw), 'worstSig',worstSig);

if opts.Print
    if fullRank
        verdict = "CONTROLLABLE";
    else
        verdict = "NOT controllable";
    end
    fprintf('--- Controllability %s ---\n', opts.Label);
    fprintf('  rank(ctrb) = %d / %d  ->  %s\n', r, n, verdict);
    fprintf('  ctrb singular values: %s\n', mat2str(s',3));
    fprintf('  worst PBH mode: lambda = %.4g%+.4gi, min sigma = %.3g\n', ...
        real(C.worstMode), imag(C.worstMode), worstSig);
end
end