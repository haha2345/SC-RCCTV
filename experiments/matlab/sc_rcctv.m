function [RecoverImg, info] = sc_rcctv(Nhsi, o)
%% SC-RCCTV: self-calibrated subspace RCCTV (band-wise Huber residual + penalised stripe).
%  Y = U*V' + S + G + T,   RCCTV regulariser on spatial gradients of U (Hu et al.),
%  lambda*|S|_1 + sum_k beta_k/2 |G_k|^2 + lambdaT*|a|_1,  T = 1_M*a'.
%  Options (struct o, all optional)
%    rk=8, lam=6/sqrt(MN), tol=1e-5, maxIter=100
%    T_mode : 'pen' (default, blind) | 'ball2' (oracle radius, every-2) | 'none'
%    G_mode : 'band' (default, blind sigma_k) | 'peng' (given sigma) | 'none'
%    beta   : stripe L1-ball radius (T_mode='ball2');  sigma: scalar/vector (G_mode='peng')
%    cT=3, sig_online=true, online_start=5, scalar_sigma=false
%  Official pieces untouched: G1/G2 unweighted SVT, FFT U-update, rho=1.5, mu schedule.
%  Self-contained: circular differences, their adjoints and soft-thresholding are local functions below
%  (bitwise identical to diff_x/diff_xT/diff_y/diff_yT of the official RCCTV code); needs psf2otf (Image Processing).
if nargin < 2, o = struct(); end
[M,N,p] = size(Nhsi); MN = M*N;
rk   = gf(o,'rk',8);          lam  = gf(o,'lam',6/sqrt(MN));
tol  = gf(o,'tol',1e-5);      maxIter = gf(o,'maxIter',100);
Tm   = gf(o,'T_mode','pen');  Gm = gf(o,'G_mode','band');
beta = gf(o,'beta',0);        sigma = gf(o,'sigma',0);
cT   = gf(o,'cT',3);          sigon = gf(o,'sig_online',true);
son  = gf(o,'online_start',5); scal = gf(o,'scalar_sigma',false);
sdir = gf(o,'stripe_dirs','vh');   % 'v' column-constant, 'h' row-constant, 'vh' both (orientation-agnostic)
rho = 1.5; max_mu = 1e6; eps = 1e-4;
D = reshape(double(Nhsi), [MN, p]);
norm_two  = norm(D, 2);
dual_norm = max(norm_two, max(abs(D(:))));
mu = 1/norm_two; normD = norm(D,'fro');
sizeU = [M,N,rk];
Eny_x = (abs(psf2otf([+1; -1], [M,N,rk]))).^2;
Eny_y = (abs(psf2otf([+1, -1], [M,N,rk]))).^2;
determ = Eny_x + Eny_y;
[u,s,v] = svd(D,'econ');
U = u(:,1:rk)*s(1:rk,1:rk);  V = v(:,1:rk);
E = zeros(MN,p); T = zeros(MN,p); G = zeros(MN,p);
M1 = zeros(MN,rk); M2 = zeros(MN,rk); M3 = D/dual_norm;
nu = 3*sqrt(2)/(sqrt(MN)+sqrt(p));          % Peng dual noise scale
switch Gm
    case 'peng', sig_k = zeros(1,p) + sigma(:)';
    case 'band'
        sig_k = noise_sigma_vdiff(Nhsi, sdir);
        if scal, sig_k = zeros(1,p) + median(sig_k); end
    otherwise,   sig_k = zeros(1,p);
end
if strcmp(Gm,'none'), beta_g = inf(1,p); else, beta_g = nu*(sqrt(MN)+sqrt(p))./((sqrt(MN)+sqrt(p))*sig_k + 1e-6); end
sa = zeros(1,N,p); sb = zeros(M,1,p);        % stripe amplitudes (column-/row-constant)
iter = 0;
while iter < maxIter
    iter = iter + 1;
    [B,C,P] = svd(reshape(cdiff_x(U,sizeU),[MN,rk]) + M1/mu, 'econ');
    G1 = B*softthre(C,1/mu)*P';
    [B,C,P] = svd(reshape(cdiff_y(U,sizeU),[MN,rk]) + M2/mu, 'econ');
    G2 = B*softthre(C,1/mu)*P';
    %% stripe
    if strcmp(Tm,'ball2')
        if beta <= 0
            T = zeros(MN,p);
        elseif rem(iter,2) == 0
            R = D - U*V' - E - G + M3/mu;
            Tmean = mean(reshape(R,[M N p]), 1);
            Tc = repmat(Tmean, [M 1 1]);
            T = reshape(project_l1_ball_oct(Tc(:), beta), [MN p]);
        end
    elseif strcmp(Tm,'pen')
        R3 = reshape(D - U*V' - E - G + M3/mu, [M N p]);
        if any(sdir == 'v')   % column-constant stripes a (1xNxp): threshold cT*nu/sqrt(M)
            sa = softthre(mean(R3 - sb, 1), cT*nu/sqrt(M)/mu);
        end
        if any(sdir == 'h')   % row-constant stripes b (Mx1xp): threshold cT*nu/sqrt(N)
            sb = softthre(mean(R3 - sa, 2), cT*nu/sqrt(N)/mu);
        end
        T = reshape(sa + sb, [MN p]);
    end
    %% U, V
    diffT_p = cdiff_xT(G1-M1/mu,sizeU) + cdiff_yT(G2-M2/mu,sizeU);
    temp    = (D-E-T-G+M3/mu)*V;
    numer1  = reshape(diffT_p + temp(:), sizeU);
    x = real(ifftn(fftn(numer1) ./ (determ + 1 + eps)));
    U = reshape(x,[MN,rk]);
    [uu,~,vv] = svd((D-E-T-G+M3/mu)'*U, 'econ');
    V = uu*vv';
    X = U*V';
    E = softthre(D-X-T-G+M3/mu, lam/mu);
    if ~strcmp(Gm,'none')
        if sigon && iter >= son
            r = D - X - T;
            md = median(r,1);
            sig_k = 1.4826*median(abs(r - md),1);
            if scal, sig_k = zeros(1,p) + median(sig_k); end
            beta_g = nu*(sqrt(MN)+sqrt(p))./((sqrt(MN)+sqrt(p))*sig_k + 1e-6);
        end
        G = (M3 + mu*(D-X-E-T)) ./ (beta_g + mu);
    end
    leq1 = reshape(cdiff_x(U,sizeU),[MN,rk]) - G1;
    leq2 = reshape(cdiff_y(U,sizeU),[MN,rk]) - G2;
    leq3 = D - X - E - T - G;
    if norm(leq1,'fro')/normD < tol && norm(leq2,'fro')/normD < tol && norm(leq3,'fro')/normD < tol
        break;
    end
    M1 = M1 + mu*leq1; M2 = M2 + mu*leq2; M3 = M3 + mu*leq3;
    mu = min(max_mu, mu*rho);
end
RecoverImg = reshape(U*V', [M,N,p]);
info = struct('iters',iter,'sig_k',sig_k,'T',T,'G',G,'E',E,'beta_g',beta_g);
end

function s = noise_sigma_vdiff(Y, dirs)
% stripe-invariant per-band sigma: MAD of vertical differences (column-constant stripes cancel);
% dirs='vh' also uses horizontal differences (row-constant stripes cancel) and takes the smaller.
[M,N,p] = size(Y);
s = mad_diff(reshape(diff(Y,1,1), [(M-1)*N, p]));
if any(dirs == 'h'), s = min(s, mad_diff(reshape(diff(Y,1,2), [M*(N-1), p]))); end
end
function s = mad_diff(d)
d = d - median(d,1);
s = 1.4826*median(abs(d),1)/sqrt(2);
end

function v = gf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

% ---- circular forward differences along rows (x) and columns (y) of each coefficient image, and their adjoints
function d = cdiff_x(x, sz)
T = reshape(x, sz); d = reshape(circshift(T, -1, 1) - T, [], 1);
end
function d = cdiff_xT(a, sz)
T = reshape(a(1:prod(sz)), sz); d = reshape(circshift(T, 1, 1) - T, [], 1);
end
function d = cdiff_y(x, sz)
T = reshape(x, sz); d = reshape(circshift(T, -1, 2) - T, [], 1);
end
function d = cdiff_yT(a, sz)
T = reshape(a(1:prod(sz)), sz); d = reshape(circshift(T, 1, 2) - T, [], 1);
end
function out = softthre(a, tau)
out = sign(a) .* max(abs(a) - tau, 0);
end
