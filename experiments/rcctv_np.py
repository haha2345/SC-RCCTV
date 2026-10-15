"""Numpy port of official Hu RCCTV_RPCA and of the generalised family used here.

Faithful to code/official_rcctv/RCCTV_RPCA.m (circular differences, FFT U-update,
column-major reshapes).  One solver, switchable terms:

  T_mode : 'none' | 'ball2' (oracle L1-ball radius, every-2, = innert2g)
           | 'pen'   (blind penalised column-sparse stripe, every iteration)
  G_mode : 'none' | 'peng' (scalar oracle sigma, = innert2g)
           | 'band' (band-wise, sigma estimated from the data)
"""
import numpy as np


def _dx(X):   # X: (h,w,d) ; circular forward difference along dim 0
    return np.roll(X, -1, axis=0) - X


def _dxT(X):
    return np.roll(X, 1, axis=0) - X


def _dy(X):
    return np.roll(X, -1, axis=1) - X


def _dyT(X):
    return np.roll(X, 1, axis=1) - X


def soft(a, tau):
    return np.sign(a) * np.maximum(np.abs(a) - tau, 0.0)


def proj_l1_ball(x, beta):
    x = x.ravel()
    if beta <= 0:
        return np.zeros_like(x)
    ax = np.abs(x)
    if ax.sum() <= beta:
        return x
    u = np.sort(ax)[::-1]
    sv = np.cumsum(u)
    idx = np.arange(1, u.size + 1)
    act = u - (sv - beta) / idx > 0
    last = np.nonzero(act)[0]
    if last.size == 0:
        return np.zeros_like(x)
    k = last[-1]
    level = (sv[k] - beta) / (k + 1)
    return np.sign(x) * np.maximum(ax - level, 0)


def noise_sigma_vdiff(Y, dirs='v'):
    """Stripe-invariant per-band Gaussian sigma: MAD of vertical differences.
    Column-constant stripes cancel in Y[i+1]-Y[i]; MAD ignores sparse impulses.
    dirs='vh': also horizontal differences (row-constant stripes cancel there); take the smaller."""
    def one(d):
        d = d - np.median(d.reshape(-1, d.shape[2]), axis=0)
        return 1.4826 * np.median(np.abs(d).reshape(-1, d.shape[2]), axis=0) / np.sqrt(2.0)
    s = one(Y[1:] - Y[:-1])
    if dirs == 'vh':
        s = np.minimum(s, one(Y[:, 1:] - Y[:, :-1]))
    return s


def rcctv(Y, rk=8, lam=None, T_mode='none', G_mode='none', beta=0.0, sigma=0.0,
          cT=3.0, sig_online=False, scalar_sigma=False, T_every=1, stripe_dirs='v', tol=1e-5, max_iter=100, rho=1.5,
          max_mu=1e6, return_info=False):
    M, N, p = Y.shape
    MN = M * N
    D = Y.reshape(MN, p, order='F').astype(float)
    if lam is None:
        lam = 6.0 / np.sqrt(MN)
    norm_two = np.linalg.norm(D, 2)
    dual_norm = max(norm_two, np.abs(D).max())
    mu = 1.0 / norm_two
    normD = np.linalg.norm(D)
    sizeU = (M, N, rk)
    kx = 2 - 2 * np.cos(2 * np.pi * np.arange(M) / M)
    ky = 2 - 2 * np.cos(2 * np.pi * np.arange(N) / N)
    determ = (kx[:, None] + ky[None, :])[:, :, None]
    eps = 1e-4

    u, s, vt = np.linalg.svd(D, full_matrices=False)
    U = u[:, :rk] * s[:rk]
    V = vt[:rk].T
    E = np.zeros_like(D)
    T = np.zeros_like(D)
    G = np.zeros_like(D)
    M1 = np.zeros((MN, rk))
    M2 = np.zeros((MN, rk))
    M3 = D / dual_norm

    # dense-term weights
    nu = 3 * np.sqrt(2) / (np.sqrt(MN) + np.sqrt(p))          # Peng dual noise scale
    if G_mode == 'peng':
        sig_k = np.broadcast_to(np.asarray(sigma, float), (p,)).copy()
    elif G_mode == 'band':
        sig_k = noise_sigma_vdiff(Y, stripe_dirs)
        if scalar_sigma:
            sig_k = np.full(p, np.median(sig_k))
    else:
        sig_k = np.zeros(p)
    beta_g = (nu * (np.sqrt(MN) + np.sqrt(p))) / ((np.sqrt(MN) + np.sqrt(p)) * sig_k + 1e-6) \
        if G_mode != 'none' else np.full(p, np.inf)
    lamT = cT * nu / np.sqrt(M)

    a = np.zeros((N, p)); b = np.zeros((M, p))
    it = 0
    sig_hist = []
    while it < max_iter:
        it += 1
        A = _dx(U.reshape(sizeU, order='F')).reshape(MN, rk, order='F') + M1 / mu
        B, C, P = np.linalg.svd(A, full_matrices=False)
        G1 = (B * soft(C, 1 / mu)) @ P
        A = _dy(U.reshape(sizeU, order='F')).reshape(MN, rk, order='F') + M2 / mu
        B, C, P = np.linalg.svd(A, full_matrices=False)
        G2 = (B * soft(C, 1 / mu)) @ P

        # stripe update
        if T_mode == 'ball2':
            if beta <= 0:
                T = np.zeros_like(D)
            elif it % 2 == 0:
                R = D - U @ V.T - E - G + M3 / mu
                Tm = R.reshape(M, N, p, order='F').mean(0)           # (N,p)
                Tc = np.broadcast_to(Tm, (M, N, p)).reshape(MN, p, order='F')
                T = proj_l1_ball(Tc.ravel(order='F'), beta).reshape(MN, p, order='F')
        elif T_mode == 'pen' and it % T_every == 0:
            R = D - U @ V.T - E - G + M3 / mu
            R3 = R.reshape(M, N, p, order='F')
            if 'v' in stripe_dirs:        # column-constant stripes a (N,p): threshold cT*nu/sqrt(M)
                a = soft((R3 - b[:, None, :]).mean(0), cT * nu / np.sqrt(M) / mu)
            if 'h' in stripe_dirs:        # row-constant stripes b (M,p): threshold cT*nu/sqrt(N)
                b = soft((R3 - a[None, :, :]).mean(1), cT * nu / np.sqrt(N) / mu)
            T = (a[None, :, :] + b[:, None, :]).reshape(MN, p, order='F')

        # U
        dT = (_dxT((G1 - M1 / mu).reshape(sizeU, order='F')) +
              _dyT((G2 - M2 / mu).reshape(sizeU, order='F')))
        temp = (D - E - T - G + M3 / mu) @ V
        numer = dT + temp.reshape(sizeU, order='F')
        x = np.real(np.fft.ifft2(np.fft.fft2(numer, axes=(0, 1)) / (determ + 1 + eps), axes=(0, 1)))
        U = x.reshape(MN, rk, order='F')
        # V
        uu, _, vv = np.linalg.svd((D - E - T - G + M3 / mu).T @ U, full_matrices=False)
        V = uu @ vv
        # E
        X = U @ V.T
        E = soft(D - X - T - G + M3 / mu, lam / mu)
        # G
        if G_mode != 'none':
            if sig_online and it >= 5:
                r = (D - X - T)
                med = np.median(r, axis=0)
                sig_k = 1.4826 * np.median(np.abs(r - med), axis=0)
                if scalar_sigma:
                    sig_k = np.full(p, np.median(sig_k))
                sig_hist.append(sig_k.copy())
                beta_g = nu * (np.sqrt(MN) + np.sqrt(p)) / ((np.sqrt(MN) + np.sqrt(p)) * sig_k + 1e-6)
            G = (M3 + mu * (D - X - E - T)) / (beta_g[None, :] + mu)
        # multipliers
        leq1 = _dx(U.reshape(sizeU, order='F')).reshape(MN, rk, order='F') - G1
        leq2 = _dy(U.reshape(sizeU, order='F')).reshape(MN, rk, order='F') - G2
        leq3 = D - X - E - T - G
        if (np.linalg.norm(leq1) / normD < tol and np.linalg.norm(leq2) / normD < tol
                and np.linalg.norm(leq3) / normD < tol):
            break
        M1 = M1 + mu * leq1
        M2 = M2 + mu * leq2
        M3 = M3 + mu * leq3
        mu = min(max_mu, mu * rho)
    out = (U @ V.T).reshape(M, N, p, order='F')
    if return_info:
        return out, dict(iters=it, sig_k=sig_k, T=T, G=G, E=E, sig_hist=sig_hist, beta_g=beta_g, M3=M3, a=a, b=b)
    return out


# ---------------- scoring (port of sc() in protocol27/M/runM.m) ----------------
def _conv_valid(a, win):
    from scipy.signal import convolve2d
    return convolve2d(a, win, mode='valid')


def score(rec, clean):
    rec = np.clip(rec, 0, 1)
    ref = np.clip(clean, 0, 1)
    recb, refb = rec[:, :, 3:-3], ref[:, :, 3:-3]
    h, w, nb = recb.shape
    winsz = max(2, int(np.floor(min(h, w) / 16)))
    co = np.arange(winsz) - (winsz - 1) / 2
    ker = np.exp(-(co ** 2) / (2 * 1.5 ** 2))
    win = np.outer(ker, ker)
    win /= win.sum()
    c1, c2 = 0.01 ** 2, 0.03 ** 2
    ps, ss = [], []
    for b in range(nb):
        a, r = recb[:, :, b], refb[:, :, b]
        ps.append(10 * np.log10(1.0 / np.mean((a - r) ** 2)))
        m1, m2 = _conv_valid(a, win), _conv_valid(r, win)
        s1 = _conv_valid(a * a, win) - m1 * m1
        s2 = _conv_valid(r * r, win) - m2 * m2
        sc = _conv_valid(a * r, win) - m1 * m2
        st = ((2 * m1 * m2 + c1) * (2 * sc + c2)) / ((m1 * m1 + m2 * m2 + c1) * (s1 + s2 + c2))
        ss.append(st.mean())
    return float(np.mean(ps)), float(np.mean(ss))
