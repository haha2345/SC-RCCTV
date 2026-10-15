%% Demo of SC-RCCTV on a synthetic hyperspectral cube (no data download needed).
% Needs MATLAB R2016b or newer with the Image Processing Toolbox (psf2otf, imgaussfilt); sc_rcctv.m has no other
% dependency. The cube is low-rank and spatially smooth, and is corrupted by band-dependent Gaussian noise,
% 5% salt-and-pepper impulses and stripes on 5% of the columns of every band. No noise parameter is passed.
addpath(fullfile(fileparts(mfilename('fullpath')), 'experiments', 'matlab'));
rng(0);
M = 100; N = 100; p = 80; r = 4;
U = zeros(M*N, r);                                       % smooth, piecewise-constant coefficient images
for i = 1:r, U(:, i) = reshape(imgaussfilt(double(rand(M, N) > 0.5), 3), [], 1); end
t = linspace(0, 1, p)';                                  % smooth spectra
V = [exp(-(t - 0.2).^2 / 0.02), exp(-(t - 0.5).^2 / 0.05), 0.5 + 0.5 * sin(4 * t), t];
X = reshape(U * V', M, N, p); X = (X - min(X(:))) / (max(X(:)) - min(X(:)));

sig = 0.02 + 0.13 * rand(1, 1, p);                       % band-dependent Gaussian noise levels
Y = X + sig .* randn(M, N, p);
imp = rand(M, N, p) < 0.05; Y(imp) = double(rand(nnz(imp), 1) > 0.5);         % impulses
for k = 1:p                                              % column stripes
  cols = randperm(N, round(0.05 * N));
  Y(:, cols, k) = Y(:, cols, k) + (rand(1, numel(cols)) - 0.5);
end

tic; [Xh, info] = sc_rcctv(Y); sec = toc;              % default settings, blind
psnr_k = @(Z, k) 10 * log10(1 / mean(reshape(min(max(Z(:, :, k), 0), 1) - X(:, :, k), [], 1) .^ 2));
mpsnr = @(Z) mean(arrayfun(@(k) psnr_k(Z, k), 1:p));
fprintf('MPSNR: noisy %.2f dB -> SC-RCCTV %.2f dB (%.1f s, %d iterations)\n', mpsnr(Y), mpsnr(Xh), sec, info.iters);
fprintf('noise levels: median |estimated sigma_k - true sigma_k| = %.4f (true range %.2f-%.2f)\n', ...
        median(abs(info.sig_k(:) - sig(:))), min(sig), max(sig));
