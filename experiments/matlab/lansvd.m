function [U, S, V] = lansvd(A, K, SIGMA)
% Octave stand-in for PROPACK lansvd. Used only by author ctv_alm_spcp.m
% line 39:  norm_two = lansvd(D, 1, 'L');
% Same closed form: singular values of A via svd (largest if SIGMA='L').
% Not a new prox. Does not change SVT / softthre / G / FFT Poisson.
% Author ctv_alm_spcp.m is unmodified.

if nargin < 1
  error('lansvd: A required');
end
if nargin < 2 || isempty(K)
  K = 1;
end
if nargin < 3 || isempty(SIGMA)
  SIGMA = 'L';
end

[Ufull, Sfull, Vfull] = svd(A, 'econ');
s = diag(Sfull);
K = min(K, numel(s));
if upper(SIGMA(1)) == 'L'
  idx = 1:K;
else
  idx = (numel(s) - K + 1):numel(s);
end

if nargout <= 1
  U = s(idx);
  return
end
U = Ufull(:, idx);
S = diag(s(idx));
if nargout >= 3
  V = Vfull(:, idx);
end
end
