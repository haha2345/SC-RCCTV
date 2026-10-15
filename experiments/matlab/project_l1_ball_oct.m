function y = project_l1_ball_oct(x, beta)
%% Euclidean projection onto the vectorized L1 ball.
%% Port of iat_trpca.s3ttv_solver.project_l1_ball (vectorized).
%% If sum(abs(x))<=beta return x; if beta==0 return 0;
%% else sort |x| descending, find level, soft-threshold.
x = x(:);
if ~isfinite(beta) || beta < 0
  error('project_l1_ball_oct: beta must be finite and nonnegative');
end
if beta == 0
  y = zeros(size(x));
  return
end
ax = abs(x);
if sum(ax) <= beta
  y = x;
  return
end
u = sort(ax, 'descend');
sv = cumsum(u);
n = numel(u);
idx = (1:n).';
active = u - (sv - beta) ./ idx > 0;
last = find(active, 1, 'last');
if isempty(last)
  y = zeros(size(x));
  return
end
level = (sv(last) - beta) / last;
y = sign(x) .* max(ax - level, 0);
end
