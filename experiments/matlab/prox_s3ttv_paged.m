function z = prox_s3ttv_paged(z, gamma, blocksize)
% Batched (pagesvd) version of Prox_S3TTV.m from external/S3TTV: identical maths -- singular-value
% soft-thresholding of every non-overlapping b1 x b2 block (all bands and both gradient directions
% stacked as columns) for every periodic shift (k,l).  Processed one shift-row k at a time to bound memory.
[v, h, c, d, s1, s2] = size(z);
b1 = blocksize(1); b2 = blocksize(2); nb1 = v/b1; nb2 = h/b2;
for k = 1:s1
    A = z(:,:,:,:,k,:);                                   % v x h x c x d x 1 x s2
    A = reshape(A, [b1, nb1, b2, nb2, c, d, s2]);
    A = permute(A, [1 3 5 6 2 4 7]);                      % b1 b2 c d nb1 nb2 s2
    A = reshape(A, [b1*b2, c*d, nb1*nb2*s2]);
    [U, S, V] = pagesvd(A, 'econ');
    s = max(0, pagediag(S) - gamma);                      % r x 1 x pages
    A = pagemtimes(U .* reshape(s, 1, [], size(s,3)), 'none', V, 'transpose');
    A = reshape(A, [b1, b2, c, d, nb1, nb2, s2]);
    A = ipermute(A, [1 3 5 6 2 4 7]);
    z(:,:,:,:,k,:) = reshape(A, [v, h, c, d, 1, s2]);
end
end

function s = pagediag(S)
r = min(size(S,1), size(S,2)); n = size(S,3);
idx = (0:r-1)*size(S,1) + (1:r);                          % linear indices of the diagonal in one page
S = reshape(S, [], n);
s = reshape(S(idx, :), [r, 1, n]);
end
