function [mpsnr, mssim, sam, ergas] = score_full(rec, clean)
% MPSNR/MSSIM exactly as score_hsi; SAM (deg) and ERGAS over all bands on clipped [0,1] cubes.
[mpsnr, mssim] = score_hsi(rec, clean);
a = min(max(double(rec),0),1); r = min(max(double(clean),0),1);
[M,N,p] = size(r); A = reshape(a,[M*N p]); B = reshape(r,[M*N p]);
cs = sum(A.*B,2) ./ max(sqrt(sum(A.^2,2)).*sqrt(sum(B.^2,2)), eps);
sam = mean(acosd(min(max(cs,-1),1)));
rm = sqrt(mean((A-B).^2,1)); mu = mean(B,1);
ergas = 100*sqrt(mean((rm./max(mu,eps)).^2));
end
