function [mpsnr, mssim] = score_hsi(rec, clean)
% Same formula as sc() in protocol27/M/runM.m: clip to [0,1], drop 3 edge bands each side,
% PSNR peak 1, Gaussian-window SSIM (win = max(2,floor(min(h,w)/16)), sigma 1.5, 'valid').
rec = min(max(double(rec),0),1); ref = min(max(double(clean),0),1);
recb = rec(:,:,4:end-3); refb = ref(:,:,4:end-3); [h,w,nb] = size(recb);
psnrs = zeros(nb,1); ssims = zeros(nb,1); winsz = max(2, floor(min(h,w)/16));
coords = (0:winsz-1) - (winsz-1)/2; kern = exp(-(coords.^2)/(2*1.5^2));
win = kern'*kern; win = win/sum(win(:)); c1 = 0.01^2; c2 = 0.03^2;
for b = 1:nb
    a = recb(:,:,b); r = refb(:,:,b); mse = mean((a(:)-r(:)).^2); psnrs(b) = 10*log10(1/mse);
    mu1 = conv2(a,win,'valid'); mu2 = conv2(r,win,'valid');
    mu1sq = mu1.*mu1; mu2sq = mu2.*mu2; mup = mu1.*mu2;
    s1 = conv2(a.*a,win,'valid')-mu1sq; s2 = conv2(r.*r,win,'valid')-mu2sq; scv = conv2(a.*r,win,'valid')-mup;
    st = ((2*mup+c1).*(2*scv+c2))./((mu1sq+mu2sq+c1).*(s1+s2+c2)); ssims(b) = mean(st(:));
end
mpsnr = mean(psnrs); mssim = mean(ssims);
end
