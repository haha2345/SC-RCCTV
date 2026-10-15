% Rank sensitivity of SC-RCCTV on real Indian Pines (OA of RBF-SVM, 10% train, 5 splits).  Report-only sensitivity.
R = setup_paths();
S = load(fullfile(R,'data','real','IndianPines','Indian_pines.mat')); Y = double(S.indian_pines);
for b = 1:size(Y,3), a = Y(:,:,b); mn = min(a(:)); mx = max(a(:)); Y(:,:,b) = (a-mn)/max(mx-mn,eps); end
G = load(fullfile(R,'data','real','IndianPines','Indian_pines_gt.mat')); gt = double(G.indian_pines_gt);
idx = find(gt>0); lab = gt(idx);
out = fullfile(res_dirs(R),'rank_sens.csv'); fid = fopen(out,'w'); fprintf(fid,'rk,OA,OA_std,AA,sec\n');
for rk = [4 8 12 16 24 32]
  t0 = tic; [X,~] = sc_rcctv(Y, struct('rk',rk,'lam',6/sqrt(145*145))); sec = toc(t0);
  F = reshape(X,[],size(X,3)); F = F(idx,:); OA = zeros(5,1); AA = zeros(5,1);
  for rep = 1:5
    rng(100+rep); tr = false(size(lab));
    for k = unique(lab)', ik = find(lab==k); n = max(5, round(0.1*numel(ik))); pr = ik(randperm(numel(ik),min(n,numel(ik)))); tr(pr) = true; end
    mdl = fitcecoc(F(tr,:), lab(tr), 'Learners', templateSVM('KernelFunction','rbf','KernelScale','auto','Standardize',true));
    pr = predict(mdl, F(~tr,:)); yt = lab(~tr); OA(rep) = mean(pr==yt); cm = confusionmat(yt,pr); AA(rep) = mean(diag(cm)./sum(cm,2));
  end
  fprintf(fid,'%d,%.3f,%.3f,%.3f,%.1f\n', rk, 100*mean(OA), 100*std(OA), 100*mean(AA), sec); fprintf('rk=%d OA=%.2f\n', rk, 100*mean(OA));
end
fclose(fid);
