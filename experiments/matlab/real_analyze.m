% Analyse real-data outputs: stripe index, level shift, residual noise, Indian Pines SVM accuracy, figures.
R = setup_paths();
base = res_dirs(R); rd = fullfile(base,'real'); fd = fullfile(base,'figs'); if ~exist(fd,'dir'), mkdir(fd); end
meths = {'Input','RCCTV','E3DTV','LRTDTV','tCTV','tCTVd','FastHyMix','L1HyMixDe','SC','SC_noT','SC_noG','SC_scalar','SC_fix'};
dss = {'indian220','urban200'};
rows = {};
for di = 1:2
  ds = dss{di}; Y = load_real(R, ds); [M,N,p] = size(Y);
  sgY = vdiff_sigma(Y); [~,ordS] = sort(sgY,'descend');
  SIy = stripe_index(Y); [~,ordI] = sort(SIy,'descend');
  for mi = 1:numel(meths)
    m = meths{mi};
    if strcmp(m,'Input'), X = Y; sec = NaN; else
      f = fullfile(rd,[ds '_' m '.mat']); if ~exist(f,'file'), continue; end
      L = load(f); X = L.X; sec = L.sec;
    end
    si = stripe_index(X); sg = vdiff_sigma(X);
    lv = squeeze(mean(mean(X - Y,1),2));                       % per-band mean shift (Y already in [0,1])
    topS = ordS(1:20);
    rows(end+1,:) = {ds, m, mean(si), mean(si(topS)), mean(sg), mean(sg(topS)), max(abs(lv)), mean(abs(lv)), sec}; %#ok<SAGROW>
    if di == 1, X1{mi} = X; end                                 %#ok<SAGROW>
    S_.(ds).(m) = X;
  end
  S_.(ds).Y = Y; S_.(ds).ordS = ordS; S_.(ds).ordI = ordI; S_.(ds).sgY = sgY;
end
T = cell2table(rows, 'VariableNames', {'ds','meth','stripeIdx_all','stripeIdx_top20','noiseMAD_all','noiseMAD_top20','shift_max','shift_mean','sec'});
writetable(T, fullfile(base,'real_metrics.csv')); disp(T);

%% Indian Pines classification (10% train/class, RBF-SVM, 5 splits)
G = load(fullfile(R,'data','real','IndianPines','Indian_pines_gt.mat')); gt = double(G.indian_pines_gt);
idx = find(gt>0); lab = gt(idx); acc = {};
for mi = 1:numel(meths)
  m = meths{mi};
  if ~isempty(getenv('SKIP_CLS')), break; end
  if ~isfield(S_.indian220, m) && ~strcmp(m,'Input'), continue; end
  if strcmp(m,'Input'), X = S_.indian220.Y; else, X = S_.indian220.(m); end
  F = reshape(X,[],size(X,3)); F = F(idx,:);
  OA = zeros(5,1); AA = zeros(5,1); KA = zeros(5,1);
  for rep = 1:5
    rng(100+rep); tr = false(size(lab));
    for k = unique(lab)', ik = find(lab==k); n = max(5, round(0.1*numel(ik))); pr = ik(randperm(numel(ik),min(n,numel(ik)))); tr(pr) = true; end
    mdl = fitcecoc(F(tr,:), lab(tr), 'Learners', templateSVM('KernelFunction','rbf','KernelScale','auto','Standardize',true));
    pr = predict(mdl, F(~tr,:)); yt = lab(~tr);
    OA(rep) = mean(pr==yt); cm = confusionmat(yt, pr); AA(rep) = mean(diag(cm)./sum(cm,2));
    pe = sum(sum(cm,1).*sum(cm,2)')/sum(cm(:))^2; KA(rep) = (OA(rep)-pe)/(1-pe);
  end
  acc(end+1,:) = {m, 100*mean(OA), 100*std(OA), 100*mean(AA), mean(KA)}; %#ok<SAGROW>
end
if ~isempty(acc), A = cell2table(acc, 'VariableNames', {'meth','OA','OA_std','AA','kappa'}); writetable(A, fullfile(base,'real_indian_cls.csv')); disp(A); end

%% Figures
cols = {'Input','RCCTV','E3DTV','LRTDTV','FastHyMix','L1HyMixDe','SC'}; ttl = {'Input','RCCTV','E-3DTV','LRTDTV','FastHyMix','L1HyMixDe','SC-RCCTV (ours)'};
% Indian Pines: noisiest band outside the water-absorption windows (104-108, 150-163, 220; the usual 200-band set).
% Urban: band 207, heavily striped but not dead. Near-dead bands (e.g. Indian 151, Urban 147) carry almost no signal
% and are discussed in the text instead (2026-10-14; previously the 10th noisiest band, which fell into these windows).
oI = S_.indian220.ordS; oI = oI(~ismember(oI, [104:108 150:163 220]));
pick = {{'indian220', oI(1)}, {'urban200', 207}};
nc = numel(cols);
wcm = 18; fig = figure('Visible','off','Units','centimeters','Position',[1 1 wcm 2*wcm/nc+0.45],'Color','w');   % physical size, screen-independent
tl = tiledlayout(2,nc,'TileSpacing','none','Padding','none');
for r = 1:2
  ds = pick{r}{1}; b = pick{r}{2};
  ref = cell2mat(cellfun(@(m) reshape(S_.(ds).(m)(:,:,b),[],1), cols(2:end), 'UniformOutput', false));   % common range from all denoised outputs
  lo = prctile(ref(:),1); hi = prctile(ref(:),99);
  for c = 1:nc
    nexttile; if strcmp(cols{c},'Input'), Z = S_.(ds).Y(:,:,b); else, Z = S_.(ds).(cols{c})(:,:,b); end
    imshow(Z,[lo hi]); if r==1, title(ttl{c},'FontSize',6.5,'FontWeight','normal'); end
    if c==1, text(4,14,sprintf('%s, band %d', strrep(strrep(ds,'220',''),'200',''), b),'Color',[1 .9 0],'FontSize',5.5,'FontWeight','bold','BackgroundColor',[0 0 0],'Margin',0.5); end
  end
end
exportgraphics(fig, fullfile(fd,'fig_real_visual.png'), 'Resolution', 300); exportgraphics(fig, fullfile(fd,'fig_real_visual.pdf'));
close(fig);
% sigma_k curve (Indian Pines): the initial estimate is the one sc_rcctv starts from, min of the vertical- and
% horizontal-difference MADs (stripe_dirs = 'vh'); before 2026-10-14 the curve showed the vertical MAD only
L = load(fullfile(rd,'indian220_SC.mat')); sk = L.info.sig_k;
Yi = S_.indian220.Y; sy = min(vdiff_sigma(Yi), vdiff_sigma(permute(Yi,[2 1 3])));
fig = figure('Visible','off','Units','centimeters','Position',[1 1 8.8 4.2],'Color','w');
semilogy(sy,'-','Color',[.6 .6 .6],'LineWidth',0.8); hold on; semilogy(sk,'-','Color',[0.85 .2 .1],'LineWidth',1.1);
xlabel('Band'); ylabel('\sigma_k'); legend({'initial estimate','after online refinement'},'Location','northwest','Box','off'); grid on; xlim([1 numel(sk)]);
set(gca,'FontSize',7); exportgraphics(fig, fullfile(fd,'fig_sigma_curve.pdf')); exportgraphics(fig, fullfile(fd,'fig_sigma_curve.png'),'Resolution',300); close(fig);
% statistics quoted in the text (make_rsl.py reads sigma_stats.csv and real_picks.csv)
wab = [104:108 150:163 220]; ok = setdiff(1:numel(sk), wab);
r = sy(ok)./sk(ok); up2 = ok(sk(ok) > 2*sy(ok)); dev = max(abs(sk(wab)./sy(wab) - 1));
fid = fopen(fullfile(base,'sigma_stats.csv'),'w');
fprintf(fid,'med_ratio,max_ratio,abs_maxdev,bands_up2\n%.4f,%.4f,%.4f,%s\n', median(r), max(r), dev, strjoin(arrayfun(@num2str, up2, 'UniformOutput', false), ' '));
fclose(fid);
fid = fopen(fullfile(base,'real_picks.csv'),'w'); fprintf(fid,'ds,band\n%s,%d\n%s,%d\n', pick{1}{1}, pick{1}{2}, pick{2}{1}, pick{2}{2}); fclose(fid);
save(fullfile(base,'real_picks.mat'), 'pick');
disp('picked bands:'); disp(pick);

function s = vdiff_sigma(Y)
[M,N,p] = size(Y); d = reshape(diff(Y,1,1),[(M-1)*N p]); d = d - median(d,1); s = 1.4826*median(abs(d),1)/sqrt(2);
end
function si = stripe_index(Y)
[M,N,p] = size(Y); c = squeeze(mean(Y,1)); r = squeeze(mean(Y,2));   % column means (N x p), row means (M x p)
si = zeros(1,p);                                                  % orientation-agnostic: roughness of both profiles
for b = 1:p, v = c(:,b); w = r(:,b); si(b) = mean(abs(v - movmedian(v,7))) + mean(abs(w - movmedian(w,7))); end
end
function Y = load_real(R, ds)
switch ds
  case 'indian220', S = load(fullfile(R,'data','real','IndianPines','Indian_pines.mat')); Y = double(S.indian_pines);
  case 'urban200',  S = load(fullfile(R,'data','real','Urban','Urban_F210.mat')); Y = double(S.Urban_F210(54:253,54:253,:));
end
for b = 1:size(Y,3), a = Y(:,:,b); mn = min(a(:)); mx = max(a(:)); Y(:,:,b) = (a-mn)/max(mx-mn,eps); end
end
