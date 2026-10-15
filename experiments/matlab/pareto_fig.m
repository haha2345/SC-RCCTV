% Accuracy vs. run time (log) -- (a) native size, case 8 seed 0, mean of both cubes;
% (b) full Pavia University 610x340x103, mean of cases 8 and 10; (c) Indian Pines OA.
% Times: idle-machine wall clock (timing.csv, timing_new.csv, large.csv).  Ours highlighted, others grey.
R = setup_paths(); base = res_dirs(R); fd = fullfile(base,'figs');
M  = {'RCCTV','LRTDTV','E3DTV','tCTVd','FastHyMix','L1HyMixDe','S3TTV','SC'};
LB = {'RCCTV','LRTDTV','E-3DTV','t-CTV','FastHyMix','L1HyMixDe','S3TTV*','SC-RCCTV'};
MK = {'o','s','d','v','^','>','p','o'};
ink = [82 81 78]/255; ours = [42 120 214]/255;
% ---- data
fl = dir(fullfile(base,'sim_*.csv')); S = table();
for i = 1:numel(fl)
  T = readtable(fullfile(base,fl(i).name)); vn = T.Properties.VariableNames;
  vn(ismember(vn,{'case','case_','xCase'})) = {'cs'}; T.Properties.VariableNames = vn;   % 'case' is a keyword; readtable renames it in some releases
  if all(ismember({'cs','meth'},vn)), S = [S; T(:,{'cube','cs','seed','meth','mpsnr'})]; end %#ok<AGROW>
end
S = S(S.cs==8 & S.seed==0,:);
T1 = readtable(fullfile(base,'timing.csv')); T2 = readtable(fullfile(base,'timing_new.csv'),'Delimiter',',');
TT = [T1(:,{'ds','meth','sec'}); T2(:,{'ds','meth','sec'})];
TT.meth(strcmp(TT.meth,'tCTV')) = {'tCTV_default'};
L  = readtable(fullfile(base,'large.csv')); L = L(startsWith(L.job,'paviaU'),:);
C  = readtable(fullfile(base,'real_indian_cls.csv'));
P = nan(numel(M),2,3);                                   % method x (time, acc) x panel
for k = 1:numel(M)
  m = M{k};
  v = S.mpsnr(strcmp(S.meth,m)); tp = TT.sec(strcmp(TT.meth,m) & ismember(TT.ds,{'pavia120_c8','belt100_c8'}));
  if numel(v)==2 && numel(tp)==2, P(k,:,1) = [mean(tp) mean(v)]; end
  r = L(strcmp(L.meth,m),:);
  if height(r)==2 && all(isfinite(r.mpsnr)), P(k,:,2) = [mean(r.sec) mean(r.mpsnr)]; end
  ti = TT.sec(strcmp(TT.meth,m) & strcmp(TT.ds,'indian220')); oa = C.OA(strcmp(C.meth,m));
  if ~isempty(ti) && ~isempty(oa), P(k,:,3) = [ti(1) oa(1)]; end
end
save(fullfile(base,'pareto_points.mat'),'P','M');
ttl = {'(a) 120\times120 / 100\times100, case 8','(b) 610\times340\times103, cases 8 & 10','(c) Indian Pines (real)'};
yl  = {'MPSNR (dB)','MPSNR (dB)','OA (%)'};
% label placement where neighbours collide (default: right of the marker)
LP = {struct('SC','below','RCCTV','below'), ...
      struct('SC','above'), ...
      struct('E3DTV','left')};
fig = figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.0 2.15]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for pnl = 1:3
  ax = nexttile; hold(ax,'on'); set(ax,'XScale','log','FontSize',7.5,'Box','off','TickDir','out', ...
      'XColor',ink,'YColor',ink,'GridColor',[.85 .85 .85],'GridAlpha',1); grid(ax,'on'); ax.XMinorGrid = 'off';
  xs = P(:,1,pnl); xs = xs(isfinite(xs)); if ~isempty(xs), xlim(ax,[min(xs)/2.2, max(xs)*6]); end
  ys = P(:,2,pnl); ys = ys(isfinite(ys)); if ~isempty(ys), ylim(ax,[min(ys)-0.6, max(ys)+0.6]); end
  dy = 0.035*diff(ylim(ax));
  for k = 1:numel(M)
    if any(isnan(P(k,:,pnl))), continue; end
    isours = strcmp(M{k},'SC'); c = ink; ms = 5; if isours, c = ours; ms = 7; end
    plot(ax, P(k,1,pnl), P(k,2,pnl), MK{k}, 'MarkerSize', ms, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
    pos = 'right'; if isfield(LP{pnl}, M{k}), pos = LP{pnl}.(M{k}); end
    switch pos
      case 'above', tx = {P(k,1,pnl), P(k,2,pnl)+dy, 'HorizontalAlignment','center','VerticalAlignment','bottom'};
      case 'below', tx = {P(k,1,pnl), P(k,2,pnl)-dy, 'HorizontalAlignment','center','VerticalAlignment','top'};
      case 'left',  tx = {P(k,1,pnl)/1.18, P(k,2,pnl), 'HorizontalAlignment','right','VerticalAlignment','middle'};
      otherwise,    tx = {P(k,1,pnl)*1.18, P(k,2,pnl), 'VerticalAlignment','middle'};
    end
    text(ax, tx{1}, tx{2}, LB{k}, tx{3:end}, 'FontSize', 6.5, 'Color', c*0.9, 'FontWeight', ternary(isours,'bold','normal'));
  end
  xlabel(ax,'Run time (s, log scale)'); ylabel(ax,yl{pnl}); title(ax,ttl{pnl},'FontWeight','normal','FontSize',8);
end
exportgraphics(fig, fullfile(fd,'fig_pareto.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile(fd,'fig_pareto.png'), 'Resolution', 300); close(fig);
disp(squeeze(P(:,:,1))); disp(squeeze(P(:,:,2))); disp(squeeze(P(:,:,3)));

function r = ternary(c,a,b), if c, r = a; else, r = b; end, end
