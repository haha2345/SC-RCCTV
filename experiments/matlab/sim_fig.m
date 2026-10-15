% Visual comparison on a simulated cube: Pavia 120x120x99, case 8 (sigma .1 + 5% impulses + stripes), seed 0 (official input).
R = setup_paths_new(); warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory');
[base, inb] = res_dirs(R); fd = fullfile(base,'figs'); if ~exist(fd,'dir'), mkdir(fd); end
C = load(fullfile(R,'data','locked_nhsi','pavia120','clean.mat')); clean = C.clean; I = load(fullfile(inb,'inputs','pavia120_c8_s0.mat')); Y = I.Nhsi; [M,N,p] = size(Y);
b = 50;
Xr = RCCTV_RPCA(Y,8,6/sqrt(M*N)); Xs = sc_rcctv(Y,struct()); rng(0); Xe = EnhancedTV(Y,0.004*sqrt(M*N),[13 13 13]); rng(0); Xl = LRTDTV(Y,1,10,[round(.8*M) round(.8*N) 10]);
Xf = reshape(run_new_method('FastHyMix', Y, struct()), [M N p]); Xm = reshape(run_new_method('L1HyMixDe', Y, struct()), [M N p]);
list = {clean, Y, Xr, Xl, Xe, Xf, Xm, Xs}; ttl = {'Clean','Noisy','RCCTV','LRTDTV','E-3DTV','FastHyMix','L1HyMixDe','SC-RCCTV (ours)'};
nc = numel(list);
% physical size (screen-independent): nc square tiles across 18 cm, two rows plus a title line
wcm = 18; fig = figure('Visible','off','Units','centimeters','Position',[1 1 wcm 2*wcm/nc+0.45],'Color','w'); tiledlayout(2,nc,'TileSpacing','none','Padding','none');
for r = 1:2
  for c = 1:nc
    nexttile; Z = list{c}(:,:,b);
    if r == 1, imshow(min(max(Z,0),1),[0 1]); title(ttl{c},'FontSize',6.5,'FontWeight','normal');
      if c >= 2, cb = clean(:,:,b); zb = min(max(Z,0),1); text(3,11,sprintf('%.2f dB', 10*log10(1/mean((zb(:)-cb(:)).^2))),'Color',[1 .9 0],'FontSize',5.5,'FontWeight','bold','BackgroundColor',[0 0 0],'Margin',0.5); end
    else, E = abs(min(max(Z,0),1) - clean(:,:,b)); imshow(E,[0 0.25]); colormap(gca,hot); end
  end
end
exportgraphics(fig, fullfile(fd,'fig_sim_visual.png'),'Resolution',300); exportgraphics(fig, fullfile(fd,'fig_sim_visual.pdf')); close(fig);
disp('sim figure done');
