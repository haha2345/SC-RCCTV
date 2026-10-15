function [HSI_restored, iteration] = s3ttv_cpu(HSI_noisy, params)
%% CPU copy of S3TTV_CPU_fast.m (Takemoto, Naganuma, Ono, JSTARS 2025; MIT licence,
%  external/S3TTV).  Only changes: the stray gpuArray(disprate) call is removed so it runs
%  without a GPU, printing is silenced, and the block-SVT prox is batched with pagesvd
%  (prox_s3ttv_paged.m, checked equal to Prox_S3TTV to 2e-6 relative).  Algorithm, step sizes
%  and stopping rule unchanged.
%  params.use_gpu = true runs on the GPU (as the authors' S3TTV_GPU_fast.m does); the batched prox is
%  used when pagesvd accepts gpuArray, otherwise the official per-block Prox_S3TTV.
use_gpu = isfield(params,'use_gpu') && params.use_gpu;
HSI_noisy = single(HSI_noisy);
[n1, n2, n3] = size(HSI_noisy);
alpha = single(params.alpha); beta = single(params.beta); epsilon = single(params.epsilon);
blocksize = double(params.blocksize); maxiter = params.maxiter; stopcri = single(params.stopcri);
U = zeros([n1, n2, n3], 'single'); S = U; T = U;
Y1 = zeros([n1, n2, n3, 2, blocksize(1), blocksize(2)], 'single');
proxf = @prox_s3ttv_paged;
if use_gpu
    HSI_noisy = gpuArray(HSI_noisy); U = gpuArray(U); S = U; T = U; Y1 = gpuArray(Y1);
    try, pagesvd(gpuArray(rand(3,3,2,'single'))); catch, proxf = @Prox_S3TTV; end
end
Y2 = U; Y3 = U; Y4 = U; Y5 = U;
D   = @(z) cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z);
Dt  = @(z) z([end,1:end-1],:,:,1) - z(:,:,:,1) + z(:,[end,1:end-1],:,2) - z(:,:,:,2);
Dv  = @(z) z([2:end, 1],:,:) - z;
Dvt = @(z) z([end,1:end-1],:,:) - z(:,:,:);
Ds  = @(z) z(:, :, [2:end, 1], :) - z;
Dst = @(z) z(:,:,[end,1:end-1],:) - z(:,:,:,:);
P  = @(z) func_PeriodicExpansion(z, blocksize);
Pt = @(z) func_PeriodicExpansionTrans(z);
gamma1_U = single(1./(prod(blocksize) * 2*2 * 2 + 1)); gamma1_S = single(1); gamma1_T = single(1/(2 + 1));
gamma2_Y1 = single(1/(2*2)); gamma2_Y2 = single(1); gamma2_Y3 = single(1); gamma2_Y4 = single(1); gamma2_Y5 = single(1/2);
for i = 1:maxiter
    U_tmp = U - gamma1_U.*(Dst(Dt(Pt(Y1))) + Y2);
    S_tmp = S - gamma1_S.*Y3;
    T_tmp = T - gamma1_T.*(Y4 + Dvt(Y5));
    Primal_sum = U_tmp + S_tmp + T_tmp;
    Primal_sum = ProjL2ball(Primal_sum, HSI_noisy, epsilon) - Primal_sum;
    U_next = U_tmp + Primal_sum/3; S_next = S_tmp + Primal_sum/3; T_next = T_tmp + Primal_sum/3;
    U_res = 2*U_next - U; S_res = 2*S_next - S; T_res = 2*T_next - T;
    Y1_tmp  = Y1 + gamma2_Y1.*(P(D(Ds(2*U_next - U))));
    Y1 = Y1_tmp - gamma2_Y1.*proxf(Y1_tmp./gamma2_Y1, 1./gamma2_Y1, blocksize);   % = Prox_S3TTV (batched on CPU)
    Y2_tmp = Y2 + gamma2_Y2.*U_res;  Y2 = Y2_tmp - gamma2_Y2*ProjBox(Y2_tmp/gamma2_Y2, 0, 1);
    Y3_tmp = Y3 + gamma2_Y3.*S_res;  Y3 = Y3_tmp - gamma2_Y3.*ProjFastL1Ball(Y3_tmp./gamma2_Y3, alpha);
    Y4_tmp = Y4 + gamma2_Y4.*T_res;  Y4 = Y4_tmp - gamma2_Y4*ProjFastL1Ball(Y4_tmp/gamma2_Y4, beta);
    Y5 = Y5 + gamma2_Y5.*Dv(T_res);
    rate = norm(U_next(:) - U(:),2)/norm(U(:),2);
    U = U_next; S = S_next; T = T_next;
    if isfield(params,'trace_clean') && mod(i, params.trace_every) == 0     % convergence probe only
        fprintf('S3TTV_TRACE it %d rate %.3g mpsnr %.3f\n', i, rate, score_full(double(gather(U)), params.trace_clean));
    end
    if i >= 2 && rate < stopcri, break; end
end
HSI_restored = double(gather(U)); iteration = i;
end
