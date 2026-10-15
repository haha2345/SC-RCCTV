# Third-party notices

`experiments/matlab/s3ttv_cpu.m` and `experiments/matlab/prox_s3ttv_paged.m` are derived from S3TTV (S. Takemoto, K. Naganuma, S. Ono, "Spatio-spectral structure tensor total variation for hyperspectral image denoising and destriping," IEEE JSTARS, 2025; https://www.mdi.c.titech.ac.jp/publications/s3ttv). Changes: CPU execution without the stray `gpuArray` call, silenced printing, and a batched (`pagesvd`) version of the block singular-value thresholding, checked equal to the official `Prox_S3TTV` to 2e-6 relative error. The S3TTV license:

```text
MIT License

Copyright (c) 2024 MDILab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

No other third-party code is included. The comparison methods (RCCTV, t-CTV, LRTDTV, E-3DTV, FastHyMix, L1HyMixDe, S3TTV, BM3D) and the data sets must be obtained from their authors; see `code/README.md`, `experiments/external/README.md` and `data/README.md`.

