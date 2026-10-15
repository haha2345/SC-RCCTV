# Data sets (not included)

| path | content | source |
|---|---|---|
| `locked_nhsi/pavia120/`, `locked_nhsi/belt100/` | `clean.mat` and the official noisy inputs `case1_in.mat` ... `case8_in.mat` of RCCTV (Pavia University 120x120x99, Beltsville 100x100x128) | RCCTV test data, https://github.com/mendy-2013/RCCTV |
| `locked_nhsi/pavia64/`, `locked_nhsi/belt64/` | 64x64 crops of the same scenes (Python checks) | as above |
| `real/IndianPines/Indian_pines.mat`, `Indian_pines_gt.mat` | AVIRIS Indian Pines, 145x145x220, and ground truth | https://huggingface.co/datasets/danaroth/indian_pines |
| `real/Urban/Urban_F210.mat` | HYDICE Urban, 307x307x210 (variable `Urban_F210`) | https://huggingface.co/datasets/danaroth/urban |
| `large/PaviaU.mat` | ROSIS Pavia University, 610x340x103 (variable `paviaU`) | Hyperspectral Remote Sensing Scenes, University of the Basque Country (EHU) |

`experiments/inputs_manifest.csv` lists the SHA-256 of these files as used for the reported results;
`run_all_main.ps1 -Phase check` verifies them.
