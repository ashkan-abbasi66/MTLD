# Multiscale Sparsifying Transform Learning for Image Denoising

This repository contains the code associated with the following paper:

Ashkan Abbasi, Amirhassan Monadjemi, Leyuan Fang, Hossein Rabbani, Neda Noormohammadi, Yi Zhang, "Multiscale Sparsifying Transform Learning for Image Denoising," [ArXiv Prepr.](https://arxiv.org/abs/2003.11265), 2020.



## Available Methods

The code for running the following methods are available in this package. All the codes needed to run the methods are included except for the sparsifying transform learning denoising. 



- Four TLD (Sparsifying Transform Learning Denoising) based methods: **TLD, MTLD, MMTLD, FMMTLD**
  - Demo scripts: `Benchmark_MTLD_for_Gaussian_denoising.m` and `Benchmark_MTLD_for_FMD.m`
  - Requirements:
    - Download sparsifying transform learning [1-3] package (`TSP2015ClosedformTL_code.zip`) from [here](http://transformlearning.csl.illinois.edu/software/).
    - Extract the package in the `METHODS` folder. So, the path should be like this: `./METHODS/TSP2015ClosedformTL_code`



- **MM K-SVD** (Multiscale Mixed K-SVD) and **Fast MM K-SVD**
  - Demo script: `Benchmark_MMKSVD_for_Gaussian_denoising.m`



- Fused K-SVD package contains **K-SVD, MS K-SVD**, and **Fused K-SVD**
  - `Benchmark_FusedKSVD_for_Gaussian_denoising.m` and `Benchmark_FusedKSVD_for_FMD.m`



- **BLS-GSM**
  - `Benchmark_BLS_GSM_for_Gaussian_denoising.m`
  - Instead of BLS-GSM, we use PURE-LET for denoising fluorescence microscopy images.



For PURE-LET, please refer to here.

All of the codes run over a Windows operating system with a proper MATLAB installation. However, just  to let you know, we carried out our experiments using MATLAB R2019a over Windows 10. 

Please, feel free to contact me or open an issue. 



## Datasets

**(1)** 12 Classic test images which are stored [here](./DATASETS/20_classic_images).  

**(2)** The image `test_011` from BSD68 dataset. 

**(3)** The `test_mix` subset of the Florescence Microscopy Denoising (FMD) dataset [1] ([GitHub page](https://github.com/yinhaoz/denoising-fluorescence), and [Download link](https://drive.google.com/drive/folders/1aygMzSDdoq63IqSk-ly8cMq0_owup8UM) (Google Drive)).



Given a folder containing the images (e.g., `./DATASETS/20_classic_images`), you can use **`Generate_synthetic_Gaussian_noise.m`** to generate noisy version of those images. For each given noise level, this script saves the noisy image and original image in separate `.mat` files and place them into a subfolder. 

By storing noisy images in `.mat` files, we can ensure that our experiments are as repeatable as possible. 



## References

[1] Y. Zhang et al., “A Poisson-Gaussian Denoising Dataset With Real Fluorescence Microscopy Images,” in IEEE Conference on Computer Vision and Pattern Recognition, 2019, pp. 11702–11710. 