# EBM_RPCA
Code of TSP2019 Nonconvex Robust Low-Rank Tensor Reconstruction via an Empirical Bayes Method <br>

![pre_face](https://github.com/wc253/EBM_RPCA/pre.png)

# Folder structure
```shell
Demo_face_desparse.m          : Facial Image Denoising (salt-and-pepper noise).
EBM_RPCA.m                    : the main function of the proposed Empirical Bayes Method for tensor RPCA
Parset.m                      : the tuned parameters
data\
├────CMUface.mat              : the CMUPIE face database (30 subjects 11 view 21 illu imgsize 32*32)
lib\                              
├───quality_assess\         
├───show_results\   
├───tensor_toolbox\         : tensor processing toolbox http://www.sandia.gov/~tgkolda/TensorToolbox/index‐2.5.html
├───compete_methods\
├───────────────────HOrpca\     : https://sites.google.com/site/tonyqin/research
├───────────────────KBR\        : http://gr.xjtu.edu.cn/web/dymeng
├───────────────────TNNTRPCA\   : https://github.com/canyilu/Tensor-Robust-Principal-Component-Analysis-TRPCA
result\ 
```

# Citation
W. Chen, X. Gong and N. Song, "Nonconvex Robust Low-Rank Tensor Reconstruction via an Empirical Bayes Method," in IEEE Transactions on Signal Processing, vol. 67, no. 22, pp. 5785-5797, 15 Nov.15, 2019, doi: 10.1109/TSP.2019.2946022.

We would like to thank those researchers for making their codes and datasets publicly available. If you have any question, please feel free to contact me via: xiaogong@bjtu.edu.cn
