### Curvature filters are efficient solvers for variational models.
These curvature filters are slightly different from the original ones that were also developed by Yuanhao Gong during his PhD. But the theory is the same. Please cite following papers if you use curvature filter in your work. Thank you!

:books: **<a href="http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7835193" target="_blank">The Paper</a>**. The general theory is in Chapter **Six** of **<a href="http://e-collection.library.ethz.ch/eserv/eth:47737/eth-47737-02.pdf" target="_blank">PhD thesis</a>** (downloaded **3500+** between June, 2015 and Dec, 2016)

:closed_book: Presentation of Gaussian Curvature Filter: **<a href="http://www.slideshare.net/YuanhaoGong/a-fast-implicit-gaussian-curvature-filter" target=" blank">LinkedIn</a>**, **<a href="https://www.dropbox.com/s/ax73park0popi4x/GCFilter_small.pdf?dl=0" target="_blank">Dropbox</a>** or **<a href="https://pan.baidu.com/s/1geS2EXH" target="_blank">Baidu</a>**. 

:blue_book: Poster of Bernstein Filter can be found **[here](images/BernsteinFilter.pdf)**. 

:gift: a short introduction in Chinese: **<a href="https://zhuanlan.zhihu.com/p/22971865" target="_blank">Zhihu(Editors' Choice)</a>**, **<a href="http://www.zhihu.com/question/35499791" target="_blank">Zhihu</a>** or this **<a href="http://blog.csdn.net/jorg_zhao/article/details/51328966" target="_blank">Blog</a>**

:trophy: **source code** in **C++** and **Java** can also be found from **<a href="http://mosaic.mpi-cbg.de/?q=downloads/curvaturefilters" target="_blank">MOSAIC group</a>**

:bell: The kernels summary and one example how to get the kernel can be found **[here](images/CF_Kernels.pdf)**.

:e-mail: gongyuanhao@gmail.com or join the **<a href="https://groups.google.com/forum/?hl=en#!forum/curvaturefilter" target="_blank">Curvature Filter Forum</a>**
***
```text
@ARTICLE{gong:cf, 
    author={Yuanhao Gong and Ivo F. Sbalzarini}, 
    journal={IEEE Transactions on Image Processing}, 
    title={Curvature filters efficiently reduce certain variational energies}, 
    year={2017}, 
    volume={26}, 
    number={4}, 
    pages={1786-1798}, 
    doi={10.1109/TIP.2017.2658954}, 
    ISSN={1057-7149}, 
    month={April},}

@phdthesis{gong:phd, 
    title={Spectrally regularized surfaces}, 
    author={Gong, Yuanhao}, 
    year={2015}, 
    school={ETH Zurich, Nr. 22616},
    note={http://dx.doi.org/10.3929/ethz-a-010438292}}

@INPROCEEDINGS{gong:Bernstein,
    author={Yuanhao Gong}, 
    booktitle={2016 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)}, 
    title={Bernstein filter: A new solver for mean curvature regularized models}, 
    year={2016}, 
    pages={1701-1705}, 
    doi={10.1109/ICASSP.2016.7471967}, 
    month={March},}
	
@article{gong:gc,
    Author = {Yuanhao Gong and Ivo F. Sbalzarini},
    Journal = {Intl. Conf. Image Proc. (ICIP)},
    Month = {September},
    Pages = {534--538},
    Title = {Local weighted {G}aussian curvature for image processing},
    Year = {2013}}
```
***
## Curvature filters' philosophy 

Traditional solvers, such as gradient descent or Euler Lagrange Euqation, start at the total energy and use diffusion scheme to carry out the minimization. When the initial condition is the original image, the data fitting energy always increases while the regularization energy always reduces during the optimization, as illustrated in the below figure. Thus, regularization energy must be the dominant part since the total energy has to decrease. 

Therefore, **Curvature filters focus on minimizing the regularization term,** whose minimizers are already known. For example, if the regularization is Gaussian curvature, the developable surfaces minimize this energy. Therefore, in curvature filter, developable surfaces are used to approximate the data. **As long as the decreased amount in the regularization part is larger than the increased amount in the data fitting energy, the total energy is reduced.**

![image](images/phs.PNG)
***
## The role of curvature filters
Curvature filters link the variational models with image filters. Meanwhile, they implicitly impose differential geometry.                                                              ![ image ](images/role.png)
***
## Features
| Theoretical  | Practical |
| ------------- |:-------------:|
| **Generality**: handle arbitrary data fitting term (BlackBox) ![ image ](images/box.png) | **Efficient**: three or four order of magnitude faster than traditional solvers ![ image ](images/fast.jpg) |
| **Convergence**: theoretically guaranteed ![ image ](images/theory.png) | **Implementation**: 40 lines in Matlab and 100 lines in C++ ![ image ](images/easy.png) |

***
## Faster and Faster 
| Filter       | Bilateral Filter | Guided Filter | Guided Filter | MC Filter | GC Filter | Bernstein Filter |
| ------------- |:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
| Lang.      | C++ | Matlab | C++ | C++ | C++| C++|
| MilliSec.      | 103 | 514 | 130 | 8 (or **32.7 MPixels/sec**) | 11| 7|

Running time with 10 iterations on 512X512 Lena image. Matlab version is R2015a and GCC version is 5.1. All tests are on a Thinkpad T410 with i7-620M core CPU (2.6GHz). We take the time for 100 iterations and divide it by 10. On average, curvature filters take 1 millisecond per iteration.

On my new taptop(Thinkpad T470p, NVIDIA GeForce 940MX), GPU version of MC filter can achieve **1143 MPixels/Second** with single precision. Another implementation of MC filter, based on ArrayFire library, can achieve **88 MPixels/Second** with single precision. 

***
## Example Applications
### 1) Only minimize the regularization 
GC = Gaussian Curvature, MC = Mean Curvature, TV = Total Variation
![image](images/curvatureFilters.png)
### 2) Minimize a variational model, showing the line profile
We show three lines' profiles during minimizing a mean curvature regularized model (MC filter used). 

| ![ image](images/Lena_three_lines.png)      | ![image ](images/MC_line1_small.gif) |
| ------------- |:-------------:|
| ![image ](images/MC_line2_small.gif)      | ![image ](images/MC_line3_small.gif) |

### 3) Denoising
![image](images/denoise.PNG)
The noise free test image can be downloaded **[here](images/developable.png)**
### 4) Cartoon Texture Decomposition
![image](images/decomposition.png)
### 5) Registration
from left to right: original reference image, distorted source image, registered results by TV filter, MC filter and GC filter.
![image](images/lena_circ.png)
***
## On Triangular Meshes (preliminary results, p.195 in the thesis)
original mesh (left) and processed mesh (right), the energy profile is shown in the middle.
![image](images/GC_mesh.jpg)
***
## FAQ:
1) Why dual mesh (DM) structure is needed?

There are two reasons. First, these four sets guarantee the convergence. Second, 
we can use the updated neighbors for current position. Therefore, it is more computational efficient.

2) What is the difference between these three filters?

In general, GC filter is better in preserving details, compared with the other two. And
TV filter is better in removing noise as well as details. MC filter is between these two.

These three filters are correspond to three types of variational models. User should decide
which prior is to be assumed about the ground truth. 

3) What is the difference between split and nosplit scheme?

In general, splitting the image into four sets and looping on them is computational faster.
However, in some cases like deconvolution, we need to merge the four sets after every iteration.
So, it is better do nosplit scheme.

These two lead to exactly the same result. The split code is just more cache friendly.
