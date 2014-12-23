GPU-SExtractor
==============

Parallel Astronomical Source Extraction tool based on SExtractor

==========

Main contribution

In this work, we propose to use the GPU (Graphics Processing Unit) to accelerate source extraction. Our work is based on SExtractor, an astronomical source extraction tool widely used in astronomy projects, and study its parallelization on the GPU. In GPU-SExtractor, we re-design and parallelize each major step in SExtractor:

1) Background Computation,

2) Multi-Threshold Object Detection,

3) Object Cleaning and

4) Object Analysis.

In particular, we identify the Multi-Threshold Object Detection step as the most complex and time-consuming, and design a parallel detection algorithm based on Connected-Component Labelling. Furthermore, we apply compaction techniques to optimize the detection algorithm to
better utilize the massive GPU thread parallelism. In the Object Analysis step, we decompose the analysis into a sequence of GPU-friendly
data-parallel primitives to compute the attributes of each extracted object. We have evaluated our GPU-SExtractor in comparison with the original
SExtractor on a desktop with an Intel i7 CPU and an NVIDIA GTX670 GPU on a set of real-world and synthetic astronomical images of different sizes.
The results show that our GPU-SExtractor outperforms the original SExtractor by a factor of 6, taking a merely 1.9 second to process a typical
4KX4K image containing 167,000 celestial objects.


==========

File description

The CUDA code is in the src/cuda directory.

cudaback.cu(.h) Parallel Background Computation

cudaclean.cu(.h) Parallel Cleaning

cudadetection.cu(.h) Parallel Raw Object Detection

cudadeblend.cu(.h) Parallel Multi-level Object Deblending

cudaanalyse.cu(.h) Parallel Object Analysis

==========

Compile

Software Dependency:

Besides the software required for SExtractor (described in "./doc/SExtractor installation - MediaWiki.html"), the following SDKs and libraries are
required

1) CUDA version 6.5 or higher

2) CUDPP version 2.2

Compile steps:

1) install the dependent softwares following the instructions in "./doc/SExtractor installation - MediaWiki.html",

2) run "./configure --with-atlas=/usr/local/atlas/lib --with-atlas-incdir=/usr/local/atlas/include --with-fftw=/usr/local/lib"

3) Replace the content in ./src/Makefile with ./src/Makefile.bak

3) run make, there should be some errors

4) execute rebuild.sh script to patch the compilation

5) run make install to install the executable


==========

Sample Usage:

sex -c param.sex image4k4k.fits

==========

Citation:
Baoxue Zhao,  Qiong Luo, Chao Wu. "Parallelizing Astronomical Source Extraction on the GPU." eScience (eScience), 2013 IEEE 9th International Conference on. IEEE, 2013.
