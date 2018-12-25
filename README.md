# Voronoi based Curve Reconstruction
An incremental labeling algorithm for curve reconstruction and medial axis extraction.

## Introduction
The repo contains the source code of the paper, Jiju Peethambaran, Amal Dev, Andrea Taglliasacchi, Ruisheng Wang and Ramanathan Muthuganapathy, "[Incremental labeling for Shape Reconstruction](https://onlinelibrary.wiley.com/doi/full/10.1111/cgf.13589#)", Computer Graphics Forum, 2018.
The algorithm is implemented using the Voronoi diagram adaptor available in [CGAL](https://www.cgal.org/project.html)(Computational Geometry Algorithms Library) 

## Abstract

An incremental Voronoi vertex labelling algorithm for approximating contours, medial axes and dominant points (high curvature points) from 2D point sets is presented. Though there exist many number of algorithms for reconstructing curves, medial axes or dominant points, a unified framework capable of approximating all the three in one place from points is missing in the literature. Our algorithm estimates the normals at each sample point through poles (farthest Voronoi vertices of a sample point) and uses the estimated normals and the corresponding tangents to determine the spatial locations (inner or outer) of the Voronoi vertices with respect to the original curve. The vertex classification helps to construct a piece‐wise linear approximation to the object boundary. We provide a theoretical analysis of the algorithm for points non‐uniformly (ε‐sampling) sampled from simple, closed, concave and smooth curves. The proposed framework has been thoroughly evaluated for its usefulness using various test data. Results indicate that even sparsely and non‐uniformly sampled curves with outliers or collection of curves are faithfully reconstructed by the proposed algorithm.

## Bibtex

@article{doi:10.1111/cgf.13589,
author = {Peethambaran, J. and Parakkat, A.D. and Tagliasacchi, A. and Wang, R. and Muthuganapathy, R.},
title = {Incremental Labelling of Voronoi Vertices for Shape Reconstruction},
journal = {Computer Graphics Forum},
volume = {0},
number = {0},
pages = {},
keywords = {curves and surfaces, modelling, computational geometry, geometric modelling, •Computing methodologies → Computer graphics; Shape analysis; •Theory of computation → Computational geometry},
doi = {10.1111/cgf.13589},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13589},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1111/cgf.13589},
}

## Teaser

![Picture](https://github.com/jijup/VoronoiCurve/blob/master/assets/teaser.jpg)

## Feedback

You may send your comments or feedback on the software to jijupnair2000@gmail.com
This is a build for unix environments. A windows build will be pushed soon.
