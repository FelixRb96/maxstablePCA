
# maxstablePCA

A package for dimensionality reduction of multivariate extremes using the idea of PCA to 
obtain a resonable compact description of the data. 

### Main functionalities

* Transform a dataset to standard margins to use well known ideas from extreme value theory
* Perform a dimensionality reduction of a dataset to a fixed number of encoding variables. For further information about the theory of this consider looking at the references.
* Evaluate the quality of this reconstruction.
* Transform the data back to the distribution of the original dataset.

### Examples on simulated and real world data

For a better feeling of what this algorithm does, please consider looking at the following repo, providing example data analyses and simulation studies
https://github.com/FelixRb96/maxstablePCA_examples. 

### References

Cran: https://cran.r-project.org/package=maxstablePCA

* Principal component analysis for max-stable distributions, Reinbott F., Janßen A. , arxiv preprint, https://arxiv.org/abs/2408.10650
* A semi-group approach to Principal Component Analysis, Schlather M., Reinbott F., arxiv preprint, https://arxiv.org/pdf/2112.04026.pdf, 2021

