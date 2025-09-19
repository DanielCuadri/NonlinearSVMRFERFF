# Nonlinear SVM RFE using RFF
Nonlinear SVM-RFE with Random Fourier Features (RFF). Extends traditional SVM-RFE beyond linear kernels by approximating nonlinear mappings with RFF, enabling scalable and efficient feature selection in high-dimensional datasets with complex relationships.

## Content

- `example.Rmd`: rmarkdown file using all main funtions as an example
- `R/`: folder containing all functions:
    - **rfe** This function performs recursive feature elimination (RFE) using random Fourier feature (RFF) transformations and SVM training with optional kernel types.
    - **rff**: This function generates Random Fourier Features (RFF) to approximate shift-invariant kernels, including Gaussian, Laplace, Cauchy, and Matérn kernels. It returns both the random projection vectors and the transformed feature matrix.
    - **orf**: This function generates orthogonal random Fourier features (ORF) to approximate shift-invariant kernels. ORF improves the variance of standard Random Fourier Features by using orthogonal projection matrices. Supports parallel computation for large feature dimensions.
    - **irff**: This function generates two sets of random vectors, `W_pos` and `W_neg`, for constructing indefinite random Fourier features (RFF).
    - **scale**: This function standardizes the numeric columns of a data frame.
    - **sigest**: This function computes heuristic estimates for the scale parameter (sigma) of a kernel based on the distances between observations. Supports Gaussian and Laplace kernels.
    - **utils**: This function computes the decision function values for a linear SVM in the original feature space, given support vector coefficients and bias.
- `data/`: datasets used for the experiment + Alon's cancer dataset + Golub's Leukemia dataset read in the example file.
