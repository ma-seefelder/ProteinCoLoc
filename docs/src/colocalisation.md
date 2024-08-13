```@meta
CurrentModule = ProteinCoLoc
```

# Computing colocalization

## Main functions und structs

### CoLocResult

The struct `CoLocResult` is used to store the results of the Bayesian inference, the input images, the identifiers of the analysed colour channels and the number of patches used for the anaylsis.

```@docs
CoLocResult
```

The following methods for the struct `CoLocResult` are defined that allow retrieving certain fields from the structure

```@docs
get_images(result::CoLocResult)
get_posterior
get_advi_result
get_config
```

### Inferring colocalization with the hierarchical Bayesian model

The function `colocalization` is one of major APIs of `ProteinCoLoc.jl`. This function performs:

- the necessary data manipulation (patching of the image)
- computes the patch-wise correlation of pixel intensities using the defined `cor_method`
- exclusion of images from the anaylsis where no signal is above background
- fitting of the hierarchical Bayesian model
- extraction of prior and posterior samples as a DataFrame

```@docs
colocalization
```

### Compute the Bayes factor

The function `compute_BayesFactor`computes the Bayes Factor for the hypothesis $H_1: \Delta \rho > \Delta \rho_0$ compared to the null hypothesis $H_0: \Delta \rho \leq \Delta \rho_0$.

Therefore, first the difference of the correlation metric $\Delta \rho$ between the sample and control images of the prior and posterior samples are computed. Afterwards, two normal distribution kernels are fitted on the samples of $\Delta \rho$ to obtain the prior $P(\Delta \rho)$ and posterior distribution $P(\Delta \rho | data)$ using the function `Distributions.kde` from `Distributions.jl`. Subsequently, the integral of probability densitiy function of the prior

```math
P(\Delta \rho > \rho_0) = 1 - \int_{-\infty}^{\Delta \rho} P(\Delta \rho > \rho_0) \,d\Delta \rho
```

and posterior distribution

```math
P(\Delta \rho > \rho_0 | data) = 1 - \int_{-\infty}^{\Delta \rho} P(\Delta \rho > \rho_0 | data) \,d\Delta \rho
```

are estimated numerically. The Bayes factor (BF) is then computed by dividing the posterior odds $ \frac{1 - P(\Delta \rho > \rho_0 | data)}{P(\Delta \rho > \rho_0 | data)} $ by the prior odds $ \frac{1 - P(\Delta \rho > \rho_0)}{P(\Delta \rho > \rho_0)} $

```math

BF = \frac{\frac{1 - P(\Delta \rho > \rho_0 | data)}{P(\Delta \rho > \rho_0 | data)}}{\frac{1 - P(\Delta \rho > \rho_0)}{P(\Delta \rho > \rho_0)}} = 

\frac{(1 - P(\Delta \rho > \rho_0 | data)) * P(\Delta \rho > \rho_0)}{P(\Delta \rho > \rho_0 | data) * (1 - P(\Delta \rho > \rho_0))}

```

Generally speaking, a $BF[H_1:H_0] = 1$ means that $H_1$ is *a posteriori* (after seeing the data) as likely as $H_0$. A $BF[H_1:H_0] = x$ means that $H_1$ is *a posteriori* (after seeing the data) x-times more likely than $H_0$. A common interpreation scale can be found in our open-access [paper](https://www.nature.com/articles/s41598-024-63884-1) (Figure 6) published in Scientific Reports.

```@docs
compute_BayesFactor
```

## Helper functions

### Dividing an image in to patches

The function `patch(img, num_patches_x, num_patches_y)` splits the image `img` into a `num_patches_x` x `num_patches_y`patches. This function returns a 4D-array with the following configuration:

1. Axis (i): Patches along the x-dimension (width of the image) with size `num_patches_x`
2. Axis (j): Patches along the y-dimension (height of the image) with size `num_patches_y`
3. Axis (k): x-dimension of the patch `patch[i,j,k,:]` with size `div(size(img, 2), num_patches_x, RoundDown)`
4. Axis (h): y-dimension of the patch `patch[i,j,:,h]` with size `div(size(img, 2), num_patches_y, RoundDown)`  

Additionally, the syntactic wrapper `patch(img::Array{Float64,2}, num_patches)` is available that splits the images into `num_patches` x `num_patches` patches.  

```@docs
patch(img::Array{Float64, 2}, num_patches::Int64)
patch(img::Array{Float64, 2}, num_patches_x::Int64, num_patches_y::Int64)
```

### Merging image patches into a single image

The function `unpatch(patches, img_size)` combines all patches given by a four-dimensional array obtained by `patch(img, num_patches_x, num_patches_y)` into a single channel image.

```@docs
unpatch
```

### Compute per-patch correlation between patched images

The function `correlation` can be used to compute correlation of pixel intensities for each image patch. Therefore, the 2D images needs first to be patched using the [patch](@ref) function. Three different correlation algorithms are available:

- `:pearson`: Pearson correlation coefficient (PCC). For details, please refer to the documentation of the function Statisctics package in Julia's standard library. Default method.
- `:spearman`: Spearman's rank correlation. For details, please refer to the [documentation of StatsBase.jl](https://juliastats.org/StatsBase.jl/stable/)
- `:kendall` Kendall rank correlation. For details, please refer to the [documentation of StatsBase.jl](https://juliastats.org/StatsBase.jl/stable/)

!!! note "Minimum number of pixels"
    A correlation score $\rho$ is calculated only when more than 5 pixels in a patch are non-masked. This requirement helps ensure the reliability of the results, as the correlation score may be unreliable if too few pixels have values above the background in colour channels A and B.

```@docs
correlation
```

### Retrieve posterior samples

The function `convert_posterior_samples` gets the parameter names from the Bayesian model, selects the necessary parameters, permutes the samples, and converts them into a DataFrame. Each column of the returning DataFrame represents a parameter of the Bayesian hierarchical model, while each row represents a sample from the posterior distribution.

```@docs
convert_posterior_samples
```
