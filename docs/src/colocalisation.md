```@meta
CurrentModule = ProteinCoLoc
```
# Computing colocalization

## Main functions

### Compute the Bayes factor 
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

# Compute per-patch correlation between patched images
The function `correlation` can be used to compute correlation of pixel intensities for each image patch. Therefore, the 2D images needs first to be patched using the [patch](@ref) function. Three different correlation algorithms are available:

- `:pearson`: Pearson correlation coefficient (PCC). For details, please refer to the documentation of the function Statisctics package in Julia's standard library. Default method.
- `:spearman`: Spearman's rank correlation. For details, please refer to the [documentation of StatsBase.jl](https://juliastats.org/StatsBase.jl/stable/)
- `:kendall` Kendall rank correlation. For details, please refer to the [documentation of StatsBase.jl](https://juliastats.org/StatsBase.jl/stable/)


!!! note "Minimum number of pixels"
    A correlation score $\rho$ is calculated only when more than 15 pixels in a patch are non-masked. This requirement helps ensure the reliability of the results, as the correlation score may be unreliable if too few pixels have values above the background in colour channels A and B.

```@docs
correlation
```



