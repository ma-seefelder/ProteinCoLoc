```@meta
CurrentModule = ProteinCoLoc
```

# Main function

The function `main` is the main API of `ProteinCoLoc.jl` and allows to perform the whole analysis including image loading, and preprocessing, and generation of plots with one function call.

The following arguments need to be specified:

- `image_path`: Path to the folder comprising images of the individual colour channels. The file names of the images need to follow the naming conventions described in our [paper](https://doi.org/10.1038/s41598-024-63884-1) and in the documentation for the function `get_images()`. In short:

  - all images belonging to the same biological condition must be in the same folder
  - the individual channels should be named as $\left[ {\text{image name}} \right]\_\left[{\text{channel ID}} \right].{\text{tiff}}$.
  - If there are additional underscores "_" in the file names, the characters following the last underscore are used as channel identifier

- `control_image_path`: Path to the directory with the control images.
- `number_patches`: number of patches used to divide the image into.
- `number_patches_loc`: number of image patches for the local correlation analysis
- `number_channels`: number of different channels.

All other arguments are optional. With the optional argument `plot_options`, it can be defined which result plot should be generated. A dicitionary with the following keys needs to be passed to the `plot_options` argument:

- patched\_correlation\_plt
- local\_correlation\_plt
- bayes\_factor\_plt
- bayes\_range\_plt
- posterior\_plt
- mask\_plt

Details on the individual plots can be found in the documention under the section *Plot* or in the methods section of our [paper](https://doi.org/10.1038/s41598-024-63884-1).

```@docs
start_analysis
```
