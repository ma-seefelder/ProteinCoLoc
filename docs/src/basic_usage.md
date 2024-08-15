# Basic use of ProteinCoLoc.jl

In this tutorial, we will demonstrate two approaches for analysing colocalization in microscopic images. The first method uses the [`start_analysis()`](@ref) function, which performs the entire analysis with a single function call. The second method highlights the use of the low-level API, offering greater customisation and flexibility. Unlike `start_analysis`, this approach does not require a specific file naming convention.

For both tutorials, we will use images of A549 cells in which the huntingtin-associated protein 40 ([HAP40](https://content.iospress.com/articles/journal-of-huntingtons-disease/jhd220543)) has been stained using two differnt antibodies. The aim of this analysis was the verfication of a new antibody against HAP40. Further details on the experimental procedure can be found in our recent [publication](https://doi.org/10.1038/s41598-024-63884-1).
As a negative control, we used images of A549 cells that do not express recombinant HAP40. Consequently, one of the antibodies (anti-Strep) is unable to detect any HAP40 within these cells. Therefore, no colocalization should be observed in the images of the negative control.
For this tutorial, unlike the analysis presented in the paper, we will use only one image per condition. All necessary image files can be found in the test directory of [ProteinCoLoc's GitHub repository](https://github.com/ma-seefelder/ProteinCoLoc).

## Standard analysis pipeline (`start_analysis`)

To run the standard anaylsis pipline, we only have to call the function `start_analysis()`(@ref). In a first step, we define the paths to the images and define the a directory to save the results.

!!! note "No log.txt or results.txt file allowed in result directory"
    The execution of the function `start_analysis` generates a `results.txt` file. Executing the function over the graphical user interface (GUI) additionally generates a `log.txt` file. To avoid overwriting of previous results, using the same result directory for a subsequent analysis will result in an error.

```julia
# Path to the images with expressed recombinant HAP40
IMAGE_PATH = "./test/test_images/positive"
# Path to the images without expression of recombinant HAP40
CONTROL_IMAGE_PATH = "./test/test_images/positive"
# Path to the folder where the results should be stored.
OUTPUT_FOLDER_PATH = "./test/test_images/results"
```

As part of the anaylsis pipline, the images are divided into `NUM_PATCHES x NUM_PATCHES` patches and a correlation score for each of these patches is computed that is then used by the hierarchical model. Additionally, we need to specify a second number of patches for the local correlation anaylsis. These smaller patches will be exclusively used to generate the local correlation plot (e.g. [Fig 1B](https://www.nature.com/articles/s41598-024-63884-1#Fig1) and [Fig 4B](https://www.nature.com/articles/s41598-024-63884-1#Fig4) of our research article). By our experience, good starting values for number of patches `NUM_PATCHES` is 16 and 200 for as the number of patches for the local correlation analysis. If there are not sufficient pixels above background, the next lowest possible number will be used.

Additionally, we have to define the number of different colour channels that have been captured.

```julia
NUM_PATCHES = 16
NUM_PATCHES_LOC = 200
NUM_CHANNELS = 3
```

Now, we can call the function `start_analysis`. Note that the execution of this function can take up to five minutes, depending on the available hardware.

```julia
start_analysis(
    IMAGE_PATH, CONTROL_IMAGE_PATH, OUTPUT_FOLDER_PATH, 
    NUM_PATCHES, NUM_PATCHES_LOC, NUMBER_CHANNELS
)
```
!!! info "Note"
    During the execution of the start_analysis function, the info message "[ADVI] Should only be seen once: optimizer created for $\theta$" will be displayed. This message appears during the fitting of the hierarchical model using Turing.jl. It is normal for this info message to be displayed multiple times, as the optimizer needs to be created for each image and each channel combination.

By default, the colococalisation will be computed for each channel combinations. This can be changed by setting the keyword arguments `channel_selection = true` and specifying the channels by using the `channel_selection_two` argument as a vector of the channel indices `[1,2]`.

```julia
start_analysis(
    IMAGE_PATH, CONTROL_IMAGE_PATH, OUTPUT_FOLDER_PATH, 
    NUM_PATCHES, NUM_PATCHES_LOC, NUM_CHANNELS,
    channel_selection = true, channel_selection_two = [2,3]
)
```

Calling this function generates the following files that are stored in the `OUTPUT_FOLDER_PATH`:

- **Patched correlation plots**: A picture where the individual color channels have been merged after normalising. The image further displays the patch division and the computed correlation metric for each patch.
- **Local correlation plot**: A local correlation plot where the correlationn metric is displayed as as colour gradient. This allows the identification of local correlation patterns.
- **Bayes factor plot**: A plot displaying the prior and posterior distribution of $\Delta \rho$, a straight line depicteing the threshold $\Delta \rho$. Additionally, it displays the Bayes Factor for the hypothesis $H_1:\Delta \rho > \Delta \rho_0$ compared to the null hypothesis $H_0:\Delta \rho \leq \Delta \rho_0$. By default a threshold of $\Delta \rho = 0.1$ is used. This value can be freely chosen by setting `ρ_threshold` to the desired argument.
- **Bayes factor range plot**: This plot displays the Bayes Factor for the hypothesis $H_1:\Delta \rho > \Delta \rho_0$ compared to the null hypothesis $H_0:\Delta \rho \leq \Delta \rho_0$ for different values of $\Delta \rho$. The horizontal dashed lined depicts a Bayes Factor of 1.
- **Plot of prior and posterior**: Density plots of prior and posterior distributions for the parameters
  - ρ: Mean correlation metric for all images.
  - ν: Mean Degrees of freedom
  - σ: Mean standard deviation of ρ between images
  - τ: Mean standard deviation of ρ within an image
  - Δρ: Differene correlation metric between sample and control image
- **mask**: displays the image mask for each channel. Each masked pixel is displayed as black and non masked pixels as white.
- **posterior_samples.csv**: csv file with all posterior samples for each parameter.
- **prior_samples.csv**: csv file with all prior samples for each parameter.
- **results.txt**: A csv file comprising summary statistics for prior and posterior samples with the following columns: ρ\_threshold, bf, mean\_prior, median\_prior, lower\_credible\_interval\_prior, upper\_credible\_interval\_prior, mean\_posterior, median\_posterior, lower\_credible\_interval\_posterior, upper\_credible\_interval\_posterior, mean\_control\_prior, median\_control\_prior, lower\_credible\_interval\_control\_prior, upper\_credible\_interval\_control\_prior, mean\_control\_posterior, median\_control\_posterior, lower\_credible\_interval\_control\_posterior, upper\_credible\_interval\_control\_posterior, mean\_Δρ, median\_Δρ, lower\_credible\_interval\_Δρ, upper\_credible\_interval\_Δρ

### Controlling output plots

If not all output plots should be generated, the desired plots can be defined by specifying the `plot_options` keyword argument. A dictionary of type `Dict{AbstractString, Bool}` with the keys "patched_correlation_plt", "local_correlation_plt", "bayes_factor_plt", "bayes_range_plt", "posterior_plt", "mask_plt" needs to be passed to `plot_options` If the respecitive plot should be generated the value needs to be set to `true`, otherwise the value should be `false`.

## Custom analysis workflow

The manual workflow can be used if the file names do not follow the above described naming convention or the user wants to reuse only parts of the analysis workflow.

First, we load the images using the function (`ProteinCoLoc.load_tiff()`)[@ref]

```julia
# sample images
files = [
    "test/test_images/positive/positive_c1.tif",
    "test/test_images/positive/positive_c2.tif",
    "test/test_images/positive/positive_c3.tif"
]

blue, green, red = ProteinCoLoc.load_tiff.(files)

# control images
files_control = [
    "test/test_images/negative/negative_c1.tif",
    "test/test_images/negative/negative_c2.tif",
    "test/test_images/negative/negative_c3.tif"
]

blue_contorl, green_control, red_control = ProteinCoLoc.load_tiff.(files)
```

Next, we have to generate the `MultiChannelImage` objects. Here it is important that channels are read in the correct order [green, blue, red]. Otherwise, the displayed colour in the patched correlation plot are not correct.

```julia
import Images

data = [blue, green, red]
img = MultiChannelImage(
    data,  
    ["DAPI", "anti-HAP4", "anti-Strep"], "recombinant A549" , 
    files, size(data[1]), Images.otsu_threshold.(data)
    )

data = [blue_control, green_control, red_control]

control = MultiChannelImage(
    data,  
    ["DAPI", "anti-HAP4", "anti-Strep"], "control A549" , 
    files, size(data[1]), Images.otsu_threshold.(data)
    )
```

Easier, the `MultiChannelImage`can be directly generated using

```julia
img = MultiChannelImage("recombinant A549", files, ["DAPI", "anti-HAP4", "anti-Strep"])
control = MultiChannelImage("A549", files_control, ["DAPI", "anti-HAP4", "anti-Strep"])
```