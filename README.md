# smBEVO
Computer vision approach to intuitive baseline drift estimation for single-molecule data series.

## [Check out the bioRxiv preprint for the associated article!](https://doi.org/10.1101/2021.11.12.468397)

[![View smBEVO on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/101904-smbevo)

---

### Script

Download `smBEVO.m` and put it in yout MATLAB path. Now you can call the smBEVO function in your scripts.

Note that for spline smoothing you will also need the [SPLINEFIT](https://www.mathworks.com/matlabcentral/fileexchange/71225-splinefit/) add-on by Jonas Lundgren (you can either install this via MATLAB's Add-On Explorer or download it from the Mathworks File Exchange).

See `example.m` for use.

### MATLAB APP

Download `smBEVO.mlappinstall` and install it via MATLAB (e.g. see the Install App button in the APPS tab). Now you can launch the app by clicking on the smBEVO app icon in MATLAB's APPS tab.

The app provides a graphical interface for loading and baselineing data series and is the recommended approach for exploratory analysis as it provides immediate visual feedback while adjusting parameters.

### Tutorial

1. Get the `data series`. In the app the `Load Data` button will prompt to select a \*.mat file containing a *data* variable (see parameters below for format). In addition to the format described below, the app allows *data* to be a cell array for multiple series, in which case a spin box will be displayed for traversing the series. The app also reads/writes data as a struct (or struct array for multiple data series) containing both the data and all smBEVO results (this is what you get with the `Save Data` button).
2. :exclamation: Estimate appropriate values for the image filter along each dimension `sigmaX` and `sigmaY`!!! These two parameters are the only required parameters and they are `critical`. Fortunately, they are `intuitively estimated` (see parameters below). In the app, the `Update Image` button will give you a preview of the effect of your choice of `sigmaX` and `sigmaY` on the image representation of the data.
3. Set any additional `optional parameters` (see parameters below).
4. `Run smBEVO` with specified parameters on the data series. In the app the `Run` button will apply smBEVO to the currently displayed data series. If multiple data series are loaded, the `Run All` button will run on all of the data series. If you are using the app, then periodically during execution smBEVO will check to see if the `Abort` button has been pressed. Although typically very fast (often nearly immediate for many single-molecule time series), speed is dependent on the size of the image representation of the data. Relatively tiny `sigmaX` and `sigmaY` can result in large images that slow computation.
5. `Visualize results`. In the app you will see both the data series and its image representaiton overlaid with the identified levels. You can optionally show the baseline-subtracted or raw data series as well as an idealized trajectory. !!! Note that the idealization is likely to be a highly filtered version of the actual data series, but may still be of use for various estimations.
6. `Adjust parameters` until you are happy with the results. It is highly recommended to start by changing `sigmaX` and `sigmaY` before you worry about anything else as these will likely have the most impact on the results.

### Parameters

* `data`: Series data in either single-colum [y] or two-column [x y] format. If not specified, x-values default to sample indices.
* `sigmaX`: Defines image filter in x. This should be approximately the longest period where we can be reasonably sure that the baseline won't change much within this period wherever we look in the data.
* `sigmaY`: Defines the image filter in y. This should be approximately 1/3 of the smallest separtation in y between levels that we wish to detect. Note that good baseline estimation does not always require detecting every level, especially if they are closely spaced.
* `pixelsPerSigmaX`: Number of pixels per sigmaX in the image representation. Default value of 4 seems to nearly always work well. Anything much higher than ~10 is likely to slow computation without helping much.
* `pixelsPerSigmaY`: Number of pixels per sigmaY in the image representation. Default value of 4 seems to nearly always work well. Anything much higher than ~10 is likely to slow computation without helping much.
* `smoothing`: Amount of spline smoothing (0: no smoothing, >0: increased smoothing).
* `minLevelSep`: Minimum allowable separation in y between ientified levels. If set to 0, a default value based on sigmaY is used.
* `snakeLevelRefinement`: Boolean indicating whether or not to refine the identified levels using active contour snakes.
* `maxSnakeIter`: Maximum number of snake refinement iterations. Default value of 50 nearly always works well.
* `alpha`: Degree of constraint on snake elasticity (1st derivative) (0: no constraint, >0: increased constraint ==> less wiggly levels). Default value of 1 seems to often works well.
* `beta`: Degree of constraint on snake curvature (2nd derivative) (0: no constraint, >0: increased constraint ==> less wiggly levels). Default value of 1 seems to often works well.
* `gamma`: Scale parameter for updating snakes on each iteration. Default value of 1 seems to nearly always work well.

### Hints

In many cases, appropriate choice of `sigmaX` and `sigmaY` are all that is needed for good results.

When necessary, additional use of `smoothing` and `minLevelSep` are usually sufficient to achieve desirable results.

For levels that more closely follow the image represenation of the data series, use `snakeLevelRefinement`.

