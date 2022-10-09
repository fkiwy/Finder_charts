# Finder_charts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7097859.svg)](https://doi.org/10.5281/zenodo.7097859)

 ```Finder_charts``` is a Python module to create multi-bands finder charts from image data of following sky surveys:
- DSS (DSS1 B, DSS1 R, DSS2 B, DSS2 R, DSS2 IR),
- 2MASS (J, H, K),
- Spitzer (IRAC1, IRAC2, IRAC3, IRAC4, MIPS24),
- WISE (W1, W2, W3, W4),
- UKIDSS (Y, J, H, K),
- UHS (J),
- VHS (Y, J, H, K),
- Pan-STARRS (g, r, i, z, y),
- DECam (g, r, i, z, Y).

It also creates a WISE time series of epochs 2010, (2013), and 2014-2021.

All images are reprojected so that north is up and east is to the left.

The ```create_finder_charts``` function only takes two mandatory arguments, RA and Dec in decimal degrees, which can be specified either as a scalar, Python sequence (list, tuple, ...), or Numpy array.

This is a minimal usage example:
```
from Finder_charts import create_finder_charts
create_finder_charts(16.9701886, 0.6992208)
```

The resulting finder charts are saved to the directory given by parameter ```directory``` and automatically opened if ```open_file``` is set to True.

### Module dependencies:
The Python Standard Library, NumPy, Matplotlib, Astropy and Pillow (PIL Fork)

### Parameter description:
- ```ra``` : right ascension in decimal degrees (type: float)
- ```dec``` : declination in decimal degrees (type: float)
- ```img_size``` : image size in arcseconds (type: int, default: 100)
- ```overlays``` : whether to plot catalog overlays (type: bool, default: False)
- ```overlay_color``` : catalog overlay color (type: str, default: 'red')
- ```dss``` : whether to create DSS image series (type: bool, default: True)
- ```twomass``` : whether to create 2MASS image series (type: bool, default: True)
- ```spitzer``` : whether to create Spitzer image series (type: bool, default: True)
- ```wise``` : whether to create WISE image series (type: bool, default: True)
- ```ukidss``` : whether to create UKIDSS image series (type: bool, default: True)
- ```vhs``` : whether to create VHS image series (type: bool, default: True)
- ```ps1``` : whether to create Pan-STARRS image series (type: bool, default: True)
- ```decam``` : whether to create DECam image series (type: bool, default: True)
- ```neowise``` : whether to create WISE time series (type: bool, default: True)
- ```neowise_contrast``` : WISE time series contrast (type: int, default: 3)
- ```chrono_order``` : whether to plot image series in chronological order (type: bool, default: True)
- ```object_info``` : whether to plot object information like coordinates, etc. (type: bool, default: True)
- ```directory``` : directory where the finder charts should be saved (type: str, default: tempfile.gettempdir())
- ```cache``` : whether to cache the downloaded files (type: bool, default: True)
- ```show_progress``` : whether to show the file download progress (type: bool, default: True)
- ```timeout``` : timeout for remote requests in seconds (type: int, default: 300)
- ```open_pdf``` : deprecated, replaced by ```open_file```
- ```open_file``` : whether to open the saved finder charts automatically (type: bool, default: True)
- ```file_format``` : output file format: pdf, png, eps, etc. (type: str, default: 'pdf')

### Example output:
![Example output](example_output.png)

#### With catalog overlays:
![Example output](example_output_with_overlays.png)
