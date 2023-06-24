import os
import sys
import math
import warnings
import requests
import tempfile
import traceback
import subprocess
import collections
import numpy as np
from io import BytesIO
from datetime import datetime
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
from matplotlib.patches import BoxStyle, Rectangle
from astropy.io import ascii
from astropy.visualization import make_lupton_rgb
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from astropy.utils.data import download_file
from astroquery.ukidss import Ukidss
from astroquery.vsa import Vsa
from astroquery.vizier import Vizier
# from astroquery.mast import Catalogs
import astropy.units as u
from astropy.coordinates import SkyCoord
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from pyvo.dal import sia


def create_finder_charts(ra, dec, img_size=100, overlays=False, overlay_color='red', dss=True, twomass=True, spitzer=True, wise=True,
                         ukidss=True, uhs=True, vhs=True, vvv=True, viking=True, ps1=True, decam=True, neowise=True, neowise_contrast=3,
                         chrono_order=True, object_info=True, directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, open_pdf=None,
                         open_file=True, file_format='pdf', save_result_tables=False, result_tables_format='ipac', result_tables_extension='dat'):
    """
    Creates multi-bands finder charts from image data of following sky surveys:
    - DSS (DSS1 B, DSS1 R, DSS2 B, DSS2 R, DSS2 IR),
    - 2MASS (J, H, K),
    - Spitzer (IRAC1, IRAC2, IRAC3, IRAC4, MIPS24),
    - WISE (W1, W2, W3, W4),
    - UKIDSS (Y, J, H, K),
    - UHS (J, K),
    - VHS (Y, J, H, K),
    - VVV (Z, Y, J, H, K),
    - VIKING (Z, Y, J, H, K),
    - Pan-STARRS (g, r, i, z, y),
    - DECam (g, r, i, z, Y).

    This function also creates a WISE time series of epochs 2010, (2013), and 2014-2021.
    Required arguments are ``ra`` and ``dec`` in decimal degrees, which can be specified either as a scalar, Python sequence (list, tuple, ...), or Numpy array.
    The resulting finder charts are saved to the directory given by parameter ``directory`` and automatically opened if ``open_file`` is set to True.

    Parameters
    ----------
    ra : float
        Right ascension in decimal degrees.
    dec : float
        Declination in decimal degrees.
    img_size : int, optional
        Image size in arcseconds. The default is 100.
    overlays : bool, optional
        Whether to plot catalog overlays. The default is False.
    overlay_color : str, optional
        Catalog overlay color. The default is 'red'.
    dss : bool, optional
        Whether to create DSS image series. The default is True.
    twomass : bool, optional
        Whether to create 2MASS image series. The default is True.
    spitzer : bool, optional
        Whether to create Spitzer image series. The default is True.
    wise : bool, optional
        Whether to create WISE image series. The default is True.
    ukidss : bool, optional
        Whether to create UKIDSS image series. The default is True.
    uhs : bool, optional
        Whether to create UHS image series. The default is True.
    vhs : bool, optional
        Whether to create VHS image series. The default is True.
    vvv : bool, optional
        Whether to create VVV image series. The default is True.
    viking : bool, optional
        Whether to create VIKING image series. The default is True.
    ps1 : bool, optional
        Whether to create Pan-STARRS image series. The default is True.
    decam : bool, optional
        Whether to create DECam image series. The default is True.
    neowise : bool, optional
        Whether to create WISE time series. The default is True.
    neowise_contrast : int, optional
        WISE time series contrast. The default is 3.
    chrono_order : bool, optional
        Whether to plot image series in chronological order. The default is True.
    object_info : bool, optional
        Whether to plot object information like coordinates, etc. The default is True.
    directory : str, optional
        Directory where the finder charts should be saved. The default is tempfile.gettempdir().
    cache : bool, optional
        Whether to cache the downloaded files. The default is True.
    show_progress : bool, optional
        Whether to show the file download progress. The default is True.
    timeout : int, optional
        Timeout for remote requests in seconds. The default is 300.
    open_pdf : bool, optional
        Deprecated, replaced by ``open_file``.
    open_file : bool, optional
        Whether to open the saved finder charts automatically. The default is True.
    file_format : str, optional
        Output file format: pdf, png, eps, etc.. The default is 'pdf'.
    save_result_tables : bool, optional
        Whether to save the catalog search result tables to the specified directory. The default is False.
    result_tables_format : str, optional
        Catalog search result tables output format. The default is 'ipac'.
    result_tables_extension : str, optional
        Catalog search result tables file extension. The default is 'dat'.

    Raises
    ------
    Exception
        - If ``ra`` and ``dec`` are neither a scalar nor a sequence or numpy array.
        - If ``ra`` and ``dec`` are either a sequence or a numpy array but don't have the same length.
    """

    class ImageBucket:
        def __init__(self, data, x, y, band, year_obs, wcs, overlay_ra=None, overlay_dec=None, overlay_phot=None):
            self.data = data
            self.x = x
            self.y = y
            self.band = band
            self.year_obs = year_obs
            self.wcs = wcs
            self.overlay_ra = overlay_ra
            self.overlay_dec = overlay_dec
            self.overlay_phot = overlay_phot

    def finder_charts(ra, dec):

        def process_image_data(hdu):
            try:
                data = hdu.data
                nanValues = np.count_nonzero(np.isnan(data))
                totValues = np.count_nonzero(data)

                if nanValues < totValues * 0.5:
                    wcs, shape = find_optimal_celestial_wcs([hdu])
                    data, _ = reproject_interp(hdu, wcs, shape_out=shape)
                    position = SkyCoord(ra*u.deg, dec*u.deg)
                    cutout = Cutout2D(data, position, img_size*u.arcsec, wcs=wcs, mode='partial')
                    data = cutout.data
                    wcs = cutout.wcs
                    x, y = wcs.world_to_pixel(position)
                    return data, x, y, wcs
                else:
                    return None, 0, 0, None
            except Exception:
                print('A problem occurred while creating an image for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                print(traceback.format_exc())
                return None, 0, 0, None

        def create_color_image(r, g, b, neowise=False):
            try:
                if r is None or g is None or b is None:
                    return None

                if neowise:
                    vmin, vmax = get_min_max(g, lo=neowise_contrast, hi=100-neowise_contrast)
                    image = Image.fromarray(make_lupton_rgb(r, g, b, minimum=vmin, stretch=vmax-vmin, Q=0))
                    image = ImageOps.invert(image)
                else:
                    xmax = max([r.shape[0], g.shape[0], b.shape[0]])
                    ymax = max([r.shape[1], g.shape[1], b.shape[1]])
                    r = Image.fromarray(create_lupton_rgb(r)).convert('L').resize((xmax, ymax), Image.ANTIALIAS)
                    g = Image.fromarray(create_lupton_rgb(g)).convert('L').resize((xmax, ymax), Image.ANTIALIAS)
                    b = Image.fromarray(create_lupton_rgb(b)).convert('L').resize((xmax, ymax), Image.ANTIALIAS)
                    image = Image.merge('RGB', (r, g, b))

                return np.array(image)
            except Exception:
                print('A problem occurred while creating a color image for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                print(traceback.format_exc())
                return None

        def plot_image(image_bucket, img_idx):
            try:
                if image_bucket.data is None:
                    return

                data = image_bucket.data
                x = image_bucket.x
                y = image_bucket.y
                band = image_bucket.band
                year_obs = int(round(image_bucket.year_obs, 0))
                wcs = image_bucket.wcs
                overlay_ra = image_bucket.overlay_ra
                overlay_dec = image_bucket.overlay_dec
                overlay_phot = image_bucket.overlay_phot

                ax = fig.add_subplot(rows, cols, img_idx, projection=wcs)
                ax.plot(x, y, 'ro', fillstyle='none', markersize=9, markeredgewidth=0.3)
                ax.plot(x, y, 'ro', fillstyle='none', markersize=0.3, markeredgewidth=0.3)
                ax.text(0.03, 0.93, band, color='black', fontsize=3.0, transform=ax.transAxes,
                        bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
                ax.text(0.03, 0.04, year_obs, color='black', fontsize=3.0, transform=ax.transAxes,
                        bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
                ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec='black', transform=ax.transAxes))

                if overlays and overlay_ra is not None and overlay_dec is not None and overlay_phot is not None:
                    ax.scatter(overlay_ra, overlay_dec, transform=ax.get_transform('icrs'), s=0/overlay_phot + 2.0,
                               edgecolor=overlay_color, facecolor='none', linewidths=0.3)

                vmin, vmax = get_min_max(data)
                ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray_r')
                ax.axis('off')
            except Exception:
                print('A problem occurred while plotting an image for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
                print(traceback.format_exc())

        def create_lupton_rgb(data):
            vmin, vmax = get_min_max(data)
            stretch = 1 if vmax-vmin == 0 else vmax-vmin
            return make_lupton_rgb(data, data, data, minimum=vmin, stretch=stretch, Q=0)

        def get_min_max(data, lo=10, hi=90):
            med = np.nanmedian(data)
            mad = np.nanmedian(abs(data - med))
            dev = np.nanpercentile(data, hi) - np.nanpercentile(data, lo)
            vmin = med - 2.0 * mad
            vmax = med + 2.0 * dev
            return vmin, vmax

        def get_neowise_image(ra, dec, epoch, band, size):
            download_url = 'http://byw.tools/cutout?ra={ra}&dec={dec}&size={size}&band={band}&epoch={epoch}'
            download_url = download_url.format(ra=ra, dec=dec, size=round(size/2.75), band=band, epoch=epoch)
            try:
                return fits.open(download_file(download_url, cache=cache, show_progress=show_progress, timeout=timeout))
            except Exception:
                return None

        def get_IRSA_image(ra, dec, survey, band, size):
            download_url = 'https://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?mode=getImage&file_type=fits&reproject=false&' + \
                'RA={ra}&DEC={dec}&subsetsize={size}&survey={survey}&{band}'
            download_url = download_url.format(ra=ra, dec=dec, size=size/60, survey=survey, band=band)
            try:
                return fits.open(download_file(download_url, cache=cache, show_progress=show_progress, timeout=timeout))
            except Exception:
                return None

        def get_UKIDSS_image(ra, dec, band, size, database):
            try:
                return Ukidss.get_images(SkyCoord(ra, dec, unit=(u.deg, u.deg)), image_width=size * u.arcsec, database=database,
                                         waveband=band, frame_type='stack', verbose=False, show_progress=show_progress)
            except Exception:
                print('A problem occurred while downloading UKIDSS images for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
                print(traceback.format_exc())
                return None

        def get_VSA_image(ra, dec, band, size, database, frame_type):
            try:
                return Vsa.get_images(SkyCoord(ra, dec, unit=(u.deg, u.deg)), image_width=size * u.arcsec, database=database,
                                      waveband=band, frame_type=frame_type, verbose=False, show_progress=show_progress)
            except Exception:
                print('A problem occurred while downloading VHS images for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
                print(traceback.format_exc())
                return None

        def get_DECam_image(ra, dec, band, size):
            try:
                fov = (size * u.arcsec).to(u.deg)
                svc = sia.SIAService('https://datalab.noirlab.edu/sia/coadd_all')
                result = svc.search((ra, dec), (fov/np.cos(dec*np.pi/180), fov), verbosity=2).to_table()
                select = np.char.startswith(result['obs_bandpass'].astype(str), band) & (result['proctype'] == 'Stack') & \
                    (result['prodtype'] == 'image') & (result['project_code'] == 'DES DR1')
                table = result[select]
                if len(table) > 0:
                    row = table[np.argmax(table['exptime'].data.data.astype('float'))]  # get image with longest exposure time
                    image = fits.open(download_file(row['access_url'], cache=cache, show_progress=show_progress, timeout=timeout))
                    instrument = row['instrument_name']
                    return image, instrument
                else:
                    return None, None
            except Exception:
                print('A problem occurred while downloading DECam images for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
                print(traceback.format_exc())
                return None, None

        def search_noirlab(ra, dec, radius):
            try:
                radius = radius.to(u.deg)
                query_url = 'https://datalab.noirlab.edu/tap/sync'
                adql = """
                    SELECT ra, dec, gmag, rmag, imag, zmag, ymag
                    FROM   nsc_dr2.object
                    WHERE  't'=q3c_radial_query(ra, dec, {ra}, {dec}, {radius})
                    """
                adql = adql.format(ra=ra, dec=dec, radius=radius.value)
                payload = {
                    'request': 'doQuery',
                    'lang': 'ADQL',
                    'format': 'csv',
                    'query': adql
                }
                response = requests.get(query_url, params=payload, timeout=timeout)
                table = ascii.read(response.text, format='csv')
                if len(table) > 0:
                    return table
                else:
                    return None
            except Exception:
                print('A problem occurred while downloading Noirlab catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                print(traceback.format_exc())
                return None

        def search_panstarrs(ra, dec, radius):
            url = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.votable'
            params = {
                'ra': ra,
                'dec': dec,
                'radius': radius.to(u.deg).value,
                'nStackDetections.gte': 2,
                'columns': ['raMean', 'decMean', 'gMeanPSFMag', 'rMeanPSFMag', 'iMeanPSFMag', 'zMeanPSFMag', 'yMeanPSFMag', 'distance'],
                'sort_by': ['distance']
            }
            response = requests.get(url, params=params, timeout=timeout)
            table = Table.read(BytesIO(response.content), format='votable')
            if len(table) > 0:
                table = table.filled(np.nan)
                return table
            else:
                return None

        def get_year_obs(hdu, date_obs_key, date_pattern):
            header = hdu.header
            # print(header)
            if date_pattern == 'MJD':
                year = get_year_from_mjd(header[date_obs_key])
            else:
                date_obs = header[date_obs_key]
                date_len = len(date_obs)
                size = date_len if date_len < 10 else 10
                time = datetime.strptime(date_obs[0:size], date_pattern)
                year = time.year
            return year

        def get_year_from_mjd(mjd):
            time = Time(mjd, scale='utc', format='mjd')
            return time.ymdhms['year']

        def create_obj_name(ra, dec, precision=6):
            ra = round(ra, precision)
            dec = round(dec, precision)
            ra_str = str(ra)
            dec_str = str(dec) if dec < 0 else '+' + str(dec)
            return ra_str + dec_str

        def save_catalog_search_results(result_table, survey, ra, dec):
            if save_result_tables:
                file_name = 'Finder_charts_' + survey + '_results_' + create_obj_name(ra, dec) + '.' + result_tables_extension
                result_table.write(file_name, format=result_tables_format, overwrite=True)

        def start_file(filename):
            if sys.platform == 'win32':
                os.startfile(filename)
            else:
                opener = 'open' if sys.platform == 'darwin' else 'evince'
                subprocess.call([opener, filename])

        # -------------------------------
        # Code for finder_charts function
        # -------------------------------
        fig = plt.figure()
        fig.set_figheight(15)
        fig.set_figwidth(15)

        coords = SkyCoord(ra*u.deg, dec*u.deg)
        radius = (img_size*math.sqrt(2)/2)*u.arcsec

        surveys = []

        # DSS
        if dss:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss1_blue', img_size)
            if image:
                data, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(data, x, y, 'DSS1 B', get_year_obs(image[0], date_obs_key, date_pattern), wcs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss1_red', img_size)
            if image:
                data, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(data, x, y, 'DSS1 R', get_year_obs(image[0], date_obs_key, date_pattern), wcs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_blue', img_size)
            if image:
                year_b = get_year_obs(image[0], date_obs_key, date_pattern)
                b, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(b, x, y, 'DSS2 B', year_b, wcs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_red', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(g, x, y, 'DSS2 R', year_g, wcs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_ir', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(r, x, y, 'DSS2 IR', year_r, wcs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'DSS2 IR-R-B', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # 2MASS
        TWO_MASS_K = '2MASS K'
        if twomass:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'ORDATE'
            date_pattern = '%y%m%d'
            survey = []

            if overlays or save_result_tables:
                try:
                    v = Vizier(columns=['+_r', 'all'])
                    tables = v.query_region(coords, radius=radius, catalog='II/246/out')
                    if tables:
                        table = tables[0]
                        overlay_ra = table['RAJ2000']
                        overlay_dec = table['DEJ2000']
                        op1 = table['Jmag']
                        op2 = table['Hmag']
                        op3 = table['Kmag']
                        save_catalog_search_results(table, '2MASS', ra, dec)
                except Exception:
                    print('A problem occurred while downloading 2MASS catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            image = get_IRSA_image(ra, dec, '2mass', 'twomass_bands=j', img_size)
            if image:
                year_b = get_year_obs(image[0], date_obs_key, date_pattern)
                b, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(b, x, y, '2MASS J', year_b, wcs, overlay_ra, overlay_dec, op1))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, '2mass', 'twomass_bands=h', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(g, x, y, '2MASS H', year_g, wcs, overlay_ra, overlay_dec, op2))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, '2mass', 'twomass_bands=k', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(r, x, y, TWO_MASS_K, year_r, wcs, overlay_ra, overlay_dec, op3))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, '2MASS K-H-J', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)
            survey.append(None)

            surveys.append(survey)

        # Spitzer
        if spitzer:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            date_obs_key = 'MJDMEAN'
            date_pattern = 'MJD'
            survey = []

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC1', img_size)
            if image:
                year_b = get_year_obs(image[0], date_obs_key, date_pattern)
                b, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(b, x, y, 'IRAC1', year_b, wcs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC2', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(g, x, y, 'IRAC2', year_g, wcs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC3', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(r, x, y, 'IRAC3', year_r, wcs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC4', img_size)
            if image:
                data, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(data, x, y, 'IRAC4', get_year_obs(image[0], date_obs_key, date_pattern), wcs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:MIPS24', img_size)
            if image:
                data, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(data, x, y, 'MIPS24', get_year_obs(image[0], date_obs_key, date_pattern), wcs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'IRAC3-2-1', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # WISE
        if wise:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'MIDOBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            if overlays or save_result_tables:
                try:
                    v = Vizier(columns=['+_r', 'all'])
                    tables = v.query_region(coords, radius=radius, catalog='II/328/allwise')
                    if tables:
                        table = tables[0]
                        overlay_ra = table['RAJ2000']
                        overlay_dec = table['DEJ2000']
                        op1 = table['W1mag']
                        op2 = table['W2mag']
                        op3 = table['W3mag']
                        op4 = table['W4mag']
                        save_catalog_search_results(table, 'AllWISE', ra, dec)
                except Exception:
                    print('A problem occurred while downloading AllWISE catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=1', img_size)
            if image:
                year_b = get_year_obs(image[0], date_obs_key, date_pattern)
                b, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(b, x, y, 'W1', year_b, wcs, overlay_ra, overlay_dec, op1))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=2', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(g, x, y, 'W2', year_g, wcs, overlay_ra, overlay_dec, op2))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=3', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(r, x, y, 'W3', year_r, wcs, overlay_ra, overlay_dec, op3))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=4', img_size)
            if image:
                data, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(data, x, y, 'W4', get_year_obs(image[0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op4))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'W3-W2-W1', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)

            surveys.append(survey)

        # UKIDSS
        if ukidss:
            x = y = x_j = y_j = 0
            r = g = b = wcs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'UKIDSSDR11PLUS'
            if overlays or save_result_tables:
                try:
                    table = Ukidss.query_region(coords, radius, database=database, programme_id='LAS')
                    if table:
                        table.sort('distance')
                        overlay_ra = table['ra']
                        overlay_dec = table['dec']
                        op1 = table['yAperMag3']
                        op2 = table['jAperMag3']
                        op3 = table['hAperMag3']
                        op4 = table['kAperMag3']
                        save_catalog_search_results(table, 'UKIDSS_DR11', ra, dec)
                except Exception:
                    print('A problem occurred while downloading UKIDSS catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            images = get_UKIDSS_image(ra, dec, 'Y', img_size, database)
            if images:
                data, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(data, x, y, 'UKIDSS Y', get_year_obs(images[0][0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op1))
            else:
                survey.append(None)

            images = get_UKIDSS_image(ra, dec, 'J', img_size, database)
            if images:
                year_b = get_year_obs(images[0][0], date_obs_key, date_pattern)
                b, x_j, y_j, wcs_j = process_image_data(images[0][1])
                survey.append(ImageBucket(b, x_j, y_j, 'UKIDSS J', year_b, wcs_j, overlay_ra, overlay_dec, op2))
            else:
                survey.append(None)

            images = get_UKIDSS_image(ra, dec, 'H', img_size, database)
            if images:
                year_g = get_year_obs(images[0][0], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(g, x, y, 'UKIDSS H', year_g, wcs, overlay_ra, overlay_dec, op3))
            else:
                survey.append(None)

            images = get_UKIDSS_image(ra, dec, 'K', img_size, database)
            if images:
                year_r = get_year_obs(images[0][0], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(r, x, y, 'UKIDSS K', year_r, wcs, overlay_ra, overlay_dec, op4))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x_j, y_j, 'UKIDSS K-H-J', mean_obs_year, wcs_j))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)

            surveys.append(survey)

        # UHS
        if uhs:
            x = y = x_j = y_j = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'UHSDR2'
            """ This is not working yet. Program id UHS is unknown.
            if overlays or save_result_tables:
                try:
                    table = Ukidss.query_region(coords, radius, database=database, programme_id='LAS')  # programme_id='UHS'
                    if table:
                        table.sort('distance')
                        table.pprint_all()
                        overlay_ra = table['ra']
                        overlay_dec = table['dec']
                        op1 = table['jAperMag3']
                        op2 = table['kAperMag3']
                        save_catalog_search_results(table, 'UHS_DR2', ra, dec)
                except Exception:
                    print(table)
                    print('A problem occurred while downloading UHS catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())
            """
            images = get_UKIDSS_image(ra, dec, 'J', img_size, database)
            if images:
                year_b = get_year_obs(images[0][0], date_obs_key, date_pattern)
                b, x_j, y_j, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(b, x_j, y_j, 'UHS J', year_b, wcs, overlay_ra, overlay_dec, op2))
            else:
                survey.append(None)

            images = get_UKIDSS_image(ra, dec, 'K', img_size, database)
            if images:
                year_r = get_year_obs(images[0][0], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(r, x, y, 'UHS K', year_r, wcs, overlay_ra, overlay_dec, op4))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, (b+r)/2, b), x_j, y_j, 'UHS K-J', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)
            survey.append(None)
            survey.append(None)

            surveys.append(survey)

        # VHS
        if vhs:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'VHSDR6'
            if overlays or save_result_tables:
                try:
                    table = Vsa.query_region(coords, radius, database=database, programme_id='VHS')
                    if table:
                        table.sort('distance')
                        overlay_ra = table['ra']
                        overlay_dec = table['dec']
                        op1 = table['yAperMag3']
                        op2 = table['jAperMag3']
                        op3 = table['hAperMag3']
                        op4 = table['ksAperMag3']
                        save_catalog_search_results(table, 'VHS_DR6', ra, dec)
                except Exception:
                    print('A problem occurred while downloading VHS catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            frame_type = 'tilestack'
            images = get_VSA_image(ra, dec, 'Y', img_size, database, frame_type)
            if images:
                data, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(data, x, y, 'VHS Y', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op1))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'J', img_size, database, frame_type)
            if images:
                year_b = get_year_obs(images[0][1], date_obs_key, date_pattern)
                b, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(b, x, y, 'VHS J', year_b, wcs, overlay_ra, overlay_dec, op2))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'H', img_size, database, frame_type)
            if images:
                year_g = get_year_obs(images[0][1], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(g, x, y, 'VHS H', year_g, wcs, overlay_ra, overlay_dec, op3))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Ks', img_size, database, frame_type)
            if images:
                year_r = get_year_obs(images[0][1], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(r, x, y, 'VHS K', year_r, wcs, overlay_ra, overlay_dec, op4))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'VHS K-H-J', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)

            surveys.append(survey)

        # VVV
        if vvv:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'VVVDR5'
            if overlays or save_result_tables:
                try:
                    table = Vsa.query_region(coords, radius, database=database, programme_id='VVV')
                    if table:
                        table.sort('distance')
                        overlay_ra = table['ra']
                        overlay_dec = table['dec']
                        op1 = table['zAperMag3']
                        op2 = table['yAperMag3']
                        op3 = table['jAperMag3']
                        op4 = table['hAperMag3']
                        op5 = table['ksAperMag3']
                        save_catalog_search_results(table, 'VVV_DR5', ra, dec)
                except Exception:
                    print('A problem occurred while downloading VVV catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            frame_type = 'deep_stack'
            images = get_VSA_image(ra, dec, 'Z', img_size, database, frame_type)
            if images:
                data, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(data, x, y, 'VVV Z', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op1))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Y', img_size, database, frame_type)
            if images:
                data, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(data, x, y, 'VVV Y', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op2))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'J', img_size, database, frame_type)
            if images:
                year_b = get_year_obs(images[0][1], date_obs_key, date_pattern)
                b, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(b, x, y, 'VVV J', year_b, wcs, overlay_ra, overlay_dec, op3))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'H', img_size, database, frame_type)
            if images:
                year_g = get_year_obs(images[0][1], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(g, x, y, 'VVV H', year_g, wcs, overlay_ra, overlay_dec, op4))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Ks', img_size, database, frame_type)
            if images:
                year_r = get_year_obs(images[0][1], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(r, x, y, 'VVV K', year_r, wcs, overlay_ra, overlay_dec, op5))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'VVV K-H-J', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # VIKING
        if viking:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'VIKINGDR5'
            if overlays or save_result_tables:
                try:
                    table = Vsa.query_region(coords, radius, database=database, programme_id='VIKING')
                    if table:
                        table.sort('distance')
                        overlay_ra = table['ra']
                        overlay_dec = table['dec']
                        op1 = table['zAperMag3']
                        op2 = table['yAperMag3']
                        op3 = table['jAperMag3']
                        op4 = table['hAperMag3']
                        op5 = table['ksAperMag3']
                        save_catalog_search_results(table, 'VIKING_DR5', ra, dec)
                except Exception:
                    print('A problem occurred while downloading VIKING catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            frame_type = 'tilestack'
            images = get_VSA_image(ra, dec, 'Z', img_size, database, frame_type)
            if images:
                data, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(data, x, y, 'VIKING Z', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op1))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Y', img_size, database, frame_type)
            if images:
                data, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(data, x, y, 'VIKING Y', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op2))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'J', img_size, database, frame_type)
            if images:
                year_b = get_year_obs(images[0][1], date_obs_key, date_pattern)
                b, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(b, x, y, 'VIKING J', year_b, wcs, overlay_ra, overlay_dec, op3))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'H', img_size, database, frame_type)
            if images:
                year_g = get_year_obs(images[0][1], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(g, x, y, 'VIKING H', year_g, wcs, overlay_ra, overlay_dec, op4))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Ks', img_size, database, frame_type)
            if images:
                year_r = get_year_obs(images[0][1], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(images[0][1])
                survey.append(ImageBucket(r, x, y, 'VIKING K', year_r, wcs, overlay_ra, overlay_dec, op5))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'VIKING K-H-J', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # PS1
        if ps1:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'MJD-OBS'
            date_pattern = 'MJD'
            survey = []

            try:
                query_url = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py'
                payload = {
                    'ra': ra,
                    'dec': dec,
                    'filters': 'grizy',
                    'sep': 'comma'
                }
                text = requests.get(query_url, params=payload, timeout=timeout).text
            except Exception:
                text = None
                print('A problem occurred while downloading Pan-STARRS image urls for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                print(traceback.format_exc())

            if text and text.count('\n') > 0:
                table = ascii.read(text)

                images = {}

                for row in table:
                    try:
                        download_url = 'http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?format=fits&red={filename}&ra={ra}&dec={dec}&size={size}'
                        download_url = download_url.format(filename=row['filename'], ra=ra, dec=dec, size=img_size*4)
                        images[row['filter']] = fits.open(download_file(download_url, cache=cache, show_progress=show_progress, timeout=timeout))
                    except Exception:
                        print('A problem occurred while downloading Pan-STARRS images for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                        print(traceback.format_exc())

                if images:
                    if overlays or save_result_tables:
                        try:
                            # table = Catalogs.query_region(coords, radius=radius, catalog='Panstarrs', data_release='dr2', table='mean',
                            #                               nStackDetections=[('gte', 2)], sort_by=[('asc', 'distance')])
                            table = search_panstarrs(ra, dec, radius)
                            if table:
                                overlay_ra = table['raMean'][:, 0]
                                overlay_dec = table['decMean'][:, 0]
                                op1 = table['gMeanPSFMag'][:, 0]
                                op2 = table['rMeanPSFMag'][:, 0]
                                op3 = table['iMeanPSFMag'][:, 0]
                                op4 = table['zMeanPSFMag'][:, 0]
                                op5 = table['yMeanPSFMag'][:, 0]
                                save_catalog_search_results(table, 'PS1_DR2', ra, dec)
                        except Exception:
                            print('A problem occurred while downloading Pan-STARRS catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                            print(traceback.format_exc())

                    if 'g' in images:
                        year_b = get_year_obs(images['g'][0], date_obs_key, date_pattern)
                        b, x, y, wcs = process_image_data(images['g'][0])
                        survey.append(ImageBucket(b, x, y, 'PS1 g', year_b, wcs, overlay_ra, overlay_dec, op1))
                    else:
                        survey.append(None)

                    if 'r' in images:
                        data, x, y, wcs = process_image_data(images['r'][0])
                        survey.append(ImageBucket(data, x, y, 'PS1 r', get_year_obs(images['r'][0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op2))
                    else:
                        survey.append(None)

                    if 'i' in images:
                        year_g = get_year_obs(images['i'][0], date_obs_key, date_pattern)
                        g, x, y, wcs = process_image_data(images['i'][0])
                        survey.append(ImageBucket(g, x, y, 'PS1 i', year_g, wcs, overlay_ra, overlay_dec, op3))
                    else:
                        survey.append(None)

                    if 'z' in images:
                        data, x, y, wcs = process_image_data(images['z'][0])
                        survey.append(ImageBucket(data, x, y, 'PS1 z', get_year_obs(images['z'][0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op4))
                    else:
                        survey.append(None)

                    if 'y' in images:
                        year_r = get_year_obs(images['y'][0], date_obs_key, date_pattern)
                        r, x, y, wcs = process_image_data(images['y'][0])
                        survey.append(ImageBucket(r, x, y, 'PS1 y', year_r, wcs, overlay_ra, overlay_dec, op5))
                    else:
                        survey.append(None)

                    if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                        mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                        survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'PS1 y-i-g', mean_obs_year, wcs))
                    else:
                        survey.insert(0, None)

                    survey.insert(0, mean_obs_year)

                    surveys.append(survey)

        # DECam
        if decam:
            x = y = 0
            r = g = b = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'MJD-OBS'
            date_pattern = 'MJD'
            survey = []

            if overlays or save_result_tables:
                """
                table = search_noirlab(ra, dec, radius)
                if table:
                    overlay_ra = table['ra']
                    overlay_dec = table['dec']
                    op1 = table['gmag']
                    op2 = table['rmag']
                    op3 = table['imag']
                    op4 = table['zmag']
                    op5 = table['ymag']
                """
                try:
                    v = Vizier(columns=['+_r', 'all'])
                    tables = v.query_region(coords, radius=radius, catalog='II/357/des_dr1')
                    if tables:
                        table = tables[0]
                        overlay_ra = table['RAJ2000']
                        overlay_dec = table['DEJ2000']
                        op1 = table['gmag']
                        op2 = table['rmag']
                        op3 = table['imag']
                        op4 = table['zmag']
                        op5 = table['Ymag']
                        save_catalog_search_results(table, 'DES_DR1', ra, dec)
                except Exception:
                    print('A problem occurred while downloading DES catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            image, instrument = get_DECam_image(ra, dec, 'g', img_size)
            if image:
                year_b = get_year_obs(image[0], date_obs_key, date_pattern)
                b, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(b, x, y, instrument + ' g', year_b, wcs, overlay_ra, overlay_dec, op1))
            else:
                survey.append(None)

            image, instrument = get_DECam_image(ra, dec, 'r', img_size)
            if image:
                data, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(data, x, y, instrument + ' r', get_year_obs(image[0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op2))
            else:
                survey.append(None)

            image, instrument = get_DECam_image(ra, dec, 'i', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(g, x, y, instrument + ' i', year_g, wcs, overlay_ra, overlay_dec, op3))
            else:
                survey.append(None)

            image, instrument = get_DECam_image(ra, dec, 'z', img_size)
            if image:
                data, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(data, x, y, instrument + ' z', get_year_obs(image[0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op4))
            else:
                survey.append(None)

            image, instrument = get_DECam_image(ra, dec, 'Y', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs = process_image_data(image[0])
                survey.append(ImageBucket(r, x, y, instrument + ' Y', year_r, wcs, overlay_ra, overlay_dec, op5))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, instrument + ' Y-i-g', mean_obs_year, wcs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # WISE time series
        if neowise:
            imageW1 = get_neowise_image(ra, dec, epoch=0, band=1, size=img_size)
            imageW2 = get_neowise_image(ra, dec, epoch=0, band=2, size=img_size)
            if imageW1 and imageW2:
                hdu = imageW1[0]
                header = hdu.header
                meanmjd = (header['MJDMIN']+header['MJDMAX'])/2
                prev_year = get_year_from_mjd(meanmjd)
                dataW1 = hdu.data
                dataW2 = imageW2[0].data
                j = 1
                images = []

                for i in range(1, 100, 1):
                    try:
                        imageW1 = get_neowise_image(ra, dec, epoch=i, band=1, size=img_size)
                        imageW2 = get_neowise_image(ra, dec, epoch=i, band=2, size=img_size)
                        if not imageW1 and not imageW2:
                            break
                        hdu = imageW1[0]
                        header = hdu.header
                        meanmjd = (header['MJDMIN']+header['MJDMAX'])/2
                        year = get_year_from_mjd(meanmjd)
                        if year == prev_year:
                            dataW1 += hdu.data
                            dataW2 += imageW2[0].data
                            j += 1
                        else:
                            hduW1 = fits.PrimaryHDU(data=dataW1/j, header=header)
                            hduW2 = fits.PrimaryHDU(data=dataW2/j, header=header)
                            images.append((hduW1, hduW2, prev_year))
                            dataW1 = hdu.data
                            dataW2 = imageW2[0].data
                            j = 1
                        prev_year = year
                    except Exception:
                        print('A problem occurred while creating WISE time series for object ra={ra}, dec={dec}, epoch={epoch}'
                              .format(ra=ra, dec=dec, epoch=i))
                        print(traceback.format_exc())

                hduW1 = fits.PrimaryHDU(data=dataW1/j, header=header)
                hduW2 = fits.PrimaryHDU(data=dataW2/j, header=header)
                images.append((hduW1, hduW2, prev_year))

                survey = []

                for image in images:
                    w1, x, y, wcs = process_image_data(image[0])
                    w2, x, y, wcs = process_image_data(image[1])
                    if w1 is not None and w2 is not None:
                        survey.append(ImageBucket(create_color_image(w1, (w1+w2)/2, w2, neowise=True), x, y, 'W1-W2', image[2], wcs))

                survey.insert(0, 9999)

                surveys.append(survey)

        # Plot image series
        if chrono_order:
            surveys.sort(key=lambda x: x[0])  # sort by mean observation year
        cols = 6
        rows = 15
        img_idx = 0
        info_idx = 0
        for survey in surveys:
            survey = survey[1:len(survey)]
            if all(image_bucket is None for image_bucket in survey):
                continue
            for image_bucket in survey:
                img_idx += 1
                if image_bucket is not None:
                    plot_image(image_bucket, img_idx)
                    if image_bucket.band == TWO_MASS_K:
                        info_idx = img_idx

        # Determine info text index
        if info_idx == 0:
            info_idx = math.ceil(img_idx / cols) * cols
        info_idx += 1

        if object_info:
            # Info text
            fontsize = 5.0
            ax = fig.add_subplot(rows, cols, info_idx)
            ax.text(0.05, 0.70, r'$\alpha$ = ' + str(round(coords.ra.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.55, r'$\delta$ = ' + str(round(coords.dec.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.40, '$l$ = ' + str(round(coords.galactic.l.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.25, '$b$ = ' + str(round(coords.galactic.b.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.axis('off')

            # Info text cont'd
            hmsdms = coords.to_string('hmsdms', sep=':', precision=2)
            hms = hmsdms[0:11]
            dms = hmsdms[12:24] if dec < 0 else hmsdms[13:24]
            ax = fig.add_subplot(rows, cols, info_idx + 1)
            ax.text(0, 0.72, '(' + hms + ')', fontsize=fontsize, transform=ax.transAxes)
            ax.text(0, 0.57, '(' + dms + ')', fontsize=fontsize, transform=ax.transAxes)
            ax.text(0, 0.42, 'Size = ' + str(int(img_size)) + ' arcsec', fontsize=fontsize, transform=ax.transAxes)
            ax.text(0, 0.27, 'North up, East left', fontsize=fontsize, transform=ax.transAxes)
            ax.axis('off')

        # Save and open the PDF file
        filename = 'Finder_charts_' + create_obj_name(ra, dec) + '.' + file_format
        plt.subplots_adjust(wspace=0, hspace=0.05, right=0.43)
        plt.savefig(filename, dpi=600, bbox_inches='tight', format=file_format)
        plt.close()

        if open_pdf if open_pdf is not None else open_file:
            start_file(filename)

    # --------------------------------------
    # Code for create_finder_charts function
    # --------------------------------------

    # Parameter deprecation warnings
    if open_pdf is not None:
        warnings.warn('Parameter ``open_pdf`` is deprecated. Please use ``open_file`` instead.', DeprecationWarning, stacklevel=2)

    warnings.simplefilter('ignore', category=AstropyWarning)
    os.chdir(directory)

    if np.isscalar(ra) and np.isscalar(dec):
        finder_charts(ra, dec)
    else:
        if not isinstance(ra, (collections.Sequence, np.ndarray)):
            raise Exception('Ra must be either a sequence or a numpy array')

        if not isinstance(dec, (collections.Sequence, np.ndarray)):
            raise Exception('Dec must be either a sequence or a numpy array')

        if len(ra) != len(dec):
            raise Exception('Ra and Dec must have the same length: len(ra)={ra}, len(dec)={dec}'.format(ra=len(ra), dec=len(dec)))

        for i in range(len(ra)):
            try:
                finder_charts(ra[i], dec[i])
            except Exception:
                print('A problem occurred while creating finder charts for object ra={ra}, dec={dec}'.format(ra=ra[i], dec=dec[i]))
                print(traceback.format_exc())
