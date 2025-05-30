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
from enum import Enum
from datetime import datetime
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
from matplotlib.patches import BoxStyle, Rectangle
from astropy.io import ascii
from astropy.visualization import make_lupton_rgb
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.time import Time
from astropy.utils.data import download_file
from astroquery.ukidss import Ukidss
from astroquery.vsa import Vsa
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
# from astroquery.mast import Catalogs
import astropy.units as u
from astropy.coordinates import SkyCoord
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from pyvo.dal import sia
from Gaia_finder_chart import create_gaia_finder_chart


class Crosshair(Enum):
    NONE = 0
    CIRCLE_DOT = 1
    MULTI_CIRCLES = 2
    CROSS_NO_CENTER = 3


class Survey(Enum):
    NONE = 'NONE'
    ALL = 'ALL'
    DSS = 'DSS'
    TWO_MASS = '2MASS'
    SPITZER = 'IRAC'
    ALLWISE = 'WISE'
    UKIDSS = 'UKIDSS'
    UHS = 'UHS'
    VHS = 'VHS'
    VVV = 'VVV'
    VIKING = 'VIKING'
    PS1 = 'PS1'
    DECAM = 'DECam'
    GAIA = 'Gaia'


class Target:
    def __init__(self, catalog, epoch, ra, dec, marker_size=3, marker_color='red', survey=Survey.NONE):
        self.catalog = catalog
        self.epoch = epoch
        self.ra = ra
        self.dec = dec
        self.marker_size = marker_size
        self.marker_color = marker_color
        self.survey = survey


def create_finder_charts(ra, dec, img_size=100, overlays=False, overlay_color='red', dss=True, twomass=True, spitzer=True, wise=True,
                         ukidss=True, uhs=True, vhs=True, vvv=True, viking=True, ps1=True, decam=True, neowise=True, neowise_contrast=3,
                         gaia_images=False, gaia_entries=False, gaia_pm_vectors=False, gaia_pm_years=10, gaia_color='green', targets=None,
                         pmra=None, pmdec=None, ref_epoch=None, crosshair_type=Crosshair.CIRCLE_DOT, propagate_gaia_positions=False,
                         chrono_order=True, object_info=True, directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300,
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
    gaia_images : bool, optional
        Whether to create simulated Gaia image series. The default is False.
    gaia_entries : bool, optional
        Whether to plot Gaia catalog overlays. The default is False.
    gaia_pm_vectors : bool, optional
        Whether to plot Gaia porper motion vectors. The default is False.
    gaia_pm_years : int, optional
        Number of years to scale the proper motion vectors. The default is 10.
    gaia_color : str, optional
        Gaia overlay color. The default is 'green'.
    targets : list of Finder_charts.Target objects, optional
        List of targets to be plotted on the specified image series/survey. The default is None.
        Example:
        targets = [
            Target(catalog='2MASS', epoch=2000.5, ra=0.123, dec=0.123, marker_size=10, marker_color='red', survey=Survey.TWO_MASS),
            Target(catalog='AllWISE', epoch=2010.6, ra=1.234, dec=1.234, marker_size=9, marker_color='blue', survey=Survey.ALLWISE),
        ]
        where ``catalog`` is the catalog label and ``survey`` is of type Finder_charts.Survey
    pmra : float, optional
        Proper motion in RA (mas/yr) used to propagate the crosshair position to the appropriate survey epoch. The default is None.
    pmdec : float, optional
        Proper motion in declination (mas/yr) used to propagate the crosshair position to the appropriate survey epoch. The default is None.
    ref_epoch : astropy.time.Time, optional
        Epoch of ``ra`` and ``dec`` used to propagate the crosshair position to the appropriate survey epoch. The default is None.
    crosshair_type : Finder_charts.Crosshair, , optional
        Shape of the crosshair. The default is Crosshair.CIRCLE_DOT
    propagate_gaia_positions : bool, optional
        Whether to propagate the Gaia catalog positions to the appropriate survey epoch. The default is False.
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
        def __init__(self, data, x, y, band, year_obs, wcs, overlay_ra=None, overlay_dec=None, overlay_phot=None, time_obs=None):
            self.data = data
            self.x = x
            self.y = y
            self.band = band
            self.year_obs = year_obs
            self.wcs = wcs
            self.overlay_ra = overlay_ra
            self.overlay_dec = overlay_dec
            self.overlay_phot = overlay_phot
            self.time_obs = time_obs

    def finder_charts(ra, dec, targets):

        def process_image_data(hdu, date_obs_key, date_pattern, hdu2=None):
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

                    date_obs = get_date_obs(hdu2 if hdu2 else hdu, date_obs_key, date_pattern)

                    if None not in (pmra, pmdec, ref_epoch):
                        # Propagate crosshair position
                        position = propagate_position(ra, dec, pmra, pmdec, ref_epoch, date_obs)

                    x, y = wcs.world_to_pixel(position)
                    return data, x, y, wcs, date_obs
                else:
                    return None, 0, 0, None, None
            except Exception:
                print('A problem occurred while creating an image for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                print(traceback.format_exc())
                return None, 0, 0, None, None

        def create_color_image(r, g, b, neowise=False):
            try:
                if r is None or b is None:
                    return None

                if g is None:
                    if b.shape == r.shape:
                        g = (b + r) / 2
                    else:
                        xmin = min([r.shape[0], b.shape[0]])
                        ymin = min([r.shape[1], b.shape[1]])
                        xmax = max([r.shape[0], b.shape[0]])
                        ymax = max([r.shape[1], b.shape[1]])
                        if xmax - xmin <= 2 and ymax - ymin <= 2:
                            r = r[0:xmin, 0:ymin]
                            b = b[0:xmin, 0:ymin]
                            g = (b + r) / 2
                        else:
                            return None

                if neowise:
                    vmin, vmax = get_min_max(g, lo=neowise_contrast, hi=100-neowise_contrast)
                    image = Image.fromarray(make_lupton_rgb(r, g, b, minimum=vmin, stretch=vmax-vmin, Q=0))
                    image = ImageOps.invert(image)
                else:
                    xmax = max([r.shape[0], g.shape[0], b.shape[0]])
                    ymax = max([r.shape[1], g.shape[1], b.shape[1]])
                    r = Image.fromarray(create_lupton_rgb(r)).convert('L').resize((xmax, ymax), Image.LANCZOS)
                    g = Image.fromarray(create_lupton_rgb(g)).convert('L').resize((xmax, ymax), Image.LANCZOS)
                    b = Image.fromarray(create_lupton_rgb(b)).convert('L').resize((xmax, ymax), Image.LANCZOS)
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

                if targets:
                    for t in targets:
                        if t.ra and t.dec and t.marker_size and t.marker_color:
                            if t.survey == Survey.ALL or t.survey.value in band:
                                ax.plot(t.ra, t.dec, 'o', fillstyle='none', markersize=t.marker_size, markeredgewidth=0.3, color=t.marker_color,
                                        transform=ax.get_transform('icrs'))

                if crosshair_type == Crosshair.CIRCLE_DOT:
                    mw = 0.4
                    mc = 'r'
                    ax.plot(x, y, 'o', fillstyle='none', markersize=9,  markeredgewidth=mw, color=mc)
                    ax.plot(x, y, 'o', fillstyle='none', markersize=mw, markeredgewidth=mw, color=mc)
                if crosshair_type == Crosshair.MULTI_CIRCLES:
                    mw = 0.4
                    mc = 'r'
                    ax.plot(x, y, 'o', fillstyle='none', markersize=9, markeredgewidth=mw, color=mc)
                    ax.plot(x, y, 'o', fillstyle='none', markersize=6, markeredgewidth=mw, color=mc)
                    ax.plot(x, y, 'o', fillstyle='none', markersize=3, markeredgewidth=mw, color=mc)
                if crosshair_type == Crosshair.CROSS_NO_CENTER:
                    z = len(data)*0.03
                    ms = 3
                    mw = 0.4
                    mc = 'r'
                    ax.plot(x-z, y, marker=0, markersize=ms, markeredgewidth=mw, color=mc)
                    ax.plot(x+z, y, marker=1, markersize=ms, markeredgewidth=mw, color=mc)
                    ax.plot(x, y+z, marker=2, markersize=ms, markeredgewidth=mw, color=mc)
                    ax.plot(x, y-z, marker=3, markersize=ms, markeredgewidth=mw, color=mc)

                if overlays and overlay_ra is not None and overlay_dec is not None and overlay_phot is not None:
                    ax.scatter(overlay_ra, overlay_dec, transform=ax.get_transform('icrs'), s=0/overlay_phot + 1.5,
                               edgecolor=overlay_color, facecolor='none', linewidths=0.2)

                gaia_overlay_ra = []
                gaia_overlay_dec = []
                gaia_overlay_pmra = []
                gaia_overlay_pmdec = []
                if gaia_results:
                    time_obs = image_bucket.time_obs
                    if propagate_gaia_positions and time_obs:
                        date_obs = Time(time_obs, format='jyear')
                        date_ref = Time(gaia_epoch, format='jyear')
                        for result in gaia_results:
                            gaia_ra, gaia_dec = result['ra'], result['dec']
                            gaia_pmra, gaia_pmdec = result['pmra'], result['pmdec']
                            gaia_position = propagate_position(gaia_ra, gaia_dec, gaia_pmra, gaia_pmdec, date_ref, date_obs)
                            gaia_overlay_ra.append(gaia_position.ra.value)
                            gaia_overlay_dec.append(gaia_position.dec.value)
                            gaia_overlay_pmra.append(gaia_pmra)
                            gaia_overlay_pmdec.append(gaia_pmdec)
                    else:
                        gaia_overlay_ra = gaia_results['ra']
                        gaia_overlay_dec = gaia_results['dec']
                        gaia_overlay_pmra = gaia_results['pmra']
                        gaia_overlay_pmdec = gaia_results['pmdec']

                if gaia_entries:
                    ms = np.log10(np.maximum(gaia_results['parallax'], 3))
                    ax.scatter(gaia_overlay_ra, gaia_overlay_dec, transform=ax.get_transform('icrs'), s=ms, edgecolor=gaia_color, facecolor='none', linewidths=0.2)

                if gaia_pm_vectors:
                    for gaia_ra, gaia_dec, gaia_pmra, gaia_pmdec in zip(gaia_overlay_ra, gaia_overlay_dec, gaia_overlay_pmra, gaia_overlay_pmdec):
                        coords1, coords2 = apply_PM(gaia_ra, gaia_dec, gaia_pmra, gaia_pmdec)
                        x1, y1 = wcs.world_to_pixel(coords1)
                        x2, y2 = wcs.world_to_pixel(coords2)
                        dx = x2 - x1
                        dy = y2 - y1
                        ax.quiver(x1, y1, dx, dy, angles='xy', scale_units='xy', headwidth=8, headlength=8, linewidths=0.2, scale=1.0, color=gaia_color)

                ax.text(0.03, 0.93, band, color='black', fontsize=3.0, transform=ax.transAxes,
                        bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
                ax.text(0.03, 0.04, year_obs, color='black', fontsize=3.0, transform=ax.transAxes,
                        bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
                ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec='black', transform=ax.transAxes))

                vmin, vmax = get_min_max(data)
                ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray_r')
                ax.axis('off')
            except Exception:
                print('A problem occurred while plotting an image for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
                print(traceback.format_exc())

        def apply_PM(ra, dec, pmra, pmdec):
            t1 = Time(gaia_epoch, format='jyear')
            t2 = Time(gaia_epoch + gaia_pm_years, format='jyear')
            coords1 = SkyCoord(
                ra=ra * u.deg,
                dec=dec * u.deg,
                pm_ra_cosdec=pmra * u.mas/u.yr,
                pm_dec=pmdec * u.mas/u.yr,
                obstime=t1,
                frame='icrs'
            )
            coords2 = coords1.apply_space_motion(t2)
            return coords1, coords2

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
            query_url = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv'
            payload = {
                'ra': ra,
                'dec': dec,
                'radius': radius.to(u.deg).value,
                'nStackDetections.gte': 2,
                'sort_by': ['distance']
            }
            response = requests.get(query_url, params=payload, timeout=timeout)
            table = ascii.read(response.text, format='csv')
            if len(table) > 0:
                return table
            else:
                return None

        def replace_table_values(table, value_to_replace, replacement):
            for colname in table.colnames:
                if table[colname].dtype.kind == 'f':
                    table[colname][table[colname] == value_to_replace] = replacement

        def propagate_position(ra, dec, pmra, pmdec, date_ref, date_obs):
            c = SkyCoord(
                ra=ra * u.deg,
                dec=dec * u.deg,
                pm_ra_cosdec=pmra * u.mas/u.yr,
                pm_dec=pmdec * u.mas/u.yr,
                obstime=date_ref,
                frame='icrs'
            )
            return c.apply_space_motion(date_obs)

        def get_date_obs(hdu, date_obs_key, date_pattern):
            header = hdu.header
            time_obs = header[date_obs_key]
            if date_pattern == 'MJD':
                date_obs = Time(time_obs, scale='utc', format='mjd')
            else:
                try:
                    date_obs = Time(time_obs, scale='utc', format='fits')
                except Exception:
                    try:
                        date_time = datetime.strptime(header[date_obs_key], date_pattern)
                        date_obs = Time(date_time.isoformat(sep='T', timespec='auto'), scale='utc', format='isot')
                    except Exception:
                        date_obs = Time(get_year_obs(hdu, date_obs_key, date_pattern), format='jyear')
            return date_obs

        def get_year_obs(hdu, date_obs_key, date_pattern):
            header = hdu.header
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
                filename = 'Finder_charts_' + survey + '_results_' + create_obj_name(ra, dec) + '.' + result_tables_extension
                filepath = os.path.join(directory, filename)
                result_table.write(filepath, format=result_tables_format, overwrite=True)

        def start_file(filename):
            if sys.platform == 'win32':
                os.startfile(filename)
            else:
                opener = 'open' if sys.platform == 'darwin' else 'evince'
                subprocess.call([opener, filename])

        # -------------------------------
        # Code for finder_charts function
        # -------------------------------
        gaia_epoch = 2016.0

        if targets is None:
            targets = []

        fig = plt.figure()
        fig.set_figheight(15)
        fig.set_figwidth(15)

        coords = SkyCoord(ra*u.deg, dec*u.deg)
        radius = (img_size*math.sqrt(2)/2)*u.arcsec

        surveys = []

        # Query row limit
        row_limit = 5000

        # Get Gaia catalog entries if requested
        gaia_results = None
        if gaia_entries or gaia_pm_vectors:
            Gaia.ROW_LIMIT = row_limit
            job = Gaia.cone_search_async(coords, radius=radius.to(u.deg))
            gaia_results = job.get_results()

        # DSS
        if dss:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss1_blue', img_size)
            if image:
                data, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'DSS1 B', get_year_obs(image[0], date_obs_key, date_pattern), wcs, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss1_red', img_size)
            if image:
                data, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'DSS1 R', get_year_obs(image[0], date_obs_key, date_pattern), wcs, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_blue', img_size)
            if image:
                year_b = get_year_obs(image[0], date_obs_key, date_pattern)
                b, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(b, x, y, 'DSS2 B', year_b, wcs, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_red', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(g, x, y, 'DSS2 R', year_g, wcs, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_ir', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(r, x, y, 'DSS2 IR', year_r, wcs, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'DSS2 IR-R-B', mean_obs_year, wcs, time_obs=time_obs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # 2MASS
        TWO_MASS_K = '2MASS K'
        if twomass:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'ORDATE'
            date_pattern = '%y%m%d'
            survey = []

            if overlays or save_result_tables:
                try:
                    v = Vizier(columns=['+_r', 'all'], row_limit=row_limit)
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
                b, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(b, x, y, '2MASS J', year_b, wcs, overlay_ra, overlay_dec, op1, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, '2mass', 'twomass_bands=h', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(g, x, y, '2MASS H', year_g, wcs, overlay_ra, overlay_dec, op2, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, '2mass', 'twomass_bands=k', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(r, x, y, TWO_MASS_K, year_r, wcs, overlay_ra, overlay_dec, op3, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, '2MASS K-H-J', mean_obs_year, wcs, time_obs=time_obs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)
            survey.append(None)

            surveys.append(survey)

        # Spitzer
        if spitzer:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            date_obs_key = 'MJDMEAN'
            date_pattern = 'MJD'
            survey = []

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC1', img_size)
            if image:
                year_b = get_year_obs(image[0], date_obs_key, date_pattern)
                b, x_j, y_j, wcs_j, time_obs_j = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(b, x_j, y_j, 'IRAC 1', year_b, wcs_j, time_obs=time_obs_j))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC2', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(g, x, y, 'IRAC 2', year_g, wcs, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC3', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(r, x, y, 'IRAC 3', year_r, wcs, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC4', img_size)
            if image:
                data, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'IRAC 4', get_year_obs(image[0], date_obs_key, date_pattern), wcs, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:MIPS24', img_size)
            if image:
                data, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'MIPS 24', get_year_obs(image[0], date_obs_key, date_pattern), wcs, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x_j, y_j, 'IRAC 3-2-1', mean_obs_year, wcs_j, time_obs=time_obs_j))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # WISE
        if wise:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'MIDOBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            if overlays or save_result_tables:
                try:
                    v = Vizier(columns=['+_r', 'all'], row_limit=row_limit)
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
                b, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(b, x, y, 'WISE 1', year_b, wcs, overlay_ra, overlay_dec, op1, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=2', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(g, x, y, 'WISE 2', year_g, wcs, overlay_ra, overlay_dec, op2, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=3', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(r, x, y, 'WISE 3', year_r, wcs, overlay_ra, overlay_dec, op3, time_obs=time_obs))
            else:
                survey.append(None)

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=4', img_size)
            if image:
                data, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'WISE 4', get_year_obs(image[0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op4, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'WISE 3-2-1', mean_obs_year, wcs, time_obs=time_obs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)

            surveys.append(survey)

        # UKIDSS
        if ukidss:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'UKIDSSDR11PLUS'
            if overlays or save_result_tables:
                try:
                    table = Ukidss.query_region(coords, radius=radius, database=database, programme_id='LAS')
                    if table:
                        replace_table_values(table, -999999500.0, np.nan)
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
                data, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern, images[0][0])
                survey.append(ImageBucket(data, x, y, 'UKIDSS Y', get_year_obs(images[0][0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op1, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_UKIDSS_image(ra, dec, 'J', img_size, database)
            if images:
                year_b = get_year_obs(images[0][0], date_obs_key, date_pattern)
                b, x_j, y_j, wcs_j, time_obs_j = process_image_data(images[0][1], date_obs_key, date_pattern, images[0][0])
                survey.append(ImageBucket(b, x_j, y_j, 'UKIDSS J', year_b, wcs_j, overlay_ra, overlay_dec, op2, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_UKIDSS_image(ra, dec, 'H', img_size, database)
            if images:
                year_g = get_year_obs(images[0][0], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern, images[0][0])
                survey.append(ImageBucket(g, x, y, 'UKIDSS H', year_g, wcs, overlay_ra, overlay_dec, op3, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_UKIDSS_image(ra, dec, 'K', img_size, database)
            if images:
                year_r = get_year_obs(images[0][0], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern, images[0][0])
                survey.append(ImageBucket(r, x, y, 'UKIDSS K', year_r, wcs, overlay_ra, overlay_dec, op4, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x_j, y_j, 'UKIDSS K-H-J', mean_obs_year, wcs_j, time_obs=time_obs_j))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)

            surveys.append(survey)

        # UHS
        if uhs:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
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
                        replace_table_values(table, -999999500.0, np.nan)
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
                b, x_j, y_j, wcs_j, time_obs_j = process_image_data(images[0][1], date_obs_key, date_pattern, images[0][0])
                survey.append(ImageBucket(b, x_j, y_j, 'UHS J', year_b, wcs_j, overlay_ra, overlay_dec, op2, time_obs=time_obs_j))
            else:
                survey.append(None)

            images = get_UKIDSS_image(ra, dec, 'K', img_size, database)
            if images:
                year_r = get_year_obs(images[0][0], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern, images[0][0])
                survey.append(ImageBucket(r, x, y, 'UHS K', year_r, wcs, overlay_ra, overlay_dec, op4, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x_j, y_j, 'UHS K-J', mean_obs_year, wcs_j, time_obs=time_obs_j))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)
            survey.append(None)
            survey.append(None)

            surveys.append(survey)

        # Override deprecated VSA access URLs
        BASE_URL = 'http://vsa.roe.ac.uk:8080/vdfs/'
        Vsa.LOGIN_URL = BASE_URL + "DBLogin"
        Vsa.IMAGE_URL = BASE_URL + "GetImage"
        Vsa.ARCHIVE_URL = BASE_URL + "ImageList"
        Vsa.REGION_URL = BASE_URL + "WSASQL"

        # VHS
        if vhs:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'VHSDR6'
            if overlays or save_result_tables:
                try:
                    table = Vsa.query_region(coords, radius=radius, database=database, programme_id='VHS')
                    if table:
                        replace_table_values(table, -999999500.0, np.nan)
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
                data, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'VHS Y', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op1, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'J', img_size, database, frame_type)
            if images:
                year_b = get_year_obs(images[0][1], date_obs_key, date_pattern)
                b, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(b, x, y, 'VHS J', year_b, wcs, overlay_ra, overlay_dec, op2, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'H', img_size, database, frame_type)
            if images:
                year_g = get_year_obs(images[0][1], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(g, x, y, 'VHS H', year_g, wcs, overlay_ra, overlay_dec, op3, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Ks', img_size, database, frame_type)
            if images:
                year_r = get_year_obs(images[0][1], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(r, x, y, 'VHS K', year_r, wcs, overlay_ra, overlay_dec, op4, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'VHS K-H-J', mean_obs_year, wcs, time_obs=time_obs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            survey.append(None)

            surveys.append(survey)

        # VVV
        if vvv:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'VVVDR5'
            if overlays or save_result_tables:
                try:
                    table = Vsa.query_region(coords, radius=radius, database=database, programme_id='VVV')
                    if table:
                        replace_table_values(table, -999999500.0, np.nan)
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
                data, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'VVV Z', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op1, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Y', img_size, database, frame_type)
            if images:
                data, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'VVV Y', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op2, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'J', img_size, database, frame_type)
            if images:
                year_b = get_year_obs(images[0][1], date_obs_key, date_pattern)
                b, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(b, x, y, 'VVV J', year_b, wcs, overlay_ra, overlay_dec, op3, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'H', img_size, database, frame_type)
            if images:
                year_g = get_year_obs(images[0][1], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(g, x, y, 'VVV H', year_g, wcs, overlay_ra, overlay_dec, op4, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Ks', img_size, database, frame_type)
            if images:
                year_r = get_year_obs(images[0][1], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(r, x, y, 'VVV K', year_r, wcs, overlay_ra, overlay_dec, op5, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'VVV K-H-J', mean_obs_year, wcs, time_obs=time_obs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # VIKING
        if viking:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
            mean_obs_year = 0
            year_r = year_g = year_b = np.nan
            overlay_ra = overlay_dec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'
            survey = []

            database = 'VIKINGDR5'
            if overlays or save_result_tables:
                try:
                    table = Vsa.query_region(coords, radius=radius, database=database, programme_id='VIKING')
                    if table:
                        replace_table_values(table, -999999500.0, np.nan)
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
                data, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'VIKING Z', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op1, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Y', img_size, database, frame_type)
            if images:
                data, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, 'VIKING Y', get_year_obs(images[0][1], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op2, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'J', img_size, database, frame_type)
            if images:
                year_b = get_year_obs(images[0][1], date_obs_key, date_pattern)
                b, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(b, x, y, 'VIKING J', year_b, wcs, overlay_ra, overlay_dec, op3, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'H', img_size, database, frame_type)
            if images:
                year_g = get_year_obs(images[0][1], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(g, x, y, 'VIKING H', year_g, wcs, overlay_ra, overlay_dec, op4, time_obs=time_obs))
            else:
                survey.append(None)

            images = get_VSA_image(ra, dec, 'Ks', img_size, database, frame_type)
            if images:
                year_r = get_year_obs(images[0][1], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(images[0][1], date_obs_key, date_pattern)
                survey.append(ImageBucket(r, x, y, 'VIKING K', year_r, wcs, overlay_ra, overlay_dec, op5, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'VIKING K-H-J', mean_obs_year, wcs, time_obs=time_obs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # PS1
        if ps1:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
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
                                replace_table_values(table, -999.0, np.nan)
                                overlay_ra = table['raMean']
                                overlay_dec = table['decMean']
                                op1 = table['gMeanPSFMag']
                                op2 = table['rMeanPSFMag']
                                op3 = table['iMeanPSFMag']
                                op4 = table['zMeanPSFMag']
                                op5 = table['yMeanPSFMag']
                                save_catalog_search_results(table, 'PS1_DR2', ra, dec)
                        except Exception:
                            print('A problem occurred while downloading Pan-STARRS catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                            print(traceback.format_exc())

                    if 'g' in images:
                        year_b = get_year_obs(images['g'][0], date_obs_key, date_pattern)
                        b, x, y, wcs, time_obs = process_image_data(images['g'][0], date_obs_key, date_pattern)
                        survey.append(ImageBucket(b, x, y, 'PS1 g', year_b, wcs, overlay_ra, overlay_dec, op1, time_obs=time_obs))
                    else:
                        survey.append(None)

                    if 'r' in images:
                        data, x, y, wcs, time_obs = process_image_data(images['r'][0], date_obs_key, date_pattern)
                        survey.append(ImageBucket(data, x, y, 'PS1 r', get_year_obs(images['r'][0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op2, time_obs=time_obs))
                    else:
                        survey.append(None)

                    if 'i' in images:
                        year_g = get_year_obs(images['i'][0], date_obs_key, date_pattern)
                        g, x, y, wcs, time_obs = process_image_data(images['i'][0], date_obs_key, date_pattern)
                        survey.append(ImageBucket(g, x, y, 'PS1 i', year_g, wcs, overlay_ra, overlay_dec, op3, time_obs=time_obs))
                    else:
                        survey.append(None)

                    if 'z' in images:
                        data, x, y, wcs, time_obs = process_image_data(images['z'][0], date_obs_key, date_pattern)
                        survey.append(ImageBucket(data, x, y, 'PS1 z', get_year_obs(images['z'][0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op4, time_obs=time_obs))
                    else:
                        survey.append(None)

                    if 'y' in images:
                        year_r = get_year_obs(images['y'][0], date_obs_key, date_pattern)
                        r, x, y, wcs, time_obs = process_image_data(images['y'][0], date_obs_key, date_pattern)
                        survey.append(ImageBucket(r, x, y, 'PS1 y', year_r, wcs, overlay_ra, overlay_dec, op5, time_obs=time_obs))
                    else:
                        survey.append(None)

                    if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                        mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                        survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'PS1 y-i-g', mean_obs_year, wcs, time_obs=time_obs))
                    else:
                        survey.insert(0, None)

                    survey.insert(0, mean_obs_year)

                    surveys.append(survey)

        # DECam
        if decam:
            x = y = x_j = y_j = 0
            r = g = b = wcs = time_obs = wcs_j = time_obs_j = None
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
                    v = Vizier(columns=['+_r', 'all'], row_limit=row_limit)
                    tables = v.query_region(coords, radius=radius, catalog='II/371/des_dr2')
                    if tables:
                        table = tables[0]
                        overlay_ra = table['RA_ICRS']
                        overlay_dec = table['DE_ICRS']
                        op1 = table['gmag']
                        op2 = table['rmag']
                        op3 = table['imag']
                        op4 = table['zmag']
                        op5 = table['Ymag']
                        for col in table.colnames:
                            if '/' in col:
                                table.rename_column(col, col.replace('/', '_'))
                        save_catalog_search_results(table, 'DES_DR2', ra, dec)
                except Exception:
                    print('A problem occurred while downloading DES catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            image, instrument = get_DECam_image(ra, dec, 'g', img_size)
            if image:
                year_b = get_year_obs(image[0], date_obs_key, date_pattern)
                b, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(b, x, y, instrument + ' g', year_b, wcs, overlay_ra, overlay_dec, op1, time_obs=time_obs))
            else:
                survey.append(None)

            image, instrument = get_DECam_image(ra, dec, 'r', img_size)
            if image:
                data, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, instrument + ' r', get_year_obs(image[0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op2, time_obs=time_obs))
            else:
                survey.append(None)

            image, instrument = get_DECam_image(ra, dec, 'i', img_size)
            if image:
                year_g = get_year_obs(image[0], date_obs_key, date_pattern)
                g, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(g, x, y, instrument + ' i', year_g, wcs, overlay_ra, overlay_dec, op3, time_obs=time_obs))
            else:
                survey.append(None)

            image, instrument = get_DECam_image(ra, dec, 'z', img_size)
            if image:
                data, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(data, x, y, instrument + ' z', get_year_obs(image[0], date_obs_key, date_pattern), wcs, overlay_ra, overlay_dec, op4, time_obs=time_obs))
            else:
                survey.append(None)

            image, instrument = get_DECam_image(ra, dec, 'Y', img_size)
            if image:
                year_r = get_year_obs(image[0], date_obs_key, date_pattern)
                r, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                survey.append(ImageBucket(r, x, y, instrument + ' Y', year_r, wcs, overlay_ra, overlay_dec, op5, time_obs=time_obs))
            else:
                survey.append(None)

            if np.isfinite(year_r) or np.isfinite(year_g) or np.isfinite(year_b):
                mean_obs_year = round(np.nanmean([year_r, year_g, year_b]), 1)
                survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, instrument + ' Y-i-g', mean_obs_year, wcs, time_obs=time_obs))
            else:
                survey.insert(0, None)

            survey.insert(0, mean_obs_year)

            surveys.append(survey)

        # Gaia simulated images
        if gaia_images:
            survey = []

            b, gaia_wcs, obs_coords = create_gaia_finder_chart(ra, dec, size=img_size, band='Bp', epoch=gaia_epoch)
            x, y = gaia_wcs.world_to_pixel(obs_coords)
            survey.append(ImageBucket(b, x, y, 'Gaia BP (simulated)', gaia_epoch, gaia_wcs, time_obs=gaia_epoch))

            g, gaia_wcs, obs_coords = create_gaia_finder_chart(ra, dec, size=img_size, band='G', epoch=gaia_epoch)
            x, y = gaia_wcs.world_to_pixel(obs_coords)
            survey.append(ImageBucket(g, x, y, 'Gaia G (simulated)', gaia_epoch, gaia_wcs, time_obs=gaia_epoch))

            r, gaia_wcs, obs_coords = create_gaia_finder_chart(ra, dec, size=img_size, band='Rp', epoch=gaia_epoch)
            x, y = gaia_wcs.world_to_pixel(obs_coords)
            survey.append(ImageBucket(r, x, y, 'Gaia RP (simulated)', gaia_epoch, gaia_wcs, time_obs=gaia_epoch))

            survey.insert(0, ImageBucket(create_color_image(r, g, b), x, y, 'Gaia RP-G-BP', gaia_epoch, gaia_wcs, time_obs=gaia_epoch))

            survey.insert(0, gaia_epoch)

            survey.append(None)
            survey.append(None)

            surveys.append(survey)

        # WISE time series
        if neowise:
            imageW1 = get_neowise_image(ra, dec, epoch=0, band=1, size=img_size)
            imageW2 = get_neowise_image(ra, dec, epoch=0, band=2, size=img_size)
            if imageW1 and imageW2:
                hdu = imageW1[0]
                header = hdu.header
                mjdmean = (header['MJDMIN']+header['MJDMAX'])/2
                header['MJDMEAN'] = mjdmean
                prev_year = get_year_from_mjd(mjdmean)
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
                        h = hdu.header
                        mjdmean = (h['MJDMIN']+h['MJDMAX'])/2
                        year = get_year_from_mjd(mjdmean)
                        if year == prev_year:
                            j += 1
                        else:
                            hduW1 = fits.PrimaryHDU(data=dataW1/j, header=header)
                            hduW2 = fits.PrimaryHDU(data=dataW2/j, header=header)
                            images.append((hduW1, hduW2, prev_year))
                            j = 1
                        header = hdu.header
                        header['MJDMEAN'] = mjdmean
                        dataW1 = hdu.data
                        dataW2 = imageW2[0].data
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
                    date_obs_key = 'MJDMEAN'
                    date_pattern = 'MJD'
                    w1, x, y, wcs, time_obs = process_image_data(image[0], date_obs_key, date_pattern)
                    w2, x, y, wcs, time_obs = process_image_data(image[1], date_obs_key, date_pattern)
                    if w1 is not None and w2 is not None:
                        wise_epoch = image[2]
                        survey.append(ImageBucket(create_color_image(w1, (w1+w2)/2, w2, neowise=True), x, y, 'W1-W2', wise_epoch, wcs, time_obs=time_obs))

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
            hmsdms = coords.to_string('hmsdms', sep=':', precision=2)
            hms = hmsdms[0:11]
            dms = hmsdms[12:24] if dec < 0 else hmsdms[13:24]
            ax.text(0.05, 0.86, r'$\alpha$ = ' + str(round(coords.ra.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.73, '(' + hms + ')', fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.60, r'$\delta$ = ' + str(round(coords.dec.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.47, '(' + dms + ')', fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.34, '$l$ = ' + str(round(coords.galactic.l.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.21, '$b$ = ' + str(round(coords.galactic.b.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, 0.08, 'Size = ' + str(int(img_size)) + ' arcsec', fontsize=fontsize, transform=ax.transAxes)
            ax.axis('off')

            # Info text cont'd
            text_y = 0.85
            fontsize = 3.0
            ax = fig.add_subplot(rows, cols, info_idx + 1)
            ax.axis('off')

            if gaia_entries:
                targets.insert(0, Target(catalog='Gaia DR3 entries', epoch='', ra=None, dec=None, marker_size=None, marker_color=gaia_color, survey=Survey.ALL))

            for t in targets:
                if t.survey == Survey.ALL:
                    ax.text(0.05, text_y, t.catalog + ' ' + str(t.epoch), fontsize=fontsize, transform=ax.transAxes, color=t.marker_color)
                    text_y -= 0.1
                if text_y < 0:
                    break

        # Save the finder charts
        filename = 'Finder_charts_' + create_obj_name(ra, dec) + '.' + file_format
        filepath = os.path.join(directory, filename)
        plt.subplots_adjust(wspace=0, hspace=0.05, right=0.43)
        plt.savefig(filepath, dpi=600, bbox_inches='tight', format=file_format)
        plt.close()

        if open_file:
            start_file(filepath)

    # --------------------------------------
    # Code for create_finder_charts function
    # --------------------------------------

    warnings.simplefilter('ignore', category=Warning)

    if np.isscalar(ra) and np.isscalar(dec):
        finder_charts(ra, dec, targets)
    else:
        if not isinstance(ra, (collections.Sequence, np.ndarray)):
            raise Exception('Ra must be either a sequence or a numpy array')

        if not isinstance(dec, (collections.Sequence, np.ndarray)):
            raise Exception('Dec must be either a sequence or a numpy array')

        if len(ra) != len(dec):
            raise Exception('Ra and Dec must have the same length: len(ra)={ra}, len(dec)={dec}'.format(ra=len(ra), dec=len(dec)))

        for i in range(len(ra)):
            try:
                finder_charts(ra[i], dec[i], targets)
            except Exception:
                print('A problem occurred while creating finder charts for object ra={ra}, dec={dec}'.format(ra=ra[i], dec=dec[i]))
                print(traceback.format_exc())
