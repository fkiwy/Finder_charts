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
from datetime import datetime
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
from matplotlib.patches import BoxStyle, Rectangle
from astropy.io import ascii
# from astropy.stats import sigma_clip
from astropy.visualization import make_lupton_rgb
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning
from astropy.utils.data import download_file
from astroquery.ukidss import Ukidss
from astroquery.vsa import Vsa
from astroquery.vizier import Vizier
from astroquery.mast import Catalogs
import astropy.units as u
from astropy.coordinates import SkyCoord
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from pyvo.dal import sia


def create_finder_charts(ra, dec, img_size=100, overlays=False, overlay_color='red', dss=True, twomass=True, spitzer=True, wise=True,
                         ukidss=True, vhs=True, ps1=True, decam=True, neowise=True, neowise_contrast=3, directory=tempfile.gettempdir(),
                         cache=True, show_progress=True, timeout=300 ,open_pdf=True, file_format='pdf'):

    def finder_charts(ra, dec):

        def create_image(hdu, img_idx, band, year_obs, ora=None, odec=None, ophot=None):
            try:
                data = hdu.data

                nanValues = np.count_nonzero(np.isnan(data))
                totValues = np.count_nonzero(data)

                if nanValues < totValues * 0.5:
                    wcs, shape = find_optimal_celestial_wcs([hdu], frame='icrs')
                    data, _ = reproject_interp(hdu, wcs, shape_out=shape)

                    position = SkyCoord(ra*u.deg, dec*u.deg)
                    cutout = Cutout2D(data, position, img_size*u.arcsec, wcs=wcs, mode='partial')
                    data = cutout.data
                    wcs = cutout.wcs

                    ax = fig.add_subplot(rows, cols, img_idx, projection=wcs)
                    x, y = wcs.world_to_pixel(position)
                    ax.plot(x, y, 'ro', fillstyle='none', markersize=7, markeredgewidth=0.2)
                    ax.plot(x, y, 'ro', fillstyle='none', markersize=0.2, markeredgewidth=0.2)
                    ax.text(0.04, 0.91, band, color='black', fontsize=1.8, transform=ax.transAxes,
                            bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
                    ax.text(0.04, 0.05, year_obs, color='black', fontsize=1.8, transform=ax.transAxes,
                            bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
                    ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec='black', transform=ax.transAxes))

                    if ora is not None and odec is not None and ophot is not None:
                        ax.scatter(ora, odec, transform=ax.get_transform('icrs'), s=0/ophot + 1.0,
                                   edgecolor=overlay_color, facecolor='none', linewidths=0.2)

                    vmin, vmax = get_min_max(data)
                    ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray_r')
                    ax.axis('off')
                    return 1, data, x, y
                else:
                    return 0, None, 0, 0
            except:
                print('A problem occurred while creating an image for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
                print(traceback.format_exc())
                return 0, None, 0, 0

        def create_color_image(r, g, b, img_idx, band, x, y, neowise=False, year_obs=None):
            try:
                if r is None or g is None or b is None:
                    return

                xmin = min([r.shape[0], g.shape[0], b.shape[0]])
                ymin = min([r.shape[1], g.shape[1], b.shape[1]])

                xmax = max([r.shape[0], g.shape[0], b.shape[0]])
                ymax = max([r.shape[1], g.shape[1], b.shape[1]])

                if xmax - xmin > 2 or ymax - ymin > 2:
                    print('Array shapes too different to create a color image for', band, '-> R shape:',
                          r.shape, 'G shape:', g.shape, 'B shape:', b.shape)
                    return

                r = r[0:xmin, 0:ymin]
                g = g[0:xmin, 0:ymin]
                b = b[0:xmin, 0:ymin]

                if neowise:
                    # _, vmin, vmax = sigma_clip(g, sigma_lower=1, sigma_upper=5, maxiters=None, return_bounds=True)
                    vmin, vmax = get_min_max(g, lo=neowise_contrast, hi=100-neowise_contrast)
                    rgb = Image.fromarray(make_lupton_rgb(r, g, b, stretch=vmax-vmin, Q=0, minimum=vmin))
                else:
                    r = Image.fromarray(create_lupton_rgb(r)).convert("L")
                    g = Image.fromarray(create_lupton_rgb(g)).convert("L")
                    b = Image.fromarray(create_lupton_rgb(b)).convert("L")
                    rgb = Image.merge("RGB", (r, g, b))

                ax = fig.add_subplot(rows, cols, img_idx)
                ax.plot(x, y, 'ro', fillstyle='none', markersize=7, markeredgewidth=0.2)
                ax.plot(x, y, 'ro', fillstyle='none', markersize=0.2, markeredgewidth=0.2)
                ax.text(0.04, 0.91, band, color='black', fontsize=1.8, transform=ax.transAxes,
                        bbox=dict(facecolor='white', alpha=0.7, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
                ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec='black', transform=ax.transAxes))

                if year_obs:
                    ax.text(0.04, 0.05, year_obs, color='black', fontsize=1.8, transform=ax.transAxes,
                            bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
                    rgb = ImageOps.invert(rgb)

                ax.imshow(rgb, origin='lower')
                ax.axis('off')
            except:
                print('A problem occurred while creating a color image for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
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
            except:
                return None

        def get_IRSA_image(ra, dec, survey, band, size):
            download_url = 'https://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?mode=getImage&file_type=fits&reproject=false&' + \
                'RA={ra}&DEC={dec}&subsetsize={size}&survey={survey}&{band}'
            download_url = download_url.format(ra=ra, dec=dec, size=size/60, survey=survey, band=band)
            try:
                return fits.open(download_file(download_url, cache=cache, show_progress=show_progress, timeout=timeout))
            except:
                return None

        def get_UKIDSS_image(ra, dec, band, size, database):
            try:
                return Ukidss.get_images(SkyCoord(ra, dec, unit=(u.deg, u.deg)), image_width=size * u.arcsec, database=database,
                                         waveband=band, frame_type='stack', verbose=False, show_progress=show_progress)
            except:
                print('A problem occurred while downloading UKIDSS images for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
                print(traceback.format_exc())
                return None

        def get_VHS_image(ra, dec, band, size):
            try:
                return Vsa.get_images(SkyCoord(ra, dec, unit=(u.deg, u.deg)), image_width=size * u.arcsec, database='VHSDR6',
                                      waveband=band, frame_type='tilestack', verbose=False, show_progress=show_progress)
            except:
                print('A problem occurred while downloading VHS images for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
                print(traceback.format_exc())
                return None

        def get_DECam_image(ra, dec, band, size):
            try:
                fov = (size * u.arcsec).to(u.deg)
                svc = sia.SIAService('https://datalab.noirlab.edu/sia')
                result = svc.search((ra, dec), (fov/np.cos(dec*np.pi/180), fov), verbosity=2).to_table()
                select = np.char.startswith(result['obs_bandpass'].astype(str), band) & \
                    (result['proctype'] == 'Stack') & (result['prodtype'] == 'image')  # & (result['instrument_name'] == 'DECam')
                table = result[select]
                if len(table) > 0:
                    row = table[np.argmax(table['exptime'].data.data.astype('float'))]  # get image with longest exposure time
                    image = fits.open(download_file(row['access_url'], cache=cache, show_progress=show_progress, timeout=timeout))
                    instrument = row['instrument_name']
                    return image, instrument
                else:
                    return None, None
            except:
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
                print('A problem occurred while downloading Noirlab catalog overlays for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                print(traceback.format_exc())
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

        def open_file(filename):
            if sys.platform == "win32":
                os.startfile(filename)
            else:
                opener = "open" if sys.platform == "darwin" else "evince"
                subprocess.call([opener, filename])

        # -------------------------------
        # Code for finder_charts function
        # -------------------------------
        fig = plt.figure()
        fig.set_figheight(5)
        fig.set_figwidth(5)
        plt.subplots_adjust(wspace=0, hspace=0.05, right=0.575)

        img_idx = 1
        coords = SkyCoord(ra*u.deg, dec*u.deg)
        radius = (img_size*math.sqrt(2)/2)*u.arcsec

        # DSS
        if dss:
            i = 0
            r = g = b = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss1_blue', img_size)
            if image:
                j, _, _, _ = create_image(image[0], img_idx, 'DSS1 B', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss1_red', img_size)
            if image:
                j, _, _, _ = create_image(image[0], img_idx, 'DSS1 R', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_blue', img_size)
            if image:
                j, b, x, y = create_image(image[0], img_idx, 'DSS2 B', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_red', img_size)
            if image:
                j, g, x, y = create_image(image[0], img_idx, 'DSS2 R', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'dss', 'dss_bands=poss2ukstu_ir', img_size)
            if image:
                j, r, x, y = create_image(image[0], img_idx, 'DSS2 IR', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            create_color_image(r, g, b, img_idx, 'IR-R-B', x, y)
            img_idx += 1

            if i == 0:
                img_idx -= cols

        # 2MASS
        if twomass:
            i = 0
            r = g = b = None
            ora = odec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'ORDATE'
            date_pattern = '%y%m%d'

            if overlays:
                try:
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'Jmag', 'Hmag',  'Kmag'])
                    tables = v.query_region(coords, radius=radius, catalog='II/246/out')
                    if tables:
                        table = tables[0]
                        ora = table['RAJ2000']
                        odec = table['DEJ2000']
                        op1 = table['Jmag']
                        op2 = table['Hmag']
                        op3 = table['Kmag']
                except:
                    print('A problem occurred while downloading 2MASS catalog overlays for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            image = get_IRSA_image(ra, dec, '2mass', 'twomass_bands=j', img_size)
            if image:
                j, b, x, y = create_image(image[0], img_idx, '2MASS J', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op1)
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, '2mass', 'twomass_bands=h', img_size)
            if image:
                j, g, x, y = create_image(image[0], img_idx, '2MASS H', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op2)
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, '2mass', 'twomass_bands=k', img_size)
            if image:
                j, r, x, y = create_image(image[0], img_idx, '2MASS K', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op3)
                i += j
            img_idx += 1

            create_color_image(r, g, b, img_idx, 'K-H-J', x, y)
            img_idx += 1

            # if i == 0:
            #    img_idx -= cols
        else:
            img_idx += 4

        # Info text
        fontsize = 2.4
        ax = fig.add_subplot(rows, cols, img_idx)
        ax.text(0.1, 0.7, r'$\alpha$ = ' + str(round(coords.ra.value, 6)), fontsize=fontsize, transform=ax.transAxes)
        ax.text(0.1, 0.6, r'$\delta$ = ' + str(round(coords.dec.value, 6)), fontsize=fontsize, transform=ax.transAxes)
        ax.text(0.1, 0.5, '$l$ = ' + str(round(coords.galactic.l.value, 6)), fontsize=fontsize, transform=ax.transAxes)
        ax.text(0.1, 0.4, '$b$ = ' + str(round(coords.galactic.b.value, 6)), fontsize=fontsize, transform=ax.transAxes)
        ax.axis('off')
        img_idx += 1

        # Info text cont'd
        hmsdms = coords.to_string('hmsdms', sep=':', precision=2)
        hms = hmsdms[0:11]
        dms = hmsdms[12:24] if dec < 0 else hmsdms[13:24]
        ax = fig.add_subplot(rows, cols, img_idx)
        ax.text(-0.1, 0.7, '(' + hms + ')', fontsize=fontsize, transform=ax.transAxes)
        ax.text(-0.1, 0.6, '(' + dms + ')', fontsize=fontsize, transform=ax.transAxes)
        ax.text(-0.1, 0.5, 'Size = ' + str(int(img_size)) + ' arcsec', fontsize=fontsize, transform=ax.transAxes)
        ax.text(-0.1, 0.4, 'North up, East left', fontsize=fontsize, transform=ax.transAxes)
        ax.axis('off')
        img_idx += 1

        # Spitzer
        if spitzer:
            i = 0
            r = g = b = None
            date_obs_key = 'MJDMEAN'
            date_pattern = 'MJD'

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC1', img_size)
            if image:
                j, b, x, y = create_image(image[0], img_idx, 'IRAC1', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC2', img_size)
            if image:
                j, g, x, y = create_image(image[0], img_idx, 'IRAC2', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC3', img_size)
            if image:
                j, r, x, y = create_image(image[0], img_idx, 'IRAC3', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:IRAC4', img_size)
            if image:
                j, _, _, _ = create_image(image[0], img_idx, 'IRAC4', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'seip', 'seip_bands=spitzer.seip_science:MIPS24', img_size)
            if image:
                j, _, _, _ = create_image(image[0], img_idx, 'MIPS24', get_year_obs(image[0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            create_color_image(r, g, b, img_idx, 'CH3-CH2-CH1', x, y)
            img_idx += 1

            if i == 0:
                img_idx -= cols

        # WISE
        if wise:
            i = 0
            r = g = b = None
            ora = odec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'MIDOBS'
            date_pattern = '%Y-%m-%d'

            if overlays:
                try:
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'W1mag', 'W2mag', 'W3mag', 'W4mag'])
                    tables = v.query_region(coords, radius=radius, catalog='II/328/allwise')
                    if tables:
                        table = tables[0]
                        ora = table['RAJ2000']
                        odec = table['DEJ2000']
                        op1 = table['W1mag']
                        op2 = table['W2mag']
                        op3 = table['W3mag']
                        op4 = table['W4mag']
                except:
                    print('A problem occurred while downloading 2MASS catalog overlays for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=1', img_size)
            if image:
                j, b, x, y = create_image(image[0], img_idx, 'W1', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op1)
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=2', img_size)
            if image:
                j, g, x, y = create_image(image[0], img_idx, 'W2', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op2)
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=3', img_size)
            if image:
                j, r, x, y = create_image(image[0], img_idx, 'W3', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op3)
                i += j
            img_idx += 1

            image = get_IRSA_image(ra, dec, 'wise', 'wise_bands=4', img_size)
            if image:
                j, _, _, _ = create_image(image[0], img_idx, 'W4', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op4)
                i += j
            img_idx += 1

            create_color_image(r, g, b, img_idx, 'W3-W2-W1', x, y)
            img_idx += 2

            if i == 0:
                img_idx -= cols

        # UKIDSS
        if ukidss:
            i = 0
            r = g = b = None
            ora = odec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'

            if overlays:
                try:
                    table = Ukidss.query_region(coords, radius, database='UKIDSSDR11PLUS', programme_id='LAS')
                    if table:
                        ora = table['ra']
                        odec = table['dec']
                        op1 = table['yAperMag3']
                        op2 = table['jAperMag3']
                        op3 = table['hAperMag3']
                        op4 = table['kAperMag3']
                except:
                    print('A problem occurred while downloading 2MASS catalog overlays for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            database = 'UKIDSSDR11PLUS'
            images = get_UKIDSS_image(ra, dec, 'Y', img_size, database)
            if images:
                j, _, _, _ = create_image(images[0][1], img_idx, 'UKIDSS Y', get_year_obs(images[0][0], date_obs_key, date_pattern), ora, odec, op1)
                i += j
            img_idx += 1

            images = get_UKIDSS_image(ra, dec, 'J', img_size, database)
            if images:
                j, _, _, _ = create_image(images[0][1], img_idx, 'UKIDSS J', get_year_obs(images[0][0], date_obs_key, date_pattern), ora, odec, op2)
                i += j
            img_idx += 1

            images = get_UKIDSS_image(ra, dec, 'H', img_size, database)
            if images:
                j, g, x, y = create_image(images[0][1], img_idx, 'UKIDSS H', get_year_obs(images[0][0], date_obs_key, date_pattern), ora, odec, op3)
                i += j
            img_idx += 1

            images = get_UKIDSS_image(ra, dec, 'K', img_size, database)
            if images:
                j, r, x, y = create_image(images[0][1], img_idx, 'UKIDSS K', get_year_obs(images[0][0], date_obs_key, date_pattern), ora, odec, op4)
                i += j
            img_idx += 1

            create_color_image(r, g, g, img_idx, 'K-H', x, y)
            img_idx += 1

            # UHS
            images = get_UKIDSS_image(ra, dec, 'J', img_size, 'UHSDR1')
            if images:
                j, _, _, _ = create_image(images[0][1], img_idx, 'UHS J', get_year_obs(images[0][0], date_obs_key, date_pattern))
                i += j
            img_idx += 1

            if i == 0:
                img_idx -= cols

        # VHS
        if vhs:
            i = 0
            r = g = b = None
            ora = odec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'

            if overlays:
                try:
                    table = Vsa.query_region(coords, radius, database='VHSDR6', programme_id='VHS')
                    if table:
                        ora = table['ra']
                        odec = table['dec']
                        op1 = table['yAperMag3']
                        op2 = table['jAperMag3']
                        op3 = table['hAperMag3']
                        op4 = table['ksAperMag3']
                except:
                    print('A problem occurred while downloading 2MASS catalog overlays for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                    print(traceback.format_exc())

            images = get_VHS_image(ra, dec, 'Y', img_size)
            if images:
                j, _, _, _ = create_image(images[0][1], img_idx, 'VHS Y', get_year_obs(images[0][1], date_obs_key, date_pattern), ora, odec, op1)
                i += j
            img_idx += 1

            images = get_VHS_image(ra, dec, 'J', img_size)
            if images:
                j, b, x, y = create_image(images[0][1], img_idx, 'VHS J', get_year_obs(images[0][1], date_obs_key, date_pattern), ora, odec, op2)
                i += j
            img_idx += 1

            images = get_VHS_image(ra, dec, 'H', img_size)
            if images:
                j, g, x, y = create_image(images[0][1], img_idx, 'VHS H', get_year_obs(images[0][1], date_obs_key, date_pattern), ora, odec, op3)
                i += j
            img_idx += 1

            images = get_VHS_image(ra, dec, 'Ks', img_size)
            if images:
                j, r, x, y = create_image(images[0][1], img_idx, 'VHS K', get_year_obs(images[0][1], date_obs_key, date_pattern), ora, odec, op4)
                i += j
            img_idx += 1

            create_color_image(r, g, b, img_idx, 'K-H-J', x, y)
            img_idx += 2

            if i == 0:
                img_idx -= cols

        # PS1
        if ps1:
            r = g = b = None
            ora = odec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'MJD-OBS'
            date_pattern = 'MJD'

            try:
                query_url = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py'
                payload = {
                    'ra': ra,
                    'dec': dec,
                    'filters': 'grizy',
                    'sep': 'comma'
                }
                text = requests.get(query_url, params=payload, timeout=timeout).text
            except:
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
                    except:
                        print('A problem occurred while downloading Pan-STARRS images for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                        print(traceback.format_exc())

                if images:
                    if overlays:
                        try:
                            columns = ['raMean', 'decMean', 'gMeanPSFMag',  'rMeanPSFMag', 'iMeanPSFMag', 'zMeanPSFMag', 'yMeanPSFMag']
                            table = Catalogs.query_region(coords, radius=radius, catalog='Panstarrs', data_release='dr2', table='mean',
                                                          nStackDetections=[("gte", 2)], columns=columns)
                            if table:
                                ora = table['raMean']
                                odec = table['decMean']
                                op1 = table['gMeanPSFMag']
                                op2 = table['rMeanPSFMag']
                                op3 = table['iMeanPSFMag']
                                op4 = table['zMeanPSFMag']
                                op5 = table['yMeanPSFMag']
                        except:
                            print('A problem occurred while downloading 2MASS catalog overlays for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
                            print(traceback.format_exc())

                    _, b, x, y = create_image(images['g'][0], img_idx, 'PS1 g', get_year_obs(images['g'][0], date_obs_key, date_pattern), ora, odec, op1)
                    img_idx += 1

                    _, _, _, _ = create_image(images['r'][0], img_idx, 'PS1 r', get_year_obs(images['r'][0], date_obs_key, date_pattern), ora, odec, op2)
                    img_idx += 1

                    _, g, x, y = create_image(images['i'][0], img_idx, 'PS1 i', get_year_obs(images['i'][0], date_obs_key, date_pattern), ora, odec, op3)
                    img_idx += 1

                    _, _, _, _ = create_image(images['z'][0], img_idx, 'PS1 z', get_year_obs(images['z'][0], date_obs_key, date_pattern), ora, odec, op4)
                    img_idx += 1

                    _, r, x, y = create_image(images['y'][0], img_idx, 'PS1 y', get_year_obs(images['y'][0], date_obs_key, date_pattern), ora, odec, op5)
                    img_idx += 1

                    create_color_image(r, g, b, img_idx, 'y-i-g', x, y)
                    img_idx += 1

        # DECam
        if decam:
            i = 0
            r = g = b = None
            ora = odec = op1 = op2 = op3 = op4 = op5 = None
            date_obs_key = 'DATE-OBS'
            date_pattern = '%Y-%m-%d'

            if overlays:
                table = search_noirlab(ra, dec, radius)
                if table:
                    ora = table['ra']
                    odec = table['dec']
                    op1 = table['gmag']
                    op2 = table['rmag']
                    op3 = table['imag']
                    op4 = table['zmag']
                    op5 = table['ymag']

            image, instrument = get_DECam_image(ra, dec, 'g', img_size)
            if image:
                j, b, x, y = create_image(image[0], img_idx, instrument + ' g', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op5)
                i += j
            img_idx += 1

            image, instrument = get_DECam_image(ra, dec, 'r', img_size)
            if image:
                j, g, x, y = create_image(image[0], img_idx, instrument + ' r', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op5)
                i += j
            img_idx += 1

            image, instrument = get_DECam_image(ra, dec, 'i', img_size)
            if image:
                j, _, _, _ = create_image(image[0], img_idx, instrument + ' i', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op5)
                i += j
            img_idx += 1

            image, instrument = get_DECam_image(ra, dec, 'z', img_size)
            if image:
                j, r, x, y = create_image(image[0], img_idx, instrument + ' z', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op5)
                i += j
            img_idx += 1

            image, instrument = get_DECam_image(ra, dec, 'Y', img_size)
            if image:
                j, _, _, _ = create_image(image[0], img_idx, instrument + ' Y', get_year_obs(image[0], date_obs_key, date_pattern), ora, odec, op5)
                i += j
            img_idx += 1

            create_color_image(r, g, b, img_idx, 'z-r-g', x, y)
            img_idx += 1

            if i == 0:
                img_idx -= cols

        # NEOWISE
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
                        if (year == prev_year):
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
                    except:
                        print('A problem occurred while creating WISE time series for object ra={ra}, dec={dec}, epoch={epoch}'.format(
                            ra=ra, dec=dec, epoch=i))
                        print(traceback.format_exc())

                hduW1 = fits.PrimaryHDU(data=dataW1/j, header=header)
                hduW2 = fits.PrimaryHDU(data=dataW2/j, header=header)
                images.append((hduW1, hduW2, prev_year))

                for image in images:
                    _, w1, x, y = create_image(image[0], img_idx, 'W1', image[2])
                    _, w2, x, y = create_image(image[1], img_idx, 'W2', image[2])
                    if w1 is not None and w2 is not None:
                        create_color_image(w1, (w1+w2)/2, w2, img_idx, 'W1-W2', x, y, True, image[2])
                        img_idx += 1

        # Save and open the PDF file
        filename = create_obj_name(ra, dec) + '.' + file_format
        plt.savefig(filename, dpi=600, bbox_inches='tight', format=file_format)
        plt.close()

        if open_pdf:
            open_file(filename)

    # --------------------------------------
    # Code for create_finder_charts function
    # --------------------------------------
    warnings.simplefilter('ignore', category=AstropyWarning)
    os.chdir(directory)

    rows = 10
    cols = 6

    if np.isscalar(ra) and np.isscalar(dec):
        finder_charts(ra, dec)
    else:
        if not isinstance(ra, (collections.Sequence, np.ndarray)):
            raise Exception('Ra must be either a sequence or a numpy array')

        if not isinstance(dec, (collections.Sequence, np.ndarray)):
            raise Exception('Dec must be either a sequence or a numpy array')

        if (len(ra) != len(dec)):
            raise Exception('Ra and Dec must have the same length: len(ra)={ra}, len(dec)={dec}'.format(ra=len(ra), dec=len(dec)))

        for i in range(len(ra)):
            try:
                finder_charts(ra[i], dec[i])
            except:
                print('A problem occurred while creating finder charts for object ra={ra}, dec={dec}'.format(ra=ra[i], dec=dec[i]))
                print(traceback.format_exc())
