import warnings
from typing import Dict, List, Tuple
import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.utils.tap.core import TapPlus
from dateutil.parser import parse
from tqdm import tqdm


def create_gaia_finder_chart(ra, dec, size=100, band='G', epoch=2016.0, pmra=0, pmdec=0, plx=0, date=None):

    FIELD_OF_VIEW = dict()
    FIELD_OF_VIEW['G'] = np.array([size, size]) * uu.arcsec

    SCALE_FACTOR = dict()
    SCALE_FACTOR['G'] = np.array([1.0, 1.0])

    RADIUS = dict()
    RADIUS['G'] = np.max(SCALE_FACTOR['G'] * FIELD_OF_VIEW['G'] / np.sqrt(2))

    PIXEL_SCALE = dict()
    PIXEL_SCALE['G'] = 0.1514 * uu.arcsec / uu.pixel

    TRANSFORM_ROTATE = dict()
    TRANSFORM_ROTATE['G'] = 180 * uu.deg + 90 * uu.deg - 37.4 * uu.deg

    FWHM = dict()
    FWHM['G'] = 2.0 * uu.arcsec

    SIGMA_LIMIT = dict()
    SIGMA_LIMIT['G'] = 25

    MAX_PM = 11 * uu.arcsec / uu.year

    FLIP_X = dict()
    FLIP_X['G'] = True

    FLIP_Y = dict()
    FLIP_Y['G'] = False

    GAIA_EPOCH = 2016.0
    GAIA_URL = 'https://gea.esac.esa.int/tap-server/tap'
    GAIA_QUERY = """
    SELECT source_id,ra,dec,parallax,pmra,pmdec,phot_g_mean_mag,phot_rp_mean_mag,
           phot_bp_mean_mag, phot_g_mean_flux, phot_rp_mean_flux, phot_bp_mean_flux
    FROM   gaiadr3.gaia_source
    WHERE  1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius}))
    """

    def propagate_coords(ra, dec, pmra, pmdec, plx, epoch, obs_time: Time) -> Tuple[SkyCoord, Time]:
        # deal with distance
        if plx <= 0:
            distance = None
        else:
            distance = Distance(parallax=plx * uu.mas)
        # need to propagate ra and dec to J2000
        coords = SkyCoord(ra=ra * uu.deg,
                          dec=dec * uu.deg,
                          distance=distance,
                          pm_ra_cosdec=pmra * uu.mas / uu.yr,
                          pm_dec=pmdec * uu.mas / uu.yr,
                          obstime=Time(epoch, format='jd'))
        # work out the delta time between epoch and J2000.0
        with warnings.catch_warnings(record=True) as _:
            jepoch = Time(epoch, format='jd')
            delta_time = (obs_time.jd - jepoch.jd) * uu.day
        # get the coordinates
        with warnings.catch_warnings(record=True) as _:
            # center the image on the current coordinates
            curr_coords = coords.apply_space_motion(dt=delta_time)
            curr_time = obs_time
        return curr_coords, curr_time

    def setup_wcs(image_shape: Tuple[int, int], cent_coords: SkyCoord,
                  pixel_scale: uu.Quantity, rotation: uu.Quantity,
                  flip_x: bool, flip_y: bool) -> WCS:
        # set up a wcs
        naxis2, naxis1 = image_shape
        pix_scale = pixel_scale.to(uu.deg / uu.pixel).value
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [naxis1 / 2, naxis2 / 2]
        if flip_y and flip_x:
            wcs.wcs.cdelt = np.array([pix_scale, -pix_scale])
        elif flip_x:
            wcs.wcs.cdelt = np.array([pix_scale, pix_scale])
        elif flip_y:
            wcs.wcs.cdelt = np.array([pix_scale, -pix_scale])
        else:
            wcs.wcs.cdelt = np.array([-pix_scale, pix_scale])
        wcs.wcs.crval = [cent_coords.ra.deg, cent_coords.dec.deg]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        # add a rotation
        wcs.wcs.pc = np.array([[np.cos(rotation), -np.sin(rotation)],
                               [np.sin(rotation), np.cos(rotation)]])
        return wcs

    def get_gaia_sources(coords: SkyCoord, obstime: Time, radius: uu.Quantity) -> Dict[str, List[float]]:
        # get the gaia time
        gaia_time = Time(GAIA_EPOCH, format='decimalyear')
        # work out the delta time between epoch and J2000.0
        with warnings.catch_warnings(record=True) as _:
            delta_time_now = (obstime.jd - gaia_time.jd) * uu.day
        # modify radius to catch all sources
        search = abs(delta_time_now * MAX_PM).to(uu.deg)
        radius = radius.to(uu.deg) + search
        # define the gaia Tap instance
        gaia = TapPlus(url=GAIA_URL)
        gaia_query = GAIA_QUERY.format(ra=coords.ra.deg, dec=coords.dec.deg,
                                       radius=radius.to(uu.deg).value)
        # launch the query
        job = gaia.launch_job(gaia_query)
        # get the query results
        table = job.get_results()
        # delete job
        del job
        # deal with query being exactly 2000 (the max size)
        if len(table) == 2000:
            print('Too many sources. Launching job asyncronously. Please wait...')
            job = gaia.launch_job_async(gaia_query)
            # get the query results
            table = job.get_results()
            # delete job
            del job
        # storage for sources
        sources = dict()
        sources['gaia_id'] = []
        sources['ra'] = []
        sources['dec'] = []
        sources['G'] = []
        sources['Rp'] = []
        sources['Bp'] = []
        sources['ra_gaia'] = []
        sources['dec_gaia'] = []
        sources['pmra'] = []
        sources['pmdec'] = []
        sources['parallax'] = []
        sources['separation'] = []
        # get parallax mask
        parallax_mask = table['parallax'].mask
        # get all stars in the region where they would be now
        for row in range(len(table)):
            # get table row
            table_row = table[row]
            # deal with distances
            if parallax_mask[row]:
                distance = None
                plx = np.nan
            elif table_row['parallax'] <= 0:
                distance = None
                plx = np.nan
            else:
                plx = table_row['parallax']
                distance = Distance(parallax=plx * uu.mas)
            # need to propagate ra and dec to J2000
            gaia_coords = SkyCoord(ra=table_row['ra'] * uu.deg,
                                   dec=table_row['dec'] * uu.deg,
                                   distance=distance,
                                   pm_ra_cosdec=table_row['pmra'] * uu.mas / uu.yr,
                                   pm_dec=table_row['pmdec'] * uu.mas / uu.yr,
                                   obstime=gaia_time, frame='icrs')
            # center the image on the current coordinates
            with warnings.catch_warnings(record=True) as _:
                curr_coords = gaia_coords.apply_space_motion(dt=delta_time_now)
            # work out the separation between this point and the center point
            separation = coords.separation(curr_coords)
            # add to sources
            sources['gaia_id'].append(table_row['source_id'])
            sources['ra'].append(curr_coords.ra.deg)
            sources['dec'].append(curr_coords.dec.deg)
            sources['G'].append(table_row['phot_g_mean_mag'])
            sources['Rp'].append(table_row['phot_rp_mean_mag'])
            sources['Bp'].append(table_row['phot_bp_mean_mag'])
            sources['ra_gaia'].append(table_row['ra'])
            sources['dec_gaia'].append(table_row['dec'])
            sources['pmra'].append(table_row['pmra'])
            sources['pmdec'].append(table_row['pmdec'])
            sources['parallax'].append(plx)
            sources['separation'].append(separation.value)
        # convert all to numpy arrays
        sources['gaia_id'] = np.array(sources['gaia_id'])
        sources['ra'] = np.array(sources['ra'])
        sources['dec'] = np.array(sources['dec'])
        # warning is for the fact some columns are masked
        with warnings.catch_warnings(record=True) as _:
            sources['G'] = np.array(sources['G'])
            sources['Rp'] = np.array(sources['Rp'])
            sources['Bp'] = np.array(sources['Bp'])
            sources['ra_gaia'] = np.array(sources['ra_gaia'])
            sources['dec_gaia'] = np.array(sources['dec_gaia'])
            sources['pmra'] = np.array(sources['pmra'])
            sources['pmdec'] = np.array(sources['pmdec'])
            sources['parallax'] = np.array(sources['parallax'])
            sources['separation'] = np.array(sources['separation'])
        return sources

    def seed_image(gaia_sources, pixel_scale, obs_coords, fwhm, field_of_view,
                   sigma_limit, band, rotation, flip_x, flip_y, scale_factor):
        # plot out to the scale factor
        field_of_view = field_of_view * scale_factor
        # number of pixels in each direction
        npixel_x = int(field_of_view[0].to(uu.deg) // (pixel_scale * uu.pixel).to(uu.deg))
        npixel_y = int(field_of_view[1].to(uu.deg) // (pixel_scale * uu.pixel).to(uu.deg))
        fwhm_pix = fwhm.to(uu.arcsec) / (pixel_scale * uu.pixel).to(uu.arcsec)
        # create wcs
        wcs = setup_wcs((npixel_y, npixel_x), obs_coords, pixel_scale, rotation, flip_x, flip_y)
        image = np.random.normal(size=(npixel_y, npixel_x), scale=1.0, loc=0)
        nsig_psf = np.array(10 ** ((sigma_limit - gaia_sources[band]) / 2.5))
        # convert all coords to pixel coords
        x_sources, y_sources = wcs.all_world2pix(gaia_sources['ra'], gaia_sources['dec'], 0)
        # create a grid of all pixels
        y, x = np.mgrid[0:npixel_y, 0:npixel_x]
        # get the width of the PSF from the fwhm
        ew = fwhm_pix.value / (2 * np.sqrt(2 * np.log(2)))
        # loop around all sources
        for i in tqdm(range(len(x_sources))):
            if np.isnan(nsig_psf[i]):
                continue
            xdiff0 = x - x_sources[i]
            ydiff0 = y - y_sources[i]
            xdiff = xdiff0 * np.cos(rotation) - ydiff0 * np.sin(rotation)
            ydiff = xdiff0 * np.sin(rotation) + ydiff0 * np.cos(rotation)
            # core of PSF
            exp1 = np.exp(-(xdiff ** 2 + ydiff ** 2) / (2 * ew ** 2))
            image += nsig_psf[i] * exp1.value
            # halo of PSF
            exp_halo = np.exp(-(xdiff ** 2 + ydiff ** 2) / (2 * (ew * 3) ** 2))
            image += nsig_psf[i] * exp_halo.value * 1e-3
            # # spike in y
            # exp_spike_y = np.exp(-np.abs(xdiff ** 2 + 0.1 * np.abs(ydiff) ** 2) / (2 * ew ** 2))
            # image += nsig_psf[i] * exp_spike_y * 1e-2
            # # spike in x
            # exp_spike_x = np.exp(-(0.1 * np.abs(xdiff) ** 2 + np.abs(ydiff) ** 2) / (2 * ew ** 2))
            # image += nsig_psf[i] * exp_spike_x * 1e-2

        # display with an arcsinh stretch
        image = np.arcsinh(image)

        # return the image and the wcs
        return image, wcs

    try:
        epoch = Time(epoch, format='decimalyear')
    except Exception:
        raise ValueError(f'Cannot parse --epoch={epoch}')
    if date:
        try:
            date = Time(parse(str(date)))
        except Exception:
            raise ValueError(f'Cannot parse --date={date}')
    else:
        date = epoch
    obs_coords, obs_time = propagate_coords(ra, dec, pmra, pmdec, plx, epoch, date)
    gaia_sources = get_gaia_sources(obs_coords, obs_time, radius=RADIUS['G'])
    kwargs_gaia = dict(pixel_scale=PIXEL_SCALE['G'], obs_coords=obs_coords,
                       fwhm=FWHM['G'], field_of_view=FIELD_OF_VIEW['G'],
                       sigma_limit=SIGMA_LIMIT['G'], band=band,
                       scale_factor=SCALE_FACTOR['G'])
    print(f'Seeding Gaia {band}-band image')
    image, wcs = seed_image(gaia_sources, flip_x=False, flip_y=False,
                            rotation=0 * uu.deg, **kwargs_gaia)
    return image, wcs, obs_coords
