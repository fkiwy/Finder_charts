from Finder_charts_VVV import create_finder_charts
from astropy.table import Table
import pandas as pd
import numpy as np
import tempfile


ra = 268.9187535  # additional NEOWISE epochs
dec = 70.8374857

ra = 133.787
dec = -7.24515

ra = 23.8559091
dec = 2.0877172

ra = 46.61254
dec = 1.9625

ra = 152.92262
dec = 28.76333

ra = 249.9235534
dec = -68.79868325

ra = 149.8377146
dec = 2.859148817

ra = 304.6036807
dec = -74.39161523

ra = 264.565964
dec = -27.955309

ra = 47.3320966
dec = -50.2706362

ra = 249.9235534
dec = -68.79868325

ra = 134.0132646  # poor neowise contrast, has all 5 DECam bands !!!!!
dec = -7.3931418

ra = 75.7761571
dec = -56.808827

ra = 1.454075
dec = -21.9557694

ra = 16.9701886
dec = 0.6992208

ra = 31.090342
dec = -2.7390737

ra = 224.3174306
dec = -21.37108003

ra = 221.0122447
dec = 8.313252505

ra = 12.36883599
dec = 4.68341893

ra = 209.2891781
dec = 55.7474398

ra = 221.5031415
dec = 0.41418472

ra = 43.53973113
dec = 2.399550061

ra = 126.3301355
dec = 21.2634308

ra = 152.7709667
dec = 38.4372530

ra = [134.0132646, 75.7761571, 1.454075, 16.9701886, 31.090342, 221.0122447, 12.36883599, 209.2891781, 221.5031415, 43.53973113]
dec = [-7.3931418, -56.808827, -21.9557694, 0.6992208, -2.7390737, 8.313252505, 4.68341893, 55.7474398, 0.41418472, 2.399550061]

# Case 1: ras and decs are each in a separate tuple
ra = (126.3301355, 152.7709667)
dec = (21.2634308, 38.4372530)

# Case 2: ras and decs are each in a separate list
ra = [126.3301355, 152.7709667]
dec = [21.2634308, 38.4372530]

# Case 3: ras and decs are each in a separate array
ra = np.array([126.3301355, 152.7709667])
dec = np.array([21.2634308, 38.4372530])

# Case 4: ra and dec are in the same tuple/list within a tuple/list
coords = [(126.3301355, 21.2634308), (152.7709667, 38.4372530)]
coords = ((126.3301355, 21.2634308), (152.7709667, 38.4372530))
coords = [[126.3301355, 21.2634308], [152.7709667, 38.4372530]]
coords = ([126.3301355, 21.2634308], [152.7709667, 38.4372530])
ra, dec = np.split(np.array(coords), 2, axis=1)
ra = ra.flatten()
dec = dec.flatten()

# Case 5: ras and decs are each in a separate column of a pandas data frame
coords = pd.DataFrame(coords, columns=['ra', 'dec'])
ra = coords['ra'].values
dec = coords['dec'].values

# Case 6: ras and decs are each in a separate column of an astropy table
coords = Table([ra, dec], names=('ra', 'dec'))
ra = np.array(coords['ra'])
dec = np.array(coords['dec'])


# Fits file has no image data for epoch 2013
ra = 337.4922649
dec = 30.40228243

# Issues
ra = 26.97838266
dec = 23.6616913
ra = 102.471579
dec = 17.37210717
ra = 280.2551196
dec = 61.34251259

"""
# Comovers
from astropy.io import ascii
test_set = ascii.read('C:/Users/wcq637/Documents/Private/BYW/Comovers/test/comovers.csv', format='csv')
test_set = test_set[test_set['component'] == 'A']
# test_set.pprint_all()
ra = np.array(test_set['nsc_ra'])
dec = np.array(test_set['nsc_dec'])
"""

ra = 8.335532
dec = -42.457337

ra = 42.915712
dec = 5.189519

ra = 356.182971
dec = 0.520115

ra = 9.973995
dec = -45.396616

ra = 10.244511
dec = -29.703435


ra = [134.0132646, 75.7761571, 1.454075, 16.9701886, 31.090342, 221.0122447, 12.36883599, 209.2891781, 221.5031415, 43.53973113]
dec = [-7.3931418, -56.808827, -21.9557694, 0.6992208, -2.7390737, 8.313252505, 4.68341893, 55.7474398, 0.41418472, 2.399550061]

"""
ra = 16.9701886
dec = 0.6992208

ra = 209.2891781
dec = 55.7474398
"""

ra = 75.7761571
dec = -56.808827

# Viking
ra = 26.9880948
dec = -35.8503841

# VVV
ra = 269.9346478
dec = -22.0947810

# Viking
ra = 221.5031415
dec = 0.41418472

# Y1 (JWST)
ra = 83.8202991
dec = -75.0068543

# Latent
ra = 22.3554008
dec = 7.3135396

# WiseView issue (case reported by Tom)
ra = 127.5804541
dec = 89.4409387

# ???
ra = 164.401
dec = 85.4526

ra = 16.9701886
dec = 0.6992208

ra = 22.2873
dec = 6.79308

# Tom's ring
ra = 66.7742411
dec = 34.251332

# UKIDSS
ra = 36.37572605
dec = 3.634077368

# VVV
ra = 211.773493
dec = -62.782159

ra = 265.094702
dec = -37.144588


create_finder_charts(ra, dec, img_size=10, overlays=False, overlay_color='red', dss=False, twomass=False, spitzer=False, wise=False,
                     ukidss=False, vhs=False, vvv=True, viking=False, ps1=False, decam=False, neowise=False, neowise_contrast=10, chrono_order=True,
                     object_info=True, directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, open_pdf=None, open_file=True,
                     file_format='pdf', save_result_tables=False, result_tables_format='ipac', result_tables_extension='dat')
