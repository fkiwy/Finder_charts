from Gaia_finder_chart import create_gaia_finder_chart


objname = "2MASS J21143529-6124082"
ra = 318.6470346
dec = -61.4024763
pmra = 20.256
pmde = -78.537
plx = 14.1171
epoch = 2016.0

create_gaia_finder_chart(objname, epoch, ra, dec, pmra=0, pmde=0, plx=0, apero=False, simbad=False)
