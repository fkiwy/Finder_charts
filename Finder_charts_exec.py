from Finder_charts import create_finder_charts

ra = 16.9701886
dec = 0.6992208

create_finder_charts(ra, dec, img_size=100, overlays=True, overlay_color='red', dss=True, twomass=True, spitzer=True, wise=True,
                     ukidss=True, uhs=True, vhs=True, vvv=True, viking=True, ps1=True, decam=True, neowise=True, neowise_contrast=10,
                     chrono_order=True, object_info=True, directory='Finder_charts_output', cache=False, show_progress=True, timeout=300,
                     open_file=False, file_format='pdf', save_result_tables=True, result_tables_format='ipac', result_tables_extension='dat')
