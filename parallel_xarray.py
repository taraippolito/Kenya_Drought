def make_8_day_avg(raster_list): 
    import os
    import rioxarray as rxr

    #open all of the files with xarray
    opened = [rxr.open_rasterio(file) for file in raster_list]
    #sub nan for nodata so that the calculations don't mess up the nodata values
    with_nan = [file.where(file.values != file.rio.nodata) for file in opened]
    #take the average over all rasters
    numerator = sum(opened)
    denominator = len(raster_list)
    avg = numerator / denominator
    
    # update attributes
    avg.rio.write_crs(opened[0].rio.crs, inplace=True)
    avg.rio.update_attrs(opened[0].attrs, inplace=True)
    avg.rio.update_encoding(opened[0].encoding, inplace=True)

    #create a new file 
    #output file name
    avg_tif_name = raster_list[0][10:-4] + "-" + raster_list[-1][15:-4] + "avg.TIF"
    #output path 
    output_path = os.path.join("F:\\CHIRTS", avg_tif_name)
    # Export data to geotiff
    avg.rio.to_raster(output_path)
    
def tmaxnc_to_tif(pair): 
    import rioxarray as rxr
    import xarray as xr
    
    file = pair[0]
    time_step = pair[1]
    #open the file
    nc_file = xr.open_dataset(file)
    # create output string
    output = file[:-7] + str(time_step)[:10] + ".tif"
    # select the time step 
    nc_file = nc_file.Tmax.sel(time = time_step)
    # set crs
    nc_file = nc_file.rio.set_crs("epsg:4326")
    # save tif 
    nc_file.rio.to_raster(output)
    
def tminnc_to_tif(pair): 
    import rioxarray as rxr
    import xarray as xr
    
    file = pair[0]
    time_step = pair[1]
    #open the file
    nc_file = xr.open_dataset(file)
    # create output string
    output = file[:-7] + str(time_step)[:10] + ".tif"
    # select the time step 
    nc_file = nc_file.Tmin.sel(time = time_step)
    # set crs
    nc_file = nc_file.rio.set_crs("epsg:4326")
    # save tif 
    nc_file.rio.to_raster(output)
    
def sample_tif(args):
    import xarray as xr
    import rioxarray as rxr
    tif_path = args[0]
    lat = args[1]
    lon = args[2]
    #open tif
    ds = xr.open_dataarray(tif_path).rio.reproject("EPSG:4326")
    #create indexing data arrays
    lats = xr.DataArray(lat, dims='z')
    lons = xr.DataArray(lon, dims='z')
    # pull data 
    data = ds.sel(x = lons, y = lats, method = 'nearest')
    #turn into a list 
    data_list = list(data.values[0])
    # return data list 
    return (tif_path, data_list)


def sample_hdf(args):
    import xarray as xr
    import rioxarray as rxr
    tif_path = args[0]
    lat = args[1]
    lon = args[2]
    #open tif
    ds = rxr.open_rasterio(tif_path).rio.reproject("EPSG:4326")
    #create indexing data arrays
    lats = xr.DataArray(lat, dims='z')
    lons = xr.DataArray(lon, dims='z')
    # pull data 
    data = ds.sel(x = lons, y = lats, method = 'nearest')
    #turn into a list 
    data_list = list(data["Npp_500m"].values[0])
    # return data list 
    return (tif_path, data_list)

    # function to turn nc files into tif for GSMap_Gauge precip data 
def pptnc_to_tif(pair):
    import rioxarray as rxr
    import xarray as xr
    
    file = pair[0]
    time_step = pair[1]
    #open the file
    nc_file = xr.open_dataset(file)
    # set spatial dimensions 
    nc_file.precip.rio.set_spatial_dims("lon", "lat", inplace=True)

    # create output string
    output =  "E:\\GSMap\\Kenya_ppt\\" + str(time_step)[:10] + ".tif"
    # select the time step & set spatial dimensions
    nc_file = nc_file.precip.sel(time = time_step).rio.set_spatial_dims("lon", "lat", inplace=True)
    # set crs
    nc_file = nc_file.rio.set_crs("epsg:4326")

    # save tif 
    nc_file.rio.to_raster(output)


def make_8_day_avg_GSM(raster_list): 
    import rioxarray as rxr
    import os 
    #open all of the files with xarray
    opened = [rxr.open_rasterio(file) for file in raster_list]
    #sub nan for nodata so that the calculations don't mess up the nodata values
    with_nan = [file.where(file.values != file.rio.nodata) for file in opened]
    #take the average over all rasters
    numerator = sum(with_nan)
    denominator = len(raster_list)
    avg = numerator / denominator
    
    # update attributes
    avg.rio.write_crs(opened[0].rio.crs, inplace=True)
    avg.rio.update_attrs(opened[0].attrs, inplace=True)
    avg.rio.update_encoding(opened[0].encoding, inplace=True)

    #create a new file 
    #output file name
    avg_tif_name = raster_list[0][-14:-4] + "_" + raster_list[-1][-14:-4] + "_avg.TIF"
    #output path 
    output_path = os.path.join("F:\\GSMap\\Kenya_ppt", avg_tif_name)
    # Export data to geotiff
    avg.rio.to_raster(output_path)

# get veg index z scores and turn into new raster for each date
# input list of files to get z score of (should be the same tile and dates within 20 day window) and the arguments for file output
# input list of files to get z score of (should be the same tile and dates within 20 day window) and the arguments for file output
def get_VI_zscore(arg):
    # pull from the argument passed
    meta_files = arg[0]
    # index_list = arg[0]
    kwargs = arg[1]
    tile = arg[2]
    date = arg[3]
    
    # import
    import rasterio as rio
    import datetime
    import os 
    import rioxarray as rxr
    import numpy as np
    
    # get the files needed which are within 20 days of date, reproject to EPSG4326 
    zscore_files = [f.rio.reproject("EPSG:4326") for f in meta_files]
    
    # pull the main file to calculate z score of 
    rxr_MAIN = zscore_files[0]
    # get the numpy array of the main file 
    MAIN = rxr_MAIN.to_numpy()[0]
    
    # make all files the same extent, pull the numpy array from the rxr file
    repr_files = [f.rio.reproject_match(rxr_MAIN) for f in zscore_files] # this works
    matched = []
    for f in repr_files: 
        matched.append(f.assign_coords({"x": rxr_MAIN.x, "y": rxr_MAIN.y})) 
    
    opened = []
    # clean up the nodata values
    for f in matched:
        opened.append(f.where(f != -255))
    print ("files opened and cleaned")
    
    #take the average over all rasters
    numerator = sum(opened)
    
#     plt.imshow(numerator)
    denominator = len(opened)
    mean = (numerator / denominator)
    
    sqd_err = [(f-mean)**2 for f in opened]
    sum_sqd_err = sum(sqd_err)
    std = np.sqrt(sum_sqd_err / len(opened))
    
    
    # z score calculation
    z_score = ((MAIN - mean) / std)
    print ("z score calculated")
    
#     return (zscore_files, opened, numerator, z_score)
#     write SAVI raster file
    new_path = "E:\\bulk_download_USGS\\Bulk_Order_Turkana\\Landsat_8-9_OLI_TIRS_C2_L2\\SAVI_zscore_new\\" + tile + "_" + date + "_zscore.TIF"
    # update attributes
    z_score.rio.write_crs(opened[0].rio.crs, inplace=True)
    z_score.rio.update_attrs(opened[0].attrs, inplace=True)
    z_score.rio.update_encoding(opened[0].encoding, inplace=True)
    z_score.rio.to_raster(new_path)
    # Image = rasterio.open(new_path, 'w', **kwargs)

    # Image.write_band(1, z_score)
    # Image.close()


    # input a date, get all date files for precip and temp 
# input a date, get all date files for precip and temp 
def get_date_files(savi_date, loc): 
    import datetime
    import os
    
    # get all dates 16 days prior to date 
    window_16 = [savi_date - datetime.timedelta(days=x) for x in range(16)]
    str_16 = [str(d)[:10] for d in window_16]
    # get all dates 32 days prior to date
    window_32 = [savi_date - datetime.timedelta(days=x) for x in range(32)]
    str_32 = [str(d)[:10] for d in window_32]
    # get all dates 48 days prior to date 
    window_48 = [savi_date - datetime.timedelta(days=x) for x in range(48)]
    str_48 = [str(d)[:10] for d in window_48]
    # get all dates 64 days prior to date 
    window_64 = [savi_date - datetime.timedelta(days=x) for x in range(64)]
    str_64 = [str(d)[:10] for d in window_64]
    
    # get precip files 
    gsmap_path = "E:\\GSMap\\Kenya_ppt"
    gsmap_files = os.listdir(gsmap_path)
    ppt_files = [os.path.join(gsmap_path, file) for file in gsmap_files if ".tif" in file]
    # get temp files
    chirts_path = "E:\\NOAA_temp\\Tavg"
    chirts_files = os.listdir(chirts_path)
    temp_files = [os.path.join(chirts_path, file) for file in chirts_files if "tavg.TIF" in file]
    # get z score file 
    if loc == "Narok":
        vi_path = "E:\\bulk_download_USGS\\Bulk_Order_Maasai_Mara\\Landsat_8-9_OLI_TIRS_C2_L2\\SAVI_zscore_new"
        vi_files = os.listdir(vi_path)
        savi_files = [os.path.join(vi_path, file) for file in vi_files]
        # og_path = "F:\\bulk_download_USGS\\Bulk_Order_Maasai_Mara\\Landsat_8-9_OLI_TIRS_C2_L2\\SAVI"
        # og_files = os.listdir(og_path)
        # savi_og = [os.path.join(og_path, file) for file in og_files]
    else: 
        vi_path = "E:\\bulk_download_USGS\\Bulk_Order_Turkana\\Landsat_8-9_OLI_TIRS_C2_L2\\SAVI_zscore_mosaic"
        vi_files = os.listdir(vi_path)
        savi_files = [os.path.join(vi_path, file) for file in vi_files]
        # og_path = "F:\\bulk_download_USGS\\Bulk_Order_Turkana\\Landsat_8-9_OLI_TIRS_C2_L2\\SAVI_mosaic"
        # og_files = os.listdir(og_path)
        # savi_og = [os.path.join(og_path, file) for file in og_files]
    
    
    # get all precip files 
    precip_16 = [file for file in ppt_files if any(substring in file for substring in str_16)]
    precip_32 = [file for file in ppt_files if any(substring in file for substring in str_32)]
    precip_48 = [file for file in ppt_files if any(substring in file for substring in str_48)]
    precip_64 = [file for file in ppt_files if any(substring in file for substring in str_64)]
    
    # get all temp files 
    temp_16 = [file for file in temp_files if any(substring in file for substring in str_16)]
    temp_32 = [file for file in temp_files if any(substring in file for substring in str_32)]
    temp_48 = [file for file in temp_files if any(substring in file for substring in str_48)]
    temp_64 = [file for file in temp_files if any(substring in file for substring in str_64)]
    
    # get z score file 
    zscore = [file for file in savi_files if str(savi_date)[:10].replace("-", "") in file][0]
    # savi = [file for file in savi_og if str(savi_date)[:10].replace("-", "") in file][0]
    
    # return (precip_16, precip_32, precip_48, precip_64, temp_16, temp_32, temp_48, temp_64, zscore, savi)
    return (precip_16, precip_32, precip_48, precip_64, temp_16, temp_32, temp_48, temp_64, zscore)

# get sum of rasters and sample 
def sum_and_sample(args): 
    import xarray as xr
    import rioxarray as rxr
    
    tif_paths = args[0]
    lat = args[1]
    lon = args[2]
    variable = args[3]
    length = 0
    if len(tif_paths) == 16 or len(tif_paths) == 15: 
        length = 16
    if len(tif_paths) == 32 or len(tif_paths) == 31: 
        length = 32
    if len(tif_paths) == 48 or len(tif_paths) == 47: 
        length = 48
    if len(tif_paths) == 64 or len(tif_paths) == 63: 
        length = 64
    #open tifs
    ds = [xr.open_dataarray(tif_path).rio.reproject("EPSG:4326") for tif_path in tif_paths]
    # sum tifs 
    tif_sum = sum(ds)
    #create indexing data arrays
    lats = xr.DataArray(lat, dims='z')
    lons = xr.DataArray(lon, dims='z')
    # pull data 
    data = tif_sum.sel(x = lons, y = lats, method = 'nearest')
    #turn into a list 
    data_list = list(data.values[0])

    # return data list 
    return ((str(length) + "_day_sum_" + variable), data_list)

# get mean of rasters and sample - arguments = paths of tifs to calculate and sample, lat, lon, variable string
def mean_and_sample(args): 
    import xarray as xr
    import rioxarray as rxr
    
    tif_paths = args[0]
    lat = args[1]
    lon = args[2]
    variable = args[3]
    length = 0
    if len(tif_paths) == 16 or len(tif_paths) == 15: 
        length = 16
    if len(tif_paths) == 32 or len(tif_paths) == 31: 
        length = 32
    if len(tif_paths) == 48 or len(tif_paths) == 47: 
        length = 48
    if len(tif_paths) == 64 or len(tif_paths) == 63: 
        length = 64
    
    #open tifs
    ds = [xr.open_dataarray(tif_path).rio.reproject("EPSG:4326") for tif_path in tif_paths]
    # sum tifs 
    tif_mean = sum(ds)/len(tif_paths)
    #create indexing data arrays
    lats = xr.DataArray(lat, dims='z')
    lons = xr.DataArray(lon, dims='z')
    # pull data 
    data = tif_mean.sel(x = lons, y = lats, method = 'nearest')
    #turn into a list 
    data_list = list(data.values[0])
    # return data list 
    return ((str(length) + "_day_mean_" + variable), data_list)

# get stdv of rasters and sample - arguments = paths of tifs to calculate and sample, lat, lon, variable string
def stdv_and_sample(args): 
    import xarray as xr
    import rioxarray as rxr
    import numpy as np
    
    tif_paths = args[0]
    lat = args[1]
    lon = args[2]
    variable = args[3]
    length = 0
    if len(tif_paths) == 16 or len(tif_paths) == 15: 
        length = 16
    if len(tif_paths) == 32 or len(tif_paths) == 31: 
        length = 32
    if len(tif_paths) == 48 or len(tif_paths) == 47: 
        length = 48
    if len(tif_paths) == 64 or len(tif_paths) == 63: 
        length = 64
    
    #open tifs
    ds = [xr.open_dataarray(tif_path).rio.reproject("EPSG:4326") for tif_path in tif_paths]
    # sum tifs 
    tif_mean = sum(ds)/len(tif_paths)
    sq_dist = [(tif - tif_mean)**2 for tif in ds]
    sum_sq_dist = sum(sq_dist)
    tif_stdv = np.sqrt(sum_sq_dist/len(tif_paths))
    
    #create indexing data arrays
    lats = xr.DataArray(lat, dims='z')
    lons = xr.DataArray(lon, dims='z')
    # pull data 
    data = tif_stdv.sel(x = lons, y = lats, method = 'nearest')
    #turn into a list 
    data_list = list(data.values[0])
    # return data list 
    return ((str(length) + "_day_stdv_" + variable), data_list)

def sample_zscore(args):
    import xarray as xr
    import rioxarray as rxr
    import numpy as np
    
    tif_path = args[0]
    lat = args[1]
    lon = args[2]
    
    # open z score tif 
    ds = xr.open_dataarray(tif_path).rio.reproject("EPSG:4326")

    #create indexing data arrays
    lats = xr.DataArray(lat, dims='z')
    lons = xr.DataArray(lon, dims='z')
    # pull data 
    data = ds.sel(x = lons, y = lats, method = 'nearest')
    #turn into a list 
    data_list = list(data.values[0])
    # return data list 
    return ("SAVI_zscore", data_list)