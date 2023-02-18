def get_date_df(arg):
    from parallel_xarray import get_date_files, sum_and_sample, mean_and_sample, stdv_and_sample, sample_zscore, sample_tif
    date = arg[0]
    pts_df = arg[1]
    loc = arg[2]

    # pull lats and lons 
    lats = list(pts_df.lat)
    lons = list(pts_df.lon)
    # drop extra columns 
    pts_df = pts_df.drop(['Unnamed: 0', 'FID', 'geometry', 'long'], axis = 1)
    pts_df['date'] = date

    # # get date files 
    # precip_16, precip_32, precip_48, precip_64, temp_16, temp_32, temp_48, temp_64, zscore, savi = get_date_files(date, loc)
    # get date files 
    precip_16, precip_32, precip_48, precip_64, temp_16, temp_32, temp_48, temp_64, zscore = get_date_files(date, loc)

    # list of precip and temp files 
    precip = [precip_16, precip_32, precip_48, precip_64]
    temp = [temp_16, temp_32, temp_48, temp_64]
    
    
    
    # instantiate lists for storing results 
#     clim_data = []
#     cols = []
    precip_sums = []
    precip_means = []
    precip_stdvs = []
    temp_means = []
    temp_stdvs = []
    # get list of data for each time window, variable, and function 
    for time in precip: 
        argument = (time, lats, lons, "ppt")
        precip_sums.append(sum_and_sample(argument))
        precip_means.append(mean_and_sample(argument))
        precip_stdvs.append(stdv_and_sample(argument))

    for time in temp: 
        argument = (time, lats, lons, "temp")
        temp_means.append(mean_and_sample(argument))
        temp_stdvs.append(stdv_and_sample(argument))
        
    # get z score 
    SAVI_zscore = sample_zscore((zscore, lats, lons))
    pts_df['SAVI_zscore'] = SAVI_zscore[1]

    # # get savi measurement 
    # savi = sample_tif((savi, lats, lons))
    # pts_df['SAVI'] = savi[1]

    # add to dataframe 
    for data in precip_sums: 
        # add a new column to the dataframe for each data retrieval
        col = data[0]
        pts_df[col] = data[1]
    # add to dataframe 
    for data in precip_means: 
        # add a new column to the dataframe for each data retrieval
        col = data[0]
        pts_df[col] = data[1]
    # add to dataframe 
    for data in precip_stdvs: 
        # add a new column to the dataframe for each data retrieval
        col = data[0]
        pts_df[col] = data[1]
        # add to dataframe 
    for data in temp_means: 
        # add a new column to the dataframe for each data retrieval
        col = data[0]
        pts_df[col] = data[1]
    # add to dataframe 
    for data in temp_stdvs: 
        # add a new column to the dataframe for each data retrieval
        col = data[0]
        pts_df[col] = data[1]
    return(pts_df)


def get_SAVI_date_df(arg):
    from parallel_xarray import sample_tif
    import os
    date = arg[0]
    pts_df = arg[1]
    loc = arg[2]

    # pull lats and lons 
    lats = list(pts_df.lat)
    lons = list(pts_df.lon)
    # drop extra columns 
    pts_df = pts_df.drop(['Unnamed: 0', 'FID', 'geometry', 'long'], axis = 1)
    pts_df['date'] = date

    # get SAVI file 
    if loc == "Narok":
        og_path = "E:\\bulk_download_USGS\\Bulk_Order_Maasai_Mara\\Landsat_8-9_OLI_TIRS_C2_L2\\SAVI"
        og_files = os.listdir(og_path)
        savi_og = [os.path.join(og_path, file) for file in og_files]
    else: 
        og_path = "E:\\bulk_download_USGS\\Bulk_Order_Turkana\\Landsat_8-9_OLI_TIRS_C2_L2\\SAVI_mosaic"
        og_files = os.listdir(og_path)
        savi_og = [os.path.join(og_path, file) for file in og_files]

    # pull the proper SAVI file 
    savi = [file for file in savi_og if str(date)[:10].replace("-", "") in file][0]

    # get savi measurement 
    savi = sample_tif((savi, lats, lons))
    pts_df['SAVI'] = savi[1]

    return(pts_df)