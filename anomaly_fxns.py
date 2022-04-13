#dictionary of dataframe name and dataframe - FULL DATA
def get_dataframes(df, date):
    date_df = df[df.Date == date]
    return (date_df)

def get_ppt_stats(row):
    ppt_sum_stats = []
    
    row_date = row['Date']
    # get string of date to match ppt columns
    ppt_date = "F" + row_date[1:5] + "_"+ row_date[5:7]+"_"+row_date[7:]+ "ppt"
    # get ppt columns list 
    ppt_col_list = [col for col in row.index if "ppt" in col]
    try: 
        # get index of date col among all ppt columns
        index_of_date = ppt_col_list.index(ppt_date) 
    except: 
        index_of_date = 0
    
    if index_of_date - 16 >= 0: 
        # get the 16 columns before that date
        cols_b4 = ppt_col_list[(index_of_date - 16): index_of_date]
        # get the cumulative ppt over the past 16 days 
        sum_ppt = row[cols_b4].sum()
        # get the mean ppt over the past 16 days 
        mean_ppt = row[cols_b4].mean()
        # get the max ppt over the past 16 days 
        max_ppt = row[cols_b4].max()
        # get the min ppt over the past 16 days 
        min_ppt = row[cols_b4].min()
        # get the median ppt over the past 16 days 
        median_ppt = row[cols_b4].median()
        # get the median ppt over the past 16 days 
        std_ppt = row[cols_b4].std()
        # number of days it rained 
        days_of_rain = row[cols_b4].astype(bool).sum()
        #upper and lower quantiles 
#         pptq1 = row[cols_b4].quantile(q=.25)
#         pptq4 = row[cols_b4].quantile(q=.75)


        row["sum_ppt"] = sum_ppt
        row["mean_ppt"] = mean_ppt
        row["max_ppt"] = max_ppt
        row["min_ppt"] = min_ppt
        row["median_ppt"] = median_ppt
        row["std_ppt"] = std_ppt
        row["days_of_rain"] = days_of_rain
#         row["pptq1"] = pptq1
#         row["pptq4"] = pptq4

    else: 
        row["sum_ppt"] = np.NaN
        row["mean_ppt"] = np.NaN
        row["max_ppt"] = np.NaN
        row["min_ppt"] = np.NaN
        row["median_ppt"] = np.NaN
        row["std_ppt"] = np.NaN
        row["days_of_rain"] = np.NaN
#         row["pptq1"] = np.NaN
#         row["pptq4"] = np.NaN
    return (row)

def make_merge_dfs(df, date): 
    #match ppt column
    ppt_date = "F" + date[1:5] + "_"+ date[5:7]+"_"+date[7:]+ "ppt"
    ppt_col_list = [col for col in df.columns if "ppt" in col]
    try: 
        # get index of date col among all ppt columns
        index_of_date = ppt_col_list.index(ppt_date) 
    except: 
        index_of_date = 0
    #check if there are a full 16 days before the date
    if index_of_date - 16 >= 0: 
        cols_b4 = ppt_col_list[(index_of_date - 16): index_of_date]
    else: 
        cols_b4 = []
    
    # pull the 16 ppt columns before the date, and the soil columns
    merge_cols = cols_b4 + ["Unnamed: 0", 'sand30',
                 'sand60','sand100','silt30','silt60','silt100','clay30','clay60','clay100',
                 'soc30','soc60','soc100','bdod30','bdod60','bdod100','cfvo30','cfvo60','cfvo100']
    merge_df = df[merge_cols]
    return (merge_df)


def merge_and_calc_ppt(df_anom, df_ppt):
	# FIRST MERGE DATAFRAMES 
	new_df = df_anom.merge(df_ppt, how = "left", left_on = "Unnamed: 0", right_on = "Unnamed: 0")

	# SECOND CALCULATE PPT STATS
	ppt_cols = [col for col in new_df.columns if "ppt" in col]
	# get the cumulative ppt over the past 16 days 
	new_df["sum_ppt"] = new_df[ppt_cols].sum(axis = 1)
	# get the mean ppt over the past 16 days 
	new_df["mean_ppt"] = new_df[ppt_cols].mean(axis = 1)
	# get the max ppt over the past 16 days 
	new_df["max_ppt"] = new_df[ppt_cols].max(axis = 1)
	# get the min ppt over the past 16 days 
	new_df["min_ppt"] = new_df[ppt_cols].min(axis = 1)
	# get the median ppt over the past 16 days 
	new_df["median_ppt"] = new_df[ppt_cols].median(axis = 1)
	# get the median ppt over the past 16 days 
	new_df["std_ppt"] = new_df[ppt_cols].std(axis = 1)
	# number of days it rained 
	new_df["days_of_rain"] = new_df[ppt_cols].astype(bool).sum(axis = 1)

	# THIRD CALCULATE ANOMALY 
	new_df['NDVI_anomaly'] = (new_df['NDVI_value'] - new_df['mean_NDVI'])/new_df['stdv_NDVI']

	DONE = new_df.drop(ppt_cols, axis = 1)

	return (DONE)


