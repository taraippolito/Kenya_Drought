# -*- coding: utf-8 -*-
"""
Generated by ArcGIS ModelBuilder on : 2021-11-22 10:08:10
"""
import arcpy

def Model5():  # Model 5

    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = False

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("spatial")

    random_samples_2_ = "random_samples"
    LC08_L2SP_170058_20140201_20200912_02_T1_NDVI_tiff = arcpy.Raster("C:\\Tara_Fall_2019\\Kenya_Drought\\Calculated_NDVI\\LC08_L2SP_170058_20140201_20200912_02_T1_NDVI.tiff")
    LC08_L2SP_170058_20140201_20200912_02_T1_NDVICLIP_tif = arcpy.Raster("C:\\Tara_Fall_2019\\Kenya_Drought\\Calculated_NDVI\\LC08_L2SP_170058_20140201_20200912_02_T1_NDVICLIP.tif")

    # Process: Extract Multi Values to Points (Extract Multi Values to Points) (sa)
    random_samples = arcpy.sa.ExtractMultiValuesToPoints(in_point_features=random_samples_2_, in_rasters=[[LC08_L2SP_170058_20140201_20200912_02_T1_NDVI_tiff, "LC08_L2SP_"], [LC08_L2SP_170058_20140201_20200912_02_T1_NDVICLIP_tif, "LC08_L2SP1"]], bilinear_interpolate_values="NONE")
    .save(Extract_Multi_Values_to_Points)


if __name__ == '__main__':
    # Global Environment settings
    with arcpy.EnvManager(scratchWorkspace=r"C:\Users\Researcher\Documents\ArcGIS\Projects\Kenya_drought\Kenya_drought.gdb", workspace=r"C:\Users\Researcher\Documents\ArcGIS\Projects\Kenya_drought\Kenya_drought.gdb"):
        Model5()
