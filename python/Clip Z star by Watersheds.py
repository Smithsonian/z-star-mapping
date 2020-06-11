# This code was written by James Holmquist on August 19, 2018
# The goal is to take maps of relative vertical vulnerability of wetlands at the national scale as inputs
# ... and clip them at the water shed level. And output the associated raster data table as the

# I repurpoused this script from the marsh meta-analysis for the Z* mapping project on 4 April 2020

# Import all necessary libraries
import arcpy
import os
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# establish scratch workspace
scratch_path = "F:/temp/ArcScratch"
arcpy.env.scratchWorkspace = scratch_path
lyr= "lyr"
poly = "in_memory/temp"

# First establish a directory where all of the large geographic data files are held
geo_directory = ""

# Huc8 input files
huc8 = geo_directory + "F:/z-star-spatial-data/derivative-maps/zStar/HUC8_Tidal_170404.shp"

# Create a list of input file names
input_names = ["F:/z-star-spatial-data/derivative-maps/zStar/zStarUncertainty_30m_maskedByCcapEEM2010_20427.tif",
               "F:/z-star-spatial-data/derivative-maps/zStar/zStar_30m_maskedByCcapEEM2010_200707.tif",
               "F:/z-star-spatial-data/input-layers/datum-layers/MHW_MSL_30m_EEM_200604.img",
               "F:/z-star-spatial-data/input-layers/datum-layers/RSLR_1983to2001_30m_EEM_200604.tif"
               ]

# Create a list of output files names
output_filepaths = ["F:/z-star-spatial-data/temp/z-star-uncertainty",
                    "F:/z-star-spatial-data/temp/z-star",
                    "F:/z-star-spatial-data/temp/MHW_MSL",
                    "F:/z-star-spatial-data/temp/RSLR"
                    ]

# Iterate through the HUC8
# set up export tags
names = arcpy.da.SearchCursor(huc8, ["Abbrev"])
name_vect=[]
for row in names:
    string1 = str(row[0])
    name_vect.append(string1)

# Iterate through the input file names
failures = []
h = 0
for input_name in input_names:

    print(input_name)

    # loop through HUC8s
    for i in range(len(name_vect)): # For every row in the polygon

        # Create output file name
        tif_export = output_filepaths[h] + "/" + name_vect[i] + ".tif"
        # dem_temp1 = (scratch_path + "temp1.tif")
        print("  " + tif_export)

        # if a huc overlaps than work with it
        if os.path.isfile(tif_export):
            print("  " + tif_export + " already done.")
        else:
            # Make a layer from the feature class
            arcpy.MakeFeatureLayer_management(huc8, lyr)
            fid_which = ' "FID" = ' + str(i) + ' '

            # select the row
            arcpy.SelectLayerByAttribute_management(lyr, selection_type="NEW_SELECTION", where_clause=fid_which)

            # export polygon to new shapefile
            arcpy.CopyFeatures_management(lyr, poly)

            arcpy.env.outputCoordinateSystem = "PROJCS['Albers_Conical_Equal_Area', GEOGCS['GCS_North_American_1983', DATUM['D_North_American_1983', SPHEROID['GRS_1980', 6378137.0, 298.257222101]], PRIMEM['Greenwich', 0.0], UNIT['Degree', 0.0174532925199433]], PROJECTION['Albers'], PARAMETER['False_Easting',0.0], PARAMETER['False_Northing',0.0], PARAMETER['central_meridian',-96.0], PARAMETER['Standard_Parallel_1',29.5], PARAMETER['Standard_Parallel_2',45.5], PARAMETER['latitude_of_origin',23.0], UNIT['Meter',1.0]]"
            arcpy.gp.ExtractByMask_sa(geo_directory + input_name, poly, tif_export)

            # if name_vect[i] == "Seal_Beac":
             #   arcpy.CopyRaster_management(dem_temp1,
             #                               output_filepaths[h] + "/" + name_vect[i] + ".tif")





            # arcpy.TableToTable_conversion(dem_temp1, output_filepaths[h], name_vect[i])

            arcpy.ClearEnvironment("outputCoordinateSystem")
            # Print that the HUC8 clip is done
            print("  Environments cleared. Done.")

    h=h+1

print(failures)
