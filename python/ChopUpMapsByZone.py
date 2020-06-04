# Need to chop up files by a few different states

# Import all necessary libraries
import arcpy
import os
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# establish scratch workspace
scratch_path = "F:/temp/ArcScratch"
arcpy.env.scratchWorkspace = scratch_path
lyr = "lyr"
poly = "in_memory/temp"

arcpy.env.extent = "MINOF"
arcpy.env.compression = "LZ77"

# First establish a directory where all of the large geographic data files are held
geo_directory = "F:/z-star-spatial-data/derivative-maps/zStar/"

# State input files
states = "F:/z-star-spatial-data/input-layers/AOI/coastal_states2_dissolved.shp"

# Create a list of input file names
input_names = ["zStarStar_lt1_30m_maskedByCcapEEM2010_200508.tif"]

output_names = ["zStarStar_lt1_30m_EEM"]

# Create a list of output files names
output_filepaths = ["F:/z-star-spatial-data/derivative-maps/zStar/zStarStarRegions"]

# Iterate through the states
# set up export tags
names = arcpy.da.SearchCursor(states, ["REGION2"])
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
        tif_import = geo_directory + input_name
        tif_export = output_filepaths[h] + "/" + output_names[h] + "_" + name_vect[i] + ".tif"

        print("  " + tif_import)
        print("  " + tif_export)

        # if a file is already done, ignore it
        if os.path.isfile(tif_export):
            print("  " + tif_export + " already done.")
        else:
            # Make a layer from the feature class
            arcpy.MakeFeatureLayer_management(states, lyr)
            fid_which = ' "FID" = ' + str(i) + ' '

            # select the row
            arcpy.SelectLayerByAttribute_management(lyr, selection_type="NEW_SELECTION", where_clause=fid_which)

            # export polygon to new shapefile
            arcpy.CopyFeatures_management(lyr, poly)

            arcpy.env.outputCoordinateSystem = "PROJCS['Albers_Conical_Equal_Area', GEOGCS['GCS_North_American_1983', DATUM['D_North_American_1983', SPHEROID['GRS_1980', 6378137.0, 298.257222101]], PRIMEM['Greenwich', 0.0], UNIT['Degree', 0.0174532925199433]], PROJECTION['Albers'], PARAMETER['False_Easting',0.0], PARAMETER['False_Northing',0.0], PARAMETER['central_meridian',-96.0], PARAMETER['Standard_Parallel_1',29.5], PARAMETER['Standard_Parallel_2',45.5], PARAMETER['latitude_of_origin',23.0], UNIT['Meter',1.0]]"
            arcpy.env.snapRaster = Raster(tif_import)

            outExtractByMask = ExtractByMask(tif_import, poly)
            outExtractByMask.save(tif_export)

            arcpy.ClearEnvironment("outputCoordinateSystem")
            # Print that the clip is done
            print("  Environments cleared. Done.")

    h=h+1

print(failures)
