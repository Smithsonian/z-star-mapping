# load dependencies
import math
import arcpy
arcpy.CheckOutExtension('Spatial')

from os import path

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *
import numpy as np

arcpy.env.compression = "LZ77"
tempWritePath = "D:/temp/ArcScratch"
arcpy.env.scratchWorkspace = tempWritePath

# set up TEMP workspaces
out_raster_folder = "D:/z-star-analyses/forDaac/"

in_raster_list = ["D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_SB.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_MidAtlantic.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Northwest.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_SouthCentral.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Southeast.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Southwest.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Northeast.tif"]

out_raster_list = ["probLowMarsh_SealBeach.tif",
                   "probLowMarsh_MidAtlantic.tif",
                   "probLowMarsh_Northeast.tif",
                   "probLowMarsh_Northwest.tif",
                   "probLowMarsh_SouthCentral.tif",
                   "probLowMarsh_Southeast.tif",
                   "probLowMarsh_Southwest.tif"]

# create all the functions
# return p-value from z-score
def z2pPct(z):
    return int((0.5 * (1 + math.erf(z / math.sqrt(2)))*100)+0.5)

# vectorize function to act on an array
z2pPct = np.vectorize(z2pPct)

# function to calculate prob. below MHHWS based on input surface and uncertainty layers + LiDAR bias and uncertainty
def zstar_prop(zRast):

    print("     ... Converting raster to array.")

    zNumpy = arcpy.RasterToNumPyArray(zRast)
    print("     Raster to Array done.")
    print("     ... Running Z to P.")
    pNumpy = z2pPct(zNumpy)
    del zNumpy

    print("      Z to P done.")
    print("     ... Converting Array to Raster.")

    lowerLeft = arcpy.Point(zRast.extent.XMin, zRast.extent.YMin)
    cellSize = zRast.meanCellWidth
    pRast = arcpy.NumPyArrayToRaster(pNumpy, lowerLeft, cellSize)

    print("     Array converted back to raster.")

    del pNumpy

    print("     .... writing final file")

    OutTemp2 = arcpy.sa.ExtractByMask(pRast, zRast)
    OutTemp2.save(out_raster_path)
    arcpy.management.BuildRasterAttributeTable(out_raster_path, "Overwrite")
    del OutTemp2

    print("     Final file done.")

#run the local task

for i in range(len(in_raster_list)):

    # change this each time
    in_raster_path = in_raster_list[i]
    out_raster_name = out_raster_list[i]

    arcpy.env.snapRaster = Raster(in_raster_path)
    arcpy.env.outputCoordinateSystem = Raster(in_raster_path).spatialReference

    # set up outpaths
    out_raster_path = out_raster_folder + "/" + out_raster_name

    if path.exists(out_raster_path):

        print("! " + out_raster_path + " already done !")

    else:
        print(in_raster_path + "...")

        zstar_prop(zRast=Raster(in_raster_path))

        print(out_raster_path + " done.")
