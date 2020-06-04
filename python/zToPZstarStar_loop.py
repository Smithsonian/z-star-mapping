
import arcpy
import numpy
import os

# load dependencies
import math
import arcpy
arcpy.CheckOutExtension('Spatial')

from os import path

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *
import numpy as np

blocksize=2500
arcpy.env.compression = "LZ77"

in_raster_list = [#"D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_SB.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_MidAtlantic.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Northwest.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_SouthCentral.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Southeast.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Southwest.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Northeast.tif"]

out_raster_list = [#"D:/z-star-analyses/forDaac/probLowMarsh_SealBeach.tif",
                   "D:/z-star-analyses/forDaac/probLowMarsh_MidAtlantic.tif",
                   "D:/z-star-analyses/forDaac/probLowMarsh_Northwest.tif",
                   "D:/z-star-analyses/forDaac/probLowMarsh_SouthCentral.tif",
                   "D:/z-star-analyses/forDaac/probLowMarsh_Southeast.tif",
                   "D:/z-star-analyses/forDaac/probLowMarsh_Southwest.tif",
                   "D:/z-star-analyses/forDaac/probLowMarsh_Northeast.tif"]

# create all the functions
# return p-value from z-score
def z2pPct(z):
    return int((0.5 * (1 + math.erf(z / math.sqrt(2)))*100)+0.5)

# vectorize function to act on an array
z2pPct = np.vectorize(z2pPct)

for i in range(len(in_raster_list)):

    myRaster = arcpy.Raster(in_raster_list[i])
    fileout = out_raster_list[i]

    if path.exists(fileout):
        print(fileout + " already done.")
    else:

        arcpy.env.mask = myRaster
        arcpy.env.snapRaster = myRaster
        arcpy.env.outputCoordinateSystem = myRaster.spatialReference

        print(myRaster)

        # Copy from
        filelist = []
        blocknum = 0
        for x in range(0, myRaster.width, blocksize):
            for y in range(0, myRaster.height, blocksize):

                print("     " + str(blocknum))

                # Lower left coordinate of block (in map units)
                mx = myRaster.extent.XMin + x * myRaster.meanCellWidth
                my = myRaster.extent.YMin + y * myRaster.meanCellHeight
                # Upper right coordinate of block (in cells)
                lx = min([x + blocksize, myRaster.width])
                ly = min([y + blocksize, myRaster.height])
                #   noting that (x, y) is the lower left coordinate (in cells)

                # Extract data block
                zNumpy = arcpy.RasterToNumPyArray(myRaster, arcpy.Point(mx, my),
                                                  lx-x, ly-y)

                pNumpy = z2pPct(zNumpy)


                # Convert data block back to raster
                myRasterBlock = arcpy.NumPyArrayToRaster(pNumpy, arcpy.Point(mx, my),
                                                         myRaster.meanCellWidth,
                                                         myRaster.meanCellHeight)

                # Save on disk temporarily as 'filename_#.ext'
                filetemp = ('_%i.' % blocknum).join(fileout.rsplit('.',1))
                myRasterBlock.save(filetemp)

                # Maintain a list of saved temporary files
                filelist.append(filetemp)
                blocknum += 1
                print("     " + filetemp, " done.")

        print(filelist)
        # Mosaic temporary files
        arcpy.Mosaic_management(';'.join(filelist[1:]), filelist[0])
        if arcpy.Exists(fileout):
            arcpy.Delete_management(fileout)
        arcpy.Rename_management(filelist[0], fileout)

        # Remove temporary files
        for fileitem in filelist:
            if arcpy.Exists(fileitem):
                arcpy.Delete_management(fileitem)

        # Release raster objects from memory
        del myRasterBlock
        del myRaster
        # ----------------------------------------------------------------------------

        arcpy.management.BuildRasterAttributeTable(arcpy.Raster(fileout), "Overwrite")
        print(fileout)
