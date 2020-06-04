
# load dependencies
import math
import numpy as np
import arcpy
arcpy.CheckOutExtension('Spatial')

import os
import glob
import sys
import time
import logging
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *
import numpy as np

# Change this each time
# zone_n = 0
zone_n = 1
cores_to_use = 5

arcpy.env.compression = "LZ77"
tempWritePath = "D:/temp/ArcScratch"
tempGDB = tempWritePath + "/temp.gdb"
arcpy.env.scratchWorkspace = tempWritePath

# set up TEMP workspaces
rasterStorage = "D:/temp/ArcScratch/local_rast_wspace"
out_fishnet_path = "D:/temp/ArcScratch/fishnets/tempFishnet.shp"
out_raster_folder = "D:/z-star-analyses/forDaac/"

in_raster_list = ["D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_SB.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_MidAtlantic.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Northeast.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Northwest.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_SouthCentral.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Southeast.tif",
                  "D:/z-star-analyses/zStarStarRegions/zStarStar_lt1_30m_EEM_Southwest.tif"]

out_raster_list = ["probLowMarsh_SealBeach.tif",
                   "probLowMarsh_MidAtlantic.tif",
                   "probLowMarsh_Northeast.tif",
                   "probLowMarsh_Northwest.tif",
                   "probLowMarsh_SouthCentral.tif",
                   "probLowMarsh_Southeast.tif",
                   "probLowMarsh_Southwest.tif"]

# change this
in_raster_path = in_raster_list[zone_n]
out_raster_name = out_raster_list[zone_n]

arcpy.env.snapRaster = Raster(in_raster_path)
arcpy.env.outputCoordinateSystem = Raster(in_raster_path).spatialReference

# set up outpaths
in_raster_name = in_raster_path.split("/")[-1]

out_raster_path = out_raster_folder + "/" + out_raster_name

# create all the functions
# return p-value from z-score
def z2pPct(z):
    return int((0.5 * (1 + math.erf(z / math.sqrt(2)))*100)+0.5)

# vectorize function to act on an array
z2pPct = np.vectorize(z2pPct)

# function to calculate prob. below MHHWS based on input surface and uncertainty layers + LiDAR bias and uncertainty
def zstar_prop(zRast):

    zNumpy = arcpy.RasterToNumPyArray(zRast)
    pNumpy = z2pPct(zNumpy)
    del zNumpy

    lowerLeft = arcpy.Point(zRast.extent.XMin, zRast.extent.YMin)
    cellSize = zRast.meanCellWidth
    pRast = arcpy.NumPyArrayToRaster(pNumpy, lowerLeft, cellSize)

    del pNumpy

    OutTemp2 = arcpy.sa.ExtractByMask(pRast, zRast)

    return(OutTemp2)

def create_fishnet(in_raster_path, out_fc_path="D:/temp/ArcScratch/fishnets/tempFishnet.shp", blockSize=2500):
    # create raster obj from in path
    ras1 = arcpy.Raster(in_raster_path)
    arcpy.env.outputCoordinateSystem = ras1.spatialReference

    # specify input parameters for fishnets
    XMin = arcpy.GetRasterProperties_management(ras1, "LEFT").getOutput(0)
    XMax = arcpy.GetRasterProperties_management(ras1, "RIGHT").getOutput(0)
    YMin = arcpy.GetRasterProperties_management(ras1, "BOTTOM").getOutput(0)
    YMax = arcpy.GetRasterProperties_management(ras1, "TOP").getOutput(0)

    origCoord = "{} {}".format(XMin, YMin)
    YAxisCord = "{} {}".format(XMin, YMax)
    CornerCord = "{} {}".format(XMax, YMax)
    cellSizeW = "0"
    cellSizeH = "0"

    rastRowCount = float(arcpy.GetRasterProperties_management(ras1, "ROWCOUNT").getOutput(0))
    rastColCount = float(arcpy.GetRasterProperties_management(ras1, "COLUMNCOUNT").getOutput(0))

    numRows = math.ceil(float(rastRowCount) / float(blockSize))
    numCols = math.ceil(float(rastColCount) / float(blockSize))

    geo_type = "POLYGON"

    # Run Fishnet Tool
    logger.info("Running fishnet creator: {} with PID {}".format(current_process().name, os.getpid()))

    arcpy.CreateFishnet_management(out_fc_path, origCoord, YAxisCord, cellSizeW, cellSizeH, numRows, numCols, CornerCord, "NO_LABELS", "", geo_type)
    arcpy.ClearEnvironment("outputCoordinateSystem")

def execute_task(in_extent_Dict):

        #start clock
        time1 = time.clock()

        #get extent count and extents
        fc_count = (in_extent_Dict[0])
        procExt = (in_extent_Dict[1])
        XMin = (procExt[0])
        YMin = (procExt[1])
        XMax = (procExt[2])
        YMax = (procExt[3])

        #set environments
        arcpy.env.extent = arcpy.Extent(XMin, YMin, XMax, YMax)

        #send process info to logger
        logger.info("Running local math task: {} with PIC {}".format(current_process().name, os.getpid()))

        #run the local task
        ras_out = zstar_prop(zRast=Raster(in_raster_path))

        #clear the extent environment
        arcpy.ClearEnvironment("extent")
        out_name = "out_pRast{}.tif".format(fc_count)
        out_path = rasterStorage + "/" + out_name
        ras_out.save(out_path)

        #end the clock
        time2 = time.clock()
        logger.info("{} with PIC {} finished in {}".format(current_process().name, os.getpid(), str(time2-time1)))# Set Geoprocessing environments

print("Starting " + out_raster_name)

# create a logger to report
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel((logging.DEBUG))
ch.setFormatter(formatter)
logger.addHandler(ch)

if __name__ == "__main__": # If this is the main file run this last. If this file is loaded by another file, don't run.

    arcpy.env.outputCoordinateSystem = Raster(in_raster_path).spatialReference

    # start clock
    time_start = time.clock()
    # call create fishnet
    logger.info("Creating fishnet features...")
    create_fishnet(in_raster_path, out_fishnet_path)

    #  extent of individual features and add it to a dictionary
    extDict = {} # empty dictionary
    count = 1 # iterating through fishnet
    for row in arcpy.da.SearchCursor(out_fishnet_path, ["SHAPE@"]): # set up a search cursor to get all of the polygon extents
        extent_curr = row[0].extent
        ls=[] # empty list
        #

        xoverlap = (extent_curr.XMax - extent_curr.XMin) / 200
        yoverlap = (extent_curr.YMax - extent_curr.YMin) / 200
        ls.append(extent_curr.XMin - xoverlap)
        ls.append(extent_curr.YMin - yoverlap)
        ls.append(extent_curr.XMax + xoverlap)
        ls.append(extent_curr.YMax + yoverlap)
        extDict[count] = ls
        count+=1

    # create a process pool
    pool = Pool(processes=cores_to_use)
    pool.map(execute_task, extDict.items())
    pool.close()
    pool.join()

    print("   Creating Mosaic dataset...")
    arcpy.CreateMosaicDataset_management(tempWritePath + "/temp.gdb", "tempMosaic", Raster(in_raster_path).spatialReference)
    # add files from temp rast workspace
    arcpy.management.AddRastersToMosaicDataset("D:/temp/ArcScratch/temp.gdb/tempMosaic", "Raster Dataset", rasterStorage)
    # Calculate Statistics
    print("   Calculating statistics...")
    arcpy.management.CalculateStatistics("D:/temp/ArcScratch/temp.gdb/tempMosaic", x_skip_factor=10, y_skip_factor=10)
    # export mosaic to raster dataset
    arcpy.management.DefineMosaicDatasetNoData("D:/temp/ArcScratch/temp.gdb/tempMosaic", 1, "BAND_1 -1")

    # writing to mosaic to temp drive
    print("   Writing mosaic to temp drive...")
    arcpy.CopyRaster_management("D:/temp/ArcScratch/temp.gdb/tempMosaic", out_raster_path)
    arcpy.management.BuildRasterAttributeTable(out_raster_path, "Overwrite")

    print("   done.")

    # clear folder of temp raster files
    time_end = time.clock()
    logger.info("Time taken in main in second(s) is: {}".format(str(time_end-time_start)))

    # delete blocks
    dirPath = 'D:/temp/ArcScratch/local_rast_wspace'
    fileList = os.listdir(dirPath)
    for fileName in fileList:
        os.remove(dirPath+"/"+fileName)
    print("   blocks deleted.")
