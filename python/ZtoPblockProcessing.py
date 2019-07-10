"""
Multiprocessing from...
http://proceedings.esri.com/library/userconf/devsummit17/papers/dev_int_39.pdf
"""

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

deleteHighRes = True

arcpy.env.scratchWorkspace = "in_memory"
tempWritePath = "D:/temp/ArcScratch"
tempGDB = tempWritePath + "/temp.gdb"

# set up TEMP workspaces
rasterStorage = "D:/temp/ArcScratch/local_rast_wspace"
out_fishnet_path = "D:/temp/ArcScratch/fishnets/tempFishnet.shp"
out_raster_folder = "D:/z-star-spatial-data/derivative-maps/pMHHWS_V2p0/"

# change this
in_raster_path = "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_MTR_STO_dems/CA_MTR_STO_dems/CA_MTR3_GCS_5m_NAVD88m.img"
min_year = 2009
max_year = 2011
hydroflattening_value = -1


# set up outpaths
in_raster_name = in_raster_path.split("/")[-1]

# Cut out other unnecessary pieces
stringsToCutOut = ["_GCS", "_5m", "_10m", "_3m", "_NAVDm", "_NAVD88m", "_dist", "_Topobathy_DEM", "JH", "_5ft", "_m", "_2m"]

out_raster_name = in_raster_name

if (('.tif' not in out_raster_name) and (".img" not in out_raster_name)):
    out_raster_name = out_raster_name + ".img"

for ii in stringsToCutOut:
        out_raster_name = out_raster_name.replace(ii, "")

if min_year == max_year:
    out_raster_name = "pMHHWS" + "_" + str(min_year) + "_"  + out_raster_name
else:
    out_raster_name = "pMHHWS" + "_"  + str(min_year) + "_" + str(max_year) + "_"  + out_raster_name

mid_raster_path = out_raster_folder + "originalRes/" + out_raster_name
out_raster_path = out_raster_folder + "30m/" + out_raster_name

# create all the functions
# return p-value from z-score
def z2p(z):
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))

# vectorize function to act on an array
z2p = np.vectorize(z2p)

# function to calculate prob. below MHHWS based on input surface and uncertainty layers + LiDAR bias and uncertainty
def mhhws_prop(in_dem,
        rmse_layer = Raster("D:/z-star-spatial-data/input-layers/datum-layers/MHHWS_NAVD_propegated_uncertainty_300m_19705.img"),
        mhhw_layer = Raster("D:/z-star-spatial-data/input-layers/datum-layers/MHHW_NAVD_300m_190531.img"),
        mhhws_layer = Raster("D:/z-star-spatial-data/input-layers/datum-layers/MHHWS_MHHW_300m_190705.img"),
        wetland_layer = Raster("D:/z-star-spatial-data/input-layers/CCAPallWetlands2010.tif"),
        wetland_offset = 0.173,
        maskLayer = Raster("D:/z-star-spatial-data/input-layers/AOI/noWater2006And2010.img")):

    # set env
    arcpy.env.outputCoordinateSystem = Raster(in_dem).spatialReference
    arcpy.env.cellSize = "MINOF"
    arcpy.env.mask = maskLayer

    # Simplify raster so that anything below
    tempRast1 = arcpy.sa.Con((Raster(in_dem) > -99) & (Raster(in_dem) != float(hydroflattening_value)), Raster(in_dem))
    tempRast2 = arcpy.sa.Con(IsNull(wetland_layer), tempRast1, (tempRast1 - float(wetland_offset)))

    zRast = (mhhw_layer + mhhws_layer - tempRast2) / rmse_layer # water surface - elevation surface / rmse
    del tempRast2

    zRastCon = arcpy.sa.Con(zRast <-4, -4, arcpy.sa.Con(zRast>4, 4, zRast))
    del zRast

    zNumpy = arcpy.RasterToNumPyArray(zRastCon)
    pNumpy = z2p(zNumpy)
    del zNumpy

    lowerLeft = arcpy.Point(zRastCon.extent.XMin, zRastCon.extent.YMin)
    cellSize = zRastCon.meanCellWidth
    pRast = arcpy.NumPyArrayToRaster(pNumpy, lowerLeft, cellSize)
    del zRastCon
    del pNumpy

    # write back in cons statement to clip out all 0's
    OutTemp = arcpy.sa.ExtractByMask(pRast, tempRast1)
    del pRast
    del tempRast1

    return(OutTemp)

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
        ras_out = mhhws_prop(in_dem=in_raster_path)

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
    logger.info("Createing fishnet features...")
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
    #print(extDict)

    # create a process pool
    pool = Pool(processes=11)
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
    arcpy.management.DefineMosaicDatasetNoData("D:/temp/ArcScratch/temp.gdb/tempMosaic", 1, "BAND_1 0")

    # writing to mosaic to temp drive
    print("   Writing mosaic to temp drive...")
    arcpy.CopyRaster_management("D:/temp/ArcScratch/temp.gdb/tempMosaic", mid_raster_path)

    # set up environments for coarser resolution resampling
    arcpy.env.snapRaster = Raster("D:/z-star-spatial-data/input-layers/AOI/noWater2006And2010.img")
    arcpy.env.outputCoordinateSystem = Raster("D:/z-star-spatial-data/input-layers/AOI/noWater2006And2010.img").spatialReference
    print("   Resampling...")

    resampledTemp = "in_memory/resampledTemp"
    arcpy.Resample_management(mid_raster_path, resampledTemp, "30")

    final_raster = arcpy.sa.Con(Raster(resampledTemp)>0, Raster(resampledTemp))
    final_raster.save(out_raster_path)

    # clear environments
    arcpy.ClearEnvironment("snapRaster")
    arcpy.ClearEnvironment("outputCoordinateSystem")

    del resampledTemp
    print("   done.")

    if deleteHighRes == True:
        del mid_raster_path
        print("   deleted Original Res.")

    # clear folder of temp raster files
    time_end = time.clock()
    logger.info("Time taken in main in second(s) is: {}".format(str(time_end-time_start)))

    # delete blocks
    dirPath = r'D:/temp/ArcScratch/local_rast_wspace'
    fileList = os.listdir(dirPath)
    for fileName in fileList:
        os.remove(dirPath+"/"+fileName)
    print("   blocks deleted.")
