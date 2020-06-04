import arcpy
arcpy.CheckOutExtension('Spatial')
import arcpy.sa
import os
import pandas as pd

taskList = pd.DataFrame.from_csv("D:/GitHub/z-star-mapping/python/zStarFullWorkOrder200427.csv")

arcpy.env.scratchWorkspace = "in_memory"
tempWritePath = "D:/temp/ArcScratch"
tempGDB = tempWritePath + "/temp.gdb"
arcpy.env.overwriteOutput = True

stringsToCutOut = ["_GCS", "_5m", "_10m", "_3m", "_NAVDm", "_NAVD88m", "_dist", "_Topobathy_DEM", "JH", "_5ft", "_m", "_2m", "_30m"]

def zStarWithUncertainty(zLayer,
                         filetag,
                         min_year,
                         max_year,
                         mslLayer = "F:/z-star-spatial-data/input-layers/datum-layers/MSL_NAVD_300m_190531.img",
                         mslSeLayer = "F:/z-star-spatial-data/input-layers/datum-layers/MSL_NAVD_propegated_uncertainty_300m_190603.img",
                         mhwLayer = "F:/z-star-spatial-data/input-layers/datum-layers/MHW_NAVD_300m_190531.img",
                         mhwSeLayer = "F:/z-star-spatial-data/input-layers/datum-layers/MHW_NAVD_propegated_uncertainty_300m_190829.img",
                         corMhwMsl = 0.873,
                         wetlandDemBias = 0.173, wetlandDemRMSE = 0.233,
                         hydroFlatteningValue = -1,
                         naValue = -9999,
                         outputPath = "F:/z-star-spatial-data/derivative-maps/zStar/national-scale",
                         wetlandLayer = "F:/z-star-spatial-data/input-layers/CCAPallWetlands2010.tif",
                         maskLayer = "F:/z-star-spatial-data/input-layers/AOI/noWater2006And2010.img",
                         mapPartialDerivatives = False
                         ):

    # set extent
    arcpy.env.extent = arcpy.sa.Raster(zLayer).extent

    # set up output's file name
    in_raster_name = zLayer.split("/")[-1]

    if '.' in in_raster_name:
        out_raster_name = in_raster_name.split(".")[-2]
    else:
        out_raster_name = in_raster_name

    for ii in stringsToCutOut:
        out_raster_name = out_raster_name.replace(ii, "")

    if min_year == max_year:
        out_raster_name = str(int(float(min_year))) + "_"  + out_raster_name
    else:
        out_raster_name = str(int(float(min_year))) + "_" + str(int(float(max_year))) + "_"  + out_raster_name

    # Simplify raster so that anything below
    # Clean out hydroflatteing

    # set all of the output filepaths
    zStarRast_path = outputPath + "/z-star/zStar_" + filetag + "_" + out_raster_name + ".img"
    totalUncertainty_path = outputPath + "/total-uncertainty/totalUncertainty_" + filetag + "_" + out_raster_name +  ".img"
    categoricalUncertainty_path =  outputPath + "/categorical-uncertainty/categoricalUncertainty_" + filetag + "_" +  out_raster_name + ".img"

    #if os.path.isfile(zStarRast_path) & os.path.isfile(totalUncertainty_path) & os.path.isfile(categoricalUncertainty_path):
    #    print(out_raster_name + " already done.")
    #else:
    print("Running " + out_raster_name)

    # Set resolution, mask and coordinates to the wetland target layer
    if os.path.isfile(wetlandLayer):

        arcpy.env.outputCoordinateSystem = arcpy.sa.Raster(maskLayer).spatialReference
        arcpy.env.cellSize = arcpy.sa.Raster(maskLayer).meanCellHeight
        arcpy.env.mask = arcpy.sa.Raster(maskLayer)
        arcpy.env.snapRaster = arcpy.sa.Raster(maskLayer)

    else:

        arcpy.env.outputCoordinateSystem = arcpy.sa.Raster(zLayer).spatialReference
        arcpy.env.cellSize = arcpy.sa.Raster(zLayer).cellSize
        arcpy.env.snapRaster = arcpy.sa.Raster(zLayer)


    print("    Parameters set.")

    print("    Removing hydroflattening, adjusting resolution, and masking out non-wetlands... ")

    tempRast1 = arcpy.sa.Con((arcpy.sa.Raster(zLayer) != float(hydroFlatteningValue)),
                             arcpy.sa.Con(arcpy.sa.Raster(zLayer) > float(naValue),
                                          arcpy.sa.Raster(zLayer)))

    print("    Accounting for " + str(wetlandDemBias) + "m of DEM bias ... ")

    # Account for DEM bias
    tempRast2 = arcpy.sa.Con(arcpy.sa.IsNull(wetlandLayer), tempRast1, (tempRast1 - float(wetlandDemBias)))
    del(tempRast1)

    print("    tempRaster1 deleted.")

    # Calculate Z-star layer (Zstar = MHW-Z/MHW-MSL)
    # water surface - elevation surface / 1/2 tide range
    print("    Calculating z-star raster...")
    # zStarRast = ((tempRast2-arcpy.sa.Raster(mslLayer))/(arcpy.sa.Raster(mhwLayer)-arcpy.sa.Raster(mslLayer)))

    # Write Z-star to file
    print("    Writing z-star raster...")
    # zStarRast.save(zStarRast_path)

    # del(zStarRast)
    print("    Z-star rater deleted.")

    # Calculate Z partial derivative

    print("    Calculating partial derivative rasters...")

    # dZ*/dZ = 1/(MHW-MSL)
    dZstarOverDZ = float(1)/(arcpy.sa.Raster(mhwLayer)-arcpy.sa.Raster(mslLayer))

    # dZ*/dMHW = -(Z-MSL)/(MHW-MSL)^2
    dZstarOverDMHW = (float(-1)*((tempRast2-arcpy.sa.Raster(mslLayer)) / arcpy.sa.Square((arcpy.sa.Raster(mhwLayer)-arcpy.sa.Raster(mslLayer)))))

    # dZ*/dMSL =  Z-MHW / (MHW-MSL)^2
    dZstarOverDMSL = ((tempRast2 - arcpy.sa.Raster(mhwLayer)) / arcpy.sa.Square((arcpy.sa.Raster(mhwLayer)-arcpy.sa.Raster(mslLayer))))

    print("    Calculating propegated uncertainties ...")

    if mapPartialDerivatives == True:
        print("    Writing partial derivative rasters ...")

        dZstarOverDZ_path = outputPath + "/partial-derivatives/dZstarOverDZ_" + filetag + "_" + out_raster_name + ".img"
        dZstarOverDZ.save(dZstarOverDZ_path)

        dZstarOverDMHW_path = outputPath + "/partial-derivatives/dZstarOverDMHW_" + filetag + "_" + out_raster_name + ".img"
        dZstarOverDMHW.save(dZstarOverDMHW_path)

        dZstarOverDMSL_path = outputPath + "/partial-derivatives/dZstarOverDMSL_" + filetag + "_" + out_raster_name + ".img"
        dZstarOverDMSL.save(dZstarOverDMSL_path)

    # Calculate propegated uncertainty sigma_z*^2 =
    propegatedUncertainty = arcpy.sa.SquareRoot((arcpy.sa.Square(dZstarOverDZ) * arcpy.sa.Square(wetlandDemRMSE)) +
                                       (arcpy.sa.Square(dZstarOverDMHW) * arcpy.sa.Square(arcpy.sa.Raster(mhwSeLayer))) +
                                       (arcpy.sa.Square(dZstarOverDMSL) * arcpy.sa.Square(arcpy.sa.Raster(mslSeLayer))) +
                                       (float(2) * arcpy.sa.Raster(mhwSeLayer) * arcpy.sa.Raster(mslSeLayer) * dZstarOverDMHW * dZstarOverDMSL * float(corMhwMsl))
                                       )

    print("    Writing propagated uncertainty raster ...")
    propegatedUncertainty.save(totalUncertainty_path)
    del(propegatedUncertainty)

    print("    Propagated uncertainty raster deleted.")

    print("    Calculating partial errors ...")

    partialZ_error = arcpy.sa.SquareRoot(arcpy.sa.Square(dZstarOverDZ) * arcpy.sa.Square(wetlandDemRMSE))

    partialDatum_error = arcpy.sa.SquareRoot((arcpy.sa.Square(dZstarOverDMHW) * arcpy.sa.Square(mhwSeLayer)) +
                                       (arcpy.sa.Square(dZstarOverDMSL) * arcpy.sa.Square(mslSeLayer)) +
                                       (float(2) * arcpy.sa.Raster(mhwSeLayer) * arcpy.sa.Raster(mslSeLayer) * dZstarOverDMHW * dZstarOverDMSL * float(corMhwMsl))
                                    )

    del(dZstarOverDZ)
    del(dZstarOverDMHW)
    del(dZstarOverDMSL)

    print("    Partial derivative rasters deleted.")

    print("    Calculating categorical uncertainty ...")
    zErrorGtDatum = arcpy.sa.Con(partialZ_error>partialDatum_error, 1, 2)

    print("    Writing categorical uncertainty raster ...")
    zErrorGtDatum.save(categoricalUncertainty_path)

    del(zErrorGtDatum)
    print("    Categorical uncertainty raster deleted.")

for index, row in taskList.iterrows():

    try:

        zStarWithUncertainty(zLayer = row["demPath"],
                             filetag = "MHW",
                             min_year = str(row["surveyYearMin"]),
                             max_year = str(row["surveyYearMax"]),
                             hydroFlatteningValue = row["hydroflatteningValue"])

    except:

        print("!!! failure !!!")
