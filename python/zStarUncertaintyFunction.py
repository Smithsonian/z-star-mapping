# Some file info started 3 June 2019. Written by James R Holmquist
# Goal is to take
# Function needs to import an NAVD88 elevation surface, a MHW or MHHW surface, a MSL surface,
# MHW/MHHW and MSL propegated uncertinaty layers,
# a correlation coefficient between MHW/MHHW and MSL, a wetland specific DEM bias, a wetland specific DEM RMSE,
# a raster specific NA or 'water flattening' value,
# a wetland yes/no layer, and an output file path.

import arcpy
arcpy.CheckOutExtension('Spatial')
import arcpy.sa
import os

arcpy.env.scratchWorkspace = "in_memory"
tempWritePath = "D:/temp/ArcScratch"
tempGDB = tempWritePath + "/temp.gdb"
arcpy.env.overwriteOutput = True

stringsToCutOut = ["_GCS", "_5m", "_10m", "_3m", "_NAVDm", "_NAVD88m", "_dist", "_Topobathy_DEM", "JH", "_5ft", "_m", "_2m"]

inputRasters = ["D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_ChannelIslands_dems/CA_ChannelIslands_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_EKA_dems/CA_EKA1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_EKA_dems/CA_EKA2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_MTR_STO_dems/CA_MTR_STO_dems/CA_MTR1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_MTR_STO_dems/CA_MTR_STO_dems/CA_MTR2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_MTR_STO_dems/CA_MTR_STO_dems/CA_MTR3_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_LOX_dems/CA_LOX_dems/CA_LOX1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_LOX_dems/CA_LOX_dems/CA_LOX2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/CA/CA_SGX_dems/CA_SGX_dems/CA_SGX_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/ALFL_MOB_TLH_dems/ALFL_MOB_TLH_dems/ALFL_A1_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/ALFL_MOB_TLH_dems/ALFL_MOB_TLH_dems/ALFL_A2_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/ALFL_MOB_TLH_dems/ALFL_MOB_TLH_dems/ALFL_A3_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/FL_JAX_dems/FL_JAX_dems/FL_JAX_1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/FL_JAX_dems/FL_JAX_dems/FL_JAX_2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/FL_MFL_dems/FL_MFL_dems/FL_MFL_1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/FL_MFL_dems/FL_MFL_dems/FL_MFL_2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/FL_MLB_dems/FL_MLB_dems/FL_MLB_1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/FL_MLB_dems/FL_MLB_dems/FL_MLB_2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/FL_TBW_dems/FL_TBW_dems/FL_TBW_1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/FL/FL_TBW_dems/FL_TBW_dems/FL_TBW_2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/GA/GA_dems/GA_dems/GA_CHS_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/GA/GA_dems/GA_dems/GA_JAX_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/MD/MD_LWX_dems/MD_LWX_GCS_10m_NAVDm_dist.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/MD/MD_LWX_dems/MD_PrinceGeorges_GCS_10m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/MD/MD_PHI_AKQ_dems/MD_PHI_AKQ_GCS_10m_NAVDm_dist.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/MS/MS_dems/MS_GCS_10m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/NC/NC_dems/NC_GCS_10m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/NewEngland/CT_dems/CT_GCS_10m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/NewEngland/MARI_BOX_dems/MARI_BOX_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/NewEngland/ME_CAR_dems/ME_CAR_GCS_10m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/NewEngland/MENH_GYX_dems/MENH_GYX_GCS_10m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/NJDEPA/NJDEPA_PHI_dems/NJDEPA_PHI_GCS_10m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/NY/NY_OKX_dems/NY_OKX1_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/NY/NY_OKX_dems/NY_OKX2_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/OR/OR_MFR_dems/OR_MFR_dems/OR_MFR_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/OR/OR_PQR_dems/OR_PQR1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/OR/OR_PQR_dems/OR_PQR2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/SC/SC_CHS_dems/SC_CHS_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/SC/SC_ILM_dems/SC_ILM_GCS_5m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/TX/TX_BRO_dems/TX_BRO_GCS_10m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/TX/TX_CRP_dems/TX_CRP_GCS_10m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/TX/TX_HGX_dems/TX_HGX_GCS_10m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/TX/TX_LCH_dems/TX_LCH_GCS_10m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/WA/WA_PQR_dems/WA_PQR_dems/WA_PQR_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs_NEW_170424/DC_dem/DC_GCS_3m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs_NEW_170424/VA_EasternShore_dem/VA_EasternShore_GCS_3m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs_NEW_170424/VA_Middle_dem/VA_Middle_GCS_3m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs_NEW_170424/VA_Northern_dem/VA_Northern_GCS_3m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs_NEW_170424/VA_Southern_dem/VA_Southern_GCS_3m_NAVDm.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs_NEW_170424/SC_Horry_dem/SC_Horry_GCS_3m_NAVDm.img",
                "D:/LIDAR DEMs/NGOM Topobathymetry/Northern_Gulf_of_Mexico_Topobathy_DEM.tif",
                "D:/LIDAR DEMs/SJ_Delta/dem_bay_delta_10m_v3_20121109_2/dem_bay_delta_10m_2JH.tif",
                "D:/LIDAR DEMs/SJ_Delta/dem_bay_delta_10m_v3_20121109_4/dem_bay_delta_10m_4JH.tif",
                "D:/LIDAR DEMs/MD_state/MD_highres/Harford_DEM_2013_1.5m/harford_5ft",
                "D:/LIDAR DEMs/MD_state/MD_highres/Baltimore_DEM_2015_2.5ft/baltimore2015.tif",
                "D:/LIDAR DEMs/MD_state/MD_highres/Calvert_DEM_2011_2m/calvert_2m",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/WA/WA_SEW_dems/WA_SEW1_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/WA/WA_SEW_dems/WA_SEW2_GCS_5m_NAVD88m.img",
                "D:/LIDAR DEMs/NOAA_OCM_SLR_Inundation_DEMs/WA/WA_SEW_dems/WA_SEW3_GCS_5m_NAVD88m.img"
                ]

def zStarWithUncertainty(zLayer,
                         mslLayer = "D:/z-star-spatial-data/input-layers/datum-layers/MSL_NAVD_300m_190531.img",
                         mslSeLayer = "D:/z-star-spatial-data/input-layers/datum-layers/MSL_NAVD_propegated_uncertainty_300m_190603.img",
                         mhwLayer = "D:/z-star-spatial-data/input-layers/datum-layers/MHW_NAVD_300m_190531.img",
                         mhwSeLayer = "D:/z-star-spatial-data/input-layers/datum-layers/MHW_NAVD_propegated_uncertainty_300m_190603.img",
                         corMhwMsl = 0.873,
                         wetlandDemBias = 0.173, wetlandDemRMSE = 0.233,
                         hydroFlatteningValue = -1,
                         naValue = -9999,
                         outputPath = "D:/z-star-spatial-data/derivative-maps/national-scale",
                         wetlandLayer = "D:/z-star-spatial-data/input-layers/CCAPallWetlands2010.tif",
                         maskLayer = "D:/z-star-spatial-data/input-layers/zStarAoiMask2_190612.tif",
                         filetag = "",
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

    print("Running " + out_raster_name)

    # Simplify raster so that anything below
    # Clean out hydroflatteing

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
    zStarRast = ((tempRast2-arcpy.sa.Raster(mslLayer))/(arcpy.sa.Raster(mhwLayer)-arcpy.sa.Raster(mslLayer)))

    # Write Z-star to file
    print("    Writing z-star raster...")
    zStarRast_path = outputPath + "/z-star/zStar_" + out_raster_name + filetag + ".img"
    zStarRast.save(zStarRast_path)

    del(zStarRast)
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

        dZstarOverDZ_path = outputPath + "/partial-derivatives/dZstarOverDZ_" + out_raster_name + filetag + ".img"
        dZstarOverDZ.save(dZstarOverDZ_path)

        dZstarOverDMHW_path = outputPath + "/partial-derivatives/dZstarOverDMHW_" + out_raster_name + filetag + ".img"
        dZstarOverDMHW.save(dZstarOverDMHW_path)

        dZstarOverDMSL_path = outputPath + "/partial-derivatives/dZstarOverDMSL_" + out_raster_name + filetag + ".img"
        dZstarOverDMSL.save(dZstarOverDMSL_path)

    # Calculate propegated uncertainty sigma_z*^2 =
    propegatedUncertainty = arcpy.sa.SquareRoot((arcpy.sa.Square(dZstarOverDZ) * arcpy.sa.Square(wetlandDemRMSE)) +
                                       (arcpy.sa.Square(dZstarOverDMHW) * arcpy.sa.Square(arcpy.sa.Raster(mhwSeLayer))) +
                                       (arcpy.sa.Square(dZstarOverDMSL) * arcpy.sa.Square(arcpy.sa.Raster(mslSeLayer))) +
                                       (float(2) * arcpy.sa.Raster(mhwSeLayer) * arcpy.sa.Raster(mslSeLayer) * dZstarOverDMHW * dZstarOverDMSL * float(corMhwMsl))
                                       )

    print("    Writing propagated uncertainty raster ...")

    totalUncertainty_path = outputPath + "/total-uncertainty/totalUncertainty_" + out_raster_name + filetag + ".img"
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
    categoricalUncertainty_path =  outputPath + "/categorical-uncertainty/categoricalUncertainty_" + out_raster_name + filetag + ".img"

    print("    Writing categorical uncertainty raster ...")
    zErrorGtDatum.save(categoricalUncertainty_path)

    del(zErrorGtDatum)
    print("    Categorical uncertainty raster deleted.")

#zStarWithUncertainty(inputRasters[0],  filetag="_MHW")
#zStarWithUncertainty(inputRasters[8],  filetag="_MHW", mapPartialDerivatives=True)
#zStarWithUncertainty(inputRasters[1],  filetag="_MHW")
#zStarWithUncertainty(inputRasters[2],  filetag="_MHW")
#zStarWithUncertainty(inputRasters[3],  filetag="_MHW")
#zStarWithUncertainty(inputRasters[4],  filetag="_MHW")
#zStarWithUncertainty(inputRasters[5],  filetag="_MHW")
#zStarWithUncertainty(inputRasters[6],  filetag="_MHW")
#zStarWithUncertainty(inputRasters[7],  filetag="_MHW")
#zStarWithUncertainty(inputRasters[9],  filetag="_MHW", hydroFlatteningValue=0, mapPartialDerivatives=True)
#zStarWithUncertainty(inputRasters[10],  filetag="_MHW", hydroFlatteningValue=0)
#zStarWithUncertainty(inputRasters[11],  filetag="_MHW", hydroFlatteningValue=0)
#zStarWithUncertainty(inputRasters[12],  filetag="_MHW", hydroFlatteningValue=0)
#zStarWithUncertainty(inputRasters[13],  filetag="_MHW", hydroFlatteningValue=-0.304801)
