''' ---------------------------------------------------------------------------
Name: TCEQ Gridded PMP Tool Python Script

Initial AWA Script Version: 1.6a
TCEQ GeoProcessing Script Version: 1.1

Python Version: 2.7

ArcGIS Version: ArcGIS Desktop 10.2

Initial Script Author: Applied Weather Associates
GeoProcessing Script Author: TCEQ Dam Safety / Morgan Dean

Usage:  The tool is designed to be executed within an the ArcMap environment
with an open MXD session.  The tool has been altered for integration through a
web app widget for users without ArcGIS.

Description:
    This tool calculates PMP depths for a given drainage basin for the
specified durations.  PMP point values are calculated (in inches) for each
grid point (spaced at 90 arc-second intervals) over the project domain. The
points are converted to gridded PMP datasets for each duration.

See TCEQ Texas Statewide PMP Study for more information.
---------------------------------------------------------------------------'''

###########################################################################
## import Python modules

import sys
import arcpy
from arcpy import env
import arcpy.analysis as an
import arcpy.management as dm
import arcpy.conversion as con
##import numpy


env.overwriteOutput = True                                                      # Set overwrite option
env.addOutputsToMap = False

###########################################################################
## get input parameters

basin = arcpy.GetParameter(0)                                                   # get AOI Basin Shapefile
home = arcpy.GetParameterAsText(1)                                              # get location of 'PMP' Project Folder
locDurations = arcpy.GetParameter(2)                                            # get local storm durations (string)
genDurations = arcpy.GetParameter(3)                                            # get general storm durations (string)
tropDurations = arcpy.GetParameter(4)                                           # get tropical storm durations (string)
weightedAve = arcpy.GetParameter(7)                                             # get option to apply weighted average (boolean)
outputLoc = arcpy.GetParameter(8)                                               # get file path for Local Storm PMP Points
outputGen = arcpy.GetParameter(9)                                               # get file path for General Storm PMP Points
outputTrop = arcpy.GetParameter(10)                                             # get file path for Tropical Storm PMP Points
outputBAT = arcpy.GetParameter(11)                                              # get file path for Basin Average Summary Table


dadGDB = home + "\\Input\\DAD_Tables.gdb"                                       # location of DAD tables
adjFactGDB = home + "\\Input\\Storm_Adj_Factors.gdb"                            # location of feature datasets containing total adjustment factors
arcpy.AddMessage("\nDAD Tables geodatabase path:  " + dadGDB)
arcpy.AddMessage("Storm Adjustment Factor geodatabase path:  " + adjFactGDB)

basAveTables = []                                                               # global list of Basin Average Summary tables

def pmpAnalysis(aoiBasin, stormType, durList):

    ###########################################################################
    ## Create PMP Point Feature Class from points within AOI basin and add fields
    def createPMPfc():

        arcpy.AddMessage("\nCreating feature class: 'PMP_Points' in Scratch.gdb...")
        dm.MakeFeatureLayer(home + "\\Input\Non_Storm_Data.gdb\Vector_Grid", "vgLayer")                 # make a feature layer of vector grid cells
        dm.SelectLayerByLocation("vgLayer", "INTERSECT", aoiBasin)                                       # select the vector grid cells that intersect the aoiBasin polygon
        dm.MakeFeatureLayer(home + "\\Input\Non_Storm_Data.gdb\Grid_Points", "gpLayer")                 # make a feature layer of grid points
        dm.SelectLayerByLocation("gpLayer", "HAVE_THEIR_CENTER_IN", "vgLayer")                           # select the grid points within the vector grid selection
        con.FeatureClassToFeatureClass("gpLayer", env.scratchGDB, "PMP_Points")                          # save feature layer as "PMP_Points" feature class
        arcpy.AddMessage("(" + str(dm.GetCount("gpLayer")) + " grid points will be analyzed)\n")

        # Add PMP Fields
        for dur in durList:
            arcpy.AddMessage("\t...adding field: PMP_" + str(dur))
            dm.AddField(env.scratchGDB + "\\PMP_Points", "PMP_" + dur, "DOUBLE")


        # Add STORM Fields (this string values identifies the driving storm by SPAS ID number)
        for dur in durList:
            arcpy.AddMessage("\t...adding field: STORM_" + str(dur))
            dm.AddField(env.scratchGDB + "\\PMP_Points", "STORM_" + dur, "TEXT", "", "", 16)

        return

    ###########################################################################
    ##  Define getAOIarea() function:
    ##  getAOIarea() calculates the area of AOI (basin outline) input shapefile/
    ##  featureclass.  The basin outline shapefile must be projected.  The area
    ##  is sqaure miles, converted from the basin layers projected units (feet
    ##  or meters).  The aoiBasin feature class should only have a single feature
    ##  (the basin outline).  If there are multiple features, the area will be stored
    ##  for the final feature only.

    def getAOIarea():
        sr = arcpy.Describe(aoiBasin).SpatialReference                                          # Determine aoiBasin spatial reference system
        srname = sr.name
        srtype = sr.type
        srunitname = sr.linearUnitName                                                          # Units
        arcpy.AddMessage("\nAOI basin spatial reference:  " + srname + "\nUnit type: " + srunitname + "\nSpatial reference type: " + srtype)

        aoiArea = 0.0
        rows = arcpy.SearchCursor(aoiBasin)
        for row in rows:
            feat = row.getValue("Shape")
            aoiArea += feat.area
        if srtype == 'Geographic':                                  # Must have a surface projection.  If one doesn't exist it projects a temporary file and uses that.
           arcpy.AddMessage("\n***The basin shapefile's spatial reference 'Geographic' is not supported.  Projecting temporary shapefile for AOI.***")
           arcpy.Project_management(aoiBasin, env.scratchGDB + "\\TempBasin", 102039)     #Projects AOI Basin (102039 = USA_Contiguous_Albers_Equal_Area_Conic_USGS_version)
           TempBasin = env.scratchGDB + "\\TempBasin"                       # Path to temporary basin created in scratch geodatabase
           sr = arcpy.Describe(TempBasin).SpatialReference                  # Determine Spatial Reference of temporary basin
           aoiArea = 0.0
           rows = arcpy.SearchCursor(TempBasin)                     # Assign area size in square meters
           for row in rows:
            feat = row.getValue("Shape")
            aoiArea += feat.area
           aoiArea = aoiArea * 0.000000386102                                           # Converts square meters to square miles
        elif srtype == 'Projected':                                  # If a projection exists, it re-projects a temporary file and uses that for data consistency.
            arcpy.AddMessage("\n***The basin shapefile's spatial reference will be reprojected to USA_Contiguous_Albers_Equal_Area_Conic_USGS_version for data consistency.  Projecting temporary shapefile for AOI.***")
            arcpy.Project_management(aoiBasin, env.scratchGDB + "\\TempBasin", 102039)     #Projects AOI Basin (102039 = USA_Contiguous_Albers_Equal_Area_Conic_USGS_version)
            TempBasin = env.scratchGDB + "\\TempBasin"                       # Path to temporary basin created in scratch geodatabase
            sr = arcpy.Describe(TempBasin).SpatialReference                  # Determine Spatial Reference of temporary basin
            aoiArea = 0.0
            rows = arcpy.SearchCursor(TempBasin)                     # Assign area size in square meters
            for row in rows:
             feat = row.getValue("Shape")
             aoiArea += feat.area
            aoiArea = aoiArea * 0.000000386102                                           # Converts square meters to square miles

        aoiArea = round(aoiArea, 3)
        arcpy.AddMessage("\nArea of interest: " + str(aoiArea) + " square miles.")

        if arcpy.GetParameter(5) == False:
            aoiArea = arcpy.GetParameter(6)                                                     # Enable a constant area size
        aoiArea = round(aoiArea, 1)
        arcpy.AddMessage("\n***Area used for PMP analysis: " + str(aoiArea) + " sqmi***")
        return aoiArea

    ###########################################################################
    ##  Define dadLookup() function:
    ##  The dadLookup() function determines the DAD value for the current storm
    ##  and duration according to the basin area size.  The DAD depth is interpolated
    ##  linearly between the two nearest areal values within the DAD table.
    def dadLookup(stormLayer, duration, area):                  # dadLookup() accepts the current storm layer name (string), the current duration (string), and AOI area size (float)
        #arcpy.AddMessage("\t\tfunction dadLookup() called.")
        durField = "H_" + duration                              # defines the name of the duration field (eg., "H_06" for 6-hour)
        dadTable = dadGDB + "\\" + stormLayer
        rows = arcpy.SearchCursor(dadTable)

        try:
            row = rows.next()                                       # Sets DAD area x1 to the value in the first row of the DAD table.
            x1 = row.AREASQMI
            y1 = row.getValue(durField)
            xFlag = "FALSE"                                         # xFlag will remain false for basins that are larger than the largest DAD area.
        except RuntimeError:                                        # return if duration does not exist in DAD table
            return

        row = rows.next()
        i = 0
        while row:                                                  # iterates through the DAD table - assiging the bounding values directly above and below the basin area size
            i += 1
            if row.AREASQMI < area:
                x1 = row.AREASQMI
                y1 = row.getValue(durField)
            else:
                xFlag = "TRUE"                                      # xFlag is switched to "TRUE" indicating area is within DAD range
                x2 = row.AREASQMI
                y2 = row.getValue(durField)
                break

            row = rows.next()
        del row, rows, i

        if xFlag == "FALSE":
            x2 = area                                           # If x2 is equal to the basin area, this means that the largest DAD area is smaller than the basin and the resulting DAD value must be extrapolated.
            arcpy.AddMessage("\t\tThe basin area size: " + str(area) + " sqmi is greater than the largest DAD area: " + str(x1) + " sqmi.\n\t\tDAD value is estimated by extrapolation.")
            y = x1 / x2 * y1                                    # y (the DAD depth) is estimated by extrapolating the DAD area to the basin area size.
            return y                                            # The extrapolated DAD depth (in inches) is returned.

        # arcpy.AddMessage("\nArea = " + str(area) + "\nx1 = " + str(x1) + "\nx2 = " + str(x2) + "\ny1 = " + str(y1) + "\ny2 = " + str(y2))

        x = area                                                # If the basin area size is within the DAD table area range, the DAD depth is interpolated
        deltax = x2 - x1                                        # to determine the DAD value (y) at area (x) based on next lower (x1) and next higher (x2) areas.
        deltay = y2 - y1
        diffx = x - x1

        y = y1 + diffx * deltay / deltax

        if x < x1:
            arcpy.AddMessage("\t\tThe basin area size: " + str(area) + " sqmi is less than the smallest DAD table area: " + str(x1) + " sqmi.\n\t\tDAD value is estimated by extrapolation.")

        return y                                                # The interpolated DAD depth (in inches) is returned.

    ###########################################################################
    ##  Define updatePMP() function:
    ##  This function updates the 'PMP_XX_' and 'STORM_XX' fields of the PMP_Points
    ##  feature class with the largest value from all analyzed storms stored in the
    ##  pmpValues list.
    def updatePMP(pmpValues, stormID, duration):                                                    # Accepts four arguments: pmpValues - largest adjusted rainfall for current duration (float list); stormID - driver storm ID for each PMP value (text list); and duration (string)
        pmpfield = "PMP_" + duration
        stormfield = "STORM_" + duration
        gridRows = arcpy.UpdateCursor(env.scratchGDB + "\\PMP_Points")                              # iterates through PMP_Points rows
        i = 0
        for row in gridRows:
            row.setValue(pmpfield, pmpValues[i])                                                    # Sets the PMP field value equal to the Max Adj. Rainfall value (if larger than existing value).
            row.setValue(stormfield, stormID[i])                                                    # Sets the storm ID field to indicate the driving storm event
            gridRows.updateRow(row)
            i += 1
        del row, gridRows, pmpfield, stormfield
        arcpy.AddMessage("\n\t" + duration + "-hour PMP values update complete. \n")
        return

    ###########################################################################
    ##  The outputPMP() function produces raster GRID files for each of the PMP durations.
    ##  Aslo, a space-delimited PMP_Distribition.txt file is created in the 'Text_Output' folder.
    def outputPMP(type, area, outPath):
        desc = arcpy.Describe(basin)
        basinName = desc.baseName
        pmpPoints = env.scratchGDB + "\\PMP_Points"                             # Location of 'PMP_Points' feature class which will provide data for output

        outType = type[:1]
        outArea = str(int(round(area,0))) + "sqmi"
        outFC = outType + "_" + outArea  #I don't think I need this.....
        arcpy.AddMessage("\nCopying PMP_Points feature class to " + outFC + "...")  #outFC might be replaced with outpath...
        dm.Merge(pmpPoints, outPath)                                            # merge the scratch feature layer(s) of vector grid cells into the outputs

        arcpy.AddMessage("\nCreating Basin Summary Table...")
        tableName = type + "_PMP_Basin_Average" + "_" + outArea
        tablePath = env.scratchGDB + "\\" + tableName
        dm.CreateTable(env.scratchGDB, tableName)                               # Create blank table
        cursor = arcpy.da.InsertCursor(tablePath, "*")                          # Create Insert cursor and add a blank row to the table
        cursor.insertRow([0])
        del cursor

        dm.AddField(tablePath, "STORM_TYPE", "TEXT", "", "", 10, "Storm Type")          # Create "Storm Type" field
        dm.CalculateField(tablePath, "STORM_TYPE", "'" + type + "'", "PYTHON_9.3")      # populate storm type field

        i = 0
        for field in arcpy.ListFields(pmpPoints, "PMP_*"):          # Add fields for each PMP duration and calculate the basin average
            fieldName = field.name
            fieldAve = basinAve(basin, fieldName)                   # Calls the basinAve() function - returns the average (weighted or not)
            dm.AddField(tablePath, fieldName, "DOUBLE", "", 2)       # Add duration field
            dm.CalculateField(tablePath, fieldName, fieldAve, "PYTHON_9.3")       # Assigns the basin average

            i += 1
        arcpy.AddMessage("\nSummary table complete.")

        basAveTables.append(tablePath)


        return


    ###########################################################################
    ##  The basin() returns the basin average PMP value for a given duration field.
    ##  If the option for a weighted average is checked in the tool parameter the script
    ##  will weight the grid point values based on proportion of area inside the basin.
    def basinAve(aoiBasin, pmpField):
        pmpPoints = env.scratchGDB + "\\PMP_Points"                                                         # Path of 'PMP_Points' scratch feature class
        if weightedAve:
            arcpy.AddMessage("\tCalculating basin average for " + pmpField + "(weighted)...")
            vectorGridClip = env.scratchGDB + "\\VectorGridClip"                                            # Path of 'PMP_Points' scratch feature class
            sumstats = env.scratchGDB + "\\SummaryStats"

            dm.MakeFeatureLayer(home + "\\Input\Non_Storm_Data.gdb\\Vector_Grid", "vgLayer")                # make a feature layer of vector grid cells
            dm.SelectLayerByLocation("vgLayer", "INTERSECT", aoiBasin)                                      # select the vector grid cells that intersect the aoiBasin polygon

            an.Clip("vgLayer", aoiBasin, vectorGridClip)                                                    # clips aoi vector grid to basin
            dm.AddField(pmpPoints, "WEIGHT", "DOUBLE")                                                      # adds 'WEIGHT' field to PMP_Points scratch feature class
            dm.MakeFeatureLayer(vectorGridClip, "vgClipLayer")                                              # make a feature layer of basin clipped vector grid cells
            dm.MakeFeatureLayer(pmpPoints, "pmpPointsLayer")                                                # make a feature layer of PMP_Points feature class

            dm.AddJoin("pmpPointsLayer", "ID", "vgClipLayer", "ID")                                         # joins PMP_Points and vectorGridBasin tables
            dm.CalculateField("pmpPointsLayer", "WEIGHT", "!vectorGridClip.Shape_Area!", "PYTHON_9.3")      # Calculates basin area proportion to use as weight for each grid cell.
            dm.RemoveJoin("pmpPointsLayer", "vectorGridClip")

            an.Statistics(pmpPoints, sumstats, [["WEIGHT", "SUM"]], "")
            stats = arcpy.SearchCursor(sumstats)
            pmpWgtAve = pmpField + "_WgtAve"

            for row in stats:
                calc = row.getValue("SUM_WEIGHT")
                express = "(!WEIGHT!/{})* !{}!".format(calc, pmpField)
                i = 0
                for field in arcpy.ListFields(pmpPoints, pmpField):
                    dm.AddField(pmpPoints, pmpWgtAve, "DOUBLE", 2)
                    dm.CalculateField(pmpPoints, pmpWgtAve, express, "PYTHON_9.3")
                    i += 1
                del stats, row

            an.Statistics(pmpPoints, sumstats, [[pmpWgtAve, "SUM"]], "")
            sumwgtave = "SUM_" + pmpWgtAve
            with arcpy.da.SearchCursor(sumstats, sumwgtave) as stats:
                for row in stats:
                    wgtAve = row [0]
                    return round(wgtAve, 2)


##            na = arcpy.da.TableToNumPyArray(pmpPoints,(pmpField, 'WEIGHT'))                                 # Assign pmpPoints values and weights to Numpy array (na)
##            wgtAve = numpy.average(na[pmpField], weights=na['WEIGHT'])                                         # Calculate weighted average with Numpy average
##            del na
##            return round(wgtAve, 2)

        else:
            arcpy.AddMessage("\tCalculating basin average for " + pmpField + "(not weighted)...")
            sumstats = env.scratchGDB + "\\SummaryStats"
            an.Statistics(pmpPoints, sumstats, [[pmpField, "MEAN"]], "")
            mean = "MEAN_" + pmpField
            with arcpy.da.SearchCursor(sumstats, mean) as stats:
                for row in stats:
                    fieldAve = row [0]
                    return round(fieldAve, 2)

##            na = arcpy.da.TableToNumPyArray(pmpPoints, pmpField)                                            # Assign pmpPoints values to Numpy array (na)
##            fieldAve = numpy.average(na[pmpField])                                                             # Calculates aritmetic mean
##            del na
##            return round(fieldAve, 2)


    ###########################################################################
    ##  This portion of the code iterates through each storm feature class in the
    ##  'Storm_Adj_Factors' geodatabase (evaluating the feature class only within
    ##  the Local, Tropical, or general feature dataset).  For each duration,
    ##  at each grid point within the aoi basin, the transpositionality is
    ##  confirmed.  Then the DAD precip depth is retrieved and applied to the
    ##  total adjustement factor to yield the total adjusted rainfall.  This
    ##  value is then sent to the updatePMP() function to update the 'PMP_Points'
    ##  feature class.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

    desc = arcpy.Describe(basin)                                                        # Check to ensure AOI input shape is a Polygon. If not - exit.
    basinShape = desc.shapeType
    if desc.shapeType == "Polygon":
        arcpy.AddMessage("\nBasin shape type: " + desc.shapeType)
    else:
        arcpy.AddMessage("\nBasin shape type: " + desc.shapeType)
        arcpy.AddMessage("\nError: Input shapefile must be a polygon!\n")
        sys.exit()

    createPMPfc()                                                                       # Call the createPMPfc() function to create the PMP_Points feature class.

    env.workspace = adjFactGDB                                                          # the workspace environment is set to the 'Storm_Adj_Factors' file geodatabase

    aoiSQMI = round(getAOIarea(),2)                                                     # Calls the getAOIarea() function to assign area of AOI shapefile to 'aoiSQMI'

    for dur in durList:
        stormList = arcpy.ListFeatureClasses("", "Point", stormType)                    # List all the total adjustment factor feature classes within the storm type feature dataset.

        arcpy.AddMessage("\n*************************************************************\nEvaluating " + dur + "-hour duration...")

        pmpList = []
        driverList = []
        gridRows = arcpy.SearchCursor(env.scratchGDB + "\\PMP_Points")
        try:
            for row in gridRows:
                pmpList.append(0.0)                                                         # creates pmpList of empty float values for each grid point to store final PMP values
                driverList.append("STORM")                                                  # creates driverList of empty text values for each grid point to store final Driver Storm IDs
            del row, gridRows
        except UnboundLocalError:
            arcpy.AddMessage("\n***Error: No data present within basin/AOI area.***\n")
            sys.exit()

        for storm in stormList:
            arcpy.AddMessage("\n\tEvaluating storm: " + storm + "...")
            dm.MakeFeatureLayer(storm, "stormLayer")                                    # creates a feature layer for the current storm
            dm.SelectLayerByLocation("stormLayer", "HAVE_THEIR_CENTER_IN", "vgLayer")   # examines only the grid points that lie within the AOI
            gridRows = arcpy.SearchCursor("stormLayer")
            pmpField = "PMP_" + dur
            i = 0
            try:
                dadPrecip = round(dadLookup(storm, dur, aoiSQMI),3)
                arcpy.AddMessage("\t\t" + dur + "-hour DAD value:  " + str(dadPrecip) + chr(34))
            except TypeError:                                                           # In no duration exists in the DAD table - move to the next storm
                arcpy.AddMessage("\t***Duration '" + str(dur) + "-hour' is not present for " + str(storm) + ".***\n")
                continue
            arcpy.AddMessage("\t\tComparing " + storm + " adjusted rainfall values against current driver values...\n")
            for row in gridRows:
                if row.TRANS == 1:                                              # Only continue if grid point is transpositionable ('1' is transpostionable, '0' is not).
                    try:                                                        # get total adj. factor if duration exists
                        adjRain = round(dadPrecip * row.TAF,1)
                        if adjRain > pmpList[i]:
                            pmpList[i] = adjRain
                            driverList[i] = storm
                    except RuntimeError:
                        arcpy.AddMessage("\t\t   *Warning*  Total Adjusted Raifnall value falied to set for row " + str(row.CNT))
                        break
                    del adjRain
                i += 1
            del row
        del storm, stormList, gridRows, dadPrecip
        updatePMP(pmpList, driverList, dur)              # calls function to update "PMP Points" feature class
    del dur, pmpList

    arcpy.AddMessage("\n'PMP_Points' Feature Class 'PMP_XX' fields update complete for all '" + stormType + "' storms.")

    outputPMP(stormType, aoiSQMI, outPath)               # calls outputPMP() function

    del aoiSQMI
    return
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

def outputBasAveTable():
    arcpy.AddMessage("\nCreating basin average summary table.\n")
    tableList = basAveTables
    for table in tableList:
        arcpy.AddMessage("\t\t...Table: " + table)
    dm.Merge(basAveTables, outputBAT)                                           # Merges scratch basin average tables into the final table output

    return


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


if locDurations:
    type = "Local"
    durations = locDurations
    outPath = outputLoc
    arcpy.AddMessage("\nRunning PMP analysis for storm type: " + type)
    pmpAnalysis(basin, type, durations)          # Calls the pmpAnalysis() function to calculate the local storm PMP
    arcpy.AddMessage("\nLocal storm analysis complete...\n*********************************************************************************************************")

if genDurations:
    type = "General"
    durations = genDurations
    outPath = outputGen
    arcpy.AddMessage("\nRunning PMP analysis for storm type: " + type)
    pmpAnalysis(basin, type, durations)          # Calls the pmpAnalysis() function to calculate the general storm PMP
    arcpy.AddMessage("\nGeneral Winter storm analysis complete...\n*********************************************************************************************************")

if tropDurations:
    type = "Tropical"
    durations = tropDurations
    outPath = outputTrop
    arcpy.AddMessage("\nRunning PMP analysis for storm type: " + type)
    pmpAnalysis(basin, type, durations)          # Calls the pmpAnalysis() function to calculate the tropical storm PMP
    arcpy.AddMessage("\nTropical storm analysis complete...\n*********************************************************************************************************")

if outputBAT:
    arcpy.AddMessage("\nVariable outputTable exists ")
    outPath = outputBAT
    outputBasAveTable()                         # Calls the outputBasAveTable() function to merge the weighted basin averages into a single table output

