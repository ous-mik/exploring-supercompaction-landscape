import numpy as np
import pandas as pd
import os
import time
import dirichletintegrate
import polyafit
import DNACompactionAnalysis as DCA

"""GO TO THE LAST LINES OF THIS DOCUMENT TO ADD INFORMATION FOR SCRIPT"""

def getListOfFilesInPlateRow(directorypath, plate_date, row_letter):
    """Create a list of filenames for all files in input directory."""
    absDirectorypath = os.path.abspath(directorypath)
    listOfFiles = os.listdir(absDirectorypath)
    allFiles = list()
    "Check if object in list is directory or not"
    for entry in listOfFiles:
        "Create full path"
        fullPath = os.path.join(absDirectorypath, entry)
        "If object is a directory then get the list of files in this directory"
        if os.path.isdir(fullPath):
            if len(os.listdir(fullPath)) == 0:
                pass
            else:
                allFiles = allFiles + getListOfFilesInPlateRow(fullPath, plate_date, row_letter)
        else:
            "Check if entry belongs to specified row of specified plate"
            entryInfo = entry.split("_")
            if (entryInfo[0] == plate_date) and (entryInfo[1][0] == row_letter):
                allFiles.append(fullPath)
            else:
                continue    
    return allFiles


def combineRowResults(directorypath, plate_date, row_letter, absParentOutputDirectorypath, cluster):

    fileList = getListOfFilesInPlateRow(directorypath, plate_date, row_letter)
    nonprocessedFiles = 0
    imageMetadataExists = False
    interestingObjectsMeasurementsExists = False
    objectMeasurementsExists = False

    for eachFile in fileList:
        if cluster:
            filePathParts = eachFile.split("/")
        else:
            filePathParts = eachFile.split("\\")
        fileInfo = filePathParts[-1].split("_")
        fileContent = None
        
        if (fileInfo[2] == "Image") or (fileInfo[2] == "ImageMetadata.csv"):
            if imageMetadataExists:
                fileContent = pd.read_csv(eachFile)
                frames = [imageMetadata,fileContent]
                imageMetadata = pd.concat(frames, axis=0, ignore_index=True)
            else:
                imageMetadata = pd.read_csv(eachFile)
                imageMetadataExists = True
        elif fileInfo[2] == "InterestingObjectsMeasurements.csv":
            if interestingObjectsMeasurementsExists:
                fileContent = pd.read_csv(eachFile, header=[0, 1])
                frames = [interestingObjectsMeasurements,fileContent]
                interestingObjectsMeasurements = pd.concat(frames, axis=0, ignore_index=True)
            else:
                interestingObjectsMeasurements = pd.read_csv(eachFile, header=[0, 1])
                interestingObjectsMeasurementsExists = True
        elif fileInfo[2] == "ObjectMeasurement.csv":
            if objectMeasurementsExists:
                fileContent = pd.read_csv(eachFile, header=[0, 1])
                frames = [objectMeasurements,fileContent]
                objectMeasurements = pd.concat(frames, axis=0, ignore_index=True)
            else:
                objectMeasurements = pd.read_csv(eachFile, header=[0, 1])
                objectMeasurementsExists = True
        else:
            nonprocessedFiles = nonprocessedFiles + 1
            print(eachFile)
            print(fileInfo)

    if nonprocessedFiles:
        print("Finished processing row " + row_letter + " in plate " + plate_date + ".\nNumber of files that were not added to the combined files: " + str(nonprocessedFiles))
    else:
        print("Finished processing row " + row_letter + " in plate " + plate_date + ".\nAll files added to the combined files")


    # row_directory = str(plate_date) + "_" + row_letter + "_row"
    # output_directory = os.path.join(absParentOutputDirectorypath, row_directory)
    output_directory = absParentOutputDirectorypath  #Check what is the best way to do this
    
    output_filename_prefix = str(plate_date) + "_" + str(row_letter) + "_RowRes_"
    filename_imageMetadata = output_filename_prefix + "ImageMetadata.csv"
    outputPath_imageMetadata = os.path.join(output_directory,filename_imageMetadata)
    filename_interestingObjects = output_filename_prefix + "InterestingObjectsMeasurements.csv"
    outputPath_interestingObjects = os.path.join(output_directory,filename_interestingObjects)
    filename_objectMeasurements = output_filename_prefix + "ObjectMeasurements.csv"
    outputPath_objectMeasurements = os.path.join(output_directory,filename_objectMeasurements)

    imageMetadata.to_csv(outputPath_imageMetadata, index=False)
    interestingObjectsMeasurements.to_csv(outputPath_interestingObjects, index=False)
    objectMeasurements.to_csv(outputPath_objectMeasurements, index=False)

    print("\nExported combined files for Row " + row_letter + " to location: " + output_directory + "\nFilenames are: \n" + filename_imageMetadata + "\n" + filename_interestingObjects + "\n" + filename_objectMeasurements)

    #return filename_imageMetadata, filename_interestingObjects, filename_objectMeasurements


def analyzeRowResults(directorypath, plate_date, row_letter, absParentOutputDirectorypath, cluster):
    
    """USE FUNCTION AnalyzeDNAinBacteria()"""
    """ENTER CORRECT INPUT VARIABLES (PREFIX OF INPUT FILENAME, OUTPUT FILENAME, OUTPUT FOLDER PATH, cluster=False)"""
    """EXAMPLE: AnalyzeDNAinBacteria("200501_F_RowRes_", "200501_F_RowAnalyzed.xlsx", "./AnalyzedPlates/200501_Plate5/output-data", cluster=False)"""
    """MAKE SURE ALL OF YOUR INPUT FILES ARE IN THE SAME FOLDER (PATH)"""
    """REMEMBER TO INCLUDE WellGeneIndexList.xlsx IN THE SAME FOLDER (PATH)"""

    combineRowResults(directorypath, plate_date, row_letter, absParentOutputDirectorypath, cluster)

    filename_prefix = str(plate_date) + "_" + str(row_letter) + "_RowRes_"
    output_filename = str(plate_date) + "_" + str(row_letter) + "_RowAnalyzed.xlsx"

    DCA.AnalyzeDNAinBacteria(filename_prefix,output_filename,absParentOutputDirectorypath, cluster)


def combinePlateResults(parent_directorypath, plate_date, cluster=False):
    """Will combine all calculation results from CellProfiler and combine them for each row
    before performing analyses of the results. The information from the analyses are then 
    combined into one file. 

    The parent_directorypath should be a default output folder in the format 
    {path}/{plate_date}/output-data."""

    row_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    used_row_letters = []

    ##Find all directories inside parent directory that match some format and make list
    absParentOutputDirectorypath = os.path.abspath(parent_directorypath)
    directoryList = os.listdir(absParentOutputDirectorypath)
    
    for eachDirectorypath in directoryList:
        fullPath = os.path.join(absParentOutputDirectorypath, eachDirectorypath)
        if os.path.isdir(fullPath):
            dirNameInfo = eachDirectorypath.split("_")
            if dirNameInfo[0] == plate_date:
                row_letter = dirNameInfo[1]
                if (row_letter in row_letters) and not (row_letter in used_row_letters):
                    analyzeRowResults(fullPath, plate_date, row_letter, absParentOutputDirectorypath, cluster)
                    used_row_letters.append(row_letter)
            else:
                continue
        else:
            continue
    
    #Combine output from all row files just created using dataframes and sort by well number
    # cwd = os.getcwd() # Could maybe reuse earlier code and specify folder for output of files
    fileList = os.listdir(absParentOutputDirectorypath)
    combinedPlateResultsExists = False

    for eachFile in fileList:
        filePath = os.path.join(absParentOutputDirectorypath, eachFile)
        if os.path.isfile(filePath):
            fileInfo = eachFile.split("_")
            if len(fileInfo) >= 3:
                if fileInfo[2] == "RowAnalyzed.xlsx":
                    if combinedPlateResultsExists:
                        fileContent = pd.read_excel(filePath, sheet_name="Combined_well_results", index_col=0)
                        frames = [combinedPlateResults,fileContent]
                        combinedPlateResults = pd.concat(frames, axis=0, ignore_index=True)
                    else: 
                        combinedPlateResults = pd.read_excel(filePath, sheet_name="Combined_well_results", index_col=0)
                        combinedPlateResultsExists = True
                else:
                    continue
            else:
                continue
        else:
            continue
    
    combinedPlateResults.sort_values(by=['Plate_date', 'Well_no'])
    output_filename = str(plate_date) + "_FullPlateAnalyzed.xlsx"
    output_filepath = os.path.join(absParentOutputDirectorypath, output_filename)

    with pd.ExcelWriter(output_filepath) as writer: # pylint: disable=abstract-class-instantiated
        combinedPlateResults.to_excel(writer, sheet_name='Combined_well_results')

    print("\nFinished analyses of the results in all wells from plate " + plate_date + "\n Combined files of all data from each row after image analysis can be found in location: " + parent_directorypath + "\nA file with results of the well analyses (" + output_filename + ") can be found in the same location.\n")


################  FILL IN INFORMATION FOR YOUR RUN BELOW THIS LINE  ################
"""USE FUNCTION combinePlateResults()"""
"""ENTER CORRECT INPUT VARIABLE (OUTPUT-DATA FOLDER PATH, PLATE DATE/NAME, cluster=False)"""
"""EXAMPLE: combinePlateResults("./AnalyzedPlates/200908_Plate1", "200908", cluster=False)"""
"""MAKE SURE ALL OF YOUR INPUT FILES ARE IN THE SAME FOLDER (PATH)"""
"""REMEMBER TO INCLUDE DNACompactionAnalysis.py, dirichletintegrate.py, polyafit.py, hypergeom.py, PlateDateNumberIndexList.xlsx and WellGeneIndexList.xlsx IN THE SAME FOLDER (PATH)"""

combinePlateResults("<OUTPUT-DATA FOLDER PATH>", "<PLATE DATE/NAME>", cluster=False)