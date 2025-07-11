import numpy as np
import pandas as pd
import math
import time
import os
import dirichletintegrate
import polyafit

# Enrichment score calculations in this script is based on the code from CellProfiler Analyst,
# Copyright (c) Broad Institute
# Licensed under the BSD 3-Clause License. 

"""GO TO THE LAST LINES OF THIS DOCUMENT TO ADD INFORMATION FOR CALCULATIONS"""

def AnalyzeDNAinBacteria(input_filename_prefix,output_filename,absParentDirectorypath,cluster):
    
    ################  POTENTIAL ERROR MESSAGES  ################

    """Check if file with output filename already exists in path"""
    output_path = os.path.join(absParentDirectorypath, output_filename)
    if os.path.exists(output_path):
        raise Exception("A file with the name " + output_filename + " already exists in path. Enter a different output_filename.")
    else: 
        pass

    """Check if important python files are in path"""

    if cluster:
        main_input_file_path = "/cluster/projects/nn9383k/kristerv/cp-mainruns"
    else:
        main_input_file_path = "."

    polyafit_path = os.path.join(main_input_file_path, "polyafit.py")
    dirichletintegrate_path = os.path.join(main_input_file_path, "dirichletintegrate.py")
    hypergeom_path = os.path.join(main_input_file_path, "hypergeom.py")

    if not os.path.exists(polyafit_path):
        raise Exception("The file polyafit.py is missing in path. Please add it to path, it is necessary to run the code.")
    elif not os.path.exists(dirichletintegrate_path):
        raise Exception("The file dirichletintegrate.py is missing in path. Please add it to path, it is necessary to run the code.")
    elif not os.path.exists(hypergeom_path):
        raise Exception("The file hypergeom.py is missing in path. Please add it to path, it is necessary to run the code.")
    else:
        pass

    wellGene_path = os.path.join(main_input_file_path, "WellGeneIndexList.xlsx")
    plateDate_path = os.path.join(main_input_file_path, "PlateDateNumberIndexList.xlsx")

    """Check if important excel files are in path"""
    if not os.path.exists(wellGene_path):
        raise Exception("The file WellGeneIndexList.xlsx is missing in path. Please add it to path, it is necessary to run the code.")
    elif not os.path.exists(plateDate_path):
        raise Exception("The file PlateDateNumberIndexList.xlsx is missing in path. Please add it to path, it is necessary to run the code.")
    else:
        pass

    filename_object_results = input_filename_prefix + "ObjectMeasurements.csv"
    filename_image_metadata = input_filename_prefix + "ImageMetadata.csv"

    filepath_object_results = os.path.join(absParentDirectorypath, filename_object_results)
    filepath_image_metadata = os.path.join(absParentDirectorypath, filename_image_metadata)

    # if not os.path.exists(filename_object_results):
    #     raise Exception("The file " + filename_object_results + " is missing in path. Please check that the input filename prefix is correct and that it is added to path.")
    # elif not os.path.exists(filename_image_metadata):
    #     raise Exception("The file " + filename_image_metadata + " is missing in path. Please check that the input filename prefix is correct and that it is added to path.")
    # else: 
    #     pass

    

    ################  LOAD/READ INPUT FILES  ################

    resultsObjects = pd.read_csv(filepath_object_results, header=[0, 1])
    imageCalculationsData = pd.read_csv(filepath_image_metadata, header=[0, 0])

    wellGeneIndexList = pd.read_excel(wellGene_path, header=[0, 1])
    plateDateNumberIndexList = pd.read_excel(plateDate_path, header=[0, 0])



    ################  EXTRACTION OF METADATA FOR CLASSIFICATION OF IMAGES  ################

    imageNumber_metadata = imageCalculationsData.loc[:,"ImageNumber"]
    imageNumber_metadata = imageNumber_metadata["ImageNumber"].tolist() 

    well_number_metadata = imageCalculationsData.loc[:,"Metadata_Well"]
    well_number_metadata = well_number_metadata["Metadata_Well"].tolist() 

    site_number_metadata = imageCalculationsData.loc[:,"Metadata_Site"]
    site_number_metadata = site_number_metadata["Metadata_Site"].tolist() 

    plate_date_metadata = imageCalculationsData.loc[:,"Metadata_Plate"]
    plate_date_metadata = plate_date_metadata["Metadata_Plate"].tolist()

    for dates in range(len(plate_date_metadata)):
        plate_date_metadata[dates] = str(plate_date_metadata[dates])

    """Lists of object counts for all bacteria"""
    cellCount_AllInterestingBacteria = imageCalculationsData.loc[:, "Count_InterestingBacteriaObjects"]
    cellCount_AllInterestingBacteria = cellCount_AllInterestingBacteria["Count_InterestingBacteriaObjects"].tolist()

    ################  EXTRACTION OF MEASUREMENTS FOR OBJECTS  ################

    """List with image numbers for objects"""
    imageNumber_forObjects = resultsObjects.loc[:, ("Image", "ImageNumber")]
    imageNumber_forObjects = pd.Series.tolist(imageNumber_forObjects)
    
    """Include all categories that the bacteria have been sorted into"""
    bacteriaCategory = ["PositiveTreated","Nontreated","NegativeTreated"]
    two_classes = len(bacteriaCategory) == 2
    

    """Lists of object numbers for objects in each image"""
    objectNumber_forObjects = resultsObjects.loc[:, (bacteriaCategory[0] + "Bacteria", "ObjectNumber")]
    objectNumber_forObjects = pd.Series.tolist(objectNumber_forObjects)


    # cellCount_NegativeTreatedBacteria = imageCalculationsData.loc[:, "Count_NegativeTreatedBacteria"]
    # cellCount_NegativeTreatedBacteria = cellCount_NegativeTreatedBacteria["Count_NegativeTreatedBacteria"].tolist()

    # cellCount_NontreatedBacteria = imageCalculationsData.loc[:, "Count_NontreatedBacteria"]
    # cellCount_NontreatedBacteria = cellCount_NontreatedBacteria["Count_NontreatedBacteria"].tolist()

    # cellCount_PositiveTreatedBacteria = imageCalculationsData.loc[:, "Count_PositiveTreatedBacteria"]
    # cellCount_PositiveTreatedBacteria = cellCount_PositiveTreatedBacteria["Count_PositiveTreatedBacteria"].tolist()

    # cellCount = {}
    cellCount = []
    bacteriaObjectIndex = []
    DNAObjectIndex = []
    parentBacteriaIndex = []
    DNAfociArea = []
    cellArea = []
    DNAcompactness = []
    DNAeccentricity = []
    DNAfociCountperBacteria = []
    cellularMeanIntensities = []
    DNAMeanIntensities = []
    cellularDNAmassDisplacement = []
    bacteriaMajorAxisLength = []
    bacteriaMinorAxisLength = []
    bacteriaMaxFeretLength = []
    bacteriaMinFeretLength = []

    distributionTotDNALevel1 = []
    distributionTotDNALevel2 = []
    distributionTotDNALevel3 = []
    distributionTotDNALevel4 = []
    distributionTotDNALevel5 = []
    distributionTotDNALevel6 = []
    distributionTotDNALevelExtra = []

    distributionConcDNALevel1 = []
    distributionConcDNALevel2 = []
    distributionConcDNALevel3 = []
    distributionConcDNALevel4 = []
    distributionConcDNALevel5 = []
    distributionConcDNALevel6 = []
    distributionConcDNALevelExtra = []


    columnLoadingVariable = None


    """Load all calculation values for every bacteria category"""
    for CellCategory in bacteriaCategory:
        # cellCount[bacteriaCategory.index(CellCategory)] = imageCalculationsData.loc[:, "Count_" + CellCategory + "Bacteria"]
        # cellCount[bacteriaCategory.index(CellCategory)] = cellCount[bacteriaCategory.index(CellCategory)]["Count_" + CellCategory + "Bacteria"].tolist()

        """Lists of object counts for different bacteria categories"""
        columnLoadingVariable = imageCalculationsData.loc[:, "Count_" + CellCategory + "Bacteria"]
        cellCount.append(columnLoadingVariable["Count_" + CellCategory + "Bacteria"].tolist())
        columnLoadingVariable = None

        """List with object number"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "Number_Object_Number")]
        bacteriaObjectIndex.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List with DNA object number"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "DNA", "Number_Object_Number")]
        DNAObjectIndex.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List for relating DNA object number to cell object number"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "DNA", "Parent_" + CellCategory + "Bacteria")]
        parentBacteriaIndex.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of DNA foci areas"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "DNA", "AreaShape_Area")]
        DNAfociArea.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of bacteria areas"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "AreaShape_Area")]
        cellArea.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of DNA compactness"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "DNA", "AreaShape_Compactness")]
        DNAcompactness.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of DNA eccentricity"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "DNA", "AreaShape_Eccentricity")]
        DNAeccentricity.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of DNA focis per bacteria"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "Children_" + CellCategory + "DNA_Count")]
        DNAfociCountperBacteria.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of average intensities inside cell objects"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "Intensity_MeanIntensity_InterestingDNAImage")]
        cellularMeanIntensities.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of average intensities inside DNA objects"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "DNA", "Intensity_MeanIntensity_InterestingDNAImage")]
        DNAMeanIntensities.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of mass displacement of DNA intensities in cells"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "Intensity_MassDisplacement_InterestingDNAImage")]
        cellularDNAmassDisplacement.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of bacteria length from Major Axis"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "AreaShape_MajorAxisLength")]
        bacteriaMajorAxisLength.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of bacteria width from Minor Axis"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "AreaShape_MinorAxisLength")]
        bacteriaMinorAxisLength.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of bacteria length from Max Feret"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "AreaShape_MaxFeretDiameter")]
        bacteriaMaxFeretLength.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """List of bacteria width from Min Feret"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "AreaShape_MinFeretDiameter")]
        bacteriaMinFeretLength.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """Lists of fractional total DNA intensities from DNA centers inside cells"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_FracAtD_EnhancedDNAImage_1of6")]
        distributionTotDNALevel1.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_FracAtD_EnhancedDNAImage_2of6")]
        distributionTotDNALevel2.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_FracAtD_EnhancedDNAImage_3of6")]
        distributionTotDNALevel3.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_FracAtD_EnhancedDNAImage_4of6")]
        distributionTotDNALevel4.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_FracAtD_EnhancedDNAImage_5of6")]
        distributionTotDNALevel5.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_FracAtD_EnhancedDNAImage_6of6")]
        distributionTotDNALevel6.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_FracAtD_EnhancedDNAImage_Overflow")]
        distributionTotDNALevelExtra.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        """Lists of fractional DNA 'concentration' (tot DNA normalized by no. of pixels)"""
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_MeanFrac_EnhancedDNAImage_1of6")]
        distributionConcDNALevel1.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_MeanFrac_EnhancedDNAImage_2of6")]
        distributionConcDNALevel2.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_MeanFrac_EnhancedDNAImage_3of6")]
        distributionConcDNALevel3.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_MeanFrac_EnhancedDNAImage_4of6")]
        distributionConcDNALevel4.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_MeanFrac_EnhancedDNAImage_5of6")]
        distributionConcDNALevel5.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_MeanFrac_EnhancedDNAImage_6of6")]
        distributionConcDNALevel6.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        
        columnLoadingVariable = resultsObjects.loc[:, (CellCategory + "Bacteria", "RadialDistribution_MeanFrac_EnhancedDNAImage_Overflow")]
        distributionConcDNALevelExtra.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None
        


    ################  FUNCTION DEFINITIONS  ################

    def find_start_end(lst):
        """This function applies to cases where not all rows are representative of the objects of interest (find correct interval)"""
        start, end = [0], []
        k = 0
        l = 0
        
        while  cellCount_AllInterestingBacteria[l] == 0:
            end.append(k)
            start.append(k)
            l = l + 1

        for bacteria1, bacteria2 in zip(lst, lst[1:]):
            if (k == 0 and bacteria1 == bacteria2) and not (imageNumber_forObjects[k + 1] == imageNumber_forObjects[k]):
                end.append(k)
                start.append(k + 1)
            elif (np.isnan(bacteria2) and not np.isnan(bacteria1)) and (imageNumber_forObjects[k + 1] == imageNumber_forObjects[k]):
                end.append(k)
            elif (np.isnan(bacteria1) and not np.isnan(bacteria2)) and not (imageNumber_forObjects[k + 1] == imageNumber_forObjects[k]):
                if len(start) == (len(end) + 1):
                    end.append(k)
                start.append(k + 1)
            elif (bacteria1 > bacteria2) and not (imageNumber_forObjects[k + 1] == imageNumber_forObjects[k]):
                end.append(k)
                start.append(k + 1)
            elif (bacteria1 == bacteria2) and not (imageNumber_forObjects[k + 1] == imageNumber_forObjects[k]):
                end.append(k)
                start.append(k + 1)
            elif (imageNumber_forObjects[k + 1] > imageNumber_forObjects[k]):
                end.append(k)
                start.append(k + 1)
            else:
                pass
            
            if ((imageNumber_forObjects[k + 1] - (imageNumber_forObjects[k])) >= 2):
                nonPresent = (imageNumber_forObjects[k + 1] - (imageNumber_forObjects[k]) - 1)
                while nonPresent:
                    end.append(np.nan)
                    start.append(k + 1)
                    nonPresent = nonPresent - 1
            else:
                pass

            k = k + 1
        if not len(end) == len(start):
            end.append(k)
        
        while len(start) < len(well_number_metadata):
            start.append(k)
            end.append(np.nan)
        
        if len(start) > len(well_number_metadata):
            raise Exception('Something went wrong when determining number of cells per well')
        else:
            pass
            
        return start, end

    def match_index(bacteria_number, image):
        """Find the index of foci belonging to the parent cells"""
        index_list = []
        for i in range(start_DNAIndex[image],end_DNAIndex[image]+1):
            if parentBacteriaIndex[CellCategoryIndex][i] == bacteria_number:
                index_list.append(DNAObjectIndex[CellCategoryIndex][i])            
        return index_list

    def convert_float_list_to_integer(lst):
        """Make sure list is integer for compatibility with different functions"""
        for i in range(0,len(lst)):
            if not np.isnan(lst[i]):
                lst[i] = round(lst[i])
        return lst

    def findNonzeroAverage(fociContent, start, end):
        """Find average and standard deviation of nonzero values in list"""
        contentInContainingCells = [] 
        allValuesAreNAN = True   
        for i in range(start,end):
            if np.isnan(fociContent[i]):
                pass
            elif fociContent[i] > 0:
                if allValuesAreNAN:
                    allValuesAreNAN = False
                contentInContainingCells.append(fociContent[i])
            else:
                pass
        if allValuesAreNAN:
            avg = np.nan
            sd = np.nan
        elif contentInContainingCells != []:
            avg = np.average(contentInContainingCells)
            sd = np.std(contentInContainingCells)
        else:
            avg = 0
            sd = 0
        return avg, sd


    ################  INDEX IMAGES WITH CORRECT PLATE  ################

    plate_date_index = plateDateNumberIndexList.loc[:, "Plate_date"]
    plate_date_index = plate_date_index["Plate_date"].tolist()

    for dates in range(len(plate_date_index)):
        plate_date_index[dates] = str(plate_date_index[dates])

    plate_number_index = plateDateNumberIndexList.loc[:, "Plate_number"]
    plate_number_index = plate_number_index["Plate_number"].tolist()

    plate_number_max = max(plate_number_index)

    dateMatchedtoPlateNumber = []
    for i in range(0,len(plate_date_metadata)):
        dateMatchedtoPlateNumber.append(plate_number_index[plate_date_index.index(str(plate_date_metadata[i]))])
    
    present_wells = [ [] for _ in range(plate_number_max) ]
    present_wells_count = 0
    
    for plate_numb in range(1,plate_number_max+1):
        for r in range(0,len(well_number_metadata)):
            if dateMatchedtoPlateNumber[r] == plate_numb:
                if well_number_metadata[r] not in present_wells[plate_numb-1]:
                    present_wells[plate_numb-1].append(well_number_metadata[r])
                    present_wells_count = present_wells_count + 1
                else:
                    pass
            else:
                pass

    
    
    ################  INDEX WELLS WITH GENES  ################
    
    well_number_index = {}
    gene_index = {}
    treatment_index = []

    for plate_num in range(1,plate_number_max+1):
        well_number_index[plate_num-1] = wellGeneIndexList.loc[:, ("Plate" + str(plate_num), "Well")]
        well_number_index[plate_num-1] = pd.Series.tolist(well_number_index[plate_num-1])

        gene_index[plate_num-1] = wellGeneIndexList.loc[:,("Plate" + str(plate_num),"Gene")]
        gene_index[plate_num-1] = pd.Series.tolist(gene_index[plate_num-1])

        if plate_num >= 14:
            treatment_index = wellGeneIndexList.loc[:, ("Plate" + str(plate_num), "Treatment")]
            treatment_index = pd.Series.tolist(treatment_index)
    
    geneMatchedtoWellNumber = []
    for i in range(0,len(well_number_metadata)):
        if plate_date_metadata[i] == "201006-1" or plate_date_metadata[i] == "201006-2":
            if well_number_metadata[i] == "F05":
                geneMatchedtoWellNumber.append("wt NT")
            elif well_number_metadata[i] == "E05":
                geneMatchedtoWellNumber.append("wt +cip")
            else:
                geneMatchedtoWellNumber.append(gene_index[dateMatchedtoPlateNumber[i]-1][well_number_index[dateMatchedtoPlateNumber[i]-1].index(well_number_metadata[i])])
        else:
            geneMatchedtoWellNumber.append(gene_index[dateMatchedtoPlateNumber[i]-1][well_number_index[dateMatchedtoPlateNumber[i]-1].index(well_number_metadata[i])])
    

    
    ###############  CREATE OUTPUT MATRIX WITH METADATA  ################

    # count0 = 0
    wellInInput = False
    tempDataSorting = np.zeros((present_wells_count,6,(2+3*len(bacteriaCategory)+23*2*(len(bacteriaCategory)+1))))
    

    plate_date_ordered = []
    plate_number_ordered = []
    well_number_ordered = []
    genes_ordered = []
    site_number_ordered = []
    site_number_combined = []
    treatment_ordered = []

    count2 = 0
    count3 = 0

    for pl in range(0,len(plate_date_index)):
        """Order by plate number and well number in final results output"""
        for k in range(0,len(well_number_index[plate_number_index[pl]-1])):
            """Create intermediate matrix that order sites and allows calculation of well averages and std"""
            for m in range(0,len(well_number_metadata)):
                if (well_number_metadata[m] == well_number_index[plate_number_index[pl]-1][k]) and (plate_date_metadata[m] == plate_date_index[pl]):
                    count2 = count2 + 1
                    if not wellInInput:
                        plate_date_ordered.append(plate_date_metadata[m])
                        plate_number_ordered.append(dateMatchedtoPlateNumber[m])
                        well_number_ordered.append(well_number_metadata[m])
                        genes_ordered.append(geneMatchedtoWellNumber[m])
                        if plate_number_index[pl] >= 14:
                            treatment_ordered.append(treatment_index[k])
                        else:
                            treatment_ordered.append(np.nan)
                        wellInInput = True
                    else:
                        pass
                else:
                    pass
            if wellInInput:
                # count0 = count0 + 1
                count3 = count3 + 1
                wellInInput = False
            else:
                pass





    # for i in range(0,len(well_number_metadata)):
    #     if i == 0:
    #         plate_number_current = dateMatchedtoPlateNumber[i]
    #         geneMatchedtoWellNumber.append(gene_index[dateMatchedtoPlateNumber[i]][well_number_index.index(well_number_metadata[i])])
    #     elif plate_date_metadata[i] != plate_date_metadata[i-1]:
    #         plate_number_current = dateMatchedtoPlateNumber[i]
    #         geneMatchedtoWellNumber.append(gene_index[well_number_index.index(well_number_metadata[i])])
    #     else:
    #         geneMatchedtoWellNumber.append(gene_index[well_number_index.index(well_number_metadata[i])])


    ################  CALCULATIONS  ################

    for CellCategoryIndex in range(0, len(bacteriaCategory)):
        """Determine the relevant index ranges for the objects of interest"""
        start_BacIndex, end_BacIndex = find_start_end(bacteriaObjectIndex[CellCategoryIndex])
        start_DNAIndex, end_DNAIndex = find_start_end(DNAObjectIndex[CellCategoryIndex])

        """Make sure object numbering is numeric and not float"""
        bacteriaObjectIndex[CellCategoryIndex] = convert_float_list_to_integer(bacteriaObjectIndex[CellCategoryIndex])
        #cellArea = convert_float_list_to_integer(cellArea)
        #DNAfociCountperBacteria = convert_float_list_to_integer(DNAfociCountperBacteria)

        """Create variables that will be used for calculations and exported"""
        avgbacteriasizelist = []
        stdbacteriasizelist = []
        avgfocicontentlist = []
        stdfocicontentlist = []
        avgAreaCoveragelist = []
        stdAreaCoveragelist = []
        avgDNAcompactionlist = []
        stdDNAcompactionlist = []
        avgEccentricitylist = []
        stdEccentricitylist = []
        avgIntensityinDNAobjectsforBaclist = []
        stdIntensityinDNAobjectsforBaclist = []
        avgIntensityinBACobjectslist = []
        stdIntensityinBACobjectslist = []
        avgMassDisplacementlist = []
        stdMassDisplacementlist = []
        avgMajorAxisLengthlist = []
        stdMajorAxisLengthlist = []
        avgMaxFeretLengthlist = []
        stdMaxFeretLengthlist = []
        numberBacteriasAnalyzedlist = []

        avgTotDNAdistributionlist1 = []
        stdTotDNAdistributionlist1 = []
        avgTotDNAdistributionlist2 = []
        stdTotDNAdistributionlist2 = []
        avgTotDNAdistributionlist3 = []
        stdTotDNAdistributionlist3 = []
        avgTotDNAdistributionlist4 = []
        stdTotDNAdistributionlist4 = []
        avgTotDNAdistributionlist5 = []
        stdTotDNAdistributionlist5 = []
        avgTotDNAdistributionlist6 = []
        stdTotDNAdistributionlist6 = []
        avgTotDNAdistributionlistExtra = []
        stdTotDNAdistributionlistExtra = []

        avgConcDNAdistributionlist1 = []
        stdConcDNAdistributionlist1 = []
        avgConcDNAdistributionlist2 = []
        stdConcDNAdistributionlist2 = []
        avgConcDNAdistributionlist3 = []
        stdConcDNAdistributionlist3 = []
        avgConcDNAdistributionlist4 = []
        stdConcDNAdistributionlist4 = []
        avgConcDNAdistributionlist5 = []
        stdConcDNAdistributionlist5 = []
        avgConcDNAdistributionlist6 = []
        stdConcDNAdistributionlist6 = []
        avgConcDNAdistributionlistExtra = []
        stdConcDNAdistributionlistExtra = []


        """Extra lists not used:"""
        #combinedfociarealist = []
        #bacteriasizelist = []
        #focitocellarearatiolist = []
        #overallDNAcompactionlist = []
        #thresholdcellsizelist = []]



        for n in range(0, len(well_number_metadata)):

            sum_cell_size_L = []
            foci_to_cell_area_ratio_L = []
            DNA_avg_compaction_L = []
            DNA_avg_eccentricity_L = []
            avg_intensities_inDNAobjects_perBac_L = []
            #DNAfociCountperBacteria_L = []
            noCellsforCategoryinImage = True

            if (cellCount[CellCategoryIndex][n] == 0) or (np.isnan(end_BacIndex[n])):
                numberBacteriasAnalyzedlist.append(0)
                    
                avgbacteriasizelist.append(np.nan)
                stdbacteriasizelist.append(np.nan)

                avgfocicontentlist.append(np.nan)
                stdfocicontentlist.append(np.nan)

                avgAreaCoveragelist.append(np.nan)
                stdAreaCoveragelist.append(np.nan)

                avgDNAcompactionlist.append(np.nan)
                stdDNAcompactionlist.append(np.nan)

                avgEccentricitylist.append(np.nan)
                stdEccentricitylist.append(np.nan)

                avgIntensityinDNAobjectsforBaclist.append(np.nan)
                stdIntensityinDNAobjectsforBaclist.append(np.nan)

                avgIntensityinBACobjectslist.append(np.nan)
                stdIntensityinBACobjectslist.append(np.nan)

                avgMassDisplacementlist.append(np.nan)
                stdMassDisplacementlist.append(np.nan)

                avgMajorAxisLengthlist.append(np.nan)
                stdMajorAxisLengthlist.append(np.nan)

                avgMaxFeretLengthlist.append(np.nan)
                stdMaxFeretLengthlist.append(np.nan)

                avgTotDNAdistributionlist1.append(np.nan)
                stdTotDNAdistributionlist1.append(np.nan)

                avgTotDNAdistributionlist2.append(np.nan)
                stdTotDNAdistributionlist2.append(np.nan)

                avgTotDNAdistributionlist3.append(np.nan)
                stdTotDNAdistributionlist3.append(np.nan)

                avgTotDNAdistributionlist4.append(np.nan)
                stdTotDNAdistributionlist4.append(np.nan)

                avgTotDNAdistributionlist5.append(np.nan)
                stdTotDNAdistributionlist5.append(np.nan)

                avgTotDNAdistributionlist6.append(np.nan)
                stdTotDNAdistributionlist6.append(np.nan)

                avgTotDNAdistributionlistExtra.append(np.nan)
                stdTotDNAdistributionlistExtra.append(np.nan)

                avgConcDNAdistributionlist1.append(np.nan)
                stdConcDNAdistributionlist1.append(np.nan)

                avgConcDNAdistributionlist2.append(np.nan)
                stdConcDNAdistributionlist2.append(np.nan)

                avgConcDNAdistributionlist3.append(np.nan)
                stdConcDNAdistributionlist3.append(np.nan)

                avgConcDNAdistributionlist4.append(np.nan)
                stdConcDNAdistributionlist4.append(np.nan)

                avgConcDNAdistributionlist5.append(np.nan)
                stdConcDNAdistributionlist5.append(np.nan)

                avgConcDNAdistributionlist6.append(np.nan)
                stdConcDNAdistributionlist6.append(np.nan)

                avgConcDNAdistributionlistExtra.append(np.nan)
                stdConcDNAdistributionlistExtra.append(np.nan)
            else:
                for j in bacteriaObjectIndex[CellCategoryIndex][start_BacIndex[n]:(end_BacIndex[n]+1)]:
                    if np.isnan(j):
                        continue
                    else:
                        noCellsforCategoryinImage = False

                        """Find the area of the bacteria"""
                        cell_size = cellArea[CellCategoryIndex][start_BacIndex[n] + j - 1]
                        
                        """Skip cell if the following values are below chosen threshold"""
                        if cell_size < 30 or (bacteriaMinorAxisLength[CellCategoryIndex][start_BacIndex[n] + j - 1] < 4 or bacteriaMinFeretLength[CellCategoryIndex][start_BacIndex[n] + j - 1] < 4):
                            """Enter correct threshold value for cell size (60 correspond to a circular cell of 1 um diameter)"""
                            """Also for cell width (6.5 correspond to a width of 0.75 um)"""
                            continue
                        #elif DNAfociCountperBacteria[start_BacIndex[n] + j - 1] == 0:
                        #    continue
                        else:
                            pass

                        sum_foci_area = 0
                        Foci_compaction_in_cell = []
                        Foci_eccentricity_in_cell = []
                        DNAobjectsinCell_Intensities = []   

                        """Obtain index of daughter foci"""
                        index_list = match_index(bacteria_number=j, image = n)
                        index_list = convert_float_list_to_integer(index_list)

                        for index in index_list:
                            """Find the sum of foci area (foci from the same parent cell is combined into a total foci area)"""
                            sum_foci_area = sum_foci_area + DNAfociArea[CellCategoryIndex][start_DNAIndex[n] + index - 1]
                            
                            """Find the compactness of DNA within the given cell (j)"""
                            Foci_compaction_in_cell.append(DNAcompactness[CellCategoryIndex][start_DNAIndex[n] + index - 1])

                            """Find the eccentricity of DNA within the given cell (j)"""
                            Foci_eccentricity_in_cell.append(DNAeccentricity[CellCategoryIndex][start_DNAIndex[n] + index - 1])

                            """Find the intensities of DNA within DNA objects in the given cell (j)"""
                            DNAobjectsinCell_Intensities.append(DNAMeanIntensities[CellCategoryIndex][start_DNAIndex[n] + index - 1])

                        """Calculate cell coverage by focis and store results for area calculations"""
                        foci_to_cell_area_ratio = sum_foci_area / cell_size
                        sum_cell_size_L.append(cell_size)
                        foci_to_cell_area_ratio_L.append(foci_to_cell_area_ratio)

                        """Number of DNA focis in bacteria"""
                        #DNAfociCountperBacteria_L.append(DNAfociCountperBacteria[start_BacIndex[n] + j - 1])
                        
                        """Find the average compactness of DNA within the cell"""
                        if Foci_compaction_in_cell == []:
                            DNA_avg_compaction_L.append(0)
                        else:
                            avgCellCompaction, _ = findNonzeroAverage(Foci_compaction_in_cell, 0, len(Foci_compaction_in_cell))
                            DNA_avg_compaction_L.append(avgCellCompaction)
                        
                        """Find the average eccentricity of DNA within the cell"""
                        if Foci_eccentricity_in_cell == []:
                            DNA_avg_eccentricity_L.append(0)
                        else:
                            avgCellEccentricity, _ = findNonzeroAverage(Foci_eccentricity_in_cell, 0, len(Foci_eccentricity_in_cell))
                            DNA_avg_eccentricity_L.append(avgCellEccentricity)
                        
                        """Find the average intensity inside DNA objects for the given cell"""
                        if DNAobjectsinCell_Intensities == []:
                            avg_intensities_inDNAobjects_perBac_L.append(0)
                        else:
                            avgDNAIntensity, _ = findNonzeroAverage(DNAobjectsinCell_Intensities, 0, len(DNAobjectsinCell_Intensities))
                            avg_intensities_inDNAobjects_perBac_L.append(avgDNAIntensity)
                
                if noCellsforCategoryinImage:
                    numberBacteriasAnalyzedlist.append(0)
                    
                    avgbacteriasizelist.append(np.nan)
                    stdbacteriasizelist.append(np.nan)

                    avgfocicontentlist.append(np.nan)
                    stdfocicontentlist.append(np.nan)

                    avgAreaCoveragelist.append(np.nan)
                    stdAreaCoveragelist.append(np.nan)

                    avgDNAcompactionlist.append(np.nan)
                    stdDNAcompactionlist.append(np.nan)

                    avgEccentricitylist.append(np.nan)
                    stdEccentricitylist.append(np.nan)

                    avgIntensityinDNAobjectsforBaclist.append(np.nan)
                    stdIntensityinDNAobjectsforBaclist.append(np.nan)

                    avgIntensityinBACobjectslist.append(np.nan)
                    stdIntensityinBACobjectslist.append(np.nan)

                    avgMassDisplacementlist.append(np.nan)
                    stdMassDisplacementlist.append(np.nan)

                    avgMajorAxisLengthlist.append(np.nan)
                    stdMajorAxisLengthlist.append(np.nan)

                    avgMaxFeretLengthlist.append(np.nan)
                    stdMaxFeretLengthlist.append(np.nan)

                    avgTotDNAdistributionlist1.append(np.nan)
                    stdTotDNAdistributionlist1.append(np.nan)

                    avgTotDNAdistributionlist2.append(np.nan)
                    stdTotDNAdistributionlist2.append(np.nan)

                    avgTotDNAdistributionlist3.append(np.nan)
                    stdTotDNAdistributionlist3.append(np.nan)

                    avgTotDNAdistributionlist4.append(np.nan)
                    stdTotDNAdistributionlist4.append(np.nan)

                    avgTotDNAdistributionlist5.append(np.nan)
                    stdTotDNAdistributionlist5.append(np.nan)

                    avgTotDNAdistributionlist6.append(np.nan)
                    stdTotDNAdistributionlist6.append(np.nan)

                    avgTotDNAdistributionlistExtra.append(np.nan)
                    stdTotDNAdistributionlistExtra.append(np.nan)

                    avgConcDNAdistributionlist1.append(np.nan)
                    stdConcDNAdistributionlist1.append(np.nan)

                    avgConcDNAdistributionlist2.append(np.nan)
                    stdConcDNAdistributionlist2.append(np.nan)

                    avgConcDNAdistributionlist3.append(np.nan)
                    stdConcDNAdistributionlist3.append(np.nan)

                    avgConcDNAdistributionlist4.append(np.nan)
                    stdConcDNAdistributionlist4.append(np.nan)

                    avgConcDNAdistributionlist5.append(np.nan)
                    stdConcDNAdistributionlist5.append(np.nan)

                    avgConcDNAdistributionlist6.append(np.nan)
                    stdConcDNAdistributionlist6.append(np.nan)

                    avgConcDNAdistributionlistExtra.append(np.nan)
                    stdConcDNAdistributionlistExtra.append(np.nan)
                else:
                    """Store number of bacterias analyzed for the image"""
                    numberBacteriasAnalyzedlist.append(len(sum_cell_size_L))

                    """Determine average size of bacteria (pixels squared), along with standard deviation"""
                    avgBacSize, stdBacSize = findNonzeroAverage(sum_cell_size_L, 0, len(sum_cell_size_L))
                    avgbacteriasizelist.append(avgBacSize)
                    stdbacteriasizelist.append(stdBacSize)
                    
                    """Determine average foci content in cells containing DNA, along with standard deviation"""
                    # avgFoci, stdFoci = findNonzeroAverage(DNAfociCountperBacteria_L, 0, len(DNAfociCountperBacteria_L))
                    avgFoci, stdFoci = findNonzeroAverage(DNAfociCountperBacteria[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgfocicontentlist.append(avgFoci)
                    stdfocicontentlist.append(stdFoci)
                    
                    """Determine average cell coverage by focis in cells containing DNA, along with standard deviation"""
                    avgAreaCoverage, stdAreaCoverage = findNonzeroAverage(foci_to_cell_area_ratio_L, 0, len(foci_to_cell_area_ratio_L))
                    avgAreaCoveragelist.append(avgAreaCoverage)
                    stdAreaCoveragelist.append(stdAreaCoverage)

                    """Determine average compaction of DNA within all cells in image, along with standard deviation"""
                    avgImageCompaction, stdImageCompaction = findNonzeroAverage(DNA_avg_compaction_L, 0, len(DNA_avg_compaction_L))
                    avgDNAcompactionlist.append(avgImageCompaction)
                    stdDNAcompactionlist.append(stdImageCompaction)

                    """Determine average eccentricity of DNA within all cells in image, along with standard deviation"""
                    avgImageEccentricity, stdImageEccentricity = findNonzeroAverage(DNA_avg_eccentricity_L, 0, len(DNA_avg_eccentricity_L))
                    avgEccentricitylist.append(avgImageEccentricity)
                    stdEccentricitylist.append(stdImageEccentricity)

                    """Determine average intensity of DNA inside DNA objects for each bacteria, along with standard deviation"""
                    avgDNAint, stdDNAint = findNonzeroAverage(avg_intensities_inDNAobjects_perBac_L, 0, len(avg_intensities_inDNAobjects_perBac_L))
                    avgIntensityinDNAobjectsforBaclist.append(avgDNAint)
                    stdIntensityinDNAobjectsforBaclist.append(stdDNAint)

                    """Determine average intensity of DNA inside whole bacteria for each cell, along with standard deviation"""
                    avgBacInt, stdBacInt = findNonzeroAverage(cellularMeanIntensities[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgIntensityinBACobjectslist.append(avgBacInt)
                    stdIntensityinBACobjectslist.append(stdBacInt)

                    """Determine average displacement of DNA intensities from cell center"""
                    avgMassDisp, stdMassDisp = findNonzeroAverage(cellularDNAmassDisplacement[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgMassDisplacementlist.append(avgMassDisp)
                    stdMassDisplacementlist.append(stdMassDisp)

                    """Determine average length of bacteria, along with standard deviation"""
                    avgBacLengthMajorAxis, stdBacLengthMajorAxis = findNonzeroAverage(bacteriaMajorAxisLength[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgMajorAxisLengthlist.append(avgBacLengthMajorAxis)
                    stdMajorAxisLengthlist.append(stdBacLengthMajorAxis)

                    """Determine average length of bacteria, along with standard deviation"""
                    avgBacLength, stdBacLength = findNonzeroAverage(bacteriaMaxFeretLength[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgMaxFeretLengthlist.append(avgBacLength)
                    stdMaxFeretLengthlist.append(stdBacLength)
                            
                    """Determine average and std of total intensity distribution of DNA around DNA center in cell"""
                    avgdistTotDNA1, stddistTotDNA1 = findNonzeroAverage(distributionTotDNALevel1[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgTotDNAdistributionlist1.append(avgdistTotDNA1)
                    stdTotDNAdistributionlist1.append(stddistTotDNA1)

                    avgdistTotDNA2, stddistTotDNA2 = findNonzeroAverage(distributionTotDNALevel2[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgTotDNAdistributionlist2.append(avgdistTotDNA2)
                    stdTotDNAdistributionlist2.append(stddistTotDNA2)

                    avgdistTotDNA3, stddistTotDNA3 = findNonzeroAverage(distributionTotDNALevel3[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgTotDNAdistributionlist3.append(avgdistTotDNA3)
                    stdTotDNAdistributionlist3.append(stddistTotDNA3)

                    avgdistTotDNA4, stddistTotDNA4 = findNonzeroAverage(distributionTotDNALevel4[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgTotDNAdistributionlist4.append(avgdistTotDNA4)
                    stdTotDNAdistributionlist4.append(stddistTotDNA4)

                    avgdistTotDNA5, stddistTotDNA5 = findNonzeroAverage(distributionTotDNALevel5[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgTotDNAdistributionlist5.append(avgdistTotDNA5)
                    stdTotDNAdistributionlist5.append(stddistTotDNA5)

                    avgdistTotDNA6, stddistTotDNA6 = findNonzeroAverage(distributionTotDNALevel6[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgTotDNAdistributionlist6.append(avgdistTotDNA6)
                    stdTotDNAdistributionlist6.append(stddistTotDNA6)

                    avgdistTotDNAExtra, stddistTotDNAExtra = findNonzeroAverage(distributionTotDNALevelExtra[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgTotDNAdistributionlistExtra.append(avgdistTotDNAExtra)
                    stdTotDNAdistributionlistExtra.append(stddistTotDNAExtra)


                    """Determine average and std of DNA intensity concentration (total intensity distribution of DNA around DNA center in cell, normalized by number of pixels)"""
                    avgdistConcDNA1, stddistConcDNA1 = findNonzeroAverage(distributionConcDNALevel1[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgConcDNAdistributionlist1.append(avgdistConcDNA1)
                    stdConcDNAdistributionlist1.append(stddistConcDNA1)

                    avgdistConcDNA2, stddistConcDNA2 = findNonzeroAverage(distributionConcDNALevel2[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgConcDNAdistributionlist2.append(avgdistConcDNA2)
                    stdConcDNAdistributionlist2.append(stddistConcDNA2)

                    avgdistConcDNA3, stddistConcDNA3 = findNonzeroAverage(distributionConcDNALevel3[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgConcDNAdistributionlist3.append(avgdistConcDNA3)
                    stdConcDNAdistributionlist3.append(stddistConcDNA3)

                    avgdistConcDNA4, stddistConcDNA4 = findNonzeroAverage(distributionConcDNALevel4[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgConcDNAdistributionlist4.append(avgdistConcDNA4)
                    stdConcDNAdistributionlist4.append(stddistConcDNA4)

                    avgdistConcDNA5, stddistConcDNA5 = findNonzeroAverage(distributionConcDNALevel5[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgConcDNAdistributionlist5.append(avgdistConcDNA5)
                    stdConcDNAdistributionlist5.append(stddistConcDNA5)

                    avgdistConcDNA6, stddistConcDNA6 = findNonzeroAverage(distributionConcDNALevel6[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgConcDNAdistributionlist6.append(avgdistConcDNA6)
                    stdConcDNAdistributionlist6.append(stddistConcDNA6)

                    avgdistConcDNAExtra, stddistConcDNAExtra = findNonzeroAverage(distributionConcDNALevelExtra[CellCategoryIndex], start_BacIndex[n], end_BacIndex[n]+1)
                    avgConcDNAdistributionlistExtra.append(avgdistConcDNAExtra)
                    stdConcDNAdistributionlistExtra.append(stddistConcDNAExtra)
                    

                
            # duration = time.process_time()
            # print("Finished Image Set " + str(n + 1) + " of " + str(len(start_BacIndex)) + ". Time: " + str(duration))


        ################  SORT RESULTS FOR EASIER OUTPUT  ################

        count0 = 0
        count1 = 0
        wellInInput = False

        for pl in range(0,len(plate_date_index)):
            """Order by plate number and well number in final results output"""
            for k in range(0,len(well_number_index[plate_number_index[pl]-1])):
                """Create intermediate matrix that order sites and allows calculation of well averages and std"""
                for m in range(0,len(well_number_metadata)):
                    if (well_number_metadata[m] == well_number_index[plate_number_index[pl]-1][k]) and (plate_date_metadata[m] == plate_date_index[pl]):
                        count1=count1+1
                        if plate_date_ordered[count0] == plate_date_metadata[m] and plate_number_ordered[count0] == dateMatchedtoPlateNumber[m] and well_number_ordered[count0] == well_number_metadata[m] and genes_ordered[count0] == geneMatchedtoWellNumber[m]:
                            wellInInput = True
                        elif pd.isnull(genes_ordered[count0]):
                            continue
                        else:
                            raise Exception('Something went wrong with sorting of results')
                        
                        

                        """For all calculated values, will have value for all categories"""
                        tempDataSorting[count0,site_number_metadata[m],0] = site_number_metadata[m]
                        tempDataSorting[count0,site_number_metadata[m],(2+3*CellCategoryIndex)] = numberBacteriasAnalyzedlist[m]
                        tempDataSorting[count0,site_number_metadata[m],(3+3*CellCategoryIndex)] = 0     #Add value for p(Enrichment score)
                        tempDataSorting[count0,site_number_metadata[m],(4+3*CellCategoryIndex)] = 0     #Add value for Enrichment score
                        tempDataSorting[count0,site_number_metadata[m],(4+3*len(bacteriaCategory)+2*CellCategoryIndex)] = avgfocicontentlist[m]
                        tempDataSorting[count0,site_number_metadata[m],(5+3*len(bacteriaCategory)+2*CellCategoryIndex)] = stdfocicontentlist[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+1*2*(len(bacteriaCategory)+1))] = avgAreaCoveragelist[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+1*2*(len(bacteriaCategory)+1))] = stdAreaCoveragelist[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+2*2*(len(bacteriaCategory)+1))] = avgDNAcompactionlist[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+2*2*(len(bacteriaCategory)+1))] = stdDNAcompactionlist[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+3*2*(len(bacteriaCategory)+1))] = avgEccentricitylist[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+3*2*(len(bacteriaCategory)+1))] = stdEccentricitylist[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+4*2*(len(bacteriaCategory)+1))] = avgIntensityinBACobjectslist[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+4*2*(len(bacteriaCategory)+1))] = stdIntensityinBACobjectslist[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+5*2*(len(bacteriaCategory)+1))] = avgIntensityinDNAobjectsforBaclist[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+5*2*(len(bacteriaCategory)+1))] = stdIntensityinDNAobjectsforBaclist[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+6*2*(len(bacteriaCategory)+1))] = avgMassDisplacementlist[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+6*2*(len(bacteriaCategory)+1))] = stdMassDisplacementlist[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+7*2*(len(bacteriaCategory)+1))] = avgMajorAxisLengthlist[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+7*2*(len(bacteriaCategory)+1))] = stdMajorAxisLengthlist[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+8*2*(len(bacteriaCategory)+1))] = avgMaxFeretLengthlist[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+8*2*(len(bacteriaCategory)+1))] = stdMaxFeretLengthlist[m]
                        
                        
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+9*2*(len(bacteriaCategory)+1))] = avgTotDNAdistributionlist1[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+9*2*(len(bacteriaCategory)+1))] = stdTotDNAdistributionlist1[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+10*2*(len(bacteriaCategory)+1))] = avgTotDNAdistributionlist2[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+10*2*(len(bacteriaCategory)+1))] = stdTotDNAdistributionlist2[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+11*2*(len(bacteriaCategory)+1))] = avgTotDNAdistributionlist3[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+11*2*(len(bacteriaCategory)+1))] = stdTotDNAdistributionlist3[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+12*2*(len(bacteriaCategory)+1))] = avgTotDNAdistributionlist4[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+12*2*(len(bacteriaCategory)+1))] = stdTotDNAdistributionlist4[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+13*2*(len(bacteriaCategory)+1))] = avgTotDNAdistributionlist5[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+13*2*(len(bacteriaCategory)+1))] = stdTotDNAdistributionlist5[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+14*2*(len(bacteriaCategory)+1))] = avgTotDNAdistributionlist6[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+14*2*(len(bacteriaCategory)+1))] = stdTotDNAdistributionlist6[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+15*2*(len(bacteriaCategory)+1))] = avgTotDNAdistributionlistExtra[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+15*2*(len(bacteriaCategory)+1))] = stdTotDNAdistributionlistExtra[m]
                        
                        
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+16*2*(len(bacteriaCategory)+1))] = avgConcDNAdistributionlist1[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+16*2*(len(bacteriaCategory)+1))] = stdConcDNAdistributionlist1[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+17*2*(len(bacteriaCategory)+1))] = avgConcDNAdistributionlist2[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+17*2*(len(bacteriaCategory)+1))] = stdConcDNAdistributionlist2[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+18*2*(len(bacteriaCategory)+1))] = avgConcDNAdistributionlist3[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+18*2*(len(bacteriaCategory)+1))] = stdConcDNAdistributionlist3[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+19*2*(len(bacteriaCategory)+1))] = avgConcDNAdistributionlist4[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+19*2*(len(bacteriaCategory)+1))] = stdConcDNAdistributionlist4[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+20*2*(len(bacteriaCategory)+1))] = avgConcDNAdistributionlist5[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+20*2*(len(bacteriaCategory)+1))] = stdConcDNAdistributionlist5[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+21*2*(len(bacteriaCategory)+1))] = avgConcDNAdistributionlist6[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+21*2*(len(bacteriaCategory)+1))] = stdConcDNAdistributionlist6[m]
                        tempDataSorting[count0,site_number_metadata[m],((4+3*len(bacteriaCategory)+2*CellCategoryIndex)+22*2*(len(bacteriaCategory)+1))] = avgConcDNAdistributionlistExtra[m]
                        tempDataSorting[count0,site_number_metadata[m],((5+3*len(bacteriaCategory)+2*CellCategoryIndex)+22*2*(len(bacteriaCategory)+1))] = stdConcDNAdistributionlistExtra[m]
                    else:
                        pass
                
                if wellInInput:
                    count0 = count0 + 1
                    wellInInput = False
                else:
                    pass

        duration = time.process_time()
        print("Finished calculations for all image sets in category: " + str(bacteriaCategory[CellCategoryIndex]) + ". Time: " + str(duration))



    """Determine values for all cells in well"""
    for count in range(0,tempDataSorting.shape[0]):
        for sites in range(1,5):
            tempDataSorting[count,sites,1] = sum(tempDataSorting[count,sites,(2+3*CategoryIndex)] for CategoryIndex in range(len(bacteriaCategory)))

            std_combined = 0
            mean_combined = 0
            number_combined = 0

            for measurement in range(0,23):
                std_combined = 0
                mean_combined = 0
                number_combined = 0
                for CategoryIndex in range(len(bacteriaCategory)):
                    if np.isnan(tempDataSorting[count,sites,(2+3*CategoryIndex)]):
                        continue
                    elif np.isnan(tempDataSorting[count,sites,((4+3*len(bacteriaCategory)+2*CategoryIndex)+measurement*2*(len(bacteriaCategory)+1))]):
                        continue
                    elif np.isnan(tempDataSorting[count,sites,((5+3*len(bacteriaCategory)+2*CategoryIndex)+measurement*2*(len(bacteriaCategory)+1))]):
                        continue
                    elif tempDataSorting[count,sites,(2+3*CategoryIndex)] <= 1:
                        continue
                    else:
                        std_combined = np.sqrt( ((number_combined-1)*(std_combined**2) + (tempDataSorting[count,sites,(2+3*CategoryIndex)]-1)*(tempDataSorting[count,sites,((5+3*len(bacteriaCategory)+2*CategoryIndex)+measurement*2*(len(bacteriaCategory)+1))]**2) + (number_combined*tempDataSorting[count,sites,(2+3*CategoryIndex)] / (number_combined + tempDataSorting[count,sites,(2+3*CategoryIndex)])) * (mean_combined**2 + tempDataSorting[count,sites,((4+3*len(bacteriaCategory)+2*CategoryIndex)+measurement*2*(len(bacteriaCategory)+1))]**2 - 2*mean_combined*tempDataSorting[count,sites,((4+3*len(bacteriaCategory)+2*CategoryIndex)+measurement*2*(len(bacteriaCategory)+1))])) / (number_combined + tempDataSorting[count,sites,(2+3*CategoryIndex)] - 1) )
                        mean_combined = ((mean_combined*number_combined + tempDataSorting[count,sites,((4+3*len(bacteriaCategory)+2*CategoryIndex)+measurement*2*(len(bacteriaCategory)+1))]*tempDataSorting[count,sites,(2+3*CategoryIndex)]) / (number_combined + tempDataSorting[count,sites,(2+3*CategoryIndex)]))
                        number_combined = number_combined + tempDataSorting[count,sites,(2+3*CategoryIndex)]
                
                if std_combined == 0 and (mean_combined == 0 and number_combined == 0):
                    tempDataSorting[count,sites,2+3*len(bacteriaCategory)+measurement*2*(len(bacteriaCategory)+1)] = np.nan
                    tempDataSorting[count,sites,3+3*len(bacteriaCategory)+measurement*2*(len(bacteriaCategory)+1)] = np.nan
                else:
                    tempDataSorting[count,sites,2+3*len(bacteriaCategory)+measurement*2*(len(bacteriaCategory)+1)] = mean_combined
                    tempDataSorting[count,sites,3+3*len(bacteriaCategory)+measurement*2*(len(bacteriaCategory)+1)] = std_combined

    """Find average of all sites present for each well"""
    for l in range(0,tempDataSorting.shape[0]): 
        """number of rows in matrix (x)"""
        site_number_ordered.append("Site_distributions")
        tempDataSorting[l,0,1], _ = findNonzeroAverage(tempDataSorting[l,:,1],1,5)

        for CellCategoryIndex in range(0,3*len(bacteriaCategory)):
            tempDataSorting[l,0,(2+CellCategoryIndex)], _ = findNonzeroAverage(tempDataSorting[l,:,(2+CellCategoryIndex)],1,5)

        for h in range(2+3*len(bacteriaCategory),tempDataSorting.shape[2],2): 
            """depth of matrix (z)"""
            tempDataSorting[l,0,h], tempDataSorting[l,0,h+1] = findNonzeroAverage(tempDataSorting[l,:,h],1,5)


    """Find combined distribution for each well for all sites combined"""
    for nn in range(0,tempDataSorting.shape[0]):
        """Number of rows in matrix (x)"""
        site_number_combined.append("Combined_distributions")
        for CellSum in range(1,(len(bacteriaCategory)+2)):
            if CellSum <= 2:
                tempDataSorting[nn,5,CellSum] = sum(tempDataSorting[nn,site,CellSum] for site in range(1,5))
            else:
                tempDataSorting[nn,5,(2+3*(CellSum-2))] = sum(tempDataSorting[nn,site,(2+3*(CellSum-2))] for site in range(1,5))

        for mm in range(2+3*len(bacteriaCategory),tempDataSorting.shape[2],8): 
            """depth of matrix (z)"""
            
            for CalculationIndex in range(0,len(bacteriaCategory)+1):
                """Calculate mean and standard deviation of entire well as if the well had never been divided into 4 sites"""
                """See 'Cochrane Handbook for Systematic Reviews of Interventions, 7.7.3.8' for calculation basis"""
                if CalculationIndex == 0:
                    CountIndex = 1
                else:
                    CountIndex = 3*CalculationIndex - 1
                std_combined = 0
                mean_combined = 0
                number_combined = 0
                abs(CalculationIndex-1)
                for sites in range(1,5):
                    if np.isnan(tempDataSorting[nn,sites,CountIndex]):
                        continue
                    elif np.isnan(tempDataSorting[nn,sites,(mm+1+2*CalculationIndex)]):
                        continue
                    elif np.isnan(tempDataSorting[nn,sites,(mm+2*CalculationIndex)]):
                        continue
                    elif tempDataSorting[nn,sites,CountIndex] <= 1:
                        continue
                    else:
                        std_combined = np.sqrt( ((number_combined-1)*(std_combined**2) + (tempDataSorting[nn,sites,CountIndex]-1)*(tempDataSorting[nn,sites,(mm+1+2*CalculationIndex)]**2) + (number_combined*tempDataSorting[nn,sites,CountIndex] / (number_combined + tempDataSorting[nn,sites,CountIndex])) * (mean_combined**2 + tempDataSorting[nn,sites,(mm+2*CalculationIndex)]**2 - 2*mean_combined*tempDataSorting[nn,sites,(mm+2*CalculationIndex)])) / (number_combined + tempDataSorting[nn,sites,CountIndex] - 1) )
                        mean_combined = ((mean_combined*number_combined + tempDataSorting[nn,sites,(mm+2*CalculationIndex)]*tempDataSorting[nn,sites,CountIndex]) / (number_combined + tempDataSorting[nn,sites,CountIndex]))
                        number_combined = number_combined + tempDataSorting[nn,sites,CountIndex]
            
                """The resulting mean and std will now represent all bacteria analyzed in the well"""
                if std_combined == 0 and (mean_combined == 0 and number_combined == 0):
                    tempDataSorting[nn,5,(mm+2*CalculationIndex)] = np.nan
                    tempDataSorting[nn,5,(mm+1+2*CalculationIndex)] = np.nan
                else:
                    tempDataSorting[nn,5,(mm+2*CalculationIndex)] = mean_combined
                    tempDataSorting[nn,5,(mm+1+2*CalculationIndex)] = std_combined


    """Determine Enrichment scores for all wells in all outputs"""

    counts_AllSites = [ [] for _ in range(len(bacteriaCategory)) ]
    for eachSite in range(1,5):
        for eachCategory in range(0,len(bacteriaCategory)):
            counts_AllSites[eachCategory] += tempDataSorting[:,eachSite,(2+3*eachCategory)].tolist()
    counts_AllSites = np.transpose(counts_AllSites)
    alpha_AllSites, converged_AllSites = polyafit.fit_betabinom_minka_alternating(counts_AllSites)

    counts_Combined = [ [] for _ in range(len(bacteriaCategory)) ]
    for eachCategory in range(0,len(bacteriaCategory)):
        counts_Combined[eachCategory] += tempDataSorting[:,5,(2+3*eachCategory)].tolist()
    counts_Combined = np.transpose(counts_Combined)
    alpha_Combined, converged_Combined = polyafit.fit_betabinom_minka_alternating(counts_Combined)
    
    if not converged_AllSites:
        raise Exception("The Pólya fit estimation did not converge for the individual sites")
    elif not converged_Combined:
        raise Exception("The Pólya fit estimation did not converge for the combined cell counts in wells")
    else:
        pass

    # if np.isnan(alpha_AllSites).any():
    #     raise Exception("The alpha value from the Pólya fit is 'Nan' for the individual sites"
    # elif np.isnan(alpha_Combined).any():
    #     raise Exception("The alpha value from the Pólya fit is 'Nan' for the combined cell counts in wells")
    # else:
    #     pass

    for eachSite in range(0,tempDataSorting.shape[1]):
        for eachWell in range(0,tempDataSorting.shape[0]):
            counts = []
            for eachCategory in range(0,len(bacteriaCategory)):
                if np.isnan(tempDataSorting[eachWell,eachSite,(2+3*eachCategory)]):
                    counts.append(0)
                else:
                    counts.append(tempDataSorting[eachWell,eachSite,(2+3*eachCategory)])
            
            scores = []
            """Determine Enrichment scores (probability) for well and/or site"""
            if eachSite == 5:
                if np.isnan(alpha_Combined).any():
                    scores = [ np.nan for _ in range(len(bacteriaCategory)) ]
                else: 
                    scores = np.array(dirichletintegrate.score(alpha_Combined, np.array(counts)))
            else:
                if np.isnan(alpha_AllSites).any():
                    scores = [ np.nan for _ in range(len(bacteriaCategory)) ]
                else:
                    scores = np.array(dirichletintegrate.score(alpha_AllSites, np.array(counts)))

            """clamp scores to [0,1]"""
            if not np.isnan(scores).any():
                scores[scores > 1.] = 1.
                scores[scores < 0.] = 0.

            """Add Enrichment score (probability) to data table"""
            for eachCategory in range(0,len(bacteriaCategory)):
                tempDataSorting[eachWell,eachSite,(3+3*eachCategory)] = scores[eachCategory]
                
            """Calculate logit Enrichment scores and add to data table"""
            """Special case: only calculate logit of "positives" for 2-classes"""
            if two_classes:
                for eachCategory in range(0,len(bacteriaCategory)):
                    tempDataSorting[eachWell,eachSite,(4+3*eachCategory)] = [np.log10(scores[0]) - (np.log10(1 - scores[0]))]
            else:
                for eachCategory in range(0,len(bacteriaCategory)):
                    x = [np.log10(scores[eachCategory]) - (np.log10(1 - scores[eachCategory]))]
                    tempDataSorting[eachWell,eachSite,(4+3*eachCategory)] += x
                    # tempDataSorting[eachWell,eachSite,(4+3*eachCategory)] = [np.log10(scores[eachCategory]) - (np.log10(1 - scores[eachCategory]))]
                
            


    
    ################  CREATE DATAFRAMES FROM SORTED RESULTS  ################

    SiteExists = [False,False,False,False,False,False]

    ExcelSheetsforPrint = {}

    """Create dataframe from averages of sites per well"""
    ExcelSheetsforPrint[0] = pd.DataFrame({'Manually_discard_data_from_analyses': False,
        'Plate_date': plate_date_ordered, 
        'Plate_number': plate_number_ordered,
        'Well_no': well_number_ordered,
        'Gene': genes_ordered,
        'Treatment': treatment_ordered,
        'Site_no': site_number_ordered,
        'Avg_All_Bacteria_count': tempDataSorting[:,0,1],
        'Avg_'+bacteriaCategory[0]+'_count': tempDataSorting[:,0,2],
        'Avg_'+bacteriaCategory[0]+'_p(Enrichment_score)': tempDataSorting[:,0,3],
        'Avg_'+bacteriaCategory[0]+'_logit_Enrichment_score': tempDataSorting[:,0,4],
        'Avg_'+bacteriaCategory[1]+'_count': tempDataSorting[:,0,5],
        'Avg_'+bacteriaCategory[1]+'_p(Enrichment_score)': tempDataSorting[:,0,6],
        'Avg_'+bacteriaCategory[1]+'_logit_Enrichment_score': tempDataSorting[:,0,7],
        'Avg_'+bacteriaCategory[2]+'_count': tempDataSorting[:,0,8],
        'Avg_'+bacteriaCategory[2]+'_p(Enrichment_score)': tempDataSorting[:,0,9],
        'Avg_'+bacteriaCategory[2]+'_logit_Enrichment_score': tempDataSorting[:,0,10],
        
        'Avg_All_focis_per_cell': tempDataSorting[:,0,11],
        'Std_All_focis_per_cell': tempDataSorting[:,0,12],
        'Avg_'+bacteriaCategory[0]+'_focis_per_cell': tempDataSorting[:,0,13],
        'Std_'+bacteriaCategory[0]+'_focis_per_cell': tempDataSorting[:,0,14],
        'Avg_'+bacteriaCategory[1]+'_focis_per_cell': tempDataSorting[:,0,15],
        'Std_'+bacteriaCategory[1]+'_focis_per_cell': tempDataSorting[:,0,16],
        'Avg_'+bacteriaCategory[2]+'_focis_per_cell': tempDataSorting[:,0,17],
        'Std_'+bacteriaCategory[2]+'_focis_per_cell': tempDataSorting[:,0,18],

        'Avg_All_DNA_area_coverage': tempDataSorting[:,0,19],
        'Std_All_DNA_area_coverage': tempDataSorting[:,0,20],
        'Avg_'+bacteriaCategory[0]+'_DNA_area_coverage': tempDataSorting[:,0,21],
        'Std_'+bacteriaCategory[0]+'_DNA_area_coverage': tempDataSorting[:,0,22],
        'Avg_'+bacteriaCategory[1]+'_DNA_area_coverage': tempDataSorting[:,0,23],
        'Std_'+bacteriaCategory[1]+'_DNA_area_coverage': tempDataSorting[:,0,24],
        'Avg_'+bacteriaCategory[2]+'_DNA_area_coverage': tempDataSorting[:,0,25],
        'Std_'+bacteriaCategory[2]+'_DNA_area_coverage': tempDataSorting[:,0,26],
        
        'Avg_All_DNA_compaction': tempDataSorting[:,0,27],
        'Std_All_DNA_compaction': tempDataSorting[:,0,28],
        'Avg_'+bacteriaCategory[0]+'_DNA_compaction': tempDataSorting[:,0,29],
        'Std_'+bacteriaCategory[0]+'_DNA_compaction': tempDataSorting[:,0,30],
        'Avg_'+bacteriaCategory[1]+'_DNA_compaction': tempDataSorting[:,0,31],
        'Std_'+bacteriaCategory[1]+'_DNA_compaction': tempDataSorting[:,0,32],
        'Avg_'+bacteriaCategory[2]+'_DNA_compaction': tempDataSorting[:,0,33],
        'Std_'+bacteriaCategory[2]+'_DNA_compaction': tempDataSorting[:,0,34],
        
        'Avg_All_DNA_eccentricity': tempDataSorting[:,0,35],
        'Std_All_DNA_eccentricity': tempDataSorting[:,0,36],
        'Avg_'+bacteriaCategory[0]+'_DNA_eccentricity': tempDataSorting[:,0,37],
        'Std_'+bacteriaCategory[0]+'_DNA_eccentricity': tempDataSorting[:,0,38],
        'Avg_'+bacteriaCategory[1]+'_DNA_eccentricity': tempDataSorting[:,0,39],
        'Std_'+bacteriaCategory[1]+'_DNA_eccentricity': tempDataSorting[:,0,40],
        'Avg_'+bacteriaCategory[2]+'_DNA_eccentricity': tempDataSorting[:,0,41],
        'Std_'+bacteriaCategory[2]+'_DNA_eccentricity': tempDataSorting[:,0,42],
        
        'Avg_All_whole_cell_intensity': tempDataSorting[:,0,43],
        'Std_All_whole_cell_intensity': tempDataSorting[:,0,44],
        'Avg_'+bacteriaCategory[0]+'_whole_cell_intensity': tempDataSorting[:,0,45],
        'Std_'+bacteriaCategory[0]+'_whole_cell_intensity': tempDataSorting[:,0,46],
        'Avg_'+bacteriaCategory[1]+'_whole_cell_intensity': tempDataSorting[:,0,47],
        'Std_'+bacteriaCategory[1]+'_whole_cell_intensity': tempDataSorting[:,0,48],
        'Avg_'+bacteriaCategory[2]+'_whole_cell_intensity': tempDataSorting[:,0,49],
        'Std_'+bacteriaCategory[2]+'_whole_cell_intensity': tempDataSorting[:,0,50],
        
        'Avg_All_DNA_objects_in_cells_intensity': tempDataSorting[:,0,51],
        'Std_All_DNA_objects_in_cells_intensity': tempDataSorting[:,0,52],
        'Avg_'+bacteriaCategory[0]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,0,53],
        'Std_'+bacteriaCategory[0]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,0,54],
        'Avg_'+bacteriaCategory[1]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,0,55],
        'Std_'+bacteriaCategory[1]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,0,56],
        'Avg_'+bacteriaCategory[2]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,0,57],
        'Std_'+bacteriaCategory[2]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,0,58],
        
        'Avg_All_Intensity_mass_displacement_in_cells': tempDataSorting[:,0,59],
        'Std_All_Intensity_mass_displacement_in_cells': tempDataSorting[:,0,60],
        'Avg_'+bacteriaCategory[0]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,0,61],
        'Std_'+bacteriaCategory[0]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,0,62],
        'Avg_'+bacteriaCategory[1]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,0,63],
        'Std_'+bacteriaCategory[1]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,0,64],
        'Avg_'+bacteriaCategory[2]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,0,65],
        'Std_'+bacteriaCategory[2]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,0,66],
        
        'Avg_All_MajorAxis_Length_of_bacteria': tempDataSorting[:,0,67],
        'Std_All_MajorAxis_Length_of_bacteria': tempDataSorting[:,0,68],
        'Avg_'+bacteriaCategory[0]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,0,69],
        'Std_'+bacteriaCategory[0]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,0,70],
        'Avg_'+bacteriaCategory[1]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,0,71],
        'Std_'+bacteriaCategory[1]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,0,72],
        'Avg_'+bacteriaCategory[2]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,0,73],
        'Std_'+bacteriaCategory[2]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,0,74],
        
        'Avg_All_MaxFeret_Length_of_bacteria': tempDataSorting[:,0,75],
        'Std_All_MaxFeret_Length_of_bacteria': tempDataSorting[:,0,76],
        'Avg_'+bacteriaCategory[0]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,0,77],
        'Std_'+bacteriaCategory[0]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,0,78],
        'Avg_'+bacteriaCategory[1]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,0,79],
        'Std_'+bacteriaCategory[1]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,0,80],
        'Avg_'+bacteriaCategory[2]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,0,81],
        'Std_'+bacteriaCategory[2]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,0,82],
        
        'Avg_All_Tot_DNA_Distribution_1of6': tempDataSorting[:,0,83],
        'Std_All_Tot_DNA_Distribution_1of6': tempDataSorting[:,0,84],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,0,85],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,0,86],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,0,87],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,0,88],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,0,89],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,0,90],
        
        'Avg_All_Tot_DNA_Distribution_2of6': tempDataSorting[:,0,91],
        'Std_All_Tot_DNA_Distribution_2of6': tempDataSorting[:,0,92],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,0,93],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,0,94],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,0,95],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,0,96],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,0,97],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,0,98],
        
        'Avg_All_Tot_DNA_Distribution_3of6': tempDataSorting[:,0,99],
        'Std_All_Tot_DNA_Distribution_3of6': tempDataSorting[:,0,100],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,0,101],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,0,102],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,0,103],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,0,104],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,0,105],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,0,106],
        
        'Avg_All_Tot_DNA_Distribution_4of6': tempDataSorting[:,0,107],
        'Std_All_Tot_DNA_Distribution_4of6': tempDataSorting[:,0,108],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,0,109],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,0,110],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,0,111],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,0,112],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,0,113],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,0,114],
        
        'Avg_All_Tot_DNA_Distribution_5of6': tempDataSorting[:,0,115],
        'Std_All_Tot_DNA_Distribution_5of6': tempDataSorting[:,0,116],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,0,117],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,0,118],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,0,119],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,0,120],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,0,121],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,0,122],
        
        'Avg_All_Tot_DNA_Distribution_6of6': tempDataSorting[:,0,123],
        'Std_All_Tot_DNA_Distribution_6of6': tempDataSorting[:,0,124],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,0,125],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,0,126],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,0,127],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,0,128],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,0,129],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,0,130],
        
        'Avg_All_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,131],
        'Std_All_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,132], 
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,133],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,134], 
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,135],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,136], 
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,137],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,138], 
        
        'Avg_All_Conc_DNA_Distribution_1of6': tempDataSorting[:,0,139],
        'Std_All_Conc_DNA_Distribution_1of6': tempDataSorting[:,0,140],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,0,141],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,0,142],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,0,143],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,0,144],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,0,145],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,0,146],
        
        'Avg_All_Conc_DNA_Distribution_2of6': tempDataSorting[:,0,147],
        'Std_All_Conc_DNA_Distribution_2of6': tempDataSorting[:,0,148],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,0,149],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,0,150],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,0,151],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,0,152],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,0,153],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,0,154],
        
        'Avg_All_Conc_DNA_Distribution_3of6': tempDataSorting[:,0,155],
        'Std_All_Conc_DNA_Distribution_3of6': tempDataSorting[:,0,156],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,0,157],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,0,158],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,0,159],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,0,160],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,0,161],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,0,162],
        
        'Avg_All_Conc_DNA_Distribution_4of6': tempDataSorting[:,0,163],
        'Std_All_Conc_DNA_Distribution_4of6': tempDataSorting[:,0,164],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,0,165],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,0,166],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,0,167],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,0,168],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,0,169],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,0,170],
        
        'Avg_All_Conc_DNA_Distribution_5of6': tempDataSorting[:,0,171],
        'Std_All_Conc_DNA_Distribution_5of6': tempDataSorting[:,0,172],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,0,173],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,0,174],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,0,175],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,0,176],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,0,177],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,0,178],
        
        'Avg_All_Conc_DNA_Distribution_6of6': tempDataSorting[:,0,179],
        'Std_All_Conc_DNA_Distribution_6of6': tempDataSorting[:,0,180],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,0,181],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,0,182],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,0,183],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,0,184],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,0,185],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,0,186],
        
        'Avg_All_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,187],
        'Std_All_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,188],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,189],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,190],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,191],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,192],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,193],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,0,194]})
    SiteExists[0] = True


    """Create dataframes with all information for each site"""
    for z in range(1,5):
        if z in site_number_metadata:
            ExcelSheetsforPrint[z] = pd.DataFrame({'Manually_discard_data_from_analyses': False,
                'Plate_date': plate_date_ordered,
                'Plate_number': plate_number_ordered,
                'Well_no': well_number_ordered,
                'Gene': genes_ordered,
                'Treatment': treatment_ordered,
                'Site_no': tempDataSorting[:,z,0],
                'Bacteria_count': tempDataSorting[:,z,1],
                bacteriaCategory[0]+'_cell_count': tempDataSorting[:,z,2],
                bacteriaCategory[0]+'_p(Enrichment_score)': tempDataSorting[:,z,3],
                bacteriaCategory[0]+'_logit_Enrichment_score': tempDataSorting[:,z,4],
                bacteriaCategory[1]+'_cell_count': tempDataSorting[:,z,5],
                bacteriaCategory[1]+'_p(Enrichment_score)': tempDataSorting[:,z,6],
                bacteriaCategory[1]+'_logit_Enrichment_score': tempDataSorting[:,z,7],
                bacteriaCategory[2]+'_cell_count': tempDataSorting[:,z,8],
                bacteriaCategory[2]+'_p(Enrichment_score)': tempDataSorting[:,z,9],
                bacteriaCategory[2]+'_logit_Enrichment_score': tempDataSorting[:,z,10],
                
                'Avg_All_focis_per_cell': tempDataSorting[:,z,11],
                'Std_All_focis_per_cell': tempDataSorting[:,z,12],
                'Avg_'+bacteriaCategory[0]+'_focis_per_cell': tempDataSorting[:,z,13],
                'Std_'+bacteriaCategory[0]+'_focis_per_cell': tempDataSorting[:,z,14],
                'Avg_'+bacteriaCategory[1]+'_focis_per_cell': tempDataSorting[:,z,15],
                'Std_'+bacteriaCategory[1]+'_focis_per_cell': tempDataSorting[:,z,16],
                'Avg_'+bacteriaCategory[2]+'_focis_per_cell': tempDataSorting[:,z,17],
                'Std_'+bacteriaCategory[2]+'_focis_per_cell': tempDataSorting[:,z,18],

                'Avg_All_DNA_area_coverage': tempDataSorting[:,z,19],
                'Std_All_DNA_area_coverage': tempDataSorting[:,z,20],
                'Avg_'+bacteriaCategory[0]+'_DNA_area_coverage': tempDataSorting[:,z,21],
                'Std_'+bacteriaCategory[0]+'_DNA_area_coverage': tempDataSorting[:,z,22],
                'Avg_'+bacteriaCategory[1]+'_DNA_area_coverage': tempDataSorting[:,z,23],
                'Std_'+bacteriaCategory[1]+'_DNA_area_coverage': tempDataSorting[:,z,24],
                'Avg_'+bacteriaCategory[2]+'_DNA_area_coverage': tempDataSorting[:,z,25],
                'Std_'+bacteriaCategory[2]+'_DNA_area_coverage': tempDataSorting[:,z,26],
                
                'Avg_All_DNA_compaction': tempDataSorting[:,z,27],
                'Std_All_DNA_compaction': tempDataSorting[:,z,28],
                'Avg_'+bacteriaCategory[0]+'_DNA_compaction': tempDataSorting[:,z,29],
                'Std_'+bacteriaCategory[0]+'_DNA_compaction': tempDataSorting[:,z,30],
                'Avg_'+bacteriaCategory[1]+'_DNA_compaction': tempDataSorting[:,z,31],
                'Std_'+bacteriaCategory[1]+'_DNA_compaction': tempDataSorting[:,z,32],
                'Avg_'+bacteriaCategory[2]+'_DNA_compaction': tempDataSorting[:,z,33],
                'Std_'+bacteriaCategory[2]+'_DNA_compaction': tempDataSorting[:,z,34],
                
                'Avg_All_DNA_eccentricity': tempDataSorting[:,z,35],
                'Std_All_DNA_eccentricity': tempDataSorting[:,z,36],
                'Avg_'+bacteriaCategory[0]+'_DNA_eccentricity': tempDataSorting[:,z,37],
                'Std_'+bacteriaCategory[0]+'_DNA_eccentricity': tempDataSorting[:,z,38],
                'Avg_'+bacteriaCategory[1]+'_DNA_eccentricity': tempDataSorting[:,z,39],
                'Std_'+bacteriaCategory[1]+'_DNA_eccentricity': tempDataSorting[:,z,40],
                'Avg_'+bacteriaCategory[2]+'_DNA_eccentricity': tempDataSorting[:,z,41],
                'Std_'+bacteriaCategory[2]+'_DNA_eccentricity': tempDataSorting[:,z,42],
                
                'Avg_All_whole_cell_intensity': tempDataSorting[:,z,43],
                'Std_All_whole_cell_intensity': tempDataSorting[:,z,44],
                'Avg_'+bacteriaCategory[0]+'_whole_cell_intensity': tempDataSorting[:,z,45],
                'Std_'+bacteriaCategory[0]+'_whole_cell_intensity': tempDataSorting[:,z,46],
                'Avg_'+bacteriaCategory[1]+'_whole_cell_intensity': tempDataSorting[:,z,47],
                'Std_'+bacteriaCategory[1]+'_whole_cell_intensity': tempDataSorting[:,z,48],
                'Avg_'+bacteriaCategory[2]+'_whole_cell_intensity': tempDataSorting[:,z,49],
                'Std_'+bacteriaCategory[2]+'_whole_cell_intensity': tempDataSorting[:,z,50],
                
                'Avg_All_DNA_objects_in_cells_intensity': tempDataSorting[:,z,51],
                'Std_All_DNA_objects_in_cells_intensity': tempDataSorting[:,z,52],
                'Avg_'+bacteriaCategory[0]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,z,53],
                'Std_'+bacteriaCategory[0]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,z,54],
                'Avg_'+bacteriaCategory[1]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,z,55],
                'Std_'+bacteriaCategory[1]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,z,56],
                'Avg_'+bacteriaCategory[2]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,z,57],
                'Std_'+bacteriaCategory[2]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,z,58],
                
                'Avg_All_Intensity_mass_displacement_in_cells': tempDataSorting[:,z,59],
                'Std_All_Intensity_mass_displacement_in_cells': tempDataSorting[:,z,60],
                'Avg_'+bacteriaCategory[0]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,z,61],
                'Std_'+bacteriaCategory[0]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,z,62],
                'Avg_'+bacteriaCategory[1]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,z,63],
                'Std_'+bacteriaCategory[1]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,z,64],
                'Avg_'+bacteriaCategory[2]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,z,65],
                'Std_'+bacteriaCategory[2]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,z,66],
                
                'Avg_All_MajorAxis_Length_of_bacteria': tempDataSorting[:,z,67],
                'Std_All_MajorAxis_Length_of_bacteria': tempDataSorting[:,z,68],
                'Avg_'+bacteriaCategory[0]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,z,69],
                'Std_'+bacteriaCategory[0]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,z,70],
                'Avg_'+bacteriaCategory[1]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,z,71],
                'Std_'+bacteriaCategory[1]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,z,72],
                'Avg_'+bacteriaCategory[2]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,z,73],
                'Std_'+bacteriaCategory[2]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,z,74],
                
                'Avg_All_MaxFeret_Length_of_bacteria': tempDataSorting[:,z,75],
                'Std_All_MaxFeret_Length_of_bacteria': tempDataSorting[:,z,76],
                'Avg_'+bacteriaCategory[0]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,z,77],
                'Std_'+bacteriaCategory[0]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,z,78],
                'Avg_'+bacteriaCategory[1]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,z,79],
                'Std_'+bacteriaCategory[1]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,z,80],
                'Avg_'+bacteriaCategory[2]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,z,81],
                'Std_'+bacteriaCategory[2]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,z,82],
                
                'Avg_All_Tot_DNA_Distribution_1of6': tempDataSorting[:,z,83],
                'Std_All_Tot_DNA_Distribution_1of6': tempDataSorting[:,z,84],
                'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,z,85],
                'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,z,86],
                'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,z,87],
                'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,z,88],
                'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,z,89],
                'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,z,90],
                
                'Avg_All_Tot_DNA_Distribution_2of6': tempDataSorting[:,z,91],
                'Std_All_Tot_DNA_Distribution_2of6': tempDataSorting[:,z,92],
                'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,z,93],
                'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,z,94],
                'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,z,95],
                'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,z,96],
                'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,z,97],
                'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,z,98],
                
                'Avg_All_Tot_DNA_Distribution_3of6': tempDataSorting[:,z,99],
                'Std_All_Tot_DNA_Distribution_3of6': tempDataSorting[:,z,100],
                'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,z,101],
                'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,z,102],
                'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,z,103],
                'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,z,104],
                'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,z,105],
                'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,z,106],
                
                'Avg_All_Tot_DNA_Distribution_4of6': tempDataSorting[:,z,107],
                'Std_All_Tot_DNA_Distribution_4of6': tempDataSorting[:,z,108],
                'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,z,109],
                'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,z,110],
                'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,z,111],
                'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,z,112],
                'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,z,113],
                'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,z,114],
                
                'Avg_All_Tot_DNA_Distribution_5of6': tempDataSorting[:,z,115],
                'Std_All_Tot_DNA_Distribution_5of6': tempDataSorting[:,z,116],
                'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,z,117],
                'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,z,118],
                'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,z,119],
                'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,z,120],
                'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,z,121],
                'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,z,122],
                
                'Avg_All_Tot_DNA_Distribution_6of6': tempDataSorting[:,z,123],
                'Std_All_Tot_DNA_Distribution_6of6': tempDataSorting[:,z,124],
                'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,z,125],
                'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,z,126],
                'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,z,127],
                'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,z,128],
                'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,z,129],
                'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,z,130],
                
                'Avg_All_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,131],
                'Std_All_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,132], 
                'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,133],
                'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,134], 
                'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,135],
                'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,136], 
                'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,137],
                'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,138], 
                
                'Avg_All_Conc_DNA_Distribution_1of6': tempDataSorting[:,z,139],
                'Std_All_Conc_DNA_Distribution_1of6': tempDataSorting[:,z,140],
                'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,z,141],
                'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,z,142],
                'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,z,143],
                'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,z,144],
                'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,z,145],
                'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,z,146],
                
                'Avg_All_Conc_DNA_Distribution_2of6': tempDataSorting[:,z,147],
                'Std_All_Conc_DNA_Distribution_2of6': tempDataSorting[:,z,148],
                'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,z,149],
                'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,z,150],
                'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,z,151],
                'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,z,152],
                'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,z,153],
                'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,z,154],
                
                'Avg_All_Conc_DNA_Distribution_3of6': tempDataSorting[:,z,155],
                'Std_All_Conc_DNA_Distribution_3of6': tempDataSorting[:,z,156],
                'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,z,157],
                'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,z,158],
                'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,z,159],
                'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,z,160],
                'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,z,161],
                'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,z,162],
                
                'Avg_All_Conc_DNA_Distribution_4of6': tempDataSorting[:,z,163],
                'Std_All_Conc_DNA_Distribution_4of6': tempDataSorting[:,z,164],
                'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,z,165],
                'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,z,166],
                'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,z,167],
                'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,z,168],
                'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,z,169],
                'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,z,170],
                
                'Avg_All_Conc_DNA_Distribution_5of6': tempDataSorting[:,z,171],
                'Std_All_Conc_DNA_Distribution_5of6': tempDataSorting[:,z,172],
                'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,z,173],
                'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,z,174],
                'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,z,175],
                'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,z,176],
                'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,z,177],
                'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,z,178],
                
                'Avg_All_Conc_DNA_Distribution_6of6': tempDataSorting[:,z,179],
                'Std_All_Conc_DNA_Distribution_6of6': tempDataSorting[:,z,180],
                'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,z,181],
                'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,z,182],
                'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,z,183],
                'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,z,184],
                'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,z,185],
                'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,z,186],
                
                'Avg_All_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,187],
                'Std_All_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,188],
                'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,189],
                'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,190],
                'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,191],
                'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,192],
                'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,193],
                'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,z,194]})
            SiteExists[z] = True
        else:
            pass


    """Create dataframe from combined results of all sites per well"""
    ExcelSheetsforPrint[5] = pd.DataFrame({'Manually_discard_data_from_analyses': False,
        'Plate_date': plate_date_ordered, 
        'Plate_number': plate_number_ordered,
        'Well_no': well_number_ordered,
        'Gene': genes_ordered,
        'Treatment': treatment_ordered,
        'Site_no': site_number_combined,
        'Total_Bacteria_count': tempDataSorting[:,5,1],
        'Total_'+bacteriaCategory[0]+'_cells': tempDataSorting[:,5,2],
        'Avg_'+bacteriaCategory[0]+'_p(Enrichment_score)': tempDataSorting[:,5,3],
        'Avg_'+bacteriaCategory[0]+'_logit_Enrichment_score': tempDataSorting[:,5,4],
        'Total_'+bacteriaCategory[1]+'_cells': tempDataSorting[:,5,5],
        'Avg_'+bacteriaCategory[1]+'_p(Enrichment_score)': tempDataSorting[:,5,6],
        'Avg_'+bacteriaCategory[1]+'_logit_Enrichment_score': tempDataSorting[:,5,7],
        'Total_'+bacteriaCategory[2]+'_cells': tempDataSorting[:,5,8],
        'Avg_'+bacteriaCategory[2]+'_p(Enrichment_score)': tempDataSorting[:,5,9],
        'Avg_'+bacteriaCategory[2]+'_logit_Enrichment_score': tempDataSorting[:,5,10],
        
        'Avg_All_focis_per_cell': tempDataSorting[:,5,11],
        'Std_All_focis_per_cell': tempDataSorting[:,5,12],
        'Avg_'+bacteriaCategory[0]+'_focis_per_cell': tempDataSorting[:,5,13],
        'Std_'+bacteriaCategory[0]+'_focis_per_cell': tempDataSorting[:,5,14],
        'Avg_'+bacteriaCategory[1]+'_focis_per_cell': tempDataSorting[:,5,15],
        'Std_'+bacteriaCategory[1]+'_focis_per_cell': tempDataSorting[:,5,16],
        'Avg_'+bacteriaCategory[2]+'_focis_per_cell': tempDataSorting[:,5,17],
        'Std_'+bacteriaCategory[2]+'_focis_per_cell': tempDataSorting[:,5,18],

        'Avg_All_DNA_area_coverage': tempDataSorting[:,5,19],
        'Std_All_DNA_area_coverage': tempDataSorting[:,5,20],
        'Avg_'+bacteriaCategory[0]+'_DNA_area_coverage': tempDataSorting[:,5,21],
        'Std_'+bacteriaCategory[0]+'_DNA_area_coverage': tempDataSorting[:,5,22],
        'Avg_'+bacteriaCategory[1]+'_DNA_area_coverage': tempDataSorting[:,5,23],
        'Std_'+bacteriaCategory[1]+'_DNA_area_coverage': tempDataSorting[:,5,24],
        'Avg_'+bacteriaCategory[2]+'_DNA_area_coverage': tempDataSorting[:,5,25],
        'Std_'+bacteriaCategory[2]+'_DNA_area_coverage': tempDataSorting[:,5,26],
        
        'Avg_All_DNA_compaction': tempDataSorting[:,5,27],
        'Std_All_DNA_compaction': tempDataSorting[:,5,28],
        'Avg_'+bacteriaCategory[0]+'_DNA_compaction': tempDataSorting[:,5,29],
        'Std_'+bacteriaCategory[0]+'_DNA_compaction': tempDataSorting[:,5,30],
        'Avg_'+bacteriaCategory[1]+'_DNA_compaction': tempDataSorting[:,5,31],
        'Std_'+bacteriaCategory[1]+'_DNA_compaction': tempDataSorting[:,5,32],
        'Avg_'+bacteriaCategory[2]+'_DNA_compaction': tempDataSorting[:,5,33],
        'Std_'+bacteriaCategory[2]+'_DNA_compaction': tempDataSorting[:,5,34],
        
        'Avg_All_DNA_eccentricity': tempDataSorting[:,5,35],
        'Std_All_DNA_eccentricity': tempDataSorting[:,5,36],
        'Avg_'+bacteriaCategory[0]+'_DNA_eccentricity': tempDataSorting[:,5,37],
        'Std_'+bacteriaCategory[0]+'_DNA_eccentricity': tempDataSorting[:,5,38],
        'Avg_'+bacteriaCategory[1]+'_DNA_eccentricity': tempDataSorting[:,5,39],
        'Std_'+bacteriaCategory[1]+'_DNA_eccentricity': tempDataSorting[:,5,40],
        'Avg_'+bacteriaCategory[2]+'_DNA_eccentricity': tempDataSorting[:,5,41],
        'Std_'+bacteriaCategory[2]+'_DNA_eccentricity': tempDataSorting[:,5,42],
        
        'Avg_All_whole_cell_intensity': tempDataSorting[:,5,43],
        'Std_All_whole_cell_intensity': tempDataSorting[:,5,44],
        'Avg_'+bacteriaCategory[0]+'_whole_cell_intensity': tempDataSorting[:,5,45],
        'Std_'+bacteriaCategory[0]+'_whole_cell_intensity': tempDataSorting[:,5,46],
        'Avg_'+bacteriaCategory[1]+'_whole_cell_intensity': tempDataSorting[:,5,47],
        'Std_'+bacteriaCategory[1]+'_whole_cell_intensity': tempDataSorting[:,5,48],
        'Avg_'+bacteriaCategory[2]+'_whole_cell_intensity': tempDataSorting[:,5,49],
        'Std_'+bacteriaCategory[2]+'_whole_cell_intensity': tempDataSorting[:,5,50],
        
        'Avg_All_DNA_objects_in_cells_intensity': tempDataSorting[:,5,51],
        'Std_All_DNA_objects_in_cells_intensity': tempDataSorting[:,5,52],
        'Avg_'+bacteriaCategory[0]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,5,53],
        'Std_'+bacteriaCategory[0]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,5,54],
        'Avg_'+bacteriaCategory[1]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,5,55],
        'Std_'+bacteriaCategory[1]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,5,56],
        'Avg_'+bacteriaCategory[2]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,5,57],
        'Std_'+bacteriaCategory[2]+'_DNA_objects_in_cells_intensity': tempDataSorting[:,5,58],
        
        'Avg_All_Intensity_mass_displacement_in_cells': tempDataSorting[:,5,59],
        'Std_All_Intensity_mass_displacement_in_cells': tempDataSorting[:,5,60],
        'Avg_'+bacteriaCategory[0]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,5,61],
        'Std_'+bacteriaCategory[0]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,5,62],
        'Avg_'+bacteriaCategory[1]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,5,63],
        'Std_'+bacteriaCategory[1]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,5,64],
        'Avg_'+bacteriaCategory[2]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,5,65],
        'Std_'+bacteriaCategory[2]+'_Intensity_mass_displacement_in_cells': tempDataSorting[:,5,66],
        
        'Avg_All_MajorAxis_Length_of_bacteria': tempDataSorting[:,5,67],
        'Std_All_MajorAxis_Length_of_bacteria': tempDataSorting[:,5,68],
        'Avg_'+bacteriaCategory[0]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,5,69],
        'Std_'+bacteriaCategory[0]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,5,70],
        'Avg_'+bacteriaCategory[1]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,5,71],
        'Std_'+bacteriaCategory[1]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,5,72],
        'Avg_'+bacteriaCategory[2]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,5,73],
        'Std_'+bacteriaCategory[2]+'_MajorAxis_Length_of_bacteria': tempDataSorting[:,5,74],
        
        'Avg_All_MaxFeret_Length_of_bacteria': tempDataSorting[:,5,75],
        'Std_All_MaxFeret_Length_of_bacteria': tempDataSorting[:,5,76],
        'Avg_'+bacteriaCategory[0]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,5,77],
        'Std_'+bacteriaCategory[0]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,5,78],
        'Avg_'+bacteriaCategory[1]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,5,79],
        'Std_'+bacteriaCategory[1]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,5,80],
        'Avg_'+bacteriaCategory[2]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,5,81],
        'Std_'+bacteriaCategory[2]+'_MaxFeret_Length_of_bacteria': tempDataSorting[:,5,82],
        
        'Avg_All_Tot_DNA_Distribution_1of6': tempDataSorting[:,5,83],
        'Std_All_Tot_DNA_Distribution_1of6': tempDataSorting[:,5,84],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,5,85],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,5,86],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,5,87],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,5,88],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,5,89],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_1of6': tempDataSorting[:,5,90],
        
        'Avg_All_Tot_DNA_Distribution_2of6': tempDataSorting[:,5,91],
        'Std_All_Tot_DNA_Distribution_2of6': tempDataSorting[:,5,92],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,5,93],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,5,94],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,5,95],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,5,96],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,5,97],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_2of6': tempDataSorting[:,5,98],
        
        'Avg_All_Tot_DNA_Distribution_3of6': tempDataSorting[:,5,99],
        'Std_All_Tot_DNA_Distribution_3of6': tempDataSorting[:,5,100],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,5,101],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,5,102],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,5,103],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,5,104],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,5,105],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_3of6': tempDataSorting[:,5,106],
        
        'Avg_All_Tot_DNA_Distribution_4of6': tempDataSorting[:,5,107],
        'Std_All_Tot_DNA_Distribution_4of6': tempDataSorting[:,5,108],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,5,109],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,5,110],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,5,111],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,5,112],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,5,113],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_4of6': tempDataSorting[:,5,114],
        
        'Avg_All_Tot_DNA_Distribution_5of6': tempDataSorting[:,5,115],
        'Std_All_Tot_DNA_Distribution_5of6': tempDataSorting[:,5,116],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,5,117],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,5,118],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,5,119],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,5,120],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,5,121],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_5of6': tempDataSorting[:,5,122],
        
        'Avg_All_Tot_DNA_Distribution_6of6': tempDataSorting[:,5,123],
        'Std_All_Tot_DNA_Distribution_6of6': tempDataSorting[:,5,124],
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,5,125],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,5,126],
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,5,127],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,5,128],
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,5,129],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_6of6': tempDataSorting[:,5,130],
        
        'Avg_All_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,131],
        'Std_All_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,132], 
        'Avg_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,133],
        'Std_'+bacteriaCategory[0]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,134], 
        'Avg_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,135],
        'Std_'+bacteriaCategory[1]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,136], 
        'Avg_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,137],
        'Std_'+bacteriaCategory[2]+'_Tot_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,138], 
        
        'Avg_All_Conc_DNA_Distribution_1of6': tempDataSorting[:,5,139],
        'Std_All_Conc_DNA_Distribution_1of6': tempDataSorting[:,5,140],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,5,141],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,5,142],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,5,143],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,5,144],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,5,145],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_1of6': tempDataSorting[:,5,146],
        
        'Avg_All_Conc_DNA_Distribution_2of6': tempDataSorting[:,5,147],
        'Std_All_Conc_DNA_Distribution_2of6': tempDataSorting[:,5,148],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,5,149],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,5,150],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,5,151],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,5,152],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,5,153],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_2of6': tempDataSorting[:,5,154],
        
        'Avg_All_Conc_DNA_Distribution_3of6': tempDataSorting[:,5,155],
        'Std_All_Conc_DNA_Distribution_3of6': tempDataSorting[:,5,156],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,5,157],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,5,158],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,5,159],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,5,160],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,5,161],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_3of6': tempDataSorting[:,5,162],
        
        'Avg_All_Conc_DNA_Distribution_4of6': tempDataSorting[:,5,163],
        'Std_All_Conc_DNA_Distribution_4of6': tempDataSorting[:,5,164],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,5,165],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,5,166],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,5,167],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,5,168],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,5,169],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_4of6': tempDataSorting[:,5,170],
        
        'Avg_All_Conc_DNA_Distribution_5of6': tempDataSorting[:,5,171],
        'Std_All_Conc_DNA_Distribution_5of6': tempDataSorting[:,5,172],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,5,173],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,5,174],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,5,175],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,5,176],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,5,177],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_5of6': tempDataSorting[:,5,178],
        
        'Avg_All_Conc_DNA_Distribution_6of6': tempDataSorting[:,5,179],
        'Std_All_Conc_DNA_Distribution_6of6': tempDataSorting[:,5,180],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,5,181],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,5,182],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,5,183],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,5,184],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,5,185],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_6of6': tempDataSorting[:,5,186],
        
        'Avg_All_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,187],
        'Std_All_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,188],
        'Avg_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,189],
        'Std_'+bacteriaCategory[0]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,190],
        'Avg_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,191],
        'Std_'+bacteriaCategory[1]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,192],
        'Avg_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,193],
        'Std_'+bacteriaCategory[2]+'_Conc_DNA_Distribution_Extra_Overflow': tempDataSorting[:,5,194]})
    SiteExists[5] = True
    

    ################  EXPORT DATA TO EXCEL FILE  ################

    with pd.ExcelWriter(output_path) as writer: # pylint: disable=abstract-class-instantiated
        if SiteExists[5]:
            ExcelSheetsforPrint[5].to_excel(writer, sheet_name='Combined_well_results')
        if SiteExists[0]:
            ExcelSheetsforPrint[0].to_excel(writer, sheet_name='Site_distribution_results')
        if SiteExists[1]:
            ExcelSheetsforPrint[1].to_excel(writer, sheet_name='Site_1')
        if SiteExists[2]:
            ExcelSheetsforPrint[2].to_excel(writer, sheet_name='Site_2')
        if SiteExists[3]:
            ExcelSheetsforPrint[3].to_excel(writer, sheet_name='Site_3')
        if SiteExists[4]:
            ExcelSheetsforPrint[4].to_excel(writer, sheet_name='Site_4')

    

    print("\nFinished analyses of all images found in files with prefix: " + input_filename_prefix + "\nExported results to excel file: " + output_filename + "\nLocation: " + absParentDirectorypath + "\n")



################  FILL IN INFORMATION FOR YOUR RUN BELOW THIS LINE  ################

"""USE FUNCTION AnalyzeDNAinBacteria()"""
"""ENTER CORRECT INPUT VARIABLES (PREFIX OF INPUT FILENAME, OUTPUT FILENAME, OUTPUT FOLDER PATH, cluster=False)"""
"""EXAMPLE: AnalyzeDNAinBacteria("200501_F_RowRes_", "200501_F_RowAnalyzed.xlsx", "./AnalyzedPlates/200501_Plate5/output-data", cluster=False)"""
"""MAKE SURE ALL OF YOUR INPUT FILES ARE IN THE SAME FOLDER (PATH)"""
"""REMEMBER TO INCLUDE WellGeneIndexList.xlsx IN THE SAME FOLDER (PATH)"""

AnalyzeDNAinBacteria("<INPUT FILENAME PREFIX>","<OUTPUT FILENAME>", "<OUTPUT FOLDER PATH>", cluster=False)
