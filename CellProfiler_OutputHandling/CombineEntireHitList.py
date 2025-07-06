import numpy as np
import pandas as pd
import os
import time
import dirichletintegrate
import polyafit
from scipy import stats

"""GO TO THE LAST LINES OF THIS DOCUMENT TO ADD INFORMATION FOR SCRIPT"""

def getListOfFiles(directorypath, indexList):
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
                allFiles = allFiles + getListOfFiles(fullPath, indexList)
        else:
            "Check if entry belongs to a plate in the index list and is a Full plate analysis"
            entryInfo = entry.split("_")
            if (entryInfo[0] in indexList) and (entryInfo[1] == "FullPlateAnalyzed.xlsx"):
                allFiles.append(fullPath)
            else:
                continue 
    return allFiles

def combineAllPlatesToHitList(main_input_filepath="."):
    """Make one output file containing classification and enrichment scores of all wells 
        from all plates in the files of the input path. 

    The main_input_filepath should be the folder containing all the FullPlateAnalyzed excel documents."""

    ##Write in the correct categories used
    bacteriaCategory = ['PositiveTreated', 'Nontreated', 'NegativeTreated']

    ##Check if output file already exists in folder
    output_filename = "HitList_AllPlatesCombined.xlsx"
    output_filepath = os.path.join(main_input_filepath, output_filename)
    
    if os.path.exists(output_filepath):
        raise Exception("A file with the name " + output_filename + " already exists in the given filepath. Rename the old file or this output file.")
    else: 
        pass

    
    ##Load plate date index to evaluate if given file belongs to a known plate
    index_input_filepath = "."
    wellGene_path = os.path.join(index_input_filepath, "WellGeneIndexList.xlsx")
    plateDate_path = os.path.join(index_input_filepath, "PlateDateNumberIndexList.xlsx")
    if not os.path.exists(wellGene_path):
        raise Exception("The file WellGeneIndexList.xlsx is missing in path. Please add it to path, it is necessary to run the code.")
    elif not os.path.exists(plateDate_path):
        raise Exception("The file PlateDateNumberIndexList.xlsx is missing in path. Please add it to path, it is necessary to run the code.")
    else:
        pass
    plateDateNumberIndexList = pd.read_excel(plateDate_path, header=[0, 0])
    plate_date_index = plateDateNumberIndexList.loc[:, "Plate_date"]
    plate_date_index = plate_date_index["Plate_date"].tolist()
    for dates in range(len(plate_date_index)):
        plate_date_index[dates] = str(plate_date_index[dates])
    
    replicate_no_index = plateDateNumberIndexList.loc[:, "Replicate_no"]
    replicate_no_index = replicate_no_index["Replicate_no"].tolist()

    wellGeneIndexList = pd.read_excel(wellGene_path, header=[0, 1])
    well_number_index = wellGeneIndexList.loc[:, ("Plate1", "Well")]
    well_number_index = pd.Series.tolist(well_number_index)


    ##Find all files inside main input filepath directory that match correct format and make list
    fileList = getListOfFiles(main_input_filepath, plate_date_index)
    combinedAllPlatesExists = False 
    # colNames = ['Plate_date', 'Plate_number', 'Well_no', 'Gene', 'Total_Bacteria_count', 'Total_'+bacteriaCategory[0]+'_cells', 'Total_'+bacteriaCategory[1]+'_cells', 'Total_'+bacteriaCategory[2]+'_cells']


    
    for eachFile in fileList:
        filePath = os.path.join(main_input_filepath, eachFile)
        if combinedAllPlatesExists:
            fileContent = pd.read_excel(filePath, sheet_name="Combined_well_results", index_col=0) #, header=[0, 0], usecols=colNames)
            frames = [combinedAllPlates,fileContent]
            combinedAllPlates = pd.concat(frames, axis=0, ignore_index=True)
        else: 
            combinedAllPlates = pd.read_excel(filePath, sheet_name="Combined_well_results", index_col=0) #, header=[0, 1], usecols=colNames)
            combinedAllPlatesExists = True

    columnLoadingVariable = None
    cellCountperCategory = []

    gene_metadata = combinedAllPlates.loc[:,"Gene"]
    # gene_metadata = gene_metadata["Gene"].tolist()
    gene_metadata = pd.Series.tolist(gene_metadata)

    plate_date_metadata = combinedAllPlates.loc[:,"Plate_date"]
    # plate_date_metadata = plate_date_metadata["Plate_date"].tolist()
    plate_date_metadata = pd.Series.tolist(plate_date_metadata)


    replicate_no_metadata = []
    
    for dates in range(len(plate_date_metadata)):
        plate_date_metadata[dates] = str(plate_date_metadata[dates])
        for indexes in range(len(plate_date_index)):
            if plate_date_metadata[dates] == plate_date_index[indexes]:
                replicate_no_metadata.append(replicate_no_index[indexes])
            else:
                pass

    plate_number_metadata = combinedAllPlates.loc[:,"Plate_number"]
    # plate_number_metadata = plate_number_metadata["Plate_number"].tolist()
    plate_number_metadata = pd.Series.tolist(plate_number_metadata)

    well_number_metadata = combinedAllPlates.loc[:,"Well_no"]
    # well_number_metadata = well_number_metadata["Well_no"].tolist()
    well_number_metadata = pd.Series.tolist(well_number_metadata) 

    total_bacteria_count = combinedAllPlates.loc[:,"Total_Bacteria_count"]
    # total_bacteria_count = total_bacteria_count["Total_Bacteria_count"].tolist()
    total_bacteria_count = pd.Series.tolist(total_bacteria_count)
    
    for eachCategory in bacteriaCategory:
        columnLoadingVariable = combinedAllPlates.loc[:, 'Total_'+eachCategory+'_cells']
        # cellCountperCategory.append(columnLoadingVariable['Total_'+eachCategory+'_cells'].tolist())
        cellCountperCategory.append(pd.Series.tolist(columnLoadingVariable))
        columnLoadingVariable = None

    # probScorePerCategory = [ [] for _ in range(len(bacteriaCategory)) ]
    # logitScorePerCategory = [ [] for _ in range(len(bacteriaCategory)) ]

    probScorePerCategory = np.zeros((len(bacteriaCategory),len(gene_metadata)))
    logitScorePerCategory = np.zeros((len(bacteriaCategory),len(gene_metadata)))

    two_classes = len(bacteriaCategory) == 2

    counts_AllPlates = [ [] for _ in range(len(bacteriaCategory)) ]
    for eachCategory in range(0,len(bacteriaCategory)):
        counts_AllPlates[eachCategory] += cellCountperCategory[eachCategory]
    counts_AllPlates = np.transpose(counts_AllPlates)
    alpha_AllPlates, converged_AllPlates = polyafit.fit_betabinom_minka_alternating(counts_AllPlates)
    
    if not converged_AllPlates:
        raise Exception("The Pólya fit estimation did not converge for the combined cell counts in wells")
    else:
        pass

    for eachSample in range(0,len(gene_metadata)):
        counts = []
        for eachCategory in range(0,len(bacteriaCategory)):
            if np.isnan(cellCountperCategory[eachCategory][eachSample]):
                counts.append(0)
            else:
                counts.append(cellCountperCategory[eachCategory][eachSample])
        
        scores = []
        """Determine Enrichment scores (probability) for well and/or site"""
        scores = np.array(dirichletintegrate.score(alpha_AllPlates, np.array(counts)))

        """clamp scores to [0,1]"""
        scores[scores > 1.] = 1.
        scores[scores < 0.] = 0.

        """Add Enrichment score (probability) to data table"""
        for eachCategory in range(0,len(bacteriaCategory)):
            probScorePerCategory[eachCategory,eachSample] = scores[eachCategory]
            
        """Calculate logit Enrichment scores and add to data table"""
        """Special case: only calculate logit of "positives" for 2-classes"""
        if two_classes:
            for eachCategory in range(0,len(bacteriaCategory)):
                logitScorePerCategory[eachCategory,eachSample] = [np.log10(scores[0]) - (np.log10(1 - scores[0]))]
        else:
            for eachCategory in range(0,len(bacteriaCategory)):
                x = [np.log10(scores[eachCategory]) - (np.log10(1 - scores[eachCategory]))]
                logitScorePerCategory[eachCategory,eachSample] += x


    ##Make list of hit indicators
    hit_indicator_perSample = []
    
    for gene in range(0,len(gene_metadata)):
        if probScorePerCategory[2][gene] > 0.75:
            hit_indicator_perSample.append("YES!")
        elif probScorePerCategory[1][gene] > 0.75:   
            #logitScorePerCategory[1][gene] > 2.5:
            hit_indicator_perSample.append("Interesting")
        elif probScorePerCategory[0][gene] > 0.75:
            hit_indicator_perSample.append("NO!")
        else:
            hit_indicator_perSample.append("Unclear")

        

    ##Make list of control indicators

    def match_control_index(control_genes, plate_date):
        """Find the index of foci belonging to the parent cells"""
        index_list = []
        for i in range(0,len(gene_metadata)):
            if (plate_date_metadata[i] == plate_date) and (gene_metadata[i] in control_genes):
                index_list.append(i) 
            else:
                pass           
        return index_list
    
    def findAverage(input_list):
        """Find average and standard deviation of nonzero values in list"""
        nonNANlistvalues = [] 
        allValuesAreNAN = True   
        for i in range(len(input_list)):
            if np.isnan(input_list[i]):
                pass
            # elif input_list[i] > 0:
            #     if allValuesAreNAN:
            #         allValuesAreNAN = False
            #     nonNANlistvalues.append(input_list[i])
            else:
                if allValuesAreNAN:
                    allValuesAreNAN = False
                nonNANlistvalues.append(input_list[i])
        if allValuesAreNAN:
            avg = np.nan
            # sd = np.nan
            sem = np.nan
        elif len(nonNANlistvalues) > 1 :
            avg = np.average(nonNANlistvalues)
            # sd = np.std(nonNANlistvalues)
            sem = stats.sem(nonNANlistvalues)
        elif len(nonNANlistvalues) == 1:
            avg = np.average(nonNANlistvalues)
            # sd = np.std(nonNANlistvalues)
            sem = 0
        else:
            avg = 0
            # sd = 0
            sem = 0
        return avg, sem
    
    
    control_indexlist = []
    control_genes = ["wt NT", "wt +cip", "recn NT", "recn +cip"]
    control_indicator_perPlate = [ [ False for _ in range(len(plate_date_index)) ] for _ in range(len(control_genes)) ]
    control_indicator_list = [ False for _ in range(len(gene_metadata)) ]
    #controlNumberCount = 0

    

    for plateIndex, eachPlate in enumerate(plate_date_index):
        #Get index values of control samples and store in variable for wt NT, recn NT, recn +cip and wt +cip
        control_indexlist = match_control_index(control_genes, eachPlate)
        
        #multipleControls = [ False for _ in range(len(control_genes)) ]
        #Change multipleControls of control gene index to True if more than one
        # for eachControlGene in range(0,control_genes):
        #     for indexValue in control_indexlist:
        #         if gene_metadata[indexValue] == control_genes[eachControlGene]:
        #             controlNumberCount = controlNumberCount + 1
        #     if controlNumberCount > 1:
        #         multipleControls[eachControlGene] = True

        # for eachControlGene in range(0,control_genes):
        #     if multipleControls[eachControlGene]:
        #         for indexValue in control_genes[]
        #     else:
        
        for indexValue in control_indexlist:
            if gene_metadata[indexValue] == "wt NT":
                #if logitScorePerCategory[1][indexValue] > 2.5:
                if probScorePerCategory[1][indexValue] > 0.75:
                    control_indicator_perPlate[0][plateIndex] = True
                else:
                    pass
            elif gene_metadata[indexValue] == "wt +cip":
                if probScorePerCategory[0][indexValue] > 0.75:
                    control_indicator_perPlate[1][plateIndex] = True
                else:
                    pass
            elif gene_metadata[indexValue] == "recn NT":
                #if logitScorePerCategory[1][indexValue] > 2.5:
                if probScorePerCategory[1][indexValue] > 0.75:
                    control_indicator_perPlate[2][plateIndex] = True
                else:
                    pass
            elif gene_metadata[indexValue] == "recn +cip":
                if probScorePerCategory[2][indexValue] > 0.75:
                    control_indicator_perPlate[3][plateIndex] = True
                else:
                    pass
            else:
                pass

        if ((control_indicator_perPlate[1][plateIndex]) and (control_indicator_perPlate[3][plateIndex])) and ((control_indicator_perPlate[0][plateIndex]) or (control_indicator_perPlate[2][plateIndex])):
            for i in range(0,len(plate_date_metadata)):
                if plate_date_metadata[i] == plate_date_index[plateIndex]:
                    control_indicator_list[i] = True
                else:
                    pass
        else:
            pass

    

    replicate_sortedIndexes = [ [] for _ in range(3) ]

    gene_metadata_perGene = []
    control_indicator_perGene = []
    plate_number_metadata_perGene = []
    well_number_metadata_perGene = []
    avg_bacteria_count_perGene = []
    sem_bacteria_count_perGene = []
    AvgProbScorePerCategory = [ [] for _ in range(len(bacteriaCategory)) ]
    AvgLogitScorePerCategory = [ [] for _ in range(len(bacteriaCategory)) ]
    SEMProbScorePerCategory = [ [] for _ in range(len(bacteriaCategory)) ]
    SEMLogitScorePerCategory = [ [] for _ in range(len(bacteriaCategory)) ]

    hit_indicator_perGene = []


    for plate_number in range(1,20):
        for sortedWells in well_number_index:
            for well_number in range(0,len(well_number_metadata)):
                if (plate_number_metadata[well_number] == plate_number) and (well_number_metadata[well_number] == sortedWells):
                    replicate_sortedIndexes[replicate_no_metadata[well_number] - 1].append(well_number)
                else:
                    pass
            while not (len(replicate_sortedIndexes[0]) == len(replicate_sortedIndexes[1]) == len(replicate_sortedIndexes[2])):
                if (len(replicate_sortedIndexes[0]) > len(replicate_sortedIndexes[1])):
                    replicate_sortedIndexes[1].append(np.nan)
                elif (len(replicate_sortedIndexes[0]) > len(replicate_sortedIndexes[2])):
                    replicate_sortedIndexes[2].append(np.nan)
                elif (len(replicate_sortedIndexes[1]) > len(replicate_sortedIndexes[0])):
                    replicate_sortedIndexes[0].append(np.nan)
                elif (len(replicate_sortedIndexes[1]) > len(replicate_sortedIndexes[2])):
                    replicate_sortedIndexes[2].append(np.nan)
                elif (len(replicate_sortedIndexes[2]) > len(replicate_sortedIndexes[0])):
                    replicate_sortedIndexes[0].append(np.nan)
                elif (len(replicate_sortedIndexes[2]) > len(replicate_sortedIndexes[1])):
                    replicate_sortedIndexes[1].append(np.nan)
                


    if not (len(replicate_sortedIndexes[0]) == len(replicate_sortedIndexes[1]) == len(replicate_sortedIndexes[2])):
        raise Exception("Issues arrised with sorting the indexes of the gene replicates")        

    for eachGene in range(len(replicate_sortedIndexes[0])):
        indexWithAllGenes = None
        if not np.isnan(replicate_sortedIndexes[0][eachGene]):
            indexWithAllGenes = 0
        elif not np.isnan(replicate_sortedIndexes[1][eachGene]):
            indexWithAllGenes = 1
        elif not np.isnan(replicate_sortedIndexes[2][eachGene]):
            indexWithAllGenes = 2
        
        gene_metadata_perGene.append(gene_metadata[replicate_sortedIndexes[indexWithAllGenes][eachGene]])
        plate_number_metadata_perGene.append(plate_number_metadata[replicate_sortedIndexes[indexWithAllGenes][eachGene]])
        well_number_metadata_perGene.append(well_number_metadata[replicate_sortedIndexes[indexWithAllGenes][eachGene]])
        
        clearYes = 0
        unclearYes = 0
        clearNo = 0
        others = 0

        gene_indexList = [replicate_sortedIndexes[0][eachGene], replicate_sortedIndexes[1][eachGene], replicate_sortedIndexes[2][eachGene]]

        if (len(gene_indexList) < 3) or (np.isnan(replicate_sortedIndexes[0][eachGene]) or np.isnan(replicate_sortedIndexes[1][eachGene]) or np.isnan(replicate_sortedIndexes[2][eachGene])):
            hit_indicator_perGene.append("Too few replicates")
        else:
            for gene_index in gene_indexList:
                if hit_indicator_perSample[gene_index] == "YES!":
                    clearYes = clearYes + 1
                elif hit_indicator_perSample[gene_index] == "Interesting":
                    unclearYes = unclearYes + 1
                elif hit_indicator_perSample[gene_index] == "NO!":
                    clearNo = clearNo + 1
                elif hit_indicator_perSample[gene_index] == "Unclear":
                    others = others + 1
                else:
                    raise Exception('Something went wrong with assigning hit indicator')
            
            if clearYes == len(gene_indexList):
                hit_indicator_perGene.append("YES! Clear")
            elif clearYes >= 1 and clearNo == 0:
                hit_indicator_perGene.append("Yes, unclear")
            elif unclearYes == len(gene_indexList):
                hit_indicator_perGene.append("INTERESTING! Clear")
            elif unclearYes >= 1 and clearNo == 0:
                hit_indicator_perGene.append("Interesting, unclear")
            elif clearNo == len(gene_indexList):
                hit_indicator_perGene.append("NO! Clear")
            elif clearNo >= 1:
                hit_indicator_perGene.append("No, unclear")
            else:
                hit_indicator_perGene.append("Unclear")
        
        allControlsOK = True

        for gene_index in gene_indexList:
            if np.isnan(gene_index):
                pass
            elif control_indicator_list[gene_index] == False:
                allControlsOK = False
            else:
                pass


        if allControlsOK:
            control_indicator_perGene.append(True)
        else:
            control_indicator_perGene.append(False)

        totalBacCountList = []
        
        probScoreList = [ [] for _ in range(len(bacteriaCategory)) ]
        logitScoreList = [ [] for _ in range(len(bacteriaCategory)) ]

        for gene_index in gene_indexList:
            if np.isnan(gene_index):
                totalBacCountList.append(np.nan)
            else:
                totalBacCountList.append(total_bacteria_count[gene_index])

            for eachCategory in range(0,len(bacteriaCategory)):
                if np.isnan(gene_index):
                    probScoreList[eachCategory].append(np.nan)
                    logitScoreList[eachCategory].append(np.nan)
                else:
                    probScoreList[eachCategory].append(probScorePerCategory[eachCategory][gene_index])
                    logitScoreList[eachCategory].append(logitScorePerCategory[eachCategory][gene_index])

        avg_totcount, sem_totcount = findAverage(totalBacCountList)
        avg_bacteria_count_perGene.append(avg_totcount) 
        sem_bacteria_count_perGene.append(sem_totcount)

        for eachCategory in range(0,len(bacteriaCategory)):
            avg_probscore, sem_probscore = findAverage(probScoreList[eachCategory])
            avg_logitscore, sem_logitscore = findAverage(logitScoreList[eachCategory])
            AvgProbScorePerCategory[eachCategory].append(avg_probscore)
            AvgLogitScorePerCategory[eachCategory].append(avg_logitscore)
            SEMProbScorePerCategory[eachCategory].append(sem_probscore)
            SEMLogitScorePerCategory[eachCategory].append(sem_logitscore)

    
    # def match_gene_index(gene, plate_number):
    #     """Find the index of foci belonging to the parent cells"""
    #     index_list = []
    #     for i in range(0,len(gene_metadata)):
    #         if (plate_number_metadata[i] == plate_number) and (gene_metadata[i] == gene):
    #             index_list.append(i) 
    #         else:
    #             pass           
    #     return index_list
    
    # hit_indicator_perGene = []        
    # # everySingleGene = []
    
    # for plate_number in range(1,14):
    #     everySingleGene = []
    #     for g in range(0,len(gene_metadata)):
    #         if (gene_metadata[g] not in everySingleGene) and (plate_number_metadata[g] == plate_number):
    #             everySingleGene.append(gene_metadata[g])

    #     for gene in range(0,len(everySingleGene)):
    #         clearYes = 0
    #         unclearYes = 0
    #         clearNo = 0
    #         others = 0

    #         gene_indexList = match_gene_index(everySingleGene[gene], plate_number)
            
    #         if len(gene_indexList) < 3:
    #             hit_indicator_perGene.append("Too few replicates")
    #         else:
    #             for gene_index in gene_indexList:
    #                 if hit_indicator_perSample[gene_index] == "YES!":
    #                     clearYes = clearYes + 1
    #                 elif hit_indicator_perSample[gene_index] == "Interesting":
    #                     unclearYes = unclearYes + 1
    #                 elif hit_indicator_perSample[gene_index] == "NO!":
    #                     clearNo = clearNo + 1
    #                 elif hit_indicator_perSample[gene_index] == "Unclear":
    #                     others = others + 1
    #                 else:
    #                     raise Exception('Something went wrong with assigning hit indicator')
                
    #             if clearYes == len(gene_indexList):
    #                 hit_indicator_perGene.append("YES! Clear")
    #             elif clearYes >= 1 and clearNo == 0:
    #                 hit_indicator_perGene.append("Yes, unclear")
    #             elif unclearYes == len(gene_indexList):
    #                 hit_indicator_perGene.append("INTERESTING! Clear")
    #             elif unclearYes >= 1 and clearNo == 0:
    #                 hit_indicator_perGene.append("Interesting, unclear")
    #             elif clearNo == len(gene_indexList):
    #                 hit_indicator_perGene.append("NO! Clear")
    #             elif clearNo >= 1:
    #                 hit_indicator_perGene.append("No, unclear")
    #             else:
    #                 hit_indicator_perGene.append("Unclear")

    ###Fix sorting of the values above

    # gene_metadata_perGene = [ [] for _ in range(3) ]
    # for plate_number in range(1,14):
    #     for well_number in range(0,len(well_number_metadata)):
    #         for accountedGenes in range(len(gene_metadata_perGene[0])):
    #             if gene_metadata_perGene[accountedGenes] == gene_metadata[well_number]:
                    
    #                 gene_metadata_perGene[0].append(gene_metadata[well_number])
    #                 gene_metadata_perGene[1].append(plate_number)
    #                 gene_metadata_perGene[2].append(well_number_metadata[well_number])
    
    
    # for g in range(0,len(gene_metadata)):
    #     if (gene_metadata[g] not in gene_metadata_perGene) and (plate_number_metadata[g] == plate_number):
    #         gene_metadata_perGene.append(gene_metadata[g])

    
    
    
    perGeneMatrixforPrint = pd.DataFrame({'Manually_discard_data_from_analyses': False,
        'Gene': gene_metadata_perGene,
        'Hit_indicator': hit_indicator_perGene,
        'Control_indicator': control_indicator_perGene,
        'Plate_number': plate_number_metadata_perGene,
        'Well_no': well_number_metadata_perGene,
        'AvgTotal_Bacteria_count': avg_bacteria_count_perGene,
        'SEMTotal_Bacteria_count': sem_bacteria_count_perGene,
        'Avg_'+bacteriaCategory[0]+'_p(Enrichment_score)': AvgProbScorePerCategory[0],
        'Avg_'+bacteriaCategory[0]+'_logit_Enrichment_score': AvgLogitScorePerCategory[0],
        'Avg_'+bacteriaCategory[1]+'_p(Enrichment_score)': AvgProbScorePerCategory[1],
        'Avg_'+bacteriaCategory[1]+'_logit_Enrichment_score': AvgLogitScorePerCategory[1],
        'Avg_'+bacteriaCategory[2]+'_p(Enrichment_score)': AvgProbScorePerCategory[2],
        'Avg_'+bacteriaCategory[2]+'_logit_Enrichment_score': AvgLogitScorePerCategory[2],
        'SEM_'+bacteriaCategory[0]+'_p(Enrichment_score)': SEMProbScorePerCategory[0],
        'SEM_'+bacteriaCategory[0]+'_logit_Enrichment_score': SEMLogitScorePerCategory[0],
        'SEM_'+bacteriaCategory[1]+'_p(Enrichment_score)': SEMProbScorePerCategory[1],
        'SEM_'+bacteriaCategory[1]+'_logit_Enrichment_score': SEMLogitScorePerCategory[1],
        'SEM_'+bacteriaCategory[2]+'_p(Enrichment_score)': SEMProbScorePerCategory[2],
        'SEM_'+bacteriaCategory[2]+'_logit_Enrichment_score': SEMLogitScorePerCategory[2]})

    perReplicateMatrixforPrint = pd.DataFrame({'Manually_discard_data_from_analyses': False,
        'Gene': gene_metadata,
        'Hit_indicator': hit_indicator_perSample,
        'Control_indicator': control_indicator_list,
        'Replicate_no': replicate_no_metadata,
        'Plate_date': plate_date_metadata, 
        'Plate_number': plate_number_metadata,
        'Well_no': well_number_metadata,
        'Total_Bacteria_count': total_bacteria_count,
        bacteriaCategory[0]+'_cells': cellCountperCategory[0],
        bacteriaCategory[0]+'_p(Enrichment_score)': probScorePerCategory[0],
        bacteriaCategory[0]+'_logit_Enrichment_score': logitScorePerCategory[0],
        bacteriaCategory[1]+'_cells': cellCountperCategory[1],
        bacteriaCategory[1]+'_p(Enrichment_score)': probScorePerCategory[1],
        bacteriaCategory[1]+'_logit_Enrichment_score': logitScorePerCategory[1],
        bacteriaCategory[2]+'_cells': cellCountperCategory[2],
        bacteriaCategory[2]+'_p(Enrichment_score)': probScorePerCategory[2],
        bacteriaCategory[2]+'_logit_Enrichment_score': logitScorePerCategory[2]})

    # perGeneMatrixforPrint.sort_values(by=['Gene', 'Plate_number'])
    # perReplicateMatrixforPrint.sort_values(by=['Gene', 'Plate_date'])


    

    with pd.ExcelWriter(output_filepath) as writer: # pylint: disable=abstract-class-instantiated
        perGeneMatrixforPrint.to_excel(writer, sheet_name='Hit_list_perGeneAveraged')
        perReplicateMatrixforPrint.to_excel(writer, sheet_name='Hit_list_per_replicate')

    #perGeneMatrixforPrint.to_excel(writer, sheet_name='Hit_list_per_geneAveraged')


    duration = time.process_time()
    print("\nFinished creating hit list for all plate analysis files in location: " + main_input_filepath + ".\nThe file containing the combined hit list (" + output_filename + ") can be found in the same location.\n Time used: " + str(duration) + ".\n")


################  FILL IN INFORMATION FOR YOUR RUN BELOW THIS LINE  ################

###FIX BELOW THIS LINE
"""USE FUNCTION combineAllPlatesToHitList()"""
"""THIS FUNCTION WILL COMBINE RESULTS FROM MULTIPLE PLATES ANALYZED USING the CombinePlateResults.py SCRIPT"""
"""ENTER CORRECT INPUT VARIABLE (PLATE RESULTS FOLDER PATH)"""
"""EXAMPLE: combineAllPlatesToHitList("./AnalysisResults")"""
"""MAKE SURE ALL OF YOUR INPUT FILES ARE IN THE SAME FOLDER (PATH)"""
"""REMEMBER TO INCLUDE DNACompactionAnalysis.py, dirichletintegrate.py, polyafit.py, hypergeom.py, PlateDateNumberIndexList.xlsx and WellGeneIndexList.xlsx IN THE SAME FOLDER (PATH)"""

combineAllPlatesToHitList("<PLATE RESULTS FOLDER PATH>")