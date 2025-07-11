# Exploring the genetic landscape of DNA supercompaction

This repository contains all pipelines, scripts, and reference documents used for the machine learning-assisted image analysis of images from the high-throughput screening in the publication below. This includes files used for batch processing on high-performance computers. Additionally, the repository contains scripts for analyzing classification parameter importance, and files from statistical analyses in R. Scripts and templates for image preprocessing and analysis with Coli-Inspector and MicrobeJ (for analyzing DNA profile widths and creating kymograph heat maps) are available in the [Zenodo repository](https://doi.org/10.5281/zenodo.14063054) of our previous paper [RecN and RecA orchestrate an ordered DNA supercompaction response following ciprofloxacin-induced DNA damage in <i>Escherichia coli</i>](https://doi.org/10.1093/nar/gkaf437).

## Publication

Vikedal K, Berges N, Riisnæs IMM, Ræder SB, Bjørnholt JV, Bjørås M, Skarstad K, Helgesen E, and Booth JA <br>
**Exploring the genetic landscape of ciprofloxacin-induced DNA supercompaction in <i>Escherichia coli</i>**<br>
bioRxiv (2025). doi: [10.1101/2025.07.xx.xxxxxx](https://doi.org/10.1101/2025.07.xx.xxxxxx)

Images and data from high-throughput imaging are available in the BioImage Archive under accession number S-BIADXXXX at [https://doi.org/10.6019/S-BIADXXXX](https://doi.org/10.6019/S-BIADXXXX).

Please see the Material and methods section and the Supplementary Material of the paper for details on screening procedure and image analysis. 

---

## Table of Contents

1. [CellProfiler pipelines](#cellprofiler-pipelines)
2. [HPC scripts](#hpc-scripts)
3. [CellProfiler output handling scripts](#cellprofiler-output-handling-scripts)
4. [Analyzing importance of classification parameters](#analyzing-importance-of-classification-parameters)
5. [R-based statistical analyses](#r-based-statistical-analyses)
6. [License Details](#license-details)
7. [Author](#author)
8. [Acknowledgements](#acknowledgements)

---

## CellProfiler pipelines

The main CellProfiler pipeline used for analysis of images from the screening, as well as pipelines used for training of the Single-Cell and Phenotype models, are located in the [CellProfiler_Pipelines](CellProfiler_Pipelines/) folder. The main CellProfiler pipeline is set up for batch processing of the screening images on a high-performance computer (cluster computer).

**Batch processing workflow:**
1. Place raw images in an `input-images` folder within the plate directory.
2. Run the pipeline ([CellProfiler_Main_KV_ScreeningPaper.cppipe](CellProfiler_Pipelines/CellProfiler_Main_KV_ScreeningPaper.cppipe)) locally to determine which images to include in the analysis and generate a `Batch_data.h5` batch file. 
3. Run CellProfiler in batch mode on the high-performance computer using the generated batch file.
4. Per-cell measurement tables (`.csv`) with data from the CellProfiler analysis of every well will be written to an `output-data` folder in the same location as the `input-images` folder.  

There are three separate files outputed for each well:
- `*_ImageMetadata.csv`: contain metadata for full images.
- `*_ObjectMeasurements.csv`: include measurements of cell and DNA features for all segmented objects prior to filtration with the Single-Cell model.
- `*_InterestingObjectMeasurements.csv`: contain measurements of cell and DNA features as well as phenotype scoring (Phenotype model) for cells remaining after filtration with the Single-Cell model. 

### Sorting rules (Single-Cell and Phenotype model)

The [sorting-rules](CellProfiler_Pipelines/sorting-rules/) folder contains the rules of the: 
- Single-Cell model ([InterestingCellFilter_MainAnalysis.txt](CellProfiler_Pipelines/sorting-rules/InterestingCellFilter_MainAnalysis.txt)) used to exclude irrelevant objects
- Phenotype model ([CompactionCategoryFilter_MainAnalysis.txt](CellProfiler_Pipelines/sorting-rules/CompactionCategoryFilter_MainAnalysis.txt)) used to classify DNA compaction phenotypes

This folder should always be included within the `input-images` directory to ensure filtering is handled correctly in the CellProfiler pipelines. 

### CellProfiler Analyst properties 
The [CPA_properties](CellProfiler_Pipelines/CPA_properties/) folder contains `.properties` files used for training the Single-Cell model (`InterestingCellFilter`) and Phenotype model (`CompactionCategoryFilter`) in CellProfiler Analyst.  The `classifier_ignore_columns` setting lists parameters excluded from training to avoid confounders (e.g. image metadata). 

## HPC scripts

The batch processing capabilities of CellProfiler were employed on high-performance computers (HPCs) to automate the large scale analysis of images from the screening. After creating a batch file in CellProfiler locally, all files needed for image analysis were uploaded to the HPC (plate directory containing `input-images` (with `sorting-rules`) and `output-data` (with `Batch_data.h5`) folders). SLURM scripts were used to manage the large-scale batch processing:

- **Environment setup:** Ensure that Anaconda or Miniconda is installed on the HPC, and that CellProfiler is available. Use [environment.yml](HPC_scripts/environment.yml) to create a Conda environment with CellProfiler and necessary dependecies. Note: this environment worked for CellProfiler v4.1.3 - newer versions may require other dependencies.
- **Job submission:** The [HPC_job-script.slurm](HPC_scripts/HPC_job-script.slurm) script was used to run the analysis. Note: edit all parameters annotated with `#CHANGE` before reusing the script.
- **Output sorting:** Analysis will produce a separate output file (`.csv`) for every well from the analyzed plate. To sort the output files into folders based on their plate row letters, we used the script [sorting-script.sh](HPC_scripts/sorting-script.sh).


## CellProfiler output handling scripts

After CellProfiler generates per-cell measurement tables for each well, a series of custom Python scripts (and functions) were used to process and summarize these data:
1. Concatenate data for all wells within a given plate row:
    - Function: `combineRowResults(<...>)` from ([CombinePlateResults.py](CellProfiler_OutputHandling/CombinePlateResults.py)).
    - Output: all three `.csv` files outputed from CellProfiler for wells within a plate row are concatenated into three files, prefixed by `<plate date>_<row>_RowRes_`.
2. Compile and summarize all per-cell measurements for every plate to get per-well results:
    - Main script: [CombinePlateResults.py](CellProfiler_OutputHandling/CombinePlateResults.py) 
    - Supporting scripts: 
        - [DNACompactionAnalysis.py](CellProfiler_OutputHandling/DNACompactionAnalysis.py)
        - [polyafit.py](CellProfiler_OutputHandling/polyafit.py)
        - [dirichletintegrate.py](CellProfiler_OutputHandling/dirichletintegrate.py) 
        - [hypergeom.py](CellProfiler_OutputHandling/hypergeom.py) 
    - Output: All concatenated row-level `.csv` files are analyzed by the [DNACompactionAnalysis.py](CellProfiler_OutputHandling/DNACompactionAnalysis.py) script to produce row-level `*_RowAnalyzed.xlsx` files with summarized results for each well and descriptive statistics. These row-level files are then combined into plate-level results-files titled `<plate date>_FullPlateAnalyzed.xlsx`.
    - Requires that the [PlateDateNumberIndexList.xlsx](CellProfiler_OutputHandling/PlateDateNumberIndexList.xlsx) and [WellGeneIndexList.xlsx](CellProfiler_OutputHandling/WellGeneIndexList.xlsx) files are up-to-date in the working folder.
3. Create Hit list which summarize results from all replicates and calculate enrichment scores for DNA supercompaction phenoytpes in all strains: 
    - Main script: [CombineEntireHitList.py](CellProfiler_OutputHandling/CombineEntireHitList.py)
    - Supporting scripts: 
        - [polyafit.py](CellProfiler_OutputHandling/polyafit.py)
        - [dirichletintegrate.py](CellProfiler_OutputHandling/dirichletintegrate.py) 
        - [hypergeom.py](CellProfiler_OutputHandling/hypergeom.py) 
    - Output: Generates the `MainHitList_AllPlatesCombined.xlsx` file, containing enrichment scores and hit status for every strain across plates and replicates, by combining and analyzing results from a folder of all `*_FullPlateAnalyzed.xlsx` files.
    - Requires that the [PlateDateNumberIndexList.xlsx](CellProfiler_OutputHandling/PlateDateNumberIndexList.xlsx) and [WellGeneIndexList.xlsx](CellProfiler_OutputHandling/WellGeneIndexList.xlsx) files are up-to-date in the working folder.

**Note:** The scripts assume the index lists ([PlateDateNumberIndexList.xlsx](CellProfiler_OutputHandling/PlateDateNumberIndexList.xlsx) and [WellGeneIndexList.xlsx](CellProfiler_OutputHandling/WellGeneIndexList.xlsx)) are complete. If using custom plate layouts, update the index files accordingly. 

## Analyzing importance of classification parameters

The [AnalyzeClassificationRules.py](AnalyzeClassificationRules/AnalyzeClassificationRules.py) script was used to calculate the importance of parameters and features for: 
- Single-Cell model ([InterestingCellFilter_MainAnalysis.txt](AnalyzeClassificationRules/InterestingCellFilter_MainAnalysis.txt))
- Phenotype model ([CompactionCategoryFilter_MainAnalysis.txt](AnalyzeClassificationRules/CompactionCategoryFilter_MainAnalysis.txt)).

## R-based statistical analyses

All scripts and data used to generate the mixed-effects model results for time- and dose-dependent survival data after CIP exposure and UV irradiation are available in the [R_StatisticalAnalyses](R_StatisticalAnalyses) folder. The rendered `.html` reports document the full statistical analyses (including our rationale for using only main effects), and the `.Rmd` files together with the accompanying `.txt` data files allow full replication of the results.

---

## License Details

This repository is released under the MIT License. See [LICENSE](LICENSE) for details. 

Enrichment scoring methods in the Python scripts are derived from the source code of CellProfiler Analyst, which is licensed under the BSD 3-Clause License. This applies to parts of [CombineEntireHitList.py](CellProfiler_OutputHandling/CombineEntireHitList.py) and [DNACompactionAnalysis.py](CellProfiler_OutputHandling/DNACompactionAnalysis.py). Additionally, the scripts [polyafit.py](CellProfiler_OutputHandling/polyafit.py), [dirichletintegrate.py](CellProfiler_OutputHandling/dirichletintegrate.py), and [hypergeom.py](CellProfiler_OutputHandling/hypergeom.py) were copied from the open-source code. The relevant license text is included in [License_BSD3](License_BSD3). 

## Author

Krister Vikedal

## Acknowledgements

Enrichment scoring methods in the Python scripts are derived from the source code of [CellProfiler Analyst](https://github.com/CellProfiler). CellProfiler and CellProfiler Analyst can be downloaded from [cellprofiler.org](https://cellprofiler.org/) and [cellprofileranalyst.org](https://cellprofileranalyst.org/), respectively. 

We thank Jon Kristen Lærdahl (Oslo University Hospital) for assistance with initial CellProfiler pipeline development. Screening image analyses were performed on resources provided by Sigma2 – the National Infrastructure for High-Performance Computing and Data Storage in Norway (projects nn5014k and nn9383k), as well as resources from the high-performance computing infrastructure at the University of Oslo (project ec100). We also acknowledge the use of ChatGPT by OpenAI for code suggestions in some of the scripts included in this repository. 
