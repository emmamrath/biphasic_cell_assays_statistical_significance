Biphasic Cell Aassays Statistical Significance
==============================================

## Overview

These R programs process cell death assay data to determine if there is a statistically significant difference between 2 cell lines.
The data may be biphasic with two inflexion points for LC50 (lethal concentration 50%) or uniphasic (one inflexion point, only one TC50).

The programs are:

* step01_choose_model.R (determine whether the data for a cell line is biphase or uniphasic)
* step02_determine_statistical_significance.R (using models chosen manually in step01, determine whether difference between TC50 of 2 cell lines is statistically significant)

The code in these programs were used to analyse cell death assays of the following paper, which can be considered the citation for these programs:  

BAMLET kills chemotherapy-resistant mesothelioma cells, holding oleic acid in an activated cytotoxic state.  
Rath EM, Cheng YY, Pinese M, Sarun KH, Hudson AL, Weir C, Wang YD, Håkansson AP, Howell VM, Liu GJ, Reid G, Knott RB, Duff AP, Church WB.  
PLoS One. 2018 Aug 29;13(8):e0203003. doi: 10.1371/journal.pone.0203003. eCollection 2018.  
PMID: 30157247 PMCID: [PMC6114908](https://www.ncbi.nlm.nih.gov/pubmed/30157247) DOI: [10.1371/journal.pone.0203003](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0203003)

## Input data format

Input data is tab-delimited, and may be the spreadsheet output of a plate reader. 
Columns are experiment details (cell type and experiment date) and increasing concentrations of treatment with the first as control.
Rows are experiments and technical replicates of experiments.

```
Organism        Line    Experiment      0.009765625     0.01953125      0.0390625       0.078125        0.15625 0.3125  0.625   1.25    2.5     5
Human   HOFM01  20151227                1.13298265      1.2518457       1.169650548     1.179535294     0.321131209     0.522045855     0.466633854     0.13501292      0.102600386
Human   HOFM01  20151227                1.133447664     1.202355076     1.130446487     1.148294016     0.386334178     0.472610525     0.404720123     0.111811275     0.104074022
Human   HOFM01  20151227                1.086991591     1.110302273     1.14235446      1.144038233     0.413532698     0.484604061     0.341553911     0.118448912     0.100895324
Human   HOFM02  20151227                1.3875474       1.4874276       1.4437425       1.4470510       0.4163004       0.5711646       0.4543922       0.1397323       0.1107171
Human   HOFM02  20151227                1.3402566       1.5280446       1.4796335       1.4625205       0.6526509       0.5887666       0.5701049       0.1751256       0.1129981
Human   HOFM02  20151227                1.3554021       1.4123655       1.4460303       1.3760835       0.4947850       0.5863324       0.5364160       0.1555323       0.1118969
Human   MM05    20151227                0.9306402       1.072803606     1.000322454     0.971798196     1.015381044     0.051270145     0.061395193     0.04015194      0.03506362
Human   MM05    20151227                0.932502068     0.98302691      0.950200394     0.917405687     0.986093263     0.05081112      0.055105287     0.04245181      0.035288504
Human   MM05    20151227                0.877522567     0.948420315     0.908850384     0.870947033     0.897884675     0.062435153     0.053259234     0.03896685      0.03922624
Human   HOFM01  20160101                        0.967966769     0.989835659     1.065876337     1.030052222     0.360388112     0.469571801     0.081465117     0.083591331
Human   HOFM01  20160101                        0.99887248      0.951469447     0.95428038      1.003324873     0.429527386     0.440047408     0.095592708     0.083782593
Human   HOFM01  20160101                        0.952357466     0.992662706     0.928494636     0.959528507     0.395327186     0.475710454     0.084401046     0.090186647
Human   HOFM02  20160101                        1.02839295      1.077715684     1.099926299     1.068857198     0.433192847     0.549805729     0.10366146      0.108999449
Human   HOFM02  20160101                        0.952455612     0.927007554     1.016159378     1.00877303      0.472461466     0.473694846     0.108300583     0.106056805
Human   HOFM02  20160101                        0.989872251     0.999210255     1.009557377     1.067603636     0.510467778     0.526365052     0.107632008     0.104853276
Human   MM05    20160101                        0.987954557     0.957794263     0.952959485     0.941152669     0.036229704     0.029511854     0.028344659     0.034761633
Human   MM05    20160101                        0.965598989     0.911933621     0.937712089     0.894310905     0.041093448     0.03371918      0.025956793     0.036116079
Human   MM05    20160101                        0.97185231      0.910456751     0.887674473     0.847337187     0.040304132     0.04323273      0.034800579     0.038541448
Human   H226    20160101                        0.869709604     0.848789789     0.76639692      0.84190623      0.060172979     0.044480836     0.058133912     0.104220456
Human   H226    20160101                        0.899923501     0.950099544     0.937156406     0.855187946     0.055260169     0.046137252     0.053929744     0.068013247
Human   H226    20160101                        0.932106895     0.944990504     0.90207903      0.943474785     0.064801593     0.04243088      0.046187829     0.067370921
```

## Input parameters

The programs can be called as: 
```
Rscript <program_name>.R <arg_1> <arg_2> <more args...>
```
The arguments of each program can be seen at the top of the program.

#### step01_choose_model.R
Rscript step01_choose_model.R <input_file> <cell_line> > <output_file>
Rscript step01_choose_model.R human_bamlet_cell_assay_data.txt HOFM01 > step01_output.txt
args_infile = args[1] # human_bamlet_cell_assay_data.txt
args_cellline = args[2] # HOFM01

#### step02_determine_statistical_significance.R
Rscript step02_determine_statistical_significance.R <input_file> <cell_line_1> <model_1> <cell_line_2> <model_2> > <output_file>
Rscript step02_determine_statistical_significance.R human_bamlet_cell_assay_data.txt MM05 1component HOFM01 2component > step02_output.txt
args_infile = args[1] # human_bamlet_cell_assay_data.txt
args_cellline_1 = args[2] # MM05
args_model_1 = args[2] # 1component # MM05 is uniphasic, only 1 TC50
args_cellline_2 = args[2] # HOFM01
args_model_2 = args[2] # 2component # HOFM01 is biphasic, has 2 TC50s

