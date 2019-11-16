# Rscript step02_determine_statistical_significance.R <input_file> <cell_line_1> <model_1> <cell_line_2> <model_2> > <output_file>
# Rscript step02_determine_statistical_significance.R human_bamlet_cell_assay_data.txt MM05 1component HOFM01 2component > step02_output.txt

library(reshape2)
library(ggplot2)
library(optimx)
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'human_bamlet_cell_assay_data.txt', 'MM05', 'HOFM01' ) # for testing

# Format of tab-delimited input data:
# Organism        Line    Experiment      0.009765625     0.01953125      0.0390625       0.078125        0.15625 0.3125  0.625   1.25    2.5     5
# Human   HOFM01  20151227                1.13298265      1.2518457       1.169650548     1.179535294     0.321131209     0.522045855     0.466633854     0.13501292      0.102600386
# Human   HOFM01  20151227                1.133447664     1.202355076     1.130446487     1.148294016     0.386334178     0.472610525     0.404720123     0.111811275     0.104074022
# Human   HOFM01  20151227                1.086991591     1.110302273     1.14235446      1.144038233     0.413532698     0.484604061     0.341553911     0.118448912     0.100895324
# Human   HOFM02  20151227                1.3875474       1.4874276       1.4437425       1.4470510       0.4163004       0.5711646       0.4543922       0.1397323       0.1107171
# Human   HOFM02  20151227                1.3402566       1.5280446       1.4796335       1.4625205       0.6526509       0.5887666       0.5701049       0.1751256       0.1129981
# Human   HOFM02  20151227                1.3554021       1.4123655       1.4460303       1.3760835       0.4947850       0.5863324       0.5364160       0.1555323       0.1118969
# Human   MM05    20151227                0.9306402       1.072803606     1.000322454     0.971798196     1.015381044     0.051270145     0.061395193     0.04015194      0.03506362
# Human   MM05    20151227                0.932502068     0.98302691      0.950200394     0.917405687     0.986093263     0.05081112      0.055105287     0.04245181      0.035288504
# Human   MM05    20151227                0.877522567     0.948420315     0.908850384     0.870947033     0.897884675     0.062435153     0.053259234     0.03896685      0.03922624
# Human   HOFM01  20160101                        0.967966769     0.989835659     1.065876337     1.030052222     0.360388112     0.469571801     0.081465117     0.083591331
# Human   HOFM01  20160101                        0.99887248      0.951469447     0.95428038      1.003324873     0.429527386     0.440047408     0.095592708     0.083782593
# Human   HOFM01  20160101                        0.952357466     0.992662706     0.928494636     0.959528507     0.395327186     0.475710454     0.084401046     0.090186647
# Human   HOFM02  20160101                        1.02839295      1.077715684     1.099926299     1.068857198     0.433192847     0.549805729     0.10366146      0.108999449
# Human   HOFM02  20160101                        0.952455612     0.927007554     1.016159378     1.00877303      0.472461466     0.473694846     0.108300583     0.106056805
# Human   HOFM02  20160101                        0.989872251     0.999210255     1.009557377     1.067603636     0.510467778     0.526365052     0.107632008     0.104853276
# Human   MM05    20160101                        0.987954557     0.957794263     0.952959485     0.941152669     0.036229704     0.029511854     0.028344659     0.034761633
# Human   MM05    20160101                        0.965598989     0.911933621     0.937712089     0.894310905     0.041093448     0.03371918      0.025956793     0.036116079
# Human   MM05    20160101                        0.97185231      0.910456751     0.887674473     0.847337187     0.040304132     0.04323273      0.034800579     0.038541448
# Human   H226    20160101                        0.869709604     0.848789789     0.76639692      0.84190623      0.060172979     0.044480836     0.058133912     0.104220456
# Human   H226    20160101                        0.899923501     0.950099544     0.937156406     0.855187946     0.055260169     0.046137252     0.053929744     0.068013247
# Human   H226    20160101                        0.932106895     0.944990504     0.90207903      0.943474785     0.064801593     0.04243088      0.046187829     0.067370921

args_infile = args[1] # human_bamlet_cell_assay_data.txt
args_cellline_1 = args[2] # MM05
args_model_1 = args[2] # 1component # MM05 is uniphasic, only 1 TC50
args_cellline_2 = args[2] # HOFM01
args_model_2 = args[2] # 2component # HOFM01 is biphasic, has 2 TC50s

human = read.table(args_infile, sep = "\t", header = TRUE)

human = melt(human, id.vars = c("Organism", "Line", "Experiment"), variable.name = "Concentration", value.name = "Viability")
human$Concentration = as.numeric(gsub("^X", "", human$Concentration))

assays = human
rm(human)
assays$Organism = as.factor(assays$Organism)
assays$Line = as.factor(assays$Line)
assays$Experiment = as.factor(as.character(assays$Experiment))
assays = assays[!is.na(assays$Viability),]
assays$log10Conc = log10(assays$Concentration)



####################################################################################################
# Is there a difference between BAMLET TC50 for mesothelioma vs TC50 of the resistant fraction of non-cancer fibroblast cells?
#
# We have already shown that the 1-component model is considered as good as the 2-component model for the following:
# mesothelioma: MM05, H226, VMC23, VMC40; immortalised mesothelial Met5A; non-cancer primary fibroblast: Humofib.
# Thus the simpler 1-component model having less parameters should be used to model those cell lines.
# We have already shown that the 2-component model is statistically significantly better than the 1-component model for the following:
# non-cancer primary fibroblast HOFM01, HOFM02, HOFF.
# Thus the biphasic 2-component model should be used to model those cell lines.
# The 2-component model has two TC50s, one for a more BAMLET resistant fraction of cells, and one for a more BAMLET sensitive fraction of cells.
#
# To determine whether the non-cancer BAMLET-resistant TC50 is statistically significantly higher than the cancer TC50
# we will create a 2-cell-line model containing both the non-cancer 2-component and cancer 1-component models.
# The xmid parameter will be the same for cancer TC50 and non-cancer BAMLET-sensitive TC50.
# Delta2 parameter will be the xshift increase of the non-cancer BAMLET-resistant TC50 compared to the xmid parameter.
# The xmid and delta2 parameters will be varied to find the optimal 2-cell-line model, 
# and then compared with a 2-cell-line model where delta2 is set to zero, which is effectively a 1-component model.
# We already know that the 1-cell-line 2-component model is better than the 1-cell-line 1-component model for HOFM01, HOFM02, and HOFF.
# However, this 2-cell-line model will allow us to test whether delta2 is statistically significantly different from zero.
# That is, to test whether the non-cancer BAMLET-resistant TC50 is statistically significantly different from cancer TC50.
# The 2 models in this test will be called 1component.2component_delta2 and 1component.1component.
#
# These models for determining statistical significance of TC50 differences will be more sophisticated than 
# the previous models for determining whether 1-component or 2-component modelling should be used for a cell line.
# These models will allow variation per experiment in the parameters whose statistical significance is not being tested,
# that is, A and B left and right asymptotes (represent cell viability at lowest and highest doses) and sensitive fraction.
# Sensitive fraction (which is 1 minus resistant fraction) in particular appears to vary per experiment for the non-cancer fibroblasts.

# A = left asymptote
# B = right asymptote
# xmid = inflection point 1 (TC50 of BAMLET-sensitive cells)
# xmid + xshift = inflection point 2 (TC50 of BAMLET-resistant cells)
# scal = factor affecting slope at inflection points
# sens = fraction of cells having inflection point 1 (fraction of BAMLET-sensitive cells)
# 1 - sens = fraction of cells having inflection point 2 (fraction of BAMLET-resistant cells)


####################################################################################################
# Define the various possible models and their functions.
# The main 2 types of models are for biphasic 
# (two LC50, defined as inflection point and shift of second inflexion point with respect to the first)
# and uniphasic (only one inflextion point).
# In case adjustments per experiment need to be made for modelling, define models permitting that.
####################################################################################################

##########
predict.1component5exp.2component3exp_delta2 = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp5_line1 = params[5]
	A_Exp1_line2 = params[6]
	A_Exp2_line2 = params[7]
	A_Exp3_line2 = params[8]
	B_Exp1_line1 = params[9]
	B_Exp2_line1 = params[10]
	B_Exp3_line1 = params[11]
	B_Exp4_line1 = params[12]
	B_Exp5_line1 = params[13]
	B_Exp1_line2 = params[14]
	B_Exp2_line2 = params[15]
	B_Exp3_line2 = params[16]
	scal_line1 = params[17]
	scal_line2 = params[18]
	xmid = params[19]
	sens_Exp1_line2 = params[20]
	sens_Exp2_line2 = params[21]
	sens_Exp3_line2 = params[22]
	delta2 = params[23]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp4_line1 * is_line1 * is_Exp4) + (A_Exp5_line1 * is_line1 * is_Exp5) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2) + (A_Exp3_line2 * is_line2 * is_Exp3)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp4_line1 * is_line1 * is_Exp4) + (B_Exp5_line1 * is_line1 * is_Exp5) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2) + (B_Exp3_line2 * is_line2 * is_Exp3)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalSens = (0 * is_line1) + (sens_Exp1_line2 * is_line2 * is_Exp1) + (sens_Exp2_line2 * is_line2 * is_Exp2) + (sens_Exp3_line2 * is_line2 * is_Exp3)
	finalDelta2 = (0 * is_line1) + (delta2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( finalSens/(1 + exp((xmid - x) / finalScal)) + (1 - finalSens)/(1 + exp((xmid + finalDelta2 - x) / finalScal)) )
	prediction_result
}
loglik.1component5exp.2component3exp_delta2 = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
{
	sigma = params[24]
	yhat = predict.1component5exp.2component3exp_delta2(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component5exp.2component3exp_delta2 = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp5_line1 = params[5]
	A_Exp1_line2 = params[6]
	A_Exp2_line2 = params[7]
	A_Exp3_line2 = params[8]
	B_Exp1_line1 = params[9]
	B_Exp2_line1 = params[10]
	B_Exp3_line1 = params[11]
	B_Exp4_line1 = params[12]
	B_Exp5_line1 = params[13]
	B_Exp1_line2 = params[14]
	B_Exp2_line2 = params[15]
	B_Exp3_line2 = params[16]
	scal_line1 = params[17]
	scal_line2 = params[18]
	xmid = params[19]
	sens_Exp1_line2 = params[20]
	sens_Exp2_line2 = params[21]
	sens_Exp3_line2 = params[22]
	delta2 = params[23]
	sensitive_fraction_average = mean(c(sens_Exp1_line2,sens_Exp2_line2,sens_Exp3_line2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1,A_Exp4_line1,A_Exp5_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2,A_Exp3_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1,B_Exp4_line1,B_Exp5_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2,B_Exp3_line2))
	c(cancer_TC50 = exp(xmid), noncancer_sensitive_fraction_TC50 = exp(xmid), noncancer_resistant_fraction_TC50 = exp(xmid + delta2), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average, A_line1_avg = A_line1_avg, A_line2_avg = A_line2_avg, B_line1_avg = B_line1_avg, B_line2_avg = B_line2_avg)
}

##########
predict.1component4exp.2component3exp_delta2 = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp1_line2 = params[5]
	A_Exp2_line2 = params[6]
	A_Exp3_line2 = params[7]
	B_Exp1_line1 = params[8]
	B_Exp2_line1 = params[9]
	B_Exp3_line1 = params[10]
	B_Exp4_line1 = params[11]
	B_Exp1_line2 = params[12]
	B_Exp2_line2 = params[13]
	B_Exp3_line2 = params[14]
	scal_line1 = params[15]
	scal_line2 = params[16]
	xmid = params[17]
	sens_Exp1_line2 = params[18]
	sens_Exp2_line2 = params[19]
	sens_Exp3_line2 = params[20]
	delta2 = params[21]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp4_line1 * is_line1 * is_Exp4) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2) + (A_Exp3_line2 * is_line2 * is_Exp3)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp4_line1 * is_line1 * is_Exp4) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2) + (B_Exp3_line2 * is_line2 * is_Exp3)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalSens = (0 * is_line1) + (sens_Exp1_line2 * is_line2 * is_Exp1) + (sens_Exp2_line2 * is_line2 * is_Exp2) + (sens_Exp3_line2 * is_line2 * is_Exp3)
	finalDelta2 = (0 * is_line1) + (delta2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( finalSens/(1 + exp((xmid - x) / finalScal)) + (1 - finalSens)/(1 + exp((xmid + finalDelta2 - x) / finalScal)) )
	prediction_result
}
loglik.1component4exp.2component3exp_delta2 = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	sigma = params[22]
	yhat = predict.1component4exp.2component3exp_delta2(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component4exp.2component3exp_delta2 = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp1_line2 = params[5]
	A_Exp2_line2 = params[6]
	A_Exp3_line2 = params[7]
	B_Exp1_line1 = params[8]
	B_Exp2_line1 = params[9]
	B_Exp3_line1 = params[10]
	B_Exp4_line1 = params[11]
	B_Exp1_line2 = params[12]
	B_Exp2_line2 = params[13]
	B_Exp3_line2 = params[14]
	xmid = params[17]
	sens_Exp1_line2 = params[18]
	sens_Exp2_line2 = params[19]
	sens_Exp3_line2 = params[20]
	delta2 = params[21]
	sensitive_fraction_average = mean(c(sens_Exp1_line2,sens_Exp2_line2,sens_Exp3_line2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1,A_Exp4_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2,A_Exp3_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1,B_Exp4_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2,B_Exp3_line2))
	c(cancer_TC50 = exp(xmid), noncancer_sensitive_fraction_TC50 = exp(xmid), noncancer_resistant_fraction_TC50 = exp(xmid + delta2), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average, A_line1_avg = A_line1_avg, A_line2_avg = A_line2_avg, B_line1_avg = B_line1_avg, B_line2_avg = B_line2_avg)
}

##########
predict.1component3exp.2component3exp_delta2 = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp1_line2 = params[4]
	A_Exp2_line2 = params[5]
	A_Exp3_line2 = params[6]
	B_Exp1_line1 = params[7]
	B_Exp2_line1 = params[8]
	B_Exp3_line1 = params[9]
	B_Exp1_line2 = params[10]
	B_Exp2_line2 = params[11]
	B_Exp3_line2 = params[12]
	scal_line1 = params[13]
	scal_line2 = params[14]
	xmid = params[15]
	sens_Exp1_line2 = params[16]
	sens_Exp2_line2 = params[17]
	sens_Exp3_line2 = params[18]
	delta2 = params[19]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2) + (A_Exp3_line2 * is_line2 * is_Exp3)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2) + (B_Exp3_line2 * is_line2 * is_Exp3)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalSens = (0 * is_line1) + (sens_Exp1_line2 * is_line2 * is_Exp1) + (sens_Exp2_line2 * is_line2 * is_Exp2) + (sens_Exp3_line2 * is_line2 * is_Exp3)
	finalDelta2 = (0 * is_line1) + (delta2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( finalSens/(1 + exp((xmid - x) / finalScal)) + (1 - finalSens)/(1 + exp((xmid + finalDelta2 - x) / finalScal)) )
	prediction_result
}
loglik.1component3exp.2component3exp_delta2 = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	sigma = params[20]
	yhat = predict.1component3exp.2component3exp_delta2(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component3exp.2component3exp_delta2 = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp1_line2 = params[4]
	A_Exp2_line2 = params[5]
	A_Exp3_line2 = params[6]
	B_Exp1_line1 = params[7]
	B_Exp2_line1 = params[8]
	B_Exp3_line1 = params[9]
	B_Exp1_line2 = params[10]
	B_Exp2_line2 = params[11]
	B_Exp3_line2 = params[12]
	scal_line1 = params[13]
	scal_line2 = params[14]
	xmid = params[15]
	sens_Exp1_line2 = params[16]
	sens_Exp2_line2 = params[17]
	sens_Exp3_line2 = params[18]
	delta2 = params[19]
	sensitive_fraction_average = mean(c(sens_Exp1_line2,sens_Exp2_line2,sens_Exp3_line2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2,A_Exp3_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2,B_Exp3_line2))
	c(cancer_TC50 = exp(xmid), noncancer_sensitive_fraction_TC50 = exp(xmid), noncancer_resistant_fraction_TC50 = exp(xmid + delta2), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average, A_line1_avg = A_line1_avg, A_line2_avg = A_line2_avg, B_line1_avg = B_line1_avg, B_line2_avg = B_line2_avg)
}

##########
predict.1component2exp.2component3exp_delta2 = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp1_line2 = params[3]
	A_Exp2_line2 = params[4]
	A_Exp3_line2 = params[5]
	B_Exp1_line1 = params[6]
	B_Exp2_line1 = params[7]
	B_Exp1_line2 = params[8]
	B_Exp2_line2 = params[9]
	B_Exp3_line2 = params[10]
	scal_line1 = params[11]
	scal_line2 = params[12]
	xmid = params[13]
	sens_Exp1_line2 = params[14]
	sens_Exp2_line2 = params[15]
	sens_Exp3_line2 = params[16]
	delta2 = params[17]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2) + (A_Exp3_line2 * is_line2 * is_Exp3)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2) + (B_Exp3_line2 * is_line2 * is_Exp3)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalSens = (0 * is_line1) + (sens_Exp1_line2 * is_line2 * is_Exp1) + (sens_Exp2_line2 * is_line2 * is_Exp2) + (sens_Exp3_line2 * is_line2 * is_Exp3)
	finalDelta2 = (0 * is_line1) + (delta2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( finalSens/(1 + exp((xmid - x) / finalScal)) + (1 - finalSens)/(1 + exp((xmid + finalDelta2 - x) / finalScal)) )
	prediction_result
}
loglik.1component2exp.2component3exp_delta2 = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	sigma = params[18]
	yhat = predict.1component2exp.2component3exp_delta2(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component2exp.2component3exp_delta2 = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp1_line2 = params[3]
	A_Exp2_line2 = params[4]
	A_Exp3_line2 = params[5]
	B_Exp1_line1 = params[6]
	B_Exp2_line1 = params[7]
	B_Exp1_line2 = params[8]
	B_Exp2_line2 = params[9]
	B_Exp3_line2 = params[10]
	scal_line1 = params[11]
	scal_line2 = params[12]
	xmid = params[13]
	sens_Exp1_line2 = params[14]
	sens_Exp2_line2 = params[15]
	sens_Exp3_line2 = params[16]
	delta2 = params[17]
	sensitive_fraction_average = mean(c(sens_Exp1_line2,sens_Exp2_line2,sens_Exp3_line2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2,A_Exp3_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2,B_Exp3_line2))
	c(cancer_TC50 = exp(xmid), noncancer_sensitive_fraction_TC50 = exp(xmid), noncancer_resistant_fraction_TC50 = exp(xmid + delta2), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average, A_line1_avg = A_line1_avg, A_line2_avg = A_line2_avg, B_line1_avg = B_line1_avg, B_line2_avg = B_line2_avg)
}

##########
predict.1component5exp.2component2exp_delta2 = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp5_line1 = params[5]
	A_Exp1_line2 = params[6]
	A_Exp2_line2 = params[7]
	B_Exp1_line1 = params[8]
	B_Exp2_line1 = params[9]
	B_Exp3_line1 = params[10]
	B_Exp4_line1 = params[11]
	B_Exp5_line1 = params[12]
	B_Exp1_line2 = params[13]
	B_Exp2_line2 = params[14]
	scal_line1 = params[15]
	scal_line2 = params[16]
	xmid = params[17]
	sens_Exp1_line2 = params[18]
	sens_Exp2_line2 = params[19]
	delta2 = params[20]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp4_line1 * is_line1 * is_Exp4) + (A_Exp5_line1 * is_line1 * is_Exp5) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp4_line1 * is_line1 * is_Exp4) + (B_Exp5_line1 * is_line1 * is_Exp5) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalSens = (0 * is_line1) + (sens_Exp1_line2 * is_line2 * is_Exp1) + (sens_Exp2_line2 * is_line2 * is_Exp2)
	finalDelta2 = (0 * is_line1) + (delta2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( finalSens/(1 + exp((xmid - x) / finalScal)) + (1 - finalSens)/(1 + exp((xmid + finalDelta2 - x) / finalScal)) )
	prediction_result
}
loglik.1component5exp.2component2exp_delta2 = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
{
	sigma = params[21]
	yhat = predict.1component5exp.2component2exp_delta2(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component5exp.2component2exp_delta2 = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp5_line1 = params[5]
	A_Exp1_line2 = params[6]
	A_Exp2_line2 = params[7]
	B_Exp1_line1 = params[8]
	B_Exp2_line1 = params[9]
	B_Exp3_line1 = params[10]
	B_Exp4_line1 = params[11]
	B_Exp5_line1 = params[12]
	B_Exp1_line2 = params[13]
	B_Exp2_line2 = params[14]
	scal_line1 = params[15]
	scal_line2 = params[16]
	xmid = params[17]
	sens_Exp1_line2 = params[18]
	sens_Exp2_line2 = params[19]
	delta2 = params[20]
	sensitive_fraction_average = mean(c(sens_Exp1_line2,sens_Exp2_line2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1,A_Exp4_line1,A_Exp5_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1,B_Exp4_line1,B_Exp5_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2))
	c(cancer_TC50 = exp(xmid), noncancer_sensitive_fraction_TC50 = exp(xmid), noncancer_resistant_fraction_TC50 = exp(xmid + delta2), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average, A_line1_avg = A_line1_avg, A_line2_avg = A_line2_avg, B_line1_avg = B_line1_avg, B_line2_avg = B_line2_avg)
}

##########
predict.1component4exp.2component2exp_delta2 = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp1_line2 = params[5]
	A_Exp2_line2 = params[6]
	B_Exp1_line1 = params[7]
	B_Exp2_line1 = params[8]
	B_Exp3_line1 = params[9]
	B_Exp4_line1 = params[10]
	B_Exp1_line2 = params[11]
	B_Exp2_line2 = params[12]
	scal_line1 = params[13]
	scal_line2 = params[14]
	xmid = params[15]
	sens_Exp1_line2 = params[16]
	sens_Exp2_line2 = params[17]
	delta2 = params[18]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp4_line1 * is_line1 * is_Exp4) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp4_line1 * is_line1 * is_Exp4) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalSens = (0 * is_line1) + (sens_Exp1_line2 * is_line2 * is_Exp1) + (sens_Exp2_line2 * is_line2 * is_Exp2)
	finalDelta2 = (0 * is_line1) + (delta2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( finalSens/(1 + exp((xmid - x) / finalScal)) + (1 - finalSens)/(1 + exp((xmid + finalDelta2 - x) / finalScal)) )
	prediction_result
}
loglik.1component4exp.2component2exp_delta2 = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	sigma = params[19]
	yhat = predict.1component4exp.2component2exp_delta2(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component4exp.2component2exp_delta2 = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp1_line2 = params[5]
	A_Exp2_line2 = params[6]
	B_Exp1_line1 = params[7]
	B_Exp2_line1 = params[8]
	B_Exp3_line1 = params[9]
	B_Exp4_line1 = params[10]
	B_Exp1_line2 = params[11]
	B_Exp2_line2 = params[12]
	scal_line1 = params[13]
	scal_line2 = params[14]
	xmid = params[15]
	sens_Exp1_line2 = params[16]
	sens_Exp2_line2 = params[17]
	delta2 = params[18]
	sensitive_fraction_average = mean(c(sens_Exp1_line2,sens_Exp2_line2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1,A_Exp4_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1,B_Exp4_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2))
	c(cancer_TC50 = exp(xmid), noncancer_sensitive_fraction_TC50 = exp(xmid), noncancer_resistant_fraction_TC50 = exp(xmid + delta2), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average, A_line1_avg = A_line1_avg, A_line2_avg = A_line2_avg, B_line1_avg = B_line1_avg, B_line2_avg = B_line2_avg)
}

##########
predict.1component3exp.2component2exp_delta2 = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp1_line2 = params[4]
	A_Exp2_line2 = params[5]
	B_Exp1_line1 = params[6]
	B_Exp2_line1 = params[7]
	B_Exp3_line1 = params[8]
	B_Exp1_line2 = params[9]
	B_Exp2_line2 = params[10]
	scal_line1 = params[11]
	scal_line2 = params[12]
	xmid = params[13]
	sens_Exp1_line2 = params[14]
	sens_Exp2_line2 = params[15]
	delta2 = params[16]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalSens = (0 * is_line1) + (sens_Exp1_line2 * is_line2 * is_Exp1) + (sens_Exp2_line2 * is_line2 * is_Exp2)
	finalDelta2 = (0 * is_line1) + (delta2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( finalSens/(1 + exp((xmid - x) / finalScal)) + (1 - finalSens)/(1 + exp((xmid + finalDelta2 - x) / finalScal)) )
	prediction_result
}
loglik.1component3exp.2component2exp_delta2 = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	sigma = params[17]
	yhat = predict.1component3exp.2component2exp_delta2(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component3exp.2component2exp_delta2 = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp1_line2 = params[4]
	A_Exp2_line2 = params[5]
	B_Exp1_line1 = params[6]
	B_Exp2_line1 = params[7]
	B_Exp3_line1 = params[8]
	B_Exp1_line2 = params[9]
	B_Exp2_line2 = params[10]
	scal_line1 = params[11]
	scal_line2 = params[12]
	xmid = params[13]
	sens_Exp1_line2 = params[14]
	sens_Exp2_line2 = params[15]
	delta2 = params[16]
	sensitive_fraction_average = mean(c(sens_Exp1_line2,sens_Exp2_line2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2))
	c(cancer_TC50 = exp(xmid), noncancer_sensitive_fraction_TC50 = exp(xmid), noncancer_resistant_fraction_TC50 = exp(xmid + delta2), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average, A_line1_avg = A_line1_avg, A_line2_avg = A_line2_avg, B_line1_avg = B_line1_avg, B_line2_avg = B_line2_avg)
}

##########
predict.1component2exp.2component2exp_delta2 = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp1_line2 = params[3]
	A_Exp2_line2 = params[4]
	B_Exp1_line1 = params[5]
	B_Exp2_line1 = params[6]
	B_Exp1_line2 = params[7]
	B_Exp2_line2 = params[8]
	scal_line1 = params[9]
	scal_line2 = params[10]
	xmid = params[11]
	sens_Exp1_line2 = params[12]
	sens_Exp2_line2 = params[13]
	delta2 = params[14]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalSens = (0 * is_line1) + (sens_Exp1_line2 * is_line2 * is_Exp1) + (sens_Exp2_line2 * is_line2 * is_Exp2)
	finalDelta2 = (0 * is_line1) + (delta2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( finalSens/(1 + exp((xmid - x) / finalScal)) + (1 - finalSens)/(1 + exp((xmid + finalDelta2 - x) / finalScal)) )
	prediction_result
}
loglik.1component2exp.2component2exp_delta2 = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2)
{
	sigma = params[15]
	yhat = predict.1component2exp.2component2exp_delta2(params, x, is_line1, is_line2, is_Exp1, is_Exp2)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component2exp.2component2exp_delta2 = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp1_line2 = params[3]
	A_Exp2_line2 = params[4]
	B_Exp1_line1 = params[5]
	B_Exp2_line1 = params[6]
	B_Exp1_line2 = params[7]
	B_Exp2_line2 = params[8]
	scal_line1 = params[9]
	scal_line2 = params[10]
	xmid = params[11]
	sens_Exp1_line2 = params[12]
	sens_Exp2_line2 = params[13]
	delta2 = params[14]
	sensitive_fraction_average = mean(c(sens_Exp1_line2,sens_Exp2_line2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2))
	c(cancer_TC50 = exp(xmid), noncancer_sensitive_fraction_TC50 = exp(xmid), noncancer_resistant_fraction_TC50 = exp(xmid + delta2), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average, A_line1_avg = A_line1_avg, A_line2_avg = A_line2_avg, B_line1_avg = B_line1_avg, B_line2_avg = B_line2_avg)
}

##########
predict.1component5exp.1component3exp = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp5_line1 = params[5]
	A_Exp1_line2 = params[6]
	A_Exp2_line2 = params[7]
	A_Exp3_line2 = params[8]
	B_Exp1_line1 = params[9]
	B_Exp2_line1 = params[10]
	B_Exp3_line1 = params[11]
	B_Exp4_line1 = params[12]
	B_Exp5_line1 = params[13]
	B_Exp1_line2 = params[14]
	B_Exp2_line2 = params[15]
	B_Exp3_line2 = params[16]
	scal_line1 = params[17]
	scal_line2 = params[18]
	xmid = params[19]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp4_line1 * is_line1 * is_Exp4) + (A_Exp5_line1 * is_line1 * is_Exp5) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2) + (A_Exp3_line2 * is_line2 * is_Exp3)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp4_line1 * is_line1 * is_Exp4) + (B_Exp5_line1 * is_line1 * is_Exp5) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2) + (B_Exp3_line2 * is_line2 * is_Exp3)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.1component5exp.1component3exp = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
{
	sigma = params[20]
	yhat = predict.1component5exp.1component3exp(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component5exp.1component3exp = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp5_line1 = params[5]
	A_Exp1_line2 = params[6]
	A_Exp2_line2 = params[7]
	A_Exp3_line2 = params[8]
	B_Exp1_line1 = params[9]
	B_Exp2_line1 = params[10]
	B_Exp3_line1 = params[11]
	B_Exp4_line1 = params[12]
	B_Exp5_line1 = params[13]
	B_Exp1_line2 = params[14]
	B_Exp2_line2 = params[15]
	B_Exp3_line2 = params[16]
	scal_line1 = params[17]
	scal_line2 = params[18]
	xmid = params[19]
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1,A_Exp4_line1,A_Exp5_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2,A_Exp3_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1,B_Exp4_line1,B_Exp5_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2,B_Exp3_line2))
	c(cancer_TC50 = exp(xmid), noncancer_TC50 = exp(xmid), A_line1_avg=A_line1_avg, A_line2_avg=A_line2_avg, B_line1_avg=B_line1_avg, B_line2_avg=B_line2_avg)
}

##########
predict.1component4exp.1component3exp = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp1_line2 = params[5]
	A_Exp2_line2 = params[6]
	A_Exp3_line2 = params[7]
	B_Exp1_line1 = params[8]
	B_Exp2_line1 = params[9]
	B_Exp3_line1 = params[10]
	B_Exp4_line1 = params[11]
	B_Exp1_line2 = params[12]
	B_Exp2_line2 = params[13]
	B_Exp3_line2 = params[14]
	scal_line1 = params[15]
	scal_line2 = params[16]
	xmid = params[17]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp4_line1 * is_line1 * is_Exp4) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2) + (A_Exp3_line2 * is_line2 * is_Exp3)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp4_line1 * is_line1 * is_Exp4) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2) + (B_Exp3_line2 * is_line2 * is_Exp3)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.1component4exp.1component3exp = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	sigma = params[18]
	yhat = predict.1component4exp.1component3exp(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component4exp.1component3exp = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp1_line2 = params[5]
	A_Exp2_line2 = params[6]
	A_Exp3_line2 = params[7]
	B_Exp1_line1 = params[8]
	B_Exp2_line1 = params[9]
	B_Exp3_line1 = params[10]
	B_Exp4_line1 = params[11]
	B_Exp1_line2 = params[12]
	B_Exp2_line2 = params[13]
	B_Exp3_line2 = params[14]
	scal_line1 = params[15]
	scal_line2 = params[16]
	xmid = params[17]
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1,A_Exp4_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2,A_Exp3_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1,B_Exp4_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2,B_Exp3_line2))
	c(cancer_TC50 = exp(xmid), noncancer_TC50 = exp(xmid), A_line1_avg=A_line1_avg, A_line2_avg=A_line2_avg, B_line1_avg=B_line1_avg, B_line2_avg=B_line2_avg)
}



##########
predict.1component3exp.1component3exp = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp1_line2 = params[4]
	A_Exp2_line2 = params[5]
	A_Exp3_line2 = params[6]
	B_Exp1_line1 = params[7]
	B_Exp2_line1 = params[8]
	B_Exp3_line1 = params[9]
	B_Exp1_line2 = params[10]
	B_Exp2_line2 = params[11]
	B_Exp3_line2 = params[12]
	scal_line1 = params[13]
	scal_line2 = params[14]
	xmid = params[15]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2) + (A_Exp3_line2 * is_line2 * is_Exp3)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2) + (B_Exp3_line2 * is_line2 * is_Exp3)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.1component3exp.1component3exp = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	sigma = params[16]
	yhat = predict.1component3exp.1component3exp(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component3exp.1component3exp = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp1_line2 = params[4]
	A_Exp2_line2 = params[5]
	A_Exp3_line2 = params[6]
	B_Exp1_line1 = params[7]
	B_Exp2_line1 = params[8]
	B_Exp3_line1 = params[9]
	B_Exp1_line2 = params[10]
	B_Exp2_line2 = params[11]
	B_Exp3_line2 = params[12]
	scal_line1 = params[13]
	scal_line2 = params[14]
	xmid = params[15]
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2,A_Exp3_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2,B_Exp3_line2))
	c(cancer_TC50 = exp(xmid), noncancer_TC50 = exp(xmid), A_line1_avg=A_line1_avg, A_line2_avg=A_line2_avg, B_line1_avg=B_line1_avg, B_line2_avg=B_line2_avg)
}

##########
predict.1component2exp.1component3exp = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp1_line2 = params[3]
	A_Exp2_line2 = params[4]
	A_Exp3_line2 = params[5]
	B_Exp1_line1 = params[6]
	B_Exp2_line1 = params[7]
	B_Exp1_line2 = params[8]
	B_Exp2_line2 = params[9]
	B_Exp3_line2 = params[10]
	scal_line1 = params[11]
	scal_line2 = params[12]
	xmid = params[13]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2) + (A_Exp3_line2 * is_line2 * is_Exp3)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2) + (B_Exp3_line2 * is_line2 * is_Exp3)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.1component2exp.1component3exp = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	sigma = params[14]
	yhat = predict.1component2exp.1component3exp(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component2exp.1component3exp = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp1_line2 = params[3]
	A_Exp2_line2 = params[4]
	A_Exp3_line2 = params[5]
	B_Exp1_line1 = params[6]
	B_Exp2_line1 = params[7]
	B_Exp1_line2 = params[8]
	B_Exp2_line2 = params[9]
	B_Exp3_line2 = params[10]
	scal_line1 = params[11]
	scal_line2 = params[12]
	xmid = params[13]
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2,A_Exp3_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2,B_Exp3_line2))
	c(cancer_TC50 = exp(xmid), noncancer_TC50 = exp(xmid), A_line1_avg=A_line1_avg, A_line2_avg=A_line2_avg, B_line1_avg=B_line1_avg, B_line2_avg=B_line2_avg)
}

##########
predict.1component5exp.1component2exp = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp5_line1 = params[5]
	A_Exp1_line2 = params[6]
	A_Exp2_line2 = params[7]
	B_Exp1_line1 = params[8]
	B_Exp2_line1 = params[9]
	B_Exp3_line1 = params[10]
	B_Exp4_line1 = params[11]
	B_Exp5_line1 = params[12]
	B_Exp1_line2 = params[13]
	B_Exp2_line2 = params[14]
	scal_line1 = params[15]
	scal_line2 = params[16]
	xmid = params[17]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp4_line1 * is_line1 * is_Exp4) + (A_Exp5_line1 * is_line1 * is_Exp5) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp4_line1 * is_line1 * is_Exp4) + (B_Exp5_line1 * is_line1 * is_Exp5) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.1component5exp.1component2exp = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
{
	sigma = params[18]
	yhat = predict.1component5exp.1component2exp(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4, is_Exp5)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component5exp.1component2exp = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp5_line1 = params[5]
	A_Exp1_line2 = params[6]
	A_Exp2_line2 = params[7]
	B_Exp1_line1 = params[8]
	B_Exp2_line1 = params[9]
	B_Exp3_line1 = params[10]
	B_Exp4_line1 = params[11]
	B_Exp5_line1 = params[12]
	B_Exp1_line2 = params[13]
	B_Exp2_line2 = params[14]
	scal_line1 = params[15]
	scal_line2 = params[16]
	xmid = params[17]
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1,A_Exp4_line1,A_Exp5_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1,B_Exp4_line1,B_Exp5_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2))
	c(cancer_TC50 = exp(xmid), noncancer_TC50 = exp(xmid), A_line1_avg=A_line1_avg, A_line2_avg=A_line2_avg, B_line1_avg=B_line1_avg, B_line2_avg=B_line2_avg)
}

##########
predict.1component4exp.1component2exp = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp1_line2 = params[5]
	A_Exp2_line2 = params[6]
	B_Exp1_line1 = params[7]
	B_Exp2_line1 = params[8]
	B_Exp3_line1 = params[9]
	B_Exp4_line1 = params[10]
	B_Exp1_line2 = params[11]
	B_Exp2_line2 = params[12]
	scal_line1 = params[13]
	scal_line2 = params[14]
	xmid = params[15]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp4_line1 * is_line1 * is_Exp4) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp4_line1 * is_line1 * is_Exp4) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.1component4exp.1component2exp = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	sigma = params[16]
	yhat = predict.1component4exp.1component2exp(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component4exp.1component2exp = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp4_line1 = params[4]
	A_Exp1_line2 = params[5]
	A_Exp2_line2 = params[6]
	B_Exp1_line1 = params[7]
	B_Exp2_line1 = params[8]
	B_Exp3_line1 = params[9]
	B_Exp4_line1 = params[10]
	B_Exp1_line2 = params[11]
	B_Exp2_line2 = params[12]
	scal_line1 = params[13]
	scal_line2 = params[14]
	xmid = params[15]
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1,A_Exp4_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1,B_Exp4_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2))
	c(cancer_TC50 = exp(xmid), noncancer_TC50 = exp(xmid), A_line1_avg=A_line1_avg, A_line2_avg=A_line2_avg, B_line1_avg=B_line1_avg, B_line2_avg=B_line2_avg)
}

##########
predict.1component3exp.1component2exp = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp1_line2 = params[4]
	A_Exp2_line2 = params[5]
	B_Exp1_line1 = params[6]
	B_Exp2_line1 = params[7]
	B_Exp3_line1 = params[8]
	B_Exp1_line2 = params[9]
	B_Exp2_line2 = params[10]
	scal_line1 = params[11]
	scal_line2 = params[12]
	xmid = params[13]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp3_line1 * is_line1 * is_Exp3) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp3_line1 * is_line1 * is_Exp3) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.1component3exp.1component2exp = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
{
	sigma = params[14]
	yhat = predict.1component3exp.1component2exp(params, x, is_line1, is_line2, is_Exp1, is_Exp2, is_Exp3)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component3exp.1component2exp = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp3_line1 = params[3]
	A_Exp1_line2 = params[4]
	A_Exp2_line2 = params[5]
	B_Exp1_line1 = params[6]
	B_Exp2_line1 = params[7]
	B_Exp3_line1 = params[8]
	B_Exp1_line2 = params[9]
	B_Exp2_line2 = params[10]
	scal_line1 = params[11]
	scal_line2 = params[12]
	xmid = params[13]
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1,A_Exp3_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1,B_Exp3_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2))
	c(cancer_TC50 = exp(xmid), noncancer_TC50 = exp(xmid), A_line1_avg=A_line1_avg, A_line2_avg=A_line2_avg, B_line1_avg=B_line1_avg, B_line2_avg=B_line2_avg)
}

##########
predict.1component2exp.1component2exp = function(params, x, is_line1, is_line2, is_Exp1, is_Exp2)
{
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp1_line2 = params[3]
	A_Exp2_line2 = params[4]
	B_Exp1_line1 = params[5]
	B_Exp2_line1 = params[6]
	B_Exp1_line2 = params[7]
	B_Exp2_line2 = params[8]
	scal_line1 = params[0]
	scal_line2 = params[10]
	xmid = params[11]
	finalA = (A_Exp1_line1 * is_line1 * is_Exp1) + (A_Exp2_line1 * is_line1 * is_Exp2) + (A_Exp1_line2 * is_line2 * is_Exp1) + (A_Exp2_line2 * is_line2 * is_Exp2)
	finalB = (B_Exp1_line1 * is_line1 * is_Exp1) + (B_Exp2_line1 * is_line1 * is_Exp2) + (B_Exp1_line2 * is_line2 * is_Exp1) + (B_Exp2_line2 * is_line2 * is_Exp2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.1component2exp.1component2exp = function(params, x, y, is_line1, is_line2, is_Exp1, is_Exp2)
{
	sigma = params[12]
	yhat = predict.1component2exp.1component2exp(params, x, is_line1, is_line2, is_Exp1, is_Exp2)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component2exp.1component2exp = function(params)
{
	params = as.numeric(params)
	A_Exp1_line1 = params[1]
	A_Exp2_line1 = params[2]
	A_Exp1_line2 = params[3]
	A_Exp2_line2 = params[4]
	B_Exp1_line1 = params[5]
	B_Exp2_line1 = params[6]
	B_Exp1_line2 = params[7]
	B_Exp2_line2 = params[8]
	scal_line1 = params[0]
	scal_line2 = params[10]
	xmid = params[11]
	A_line1_avg = mean(c(A_Exp1_line1,A_Exp2_line1))
	A_line2_avg = mean(c(A_Exp1_line2,A_Exp2_line2))
	B_line1_avg = mean(c(B_Exp1_line1,B_Exp2_line1))
	B_line2_avg = mean(c(B_Exp1_line2,B_Exp2_line2))
	c(cancer_TC50 = exp(xmid), noncancer_TC50 = exp(xmid), A_line1_avg=A_line1_avg, A_line2_avg=A_line2_avg, B_line1_avg=B_line1_avg, B_line2_avg=B_line2_avg)
}

####################################################################################################

##########
A_start = 1;            A_low = -Inf;            A_high = Inf;
B_start = 0.01;         B_low = 0;               B_high = 0.3;
scal_start = 0.1;       scal_low = 0.02;         scal_high = Inf;
xmid_start = 0;         xmid_low = -Inf;         xmid_high = Inf;
sens_start = 0.5;       sens_low = 0;            sens_high = 1;
xshift_start = 1;       xshift_low = 0;          xshift_high = Inf;
sigma_start = 0.1;      sigma_low = 1e-4;        sigma_high = Inf;
delta_start = 1;        delta_low = 0;           delta_high = Inf;

##########
start.1component5exp.2component3exp_delta2 = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp4_line1=A_start,  A_Exp5_line1=A_start,    A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   A_Exp3_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp4_line1=B_start,   B_Exp5_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   B_Exp3_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sens_Exp1_line2=sens_start,   sens_Exp2_line2=sens_start,   sens_Exp3_line2=sens_start,  delta2=delta_start,   sigma=sigma_start)

start.1component4exp.2component3exp_delta2 = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp4_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   A_Exp3_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp4_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   B_Exp3_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sens_Exp1_line2=sens_start,   sens_Exp2_line2=sens_start,   sens_Exp3_line2=sens_start,  delta2=delta_start,   sigma=sigma_start)

start.1component3exp.2component3exp_delta2 = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   A_Exp3_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   B_Exp3_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sens_Exp1_line2=sens_start,   sens_Exp2_line2=sens_start,   sens_Exp3_line2=sens_start,  delta2=delta_start,   sigma=sigma_start)

start.1component2exp.2component3exp_delta2 = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   A_Exp3_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   B_Exp3_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sens_Exp1_line2=sens_start,   sens_Exp2_line2=sens_start,   sens_Exp3_line2=sens_start,  delta2=delta_start,   sigma=sigma_start)

start.1component5exp.2component2exp_delta2 = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp4_line1=A_start,  A_Exp5_line1=A_start,    A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp4_line1=B_start,   B_Exp5_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sens_Exp1_line2=sens_start,   sens_Exp2_line2=sens_start,   delta2=delta_start,   sigma=sigma_start)

start.1component4exp.2component2exp_delta2 = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp4_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp4_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sens_Exp1_line2=sens_start,   sens_Exp2_line2=sens_start,   delta2=delta_start,   sigma=sigma_start)

start.1component3exp.2component2exp_delta2 = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sens_Exp1_line2=sens_start,   sens_Exp2_line2=sens_start,  delta2=delta_start,   sigma=sigma_start)

start.1component2exp.2component2exp_delta2 = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sens_Exp1_line2=sens_start,   sens_Exp2_line2=sens_start,   delta2=delta_start,   sigma=sigma_start)

##########
lower.1component5exp.2component3exp_delta2 = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp4_line1=A_low,   A_Exp5_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   A_Exp3_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp4_line1=B_low,   B_Exp5_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   B_Exp3_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sens_Exp1_line2=sens_low,   sens_Exp2_line2=sens_low,   sens_Exp3_line2=sens_low,   delta2=delta_low,   sigma=sigma_low)

lower.1component4exp.2component3exp_delta2 = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp4_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   A_Exp3_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp4_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   B_Exp3_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sens_Exp1_line2=sens_low,   sens_Exp2_line2=sens_low,   sens_Exp3_line2=sens_low,   delta2=delta_low,   sigma=sigma_low)

lower.1component3exp.2component3exp_delta2 = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   A_Exp3_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   B_Exp3_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sens_Exp1_line2=sens_low,   sens_Exp2_line2=sens_low,   sens_Exp3_line2=sens_low,   delta2=delta_low,   sigma=sigma_low)

lower.1component2exp.2component3exp_delta2 = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   A_Exp3_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   B_Exp3_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sens_Exp1_line2=sens_low,   sens_Exp2_line2=sens_low,   sens_Exp3_line2=sens_low,   delta2=delta_low,   sigma=sigma_low)

lower.1component5exp.2component2exp_delta2 = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp4_line1=A_low,   A_Exp5_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp4_line1=B_low,   B_Exp5_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sens_Exp1_line2=sens_low,   sens_Exp2_line2=sens_low,   delta2=delta_low,   sigma=sigma_low)

lower.1component4exp.2component2exp_delta2 = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp4_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp4_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sens_Exp1_line2=sens_low,   sens_Exp2_line2=sens_low,   delta2=delta_low,   sigma=sigma_low)

lower.1component3exp.2component2exp_delta2 = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sens_Exp1_line2=sens_low,   sens_Exp2_line2=sens_low,   delta2=delta_low,   sigma=sigma_low)

lower.1component2exp.2component2exp_delta2 = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sens_Exp1_line2=sens_low,   sens_Exp2_line2=sens_low,   delta2=delta_low,   sigma=sigma_low)

##########
upper.1component5exp.2component3exp_delta2 = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp4_line1=A_high,  A_Exp5_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  A_Exp3_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp4_line1=B_high,  B_Exp5_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  B_Exp3_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sens_Exp1_line2=sens_high,  sens_Exp2_line2=sens_high,  sens_Exp3_line2=sens_high,  delta2=delta_high,   sigma=sigma_high)

upper.1component4exp.2component3exp_delta2 = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp4_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  A_Exp3_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp4_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  B_Exp3_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sens_Exp1_line2=sens_high,  sens_Exp2_line2=sens_high,  sens_Exp3_line2=sens_high,  delta2=delta_high,   sigma=sigma_high)

upper.1component3exp.2component3exp_delta2 = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high, A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  A_Exp3_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  B_Exp3_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sens_Exp1_line2=sens_high,  sens_Exp2_line2=sens_high,  sens_Exp3_line2=sens_high,  delta2=delta_high,   sigma=sigma_high)

upper.1component2exp.2component3exp_delta2 = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  A_Exp3_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  B_Exp3_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sens_Exp1_line2=sens_high,  sens_Exp2_line2=sens_high,  sens_Exp3_line2=sens_high,  delta2=delta_high,   sigma=sigma_high)

upper.1component5exp.2component2exp_delta2 = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp4_line1=A_high,  A_Exp5_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp4_line1=B_high,  B_Exp5_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sens_Exp1_line2=sens_high,  sens_Exp2_line2=sens_high,  delta2=delta_high,   sigma=sigma_high)

upper.1component4exp.2component2exp_delta2 = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp4_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp4_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sens_Exp1_line2=sens_high,  sens_Exp2_line2=sens_high,  delta2=delta_high,   sigma=sigma_high)

upper.1component3exp.2component2exp_delta2 = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high, A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sens_Exp1_line2=sens_high,  sens_Exp2_line2=sens_high,  delta2=delta_high,   sigma=sigma_high)

upper.1component2exp.2component2exp_delta2 = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sens_Exp1_line2=sens_high,  sens_Exp2_line2=sens_high,  delta2=delta_high,   sigma=sigma_high)

##########
start.1component5exp.1component3exp = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp4_line1=A_start,   A_Exp5_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   A_Exp3_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp4_line1=B_start,   B_Exp5_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   B_Exp3_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)

start.1component4exp.1component3exp = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp4_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   A_Exp3_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp4_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   B_Exp3_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)

start.1component3exp.1component3exp = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   A_Exp3_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   B_Exp3_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)

start.1component2exp.1component3exp = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   A_Exp3_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   B_Exp3_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)

start.1component5exp.1component2exp = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp4_line1=A_start,   A_Exp5_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp4_line1=B_start,   B_Exp5_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)

start.1component4exp.1component2exp = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp4_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp4_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)

start.1component3exp.1component2exp = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp3_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp3_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)

start.1component2exp.1component2exp = c(A_Exp1_line1=A_start,   A_Exp2_line1=A_start,   A_Exp1_line2=A_start,   A_Exp2_line2=A_start,   B_Exp1_line1=B_start,   B_Exp2_line1=B_start,   B_Exp1_line2=B_start,   B_Exp2_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)


##########
lower.1component5exp.1component3exp = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp4_line1=A_low,   A_Exp5_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   A_Exp3_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp4_line1=B_low,   B_Exp5_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   B_Exp3_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sigma=sigma_low)

lower.1component4exp.1component3exp = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp4_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   A_Exp3_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp4_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   B_Exp3_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sigma=sigma_low)

lower.1component3exp.1component3exp = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   A_Exp3_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   B_Exp3_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sigma=sigma_low)

lower.1component2exp.1component3exp = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   A_Exp3_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   B_Exp3_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sigma=sigma_low)

lower.1component5exp.1component2exp = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp4_line1=A_low,   A_Exp5_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp4_line1=B_low,   B_Exp5_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sigma=sigma_low)

lower.1component4exp.1component2exp = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp4_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp4_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sigma=sigma_low)

lower.1component3exp.1component2exp = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp3_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp3_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sigma=sigma_low)

lower.1component2exp.1component2exp = c(A_Exp1_line1=A_low,   A_Exp2_line1=A_low,   A_Exp1_line2=A_low,   A_Exp2_line2=A_low,   B_Exp1_line1=B_low,   B_Exp2_line1=B_low,   B_Exp1_line2=B_low,   B_Exp2_line2=B_low,   scal_line1=scal_low,   scal_line2=scal_low,   xmid=xmid_low,   sigma=sigma_low)

##########
upper.1component5exp.1component3exp = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp4_line1=A_high,  A_Exp5_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  A_Exp3_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp4_line1=B_high,  B_Exp5_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  B_Exp3_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sigma=sigma_high)

upper.1component4exp.1component3exp = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp4_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  A_Exp3_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp4_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  B_Exp3_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sigma=sigma_high)

upper.1component3exp.1component3exp = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  A_Exp3_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  B_Exp3_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sigma=sigma_high)

upper.1component2exp.1component3exp = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  A_Exp3_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  B_Exp3_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sigma=sigma_high)

upper.1component5exp.1component2exp = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp4_line1=A_high,  A_Exp5_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp4_line1=B_high,  B_Exp5_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sigma=sigma_high)

upper.1component4exp.1component2exp = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp4_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp4_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sigma=sigma_high)

upper.1component3exp.1component2exp = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp3_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp3_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sigma=sigma_high)

upper.1component2exp.1component2exp = c(A_Exp1_line1=A_high,  A_Exp2_line1=A_high,  A_Exp1_line2=A_high,  A_Exp2_line2=A_high,  B_Exp1_line1=B_high,  B_Exp2_line1=B_high,  B_Exp1_line2=B_high,  B_Exp2_line2=B_high,  scal_line1=scal_high,  scal_line2=scal_high,   xmid=xmid_high,   sigma=sigma_high)

####################################################################################################
# Compare 2 cell-lines using the models provided on input.
# Those models were manually chosen by manually inspecting the output of step01_choose_model.R for each cell line.
####################################################################################################

cellline = assays[ ((assays$Line==args_cellline_1) | (assays$Line==args_cellline_2)), ]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_1), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_2), 1, 0)
cellline_1 = assays[ (assays$Line==args_cellline_1), ]
cellline_2 = assays[ (assays$Line==args_cellline_2), ]
exp_1 = unique(cellline_1$Experiment)
exp_2 = unique(cellline_2$Experiment)
# Actually, extra coding is required here when the 2 cell lines do not have the same number of experiments, 
# to ensure that the same experiment identifier appears with the same index in exp_1 and exp_2,
# and that the experiment identifiers existing in both cell lines are placed in the first positions of exp_1 and exp_2.

if ((args_model_1 == '1component') & (args_model_2 == '2component')) {

	if ((length(exp_1) == 5) & (length(exp_2) == 3)) {
		cellline$is_Exp1 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[1])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[1])), 1, 0)
		cellline$is_Exp2 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[2])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[2])), 1, 0)
		cellline$is_Exp3 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[3])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[3])), 1, 0)
		cellline$is_Exp4 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[4])), 1, 0)
		cellline$is_Exp5 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[5])), 1, 0)

		cellline1_cellline2_fit.1component5exp.2component3exp_delta2 = optimx(start.1component5exp.2component3exp_delta2, function(par) loglik.1component5exp.2component3exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component5exp.2component3exp_delta2, upper = upper.1component5exp.2component3exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))

		cellline1_cellline2_fit.1component5exp.1component3exp = optimx(start.1component5exp.1component3exp, function(par) loglik.1component5exp.1component3exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component5exp.1component3exp, upper = upper.1component5exp.1component3exp, method = c("L-BFGS-B"), control = list(fnscale = -1))

		degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
		cellline1_cellline2.2component_delta2_vs_1component.pvalue = pchisq(2*(cellline1_cellline2_fit.1component5exp.2component3exp_delta2$value - cellline1_cellline2_fit.1component5exp.1component3exp$value), degrees_of_freedom, lower.tail = FALSE)

		cellline1_cellline2_fit.1component5exp.2component3exp_delta2
		cellline1_cellline2_fit.1component5exp.1component3exp
		convert_values_back_to_interpretable.1component5exp.2component3exp_delta2(cellline1_cellline2_fit.1component5exp.2component3exp_delta2)
		convert_values_back_to_interpretable.1component5exp.1component3exp(cellline1_cellline2_fit.1component5exp.1component3exp)
		cellline1_cellline2.2component_delta2_vs_1component.pvalue
	}

	if ((length(exp_1) == 4) & (length(exp_2) == 3)) {
		cellline$is_Exp1 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[1])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[1])), 1, 0)
		cellline$is_Exp2 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[2])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[2])), 1, 0)
		cellline$is_Exp3 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[3])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[3])), 1, 0)
		cellline$is_Exp4 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[4])), 1, 0)

		cellline1_cellline2_fit.1component4exp.2component3exp_delta2 = optimx(start.1component4exp.2component3exp_delta2, function(par) loglik.1component4exp.2component3exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component4exp.2component3exp_delta2, upper = upper.1component4exp.2component3exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))

		cellline1_cellline2_fit.1component4exp.1component3exp = optimx(start.1component4exp.1component3exp, function(par) loglik.1component4exp.1component3exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component4exp.1component3exp, upper = upper.1component4exp.1component3exp, method = c("L-BFGS-B"), control = list(fnscale = -1))

		degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
		cellline1_cellline2.2component_delta2_vs_1component.pvalue = pchisq(2*(cellline1_cellline2_fit.1component4exp.2component3exp_delta2$value - cellline1_cellline2_fit.1component4exp.1component3exp$value), degrees_of_freedom, lower.tail = FALSE)

		cellline1_cellline2_fit.1component4exp.2component3exp_delta2
		cellline1_cellline2_fit.1component4exp.1component3exp
		convert_values_back_to_interpretable.1component4exp.2component3exp_delta2(cellline1_cellline2_fit.1component4exp.2component3exp_delta2)
		convert_values_back_to_interpretable.1component4exp.1component3exp(cellline1_cellline2_fit.1component4exp.1component3exp)
		cellline1_cellline2.2component_delta2_vs_1component.pvalue
	}

	if ((length(exp_1) == 3) & (length(exp_2) == 3)) {
		cellline$is_Exp1 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[1])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[1])), 1, 0)
		cellline$is_Exp2 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[2])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[2])), 1, 0)
		cellline$is_Exp3 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[3])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[3])), 1, 0)

		cellline1_cellline2_fit.1component3exp.2component3exp_delta2 = optimx(start.1component3exp.2component3exp_delta2, function(par) loglik.1component3exp.2component3exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component3exp.2component3exp_delta2, upper = upper.1component3exp.2component3exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))

		cellline1_cellline2_fit.1component3exp.1component3exp = optimx(start.1component3exp.1component3exp, function(par) loglik.1component3exp.1component3exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component3exp.1component3exp, upper = upper.1component3exp.1component3exp, method = c("L-BFGS-B"), control = list(fnscale = -1))

		degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
		cellline1_cellline2.2component_delta2_vs_1component.pvalue = pchisq(2*(cellline1_cellline2_fit.1component3exp.2component3exp_delta2$value - cellline1_cellline2_fit.1component3exp.1component3exp$value), degrees_of_freedom, lower.tail = FALSE)

		cellline1_cellline2_fit.1component3exp.2component3exp_delta2
		cellline1_cellline2_fit.1component3exp.1component3exp
		convert_values_back_to_interpretable.1component3exp.2component3exp_delta2(cellline1_cellline2_fit.1component3exp.2component3exp_delta2)
		convert_values_back_to_interpretable.1component3exp.1component3exp(cellline1_cellline2_fit.1component3exp.1component3exp)
		cellline1_cellline2.2component_delta2_vs_1component.pvalue
	}

	if ((length(exp_1) == 2) & (length(exp_2) == 3)) {
		cellline$is_Exp1 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[1])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[1])), 1, 0)
		cellline$is_Exp2 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[2])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[2])), 1, 0)
		cellline$is_Exp3 <- ifelse( ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[3])), 1, 0)

		cellline1_cellline2_fit.1component2exp.2component3exp_delta2 = optimx(start.1component2exp.2component3exp_delta2, function(par) loglik.1component2exp.2component3exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component3exp.2component3exp_delta2, upper = upper.1component3exp.2component3exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))

		cellline1_cellline2_fit.1component2exp.1component3exp = optimx(start.1component2exp.1component3exp, function(par) loglik.1component2exp.1component3exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component2exp.1component3exp, upper = upper.1component2exp.1component3exp, method = c("L-BFGS-B"), control = list(fnscale = -1))

		degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
		cellline1_cellline2.2component_delta2_vs_1component.pvalue = pchisq(2*(cellline1_cellline2_fit.1component2exp.2component3exp_delta2$value - cellline1_cellline2_fit.1component2exp.1component3exp$value), degrees_of_freedom, lower.tail = FALSE)

		cellline1_cellline2_fit.1component2exp.2component3exp_delta2
		cellline1_cellline2_fit.1component2exp.1component3exp
		convert_values_back_to_interpretable.1component2exp.2component3exp_delta2(cellline1_cellline2_fit.1component2exp.2component3exp_delta2)
		convert_values_back_to_interpretable.1component2exp.1component3exp(cellline1_cellline2_fit.1component2exp.1component3exp)
		cellline1_cellline2.2component_delta2_vs_1component.pvalue
	}

	if ((length(exp_1) == 5) & (length(exp_2) == 2)) {
		cellline$is_Exp1 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[1])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[1])), 1, 0)
		cellline$is_Exp2 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[2])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[2])), 1, 0)
		cellline$is_Exp3 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[3])), 1, 0)
		cellline$is_Exp4 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[4])), 1, 0)
		cellline$is_Exp5 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[5])), 1, 0)

		cellline1_cellline2_fit.1component5exp.2component2exp_delta2 = optimx(start.1component5exp.2component2exp_delta2, function(par) loglik.1component5exp.2component2exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component5exp.2component2exp_delta2, upper = upper.1component5exp.2component2exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))

		cellline1_cellline2_fit.1component5exp.1component2exp = optimx(start.1component5exp.1component2exp, function(par) loglik.1component5exp.1component2exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component5exp.1component2exp, upper = upper.1component5exp.1component2exp, method = c("L-BFGS-B"), control = list(fnscale = -1))

		degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
		cellline1_cellline2.2component_delta2_vs_1component.pvalue = pchisq(2*(cellline1_cellline2_fit.1component5exp.2component2exp_delta2$value - cellline1_cellline2_fit.1component5exp.1component2exp$value), degrees_of_freedom, lower.tail = FALSE)

		cellline1_cellline2_fit.1component5exp.2component2exp_delta2
		cellline1_cellline2_fit.1component5exp.1component2exp
		convert_values_back_to_interpretable.1component5exp.2component2exp_delta2(cellline1_cellline2_fit.1component5exp.2component2exp_delta2)
		convert_values_back_to_interpretable.1component5exp.1component2exp(cellline1_cellline2_fit.1component5exp.1component2exp)
		cellline1_cellline2.2component_delta2_vs_1component.pvalue
	}

	if ((length(exp_1) == 4) & (length(exp_2) == 2)) {
		cellline$is_Exp1 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[1])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[1])), 1, 0)
		cellline$is_Exp2 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[2])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[2])), 1, 0)
		cellline$is_Exp3 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[3])), 1, 0)
		cellline$is_Exp4 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[4])), 1, 0)

		cellline1_cellline2_fit.1component4exp.2component2exp_delta2 = optimx(start.1component4exp.2component2exp_delta2, function(par) loglik.1component4exp.2component2exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component4exp.2component2exp_delta2, upper = upper.1component4exp.2component2exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))

		cellline1_cellline2_fit.1component4exp.1component2exp = optimx(start.1component4exp.1component2exp, function(par) loglik.1component4exp.1component2exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component4exp.1component2exp, upper = upper.1component4exp.1component2exp, method = c("L-BFGS-B"), control = list(fnscale = -1))

		degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
		cellline1_cellline2.2component_delta2_vs_1component.pvalue = pchisq(2*(cellline1_cellline2_fit.1component4exp.2component2exp_delta2$value - cellline1_cellline2_fit.1component4exp.1component2exp$value), degrees_of_freedom, lower.tail = FALSE)

		cellline1_cellline2_fit.1component4exp.2component2exp_delta2
		cellline1_cellline2_fit.1component4exp.1component2exp
		convert_values_back_to_interpretable.1component4exp.2component2exp_delta2(cellline1_cellline2_fit.1component4exp.2component2exp_delta2)
		convert_values_back_to_interpretable.1component4exp.1component2exp(cellline1_cellline2_fit.1component4exp.1component2exp)
		cellline1_cellline2.2component_delta2_vs_1component.pvalue
	}

	if ((length(exp_1) == 3) & (length(exp_2) == 2)) {
		cellline$is_Exp1 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[1])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[1])), 1, 0)
		cellline$is_Exp2 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[2])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[2])), 1, 0)
		cellline$is_Exp3 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[3])), 1, 0)

		cellline1_cellline2_fit.1component3exp.2component2exp_delta2 = optimx(start.1component3exp.2component2exp_delta2, function(par) loglik.1component3exp.2component2exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component3exp.2component2exp_delta2, upper = upper.1component3exp.2component2exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))

		cellline1_cellline2_fit.1component3exp.1component2exp = optimx(start.1component3exp.1component2exp, function(par) loglik.1component3exp.1component2exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component3exp.1component2exp, upper = upper.1component3exp.1component2exp, method = c("L-BFGS-B"), control = list(fnscale = -1))

		degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
		cellline1_cellline2.2component_delta2_vs_1component.pvalue = pchisq(2*(cellline1_cellline2_fit.1component3exp.2component2exp_delta2$value - cellline1_cellline2_fit.1component3exp.1component2exp$value), degrees_of_freedom, lower.tail = FALSE)

		cellline1_cellline2_fit.1component3exp.2component2exp_delta2
		cellline1_cellline2_fit.1component3exp.1component2exp
		convert_values_back_to_interpretable.1component3exp.2component2exp_delta2(cellline1_cellline2_fit.1component3exp.2component2exp_delta2)
		convert_values_back_to_interpretable.1component3exp.1component2exp(cellline1_cellline2_fit.1component3exp.1component2exp)
		cellline1_cellline2.2component_delta2_vs_1component.pvalue
	}

	if ((length(exp_1) == 2) & (length(exp_2) == 2)) {
		cellline$is_Exp1 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[1])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[1])), 1, 0)
		cellline$is_Exp2 <- ifelse( ((cellline$Line==args_cellline_1)&(cellline$Experiment==exp[2])) | ((cellline$Line==args_cellline_2)&(cellline$Experiment==exp[2])), 1, 0)

		cellline1_cellline2_fit.1component2exp.2component2exp_delta2 = optimx(start.1component2exp.2component2exp_delta2, function(par) loglik.1component2exp.2component2exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component2exp.2component2exp_delta2, upper = upper.1component2exp.2component2exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))

		cellline1_cellline2_fit.1component2exp.1component2exp = optimx(start.1component2exp.1component2exp, function(par) loglik.1component2exp.1component2exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component2exp.1component2exp, upper = upper.1component2exp.1component2exp, method = c("L-BFGS-B"), control = list(fnscale = -1))

		degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
		cellline1_cellline2.2component_delta2_vs_1component.pvalue = pchisq(2*(cellline1_cellline2_fit.1component2exp.2component2exp_delta2$value - cellline1_cellline2_fit.1component2exp.1component2exp$value), degrees_of_freedom, lower.tail = FALSE)

		cellline1_cellline2_fit.1component2exp.2component2exp_delta2
		cellline1_cellline2_fit.1component2exp.1component2exp
		convert_values_back_to_interpretable.1component2exp.2component2exp_delta2(cellline1_cellline2_fit.1component2exp.2component2exp_delta2)
		convert_values_back_to_interpretable.1component42exp.1component2exp(cellline1_cellline2_fit.1component2exp.1component2exp)
		cellline1_cellline2.2component_delta2_vs_1component.pvalue
	}
}

if ((args_model_1 == '1component') & (args_model_2 == '1component')) {

	if ((length(exp_1) == 5) & (length(exp_2) == 3)) {
	}

	if ((length(exp_1) == 4) & (length(exp_2) == 3)) {
	}

	if ((length(exp_1) == 3) & (length(exp_2) == 3)) {
	}

	if ((length(exp_1) == 2) & (length(exp_2) == 3)) {
	}

	if ((length(exp_1) == 5) & (length(exp_2) == 2)) {
	}

	if ((length(exp_1) == 4) & (length(exp_2) == 2)) {
	}

	if ((length(exp_1) == 3) & (length(exp_2) == 2)) {
	}

	if ((length(exp_1) == 2) & (length(exp_2) == 2)) {
	}
}

# Example outputs

####################################################################################################
##### Compare MM05 mesothelioma and HOFM01 non-cancer primary fibroblasts
# 
# cellline = assays[ ((assays$Line=="MM05") | (assays$Line=="HOFM01")), ]
# cellline$is_line1 <- ifelse( (cellline$Line=="MM05"), 1, 0)
# cellline$is_line2 <- ifelse( (cellline$Line=="HOFM01"), 1, 0)
# cellline$is_Exp1 <- ifelse( ((cellline$Line=="MM05")&(cellline$Experiment=="20151227")) | ((cellline$Line=="HOFM01")&(cellline$Experiment=="20151227")), 1, 0)
# cellline$is_Exp2 <- ifelse( ((cellline$Line=="MM05")&(cellline$Experiment=="20160101")) | ((cellline$Line=="HOFM01")&(cellline$Experiment=="20160101")), 1, 0)
# cellline$is_Exp3 <- ifelse( ((cellline$Line=="MM05")&(cellline$Experiment=="20160108")) | ((cellline$Line=="HOFM01")&(cellline$Experiment=="20160108")), 1, 0)
# cellline$is_Exp4 <- ifelse( ((cellline$Line=="MM05")&(cellline$Experiment=="20160206")), 1, 0)
# 
# MM05_HOFM01_fit.1component4exp.2component3exp_delta2 = optimx(start.1component4exp.2component3exp_delta2, function(par) loglik.1component4exp.2component3exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component4exp.2component3exp_delta2, upper = upper.1component4exp.2component3exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))
# 
# MM05_HOFM01_fit.1component4exp.1component3exp = optimx(start.1component4exp.1component3exp, function(par) loglik.1component4exp.1component3exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4), lower = lower.1component4exp.1component3exp, upper = upper.1component4exp.1component3exp, method = c("L-BFGS-B"), control = list(fnscale = -1))
# 
# degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
# MM05_HOFM01.2component_delta2_vs_1component.pvalue = pchisq(2*(MM05_HOFM01_fit.1component4exp.2component3exp_delta2$value - MM05_HOFM01_fit.1component4exp.1component3exp$value), degrees_of_freedom, lower.tail = FALSE)
# 
# MM05_HOFM01_fit.1component4exp.2component3exp_delta2
# MM05_HOFM01_fit.1component4exp.1component3exp
# convert_values_back_to_interpretable.1component4exp.2component3exp_delta2(MM05_HOFM01_fit.1component4exp.2component3exp_delta2)
# convert_values_back_to_interpretable.1component4exp.1component3exp(MM05_HOFM01_fit.1component4exp.1component3exp)
# MM05_HOFM01.2component_delta2_vs_1component.pvalue
# 
# 	> MM05_HOFM01_fit.1component4exp.2component3exp_delta2
# 		 A_Exp1_line1 A_Exp2_line1 A_Exp3_line1 A_Exp4_line1 A_Exp1_line2
# 	L-BFGS-B    0.9889134    0.9831092    0.9675529    0.8781448     1.098271
# 		 A_Exp2_line2 A_Exp3_line2 B_Exp1_line1 B_Exp2_line1 B_Exp3_line1
# 	L-BFGS-B     1.019177     1.232038   0.06285519   0.04661661   0.09947674
# 		 B_Exp4_line1 B_Exp1_line2 B_Exp2_line2 B_Exp3_line2 scal_line1
# 	L-BFGS-B  0.003711902    0.0980046   0.05922539    0.2868314       0.02
# 		 scal_line2       xmid sens_Exp1_line2 sens_Exp2_line2 sens_Exp3_line2
# 	L-BFGS-B 0.07375564 -0.4855751       0.7374322       0.5541991       0.3698963
# 		    delta2     sigma    value fevals gevals niter convcode  kkt1  kkt2
# 	L-BFGS-B 0.7478481 0.1215929 125.9423    134    134    NA        0 FALSE FALSE
# 		 xtimes
# 	L-BFGS-B   1.52
# 	> MM05_HOFM01_fit.1component4exp.1component3exp
# 		 A_Exp1_line1 A_Exp2_line1 A_Exp3_line1 A_Exp4_line1 A_Exp1_line2
# 	L-BFGS-B    0.9818116    0.9725007    0.9622348    0.8618916      1.27281
# 		 A_Exp2_line2 A_Exp3_line2 B_Exp1_line1 B_Exp2_line1 B_Exp3_line1
# 	L-BFGS-B     1.230804     1.395433   0.03706417   0.02201568   0.06714598
# 		 B_Exp4_line1 B_Exp1_line2 B_Exp2_line2 B_Exp3_line2 scal_line1
# 	L-BFGS-B            0            0   0.04478798          0.3 0.07898562
# 		 scal_line2       xmid     sigma    value fevals gevals niter convcode
# 	L-BFGS-B  0.4369826 -0.3952124 0.1433815 95.77233     51     51    NA        0
# 		  kkt1  kkt2 xtimes
# 	L-BFGS-B FALSE FALSE  0.396
# 	> convert_values_back_to_interpretable.1component4exp.2component3exp_delta2(MM05_HOFM01_fit.1component4exp.2component3exp_delta2)
# 		              cancer_TC50 noncancer_sensitive_fraction_TC50 
# 		               0.61534323                        0.61534323 
# 	noncancer_resistant_fraction_TC50        sensitive_fraction_average 
# 		               1.29988134                        0.55384252 
# 	       resistant_fraction_average                       A_line1_avg 
# 		               0.44615748                        0.95443007 
# 		              A_line2_avg                       B_line1_avg 
# 		               1.11649535                        0.05316511 
# 		              B_line2_avg 
# 		               0.14802047 
# 	> convert_values_back_to_interpretable.1component4exp.1component3exp(MM05_HOFM01_fit.1component4exp.1component3exp)
# 	   cancer_TC50 noncancer_TC50    A_line1_avg    A_line2_avg    B_line1_avg 
# 	    0.67353695     0.67353695     0.94460967     1.29968225     0.03155646 
# 	   B_line2_avg 
# 	    0.11492933 
# 	> MM05_HOFM01.2component_delta2_vs_1component.pvalue
# 	[1] 7.894963e-14
# 
# 	# Non-cancer HOFM01 2component_delta2 model BAMLET-resistant TC50 = 1.3 mg/ml and mesothelioma cancer TC50 is 0.6 mg/ml.
# 	# On average 45% (avg. of 26%, 45%, and 63%) of HOFM01 cells are BAMLET-resistant.
# 	# This BAMLET-resistant TC50 is statistically significantly higher than the mesothelioma MM05 TC50, with p-value of 8e-14.
# 	# The difference in cancer and non-cancer BAMLET-resistant TC50 is considered to be "clinically significant"
# 	# because the BAMLET-resistant TC50 is more than twice the cancer TC50.

####################################################################################################
##### Compare H226 mesothelioma and HOFM02 non-cancer primary fibroblasts
# 
# cellline = assays[ ((assays$Line=="H226") | (assays$Line=="HOFM02")), ]
# cellline$is_line1 <- ifelse( (cellline$Line=="H226"), 1, 0)
# cellline$is_line2 <- ifelse( (cellline$Line=="HOFM02"), 1, 0)
# cellline$is_Exp1 <- ifelse( ((cellline$Line=="H226")&(cellline$Experiment=="20151218")) | ((cellline$Line=="HOFM02")&(cellline$Experiment=="20151227")), 1, 0)
# cellline$is_Exp2 <- ifelse( ((cellline$Line=="H226")&(cellline$Experiment=="20151227")) | ((cellline$Line=="HOFM02")&(cellline$Experiment=="20160101")), 1, 0)
# cellline$is_Exp3 <- ifelse( ((cellline$Line=="H226")&(cellline$Experiment=="20160101")) | ((cellline$Line=="HOFM02")&(cellline$Experiment=="20160108")), 1, 0)
# cellline$is_Exp4 <- ifelse( ((cellline$Line=="H226")&(cellline$Experiment=="20160108")), 1, 0)
# cellline$is_Exp5 <- ifelse( ((cellline$Line=="H226")&(cellline$Experiment=="20160206")), 1, 0)
# 
# H226_HOFM02_fit.1component5exp.2component3exp_delta2 = optimx(start.1component5exp.2component3exp_delta2, function(par) loglik.1component5exp.2component3exp_delta2(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4, is_Exp5=cellline$is_Exp5), lower = lower.1component5exp.2component3exp_delta2, upper = upper.1component5exp.2component3exp_delta2, method = c("L-BFGS-B"), control = list(fnscale = -1))
# 
# H226_HOFM02_fit.1component5exp.1component3exp = optimx(start.1component5exp.1component3exp, function(par) loglik.1component5exp.1component3exp(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3, is_Exp4=cellline$is_Exp4, is_Exp5=cellline$is_Exp5), lower = lower.1component5exp.1component3exp, upper = upper.1component5exp.1component3exp, method = c("L-BFGS-B"), control = list(fnscale = -1))
# 
# degrees_of_freedom=2 # 2component_delta2 has 2 additional variables (delta2 and sens) that 1component doesn't have
# H226_HOFM02.2component_delta2_vs_1component.pvalue = pchisq(2*(H226_HOFM02_fit.1component5exp.2component3exp_delta2$value - H226_HOFM02_fit.1component5exp.1component3exp$value), degrees_of_freedom, lower.tail = FALSE)
# 
# H226_HOFM02_fit.1component5exp.2component3exp_delta2
# H226_HOFM02_fit.1component5exp.1component3exp
# convert_values_back_to_interpretable.1component5exp.2component3exp_delta2(H226_HOFM02_fit.1component5exp.2component3exp_delta2)
# convert_values_back_to_interpretable.1component5exp.1component3exp(H226_HOFM02_fit.1component5exp.1component3exp)
# H226_HOFM02.2component_delta2_vs_1component.pvalue
# 
# 	> H226_HOFM02_fit.1component5exp.2component3exp_delta2
# 		 A_Exp1_line1 A_Exp2_line1 A_Exp3_line1 A_Exp4_line1 A_Exp5_line1
# 	L-BFGS-B     1.080497     1.071855    0.9400479    0.9890644    0.8069672
# 		 A_Exp1_line2 A_Exp2_line2 A_Exp3_line2 B_Exp1_line1 B_Exp2_line1
# 	L-BFGS-B     1.372845     1.052932     1.133171   0.08928352   0.07627034
# 		 B_Exp3_line1 B_Exp4_line1 B_Exp5_line1 B_Exp1_line2 B_Exp2_line2
# 	L-BFGS-B    0.0732545   0.05800311            0     0.112372   0.07004363
# 		 B_Exp3_line2 scal_line1 scal_line2       xmid sens_Exp1_line2
# 	L-BFGS-B    0.2101461       0.02 0.07727609 -0.4873377       0.7309922
# 		 sens_Exp2_line2 sens_Exp3_line2    delta2     sigma    value fevals
# 	L-BFGS-B       0.5040676       0.4783274 0.7570214 0.1163585 153.7624    112
# 		 gevals niter convcode  kkt1  kkt2 xtimes
# 	L-BFGS-B    112    NA        0 FALSE FALSE  1.524
# 	> H226_HOFM02_fit.1component5exp.1component3exp
# 		 A_Exp1_line1 A_Exp2_line1 A_Exp3_line1 A_Exp4_line1 A_Exp5_line1
# 	L-BFGS-B      1.04918     1.036774    0.9096231    0.9545356    0.7311433
# 		 A_Exp1_line2 A_Exp2_line2 A_Exp3_line2 B_Exp1_line1 B_Exp2_line1
# 	L-BFGS-B      1.36731     1.083081     1.216204            0            0
# 		 B_Exp3_line1 B_Exp4_line1 B_Exp5_line1 B_Exp1_line2 B_Exp2_line2
# 	L-BFGS-B            0            0            0    0.1143339    0.1188122
# 		 B_Exp3_line2 scal_line1 scal_line2       xmid     sigma  value fevals
# 	L-BFGS-B    0.1920026       0.02  0.3118619 -0.5224716 0.1643517 23.288     13
# 		 gevals niter convcode  kkt1  kkt2 xtimes
# 	L-BFGS-B     13    NA        0 FALSE FALSE  0.124
# 	> convert_values_back_to_interpretable.1component5exp.2component3exp_delta2(H226_HOFM02_fit.1component5exp.2component3exp_delta2)
# 		              cancer_TC50 noncancer_sensitive_fraction_TC50 
# 		               0.61425959                        0.61425959 
# 	noncancer_resistant_fraction_TC50        sensitive_fraction_average 
# 		               1.30955024                        0.57112909 
# 	       resistant_fraction_average                       A_line1_avg 
# 		               0.42887091                        0.97768636 
# 		              A_line2_avg                       B_line1_avg 
# 		               1.18631579                        0.05936229 
# 		              B_line2_avg 
# 		               0.13085392 
# 	> convert_values_back_to_interpretable.1component5exp.1component3exp(H226_HOFM02_fit.1component5exp.1component3exp)
# 	   cancer_TC50 noncancer_TC50    A_line1_avg    A_line2_avg    B_line1_avg 
# 	     0.5930530      0.5930530      0.9362513      1.2221984      0.0000000 
# 	   B_line2_avg 
# 	     0.1417162 
# 	> H226_HOFM02.2component_delta2_vs_1component.pvalue
# 	[1] 2.166151e-57
# 
# 	# Non-cancer HOFM02 2component_delta2 model BAMLET-resistant TC50 = 1.3 mg/ml and mesothelioma cancer TC50 is 0.6 mg/ml.
# 	# On average 43% (avg. of 27%, 50%, and 52%) of HOFM02 cells are BAMLET-resistant.
# 	# This BAMLET-resistant TC50 is statistically significantly higher than the mesothelioma H226 TC50, with p-value of 2e-57.
# 	# The difference in cancer and non-cancer BAMLET-resistant TC50 is considered to be "clinically significant"
# 	# because the BAMLET-resistant TC50 is more than twice the cancer TC50.


