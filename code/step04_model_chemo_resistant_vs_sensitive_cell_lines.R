# Rscript step04_model_chemo_resistant_vs_sensitive_cell_lines.R <input_file_1> <input_file_2> <cell_line> > <output_file>
# Rscript step04_model_chemo_resistant_vs_sensitive_cell_lines.R rat_chemo_cell_assay_data_1.txt rat_chemo_cell_assay_data_2.txt rat_bamlet_cell_assay_data.txt IL45CisR_highDose IL45_high_passage IL45CisR_lowDose IL45_low_passage > step04_output.txt

# Hypothesis for step04: That chemo-resistant cells are just as sensitive to BAMLET as chemo-sensitive cells.
# First establish that chemo-resistant cells are in fact less sensitive to chemo (cisplatin) than chemo-sensitive cells.
# Also establish that there is no difference in sensitivity (to chemo or BAMLET) between low-passage and high-passage cells,
# so counter the possibility that chemo-resistant cells are actually resistant because they are usually high-passage cells.
# Additionally, the chemo-resistant cells had their chemo-resistance maintained by regular chemotherapy treatment that was either high-dose or low-lose.
# These cell lines have a standard uniphasic response, not a biphasic response. No need to model 2component models for different TC50 in same cell line.

library(reshape2)
library(ggplot2)
library(optimx)
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'rat_chemo_cell_assay_data_1.txt', 'rat_chemo_cell_assay_data_2.txt', 'rat_bamlet_cell_assay_data.txt', 'IL45CisR_highDose', 'IL45_high_passage', 'IL45CisR_lowDose', 'IL45_low_passage' ) # for testing

args_infile_chemo1_1 = args[1] # rat_chemo_cell_assay_data_1.txt # Can be either the higher or lower concentrations,
args_infile_chemo1_2 = args[2] # rat_chemo_cell_assay_data_2.txt # because will all go into one data frame.
args_infile_treatment = args[3] # rat_bamlet_cell_assay_data.txt # 
args_cellline_chemoResistant_highDose = args[4] # IL45CisR_highDose
args_cellline_chemoResistant_lowDose = args[5] # IL45_high_passage
args_cellline_chemoSensitive_highPassage = args[6] # IL45CisR_lowDose
args_cellline_chemoSensitive_lowPassage = args[7] # IL45_low_passage

# Read in 2 data files. They have different drug concentration columns
rat = read.table(args_infile_chemo1_1, sep = "\t", header = TRUE)
rat = melt(rat, id.vars = c("Drug", "Line", "Experiment"), variable.name = "Concentration", value.name = "Viability")
rat$Concentration = as.numeric(gsub("^X", "", rat$Concentration))
rat = rat[!is.na(rat$Viability),]
rat_chemo = rat
rat = read.table(args_infile_chemo1_2, sep = "\t", header = TRUE)
rat = melt(rat, id.vars = c("Drug", "Line", "Experiment"), variable.name = "Concentration", value.name = "Viability")
rat$Concentration = as.numeric(gsub("^X", "", rat$Concentration))
rat = rat[!is.na(rat$Viability),]
rat_chemo = rbind( rat_chemo, rat )
rat_chemo$Organism = 'Rat'
rat_chemo = rat_chemo[ (rat_chemo$Concentration != 0), ]

assays = rat_chemo
assays$Drug = as.factor(assays$Drug)
assays$Organism = as.factor(assays$Organism)
assays$Line = as.factor(assays$Line)
assays$Experiment = as.factor(as.character(assays$Experiment))
assays = assays[!is.na(assays$Viability),]
assays$log10Conc = log10(assays$Concentration)
rat_chemo = assays

rat_Cisplatin = rat_chemo[ (rat_chemo$Drug=="Cisplatin"), ]

rat = read.table(args_infile_treatment, sep = "\t", header = TRUE)

rat = melt(rat, id.vars = c("Drug", "Line", "Experiment"), variable.name = "Concentration", value.name = "Viability")
rat$Concentration = as.numeric(gsub("^X", "", rat$Concentration))
rat_bamlet = rat
rat_bamlet$Drug = 'BAMLET'
rat_bamlet$Organism = 'Rat'
rat_bamlet = rat_bamlet[ (rat_bamlet$Concentration != 0), ]

assays = rat_bamlet
assays$Drug = as.factor(assays$Drug)
assays$Organism = as.factor(assays$Organism)
assays$Line = as.factor(assays$Line)
assays$Experiment = as.factor(as.character(assays$Experiment))
assays = assays[!is.na(assays$Viability),]
assays$log10Conc = log10(assays$Concentration)
rat_bamlet = assays

assays = rat_bamlet
ggplot(assays, aes(x = log10(Concentration), y = Viability, col = Experiment)) + geom_point() + facet_wrap(~ Line)


####################################################################################################
####################################################################################################

##########
A_start = 1;			A_low = -Inf;			A_high = Inf;
B_start = 0;			B_low = -Inf;			B_high = Inf;
scal_start = 0.1;		scal_low = 0.02;		scal_high = Inf;
xmid_start = 0;			xmid_low = -Inf;		xmid_high = Inf;
sens_start = 0.5;		sens_low = 0;			sens_high = 1;
xshift_start = 1;		xshift_low = 0;			xshift_high = Inf;
sigma_start = 0.1;		sigma_low = 1e-4;		sigma_high = Inf;
ResistanceDelta_start = 1;	ResistanceDelta_low = 0;	ResistanceDelta_high = Inf;
ExperimentDelta_start = 1;	ExperimentDelta_low = 0;	ExperimentDelta_high = Inf;
NotReachedDelta_start = 1;	NotReachedDelta_low = 0;	NotReachedDelta_high = Inf;

##########
# Sigmoidal model
# A = left asymptote
# B = right asymptote
# xmid = inflection point (TC50)
# scal = factor affecting slope at inflection point
predict.1component = function(params, x)
{
	A = params[1]
	B = params[2]
	scal = params[3]
	xmid = params[4]
	prediction_result = A + (B - A) * (1/(1 + exp((xmid - x) / scal)))
	prediction_result
}
loglik.1component = function(params, x, y)
{
	sigma = params[5]
	yhat = predict.1component(params, x)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component = function(params)
{
	names(params) <- NULL 
	xmid = params[4]
	c(TC50 = exp(xmid))
}

start.1component  = c(A=A_start,   B=B_start,   scal=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.1component  = c(A=A_low,     B=B_low,     scal=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.1component  = c(A=A_high,    B=B_high,    scal=scal_high,    xmid=xmid_high,    sigma=sigma_high)

##########
# A sigmoidal model where line1 and line2 experiments have the same TC50.
predict.model2lines = function(params, x, is_line1, is_line2)
{
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	finalA = (A_line1 * is_line1) + (A_line2 * is_line2)
	finalB = (B_line1 * is_line1) + (B_line2 * is_line2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp((xmid - x) / finalScal)) ) 
	prediction_result
}
loglik.model2lines = function(params, x, y, is_line1, is_line2)
{
	sigma = params[8]
	yhat = predict.model2lines(params, x, is_line1, is_line2)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.model2lines = function(params)
{
	params = as.numeric(params)
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	c(TC50 = exp(xmid))
}

start.model2lines  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)

##########
# A sigmoidal model where line1 experiments have a resistance delta added to the line2 TC50 value.
predict.ResistanceDelta = function(params, x, is_line1, is_line2)
{
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	ResistanceDelta_line1 = params[8]
	finalA = (A_line1 * is_line1) + (A_line2 * is_line2)
	finalB = (B_line1 * is_line1) + (B_line2 * is_line2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalResistanceDelta = (ResistanceDelta_line1 * is_line1)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp(( (xmid+finalResistanceDelta) - x) / finalScal)) ) 
	prediction_result
}
loglik.ResistanceDelta = function(params, x, y, is_line1, is_line2)
{
	sigma = params[9]
	yhat = predict.ResistanceDelta(params, x, is_line1, is_line2)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.ResistanceDelta = function(params)
{
	params = as.numeric(params)
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	ResistanceDelta_line1 = params[8]
	c(higher_TC50 = exp(xmid+ResistanceDelta_line1), lower_TC50 = exp(xmid) )
}

start.ResistanceDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)


##########
# A sigmoidal model where line2 sometimes has a delta compared to its baseline, 
# and line1 has a delta compared to whatever value line2 has (ie, regardless of whether line2 has its little delta or not),
# and line1 sometimes has an additional delta.
# Although line1 and line2 can both change, line1 always has a delta compared to line2.
predict.ResistanceDelta_ExperimentDelta = function(params, x, is_line1, is_line2, has_line1_ExperimentDelta, has_line2_ExperimentDelta)
{
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	ResistanceDelta_line1 = params[8]
	ExperimentDelta_line1 = params[9]
	ExperimentDelta_line2 = params[10]
	finalA = (A_line1 * is_line1) + (A_line2 * is_line2)
	finalB = (B_line1 * is_line1) + (B_line2 * is_line2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalDelta = (ResistanceDelta_line1 * is_line1) + (ExperimentDelta_line1 * has_line1_ExperimentDelta * is_line1) + (ExperimentDelta_line2 * has_line2_ExperimentDelta * is_line2)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp(( (xmid+finalDelta) - x) / finalScal)) ) 
	prediction_result
}
loglik.ResistanceDelta_ExperimentDelta = function(params, x, y, is_line1, is_line2, has_line1_ExperimentDelta, has_line2_ExperimentDelta)
{
	sigma = params[11]
	yhat = predict.ResistanceDelta_ExperimentDelta(params, x, is_line1, is_line2, has_line1_ExperimentDelta, has_line2_ExperimentDelta)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.ResistanceDelta_ExperimentDelta = function(params)
{
	params = as.numeric(params)
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	ResistanceDelta_line1 = params[8]
	ExperimentDelta_line1 = params[9]
	ExperimentDelta_line2 = params[10]
	c(higher_TC50_for_resistant_line = exp(xmid+ResistanceDelta_line1+ExperimentDelta_line1), lower_TC50_for_resistant_line = exp(xmid+ResistanceDelta_line1), higher_TC50_for_sensitive_line = exp(xmid+ExperimentDelta_line2), lower_TC50_for_sensitive_line = exp(xmid) )
}

start.ResistanceDelta_ExperimentDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   ExperimentDelta_line1=ExperimentDelta_start,  ExperimentDelta_line2=ExperimentDelta_start,  sigma=sigma_start)
lower.ResistanceDelta_ExperimentDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     ExperimentDelta_line1=ExperimentDelta_low,    ExperimentDelta_line2=ExperimentDelta_low,    sigma=sigma_low)
upper.ResistanceDelta_ExperimentDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    ExperimentDelta_line1=ExperimentDelta_high,   ExperimentDelta_line2=ExperimentDelta_high,   sigma=sigma_high)


##########
# A sigmoidal model where line2 sometimes has a delta compared to its baseline, 
# and line1 has a delta compared to whatever value line2 has (ie, regardless of whether line2 has its little delta or not),
# and line1 sometimes has an additional delta,
# and line1 sometimes has an additional delta that is so large that TC50 was not reached.
# Although line1 and line2 can both change, line1 always has a delta compared to line2.
predict.ResistanceDelta_ExperimentDelta_NotReachedDelta = function(params, x, is_line1, is_line2, has_line1_ExperimentDelta, has_line2_ExperimentDelta, has_line1_NotReachedDelta)
{
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	ResistanceDelta_line1 = params[8]
	ExperimentDelta_line1 = params[9]
	ExperimentDelta_line2 = params[10]
	NotReachedDelta_line1 = params[11]
	finalA = (A_line1 * is_line1) + (A_line2 * is_line2)
	finalB = (B_line1 * is_line1) + (B_line2 * is_line2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalDelta = (ResistanceDelta_line1 * is_line1) + (ExperimentDelta_line1 * has_line1_ExperimentDelta * is_line1) + (ExperimentDelta_line2 * has_line2_ExperimentDelta * is_line2) + (NotReachedDelta_line1 * has_line1_NotReachedDelta * is_line1)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp(( (xmid+finalDelta) - x) / finalScal)) ) 
	prediction_result
}
loglik.ResistanceDelta_ExperimentDelta_NotReachedDelta = function(params, x, y, is_line1, is_line2, has_line1_ExperimentDelta, has_line2_ExperimentDelta, has_line1_NotReachedDelta)
{
	sigma = params[12]
	yhat = predict.ResistanceDelta_ExperimentDelta_NotReachedDelta(params, x, is_line1, is_line2, has_line1_ExperimentDelta, has_line2_ExperimentDelta, has_line1_NotReachedDelta)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.ResistanceDelta_ExperimentDelta_NotReachedDelta = function(params)
{
	params = as.numeric(params)
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	ResistanceDelta_line1 = params[8]
	ExperimentDelta_line1 = params[9]
	ExperimentDelta_line2 = params[10]
	NotReachedDelta_line1 = params[11]
	c(higher_TC50_for_resistant_line = exp(xmid+ResistanceDelta_line1+ExperimentDelta_line1), lower_TC50_for_resistant_line = exp(xmid+ResistanceDelta_line1), higher_TC50_for_sensitive_line = exp(xmid+ExperimentDelta_line2), lower_TC50_for_sensitive_line = exp(xmid), additional_TC50_when_TC50_not_reached = exp(NotReachedDelta_line1) )
}

start.ResistanceDelta_ExperimentDelta_NotReachedDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   ExperimentDelta_line1=ExperimentDelta_start,  ExperimentDelta_line2=ExperimentDelta_start,  NotReachedDelta_line1=NotReachedDelta_start,  sigma=sigma_start)
lower.ResistanceDelta_ExperimentDelta_NotReachedDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     ExperimentDelta_line1=ExperimentDelta_low,    ExperimentDelta_line2=ExperimentDelta_low,    NotReachedDelta_line1=NotReachedDelta_low,    sigma=sigma_low)
upper.ResistanceDelta_ExperimentDelta_NotReachedDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    ExperimentDelta_line1=ExperimentDelta_high,   ExperimentDelta_line2=ExperimentDelta_high,   NotReachedDelta_line1=NotReachedDelta_high,   sigma=sigma_high)


##########
# A sigmoidal model where line1 and/or line2 experiments sometimes have a resistance delta added to the line2 TC50 value.
predict.SometimesResistanceDelta = function(params, x, is_line1, is_line2, has_line1_ResistanceDelta, has_line2_ResistanceDelta)
{
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	ResistanceDelta = params[8]
	finalA = (A_line1 * is_line1) + (A_line2 * is_line2)
	finalB = (B_line1 * is_line1) + (B_line2 * is_line2)
	finalScal = (scal_line1 * is_line1) + (scal_line2 * is_line2)
	finalResistanceDelta = (is_line1 * has_line1_ResistanceDelta * ResistanceDelta) + (is_line2 * has_line2_ResistanceDelta * ResistanceDelta)
	prediction_result = finalA + (finalB - finalA) * ( 1/(1 + exp(( (xmid+finalResistanceDelta) - x) / finalScal)) ) 
	prediction_result
}
loglik.SometimesResistanceDelta = function(params, x, y, is_line1, is_line2, has_line1_ResistanceDelta, has_line2_ResistanceDelta)
{
	sigma = params[9]
	yhat = predict.SometimesResistanceDelta(params, x, is_line1, is_line2, has_line1_ResistanceDelta, has_line2_ResistanceDelta)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.SometimesResistanceDelta = function(params)
{
	params = as.numeric(params)
	A_line1 = params[1]
	A_line2 = params[2]
	B_line1 = params[3]
	B_line2 = params[4]
	scal_line1 = params[5]
	scal_line2 = params[6]
	xmid = params[7]
	ResistanceDelta = params[8]
	c(higher_TC50 = exp(xmid+ResistanceDelta), lower_TC50 = exp(xmid) )
}

start.SometimesResistanceDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   sigma=sigma_start)
lower.SometimesResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.SometimesResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)


##########
# A linear model where line1 and line2 experiments are the same linear models
# y = mx + b
# y = slope * x + yintercept
predict.model1straightline = function(params, x, is_line1, is_line2)
{
	slope = params[1]
	yintercept = params[2]
	final_slope = slope
	final_yintercept = yintercept
	prediction_result = final_slope * x + final_yintercept
	prediction_result
}
loglik.model1straightline = function(params, x, y, is_line1, is_line2)
{
	sigma = params[3]
	yhat = predict.model1straightline(params, x, is_line1, is_line2)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.model1straightline = function(params)
{
	params = as.numeric(params)
	slope = params[1]
	yintercept = params[2]
	c(difference_in_y_between_line1_and_line2 = 0)
}

start.model1straightline  = c(slope=0.1,   yintercept=1,      sigma=0.1)
lower.model1straightline  = c(slope=-2,    yintercept=1e-4,   sigma=1e-4)
upper.model1straightline  = c(slope=2,     yintercept=5,      ysigma=Inf)


##########
# A linear model where line1 and line2 experiments are linear models and line1 has a delta for slope and/or for y-intercept.
# Thus line1 is forced to always have higher y values for viability than line2.
# y = mx + b
# y = slope * x + yintercept
predict.model2straightlines = function(params, x, is_line1, is_line2)
{
	slope = params[1]
	yintercept = params[2]
	delta_slope = params[3]
	delta_yintercept = params[4]
	final_slope = (is_line2 * slope) + (is_line1 * (slope + delta_slope))
	final_yintercept = (is_line2 * yintercept) + (is_line1 * (yintercept + delta_yintercept))
	prediction_result = final_slope * x + final_yintercept
	prediction_result
}
loglik.model2straightlines = function(params, x, y, is_line1, is_line2)
{
	sigma = params[5]
	yhat = predict.model2straightlines(params, x, is_line1, is_line2)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.model2straightlines = function(params)
{
	params = as.numeric(params)
	slope = params[1]
	yintercept = params[2]
	delta_slope = params[3]
	delta_yintercept = params[4]
	c(difference_in_slope_between_line1_and_line2 = delta_slope, difference_in_yintercept_between_line1_and_line2 = delta_yintercept)
}

start.model2straightlines  = c(slope=0.1,   yintercept=1,      delta_slope=-0.1,   delta_yintercept=0.1,   sigma=0.1)
lower.model2straightlines  = c(slope=-2,    yintercept=1e-4,   delta_slope=-0.2,   delta_yintercept=0,     sigma=1e-4)
upper.model2straightlines  = c(slope=2,     yintercept=5,      delta_slope=0,      delta_yintercept=0.2,   sigma=Inf)


#################################################################################################### Cisplatin resistant cell lines vs parental cell lines.
####################################################################################################

################################################################################ args_cellline_chemoResistant_highDose vs args_cellline_chemoResistant_lowDose.
#################### Cisplatin. IL45CisR_highDose vs IL45_high_passage.
assays = rat_Cisplatin; cellline = assays[( ((assays$Line==args_cellline_chemoResistant_highDose)|(assays$Line==args_cellline_chemoResistant_lowDose)) ),]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_chemoResistant_highDose), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_chemoResistant_lowDose), 1, 0)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line()
#
start.model2lines  = c(A_line1=0.6,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)
#
fit.model2lines = optimx(start.model2lines, function(par) loglik.model2lines(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.model2lines, upper=upper.model2lines, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_IL45_high_passage_noResistanceDelta = fit.model2lines$value
fit.model2lines; convert_values_back_to_interpretable.model2lines(fit.model2lines)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.model2lines(coef(fit.model2lines), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.model2lines(coef(fit.model2lines), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1   A_line2   B_line1   B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9965714 0.9577136 0.3623274 0.1364834  0.2884864  0.3236071
	# 	      xmid     sigma   value fevals gevals niter convcode  kkt1  kkt2
	# L-BFGS-B 0.2997228 0.1354483 213.548     42     42    NA        0 FALSE FALSE
	# 	 xtimes
	# L-BFGS-B  0.116
	#     TC50 
	# 1.349485 
#
start.ResistanceDelta  = c(A_line1=0.6,   A_line2=A_start,   B_line1=0.3,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=0.5,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=0.2,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)
#
fit.ResistanceDelta = optimx(start.ResistanceDelta, function(par) loglik.ResistanceDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.ResistanceDelta, upper=upper.ResistanceDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
fit.ResistanceDelta; convert_values_back_to_interpretable.ResistanceDelta(fit.ResistanceDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1  A_line2   B_line1   B_line2 scal_line1 scal_line2      xmid
	# L-BFGS-B 0.9682628 0.978282 0.3081784 0.1773617   0.295392   0.303016 0.1941367
	# 	 ResistanceDelta     sigma   value fevals gevals niter convcode kkt1
	# L-BFGS-B       0.2761193 0.1328767 220.592     71     71    NA        0 TRUE
	# 	  kkt2 xtimes
	# L-BFGS-B FALSE   0.22
	# higher_TC50  lower_TC50 
	#    1.600404    1.214262 
#
degrees_of_freedom=1 # withResistanceDelta has 1 additional variable (ResistanceDelta) that noResistanceDelta doesn't have
pvalue = pchisq(2*(fit.ResistanceDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
IL45CisR_highDose_IL45_high_passage_ResistanceDelta_vs_noResistanceDelta_pvalue = pvalue; pvalue
	# [1] 0.0001744485
	# There is a significant difference between IL45CisR_highDose and IL45_high_passage sensitivity to Cisplatin.

#################### BAMLET. args_cellline_chemoResistant_highDose vs args_cellline_chemoResistant_lowDose.
assays = rat_bamlet; cellline = assays[( ((assays$Line==args_cellline_chemoResistant_highDose)|(assays$Line==args_cellline_chemoResistant_lowDose)) ),]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_chemoResistant_highDose), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_chemoResistant_lowDose), 1, 0)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line()
#
start.model2lines  = c(A_line1=A_start,   A_line2=A_start,   B_line1=0.1,   B_line2=0.1,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=0,     B_line2=0,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)
#
fit.model2lines = optimx(start.model2lines, function(par) loglik.model2lines(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.model2lines, upper=upper.model2lines, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_IL45_high_passage_noResistanceDelta = fit.model2lines$value
fit.model2lines; convert_values_back_to_interpretable.model2lines(fit.model2lines)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.model2lines(coef(fit.model2lines), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.model2lines(coef(fit.model2lines), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1 A_line2    B_line1    B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9954771 1.04633 0.06143478 0.06824846  0.1334612 0.09406978
	# 	      xmid     sigma    value fevals gevals niter convcode  kkt1  kkt2
	# L-BFGS-B -0.308639 0.1605796 141.4707     30     30    NA        0 FALSE FALSE
	# 	 xtimes
	# L-BFGS-B   0.08
	#      TC50 
	# 0.7344459 
#
start.ResistanceDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=0.1,   B_line2=0.1,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=0,     B_line2=0,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)
#
fit.ResistanceDelta = optimx(start.ResistanceDelta, function(par) loglik.ResistanceDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.ResistanceDelta, upper=upper.ResistanceDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
fit.ResistanceDelta; convert_values_back_to_interpretable.ResistanceDelta(fit.ResistanceDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y))
	# 	   A_line1  A_line2    B_line1    B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9954913 1.046325 0.06144025 0.06826086  0.1334617 0.09406941
	# 	       xmid ResistanceDelta     sigma    value fevals gevals niter
	# L-BFGS-B -0.3086414               0 0.1605772 141.4707     46     46    NA
	# 	 convcode  kkt1  kkt2 xtimes
	# L-BFGS-B        0 FALSE FALSE  0.144
	# higher_TC50  lower_TC50 
	#   0.7344441   0.7344441 
#
degrees_of_freedom=1 # withResistanceDelta has 1 additional variable (ResistanceDelta) that noResistanceDelta doesn't have
pvalue = pchisq(2*(fit.ResistanceDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
IL45CisR_highDose_IL45_high_passage_ResistanceDelta_vs_noResistanceDelta_pvalue = pvalue; pvalue
	# [1] 0.9995776
	# There is no significant difference between IL45CisR_highDose and IL45_high_passage sensitivity to BAMLET.

################################################################################ args_cellline_chemoResistant_highDose vs args_cellline_chemoSensitive_lowPassage.
#################### Cisplatin. IL45CisR_highDose vs IL45_low_passage.
assays = rat_Cisplatin; cellline = assays[( ((assays$Line==args_cellline_chemoResistant_highDose)|(assays$Line==args_cellline_chemoSensitive_lowPassage)) ),]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_chemoResistant_highDose), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_chemoSensitive_lowPassage), 1, 0)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line()
#
start.model2lines  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)
#
fit.model2lines = optimx(start.model2lines, function(par) loglik.model2lines(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.model2lines, upper=upper.model2lines, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_IL45_low_passage_noResistanceDelta = fit.model2lines$value
fit.model2lines; convert_values_back_to_interpretable.model2lines(fit.model2lines)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.model2lines(coef(fit.model2lines), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.model2lines(coef(fit.model2lines), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1 A_line2   B_line1   B_line2 scal_line1 scal_line2     xmid
	# L-BFGS-B 0.9754008 0.96539 0.3310762 0.1122777  0.2847296  0.4003645 0.410469
	# 	     sigma    value fevals gevals niter convcode  kkt1  kkt2 xtimes
	# L-BFGS-B 0.1394843 198.3274     47     47    NA        0 FALSE FALSE  0.124
	#     TC50 
	# 1.507525
#
start.ResistanceDelta  = c(A_line1=0.6,   A_line2=A_start,   B_line1=0.3,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=0.5,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=0.2,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)
#
fit.ResistanceDelta = optimx(start.ResistanceDelta, function(par) loglik.ResistanceDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.ResistanceDelta, upper=upper.ResistanceDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45GemR_highDose_IL45_low_passage_20151204_withResistanceDelta = fit.ResistanceDelta$value;
IL45GemR_highDose_IL45_low_passage_20151204_TC50ResistanceDelta = exp(fit.ResistanceDelta$xmid + fit.ResistanceDelta$ResistanceDelta) - exp(fit.ResistanceDelta$xmid)
fit.ResistanceDelta; convert_values_back_to_interpretable.ResistanceDelta(fit.ResistanceDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) + facet_wrap(~ Experiment)
	# 	   A_line1   A_line2   B_line1   B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9682374 0.9718665 0.3082464 0.1471883  0.2952946  0.3781773
	# 	      xmid ResistanceDelta     sigma    value fevals gevals niter
	# L-BFGS-B 0.3412967       0.1288745 0.1390795 199.3701     79     79    NA
	# 	 convcode kkt1  kkt2 xtimes
	# L-BFGS-B        0 TRUE FALSE  0.292
	# higher_TC50  lower_TC50 
	#    1.600268    1.406771 
#
degrees_of_freedom=1 # withResistanceDelta has 1 additional variable (ResistanceDelta) that noResistanceDelta doesn't have
pvalue = pchisq(2*(fit.ResistanceDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
IL45CisR_highDose_IL45_low_passage_ResistanceDelta_vs_noResistanceDelta_pvalue = pvalue; pvalue
	# [1] 0.1487169
	# There is a bit of a difference, however, there is no significant difference between the simpler models for IL45CisR_highDose and IL45_low_passage sensitivity to Cisplatin.
	# The additional plots below show that IL45_low_passage changes in a small way from experiment to experiment,
	# and IL45CisR_highDose changes in a larger way from experiment to experiment,
	# and there is usually a difference between them in 7 out of 8 experiment comparisons.
	# The difference in TC50s changes for each experiment, ie. the size of delta changes for each experiment.
#
cellline2 = cellline[ (cellline$Line==args_cellline_chemoResistant_highDose), ]; ggplot(cellline2, aes(x=log10(Concentration), y=Viability, col=Experiment)) + geom_line()
cellline2 = cellline[ (cellline$Line==args_cellline_chemoSensitive_lowPassage), ]; ggplot(cellline2, aes(x=log10(Concentration), y=Viability, col=Experiment)) + geom_line()
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + facet_wrap(~ Experiment)
#
start.ResistanceDelta_ExperimentDelta  = c(A_line1=0.6,   A_line2=A_start,   B_line1=0.3,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=0.1,    ExperimentDelta_line1=0.1,      ExperimentDelta_line2=0.1,   sigma=sigma_start)
lower.ResistanceDelta_ExperimentDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=0.2,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=-0.5,   ExperimentDelta_line1=-0.5,     ExperimentDelta_line2=-0.5,  sigma=sigma_low)
upper.ResistanceDelta_ExperimentDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=0.5,    ExperimentDelta_line1=0.5,      ExperimentDelta_line2=0.5,   sigma=sigma_high)
#
cellline$has_line1_ExperimentDelta <- ifelse( ((cellline$is_line1==1) & (cellline$Experiment %in% c('20151204','20151211','20151218','20160101','20160122','20160129','20160206'))), 1, 0)
cellline$has_line2_ExperimentDelta <- ifelse( ((cellline$is_line2==1) & (cellline$Experiment %in% c('20151204','20151211','20151218','20160101','20160122','20160129'))), 1, 0)
#
fit.ResistanceDelta_ExperimentDelta = optimx(start.ResistanceDelta_ExperimentDelta, function(par) loglik.ResistanceDelta_ExperimentDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2, has_line1_ExperimentDelta=cellline$has_line1_ExperimentDelta, has_line2_ExperimentDelta=cellline$has_line2_ExperimentDelta), lower=lower.ResistanceDelta_ExperimentDelta, upper=upper.ResistanceDelta_ExperimentDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
#
fit.ResistanceDelta_ExperimentDelta; convert_values_back_to_interpretable.ResistanceDelta_ExperimentDelta(fit.ResistanceDelta_ExperimentDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta_ExperimentDelta(coef(fit.ResistanceDelta_ExperimentDelta), fit1_x, 1, 0, 0, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta_ExperimentDelta(coef(fit.ResistanceDelta_ExperimentDelta), fit2_x, 0, 1, 0, 0); fit2=data.frame(x=fit2_x,y=fit2_y)
fit3_x = seq(-2, 2, 0.01); fit3_y = predict.ResistanceDelta_ExperimentDelta(coef(fit.ResistanceDelta_ExperimentDelta), fit3_x, 1, 0, 1, 0); fit3=data.frame(x=fit3_x,y=fit3_y)
fit4_x = seq(-2, 2, 0.01); fit4_y = predict.ResistanceDelta_ExperimentDelta(coef(fit.ResistanceDelta_ExperimentDelta), fit4_x, 0, 1, 0, 1); fit4=data.frame(x=fit4_x,y=fit4_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="orange",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) + geom_line(data=fit3,color="red",aes(x=x,y=y)) + geom_line(data=fit4,color="blue",aes(x=x,y=y)) + facet_wrap(~ Experiment)
	# 	   A_line1   A_line2   B_line1   B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9856899 0.9625537 0.3032535 0.1602095  0.3060391  0.3174066
	# 	        xmid ResistanceDelta ExperimentDelta_line1
	# L-BFGS-B -0.01708188     -0.01798992                   0.5
	# 	 ExperimentDelta_line2     sigma    value fevals gevals niter convcode
	# L-BFGS-B             0.4475287 0.1283659 228.2325    104    104    NA        0
	# 	  kkt1  kkt2 xtimes
	# L-BFGS-B FALSE FALSE   0.68
	# higher_TC50_for_resistant_line  lower_TC50_for_resistant_line 
	# 	             1.5918999                      0.9655361 
	# higher_TC50_for_sensitive_line  lower_TC50_for_sensitive_line 
	# 	             1.5379446                      0.9830632 
#
degrees_of_freedom=2 # Their are 2 extra parameters - the line1_experiment_delta and line2_experiment_delta
pvalue = pchisq(2*(fit.ResistanceDelta_ExperimentDelta$value - fit.ResistanceDelta$value), degrees_of_freedom, lower.tail = FALSE)
pvalue
	# [1] 2.918689e-13
degrees_of_freedom=3 # Their are 3 extra parameters - the resistance_delta, line1_experiment_delta and line2_experiment_delta
pvalue = pchisq(2*(fit.ResistanceDelta_ExperimentDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
pvalue
	# [1] 6.45316e-13
	# The model that allows differences in TC50 for different experiments, 
	# yet always maintains a difference in TC50 between the resistant and parental cell lines,
	# is statistically significantly better than the other models.
	# There is a statistically significant difference between IL45CisR_highDose and IL45_low_passage sensitivity to Cisplatin.

#################### BAMLET. args_cellline_chemoResistant_highDose vs args_cellline_chemoSensitive_lowPassage.
assays = rat_bamlet; cellline = assays[( ((assays$Line==args_cellline_chemoResistant_highDose)|(assays$Line==args_cellline_chemoSensitive_lowPassage)) ),]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_chemoResistant_highDose), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_chemoSensitive_lowPassage), 1, 0)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line()
#
start.model2lines  = c(A_line1=A_start,   A_line2=A_start,   B_line1=0.2,   B_line2=0.2,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=0,     B_line2=0,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)
#
fit.model2lines = optimx(start.model2lines, function(par) loglik.model2lines(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.model2lines, upper=upper.model2lines, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_IL45_low_passage_noResistanceDelta = fit.model2lines$value
fit.model2lines; convert_values_back_to_interpretable.model2lines(fit.model2lines)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.model2lines(coef(fit.model2lines), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.model2lines(coef(fit.model2lines), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1  A_line2    B_line1    B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9904435 1.005135 0.05435726 0.06773019  0.1355122 0.07555994
	# 	       xmid     sigma    value fevals gevals niter convcode kkt1  kkt2
	# L-BFGS-B -0.2928708 0.1497447 165.5722     32     32    NA        0 TRUE FALSE
	# 	 xtimes
	# L-BFGS-B  0.088
	#      TC50 
	# 0.7461186 
#
start.ResistanceDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=0.1,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=0,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=0.25,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)
#
fit.ResistanceDelta = optimx(start.ResistanceDelta, function(par) loglik.ResistanceDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.ResistanceDelta, upper=upper.ResistanceDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45GemR_highDose_IL45_low_passage_20151204_withResistanceDelta = fit.ResistanceDelta$value;
IL45GemR_highDose_IL45_low_passage_20151204_TC50ResistanceDelta = exp(fit.ResistanceDelta$xmid + fit.ResistanceDelta$ResistanceDelta) - exp(fit.ResistanceDelta$xmid)
fit.ResistanceDelta; convert_values_back_to_interpretable.ResistanceDelta(fit.ResistanceDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1  A_line2    B_line1    B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9904475 1.005137 0.05434792 0.06773486  0.1355167 0.07555857
	# 	       xmid ResistanceDelta     sigma    value fevals gevals niter
	# L-BFGS-B -0.2928686               0 0.1497437 165.5722     77     77    NA
	# 	 convcode  kkt1  kkt2 xtimes
	# L-BFGS-B        0 FALSE FALSE  0.232
	# higher_TC50  lower_TC50 
	#   0.7461202   0.7461202 
#
degrees_of_freedom=1 # withResistanceDelta has 1 additional variable (ResistanceDelta) that noResistanceDelta doesn't have
pvalue = pchisq(2*(fit.ResistanceDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
IL45CisR_highDose_IL45_low_passage_ResistanceDelta_vs_noResistanceDelta_pvalue = pvalue; pvalue
	# [1] 0.9994573
	# There is no significant difference between IL45CisR_highDose and IL45_low_passage sensitivity to BAMLET.

################################################################################ args_cellline_chemoSensitive_highPassage vs args_cellline_chemoResistant_lowDose.
#################### Cisplatin. IL45CisR_lowDose vs IL45_high_passage.
assays = rat_Cisplatin; cellline = assays[( ((assays$Line==args_cellline_chemoSensitive_highPassage)|(assays$Line==args_cellline_chemoResistant_lowDose)) ),]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_chemoSensitive_highPassage), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_chemoResistant_lowDose), 1, 0)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line()
#
start.model2lines  = c(A_line1=A_start,   A_line2=A_start,   B_line1=0.3,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=0.5,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)
#
fit.model2lines = optimx(start.model2lines, function(par) loglik.model2lines(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.model2lines, upper=upper.model2lines, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_IL45_high_passage_noResistanceDelta = fit.model2lines$value
fit.model2lines; convert_values_back_to_interpretable.model2lines(fit.model2lines)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.model2lines(coef(fit.model2lines), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.model2lines(coef(fit.model2lines), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1  A_line2   B_line1   B_line2 scal_line1 scal_line2    xmid
	# L-BFGS-B 0.9769283 0.954067 0.4112523 0.1161993  0.3410254   0.340072 0.33937
	# 	     sigma    value fevals gevals niter convcode  kkt1  kkt2 xtimes
	# L-BFGS-B 0.1408826 209.8874     47     47    NA        0 FALSE FALSE  0.136
	#     TC50 
	# 1.404063 
#
start.ResistanceDelta  = c(A_line1=1.1,   A_line2=A_start,   B_line1=0.3,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=0.85,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)
#
fit.ResistanceDelta = optimx(start.ResistanceDelta, function(par) loglik.ResistanceDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.ResistanceDelta, upper=upper.ResistanceDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
fit.ResistanceDelta; convert_values_back_to_interpretable.ResistanceDelta(fit.ResistanceDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1   A_line2   B_line1   B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9338153 0.9782644 0.2550204 0.1773719  0.3763126  0.3029838
	# 	      xmid ResistanceDelta    sigma    value fevals gevals niter
	# L-BFGS-B 0.1941506       0.5574779 0.134879 226.7821     75     75    NA
	# 	 convcode  kkt1  kkt2 xtimes
	# L-BFGS-B        0 FALSE FALSE   0.24
	# higher_TC50  lower_TC50 
	#    2.120451    1.214279 
#
degrees_of_freedom=1 # withResistanceDelta has 1 additional variable (ResistanceDelta) that noResistanceDelta doesn't have
pvalue = pchisq(2*(fit.ResistanceDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
IL45CisR_lowDose_IL45_high_passage_ResistanceDelta_vs_noResistanceDelta_pvalue = pvalue; pvalue
	# [1] 6.14159e-09
	# There is a significant difference between IL45CisR_lowDose and IL45_high_passage sensitivity to Cisplatin.

#################### BAMLET. args_cellline_chemoSensitive_highPassage vs args_cellline_chemoResistant_lowDose.
assays = rat_bamlet; cellline = assays[( ((assays$Line==args_cellline_chemoSensitive_highPassage)|(assays$Line==args_cellline_chemoResistant_lowDose)) ),]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_chemoSensitive_highPassage), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_chemoResistant_lowDose), 1, 0)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line()
#
start.model2lines  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=0.2,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=0,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)
#
fit.model2lines = optimx(start.model2lines, function(par) loglik.model2lines(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.model2lines, upper=upper.model2lines, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_IL45_high_passage_noResistanceDelta = fit.model2lines$value
fit.model2lines; convert_values_back_to_interpretable.model2lines(fit.model2lines)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.model2lines(coef(fit.model2lines), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.model2lines(coef(fit.model2lines), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	  A_line1  A_line2   B_line1   B_line2 scal_line1 scal_line2      xmid
	# L-BFGS-B 1.031595 1.031858 0.1060116 0.0590163  0.1422762 0.08231622 -0.275298
	# 	     sigma    value fevals gevals niter convcode kkt1  kkt2 xtimes
	# L-BFGS-B 0.1647676 146.0381     33     33    NA        0 TRUE FALSE  0.096
	#      TC50 
	# 0.7593458 
#
start.ResistanceDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)
#
fit.ResistanceDelta = optimx(start.ResistanceDelta, function(par) loglik.ResistanceDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.ResistanceDelta, upper=upper.ResistanceDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
fit.ResistanceDelta; convert_values_back_to_interpretable.ResistanceDelta(fit.ResistanceDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	  A_line1  A_line2    B_line1    B_line2 scal_line1 scal_line2
	# L-BFGS-B 1.025336 1.037817 0.09529955 0.06177755  0.1446478 0.08838678
	# 	       xmid ResistanceDelta     sigma    value fevals gevals niter
	# L-BFGS-B -0.2881194       0.0342574 0.1644785 146.7074     65     65    NA
	# 	 convcode  kkt1  kkt2 xtimes
	# L-BFGS-B        0 FALSE FALSE  0.236
	# higher_TC50  lower_TC50 
	#   0.7757988   0.7496721 
#
degrees_of_freedom=1 # withResistanceDelta has 1 additional variable (ResistanceDelta) that noResistanceDelta doesn't have
pvalue = pchisq(2*(fit.ResistanceDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
IL45CisR_lowDose_IL45_high_passage_ResistanceDelta_vs_noResistanceDelta_pvalue = pvalue; pvalue
	# [1] 0.2472679
	# There is no significant difference between IL45CisR_lowDose and IL45_high_passage sensitivity to BAMLET.
#
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + facet_wrap(~ Experiment)

################################################################################ args_cellline_chemoSensitive_highPassage vs args_cellline_chemoSensitive_lowPassage.
#################### Cisplatin. IL45CisR_lowDose vs IL45_low_passage.
assays = rat_Cisplatin; cellline = assays[( ((assays$Line==args_cellline_chemoSensitive_highPassage)|(assays$Line==args_cellline_chemoSensitive_lowPassage)) ),]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_chemoSensitive_highPassage), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_chemoSensitive_lowPassage), 1, 0)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line()
#
start.model2lines  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)
#
fit.model2lines = optimx(start.model2lines, function(par) loglik.model2lines(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.model2lines, upper=upper.model2lines, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_IL45_low_passage_noResistanceDelta = fit.model2lines$value
fit.model2lines; convert_values_back_to_interpretable.model2lines(fit.model2lines)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.model2lines(coef(fit.model2lines), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.model2lines(coef(fit.model2lines), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1   A_line2   B_line1    B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9462757 0.9725252 0.3559036 0.02227172  0.3257688   0.479169
	# 	      xmid     sigma    value fevals gevals niter convcode  kkt1  kkt2
	# L-BFGS-B 0.5368871 0.1430855 199.6604     57     57    NA        0 FALSE FALSE
	# 	 xtimes
	# L-BFGS-B  0.164
	#     TC50 
	# 1.710673 
#
start.ResistanceDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=0.5,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)
#
fit.ResistanceDelta = optimx(start.ResistanceDelta, function(par) loglik.ResistanceDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.ResistanceDelta, upper=upper.ResistanceDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45GemR_highDose_IL45_low_passage_20151204_withResistanceDelta = fit.ResistanceDelta$value;
IL45GemR_highDose_IL45_low_passage_20151204_TC50ResistanceDelta = exp(fit.ResistanceDelta$xmid + fit.ResistanceDelta$ResistanceDelta) - exp(fit.ResistanceDelta$xmid)
fit.ResistanceDelta; convert_values_back_to_interpretable.ResistanceDelta(fit.ResistanceDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	   A_line1   A_line2   B_line1   B_line2 scal_line1 scal_line2
	# L-BFGS-B 0.9338046 0.9718535 0.2549534 0.1472879   0.376324  0.3780898
	# 	      xmid ResistanceDelta     sigma   value fevals gevals niter
	# L-BFGS-B 0.3411634       0.4105923 0.1407176 205.998     73     73    NA
	# 	 convcode  kkt1  kkt2 xtimes
	# L-BFGS-B        0 FALSE FALSE  0.232
	# higher_TC50  lower_TC50 
	#    2.120720    1.406583 
#
degrees_of_freedom=1 # withResistanceDelta has 1 additional variable (ResistanceDelta) that noResistanceDelta doesn't have
pvalue = pchisq(2*(fit.ResistanceDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
IL45CisR_lowDose_IL45_low_passage_ResistanceDelta_vs_noResistanceDelta_pvalue = pvalue; pvalue
	# [1] 0.0003705024
	# There is a significant difference between IL45CisR_lowDose and IL45_low_passage sensitivity to Cisplatin.

#################### BAMLET. args_cellline_chemoSensitive_highPassage vs args_cellline_chemoSensitive_lowPassage.
assays = rat_bamlet; cellline = assays[( ((assays$Line==args_cellline_chemoSensitive_highPassage)|(assays$Line==args_cellline_chemoSensitive_lowPassage)) ),]
cellline$is_line1 <- ifelse( (cellline$Line==args_cellline_chemoSensitive_highPassage), 1, 0)
cellline$is_line2 <- ifelse( (cellline$Line==args_cellline_chemoSensitive_lowPassage), 1, 0)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line()
#
start.model2lines  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=0.2,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.model2lines  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=0,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.model2lines  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    sigma=sigma_high)
#
fit.model2lines = optimx(start.model2lines, function(par) loglik.model2lines(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.model2lines, upper=upper.model2lines, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_IL45_low_passage_noResistanceDelta = fit.model2lines$value
fit.model2lines; convert_values_back_to_interpretable.model2lines(fit.model2lines)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.model2lines(coef(fit.model2lines), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.model2lines(coef(fit.model2lines), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	  A_line1  A_line2    B_line1    B_line2 scal_line1 scal_line2
	# L-BFGS-B 1.059285 1.011005 0.09445872 0.05195412  0.1618677       0.02
	# 	       xmid     sigma   value fevals gevals niter convcode  kkt1  kkt2
	# L-BFGS-B -0.2487896 0.1518233 153.042     16     16    NA        0 FALSE FALSE
	# 	 xtimes
	# L-BFGS-B  0.044
	#      TC50 
	# 0.7797441 
#
start.ResistanceDelta  = c(A_line1=A_start,   A_line2=A_start,   B_line1=B_start,   B_line2=B_start,   scal_line1=scal_start,   scal_line2=scal_start,   xmid=xmid_start,   ResistanceDelta=ResistanceDelta_start,   sigma=sigma_start)
lower.ResistanceDelta  = c(A_line1=A_low,     A_line2=A_low,     B_line1=B_low,     B_line2=B_low,     scal_line1=scal_low,     scal_line2=scal_low,     xmid=xmid_low,     ResistanceDelta=ResistanceDelta_low,     sigma=sigma_low)
upper.ResistanceDelta  = c(A_line1=A_high,    A_line2=A_high,    B_line1=B_high,    B_line2=B_high,    scal_line1=scal_high,    scal_line2=scal_high,   xmid=xmid_high,    ResistanceDelta=ResistanceDelta_high,    sigma=sigma_high)
#
fit.ResistanceDelta = optimx(start.ResistanceDelta, function(par) loglik.ResistanceDelta(par, x=cellline$log10Conc, y=cellline$Viability, is_line1=cellline$is_line1, is_line2=cellline$is_line2), lower=lower.ResistanceDelta, upper=upper.ResistanceDelta, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45GemR_highDose_IL45_low_passage_20151204_withResistanceDelta = fit.ResistanceDelta$value;
IL45GemR_highDose_IL45_low_passage_20151204_TC50ResistanceDelta = exp(fit.ResistanceDelta$xmid + fit.ResistanceDelta$ResistanceDelta) - exp(fit.ResistanceDelta$xmid)
fit.ResistanceDelta; convert_values_back_to_interpretable.ResistanceDelta(fit.ResistanceDelta)
fit1_x = seq(-2, 2, 0.01); fit1_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit1_x, 1, 0); fit1=data.frame(x=fit1_x,y=fit1_y)
fit2_x = seq(-2, 2, 0.01); fit2_y = predict.ResistanceDelta(coef(fit.ResistanceDelta), fit2_x, 0, 1); fit2=data.frame(x=fit2_x,y=fit2_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit1,color="red",aes(x=x,y=y)) + geom_line(data=fit2,color="cyan",aes(x=x,y=y)) 
	# 	  A_line1   A_line2    B_line1    B_line2 scal_line1 scal_line2
	# L-BFGS-B 1.024872 0.9926535 0.09432556 0.06122137  0.1449987 0.04947704
	# 	       xmid ResistanceDelta     sigma    value fevals gevals niter
	# L-BFGS-B -0.2520602               0 0.1539756 171.7817     71     71    NA
	# 	 convcode  kkt1  kkt2 xtimes
	# L-BFGS-B        0 FALSE FALSE  0.256
	# higher_TC50  lower_TC50 
	#    0.777198    0.777198 
#
degrees_of_freedom=1 # withResistanceDelta has 1 additional variable (ResistanceDelta) that noResistanceDelta doesn't have
pvalue = pchisq(2*(fit.ResistanceDelta$value - fit.model2lines$value), degrees_of_freedom, lower.tail = FALSE)
IL45CisR_lowDose_IL45_low_passage_ResistanceDelta_vs_noResistanceDelta_pvalue = pvalue; pvalue
	# [1] 9.238469e-10
	# The model where the TC50 of IL45CisR_lowDose is different to and higher than TC50 of IL45_low_passage,
	# ie., IL45CisR_lowDose's TC50 = IL45_low_passage's TC50 + a delta value
	# is significantly better than modelling the TC50 as the same.
	# The difference (ie, the delta) is 0.
	# Thus, there is no significant difference between IL45CisR_lowDose and IL45_low_passage sensitivity to BAMLET.


