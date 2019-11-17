# Rscript step05_TC50_and_plots_chemo_resistant_vs_sensitive.R <input_file_1> <input_file_2> <cell_line> > <output_file>
# Rscript step05_TC50_and_plots_chemo_resistant_vs_sensitive.R rat_chemo_cell_assay_data_1.txt rat_chemo_cell_assay_data_2.txt rat_bamlet_cell_assay_data.txt IL45CisR_highDose IL45_high_passage IL45CisR_lowDose IL45_low_passage > step04_output.txt

# Hypothesis for step04: That chemo-resistant cells are just as sensitive to BAMLET as chemo-sensitive cells.
# Step05 are the plots and TC50 calculations for the step04 work.
# First establish that chemo-resistant cells are in fact less sensitive to chemo (cisplatin) than chemo-sensitive cells.
# Also establish that there is no difference in sensitivity (to chemo or BAMLET) between low-passage and high-passage cells,
# so counter the possibility that chemo-resistant cells are actually resistant because they are usually high-passage cells.
# Additionally, the chemo-resistant cells had their chemo-resistance maintained by regular chemotherapy treatment that was either high-dose or low-lose.
# These cell lines have a standard uniphasic response, not a biphasic response. No need to model 2component models for different TC50 in same cell line.

# Please note, the code for plots and TC50 calculations in this program need to be parameterised 
# to use the input cell lines and the list of unique experiments found in the data.

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

rat_cisplatin = rat_chemo[ (rat_chemo$Drug=="Cisplatin"), ]
rat_pemetrexed = rat_chemo[ (rat_chemo$Drug=="Pemetrexed"), ]
rat_gemcitabine = rat_chemo[ (rat_chemo$Drug=="Gemcitabine"), ]
rat_vinorelbine = rat_chemo[ (rat_chemo$Drug=="Vinorelbine"), ]

assays = rat_cisplatin
rat_cisplatin_experiments = unique(rat_cisplatin$Experiment)
ggplot(assays, aes(x = log10(Concentration), y = Viability, col = Experiment)) + geom_line() + facet_wrap(~ Line)

rat = read.table(args_infile_treatment, sep = "\t", header = TRUE)

rat = melt(rat, id.vars = c("Drug", "Line", "Experiment"), variable.name = "Concentration", value.name = "Viability")
rat$Concentration = as.numeric(gsub("^X", "", rat$Concentration)
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
rat_bamlet_experiments = unique(rat_bamlet$Experiment)
ggplot(assays, aes(x = log10(Concentration), y = Viability, col = Experiment)) + geom_point() + facet_wrap(~ Line)


####################################################################################################
####################################################################################################

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
extract_TC50 = function(params)
{
	A = params[1]
	B = params[2]
	scal = params[3]
	xmid = params[4]
	TC50 = exp(xmid - log( (B - A)/(0.5 - A) - 1 ))
	TC50
}
extract_TC50_stdev = function(params, sigma)
{
	A = params[1]
	B = params[2]
	scal = params[3]
	xmid = params[4]
	sigma = params[5]
	xmid_1 = xmid - sigma
	xmid_2 = xmid + sigma
	TC50_1 = exp(xmid_1 - log( (B - A)/(0.5 - A) - 1 ))
	TC50_2 = exp(xmid_2 - log( (B - A)/(0.5 - A) - 1 ))
	TC50_stdev = abs(TC50_2 - TC50_1)/2
	TC50_stdev
}

##########
A_start = 1;            A_low = -Inf;            A_high = Inf;
B_start = 0;            B_low = -Inf;            B_high = Inf;
scal_start = 0.1;       scal_low = 0.02;         scal_high = Inf;
xmid_start = 0;         xmid_low = -Inf;         xmid_high = Inf;
sens_start = 0.5;       sens_low = 0;            sens_high = 1;
xshift_start = 1;       xshift_low = 0;          xshift_high = Inf;
sigma_start = 0.1;      sigma_low = 1e-4;        sigma_high = Inf;
delta_start = 1;        delta_low = 0;           delta_high = Inf;

start.1component  = c(A=A_start,   B=B_start,   scal=scal_start,   xmid=xmid_start,   sigma=sigma_start)
lower.1component  = c(A=A_low,     B=B_low,     scal=scal_low,     xmid=xmid_low,     sigma=sigma_low)
upper.1component  = c(A=A_high,    B=B_high,    scal=scal_high,    xmid=xmid_high,    sigma=sigma_high)

####################################################################################################
extract_TC50_from_lm = function(model)
{
	model_yintercept = model$coef[[1]]
	model_slope = model$coef[[2]]
	TC50 = exp((0.5 - model_yintercept) / model_slope)
	TC50
}

extract_TC50_stdev_from_lm = function(model)
{
	model_yintercept = model$coef[[1]]
	model_slope = model$coef[[2]]
	model_deviance = deviance(model)
	yintercept_1 = model_yintercept - model_deviance
	yintercept_2 = model_yintercept + model_deviance
	TC50_1 = exp((0.5 - yintercept_1) / model_slope)
	TC50_2 = exp((0.5 - yintercept_2) / model_slope)
	TC50_stdev = abs(TC50_2 - TC50_1)/2
	TC50_stdev
}
####################################################################################################

################################################################################ Cisplatin
######################################## Cisplatin. Experiment 20151204
#################### Cisplatin. IL45CisR_highDose. Experiment 20151204
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20151204")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20151204_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20151204_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45CisR_lowDose. Experiment 20151204
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20151204")),]
start.1component = c(A=0.9, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=0.85,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=0.95,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20151204_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20151204_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_high_passage. Experiment 20151204
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20151204")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20151204_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20151204_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_low_passage. Experiment 20151204
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20151204")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20151204_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20151204_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

######################################## Cisplatin. Experiment 20151211
#################### Cisplatin. IL45CisR_highDose. Experiment 20151211
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20151211")),]
start.1component = c(A=0.95, B=B_start, scal=scal_start, xmid=0.8, sigma=sigma_start)
lower.1component = c(A=0.8,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=1.1,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20151211_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20151211_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45CisR_lowDose. Experiment 20151211
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20151211")),]
start.1component = c(A=1, B=B_start, scal=scal_start, xmid=0.7, sigma=sigma_start)
lower.1component = c(A=0.85,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=1.15,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20151211_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20151211_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_high_passage. Experiment 20151211
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20151211")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20151211_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20151211_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_low_passage. Experiment 20151211
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20151211")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20151211_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20151211_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ Cisplatin
######################################## Cisplatin. Experiment 20151218
#################### Cisplatin. IL45CisR_highDose. Experiment 20151218
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20151218")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20151218_Cisplatin_TC50 = 100 # TC50 not reached in this experiment
IL45CisR_highDose_20151218_Cisplatin_TC50_stdev = 0
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45CisR_lowDose. Experiment 20151218
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20151218")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20151218_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20151218_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_high_passage. Experiment 20151218
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20151218")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20151218_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20151218_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_low_passage. Experiment 20151218
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20151218")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20151218_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20151218_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ Cisplatin
######################################## Cisplatin. Experiment 20151227
#################### Cisplatin. IL45CisR_highDose. Experiment 20151227
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20151227")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20151227_Cisplatin_TC50 = 100 # TC50 not reached in this experiment
IL45CisR_highDose_20151227_Cisplatin_TC50_stdev = 0
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45CisR_lowDose. Experiment 20151227
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20151227")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20151227_Cisplatin_TC50 = 100 # TC50 not reached in this experiment
IL45CisR_lowDose_20151227_Cisplatin_TC50_stdev = 0
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_high_passage. Experiment 20151227
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20151227")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20151227_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20151227_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_low_passage. Experiment 20151227
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20151227")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20151227_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20151227_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ Cisplatin
######################################## Cisplatin. Experiment 20160101
#################### Cisplatin. IL45CisR_highDose. Experiment 20160101
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160101")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160101_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160101_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45CisR_lowDose. Experiment 20160101
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160101")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160101_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20160101_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_high_passage. Experiment 20160101
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20160101")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160101_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160101_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_low_passage. Experiment 20160101
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160101")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160101_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160101_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ Cisplatin
######################################## Cisplatin. Experiment20160122
#################### Cisplatin. IL45CisR_highDose. Experiment20160122
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160122")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160122_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160122_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45CisR_lowDose. Experiment20160122
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160122")),]
start.1component = c(A=1.1, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=1.0,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=1.2,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160122_Cisplatin_TC50 = 100 # TC50 not reached in this experiment
IL45CisR_lowDose_20160122_Cisplatin_TC50_stdev = 0
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_high_passage. Experiment20160122
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20160122")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=0.5, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160122_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160122_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_low_passage. Experiment20160122
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160122")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160122_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160122_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ Cisplatin
######################################## Cisplatin. Experiment 20160129
#################### Cisplatin. IL45CisR_highDose. Experiment 20160129
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160129")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=0.8, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160129_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160129_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45CisR_lowDose. Experiment 20160129
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160129")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160129_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20160129_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_high_passage. Experiment 20160129
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20160129")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160129_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160129_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_low_passage. Experiment 20160129
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160129")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160129_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160129_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ Cisplatin
######################################## Cisplatin. Experiment 20160206
#################### Cisplatin. IL45CisR_highDose. Experiment 20160206
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160206")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160206_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160206_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45CisR_lowDose. Experiment 20160206
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160206")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160206_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20160206_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_high_passage. Experiment 20160206
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20160206")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160206_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160206_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### Cisplatin. IL45_low_passage. Experiment 20160206
assays = rat_cisplatin; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160206")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160206_Cisplatin_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160206_Cisplatin_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ BAMLET
######################################## BAMLET. Experiment 20151211
#################### BAMLET. IL45CisR_highDose. Experiment 20151211
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20151211")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20151211_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20151211_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45CisR_lowDose. Experiment 20151211
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20151211")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20151211_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20151211_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_high_passage. Experiment 20151211
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20151211")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.3, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20151211_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20151211_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_low_passage. Experiment 20151211
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20151211")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.3, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20151211_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20151211_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ BAMLET
######################################## BAMLET. Experiment 20151218
#################### BAMLET. IL45CisR_lowDose. Experiment 20151218
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20151218")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20151218_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20151218_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_high_passage. Experiment 20151218
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20151218")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=0, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20151218_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20151218_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_low_passage. Experiment 20151218
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20151218")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=0, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20151218_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20151218_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ BAMLET
######################################## BAMLET. Experiment 20151227
#################### BAMLET. IL45CisR_highDose. Experiment 20151227
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20151227")),]
start.1component = c(A=A_start, B=B_start, scal=0.1, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20151227_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20151227_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45CisR_lowDose. Experiment 20151227
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20151227")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20151227_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20151227_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_high_passage. Experiment 20151227
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20151227")),]
start.1component = c(A=A_start, B=B_start, scal=0.1, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20151227_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20151227_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_low_passage. Experiment 20151227
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20151227")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20151227_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20151227_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ BAMLET
######################################## BAMLET. Experiment 20160101
#################### BAMLET. IL45CisR_highDose. Experiment 20160101
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160101")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160101_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160101_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45CisR_lowDose. Experiment 20160101
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160101")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=0,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160101_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20160101_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_high_passage. Experiment 20160101
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20160101")),]
start.1component = c(A=A_start, B=B_start, scal=-0.4, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160101_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160101_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_low_passage. Experiment 20160101
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160101")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160101_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160101_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ BAMLET
######################################## BAMLET. Experiment 20160108
#################### BAMLET. IL45CisR_highDose. Experiment 20160108
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160108")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.6, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160108_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160108_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45CisR_lowDose. Experiment 20160108
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160108")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=-0.5, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160108_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20160108_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_high_passage. Experiment 20160108
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20160108")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160108_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160108_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_low_passage. Experiment 20160108
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160108")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160108_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160108_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ BAMLET
######################################## BAMLET. Experiment 20160122
#################### BAMLET. IL45CisR_highDose. Experiment 20160122
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160122")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160122_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160122_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45CisR_lowDose. Experiment 20160122
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160122")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=-0.3, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160122_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20160122_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
t = c(A=A_start, B=B_start, scal=0.5, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160122_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160122_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_low_passage. Experiment 20160122
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160122")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160122_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160122_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ BAMLET
######################################## BAMLET. Experiment 20160129
#################### BAMLET. IL45CisR_highDose. Experiment 20160129
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160129")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160129_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160129_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45CisR_lowDose. Experiment 20160129
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160129")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160129_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20160129_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_high_passage. Experiment 20160129
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20160129")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.3, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160129_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160129_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_low_passage. Experiment 20160129
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160129")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160129_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160129_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ BAMLET
######################################## BAMLET. Experiment 20160206
#################### BAMLET. IL45CisR_highDose. Experiment 20160206
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_highDose") & (assays$Experiment=="20160206")),]
start.1component = c(A=A_start, B=B_start, scal=scal_start, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_highDose_20160206_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_highDose_20160206_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45CisR_lowDose. Experiment 20160206
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45CisR_lowDose") & (assays$Experiment=="20160206")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=xmid_start, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45CisR_lowDose_20160206_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45CisR_lowDose_20160206_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_high_passage. Experiment 20160206
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_high_passage") & (assays$Experiment=="20160206")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_high_passage_20160206_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_high_passage_20160206_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))
#################### BAMLET. IL45_low_passage. Experiment 20160206
assays = rat_bamlet; cellline = assays[((assays$Line=="IL45_low_passage") & (assays$Experiment=="20160206")),]
start.1component = c(A=A_start, B=B_start, scal=0.5, xmid=-0.4, sigma=sigma_start)
lower.1component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sigma=sigma_low)
upper.1component = c(A=A_high,  B=B_high,  scal=scal_high,  xmid=xmid_high,  sigma=sigma_high)
fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
IL45_low_passage_20160206_BAMLET_TC50 = unname(unlist(extract_TC50(fit.1component)))
IL45_low_passage_20160206_BAMLET_TC50_stdev = unname(unlist(extract_TC50_stdev(fit.1component)))
fit_x = seq(-1, 2, 0.01); fit_y = predict.1component(coef(fit.1component), fit_x); fit=data.frame(x=fit_x,y=fit_y)
ggplot(cellline, aes(x=log10(Concentration), y=Viability, col=Line)) + geom_line() + geom_line(data=fit,color="red",aes(x=x,y=y))

################################################################################ Data frames for Plots
##### Revived 20151002

cellline_TC50 = c(IL45CisR_highDose_20151204_Cisplatin_TC50, IL45CisR_highDose_20151211_Cisplatin_TC50, IL45CisR_highDose_20151218_Cisplatin_TC50, IL45CisR_highDose_20151227_Cisplatin_TC50, IL45CisR_highDose_20160101_Cisplatin_TC50, IL45CisR_highDose_20160122_Cisplatin_TC50, IL45CisR_highDose_20160129_Cisplatin_TC50, IL45CisR_highDose_20160206_Cisplatin_TC50)
cellline_TC50_stdev = c(IL45CisR_highDose_20151204_Cisplatin_TC50_stdev, IL45CisR_highDose_20151211_Cisplatin_TC50_stdev, IL45CisR_highDose_20151218_Cisplatin_TC50_stdev, IL45CisR_highDose_20151227_Cisplatin_TC50_stdev, IL45CisR_highDose_20160101_Cisplatin_TC50_stdev, IL45CisR_highDose_20160122_Cisplatin_TC50_stdev, IL45CisR_highDose_20160129_Cisplatin_TC50_stdev, IL45CisR_highDose_20160206_Cisplatin_TC50_stdev)
cellline_weeks_since_revival = c(9, 10, 11, 12, 13, 15, 16, 17)
cellline_weeks_since_revival_TC50 = as.data.frame(cbind(cellline_TC50, cellline_TC50_stdev, cellline_weeks_since_revival))
colnames(cellline_weeks_since_revival_TC50) = c("TC50", "TC50_stdev", "weeks_since_revival")
cellline_weeks_since_revival_TC50$Line = "IL45CisR_highDose"
cellline_weeks_since_revival_TC50$Treatment = "Cisplatin"
Cisplatin_IL45CisR_highDose_weeks_since_revival_TC50 = cellline_weeks_since_revival_TC50

cellline_TC50 = c(IL45CisR_lowDose_20151204_Cisplatin_TC50, IL45CisR_lowDose_20151211_Cisplatin_TC50, IL45CisR_lowDose_20151218_Cisplatin_TC50, IL45CisR_lowDose_20151227_Cisplatin_TC50, IL45CisR_lowDose_20160101_Cisplatin_TC50, IL45CisR_lowDose_20160122_Cisplatin_TC50, IL45CisR_lowDose_20160129_Cisplatin_TC50, IL45CisR_lowDose_20160206_Cisplatin_TC50)
cellline_TC50_stdev = c(IL45CisR_lowDose_20151204_Cisplatin_TC50_stdev, IL45CisR_lowDose_20151211_Cisplatin_TC50_stdev, IL45CisR_lowDose_20151218_Cisplatin_TC50_stdev, IL45CisR_lowDose_20151227_Cisplatin_TC50_stdev, IL45CisR_lowDose_20160101_Cisplatin_TC50_stdev, IL45CisR_lowDose_20160122_Cisplatin_TC50_stdev, IL45CisR_lowDose_20160129_Cisplatin_TC50_stdev, IL45CisR_lowDose_20160206_Cisplatin_TC50_stdev)
cellline_weeks_since_revival = c(9, 10, 11, 12, 13, 15, 16, 17)
cellline_weeks_since_revival_TC50 = as.data.frame(cbind(cellline_TC50, cellline_TC50_stdev, cellline_weeks_since_revival))
colnames(cellline_weeks_since_revival_TC50) = c("TC50", "TC50_stdev", "weeks_since_revival")
cellline_weeks_since_revival_TC50$Line = "IL45CisR_lowDose"
cellline_weeks_since_revival_TC50$Treatment = "Cisplatin"
Cisplatin_IL45CisR_lowDose_weeks_since_revival_TC50 = cellline_weeks_since_revival_TC50

cellline_TC50 = c(IL45_high_passage_20151204_Cisplatin_TC50, IL45_high_passage_20151211_Cisplatin_TC50, IL45_high_passage_20151218_Cisplatin_TC50, IL45_high_passage_20151227_Cisplatin_TC50, IL45_high_passage_20160101_Cisplatin_TC50, IL45_high_passage_20160122_Cisplatin_TC50, IL45_high_passage_20160129_Cisplatin_TC50, IL45_high_passage_20160206_Cisplatin_TC50)
cellline_TC50_stdev = c(IL45_high_passage_20151204_Cisplatin_TC50_stdev, IL45_high_passage_20151211_Cisplatin_TC50_stdev, IL45_high_passage_20151218_Cisplatin_TC50_stdev, IL45_high_passage_20151227_Cisplatin_TC50_stdev, IL45_high_passage_20160101_Cisplatin_TC50_stdev, IL45_high_passage_20160122_Cisplatin_TC50_stdev, IL45_high_passage_20160129_Cisplatin_TC50_stdev, IL45_high_passage_20160206_Cisplatin_TC50_stdev)
cellline_weeks_since_revival = c(9, 10, 11, 12, 13, 15, 16, 17)
cellline_weeks_since_revival_TC50 = as.data.frame(cbind(cellline_TC50, cellline_TC50_stdev, cellline_weeks_since_revival))
colnames(cellline_weeks_since_revival_TC50) = c("TC50", "TC50_stdev", "weeks_since_revival")
cellline_weeks_since_revival_TC50$Line = "IL45_high_passage"
cellline_weeks_since_revival_TC50$Treatment = "Cisplatin"
Cisplatin_IL45_high_passage_weeks_since_revival_TC50 = cellline_weeks_since_revival_TC50

cellline_TC50 = c(IL45_low_passage_20151204_Cisplatin_TC50, IL45_low_passage_20151211_Cisplatin_TC50, IL45_low_passage_20151218_Cisplatin_TC50, IL45_low_passage_20151227_Cisplatin_TC50, IL45_low_passage_20160101_Cisplatin_TC50, IL45_low_passage_20160122_Cisplatin_TC50, IL45_low_passage_20160129_Cisplatin_TC50, IL45_low_passage_20160206_Cisplatin_TC50)
cellline_TC50_stdev = c(IL45_low_passage_20151204_Cisplatin_TC50_stdev, IL45_low_passage_20151211_Cisplatin_TC50_stdev, IL45_low_passage_20151218_Cisplatin_TC50_stdev, IL45_low_passage_20151227_Cisplatin_TC50_stdev, IL45_low_passage_20160101_Cisplatin_TC50_stdev, IL45_low_passage_20160122_Cisplatin_TC50_stdev, IL45_low_passage_20160129_Cisplatin_TC50_stdev, IL45_low_passage_20160206_Cisplatin_TC50_stdev)
cellline_weeks_since_revival = c(9, 10, 11, 12, 13, 15, 16, 17)
cellline_weeks_since_revival_TC50 = as.data.frame(cbind(cellline_TC50, cellline_TC50_stdev, cellline_weeks_since_revival))
colnames(cellline_weeks_since_revival_TC50) = c("TC50", "TC50_stdev", "weeks_since_revival")
cellline_weeks_since_revival_TC50$Line = "IL45_low_passage"
cellline_weeks_since_revival_TC50$Treatment = "Cisplatin"
Cisplatin_IL45_low_passage_weeks_since_revival_TC50 = cellline_weeks_since_revival_TC50

cellline_TC50 = c(IL45CisR_highDose_20151211_BAMLET_TC50, IL45CisR_highDose_20151227_BAMLET_TC50, IL45CisR_highDose_20160101_BAMLET_TC50, IL45CisR_highDose_20160108_BAMLET_TC50, IL45CisR_highDose_20160122_BAMLET_TC50, IL45CisR_highDose_20160129_BAMLET_TC50, IL45CisR_highDose_20160206_BAMLET_TC50)
cellline_TC50_stdev = c(IL45CisR_highDose_20151211_BAMLET_TC50_stdev, IL45CisR_highDose_20151227_BAMLET_TC50_stdev, IL45CisR_highDose_20160101_BAMLET_TC50_stdev, IL45CisR_highDose_20160108_BAMLET_TC50_stdev, IL45CisR_highDose_20160122_BAMLET_TC50_stdev, IL45CisR_highDose_20160129_BAMLET_TC50_stdev, IL45CisR_highDose_20160206_BAMLET_TC50_stdev)
cellline_weeks_since_revival = c(10, 12, 13, 14, 15, 16, 17)
cellline_weeks_since_revival_TC50 = as.data.frame(cbind(cellline_TC50, cellline_TC50_stdev, cellline_weeks_since_revival))
colnames(cellline_weeks_since_revival_TC50) = c("TC50", "TC50_stdev", "weeks_since_revival")
cellline_weeks_since_revival_TC50$Line = "IL45CisR_highDose"
cellline_weeks_since_revival_TC50$Treatment = "BAMLET"
BAMLET_IL45CisR_highDose_weeks_since_revival_TC50 = cellline_weeks_since_revival_TC50

cellline_TC50 = c(IL45CisR_lowDose_20151211_BAMLET_TC50, IL45CisR_lowDose_20151218_BAMLET_TC50, IL45CisR_lowDose_20151227_BAMLET_TC50, IL45CisR_lowDose_20160101_BAMLET_TC50, IL45CisR_lowDose_20160108_BAMLET_TC50, IL45CisR_lowDose_20160122_BAMLET_TC50, IL45CisR_lowDose_20160129_BAMLET_TC50, IL45CisR_lowDose_20160206_BAMLET_TC50)
cellline_TC50_stdev = c(IL45CisR_lowDose_20151211_BAMLET_TC50_stdev, IL45CisR_lowDose_20151218_BAMLET_TC50_stdev, IL45CisR_lowDose_20151227_BAMLET_TC50_stdev, IL45CisR_lowDose_20160101_BAMLET_TC50_stdev, IL45CisR_lowDose_20160108_BAMLET_TC50_stdev, IL45CisR_lowDose_20160122_BAMLET_TC50_stdev, IL45CisR_lowDose_20160129_BAMLET_TC50_stdev, IL45CisR_lowDose_20160206_BAMLET_TC50_stdev)
cellline_weeks_since_revival = c(10, 11, 12, 13, 14, 15, 16, 17)
cellline_weeks_since_revival_TC50 = as.data.frame(cbind(cellline_TC50, cellline_TC50_stdev, cellline_weeks_since_revival))
colnames(cellline_weeks_since_revival_TC50) = c("TC50", "TC50_stdev", "weeks_since_revival")
cellline_weeks_since_revival_TC50$Line = "IL45CisR_lowDose"
cellline_weeks_since_revival_TC50$Treatment = "BAMLET"
BAMLET_IL45CisR_lowDose_weeks_since_revival_TC50 = cellline_weeks_since_revival_TC50

cellline_TC50 = c(IL45_high_passage_20151211_BAMLET_TC50, IL45_high_passage_20151218_BAMLET_TC50, IL45_high_passage_20151227_BAMLET_TC50, IL45_high_passage_20160101_BAMLET_TC50, IL45_high_passage_20160108_BAMLET_TC50, IL45_high_passage_20160122_BAMLET_TC50, IL45_high_passage_20160129_BAMLET_TC50, IL45_high_passage_20160206_BAMLET_TC50)
cellline_TC50_stdev = c(IL45_high_passage_20151211_BAMLET_TC50_stdev, IL45_high_passage_20151218_BAMLET_TC50_stdev, IL45_high_passage_20151227_BAMLET_TC50_stdev, IL45_high_passage_20160101_BAMLET_TC50_stdev, IL45_high_passage_20160108_BAMLET_TC50_stdev, IL45_high_passage_20160122_BAMLET_TC50_stdev, IL45_high_passage_20160129_BAMLET_TC50_stdev, IL45_high_passage_20160206_BAMLET_TC50_stdev)
cellline_weeks_since_revival = c(10, 11, 12, 13, 14, 15, 16, 17)
cellline_weeks_since_revival_TC50 = as.data.frame(cbind(cellline_TC50, cellline_TC50_stdev, cellline_weeks_since_revival))
colnames(cellline_weeks_since_revival_TC50) = c("TC50", "TC50_stdev", "weeks_since_revival")
cellline_weeks_since_revival_TC50$Line = "IL45_high_passage"
cellline_weeks_since_revival_TC50$Treatment = "BAMLET"
BAMLET_IL45_high_passage_weeks_since_revival_TC50 = cellline_weeks_since_revival_TC50

cellline_TC50 = c(IL45_low_passage_20151211_BAMLET_TC50, IL45_low_passage_20151218_BAMLET_TC50, IL45_low_passage_20151227_BAMLET_TC50, IL45_low_passage_20160101_BAMLET_TC50, IL45_low_passage_20160108_BAMLET_TC50, IL45_low_passage_20160122_BAMLET_TC50, IL45_low_passage_20160129_BAMLET_TC50, IL45_low_passage_20160206_BAMLET_TC50)
cellline_TC50_stdev = c(IL45_low_passage_20151211_BAMLET_TC50_stdev, IL45_low_passage_20151218_BAMLET_TC50_stdev, IL45_low_passage_20151227_BAMLET_TC50_stdev, IL45_low_passage_20160101_BAMLET_TC50_stdev, IL45_low_passage_20160108_BAMLET_TC50_stdev, IL45_low_passage_20160122_BAMLET_TC50_stdev, IL45_low_passage_20160129_BAMLET_TC50_stdev, IL45_low_passage_20160206_BAMLET_TC50_stdev)
cellline_weeks_since_revival = c(10, 11, 12, 13, 14, 15, 16, 17)
cellline_weeks_since_revival_TC50 = as.data.frame(cbind(cellline_TC50, cellline_TC50_stdev, cellline_weeks_since_revival))
colnames(cellline_weeks_since_revival_TC50) = c("TC50", "TC50_stdev", "weeks_since_revival")
cellline_weeks_since_revival_TC50$Line = "IL45_low_passage"
cellline_weeks_since_revival_TC50$Treatment = "BAMLET"
BAMLET_IL45_low_passage_weeks_since_revival_TC50 = cellline_weeks_since_revival_TC50

all_TC50 = rbind(Cisplatin_IL45CisR_highDose_weeks_since_revival_TC50, Cisplatin_IL45CisR_lowDose_weeks_since_revival_TC50, Pemetrexed_IL45PemR_highDose_weeks_since_revival_TC50, Pemetrexed_IL45PemR_lowDose_weeks_since_revival_TC50, Gemcitabine_IL45GemR_highDose_weeks_since_revival_TC50, Gemcitabine_IL45GemR_lowDose_weeks_since_revival_TC50, Cisplatin_IL45_high_passage_weeks_since_revival_TC50, Cisplatin_IL45_low_passage_weeks_since_revival_TC50, Pemetrexed_IL45_high_passage_weeks_since_revival_TC50, Pemetrexed_IL45_low_passage_weeks_since_revival_TC50, Gemcitabine_IL45_high_passage_weeks_since_revival_TC50, Gemcitabine_IL45_low_passage_weeks_since_revival_TC50, BAMLET_IL45CisR_highDose_weeks_since_revival_TC50, BAMLET_IL45CisR_lowDose_weeks_since_revival_TC50, BAMLET_IL45PemR_highDose_weeks_since_revival_TC50, BAMLET_IL45PemR_lowDose_weeks_since_revival_TC50, BAMLET_IL45GemR_highDose_weeks_since_revival_TC50, BAMLET_IL45GemR_lowDose_weeks_since_revival_TC50, BAMLET_IL45VLBR_highDose_weeks_since_revival_TC50, BAMLET_IL45VLBR_lowDose_weeks_since_revival_TC50, BAMLET_IL45_high_passage_weeks_since_revival_TC50, BAMLET_IL45_low_passage_weeks_since_revival_TC50)

write.table(all_TC50, file = "cell_assays_TC50.tsv", row.names=FALSE, sep="\t", quote=FALSE)

################################################################################ Plots
# highest BAMLET dose tried was 5 mg/ml

Cisplatin_TC50 = rbind(Cisplatin_IL45CisR_highDose_weeks_since_revival_TC50, Cisplatin_IL45CisR_lowDose_weeks_since_revival_TC50, Cisplatin_IL45_high_passage_weeks_since_revival_TC50, Cisplatin_IL45_low_passage_weeks_since_revival_TC50)
Cisplatin_TC50$TC50 = ifelse((Cisplatin_TC50$TC50==100), 10, Cisplatin_TC50$TC50) # highest dose tried was 66.67 uM, none of the reached TC50 are higher than 10 uM
ggplot(Cisplatin_TC50, aes(x=weeks_since_revival, y=TC50, col=Line, size=1)) + geom_line() + scale_x_continuous(breaks=c(9,10,11,12,13,14,15,16,17)) + scale_colour_manual(values = c('red','orange','blue','cyan')) + ggtitle("Cisplatin treatment") + xlab("time since revival (weeks)") + ylab("TC50 (µM)") + theme(text=element_text(size=20)) + expand_limits(y = 0)
#dev.copy(png,'Cisplatin_TC50.png', width=800, height=600); dev.off()

BAMLET_IL45CisR_TC50 = rbind(BAMLET_IL45CisR_highDose_weeks_since_revival_TC50, BAMLET_IL45CisR_lowDose_weeks_since_revival_TC50, BAMLET_IL45_high_passage_weeks_since_revival_TC50, BAMLET_IL45_low_passage_weeks_since_revival_TC50)
ggplot(BAMLET_IL45CisR_TC50, aes(x=weeks_since_revival, y=TC50, col=Line, size=1)) + geom_line() + scale_x_continuous(breaks=c(9,10,11,12,13,14,15,16,17)) + scale_colour_manual(values = c('red','orange','blue','cyan')) + ggtitle("BAMLET treatment") + xlab("time since revival (weeks)") + ylab("TC50 (mg/ml)") + theme(text=element_text(size=20)) + expand_limits(y = 0)
#dev.copy(png,'BAMLET_IL45CisR_TC50.png', width=800, height=600); dev.off()

################################################################################ TC50 ± stdev

Cisplatin_IL45CisR_highDose_TC50 = mean(c(IL45CisR_highDose_20151204_Cisplatin_TC50, IL45CisR_highDose_20151211_Cisplatin_TC50, IL45CisR_highDose_20160101_Cisplatin_TC50, IL45CisR_highDose_20160122_Cisplatin_TC50, IL45CisR_highDose_20160129_Cisplatin_TC50, IL45CisR_highDose_20160206_Cisplatin_TC50))
Cisplatin_IL45CisR_highDose_TC50_stdev = sd(c(IL45CisR_highDose_20151204_Cisplatin_TC50, IL45CisR_highDose_20151211_Cisplatin_TC50, IL45CisR_highDose_20160101_Cisplatin_TC50, IL45CisR_highDose_20160122_Cisplatin_TC50, IL45CisR_highDose_20160129_Cisplatin_TC50, IL45CisR_highDose_20160206_Cisplatin_TC50))

Cisplatin_IL45CisR_lowDose_TC50 = mean(c(IL45CisR_lowDose_20151204_Cisplatin_TC50, IL45CisR_lowDose_20151211_Cisplatin_TC50, IL45CisR_lowDose_20151218_Cisplatin_TC50, IL45CisR_lowDose_20160101_Cisplatin_TC50, IL45CisR_lowDose_20160129_Cisplatin_TC50, IL45CisR_lowDose_20160206_Cisplatin_TC50))
Cisplatin_IL45CisR_lowDose_TC50_stdev = sd(c(IL45CisR_lowDose_20151204_Cisplatin_TC50, IL45CisR_lowDose_20151211_Cisplatin_TC50, IL45CisR_lowDose_20151218_Cisplatin_TC50, IL45CisR_lowDose_20160101_Cisplatin_TC50, IL45CisR_lowDose_20160129_Cisplatin_TC50, IL45CisR_lowDose_20160206_Cisplatin_TC50))

Cisplatin_IL45_high_passage_TC50 = mean(c(IL45_high_passage_20151204_Cisplatin_TC50, IL45_high_passage_20151211_Cisplatin_TC50, IL45_high_passage_20151218_Cisplatin_TC50, IL45_high_passage_20151227_Cisplatin_TC50, IL45_high_passage_20160101_Cisplatin_TC50, IL45_high_passage_20160122_Cisplatin_TC50, IL45_high_passage_20160129_Cisplatin_TC50, IL45_high_passage_20160206_Cisplatin_TC50))
Cisplatin_IL45_high_passage_TC50_stdev = sd(c(IL45_high_passage_20151204_Cisplatin_TC50, IL45_high_passage_20151211_Cisplatin_TC50, IL45_high_passage_20151218_Cisplatin_TC50, IL45_high_passage_20151227_Cisplatin_TC50, IL45_high_passage_20160101_Cisplatin_TC50, IL45_high_passage_20160122_Cisplatin_TC50, IL45_high_passage_20160129_Cisplatin_TC50, IL45_high_passage_20160206_Cisplatin_TC50))

Cisplatin_IL45_low_passage_TC50 = mean(c(IL45_low_passage_20151204_Cisplatin_TC50, IL45_low_passage_20151211_Cisplatin_TC50, IL45_low_passage_20151218_Cisplatin_TC50, IL45_low_passage_20151227_Cisplatin_TC50, IL45_low_passage_20160101_Cisplatin_TC50, IL45_low_passage_20160122_Cisplatin_TC50, IL45_low_passage_20160129_Cisplatin_TC50, IL45_low_passage_20160206_Cisplatin_TC50))
Cisplatin_IL45_low_passage_TC50_stdev = sd(c(IL45_low_passage_20151204_Cisplatin_TC50, IL45_low_passage_20151211_Cisplatin_TC50, IL45_low_passage_20151218_Cisplatin_TC50, IL45_low_passage_20151227_Cisplatin_TC50, IL45_low_passage_20160101_Cisplatin_TC50, IL45_low_passage_20160122_Cisplatin_TC50, IL45_low_passage_20160129_Cisplatin_TC50, IL45_low_passage_20160206_Cisplatin_TC50))

BAMLET_IL45CisR_highDose_TC50 = mean(c(IL45CisR_highDose_20151211_BAMLET_TC50, IL45CisR_highDose_20151227_BAMLET_TC50, IL45CisR_highDose_20160101_BAMLET_TC50, IL45CisR_highDose_20160108_BAMLET_TC50, IL45CisR_highDose_20160122_BAMLET_TC50, IL45CisR_highDose_20160129_BAMLET_TC50, IL45CisR_highDose_20160206_BAMLET_TC50))
BAMLET_IL45CisR_highDose_TC50_stdev = sd(c(IL45CisR_highDose_20151211_BAMLET_TC50, IL45CisR_highDose_20151227_BAMLET_TC50, IL45CisR_highDose_20160101_BAMLET_TC50, IL45CisR_highDose_20160108_BAMLET_TC50, IL45CisR_highDose_20160122_BAMLET_TC50, IL45CisR_highDose_20160129_BAMLET_TC50, IL45CisR_highDose_20160206_BAMLET_TC50))

BAMLET_IL45CisR_lowDose_TC50 = mean(c(IL45CisR_lowDose_20151211_BAMLET_TC50, IL45CisR_lowDose_20151218_BAMLET_TC50, IL45CisR_lowDose_20151227_BAMLET_TC50, IL45CisR_lowDose_20160101_BAMLET_TC50, IL45CisR_lowDose_20160108_BAMLET_TC50, IL45CisR_lowDose_20160122_BAMLET_TC50, IL45CisR_lowDose_20160129_BAMLET_TC50, IL45CisR_lowDose_20160206_BAMLET_TC50))
BAMLET_IL45CisR_lowDose_TC50_stdev = sd(c(IL45CisR_lowDose_20151211_BAMLET_TC50, IL45CisR_lowDose_20151218_BAMLET_TC50, IL45CisR_lowDose_20151227_BAMLET_TC50, IL45CisR_lowDose_20160101_BAMLET_TC50, IL45CisR_lowDose_20160108_BAMLET_TC50, IL45CisR_lowDose_20160122_BAMLET_TC50, IL45CisR_lowDose_20160129_BAMLET_TC50, IL45CisR_lowDose_20160206_BAMLET_TC50))

BAMLET_IL45_high_passage_TC50 = mean(c(IL45_high_passage_20151211_BAMLET_TC50, IL45_high_passage_20151218_BAMLET_TC50, IL45_high_passage_20151227_BAMLET_TC50, IL45_high_passage_20160101_BAMLET_TC50, IL45_high_passage_20160108_BAMLET_TC50, IL45_high_passage_20160122_BAMLET_TC50, IL45_high_passage_20160129_BAMLET_TC50, IL45_high_passage_20160206_BAMLET_TC50))
BAMLET_IL45_high_passage_TC50_stdev = sd(c(IL45_high_passage_20151211_BAMLET_TC50, IL45_high_passage_20151218_BAMLET_TC50, IL45_high_passage_20151227_BAMLET_TC50, IL45_high_passage_20160101_BAMLET_TC50, IL45_high_passage_20160108_BAMLET_TC50, IL45_high_passage_20160122_BAMLET_TC50, IL45_high_passage_20160129_BAMLET_TC50, IL45_high_passage_20160206_BAMLET_TC50))

BAMLET_IL45_low_passage_TC50 = mean(c(IL45_low_passage_20151211_BAMLET_TC50, IL45_low_passage_20151218_BAMLET_TC50, IL45_low_passage_20151227_BAMLET_TC50, IL45_low_passage_20160101_BAMLET_TC50, IL45_low_passage_20160108_BAMLET_TC50, IL45_low_passage_20160122_BAMLET_TC50, IL45_low_passage_20160129_BAMLET_TC50, IL45_low_passage_20160206_BAMLET_TC50))
BAMLET_IL45_low_passage_TC50_stdev = sd(c(IL45_low_passage_20151211_BAMLET_TC50, IL45_low_passage_20151218_BAMLET_TC50, IL45_low_passage_20151227_BAMLET_TC50, IL45_low_passage_20160101_BAMLET_TC50, IL45_low_passage_20160108_BAMLET_TC50, IL45_low_passage_20160122_BAMLET_TC50, IL45_low_passage_20160129_BAMLET_TC50, IL45_low_passage_20160206_BAMLET_TC50))



Cisplatin_IL45CisR_highDose_TC50; Cisplatin_IL45CisR_highDose_TC50_stdev
Cisplatin_IL45CisR_lowDose_TC50; Cisplatin_IL45CisR_lowDose_TC50_stdev
Cisplatin_IL45_high_passage_TC50; Cisplatin_IL45_high_passage_TC50_stdev
Cisplatin_IL45_low_passage_TC50; Cisplatin_IL45_low_passage_TC50_stdev
BAMLET_IL45CisR_highDose_TC50; BAMLET_IL45CisR_highDose_TC50_stdev
BAMLET_IL45CisR_lowDose_TC50; BAMLET_IL45CisR_lowDose_TC50_stdev
BAMLET_IL45_high_passage_TC50; BAMLET_IL45_high_passage_TC50_stdev
BAMLET_IL45_low_passage_TC50; BAMLET_IL45_low_passage_TC50_stdev

	# > Cisplatin_IL45CisR_highDose_TC50; Cisplatin_IL45CisR_highDose_TC50_stdev
	# [1] 3.178676
	# [1] 1.896892
	# > Cisplatin_IL45CisR_lowDose_TC50; Cisplatin_IL45CisR_lowDose_TC50_stdev
	# [1] 4.222237
	# [1] 2.834081
	# > Cisplatin_IL45_high_passage_TC50; Cisplatin_IL45_high_passage_TC50_stdev
	# [1] 1.73854
	# [1] 0.3079673
	# > Cisplatin_IL45_low_passage_TC50; Cisplatin_IL45_low_passage_TC50_stdev
	# [1] 1.90977
	# [1] 0.837261
	# > BAMLET_IL45CisR_highDose_TC50; BAMLET_IL45CisR_highDose_TC50_stdev
	# [1] 0.8734641
	# [1] 0.2089909
	# > BAMLET_IL45CisR_lowDose_TC50; BAMLET_IL45CisR_lowDose_TC50_stdev
	# [1] 1.062096
	# [1] 0.3025256
	# > BAMLET_IL45_high_passage_TC50; BAMLET_IL45_high_passage_TC50_stdev
	# [1] 0.8765398
	# [1] 0.1777053
	# > BAMLET_IL45_low_passage_TC50; BAMLET_IL45_low_passage_TC50_stdev
	# [1] 0.873435
	# [1] 0.1820021

c(IL45CisR_highDose_20151204_Cisplatin_TC50, IL45CisR_highDose_20151211_Cisplatin_TC50, IL45CisR_highDose_20151218_Cisplatin_TC50, IL45CisR_highDose_20151227_Cisplatin_TC50, IL45CisR_highDose_20160101_Cisplatin_TC50, IL45CisR_highDose_20160122_Cisplatin_TC50, IL45CisR_highDose_20160129_Cisplatin_TC50, IL45CisR_highDose_20160206_Cisplatin_TC50)
c(IL45CisR_lowDose_20151204_Cisplatin_TC50, IL45CisR_lowDose_20151211_Cisplatin_TC50, IL45CisR_lowDose_20151218_Cisplatin_TC50, IL45CisR_lowDose_20151227_Cisplatin_TC50, IL45CisR_lowDose_20160101_Cisplatin_TC50, IL45CisR_lowDose_20160122_Cisplatin_TC50, IL45CisR_lowDose_20160129_Cisplatin_TC50, IL45CisR_lowDose_20160206_Cisplatin_TC50)
c(IL45_low_passage_20151204_Cisplatin_TC50, IL45_low_passage_20151211_Cisplatin_TC50, IL45_low_passage_20151218_Cisplatin_TC50, IL45_low_passage_20151227_Cisplatin_TC50, IL45_low_passage_20160101_Cisplatin_TC50, IL45_low_passage_20160122_Cisplatin_TC50, IL45_low_passage_20160129_Cisplatin_TC50, IL45_low_passage_20160206_Cisplatin_TC50)
c(IL45_high_passage_20151204_Pemetrexed_TC50, IL45_high_passage_20151211_Pemetrexed_TC50, IL45_high_passage_20151227_Pemetrexed_TC50, IL45_high_passage_20160101_Pemetrexed_TC50, IL45_high_passage_20160108_Pemetrexed_TC50, IL45_high_passage_20160122_Pemetrexed_TC50, IL45_high_passage_20160206_Pemetrexed_TC50)
c(IL45_low_passage_20151211_Pemetrexed_TC50, IL45_low_passage_20160101_Pemetrexed_TC50, IL45_low_passage_20160108_Pemetrexed_TC50, IL45_low_passage_20160122_Pemetrexed_TC50, IL45_low_passage_20160206_Pemetrexed_TC50)
c(IL45_high_passage_20151204_Gemcitabine_TC50, IL45_high_passage_20151211_Gemcitabine_TC50, IL45_high_passage_20151218_Gemcitabine_TC50, IL45_high_passage_20151227_Gemcitabine_TC50, IL45_high_passage_20160101_Gemcitabine_TC50, IL45_high_passage_20160108_Gemcitabine_TC50, IL45_high_passage_20160122_Gemcitabine_TC50, IL45_high_passage_20160129_Gemcitabine_TC50, IL45_high_passage_20160206_Gemcitabine_TC50)
c(IL45_low_passage_20151204_Gemcitabine_TC50, IL45_low_passage_20151211_Gemcitabine_TC50, IL45_low_passage_20151218_Gemcitabine_TC50, IL45_low_passage_20151227_Gemcitabine_TC50, IL45_low_passage_20160101_Gemcitabine_TC50, IL45_low_passage_20160108_Gemcitabine_TC50, IL45_low_passage_20160122_Gemcitabine_TC50, IL45_low_passage_20160129_Gemcitabine_TC50, IL45_low_passage_20160206_Gemcitabine_TC50)
c(IL45CisR_highDose_20151211_BAMLET_TC50, IL45CisR_highDose_20151227_BAMLET_TC50, IL45CisR_highDose_20160101_BAMLET_TC50, IL45CisR_highDose_20160108_BAMLET_TC50, IL45CisR_highDose_20160122_BAMLET_TC50, IL45CisR_highDose_20160129_BAMLET_TC50, IL45CisR_highDose_20160206_BAMLET_TC50)
c(IL45CisR_lowDose_20151211_BAMLET_TC50, IL45CisR_lowDose_20151218_BAMLET_TC50, IL45CisR_lowDose_20151227_BAMLET_TC50, IL45CisR_lowDose_20160101_BAMLET_TC50, IL45CisR_lowDose_20160108_BAMLET_TC50, IL45CisR_lowDose_20160122_BAMLET_TC50, IL45CisR_lowDose_20160129_BAMLET_TC50, IL45CisR_lowDose_20160206_BAMLET_TC50)
c(IL45_high_passage_20151211_BAMLET_TC50, IL45_high_passage_20151218_BAMLET_TC50, IL45_high_passage_20151227_BAMLET_TC50, IL45_high_passage_20160101_BAMLET_TC50, IL45_high_passage_20160108_BAMLET_TC50, IL45_high_passage_20160122_BAMLET_TC50, IL45_high_passage_20160129_BAMLET_TC50, IL45_high_passage_20160206_BAMLET_TC50)
c(IL45_low_passage_20151211_BAMLET_TC50, IL45_low_passage_20151218_BAMLET_TC50, IL45_low_passage_20151227_BAMLET_TC50, IL45_low_passage_20160101_BAMLET_TC50, IL45_low_passage_20160108_BAMLET_TC50, IL45_low_passage_20160122_BAMLET_TC50, IL45_low_passage_20160129_BAMLET_TC50, IL45_low_passage_20160206_BAMLET_TC50)

	# > c(IL45CisR_highDose_20151204_Cisplatin_TC50, IL45CisR_highDose_20151211_Cisplatin_TC50, IL45CisR_highDose_20151218_Cisplatin_TC50, IL45CisR_highDose_20151227_Cisplatin_TC50, IL45CisR_highDose_20160101_Cisplatin_TC50, IL45CisR_highDose_20160122_Cisplatin_TC50, IL45CisR_highDose_20160129_Cisplatin_TC50, IL45CisR_highDose_20160206_Cisplatin_TC50)
	# [1]   6.983844   2.905043 100.000000 100.000000   2.412902   2.697545   2.104016
	# [8]   1.968707
	# > c(IL45CisR_lowDose_20151204_Cisplatin_TC50, IL45CisR_lowDose_20151211_Cisplatin_TC50, IL45CisR_lowDose_20151218_Cisplatin_TC50, IL45CisR_lowDose_20151227_Cisplatin_TC50, IL45CisR_lowDose_20160101_Cisplatin_TC50, IL45CisR_lowDose_20160122_Cisplatin_TC50, IL45CisR_lowDose_20160129_Cisplatin_TC50, IL45CisR_lowDose_20160206_Cisplatin_TC50)
	# [1]   8.5117833   3.9682222   3.3946562 100.0000000   2.5990323 100.0000000
	# [7]   0.5041196   6.3556055
	# > c(IL45_low_passage_20151204_Cisplatin_TC50, IL45_low_passage_20151211_Cisplatin_TC50, IL45_low_passage_20151218_Cisplatin_TC50, IL45_low_passage_20151227_Cisplatin_TC50, IL45_low_passage_20160101_Cisplatin_TC50, IL45_low_passage_20160122_Cisplatin_TC50, IL45_low_passage_20160129_Cisplatin_TC50, IL45_low_passage_20160206_Cisplatin_TC50)
	# [1] 2.2248707 1.8916575 1.0050285 0.6226676 2.0247632 3.3486285 2.3889451
	# [8] 1.7716026
	# > c(IL45_high_passage_20151204_Pemetrexed_TC50, IL45_high_passage_20151211_Pemetrexed_TC50, IL45_high_passage_20151227_Pemetrexed_TC50, IL45_high_passage_20160101_Pemetrexed_TC50, IL45_high_passage_20160108_Pemetrexed_TC50, IL45_high_passage_20160122_Pemetrexed_TC50, IL45_high_passage_20160206_Pemetrexed_TC50)
	# [1]     3.399207     2.272662     9.557566    19.639390 10000.000000
	# [6] 10000.000000     5.150058
	# > c(IL45_low_passage_20151211_Pemetrexed_TC50, IL45_low_passage_20160101_Pemetrexed_TC50, IL45_low_passage_20160108_Pemetrexed_TC50, IL45_low_passage_20160122_Pemetrexed_TC50, IL45_low_passage_20160206_Pemetrexed_TC50)
	# [1]     2.164996     2.735298 10000.000000 10000.000000    10.058034
	# > c(IL45_high_passage_20151204_Gemcitabine_TC50, IL45_high_passage_20151211_Gemcitabine_TC50, IL45_high_passage_20151218_Gemcitabine_TC50, IL45_high_passage_20151227_Gemcitabine_TC50, IL45_high_passage_20160101_Gemcitabine_TC50, IL45_high_passage_20160108_Gemcitabine_TC50, IL45_high_passage_20160122_Gemcitabine_TC50, IL45_high_passage_20160129_Gemcitabine_TC50, IL45_high_passage_20160206_Gemcitabine_TC50)
	# [1] 0.3106706 0.2305606 1.2622178 0.2671923 0.2386947 0.1825922 0.4362875
	# [8] 0.4849251 0.7429705
	# > c(IL45_low_passage_20151204_Gemcitabine_TC50, IL45_low_passage_20151211_Gemcitabine_TC50, IL45_low_passage_20151218_Gemcitabine_TC50, IL45_low_passage_20151227_Gemcitabine_TC50, IL45_low_passage_20160101_Gemcitabine_TC50, IL45_low_passage_20160108_Gemcitabine_TC50, IL45_low_passage_20160122_Gemcitabine_TC50, IL45_low_passage_20160129_Gemcitabine_TC50, IL45_low_passage_20160206_Gemcitabine_TC50)
	# [1] 0.5150911 0.3151779 0.1000963 0.2483025 0.4104127 0.2316589 0.4622327
	# [8] 0.5433966 0.3752460
	# > c(IL45CisR_highDose_20151211_BAMLET_TC50, IL45CisR_highDose_20151227_BAMLET_TC50, IL45CisR_highDose_20160101_BAMLET_TC50, IL45CisR_highDose_20160108_BAMLET_TC50, IL45CisR_highDose_20160122_BAMLET_TC50, IL45CisR_highDose_20160129_BAMLET_TC50, IL45CisR_highDose_20160206_BAMLET_TC50)
	# [1] 0.7479777 1.2596414 0.8138849 0.6769419 0.7117944 0.8610682 1.0429402
	# > c(IL45CisR_lowDose_20151211_BAMLET_TC50, IL45CisR_lowDose_20151218_BAMLET_TC50, IL45CisR_lowDose_20151227_BAMLET_TC50, IL45CisR_lowDose_20160101_BAMLET_TC50, IL45CisR_lowDose_20160108_BAMLET_TC50, IL45CisR_lowDose_20160122_BAMLET_TC50, IL45CisR_lowDose_20160129_BAMLET_TC50, IL45CisR_lowDose_20160206_BAMLET_TC50)
	# [1] 1.3093751 1.2756277 0.9792075 0.5755538 0.7493471 1.4755103 0.9622455
	# [8] 1.1699011
	# > c(IL45_high_passage_20151211_BAMLET_TC50, IL45_high_passage_20151218_BAMLET_TC50, IL45_high_passage_20151227_BAMLET_TC50, IL45_high_passage_20160101_BAMLET_TC50, IL45_high_passage_20160108_BAMLET_TC50, IL45_high_passage_20160122_BAMLET_TC50, IL45_high_passage_20160129_BAMLET_TC50, IL45_high_passage_20160206_BAMLET_TC50)
	# [1] 0.8965339 0.8749830 1.2874866 0.8130883 0.8808707 0.8180858 0.7442797
	# [8] 0.7807467
	# > c(IL45_low_passage_20151211_BAMLET_TC50, IL45_low_passage_20151218_BAMLET_TC50, IL45_low_passage_20151227_BAMLET_TC50, IL45_low_passage_20160101_BAMLET_TC50, IL45_low_passage_20160108_BAMLET_TC50, IL45_low_passage_20160122_BAMLET_TC50, IL45_low_passage_20160129_BAMLET_TC50, IL45_low_passage_20160206_BAMLET_TC50)
	# [1] 0.6358328 1.1359619 1.0720092 0.7293322 0.7515057 0.8259691 1.0305305
	# [8] 0.8063385


