# Rscript step03_determine_TC50_and_stdev.R <input_file> <cell_line> <model> > <output_file>
# Rscript step03_determine_TC50_and_stdev.R human_bamlet_cell_assay_data.txt HOFM01 2component > step01_output.txt
# Rscript step03_determine_TC50_and_stdev.R human_bamlet_cell_assay_data.txt MM05 1component > step01_output.txt

library(reshape2)
library(ggplot2)
library(optimx)
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'human_bamlet_cell_assay_data.txt', 'HOFM01' ) # for testing

args_infile = args[1] # human_bamlet_cell_assay_data.txt
args_cellline = args[2] # HOFM01
args_model = args[3] # 1component or 2component

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
# We have determined that the appropriate model for HOFM01, HOFM02, and HOFF is a biphasic 2-component model,
# and the appropriate model for the other cell lines (all human mesothelioma cell lines, Met5A, and Humofib) is the 1-component model.
# We have determined that HOFM01, HOFM02, and HOFF have BAMLET-resistant cell fractions whose BAMLET-treatment TC50s
# are statistically significantly higher (more than double) than the mesothelioma TC50s.
# The statistical significances were determined by comparing pairs of models and the models produced TC50s for each pair.
# Now we will determine the TC50 and its standard deviation for BAMLET treatment of these cell lines using these models
# for each cell line individually.
# We will determine the TC50 for each experiment and calculate the stdev of TC50 from the multiple experiments.

####################################################################################################
# Define the various possible models and their functions.
# The main 2 types of models are for biphasic (2component; two LC50, defined as inflection point and shift of second inflexion point with respect to the first)
# and uniphasic (1component; only one inflextion point).
# 
# A = left asymptote
# B = right asymptote
# xmid = inflection point (TC50)
# scal = factor affecting slope at inflection point, Hill's slope, steepness of curve
####################################################################################################

##########
predict.1component = function(params, x)
{
	A = as.numeric(params[1])
	B = as.numeric(params[2])
	scal = as.numeric(params[3])
	xmid = as.numeric(params[4])
	prediction_result = A + (B - A) * (1/(1 + exp((xmid - x) / scal)))
	prediction_result
}
loglik.1component = function(params, x, y)
{
	sigma = as.numeric(params[5])
	yhat = predict.1component(params, x)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.1component = function(params)
{
	names(params) <- NULL 
	xmid = as.numeric(params[4])
	c(TC50 = exp(xmid))
}

##############################
# A = left asymptote
# B = right asymptote
# xmid = inflection point 1 (TC50 of BAMLET-sensitive cells)
# xmid + xshift = inflection point 2 (TC50 of BAMLET-resistant cells)
# scal = factor affecting slope at inflection points, Hill's slope, steepness of curves
# sens = fraction of cells having inflection point 1 (fraction of BAMLET-sensitive cells)
# 1 - sens = fraction of cells having inflection point 2 (fraction of BAMLET-resistant cells)
# sens can vary per experiment
##############################
predict.2component = function(params, x)
{
	A = as.numeric(params[1])
	B = as.numeric(params[2])
	scal = as.numeric(params[3])
	xmid = as.numeric(params[4])
	sens = as.numeric(params[5])
	xshift = as.numeric(params[6])
	finalSens = sens
	prediction_result = A + (B - A) * (finalSens / (1 + exp((xmid - x) / scal)) + (1 - finalSens) / (1 + exp((xmid + xshift - x) / scal)))
	prediction_result
}
loglik.2component = function(params, x, y)
{
	sigma = as.numeric(params[7])
	yhat = predict.2component(params, x)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.2component = function(params)
{
	names(params) <- NULL 
	xmid = as.numeric(params[4])
	sens = as.numeric(params[5])
	xshift = as.numeric(params[6])
	sensitive_fraction_average = sens
	resistant_fraction_average = 1 - sensitive_fraction_average
	c(sensitive_fraction_TC50 = exp(xmid), resistant_fraction_TC50 = exp(xmid + xshift), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average)
}

##############################
assign_start_lower_upper_values = function()
{
	A_start = 1;            A_low = -Inf;            A_high = Inf;
	B_start = 0;            B_low = -Inf;            B_high = Inf;
	scal_start = 0.1;       scal_low = 0.02;         scal_high = Inf;
	xmid_start = 0;         xmid_low = -Inf;         xmid_high = Inf;
	sens_start = 0.5;       sens_low = 0;            sens_high = 1;
	xshift_start = 1;       xshift_low = 0;          xshift_high = Inf;
	sigma_start = 0.1;      sigma_low = 1e-4;        sigma_high = Inf;

	start.1component  = c(A=A_start,   B=B_start,   scal=scal_start,   xmid=xmid_start,   sigma=sigma_start)
	start2.1component = c(A=A_start,   B=0.5,       scal=scal_start,   xmid=xmid_start,   sigma=sigma_start)
	lower.1component  = c(A=A_low,     B=B_low,     scal=scal_low,     xmid=xmid_low,     sigma=sigma_low)
	lower2.1component = c(A=A_low,     B=0,         scal=scal_low,     xmid=xmid_low,     sigma=sigma_low)
	upper.1component  = c(A=A_high,    B=B_high,    scal=scal_high,    xmid=xmid_high,    sigma=sigma_high)
	upper2.1component = c(A=A_high,    B=1,         scal=scal_high,    xmid=xmid_high,    sigma=sigma_high)
	assign("start.1component", start.1component, envir = .GlobalEnv)
	assign("lower.1component", lower.1component, envir = .GlobalEnv)
	assign("upper.1component", upper.1component, envir = .GlobalEnv)
	assign("start2.1component", start2.1component, envir = .GlobalEnv)
	assign("lower2.1component", lower2.1component, envir = .GlobalEnv)
	assign("upper2.1component", upper2.1component, envir = .GlobalEnv)

	start.2component = c(A=A_start,   B=B_start,   scal=scal_start,   xmid=xmid_start,   sens=sens_start,   xshift=xshift_start,   sigma=sigma_start)
	lower.2component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sens=sens_low,   xshift=xshift_low,   sigma=sigma_low)
	upper.2component = c(A=A_high,   B=B_high,   scal=scal_high,   xmid=xmid_high,   sens=sens_high,   xshift=xshift_high,   sigma=sigma_high)
	assign("start.2component", start.2component, envir = .GlobalEnv)
	assign("lower.2component", lower.2component, envir = .GlobalEnv)
	assign("upper.2component", upper.2component, envir = .GlobalEnv)
}
assign_start_lower_upper_values()

####################################################################################################
##### Find the TC50 and stdev for input cell line, using the input model
####################################################################################################

cellline_0 = assays[ (assays$Line == args_cellline), ]
exps = unique(cellline_0$Experiment)

if (args_model == '2component') {

	cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[1])), ]
	fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
	TC50_sensitive = exp(fit.2component$xmid)
	TC50_resistant = exp(fit.2component$xmid + fit.2component$xshift)

	list_TC50_sensitive = c(TC50_sensitive)
	list_TC50_resistant = c(TC50_resistant)

	if (length(exps) >= 2) {
		cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[2])), ]
		fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
		TC50_sensitive = exp(fit.2component$xmid)
		TC50_resistant = exp(fit.2component$xmid + fit.2component$xshift)

		list_TC50_sensitive = append(list_TC50_sensitive, c(TC50_sensitive))
		list_TC50_resistant = append(list_TC50_resistant, c(TC50_resistant))
	}

	if (length(exps) >= 3) {
		cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[3])), ]
		fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
		TC50_sensitive = exp(fit.2component$xmid)
		TC50_resistant = exp(fit.2component$xmid + fit.2component$xshift)

		list_TC50_sensitive = append(list_TC50_sensitive, c(TC50_sensitive))
		list_TC50_resistant = append(list_TC50_resistant, c(TC50_resistant))
	}

	if (length(exps) >= 4) {
		cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[4])), ]
		fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
		TC50_sensitive = exp(fit.2component$xmid)
		TC50_resistant = exp(fit.2component$xmid + fit.2component$xshift)

		list_TC50_sensitive = append(list_TC50_sensitive, c(TC50_sensitive))
		list_TC50_resistant = append(list_TC50_resistant, c(TC50_resistant))
	}

	if (length(exps) >= 5) {
		cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[5])), ]
		fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
		TC50_sensitive = exp(fit.2component$xmid)
		TC50_resistant = exp(fit.2component$xmid + fit.2component$xshift)

		list_TC50_sensitive = append(list_TC50_sensitive, c(TC50_sensitive))
		list_TC50_resistant = append(list_TC50_resistant, c(TC50_resistant))
	}

	TC50_sensitive_mean = mean(list_TC50_sensitive)
	TC50_sensitive_sd = sd(list_TC50_sensitive)
	TC50_resistant_mean = mean(list_TC50_resistant)
	TC50_resistant_sd = sd(list_TC50_resistant)
	TC50_sensitive
	TC50_sensitive_sd
	TC50_resistant
	TC50_resistant_sd
}

if (args_model == '1component') {

	cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[1])), ]
	fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
	TC50_sensitive = exp(fit.1component$xmid)

	list_TC50_sensitive = c(TC50_sensitive)

	if (length(exps) >= 2) {
		cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[2])), ]
		fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
		TC50_sensitive = exp(fit.1component$xmid)
	}

	if (length(exps) >= 3) {
		cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[3])), ]
		fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
		TC50_sensitive = exp(fit.1component$xmid)
	}

	if (length(exps) >= 4) {
		cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[4])), ]
		fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
		TC50_sensitive = exp(fit.1component$xmid)
}

	if (length(exps) >= 5) {
		cellline = assays[ ((assays$Line == args_cellline) & (assays$Experiment == exps[5])), ]
		fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
		TC50_sensitive = exp(fit.1component$xmid)
	}

	TC50_sensitive_mean = mean(list_TC50_sensitive)
	TC50_sensitive_sd = sd(list_TC50_sensitive)
	TC50_sensitive
	TC50_sensitive_sd
}

# Example output
# 
# cellline = assays[ ((assays$Line == "HOFM01") & (assays$Experiment == "20151227")), ]
# fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# TC50_sensitive_1 = exp(fit.2component$xmid)
# TC50_resistant_1 = exp(fit.2component$xmid + fit.2component$xshift)
# 
# cellline = assays[ ((assays$Line == "HOFM01") & (assays$Experiment == "20160101")), ]
# fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# TC50_sensitive_2 = exp(fit.2component$xmid)
# TC50_resistant_2 = exp(fit.2component$xmid + fit.2component$xshift)
# 
# cellline = assays[ ((assays$Line == "HOFM01") & (assays$Experiment == "20160108")), ]
# fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# TC50_sensitive_3 = exp(fit.2component$xmid)
# TC50_resistant_3 = exp(fit.2component$xmid + fit.2component$xshift)
# 
# TC50_sensitive = mean(c(TC50_sensitive_1,TC50_sensitive_2,TC50_sensitive_3))
# TC50_sensitive_sd = sd(c(TC50_sensitive_1,TC50_sensitive_2,TC50_sensitive_3))
# TC50_resistant = mean(c(TC50_resistant_1,TC50_resistant_2,TC50_resistant_3))
# TC50_resistant_sd = sd(c(TC50_resistant_1,TC50_resistant_2,TC50_resistant_3))
# TC50_sensitive
# TC50_sensitive_sd
# TC50_resistant
# TC50_resistant_sd
# 
# 	> TC50_sensitive
# 	[1] 0.7222377
# 	> TC50_sensitive_sd
# 	[1] 0.2288617
# 	> TC50_resistant
# 	[1] 1.786305
# 	> TC50_resistant_sd
# 	[1] 0.7121097

# Example output
# 
# cellline = assays[ ((assays$Line == "MM05") & (assays$Experiment == "20151227")), ]
# fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# TC50_sensitive_1 = exp(fit.1component$xmid)
# 
# cellline = assays[ ((assays$Line == "MM05") & (assays$Experiment == "20160101")), ]
# fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# TC50_sensitive_2 = exp(fit.1component$xmid)
# 
# cellline = assays[ ((assays$Line == "MM05") & (assays$Experiment == "20160108")), ]
# fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# TC50_sensitive_3 = exp(fit.1component$xmid)
# 
# cellline = assays[ ((assays$Line == "MM05") & (assays$Experiment == "20160206")), ]
# fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# TC50_sensitive_4 = exp(fit.1component$xmid)
# 
# TC50_sensitive = mean(c(TC50_sensitive_1,TC50_sensitive_2,TC50_sensitive_3,TC50_sensitive_4))
# TC50_sensitive_sd = sd(c(TC50_sensitive_1,TC50_sensitive_2,TC50_sensitive_3,TC50_sensitive_4))
# TC50_sensitive
# TC50_sensitive_sd
# 
# 	> TC50_sensitive
# 	[1] 0.6810862
# 	> TC50_sensitive_sd
# 	[1] 0.1001573


