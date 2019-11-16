# Rscript step01_choose_model.R <input_file> <cell_line> > <output_file>
# Rscript step01_choose_model.R human_bamlet_cell_assay_data.txt HOFM01 > step01_output.txt

library(reshape2)
library(ggplot2)
library(optimx)
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'human_bamlet_cell_assay_data.txt', 'HOFM01' ) # for testing

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
args_cellline = args[2] # HOFM01

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

# A first look at data of all cell lines.
# ggplot(assays, aes(x = log10(Concentration), y = Viability, col = Experiment)) + geom_point() + facet_wrap(~ Line)

####################################################################################################
# Define the various possible models and their functions.
# The main 2 types of models are for biphasic 
# (two LC50, defined as inflection point and shift of second inflexion point with respect to the first)
# and uniphasic (only one inflextion point).
# In case adjustments per experiment need to be made for modelling, define models permitting that.
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

# A = left asymptote
# B = right asymptote
# xmid = inflection point 1 (TC50 of BAMLET-sensitive cells)
# xmid + xshift = inflection point 2 (TC50 of BAMLET-resistant cells)
# scal = factor affecting slope at inflection points
# sens = fraction of cells having inflection point 1 (fraction of BAMLET-sensitive cells)
# 1 - sens = fraction of cells having inflection point 2 (fraction of BAMLET-resistant cells)
predict.2component = function(params, x)
{
	A = params[1]
	B = params[2]
	scal = params[3]
	xmid = params[4]
	sens = params[5]
	xshift = params[6]
	finalSens = sens
	prediction_result = A + (B - A) * (finalSens / (1 + exp((xmid - x) / scal)) + (1 - finalSens) / (1 + exp((xmid + xshift - x) / scal)))
	prediction_result
}
loglik.2component = function(params, x, y)
{
	sigma = params[7]
	yhat = predict.2component(params, x)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.2component = function(params)
{
	names(params) <- NULL 
	xmid = params[4]
	sens = params[5]
	xshift = params[6]
	sensitive_fraction_average = sens
	resistant_fraction_average = 1 - sensitive_fraction_average
	c(sensitive_fraction_TC50 = exp(xmid), resistant_fraction_TC50 = exp(xmid + xshift), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average)
}

predict.2component.2exp = function(params, x, is_Exp1, is_Exp2)
{
	A = params[1]
	B = params[2]
	scal = params[3]
	xmid = params[4]
	sens_Exp1 = params[5]
	sens_Exp2 = params[6]
	xshift = params[7]
	finalSens = (sens_Exp1 * is_Exp1) + (sens_Exp2 * is_Exp2)
	prediction_result = A + (B - A) * (finalSens / (1 + exp((xmid - x) / scal)) + (1 - finalSens) / (1 + exp((xmid + xshift - x) / scal)))
	prediction_result
}
loglik.2component.2exp = function(params, x, y, is_Exp1, is_Exp2)
{
	sigma = params[8]
	yhat = predict.2component.2exp(params, x, is_Exp1, is_Exp2)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.2component.2exp = function(params)
{
	#names(params) <- NULL 
	xmid = params[4]
	sens_Exp1 = as.numeric(params[5])
	sens_Exp2 = as.numeric(params[6])
	xshift = params[7]
	sensitive_fraction_average = mean(c(sens_Exp1,sens_Exp2))
	resistant_fraction_average = 1 - sensitive_fraction_average
	c(sensitive_fraction_TC50 = exp(xmid), resistant_fraction_TC50 = exp(xmid + xshift), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average)
}

predict.2component.3exp = function(params, x, is_Exp1, is_Exp2, is_Exp3)
{
	A = params[1]
	B = params[2]
	scal = params[3]
	xmid = params[4]
	sens_Exp1 = params[5]
	sens_Exp2 = params[6]
	sens_Exp3 = params[7]
	xshift = params[8]
	finalSens = (sens_Exp1 * is_Exp1) + (sens_Exp2 * is_Exp2) + (sens_Exp3 * is_Exp3)
	prediction_result = A + (B - A) * (finalSens / (1 + exp((xmid - x) / scal)) + (1 - finalSens) / (1 + exp((xmid + xshift - x) / scal)))
	prediction_result
}
loglik.2component.3exp = function(params, x, y, is_Exp1, is_Exp2, is_Exp3)
{
	sigma = params[9]
	yhat = predict.2component.3exp(params, x, is_Exp1, is_Exp2, is_Exp3)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.2component.3exp = function(params)
{
	#names(params) <- NULL 
	xmid = params[4]
	sens_Exp1 = as.numeric(params[5])
	sens_Exp2 = as.numeric(params[6])
	sens_Exp3 = as.numeric(params[7])
	xshift = params[8]
	sensitive_fraction_average = mean(c(sens_Exp1,sens_Exp2,sens_Exp3))
	resistant_fraction_average = 1 - sensitive_fraction_average
	c(sensitive_fraction_TC50 = exp(xmid), resistant_fraction_TC50 = exp(xmid + xshift), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average)
}

predict.2component.4exp = function(params, x, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	A = params[1]
	B = params[2]
	scal = params[3]
	xmid = params[4]
	sens_Exp1 = params[5]
	sens_Exp2 = params[6]
	sens_Exp3 = params[7]
	sens_Exp4 = params[8]
	xshift = params[9]
	finalSens = (sens_Exp1 * is_Exp1) + (sens_Exp2 * is_Exp2) + (sens_Exp3 * is_Exp3) + (sens_Exp4 * is_Exp4)
	prediction_result = A + (B - A) * (finalSens / (1 + exp((xmid - x) / scal)) + (1 - finalSens) / (1 + exp((xmid + xshift - x) / scal)))
	prediction_result
}
loglik.2component.4exp = function(params, x, y, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
{
	sigma = params[10]
	yhat = predict.2component.4exp(params, x, is_Exp1, is_Exp2, is_Exp3, is_Exp4)
	residuals = y - yhat
	loglik_result = sum(dnorm(residuals, sd = sigma, log = TRUE))
	loglik_result
}
convert_values_back_to_interpretable.2component.4exp = function(params)
{
	names(params) <- NULL 
	xmid = params[4]
	sens_Exp1 = params[5]
	sens_Exp2 = params[6]
	sens_Exp3 = params[7]
	sens_Exp4 = params[8]
	xshift = params[9]
	sensitive_fraction_average = mean(c(sens_Exp1,sens_Exp2,sens_Exp3,sens_Exp4))
	resistant_fraction_average = 1 - sensitive_fraction_average
	c(sensitive_fraction_TC50 = exp(xmid), resistant_fraction_TC50 = exp(xmid + xshift), sensitive_fraction_average = sensitive_fraction_average, resistant_fraction_average = resistant_fraction_average)
}

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
	upper.1component  = c(A=A_high,    B=B_high,    scal=scal_high,    xmid=xmid_high,    sigma=sigma_high)
	assign("start.1component", start.1component, envir = .GlobalEnv)
	assign("lower.1component", lower.1component, envir = .GlobalEnv)
	assign("upper.1component", upper.1component, envir = .GlobalEnv)
	assign("start2.1component", start2.1component, envir = .GlobalEnv)

	start.2component = c(A=A_start,   B=B_start,   scal=scal_start,   xmid=xmid_start,   sens=sens_start,   xshift=xshift_start,   sigma=sigma_start)
	lower.2component = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sens=sens_low,   xshift=xshift_low,   sigma=sigma_low)
	upper.2component = c(A=A_high,   B=B_high,   scal=scal_high,   xmid=xmid_high,   sens=sens_high,   xshift=xshift_high,   sigma=sigma_high)
	assign("start.2component", start.2component, envir = .GlobalEnv)
	assign("lower.2component", lower.2component, envir = .GlobalEnv)
	assign("upper.2component", upper.2component, envir = .GlobalEnv)

	start.2component.2exp = c(A=A_start,   B=B_start,   scal=scal_start,   xmid=xmid_start,   sens_Exp1=sens_start,   sens_Exp2=sens_start,   xshift=xshift_start,   sigma=sigma_start)
	lower.2component.2exp = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sens_Exp1=sens_low,   sens_Exp2=sens_low,   xshift=xshift_low,   sigma=sigma_low)
	upper.2component.2exp = c(A=A_high,   B=B_high,   scal=scal_high,   xmid=xmid_high,   sens_Exp1=sens_high,   sens_Exp2=sens_high,   xshift=xshift_high,   sigma=sigma_high)
	assign("start.2component.2exp", start.2component.2exp, envir = .GlobalEnv)
	assign("lower.2component.2exp", lower.2component.2exp, envir = .GlobalEnv)
	assign("upper.2component.2exp", upper.2component.2exp, envir = .GlobalEnv)

	start.2component.3exp = c(A=A_start,   B=B_start,   scal=scal_start,   xmid=xmid_start,   sens_Exp1=sens_start,   sens_Exp2=sens_start,   sens_Exp3=sens_start,   xshift=xshift_start,   sigma=sigma_start)
	lower.2component.3exp = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sens_Exp1=sens_low,   sens_Exp2=sens_low,   sens_Exp3=sens_low,   xshift=xshift_low,   sigma=sigma_low)
	upper.2component.3exp = c(A=A_high,   B=B_high,   scal=scal_high,   xmid=xmid_high,   sens_Exp1=sens_high,   sens_Exp2=sens_high,   sens_Exp3=sens_high,   xshift=xshift_high,   sigma=sigma_high)
	assign("start.2component.3exp", start.2component.3exp, envir = .GlobalEnv)
	assign("lower.2component.3exp", lower.2component.3exp, envir = .GlobalEnv)
	assign("upper.2component.3exp", upper.2component.3exp, envir = .GlobalEnv)

	start.2component.4exp = c(A=A_start,   B=B_start,   scal=scal_start,   xmid=xmid_start,   sens_Exp1=sens_start,   sens_Exp2=sens_start,   sens_Exp3=sens_start,   sens_Exp4=sens_start,  xshift=xshift_start,   sigma=sigma_start)
	lower.2component.4exp = c(A=A_low,   B=B_low,   scal=scal_low,   xmid=xmid_low,   sens_Exp1=sens_low,   sens_Exp2=sens_low,   sens_Exp3=sens_low,   sens_Exp4=sens_low,  xshift=xshift_low,   sigma=sigma_low)
	upper.2component.4exp = c(A=A_high,   B=B_high,   scal=scal_high,   xmid=xmid_high,   sens_Exp1=sens_high,   sens_Exp2=sens_high,   sens_Exp3=sens_high,   sens_Exp4=sens_high,  xshift=xshift_high,   sigma=sigma_high)
	assign("start.2component.4exp", start.2component.4exp, envir = .GlobalEnv)
	assign("lower.2component.4exp", lower.2component.4exp, envir = .GlobalEnv)
	assign("upper.2component.4exp", upper.2component.4exp, envir = .GlobalEnv)
}
assign_start_lower_upper_values()

####################################################################################################
# Try each model and output the results for manual inspection.
####################################################################################################

cellline = assays[assays$Line == args_cellline,]

####################################################################################################
# Try 1component (uniphasic) and 2component (biphasic) models.
# Compare the two models. What is the p-value.
# If the p-value is significant, then there really is a difference between the two models, 
# and thus the 2component model is better than the 1component model and the data really is biphasic.
this_cellline_fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
this_cellline_fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
this_cellline.1component_vs_2component.pvalue = pchisq(2*(this_cellline_fit.2component$value - this_cellline_fit.1component$value), 2, lower.tail = FALSE)
this_cellline_fit.1component
this_cellline_fit.2component
this_cellline.1component_vs_2component.pvalue
convert_values_back_to_interpretable.1component(this_cellline_fit.1component)

# Example output for H226.
# Which model is appropriate for H226 human mesothelioma? Does it have one TC50 or two?
# 
# cellline = assays[assays$Line == "H226",]
# H226_fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# H226_fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# H226.1component_vs_2component.pvalue = pchisq(2*(H226_fit.2component$value - H226_fit.1component$value), 2, lower.tail = FALSE)
# H226_fit.1component
# H226_fit.2component
# H226.1component_vs_2component.pvalue
# convert_values_back_to_interpretable.1component(H226_fit.1component)
# 
# > H226_fit.1component
# 	        A          B       scal       xmid     sigma    value fevals
# L-BFGS-B 0.972522 0.05654325 0.04805096 -0.4618924 0.1379889 74.14672     50
# 	 gevals niter convcode  kkt1  kkt2 xtimes
# L-BFGS-B     50    NA        0 FALSE FALSE  0.056
# > H226_fit.2component
# 	         A          B       scal       xmid sens   xshift     sigma
# L-BFGS-B 0.9724325 0.05674729 0.04756713 -0.4623562    1 1.194868 0.1379909
# 	    value fevals gevals niter convcode  kkt1  kkt2 xtimes
# L-BFGS-B 74.14676     73     73    NA        0 FALSE FALSE  0.096
# > H226.1component_vs_2component.pvalue
# [1] 0.9999608
# > convert_values_back_to_interpretable.1component(H226_fit.1component)
# $TC50
# [1] 0.6300901

##### p-value = 1.00. 2component model is not better than 1component model for H226. H226 is not biphasic.
##### H226 TC50 = 0.630 mg/ml

####################################################################################################
# Now try modelling with adjustments for each experiment, for those cases needing this extra modelling.

exp = unique(cellline$Experiment)

if (length(exp) == 2) { # There are 2 different experiments.
	# Try 1component (uniphasic) and 2component (biphasic) models.
	# Compare the two models. What is the p-value.
	# If the p-value is significant, then there really is a difference between the two models, 
	# and thus the 2component model is better than the 1component model and the data really is biphasic.
	cellline$is_Exp1 <- ifelse(cellline$Experiment == exp[1], 1, 0)
	cellline$is_Exp2 <- ifelse(cellline$Experiment == exp[2], 1, 0)
	this_cellline_fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
	this_cellline_fit.2component.2exp = optimx(start.2component.2exp, function(par) loglik.2component.2exp(par, x=cellline$log10Conc, y=cellline$Viability, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3), lower=lower.2component.2exp, upper=upper.2component.2exp, method=c("L-BFGS-B"), control=list(fnscale = -1))
	this_cellline.1component_vs_2component.2exp.pvalue = pchisq(2*(HOFM01_fit.2component.2exp$value - HOFM01_fit.1component$value), 2, lower.tail = FALSE)
	this_cellline_fit.1component
	this_cellline_fit.2component.2exp
	this_cellline.1component_vs_2component.2exp.pvalue
	convert_values_back_to_interpretable.2component.2exp(HOFM01_fit.2component.2exp)
}

if (length(exp) == 3) { # There are 3 different experiments.
	# Try 1component (uniphasic) and 2component (biphasic) models.
	# Compare the two models. What is the p-value.
	# If the p-value is significant, then there really is a difference between the two models, 
	# and thus the 2component model is better than the 1component model and the data really is biphasic.
	cellline$is_Exp1 <- ifelse(cellline$Experiment == exp[1], 1, 0)
	cellline$is_Exp2 <- ifelse(cellline$Experiment == exp[2], 1, 0)
	cellline$is_Exp3 <- ifelse(cellline$Experiment == exp[3], 1, 0)
	this_cellline_fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
	this_cellline_fit.2component.3exp = optimx(start.2component.3exp, function(par) loglik.2component.3exp(par, x=cellline$log10Conc, y=cellline$Viability, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3), lower=lower.2component.3exp, upper=upper.2component.3exp, method=c("L-BFGS-B"), control=list(fnscale = -1))
	this_cellline.1component_vs_2component.3exp.pvalue = pchisq(2*(HOFM01_fit.2component.3exp$value - HOFM01_fit.1component$value), 2, lower.tail = FALSE)
	this_cellline_fit.1component
	this_cellline_fit.2component.3exp
	this_cellline.1component_vs_2component.3exp.pvalue
	convert_values_back_to_interpretable.2component.3exp(HOFM01_fit.2component.3exp)
}

if (length(exp) == 4) { # There are 4 different experiments.
	# Try 1component (uniphasic) and 2component (biphasic) models.
	# Compare the two models. What is the p-value.
	# If the p-value is significant, then there really is a difference between the two models, 
	# and thus the 2component model is better than the 1component model and the data really is biphasic.
	cellline$is_Exp1 <- ifelse(cellline$Experiment == exp[1], 1, 0)
	cellline$is_Exp2 <- ifelse(cellline$Experiment == exp[2], 1, 0)
	cellline$is_Exp3 <- ifelse(cellline$Experiment == exp[3], 1, 0)
	cellline$is_Exp4 <- ifelse(cellline$Experiment == exp[4], 1, 0)
	this_cellline_fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
	this_cellline_fit.2component.4exp = optimx(start.2component.4exp, function(par) loglik.2component.4exp(par, x=cellline$log10Conc, y=cellline$Viability, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3), lower=lower.2component.4exp, upper=upper.2component.4exp, method=c("L-BFGS-B"), control=list(fnscale = -1))
	this_cellline.1component_vs_2component.4exp.pvalue = pchisq(2*(HOFM01_fit.2component.4exp$value - HOFM01_fit.1component$value), 2, lower.tail = FALSE)
	this_cellline_fit.1component
	this_cellline_fit.2component.4exp
	this_cellline.1component_vs_2component.4exp.pvalue
	convert_values_back_to_interpretable.2component.4exp(HOFM01_fit.2component.4exp)
}

# Example output for HOFM01.
# Which model is appropriate for HOFM01 human non-cancer fibroblasts? Does it have one TC50 or two?
# 
# cellline = assays[assays$Line == "HOFM01",]
# cellline$is_Exp1 <- ifelse(cellline$Experiment == "20151227", 1, 0)
# cellline$is_Exp2 <- ifelse(cellline$Experiment == "20160101", 1, 0)
# cellline$is_Exp3 <- ifelse(cellline$Experiment == "20160108", 1, 0)
# HOFM01_fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# HOFM01_fit.2component.3exp = optimx(start.2component.3exp, function(par) loglik.2component.3exp(par, x=cellline$log10Conc, y=cellline$Viability, is_Exp1=cellline$is_Exp1, is_Exp2=cellline$is_Exp2, is_Exp3=cellline$is_Exp3), lower=lower.2component.3exp, upper=upper.2component.3exp, method=c("L-BFGS-B"), control=list(fnscale = -1))
# HOFM01.1component_vs_2component.3exp.pvalue = pchisq(2*(HOFM01_fit.2component.3exp$value - HOFM01_fit.1component$value), 2, lower.tail = FALSE)
# HOFM01_fit.1component
# HOFM01_fit.2component.3exp
# HOFM01.1component_vs_2component.3exp.pvalue
# convert_values_back_to_interpretable.2component.3exp(HOFM01_fit.2component.3exp)
# 
# > HOFM01_fit.1component
# 	        A           B      scal        xmid     sigma    value fevals
# L-BFGS-B 1.200066 -0.08267779 0.4032885 -0.04799031 0.1864056 20.35132     44
# 	 gevals niter convcode kkt1  kkt2 xtimes
# L-BFGS-B     44    NA        0 TRUE FALSE  0.044
# > HOFM01_fit.2component.3exp
# 	        A          B      scal       xmid sens_Exp1 sens_Exp2 sens_Exp3
# L-BFGS-B 1.155652 0.06401558 0.1127083 -0.5811668 0.7213794 0.5777401 0.1591135
# 	    xshift     sigma    value fevals gevals niter convcode kkt1  kkt2
# L-BFGS-B 0.8666875 0.1235078 52.45802     63     63    NA        0 TRUE FALSE
# 	 xtimes
# L-BFGS-B  0.144
# 
# > HOFM01.1component_vs_2component.3exp.pvalue
# [1] 1.138244e-14
# > convert_values_back_to_interpretable.2component.3exp(HOFM01_fit.2component.3exp)
# $sensitive_fraction_TC50.xmid
# [1] 0.5592455
# 
# $resistant_fraction_TC50.xmid
# [1] 1.330455
# 
# $sensitive_fraction_average
# [1] 0.4860777
# 
# $resistant_fraction_average
# [1] 0.5139223
# 
##### p-value = very small. 2component model is much better than 1component model for HOFM01. HOFM01 is biphasic.
##### HOFM01 resistant fraction's TC50 = 1.330 mg/ml
##### HOFM01 sensitive fraction's TC50 = 0.486 mg/ml
##### HOFM01 % of cells that are resistant = 51.3%

# Example output for MM05.
# Which model is appropriate for MM05 human mesothelioma? Does it have one TC50 or two?
# 
# cellline = assays[assays$Line == "MM05",]
# MM05_fit.1component = optimx(start.1component, function(par) loglik.1component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.1component, upper=upper.1component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# MM05_fit.2component = optimx(start.2component, function(par) loglik.2component(par, x=cellline$log10Conc, y=cellline$Viability), lower=lower.2component, upper=upper.2component, method=c("L-BFGS-B"), control=list(fnscale = -1))
# MM05.1component_vs_2component.pvalue = pchisq(2*(MM05_fit.2component$value - MM05_fit.1component$value), 2, lower.tail = FALSE)
# MM05_fit.1component
# MM05_fit.2component
# MM05.1component_vs_2component.pvalue
# convert_values_back_to_interpretable.1component(MM05_fit.1component)
# 
# > MM05_fit.1component
# 	         A          B       scal      xmid     sigma    value fevals
# L-BFGS-B 0.9482784 0.05017752 0.05390604 -0.448927 0.1340998 61.97863     48
# 	 gevals niter convcode  kkt1  kkt2 xtimes
# L-BFGS-B     48    NA        0 FALSE FALSE  0.044
# > MM05_fit.2component
# 	         A         B       scal       xmid      sens   xshift     sigma
# L-BFGS-B 0.9481795 0.0180738 0.05020768 -0.4530067 0.9618023 1.155636 0.1340469
# 	    value fevals gevals niter convcode  kkt1  kkt2 xtimes
# L-BFGS-B 62.01836     82     82    NA        0 FALSE FALSE  0.116
# > MM05.1component_vs_2component.pvalue
# [1] 0.9610489
# > convert_values_back_to_interpretable.1component(MM05_fit.1component)
# $TC50
# [1] 0.6383127
# 
##### p-value = 0.96. 2component model is not better than 1component model for MM05. MM05 is not biphasic.
##### MM05 TC50 = 0.638 mg/ml



