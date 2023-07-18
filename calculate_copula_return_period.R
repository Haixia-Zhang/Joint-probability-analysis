library(goft)
library(stats4)
library(MASS)
library(evd)
library(gumbel)
library(evir)
library(copula)
library(scatterplot3d)
library(plot3D)
library(tcltk)
library(Rdonlp2)
source("http://www.datall-analyse.nl/R/eva_max.R")

setwd("F://4_Paper//copula")

##################### 1.Copula fitting
para_sur <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_surge.csv")
para_wav <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_wave.csv")

# Generate the empty matrix
para_cop <- matrix(data = NA, nrow = 1665, ncol = 15)

# Fitting copula parameters and test parameters
for(i in 1:1665){
	print(i)
	# Input data
	nod = para_sur[i, 'nod']
	lon = para_sur[i, 'lon']
	lat = para_sur[i, 'lat']

	infile_mid = paste("mid/rect_year_max_select/", nod, ".csv", sep = "")
	data_sur_wav <- read.csv(infile_mid)
	data_sur_wav <- data_sur_wav[, 1:2]

	# Fitting copula
	myCop.clayton <- copula::archmCopula(family = "clayton", dim = 2, param = NA_real_)
	myCop.frank <- copula::archmCopula(family = "frank", dim = 2, param = NA_real_)
	myCop.gumbel <- copula::archmCopula(family = "gumbel", dim = 2, param = NA_real_)

	myMvd_clayton <- copula::mvdc(copula = myCop.clayton, margins = c("gev", "gev"), 
								  paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']),
                                                 list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))
    myMvd_frank <- copula::mvdc(copula = myCop.frank, margins = c("gev", "gev"), 
                                paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']),
                                               list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))
    myMvd_gumbel <- copula::mvdc(copula = myCop.gumbel, margins = c("gev", "gev"), 
                                 paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']),
                                                list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))

  # Fitting copula
	fit.Cop.clayton <- copula::fitCopula(myMvd_clayton@copula, pobs(data_sur_wav), method = "ml")
	fit.clayton.alpha <- coef(fit.Cop.clayton)
	aic_clayton <- round(AIC(fit.Cop.clayton), 3)
	bic_clayton <- round(BIC(fit.Cop.clayton), 3)
	p.clayton <- tryCatch(
		{gof.clayton <- copula::gofCopula(myMvd_frank@copula, pobs(data_sur_wav), method = c("Sn"), estim.method = c("itau"))
		p.clayton <- gof.clayton$p.value}, error = function(e){
		p.clayton <- 0
		}
	)

	fit.Cop.frank <- copula::fitCopula(myMvd_frank@copula, pobs(data_sur_wav), method="ml")
	fit.frank.alpha <- coef(fit.Cop.frank)
	aic_frank <- round(AIC(fit.Cop.frank), 3)
	bic_frank <- round(BIC(fit.Cop.frank), 3)
	p.frank <- tryCatch(
		{gof.frank <- copula::gofCopula(myMvd_frank@copula, pobs(data_sur_wav), method = c("Sn"), estim.method = c("itau"))
		p.frank <- gof.frank$p.value}, error = function(e){
		p.frank <- 0
		}
	)

	fit.Cop.gumbel <- copula::fitCopula(myMvd_gumbel@copula, pobs(data_sur_wav), method="ml")
	fit.gumbel.alpha <- coef(fit.Cop.gumbel)
	aic_gumbel <- round(AIC(fit.Cop.gumbel), 3)
	bic_gumbel <- round(BIC(fit.Cop.gumbel), 3)
	p.gumbel <- tryCatch(
		{gof.gumbel <- copula::gofCopula(myMvd_frank@copula,pobs(data_sur_wav),method = c("Sn"),estim.method = c("itau"))
		p.gumbel <- gof.gumbel$p.value}, error = function(e){
		p.gumbel <- 0
		}
	)

	# Save parameters
	para_cop[i, ] = c(nod, lon, lat,
					fit.clayton.alpha, aic_clayton, bic_clayton, p.clayton,
					fit.frank.alpha, aic_frank, bic_frank, p.frank,
                    fit.gumbel.alpha, aic_gumbel, bic_gumbel, p.gumbel)
}

# Save file
write.table(para_cop, "rst/bivariate_surge_wave/rect_bivariate_fitting_all_params.csv",
            row.names = FALSE,
            col.names = c("nod", "lon", "lat",
                          "clayton_theta", "clayton_aic", "clayton_bic", "clayton_p",
                          "frank_theta", "frank_aic", "frank_bic", "frank_p",
                          "gumbel_theta", "gumbel_aic", "gumbel_bic", "gumbel_p"), sep = ",")


##################### 2.Copula drafting
para_sur <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_surge.csv")
para_wav <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_wave.csv")
para_cop <- read.csv("rst/bivariate_surge_wave/rect_bivariate_fitting_all_params.csv")

# Drafting
for(i in 1:1665){
	print(i)

	# Input data
	nod = para_sur[i, 'nod']

	infile_mid = paste("mid/rect_year_max_select/", nod, ".csv", sep = "")
	data_sur_wav <- read.csv(infile_mid)
	data_sur_wav <- data_sur_wav[, 1:2]
  
	# Copula fitting
	sur_min <- qgev(0.001, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	sur_max <- qgev(0.999, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	wav_min <- qgev(0.001, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
	wav_max <- qgev(0.999, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
	myCop.fit <- copula::archmCopula(family = "gumbel", dim = 2, param = para_cop[i, "gumbel_theta"])
	myMvd.fit <- copula::mvdc(copula = myCop.fit, margins = c("gev", "gev"), 
                              paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']), 
                                             list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))

	# Mariginal fittiing
	sur <- dgev(seq(sur_min, sur_max, length=100), xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	wav <- dgev(seq(wav_min, wav_max, length=100), xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])

	# Save data
	para_pdf_cdf <- matrix(data = NA, nrow = n * n, ncol = 4)
	para_pdf_cdf[, 1] <- X[, 1]
	para_pdf_cdf[, 2] <- X[, 2]

	# pdf/cdf
	n <- 50
	X <- matrix(data = NA, nrow = n * n, ncol = 2)
	X[, 1] <- rep(seq(sur_min, sur_max, length = n), n)
	X[, 2] <- rep(seq(wav_min, wav_max, length = n), each = n)

	dmyMvd <- dMvdc(X, myMvd.fit)
	dZ <- matrix(data = NA, nrow = n * n, ncol = 1)
	dZ[, 1] <- dmyMvd
	para_pdf_cdf[, 3] <- dZ[, 1]  
	dZ <- matrix(dmyMvd, ncol = n)

	pmyMvd <- pMvdc(X, myMvd.fit)
	pZ <- matrix(data = NA, nrow = n * n, ncol = 1)
	pZ[, 1] <- pmyMvd
	para_pdf_cdf[, 4] <- pZ[, 1]
	pZ <- matrix(pmyMvd, ncol = n)

	pdf_cdf_file = paste("rst/bivariate_surge_wave_PDF_CDF/copula_pdf_cdf_", nod, ".csv", sep = "")
	write.table(para_pdf_cdf, pdf_cdf_file, row.names = FALSE,
                col.names = c("surge_height", "wave_height", "copula_PDF", "copula_CDF"), sep = ",")

	# Drawing pdf and cdf
	png_a = paste("jpg/copula_pdf_cdf/copula_pdf_cdf_", nod, ".png", sep = "")
	png(file = png_a, width = 3800, height = 2500, res = 72 * 2)
	op <- par(mfcol = c(1, 2))

	dp <- graphics::persp (seq(sur_min, sur_max, length=n),
							seq(wav_min, wav_max, length=n),
							dZ,
							#xlim=c(0, max(data_sur_wav[,1])),
							#ylim=range(data_sur_wav[,2], finite = TRUE),
							#zlim=c(0, 1),
							theta = 330, phi = 15,
							expand = 0.6, col = "white",
							#ltheta = 120, lphi=0, shade = 0.3,
							ticktype = "detailed",
							xlab = "Storm Surge Height (m)", ylab = "Significant Wave Height (m)", zlab = "PDF",
							cex.lab = 2.3, cex.axis = 2.3)

	pp <- graphics::persp (seq(sur_min, sur_max, length=n),
							seq(wav_min, wav_max, length=n),
							pZ,
							#xlim = c(0, 81),
							#ylim=c(0, 1277),
							zlim=c(0, 1),
							theta = 330, phi = 15,
							expand = 0.6, col = "white",
							#ltheta = 120, lphi=0, shade = 0.3,
							ticktype = "detailed",
							xlab = "Storm Surge Height (m)", ylab = "Significant Wave Height (m)", zlab = "CDF",
							cex.lab = 2.3, cex.axis = 2.3)

	par(op)
	dev.off( )

	# Plotting random number scatter
	x.samp <- rMvdc(1000, myMvd.fit)
	plot(x.samp)

	sample_file = paste("rst/bivariate_surge_wave_sample/copula_samples_", nod, ".csv", sep = "")
	write.table(x.samp, sample_file, row.names = FALSE, col.names = c("surge_height", "wave_height"), sep = ",")

	png_b = paste("jpg/copula_sample/copula_samples_", nod, ".png", sep = "")
	png(file = png_b, width = 1000, height = 600, res = 54)
	op <- par(mfcol = c(1, 1))

	dp <- plot(x.samp[, 1], x.samp[, 2],
				xlim = c(0, max(x.samp[, 1]) + 0.5),
				ylim = c(0, max(x.samp[, 2]) + 0.5),
				xlab = "Surge height (m)",
				ylab = "Significant wave height (m)",
				cex.lab = 2.5, cex.axis = 2.5)

	par(op)
	dev.off( )
}


##################### 3.Calculate return period
para_sur <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_surge.csv")
para_wav <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_wave.csv")
para_cop <- read.csv("rst/bivariate_surge_wave/rect_bivariate_fitting_all_params.csv")

for(i in 1:1665){
	print(i)
	# Input data
	nod = para_sur[i, 'nod']

	infile_mid = paste("mid/rect_year_max_select/", nod, ".csv", sep = "")
	data_sur_wav <- read.csv(infile_mid)
	data_sur_wav <- data_sur_wav[, 1:2]

	# Copula fitting
	sur_min <- qgev(0.001, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	sur_max <- qgev(0.999, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	wav_min <- qgev(0.001, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
	wav_max <- qgev(0.999, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
	myCop.fit <- copula::archmCopula(family = "gumbel", dim = 2, param = para_cop[i, "gumbel_theta"])
	myMvd.fit <- copula::mvdc(copula = myCop.fit, margins = c("gev", "gev"),
							  paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']),
											 list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))

	n <- 50
	X <- matrix(data = NA, nrow = n * n, ncol = 2)
	X[, 1] <- rep(seq(sur_min, sur_max, length = n), n)
	X[, 2] <- rep(seq(wav_min, wav_max, length = n), each = n)

	p_s <- pgev(X[, 1], xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	p_s <- matrix(p_s, ncol = n)
	p_w <- pgev(X[, 2], xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
	p_w <- matrix(p_w, ncol = n)

	p_copula <- pMvdc(X, myMvd.fit)
	p_copula <- matrix(p_copula, ncol = n)

	ReturnPeriod_and <- 1 / (1 - p_s - p_w + p_copula)
	ReturnPeriod_or <- 1 / (1 - p_copula)

	# Save data
	data_out <- matrix(data = NA, nrow = n * n, ncol = 4)
	data_out[, 1] <- X[, 1]
	data_out[, 2] <- X[, 2]
	data_out[, 3] <- ReturnPeriod_and
	data_out[, 4] <- ReturnPeriod_or

	rstPath = paste("rst/bivariate_surge_wave_RP/bivariate_fitting_RP_", nod, ".csv", sep = "")
	write.table(data_out, rstPath, row.names = FALSE,
                col.names = c("surge_height", "wave_height", "RP_AND", "RP_OR"), sep = ",")

	# return period or / and
	png_c = paste("jpg/copula3d_RP_/copula_RP_or_and_", nod, ".png", sep = "")
	png(file = png_c, width = 1000, height = 1000, res = 72 * 2)

	par(mfrow = c(1, 1), mgp=c(3,1,0))
	print(sur_min)
	print(sur_max)
	contour(seq(0, 4, length = n),
                seq(0, 10, length = n),
                ReturnPeriod_and,
                nlevels = 10,
                col = 'black',
                xlab = "Surge height (m)",
                ylab = "Significant wave height (m)",
                cex.lab = 1, cex.axis = 1)
	par(new=TRUE)
  
	contour(seq(0, 4, length = n),
                seq(0, 10, length = n),
                ReturnPeriod_or,
                nlevels = 10,
                col = 'red',
                xlab = "Surge height (m)",
                ylab = "Significant wave height (m)",
                cex.lab = 1, cex.axis = 1)

  dev.off( )
}


##################### 4.Calculate simultaneous probability/joint probability

para_sur <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_surge.csv")
para_wav <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_wave.csv")
para_cop <- read.csv("rst/bivariate_surge_wave/rect_bivariate_fitting_all_params.csv")

# Generate the empty matrix
sur_rp_out <- matrix(data = NA, nrow = 1665, ncol = 1002)
wav_rp_out <- matrix(data = NA, nrow = 1665, ncol = 1002)

p <- c(0.5)

for (n in 1:1665){
	p <- append(p, 1 - 1 / n)
}

# Calculate return period of marginal function
for(i in 1:1665){
	print(i)
	# Input data
	nod = para_sur[i, 'nod']
	lon = para_sur[i, 'lon']
	lat = para_sur[i, 'lat']

	sur_rp <- qgev(p, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	wav_rp <- qgev(p, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])

	sur_rp_out[i, 1] = nod
	sur_rp_out[i, 2] = lon
	sur_rp_out[i, 3] = lat
	sur_rp_out[i, 4:1002] <- sur_rp

	wav_rp_out[i, 1] = nod
	wav_rp_out[i, 2] = lon
	wav_rp_out[i, 3] = lat
	wav_rp_out[i, 4:1002] <- wav_rp
}

write.table(sur_rp_out, "rst/marginal_surge_wave_RP/rect_marginal_fitting_GEV_RP_surge.csv",
            row.names = FALSE, col.names = c("nod", "lon", "lat", 2:1000), sep = ",")

write.table(wav_rp_out, "rst/marginal_surge_wave_RP/rect_marginal_fitting_GEV_RP_wave.csv",
            row.names = FALSE, col.names = c("nod", "lon", "lat", 2:1000), sep = ",")

# Calculation probability
P_and_rp <- matrix(data = NA, nrow = 1665, ncol = 39)
P_or_rp <- matrix(data = NA, nrow = 1665, ncol = 39)

p_rp <- c()

for (m in c(5, 10, 20, 50, 100, 200)){
	p_rp <- append(p_rp, 1 - 1 / m)
}

for(i in 1:1665){
	print(i)
	# Input data
	nod = para_sur[i, 'nod']
	lon = para_sur[i, 'lon']
	lat = para_sur[i, 'lat']

	infile_mid = paste("mid/rect_year_max_select/", nod, ".csv", sep = "")
	data_sur_wav <- read.csv(infile_mid)
	data_sur_wav <- data_sur_wav[, 1:2]

	sur_rp <- qgev(p_rp, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	wav_rp <- qgev(p_rp, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])

	myCop.fit <- copula::archmCopula(family = "gumbel", dim = 2, param = para_cop[i, "gumbel_theta"])
	myMvd.fit <- copula::mvdc(copula = myCop.fit, margins = c("gev", "gev"), 
							  paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']), 
                                             list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))

	n <- 6
	X <- matrix(data = NA, nrow = n * n, ncol = 2)
	X[, 1] <- rep(sur_rp, n)
	X[, 2] <- rep(wav_rp, each = n)

	p_s <- pgev(X[, 1], xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	p_s <- matrix(p_s, ncol = n)
	p_w <- pgev(X[, 2], xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
	p_w <- matrix(p_w, ncol = n)

	p_copula <- pMvdc(X, myMvd.fit)
	p_copula <- matrix(p_copula, ncol = n)

	probability_and <- (1 - p_s - p_w + p_copula)
	probability_or <- (1 - p_copula)

	probability_and <- matrix(probability_and, ncol = 36)
	probability_or <- matrix(probability_or, ncol = 36)

	# Save data
	P_and_rp[i, 1] = nod
	P_and_rp[i, 2] = lon
	P_and_rp[i, 3] = lat
	P_and_rp[i, 4:39] <- probability_and

	P_or_rp[i, 1] = nod
	P_or_rp[i, 2] = lon
	P_or_rp[i, 3] = lat
	P_or_rp[i, 4:39] <- probability_or
}

write.table(P_and_rp, "rst/bivariate_surge_wave/rect_bivariate_fitting_probability_rp_and.csv",
            row.names = FALSE, col.names = c("nod", "lon", "lat",
                                             "s5w5", "s10w5", "s20w5", "s50w5", "s100w5", "s200w5",
                                             "s5w10", "s10w10", "s20w10", "s50w10", "s100w10", "s200w10",
                                             "s5w20", "s10w20", "s20w20", "s50w20", "s100w20", "s200w20",
                                             "s5w50", "s10w50", "s20w50", "s50w50", "s100w50", "s200w50",
                                             "s5w100", "s10w100", "s20w100", "s50w100", "s100w100", "s200w100",
                                             "s5w200", "s10w200", "s20w200", "s50w200", "s100w200", "s200w200"), sep = ",")

write.table(P_or_rp, "rst/bivariate_surge_wave/rect_bivariate_fitting_probability_rp_or.csv",
            row.names = FALSE, col.names = c("nod", "lon", "lat",
                                             "s5w5", "s10w5", "s20w5", "s50w5", "s100w5", "s200w5",
                                             "s5w10", "s10w10", "s20w10", "s50w10", "s100w10", "s200w10",
                                             "s5w20", "s10w20", "s20w20", "s50w20", "s100w20", "s200w20",
                                             "s5w50", "s10w50", "s20w50", "s50w50", "s100w50", "s200w50",
                                             "s5w100", "s10w100", "s20w100", "s50w100", "s100w100", "s200w100",
                                             "s5w200", "s10w200", "s20w200", "s50w200", "s100w200", "s200w200"), sep = ",")


##################### 5.Surge and Wave Division
para_sur <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_surge.csv")
para_wav <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_wave.csv")
para_cop <- read.csv("rst/bivariate_surge_wave/rect_bivariate_fitting_all_params.csv")

# Generate the empty matrix
P_level <- matrix(data = NA, nrow = 1665, ncol = 28)

# Calculation probability
for(i in 1:1665){
	print(i)
	# Input data
	nod = para_sur[i, 'nod']
	lon = para_sur[i, 'lon']
	lat = para_sur[i, 'lat']

	P_level[i, 1] = nod
	P_level[i, 2] = lon
	P_level[i, 3] = lat

	# Set class of surges and waves
	sur_min <- qgev(0.0001, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	sur_max <- qgev(0.9999, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	wav_min <- qgev(0.0001, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
	wav_max <- qgev(0.9999, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])

	sur_level <- c(0, 1.0, 1.5, 2.0, 2.5, 40)
	# wav_level <- c(0, 2.5, 3.5, 4.5, 6.0, 40)
	wav_level <- c(0, 4.0, 6.0, 9.0, 14.0, 40)
	# Copula fitting
	myCop.fit <- copula::archmCopula(family = "gumbel", dim = 2, param = para_cop[i, "gumbel_theta"])
	myMvd.fit <- copula::mvdc(copula = myCop.fit, margins = c("gev", "gev"), 
                              paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']), 
                                             list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))

	# Calculate segmentation probability
	for (j in 1:5){
		x1 <- sur_level[j]
		x2 <- sur_level[j + 1]

		if (x2 <= sur_max){
			for (k in 1:5){
				y1 <- wav_level[k]
				y2 <- wav_level[k + 1]

				if (y2 <= wav_max){
					xy <- c(pMvdc(c(x2, y2), myMvd.fit), pMvdc(c(x2, y1), myMvd.fit),
							pMvdc(c(x1, y2), myMvd.fit), pMvdc(c(x1, y1), myMvd.fit))
				xy[is.nan(xy)] <- 0
				P_level[i, 3 + 5 * (j - 1) + k] <- xy[1] - xy[2] - xy[3] + xy[4]
				}
				else{
					y2 <- wav_max
					xy <- c(pMvdc(c(x2, y2), myMvd.fit), pMvdc(c(x2, y1), myMvd.fit),
							pMvdc(c(x1, y2), myMvd.fit), pMvdc(c(x1, y1), myMvd.fit))
					xy[is.nan(xy)] <- 0
					P_level[i, 3 + 5 * (j - 1) + k] <- xy[1] - xy[2] - xy[3] + xy[4]
					break;
				}
			}
		}
		else{
			x2 <- sur_max
			for (k in 1:5){
				y1 <- wav_level[k]
				y2 <- wav_level[k + 1]

				if (y2 <= wav_max){
					xy <- c(pMvdc(c(x2, y2), myMvd.fit), pMvdc(c(x2, y1), myMvd.fit),
							pMvdc(c(x1, y2), myMvd.fit), pMvdc(c(x1, y1), myMvd.fit))
					xy[is.nan(xy)] <- 0
					P_level[i, 3 + 5 * (j - 1) + k] <- xy[1] - xy[2] - xy[3] + xy[4]
					}
				else{
					y2 <- wav_max
					xy <- c(pMvdc(c(x2, y2), myMvd.fit), pMvdc(c(x2, y1), myMvd.fit),
							pMvdc(c(x1, y2), myMvd.fit), pMvdc(c(x1, y1), myMvd.fit))
					xy[is.nan(xy)] <- 0
					P_level[i, 3 + 5 * (j - 1) + k] <- xy[1] - xy[2] - xy[3] + xy[4]
					break;
					}
			}
			break;
		}
	}
}

write.table(P_level, "rst/bivariate_surge_wave/rect_bivariate_fitting_zonation.csv",
            row.names = FALSE, col.names = c("nod", "lon", "lat",
                                             "s1w1", "s1w2", "s1w3", "s1w4", "s1w5",
                                             "s2w1", "s2w2", "s2w3", "s2w4", "s2w5",
                                             "s3w1", "s3w2", "s3w3", "s3w4", "s3w5",
                                             "s4w1", "s4w2", "s4w3", "s4w4", "s4w5",
                                             "s5w1", "s5w2", "s5w3", "s5w4", "s5w5"), sep = ",")


##################### 6.Design values for surge and wave
para_sur <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_surge.csv")
para_wav <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_wave.csv")
para_cop <- read.csv("rst/bivariate_surge_wave/rect_bivariate_fitting_all_params.csv")

# Generate the empty matrix
value_rp_or <- matrix(data = NA, nrow = 1665, ncol = 21)

# Calculate return periods of marginal function
for(i in 1:1665){
	print(i)
	# Input data
	nod = para_sur[i, 'nod']
	lon = para_sur[i, 'lon']
	lat = para_sur[i, 'lat']

	value_rp_or[i, 1] = nod
	value_rp_or[i, 2] = lon
	value_rp_or[i, 3] = lat

	myCop.fit <- copula::archmCopula(family = "gumbel", dim = 2, param = para_cop[i, "gumbel_theta"])
	myMvd.fit <- copula::mvdc(copula = myCop.fit, margins = c("gev", "gev"), 
							  paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']), 
												  list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))
	j = 4
	for (n in c(5, 10, 20, 50, 100, 200)){
		obj <- function(X){
		1 / (1 - pgev(X[1], xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
			- pgev(X[2], xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
			+ pMvdc(c(X[1], X[2]), myMvd.fit))
		}
		par.l <- c(0, 0)
		par.u <- c(40, 40)

		nlcon1 <- function(X){
		  1 / (1 - pMvdc(c(X[1], X[2]), myMvd.fit))
		}
		nlin.l <- c(n)
		nlin.u <- c(n)

		p <- c(1, 1)
		ret <- donlp2(p, obj, par.u=par.u, par.l=par.l, nlin=list(nlcon1), nlin.u=nlin.u, nlin.l=nlin.l)

		if (ret$message == "KT-conditions satisfied, no further correction computed"){
			# save data
			cat(ret$par, ret$fx, '\n')

			value_rp_or[i, j] = ret$par[1]
			value_rp_or[i, j + 1] = ret$par[2]
			value_rp_or[i, j + 2] = ret$fx[1]
		}

		else{
			value_rp_or[i, j] = NaN
			value_rp_or[i, j + 1] = NaN
			value_rp_or[i, j + 2] = NaN
		}

		j <- j + 3
	}
}

write.table(value_rp_or, "rst/bivariate_surge_wave/rect_bivariate_fitting_joint_return_period_design_value.csv",
            row.names = FALSE, col.names = c("nod", "lon", "lat",
                                             "surge_5a", "wave_5a", "RP_and_5a",
                                             "surge_10a", "wave_10a", "RP_and_10a",
                                             "surge_20a", "wave_20a", "RP_and_20a",
                                             "surge_50a", "wave_50a", "RP_and_50a",
                                             "surge_100a", "wave_100a", "RP_and_100a",
                                             "surge_200a", "wave_200a",  "RP_and_200a"), sep = ",")


## Given RP_or, solve for the minimum value of RP_and
myCop.fit <- copula::archmCopula(family = "gumbel", dim = 2, param = 1.897174029)
myMvd.fit <- copula::mvdc(copula = myCop.fit, margins = c("gev", "gev"), 
                          paramMargins = list(list(xi = 0.480647084, mu = 0.25654794, sigma = 0.173433132), 
                                              list(xi = 0.097571295, mu = 2.094932497, sigma = 1.027092584)))
obj <- function(X){
	1 / (1 - (pgev(X[1], xi = 0.480647084, mu = 0.25654794, sigma = 0.173433132)) -
	(pgev(X[2], xi = 0.097571295, mu = 2.094932497, sigma = 1.027092584)) + pMvdc(c(X[1], X[2]), myMvd.fit))
}
par.l <- c(0, 0)
par.u <- c(10, 14)

nlcon1 <- function(X){
	1 / (1 - pMvdc(c(X[1], X[2]), myMvd.fit))
}
nlin.l <- c(20)
nlin.u <- c(20)

p <- c(0, 0)
ret <- donlp2(p, obj, par.u=par.u, par.l=par.l, nlin=list(nlcon1), nlin.u=nlin.u, nlin.l=nlin.l)
ret


## Given RP_and, solve for the maximum value of RP_or
myCop.fit <- copula::archmCopula(family = "gumbel", dim = 2, param = 1.897174029)
myMvd.fit <- copula::mvdc(copula = myCop.fit, margins = c("gev", "gev"), 
                          paramMargins = list(list(xi = 0.480647084, mu = 0.25654794, sigma = 0.173433132), 
                                              list(xi = 0.097571295, mu = 2.094932497, sigma = 1.027092584)))
obj <- function(X){
		(1 - pMvdc(c(X[1], X[2]), myMvd.fit))
	}
par.l <- c(0, 0)
par.u <- c(10, 14)

nlcon1 <- function(X){
	1 / (1 - (pgev(X[1], xi = 0.480647084, mu = 0.25654794, sigma = 0.173433132)) -
	(pgev(X[2], xi = 0.097571295, mu = 2.094932497, sigma = 1.027092584)) + pMvdc(c(X[1], X[2]), myMvd.fit))
}
nlin.l <- c(50)
nlin.u <- c(50)

p <- c(0, 0)
ret <- donlp2(p, obj, par.u=par.u, par.l=par.l, nlin=list(nlcon1), nlin.u=nlin.u, nlin.l=nlin.l)
ret$par


##################### 7.Calculate conditional probability

para_sur <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_surge.csv")
para_wav <- read.csv("rst/bivariate_surge_wave/rect_marginal_fitting_GEV_params_wave.csv")
para_cop <- read.csv("rst/bivariate_surge_wave/rect_bivariate_fitting_all_params.csv")

# Calculate conditional probability under different return periods
P_sur_rp <- matrix(data = NA, nrow = 1665, ncol = 39)
P_wav_rp <- matrix(data = NA, nrow = 1665, ncol = 39)

p_rp <- c()

for (m in c(5, 10, 20, 50, 100, 200)){
	p_rp <- append(p_rp, 1 - 1 / m)
}

for(i in 1:1665){
	print(i)
	# Input data
	nod = para_sur[i, 'nod']
	lon = para_sur[i, 'lon']
	lat = para_sur[i, 'lat']

	sur_rp <- qgev(p_rp, xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	wav_rp <- qgev(p_rp, xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
  
	myCop.fit <- copula::archmCopula(family = "gumbel", dim = 2, param = para_cop[i, "gumbel_theta"])
	myMvd.fit <- copula::mvdc(copula = myCop.fit, margins = c("gev", "gev"), 
							  paramMargins = list(list(xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c']), 
                                             list(xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])))

	n <- 6
	X <- matrix(data = NA, nrow = n * n, ncol = 2)
	X[, 1] <- rep(sur_rp, n)
	X[, 2] <- rep(wav_rp, each = n)

	p_s_ <- pgev(X[, 1], xi = - para_sur[i, 'a'], mu = para_sur[i, 'b'], sigma = para_sur[i, 'c'])
	p_s <- matrix(p_s_, ncol = n)
	p_w_ <- pgev(X[, 2], xi = - para_wav[i, 'a'], mu = para_wav[i, 'b'], sigma = para_wav[i, 'c'])
	p_w <- matrix(p_w_, ncol = n)

	p_copula <- pMvdc(X, myMvd.fit)
	p_copula <- matrix(p_copula, ncol = n)

	probability_sur <- (1 - p_s - p_w + p_copula) / (1 - p_s)
	probability_wav <- (1 - p_s - p_w + p_copula) / (1 - p_w)

	probability_sur <- matrix(probability_sur, ncol = 36)
	probability_wav <- matrix(probability_wav, ncol = 36)

	# Save data
	P_sur_rp[i, 1] = nod
	P_sur_rp[i, 2] = lon
	P_sur_rp[i, 3] = lat
	P_sur_rp[i, 4:39] <- probability_sur

	P_wav_rp[i, 1] = nod
	P_wav_rp[i, 2] = lon
	P_wav_rp[i, 3] = lat
	P_wav_rp[i, 4:39] <- probability_wav
}

write.table(P_sur_rp, "rst/bivariate_surge_wave/rect_bivariate_fitting_conditional_probability_surge.csv",
            row.names = FALSE, col.names = c("nod", "lon", "lat",
                                             "w5_s5", "w5_s10", "w5_s20", "w5_s50", "w5_s100", "w5_s200",
                                             "w10_s5", "w10_s10", "w10_s20", "w10_s50", "w10_s100", "w10_s200",
                                             "w20_s5", "w20_s10", "w20_s20", "w20_s50", "w20_s100", "w20_s200",
                                             "w50_s5", "w50_s10", "w50_s20", "w50_s50", "w50_s100", "w50_s200",
                                             "w100_s5", "w100_s10", "w100_s20", "w100_s50", "w100_s100", "w100_s200",
                                             "w200_s5", "w200_s10", "w200_s20", "w200_s50", "w200_s100", "w200_s200"), sep = ",")

write.table(P_wav_rp, "rst/bivariate_surge_wave/rect_bivariate_fitting_conditional_probability_wave.csv",
            row.names = FALSE, col.names = c("nod", "lon", "lat",
                                             "s5_w5", "s10_w5", "s20_w5", "s50_w5", "s100_w5", "s200_w5",
                                             "s5_w10", "s10_w10", "s20_w10", "s50_w10", "s100_w10", "s200_w10",
                                             "s5_w20", "s10_w20", "s20_w20", "s50_w20", "s100_w20", "s200_w20",
                                             "s5_w50", "s10_w50", "s20_w50", "s50_w50", "s100_w50", "s200_w50",
                                             "s5_w100", "s10_w100", "s20_w100", "s50_w100", "s100_w100", "s200_w100",
                                             "s5_w200", "s10_w200", "s20_w200", "s50_w200", "s100_w200", "s200_w200"), sep = ",")
