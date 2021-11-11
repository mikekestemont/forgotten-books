
#############################################################################
#                                                                           #
# Copia.R: R script for inferring culture diversity based on abundance data #
#                                                                           #
#############################################################################

# Copia.R represents the R version of python package Copia (available from PyPI) developed by Folgert Karsdorp. Copia means "abundance" in Latin. 
# Copia.R facilitates computation and graphics presented in the paper "Forgotten Books: The Application of Unseen Species Models
# to the Survival of Culture (2021) by Mike Kestemont et al. (2021). 
# Copia.R computes the Chao1 species richness estimator, minimum sample size required to reach the Chao1 estimate, and plots diversity profile,
# and evenness profiles based on data of 6 vernaculars/corpora used in Kestemont et al. (2021).
# If you use Copia.R for publishing papers, please cite Kestemoent et al. (2021).
 
#----------------------------------
# Packages needed to run Copia.R
#----------------------------------
# You must first install two packages (tidyverse and devtools) from CRAN; and two packages (iNEXT.3D and iNEXT.4steps) from Anne Chao's github
install.packages("tidyverse")
install.packages("devtools")
install_github("AnneChao/iNEXT.3D")  ## Press 'Enter' to skip number selection.
install_github("AnneChao/iNEXT.4steps")  ## Press 'Enter' to skip number selection.

# Import packages
library(iNEXT.3D)
library(iNEXT.4steps)
library(tidyverse)
library(devtools)

#-------------------------
# Load data of 6 corpora
#-------------------------
abundance <- readRDS("abundance.rds")

#------------------------------------
# The Chao1 species richness estimate
#------------------------------------
# Chao_1.est is a function to obtain the Chao1 species richness estimate 
# Chao_1.est.table is a function to construct a table including f1, f2, S.obs, n, Chao1, CH1(= S.obs/Chao1), CH1-lCI, CH1-uCI
# parameter y a vector of class/species sample frequencies. 
# parameter nboot a number specifying the number of bootstrap replications; default is 1000; change it to 10000 for more accurate results, 
# but it may be time-consuming.
# return a table

Chao_1.est = function(x){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  S_obs = sum(x>0)
  S_Chao1 = S_obs + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
  S_Chao1
}

Chao_1.est.table <- function(y, nboot = 1000){
  y = y[y>0]
  n = sum(y)      
  f1 = sum(y==1)  
  f2 = sum(y==2)
  S.obs = sum(y>0)
  Chao1 = Chao_1.est(y)
  if(nboot >1){
  Bt = rmultinom(nboot, n, iNEXT.3D:::EstiBootComm.Ind(y))
  chao1s = apply(Bt, 2, function(x) Chao_1.est(x))
  chao1s.UCL = chao1s %>% quantile(., 0.975)
  chao1s.LCL = chao1s %>% quantile(., 0.025)
  ch1_LCL = S.obs/(Chao1 + (chao1s.UCL - mean(chao1s)))
  ch1_UCL = S.obs/(Chao1 - (mean(chao1s) - chao1s.LCL))
  table = matrix(c(f1, f2, S.obs, n, round(Chao1, 3), round(S.obs/Chao1, 3), round(ch1_LCL, 3),
                   round(ch1_UCL, 3)), ncol = 8) %>% as.data.frame()
  colnames(table) <- c("f1", "f2", "S.obs", "n", "S.Chao1", "CH1", "CH1-lCI", "CH1-uCI")
  }else{
    table = matrix(c(f1, f2, S.obs, n, round(Chao1, 3), round(S.obs/Chao1, 3), NA,NA), ncol = 8) %>% as.data.frame()
    colnames(table) <- c("f1", "f2", "S.obs", "n", "S.Chao1", "CH1", "CH1-lCI", "CH1-uCI")
  }
  return(table)
}

#----------------------------------------------------------
# Minimum sample size required to reach the Chao1 estimate
#----------------------------------------------------------
# MS.x is a function to solve equation 2*f1*(1+x)=exp(x*(2*f2/f1))
# Msize : minimum sample size required to reach the Chao1 estimate
# MS.est.table is a function to construct a table including n, Msize(=(1+x)*n), MS(= n/Msize), MS-lCI,MS-uCI
# parameter y a vector of species sample frequencies. 
# parameter nboot a number specifying the number of bootstrap replications; default is 1000; change it to 10000 for more accurate results, 
# but it may be time-consuming.
# return a table

MS.x <- function(x){
  f1 = sum(x==1)
  f2 = sum(x==2)
  n = sum(x)
  if (f1 == 0 | f2 == 0) 
  stop("singlton and dobleton must be non-zero values")
  f <- function(y){2*f1*(1+ y)-exp(y*(2*f2/f1))}
  root <- uniroot(f, c(0,10000*n), tol = 0.000001)
  xs <- root$root
  MS = xs
  MS
}


MS.est.table <- function(y, nboot = 1000){
  n <- sum(y)
  f1 = sum(y==1)  
  f2 = sum(y==2)
  xs = MS.x(y)
  Msize = (1+xs)*n
  MS = n/(n*(1 + xs))
  if(nboot > 1){
  Bt = rmultinom(nboot, n, iNEXT.3D:::EstiBootComm.Ind(y))
  ms = apply(Bt, 2, function(y) MS.x(y))
  ms.UCL =  quantile(ms, 0.975)
  ms.LCL =  quantile(ms, 0.025)
  MS_lCI = 1/(1 + (xs + (ms.UCL - mean(ms))))
  MS_uCI = 1/(1 + (xs - (mean(ms) - ms.LCL)))
  table = matrix(c(n, Msize, round(MS, 3), round(MS_lCI, 3), round(MS_uCI, 3)), ncol = 5) %>% as.data.frame()
  colnames(table) = c("n", "Msize", "MS", "MS-lCI", "MS-uCI")
  }else{
    table = matrix(c(n, Msize, round(MS, 3), NA, NA), ncol = 5) %>% as.data.frame()
    colnames(table) = c("n", "Msize", "MS", "MS-lCI", "MS-uCI")
  }
  return(table)
}

#----------------------------------
# Combine Chao1 table and MS table
#----------------------------------
# combine_table is a function to combine Chao1 table and MS table
# parameter y a vector of species sample frequencies. 
# parameter nboot a number specifying the number of bootstrap replications.
# return a table
combine_table <- function(y, nboot=1000){
  Chao1.table = Chao_1.est.table(y, nboot)
  MS.table = MS.est.table(y, nboot)
  table = cbind(Chao1.table, MS.table[-1])
  return(table)
}

#------------------------------------------------------------
# Output table for 6 corpora data (Table 1 in Kestemont et al.)
#------------------------------------------------------------

table = t(sapply(abundance, combine_table, nboot = 1000))
table

#--------------------------------------------------------
# Diversity profile plot (Figure 3A in Kestemont et al.)
#--------------------------------------------------------

outcome = rbind(iNEXT.3D::asy3D(unlist(abundance), q = seq(0,3,0.2)), iNEXT.3D::obs3D(unlist(abundance), q = seq(0,3,0.2)))

ggplot(outcome, aes(x = Order.q, y = qD, colour = Method, 
                    fill = Method))+ 
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, fill = Method), linetype = 0, alpha = 0.2)+
  labs(x = "Order q", y = "Hill numbers") + theme_bw() + 
  theme(legend.position = "bottom", legend.box = "vertical", legend.key.width = unit(1.2, "cm"), 
        legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(-10, -10, -5, -10), text = element_text(size = 16), 
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) + 
  guides(linetype = guide_legend(keywidth = 2.5))




#------------------------------------------------------------
# Evenness plots (Figure 4 and Figure S1 in Kestemont et al.)
#------------------------------------------------------------
# Plot 5 classes of evenness measures
 
even = iNEXT.4steps::Evenness(abundance, nboot = 0,datatype = "abundance", q = seq(0,3,0.01), method = "Estimated",C = 1)
even = even[-1]
cbPalette <- rev(c("#D55E00", "#0072B2", "grey22", "purple1" , "red1", "orange1", 
                   "seagreen3", "royalblue4"))
classdata = cbind(do.call(rbind, even), class = rep(names(even), 
                                                    each = nrow(even[[1]])))
fig = ggplot(classdata, aes(x = Order.q, y = Evenness, colour = Assemblage)) + 
  geom_line(size = 1.05) + 
  scale_colour_manual(values = cbPalette) + scale_fill_manual(values = cbPalette) + 
  labs(x = "Order q", y = "Evenness") + theme_bw() + 
  theme(legend.position = "bottom", legend.box = "vertical", legend.key.width = unit(1.2, "cm"), 
        legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(-10, -10, -5, -10), 
        text = element_text(size = 16), 
        plot.margin = unit(c(5.5,5.5, 5.5, 5.5), "pt")) + 
  guides(linetype = guide_legend(keywidth = 2.5))+
  facet_wrap(~class)

fig

#-----------------------------------------------------------------------------
# Distribution of minimum required sample size (Figure 3C in Kestemont et al.)
#-----------------------------------------------------------------------------
# MSS.est is a function to obtain minimum sample size. 
# parameter y a vector of species sample frequencies. 
# parameter nboot a number specifying the number of bootstrap replications.

MSS.est <- function(y, nboot=1000){
  n <- sum(y)
  Abun.Mat <- rmultinom(nboot, n, iNEXT.3D:::EstiBootComm.Ind(y))
  Ms = n*(1+MS.x(y))
  pool = apply(Abun.Mat, 2, function(x) n*(1+MS.x(x)))
  pro.mean = mean(pool)
  pool = pro.mean - pool
  pool = Ms - pool 
  pool
}


minsample = MSS.est(unlist(abundance), nboot = 10000) %>% as.data.frame()


ggplot(minsample)+
  geom_histogram(aes(x= .,  y = ..density..), color="#e9ecef", bins = 30)+
  geom_density(aes(x = .), color = "red")+
  xlab("documents")+
  ylab("works density")+
  theme_bw()


