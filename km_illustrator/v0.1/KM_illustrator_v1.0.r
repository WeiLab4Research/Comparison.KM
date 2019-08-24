

# check packages and load them; define theme
init = function(journal_theme = "jco"){
  if(!require(pacman)) #install.packages("pacman")
  
  #pacman::p_load(survival, survminer, ggplot2, gridExtra, customLayout, magrittr,
  #               magrittr, ggpubr, ggsci, KMsurv, survMisc, TSHRC, surv2sample)
  
  dyn.load("d:/surv2sample.dll")
  
  # default journal theme, global scope variable
  C_JOURNAL <<- journal_theme
  journal_theme <<- journal_theme
  if(tolower(journal_theme) == "nejm"){
    C_COLOR <<- scale_color_nejm()
    C_FILL <<- scale_fill_nejm()
  }else if(tolower(journal_theme) == "lancet"){
    C_COLOR <<- scale_color_lancet()
    C_FILL <<- scale_fill_lancet()
  }else if(tolower(journal_theme) == "jama"){
    C_COLOR <<- scale_color_jama()
    C_FILL <<- scale_fill_jama()
  }else if(tolower(journal_theme) == "jco"){
    C_COLOR <<- scale_color_jco()
    C_FILL <<- scale_fill_jco()
  }else{
    C_COLOR <<- scale_color_aaas()
    C_FILL <<- scale_fill_aaas()
  }
}

 km_compare = function(fit, ds){
   # I believe ds can be retrived from the fitted model, which is not necessary to be an independent parameter
   
   terms = terms.formula(fit$call$formula, keep.order = TRUE)
   surv = terms[[2]]
   
   time_var = as.character(surv)[2]
   event_var = as.character(surv)[3]
   by_var = as.character(terms)[3]
   
   time = ds[, time_var]
   event = ds[, event_var]
   by = ds[, by_var]
   by_num = match(by, names(table(by))) - 1
   
   time[time <= 0] = 0.1
   
   # log ranks and renyi
   rlt_logranks_renyi = ten(fit)
   comp(rlt_logranks_renyi, p = c(1, 0, 1), q = c(0, 1, 1))
   
   rlt_logranks = as.data.frame(attributes(rlt_logranks_renyi)$lrt)

   # two stage, KS and Neyman only fit to 2 groups   
   if(length(table(by_num)) == 2){
     rlt_renyi = as.data.frame(attributes(rlt_logranks_renyi)$sup)

     # two stage
     rlt_ts = twostage(time = time, delta = event, group = by_num, 
                      nboot = 2000, alpha = 0.05, eps = 0.1)

     # Kolmogorov-Smirnov test
     # the group coding only can be 1 or 2 (actually very wiered), (by_num + 1) is adaption to this
     res_ks = surv2.ks(Surv(time, event), group = by_num + 1,
                       process = "w", approx = "boot", 
                       nsim = 2000, nsim.plot = 500)
     
     if(res_ks$pval.ad == 0){
       cat("\n Too significant to get the p for AD and CVM by 2000 simulations...\n")
       cat("\n We are increasing the nsim for AD and CVM methods to obtain the non zero pvalue...\n")
       res_ks = surv2.ks(Surv(time, event), group = by_num + 1,
                         process = "w", approx = "boot", 
                         nsim = 1000000, nsim.plot = 500)
     }
     
     if(res_ks$pval.ad == 0){
       cat("\n Too significant to get the p for AD and CVM by 1000000 simulations. P was set to 0.000001..\n")
       res_ks$pval.ad  = 1/1000000
       res_ks$pval.cm  = 1/1000000
     }
   
    # Neyman test
     res_neyman = surv2.neyman(Surv(time, event), group = by_num + 1,
                               data.driven = TRUE)
   }

   ## construct p table
   pvalues = c()
   methods = c("LR","GW","TW","PP","mPP","HF(1,0)","HF(0,1)","HF(1,1)",
               "LR*","GW*","TW*","PP*","PPb*","HF(1,0)*","HF(0,1)*","HF(1,1)*",
               "TS","KS","CVM","AD","NM")

   pvalues = unlist(rlt_logranks[, ncol(rlt_logranks)])
   
   if(length(table(by)) == 2){
      pvalues = c(pvalues, 
                  unlist(rlt_renyi[, ncol(rlt_renyi)]),
                  rlt_ts[3],
                  as.numeric(res_ks$pval.ks.asympt), 
                  as.numeric(res_ks$pval.cm),
                  as.numeric(res_ks$pval.ad),
                  as.numeric(res_neyman$pval.asympt))
   }
   
   p_csv = data.frame(p_value = pvalues, Methods = methods[1:length(pvalues)])
   ptable = data.frame(p_value = pvalues, Methods = methods[1:length(pvalues)])   
   ptable$group = "A_Original"
   
   ptable$p_value_log10 = -log10(ptable$p_value)
   
   ptable$p_value_log10 = ifelse(is.infinite(ptable$p_value_log10), NA, ptable$p_value_log10)
   
   ptable$group = ifelse(ptable$p_value > 0.05, "insignificance", ptable$group)
   
   # define p labels
   plabels = c(0, 2, 4)
   maxp = max(ptable$p_value_log10, na.rm = TRUE)
   if(maxp > 4){
     plabels = c(plabels, seq(from = 4, to = ceiling(maxp) - 1, length.out = 2))
   }
   
   ptable$Significance = ptable$p_value_log10
   
   # P-value plot
   pplot = ggdotchart(ptable,  
                    x = "Methods", 
                    y = "p_value_log10",
                    main = "Kaplan-Meier\nComparison",
                    sorting = "descending",
                    rotate = TRUE,
                    font.main = 14,
                    add = "segments",
                    xlab = "Methods",
                    ylab = expression(-log[10](italic(P))),
                    font.xtickslab = c(12, "black"),
                    font.ytickslab = c(10, "black"),
                    font.x = c(14, "black"),
                    font.y = c(14, "black"),
                    y.text.angle = 0,
                    color = "group",
                    sort.by.groups = FALSE,
                    position = position_dodge(0.7)) + 
     labs(color = "", hjust = 1, angle = 90) +
     geom_hline(yintercept = -log10(0.05), color = "grey60", linetype = "dashed", size = 0.5) + 
     scale_y_continuous(breaks = plabels) + 
     guides(color = FALSE) + 
     scale_size(guide = FALSE) +
     geom_point(aes(x = Methods, size = Significance, y = p_value_log10),
                 data = subset(ptable, group == "insignificance"), color = "lightgrey") +
     geom_point(aes(x = Methods, size = Significance,color = group, y = p_value_log10), 
                data = subset(ptable, group == "A_Original")) +
     C_COLOR + C_FILL
   
    return(list(ptable = ptable, p_csv = p_csv, pplot = pplot))
 }
 
 
km_plot = function(fit, x_step = 24, risk_table = TRUE,
                  conf = TRUE, xlab = "Overall Survival (months)",
                  by_labels = unique(by)
                  ){

   
  ggsurv = ggsurvplot(fit, size = 1, # change line size
                       conf.int = conf, # Add confidence interval
                       xlab = xlab,
                       log.rank.weights = "1",
                       font.legend = 14,
                       break.x.by = ifelse(x_step >= median(fit$time), ceiling(quantile(fit$time, probs = 0.5)), x_step),
                       legend.labs = by_labels[1:length(fit$strata)],
                       font.xtickslab = c(12,"black"),
                       font.ytickslab = c(12,"black"),
                       font.x = c(14,"black"),
                       font.y = c(14,"black"),
                       legend.= "",
                       palette = journal_theme,
                       risk.table.y.text.col = TRUE,
                       risk.table = risk_table, # Add risk table
                       risk.table.fontsize = 4.5)
   
   ggsurv$plot = ggsurv$plot
   
   if(risk_table){
     ggsurv$table = ggsurv$table +  
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14,angle = 90))
   }
   
   return(ggsurv)
 }

km_illustrator = function(km_compare, km_plot){
   
  lay  <- lay_new(matrix(1:2, ncol = 1), heights = c(2.5, 1))
  lay2 <- lay_new(matrix(1:1))
  cl   <- lay_bind_col(lay, lay2, widths = c(2, 1))
  
  p = lay_grid(list(km_plot$plot, km_plot$table, km_compare$pplot), cl) 
  
  return(p)
} 


## example
# check packages and set journal theme
init(journal_theme = "jco")

# load testing dataset
ds = survminer::BRCAOV.survInfo
set.seed(20190822)
ds = ds[sample(1:nrow(ds), replace = TRUE, size = 200), ]
ds$times = ds$times/365.25

#ds[1:20, "group"] = 3
fit = survfit(Surv(times, patient.vital_status) ~ admin.disease_code, data = ds) 
km_pval = km_compare(fit, ds)
km_fig = km_plot(fit, x_step = 2, conf = TRUE,xlab = "OS (year)", by_labels = c("BRCA","OV"))
km_illustrator(km_pval, km_fig)


