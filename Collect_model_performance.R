
### Requited Packages####

required_packages <- c(
  "caret",           # Main ML framework
  "pls",             # PLS method
  "randomForest",    # Random Forest method
  "kernlab",         # SVM method
  "glmnet",          # GLM method
  "nnet",            # Neural Network method
  "doParallel",      # Parallel processing
  "foreach",          # Parallel loops
  "gmodels",
  "ggplot2",
  "reshape2",
  "rstatix",
  "ggpubr",
  "dplyr"
)

# Install missing packages
new.packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
lapply(required_packages, library, character.only = TRUE)

## Initial set up

wdir <- getwd()
gene.list =  c( "BBK32","OspC","DbpA","RevA","P66" )
sslist = c("fs055", "ImpQ1", "ImpQ2", "ImpQ3", "VarQ1", "VarQ2", "VarQ3")
models =c('PLS', "RF", "SVM", 'GLM',"PCANN")

inlist = ls()

inlist = ls()

#### Model Performance  #####
for (gene in gene.listS) {
  data_path <- file.path(wdir, "data", gene)
  results_path <- file.path(wdir, "results", gene)
  setwd(results_path)
  
  for (ss in c(sslist)) {
    load(paste0(ss, "_results.RData"))
    
    ####
    metrics = colnames(get(paste0(ss, "_", models[1], "_O_LOOfit"))[[1]]$Test.Perf)
    for (m in seq(along = models)) {
      LOOTrain_Perf = data.frame(matrix(NA, nrow = 100, ncol = 11))
      LOOTest_Perf = data.frame(matrix(NA, nrow = 100, ncol = 11))
      bsTrain_Perf = data.frame(matrix(NA, nrow = 100, ncol = 11))
      bsTest_Perf = data.frame(matrix(NA, nrow = 100, ncol = 11))
      for (i in seq(100)) {
        LOOTrain_Perf[i, ] = get(paste0(ss, "_", models[m], "_O_LOOfit"))[[i]][["Train.Perf"]][1:11]
        LOOTest_Perf[i, ] = t(get(paste0(ss, "_", models[m], "_O_LOOfit"))[[i]][["Test.Perf"]])
        bsTrain_Perf[i, ] = get(paste0(ss, "_", models[m], "_O_bsfit"))[[i]][["Train.Perf"]][1:11]
        bsTest_Perf[i, ] = t(get(paste0(ss, "_", models[m], "_O_bsfit"))[[i]][["Test.Perf"]])
      }
      colnames(LOOTrain_Perf) = colnames(LOOTest_Perf) = colnames(bsTrain_Perf) =
        colnames(bsTest_Perf) = metrics
      Perf = cbind(LOOTrain_Perf, LOOTest_Perf, bsTrain_Perf, bsTest_Perf)
      colnames(Perf) = c(
        paste0("LOOTrain.", metrics),
        paste0("LOOTest.", metrics),
        paste0("bsTrain.", metrics),
        paste0("bsTest.", metrics)
      )
      assign(paste0(ss, "_", models[m], "_LOOTrainPerf"), LOOTrain_Perf)
      assign(paste0(ss, "_", models[m], "_LOOTestPerf"), LOOTest_Perf)
      assign(paste0(ss, "_", models[m], "_bsTrainPerf"), bsTrain_Perf)
      assign(paste0(ss, "_", models[m], "_bsTestPerf"), bsTest_Perf)
      df = sapply(Perf, ci)
      rownames(df) = paste0(models[m], "_", rownames(df))
      Perf$model = models[m]
      if (m == 1) {
        ciPerf = df
        fullPerf = Perf
      } else {
        ciPerf = rbind(ciPerf, df)
        fullPerf = rbind(fullPerf, Perf)
      }
      rm (LOOTrain_Perf,
          LOOTest_Perf,
          bsTrain_Perf,
          bsTest_Perf,
          df,
          Perf)
    }
    write.csv(fullPerf, paste(gene, ss, "fullPerf.csv", sep = "_"))
    write.csv(ciPerf, paste(gene, ss, "ciPerf.csv", sep = "_"))
    m.ciPerf = round(ciPerf, digits = 3)
    sumciPerf = data.frame(matrix(NA, nrow = 5, ncol = 20))
    
    col.idx = c(3, 1, 4, 5, 10,
                c(3, 1, 4, 5, 10) + 11,
                c(3, 1, 4, 5, 10) + 22,
                c(3, 1, 4, 5, 10) + 33)
    m.ciPerf = m.ciPerf[, col.idx]
    for (i in seq(5)) {
      for (j in seq(20)) {
        sumciPerf[i, j] = paste0(m.ciPerf[1 + 4 * (i - 1), j], " (", m.ciPerf[2 +
                                                                                4 * (i - 1), j], " - ", m.ciPerf[3 + 4 * (i - 1), j], ")")
      }
    }
    rownames(sumciPerf) = paste(gene, ss, models, sep = "_")
    colnames(sumciPerf) = colnames(m.ciPerf)
    sumciPerf
    write.csv(sumciPerf, paste(gene, ss, "sumciPerf", sep = "_"))
    if (ss == sslist[1]) {
      dfPerf = sumciPerf
    } else {
      dfPerf = rbind(dfPerf, sumciPerf)
    }
   
    ##### Compare ROC value ####
    oj = c("LOOTrain" , "LOOTest", "bsTrain", "bsTest")
    for (o in oj) {
      Perf_ROC = data.frame(matrix(NA, nrow = 100, ncol = 5))
      colnames(Perf_ROC) = models
      for (i in seq(100)) {
        for (j in seq(along = models)) {
          Perf_ROC[i, j] = get(paste0(ss, "_", models[j], "_", o, "Perf"))[i, 3]
        }
      }
      test_ROC = melt(Perf_ROC, variable.name = "model")
      test_ROC_stat <- test_ROC %>%
        t_test(value ~ model, p.adjust.method = "bonferroni")
      bp <- ggplot(test_ROC, aes(x = model, y = value)) +
        scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, by = 0.2)) +
        stat_boxplot(
          geom = "errorbar",
          # Boxplot with error bars
          width = 0.2
        ) +
        geom_boxplot(
          aes(fill = model),
          colour = "#1F3552",
          # Colors
          alpha = 0.9 #, outlier.colour = "red"
        ) +
        geom_jitter(height = 0, size = rel(0.3)) +
        labs(y = "AUC-ROC",
             title = paste0(gene, "_", ss, "_", o),
             x = NULL) +
        
        theme(
          axis.line = element_line(
            colour = "black",
            # Theme customization
            linewidth = 0.25
          ),
          legend.position = "none",
          plot.title = element_text(size = rel(1.2)),
          # plot.background = element_rect(fill="white"),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(
            colour = "gray",
            # Theme customization
            linewidth = 0.15
          ),
          axis.text.y  = element_text(size = rel(1), color = "black"),
          axis.text.x = element_text(size = rel(1), color = "black"),
          axis.title = element_text(size = rel(1))
        )
      
      test_ROC_stat  <-  test_ROC_stat %>% add_xy_position(x = "model")
      bp <- bp + stat_pvalue_manual(
        test_ROC_stat ,
        vjust = 0,
        label = "p.adj.signif",
        bracket.nudge.y = 0.01,
        #Vertical adjustment to nudge brackets by.positive value, brackets will be moved up; if negative value, brackets are moved down.
        step.increase = 0.04,
        remove.bracket = F,
        tip.length = 0.02,
        hide.ns = T
      )
      test_ROC_stat$set = o
      if (o == oj[1])
      {
        paired_ROC = test_ROC_stat
      } else {
        paired_ROC = rbind(paired_ROC, test_ROC_stat)
      }
      assign(paste0(o, "Perf_ROC"), Perf_ROC)
      assign(paste0(o, "_ROC_bp"), bp)
      rm(test_ROC_stat, test_ROC, bp, Perf_ROC)
    }
    write.csv(sapply(paired_ROC, as.character),
              paste(gene, ss, "paired_ROC.csv", sep = "_"))
    
    gplot = ggarrange(
      LOOTrain_ROC_bp ,
      LOOTest_ROC_bp ,
      bsTrain_ROC_bp  ,
      bsTest_ROC_bp,
      ncol = 1,
      nrow = 4
    )
    assign(paste0(ss, "_gplot"), gplot)
    #Optional: Save plot
    # ggsave(
    #   filename = paste(gene,ss, "ROC.png", sep = "_"),
    #   plot = gplot,
    #   width = 8,
    #   height = 12,
    #   dpi = 300
    # )
    
    ##### VIP scores ####
    for (m in seq(5)){
      for (i in seq(100)) {
        bst = get(paste0(ss,"_",models[m],"_O_bsfit"))[[i]]$VIP_score[1]
        LOOt =  get(paste0(ss,"_",models[m],"_O_LOOfit"))[[i]]$VIP_score[1]
        LOOt$Feature = bst$Feature = rownames(bst)
        colnames(LOOt)[1] = colnames(bst)[1]=paste0(models[m],i)
        if (i == 1) {
          bsscore <- bst
          LOOscore <- LOOt
        } 
        else {
          bsscore <- merge(bsscore,bst, by ="Feature",sort = FALSE, all=TRUE)
          LOOscore <- merge(LOOscore,LOOt, by ="Feature",sort = FALSE, all=TRUE)
        }
      }
      #colnames(score) <- c("Feature",paste0(models[j],seq(100))) 
      if (m == 1) {
        bsVIP_scores <- bsscore
        LOOVIP_scores <- LOOscore
      } else {
        bsVIP_scores <- merge(bsVIP_scores,bsscore,by ="Feature",sort = FALSE, all=TRUE) 
        LOOVIP_scores <- merge(LOOVIP_scores,LOOscore,by ="Feature",sort = FALSE, all=TRUE) 
      }
      rm(bst,LOOt,bsscore, LOOscore)
    }
    
    write.csv(bsVIP_scores, paste(gene,ss,"bsVIP_Scores.csv", sep ="_"))
    write.csv(LOOVIP_scores, paste(gene,ss,"LOOVIP_Scores.csv", sep ="_"))
    
    
    #### Top20 features #####
    for (str in c("bs","LOO")) {
      Top20_var = foreach( i =2:501) %do% {
        df = get(paste0(str,"VIP_scores"))[c(1,i)]
        df %>% top_n(20)
      }
      Top20 =  Top20_var  
      for (i in seq(500)){ 
        Top20_var[[i]] =Top20_var[[i]][1]
      }   
      PLS = data.frame(table(unlist(Top20_var[1:100])))
      RF = data.frame(table(unlist(Top20_var[101:200])))
      SVM = data.frame(table(unlist(Top20_var[201:300])))
      GLM = data.frame(table(unlist(Top20_var[301:400])))
      PCANN = data.frame(table(unlist(Top20_var[401:500])))
      ImpTop20 = merge( PLS, RF,by ="Var1", sort=FALSE, all= TRUE)
      ImpTop20= merge( ImpTop20, SVM,by ="Var1", sort=FALSE, all= TRUE)
      ImpTop20= merge( ImpTop20, GLM,by ="Var1", sort=FALSE, all= TRUE)
      ImpTop20 = merge( ImpTop20, PCANN,by ="Var1", sort=FALSE, all= TRUE)
      colnames(ImpTop20) = c("Feature", paste0("Freq_", models) )
      # assign(paste0(str,"_Top20"),  Top20)
      write.csv(ImpTop20, paste(gene,ss, str,"ImpTop20.csv",sep="_"))
      rm( PLS, RF, SVM, PCANN, Top20_var)
    } 
    
    ###################################  
  }
  ### For Figures S1-S5 #####
  ROCplot = ggarrange(
    fs055_gplot,
    ImpQ1_gplot,
    ImpQ2_gplot,
    ImpQ3_gplot,
    VarQ1_gplot,
    VarQ2_gplot,
    VarQ3_gplot,
    ncol = 7,
    nrow = 1
  )
  print(ROCplot)
  ggsave(
    filename = paste(gene, "ROC.png", sep = "_"),
    plot = ROCplot,
    width = 45,
    height = 15,
    dpi = 300
  )
  
  ### For Table S3 ####
  if (gene == gene.list[1]) {
    sumPerf = dfPerf
  } else {
    sumPerf = rbind(sumPerf, dfPerf)
  }
}


write.csv(sumPerf, paste0(wdir, "/results/", "sumPerf.csv")). ###Table S3
setwd(wdir)





