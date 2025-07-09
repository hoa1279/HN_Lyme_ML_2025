###Data setup ####
cwd <- getwd()

gene.list = c( "OspC","DbpA","BBK32","RevA","P66" )
for (g in seq(7)) {
gene = gene.list[g]
setwd(paste0(cwd,"/data/",gene))


data <- readxl::read_excel(paste0(gene,"_unq_seqs.xlsx"),sheet="use")
# Identify and remove NA dissemination records
na_indices <- which(data$Dissemination == "NA")
data <- if(length(na_indices) > 0) data[-na_indices, ] else data

ids <- data$seq_name
Y <- factor(data$Dissemination, level=c("YES","NO"))
OHE <- read.csv(paste0(gene,"_OHE.csv"))
row.names(OHE) <- OHE$X
X_OHE <- OHE[ids,-1]

# Model independent VarImP
O_roc_imp <- filterVarImp(X_OHE,Y) 
summary(O_roc_imp)
O_imp_vars <- as_tibble(O_roc_imp,rownames='variables') %>% 
  arrange(desc(YES)) %>% 
  filter(YES>=0.55) %>% 
  pull(variables)
fs_X_OHE <- X_OHE %>% select( all_of(O_imp_vars))

quantile(O_roc_imp$YES)
O_ImpQ1 <-  as_tibble(O_roc_imp,rownames='variables') %>% 
  arrange(desc(YES)) %>% 
  filter(YES>=quantile(YES)[2])  %>% # dbpA use >
  pull(variables)
X_O_ImpQ1 = X_OHE %>% select( all_of(O_ImpQ1))
O_ImpQ2 <-  as_tibble(O_roc_imp,rownames='variables') %>% 
  arrange(desc(YES)) %>% 
  filter(YES>quantile(O_roc_imp$YES)[3]) %>% 
  pull(variables)
X_O_ImpQ2 = X_OHE %>% select( all_of(O_ImpQ2))
O_ImpQ3 <-  as_tibble(O_roc_imp,rownames='variables') %>% 
  arrange(desc(YES)) %>% 
  filter(YES>=quantile(O_roc_imp$YES)[4]) %>% 
  pull(variables)
X_O_ImpQ3 = X_OHE %>% select( all_of(O_ImpQ3))

#feature filter by variation
O_vars = apply(X_OHE,2,var)
quantile(O_vars)
O_VarQ1 = X_OHE[,O_vars >=quantile(O_vars)[2]]
O_VarQ2 = X_OHE[,O_vars >quantile(O_vars)[3]]
O_VarQ3 = X_OHE[,O_vars >=quantile(O_vars)[4]]
save.image(paste0(cwd,"/data/",gene,"/",
                  gene,"_fs_ML_input_data.RData"))
}