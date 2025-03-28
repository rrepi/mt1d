### Methylation score
#INPUT:
# Methylation M-values of differentially methylated CpGs (e.g., test_mvalues.txt)
# Metadata file including the subject_ID and covariables (e.g., test_metadata.txt)
# List of CpGs of interest (e.g. CpGs related to type 1 diabetes risk genes)
# (e.g., test_CpGs.txt)

#load packages
library(tidymodels)
library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(faux)
library(caret)
library(randomForest)

# set path of your data
setwd("")
set.seed(42)

#######################
# Prepare the data ----
#######################

# read the M-values of differentially methylated CpGs (columns=subject_ID, rownames=ProbeID)
m_values<- 
  read.csv("test_mvalues.txt",
           sep=""
           )

m_values <-
  rownames_to_column(m_values)

m_values<- 
  m_values %>% 
  dplyr::rename(cpg="rowname")

names(m_values)<-sapply(str_remove_all(colnames(m_values),"X"),"[") #in case sample_ID starts with X

# read the cpgs to be included (e.g. CpGs associated with T1D susceptibility genes)
cpgs <- 
  read.table("test_CpGs.txt", 
             header = TRUE) %>%
  filter(CpGT1Dgene == "Yes") %>%
  pull(CpG)

# filter methylation values of selected CpGs
m_values <-
  m_values %>%
  filter(cpg %in% cpgs)

# and transpose the matrix
m_values <-
  m_values %>%
  pivot_longer(!cpg, names_to = "subject_ID") %>%
  pivot_wider(names_from = cpg)

# read the metadata
metadata <- 
  read.table("test_metadata.txt", 
             header = TRUE)

# keep only CpGs that are uncorrelated (< 0.8)
cor_mat <- 
  cor(m_values[-1], use = "pairwise.complete.obs")
cols_to_keep <- 
  c()
for (i in 1:ncol(cor_mat)) {
  if (all(abs(cor_mat[i, -c(1:i)]) < 0.8)) {
    cols_to_keep <- c(cols_to_keep, i)
  }
}

m_values <- 
  m_values[c(1, cols_to_keep + 1)]

# merge the metadata and the data
data_full <-
  m_values %>%
  left_join(metadata, by = "subject_ID") %>%
  dplyr::select(all_of(names(metadata)), everything())


#####################################################
# Random forest selection to get predictive CpGs ----
#####################################################

# keep only complete observations of the included variables
data_reg <-
  data_full %>%
  dplyr::select(subject_ID, maternal_t1d, starts_with("cg")) %>%
  drop_na() %>%
  mutate(across(maternal_t1d, as.factor))

# Define the control using a random forest selection function
control <- 
  rfeControl(functions = rfFuncs, # random forest
             method = "repeatedcv", # repeated cv
             repeats = 5, # number of repeats
             number = 10) # number of folds


# Features
x <- 
  data_reg %>%
  dplyr::select(starts_with("cg")) %>%
  as.data.frame()

# Target variable
y <- 
  data_reg$maternal_t1d

# Training: 80%; Test: 20%
inTrain <- 
  createDataPartition(y, p = .80, list = FALSE)[,1]

x_train <- 
  x[ inTrain, ]

x_test  <- 
  x[-inTrain, ]

y_train <- 
  y[ inTrain]

y_test  <- 
  y[-inTrain]

# Run RFE
result_rfe1 <- 
  rfe(
    x = x_train, 
    y = y_train, 
    rfeControl = control)

# Print the results
result_rfe1

# Print the selected features
cpg_keep <-
  predictors(result_rfe1)

# Print the results visually
ggplot(
  data = result_rfe1, 
  metric = "Accuracy") + 
  theme_bw()

ggplot(data = result_rfe1,
       metric = "Kappa") +
  theme_bw()

varimp_data <- 
  data.frame(feature = row.names(varImp(result_rfe1))[1:5],
                          importance = varImp(result_rfe1)[1:5, 1])

ggplot(data = varimp_data, 
       aes(x = reorder(feature, -importance),
           y = importance, 
           fill = feature)) +
  geom_bar(stat = "identity") + 
  scale_fill_brewer(palette = "Paired") +
  geom_text(aes(label = round(importance, 1)),
            vjust = 1.6, 
            color = "white", 
            size = 4) +
  labs(x = "",
       y = "") +
  theme_bw() + 
  theme(legend.position = "none", axis.text = element_text(size = 18))

### use the predicted CpGs for the score
# merge the metadata and the data
m_values_pred <-
  m_values[,colnames(m_values) %in% cpg_keep]
m_values_pred$subject_ID <-
  m_values$subject_ID

data_full <-
  m_values_pred %>%
  left_join(metadata, by = "subject_ID") %>%
  dplyr::select(all_of(names(metadata)), everything())

data_reg <-
  data_full %>%
  dplyr::select(subject_ID, maternal_t1d, Study, starts_with("cg")) %>%
  drop_na() %>%
  mutate(across(maternal_t1d, as.factor))

# split the data
train_data <- data_reg %>%
  filter(Study=="0")
train_data <- train_data %>%
  dplyr::select(-Study)

test_data <- data_reg %>%
  filter(Study=="1")
test_data <- test_data %>%
  dplyr::select(-Study)

# create the classifier using the train data
the_recipe <-
  recipe(maternal_t1d ~ ., 
         data = train_data) %>%
  update_role(subject_ID, 
              new_role = "ID")

lr_mod <-
  logistic_reg() %>%
  set_engine("glm")

the_workflow <-
  workflow() %>%
  add_model(lr_mod) %>%
  add_recipe(the_recipe)

the_fit <-
  the_workflow %>%
  fit(data = train_data)

the_fit %>%
  extract_fit_parsnip() %>%
  tidy() %>%
  print(n = Inf)

# ROC curve using the train data
augment(the_fit, train_data) %>%
  roc_curve(truth = maternal_t1d, .pred_No) %>%
  autoplot()

augment(the_fit, train_data) %>%
  roc_auc(truth = maternal_t1d, .pred_No)

# check the classifier using the test data
predict(the_fit, test_data)

the_aug <-
  augment(the_fit, test_data)

the_aug %>%
  dplyr::select(maternal_t1d, subject_ID, .pred_class, .pred_Yes)

the_aug %>%
  roc_curve(truth = maternal_t1d, .pred_No) %>%
  autoplot()

the_aug %>%
  roc_auc(truth = maternal_t1d, .pred_No)

# calculate the scores for all samples
intercept <-
  the_fit %>%
  extract_fit_parsnip() %>%
  tidy() %>%
  dplyr::select(term, estimate) %>%
  dplyr::slice(1) %>%
  pull(estimate)

cpg_weights <-
  the_fit %>%
  extract_fit_parsnip() %>%
  tidy() %>%
  dplyr::select(term, estimate) %>%
  dplyr::slice(-1)

# make sure the column order and the weights order is the same
cpgs_matrix <-
  dplyr::select(data_full, starts_with("cg")) %>%
  as.matrix()

cpg_vector <-
  cpg_weights[
    terms = names(cpgs_matrix),
  ] %>%
  pull(estimate)

# matrix multiplication of the weights by the CpGs values
samples_scores <-
  cpgs_matrix %*% cpg_vector

# add methylation score to the samples
data_full$risk_score <- samples_scores

# get calculated weights per CpG integrated in the methylation score
cpg_weights <- cpg_weights %>% 
  dplyr::rename(CpG = term, weight = estimate)

# end
