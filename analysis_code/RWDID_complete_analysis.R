# =============================================================================
# RWDID COMPLETE ANALYSIS CODE
# Single consolidated source file for the RWDID research pipeline
# =============================================================================
#
# PURPOSE
#   This single file brings together the supplied RWDID analysis source code:
#   dietary item preparation, frozen 70/30 cohort split, Model 2 rotating
#   lowest-risk references, final-15 item/portion audit, RCS20 interval
#   derivation, final binary scoring, continuous RWDID-C scoring, training-only
#   beta estimation, and held-out test validation.
#
# WEB CALCULATOR VS. RESEARCH ANALYSIS
#   - The GitHub Pages calculator is a static browser application. It applies
#     frozen score rules and coefficients; it does not fit Cox/RCS models.
#   - This file is the research-analysis source code used to derive and validate
#     the score. The machine-readable parameters used by the website are
#     separately published under ./data/.
#
# EXECUTION ORDER
#   01. Step_02_1_1: dietary food-group preparation.
#   02. Step_03_1_0_SPLIT: frozen 70/30 split manifest and cohort RDS files.
#   03. Step_03_1_0_MODEL2_REFERENCE: Model 2 rotating lowest-risk reference.
#   04. Step_03_2_7_FINAL15_AUDIT: original Table S1 item/portion audit.
#   05. Step_04_1_1_1_RCS_BASE: initial RCS analysis and anchor outputs.
#   06. Step_04_1_1_2_INTERVAL_RECONCILIATION: data-derived interval audit.
#   07. Step_04_1_1_3_RCS20_OVERLAP: frozen 5%/10%/20% overlap variants.
#   08. Step_04_1_1_4_FINAL_RCS_BINARY: final RCS diagnostics and binary score.
#   09. Step_04_1_2_RWDID_C: continuous equal-weight RWDID-C (0-15).
#   10. Step_04_2_1_BETA_WEIGHTED_TRAINING: training-only binary beta weights
#       and beta-weighted RWDID (0-1).
#   11. Step_04_2_2_HELDOUT_VALIDATION: frozen-parameter test-set validation.
#
# RUNNING THE CONSOLIDATED FILE
#   - Run from the UK Biobank project root, not from this GitHub Pages folder.
#   - The source blocks retain their own rm(list = ls()) commands and original
#     file-based hand-offs. Each stage reloads the required upstream files.
#   - Execute this file from top to bottom. Do not run a later section before
#     its required upstream outputs have been generated.
#
# REQUIRED AUTHORISED PROJECT INPUTS
#   - UK Biobank individual-level source files and the project-specific
#     preprocessed objects:
#       ./data/ukb_data_dietary_org.RData
#       ./data/ukb_data_outcomes.RData
#       ./data/ukb_data_cancer_outcomes.RData
#       ./data/ukb_data_covariates.RData
#   - Dietary metadata / Table S1 mapping audit files:
#       ./data/Step_01_1_file_id_and_description_dietary.csv
#       ./data/Step_01_1_file_id_and_description_dietary_portion_data.csv
#       ./data/field.csv
#   - UKB access and project approvals are required. Individual-level UKB data
#     and participant-level score files are not included in this public site.
#
# PATH NORMALISATION
#   The 2026-11-23 split-manifest source originally wrote to
#   ./results/ukb_70_30_split. In this consolidated copy, only SPLIT_DIR is
#   normalised to ./results/Step_03/ukb_70_30_split so it matches the exact
#   downstream paths used by the supplied Step_03 / Step_04 analysis scripts.
#   The split seed, 70/30 allocation, ID handling, and split algorithm are
#   unchanged.
#
# FROZEN ANALYSIS DESIGN
#   - 15 final dietary components; no green leafy vegetables or high-fat dairy.
#   - Binary RWDID: 15 frozen 0/1 component scores weighted by the 15
#     mutually adjusted Model 3 Cox beta coefficients estimated in training.
#   - Continuous RWDID-C: 15 proportional component scores (0-1), equally
#     weighted and summed to 0-15.
#   - The held-out test cohort does not re-estimate intervals, beta weights,
#     tertile cutpoints, or per-20% scaling parameters.
#   - Missing component intake is not imputed.
#
# SOURCE INTEGRITY
#   The source blocks below are delimited and labelled with their original
#   filenames. Except for the explicitly documented SPLIT_DIR normalisation,
#   source code is included as supplied. The consolidated file does not
#   silently change model covariates, scoring intervals, beta estimates, or
#   statistical specifications.
#
# =============================================================================



########################################################################################
# BEGIN ANALYSIS BLOCK 01/11: Step_02_1_1
# Source: Step_02_1_1_data_preparation_and_wrangling_food_group(5).R
# Provenance: local current source
# SHA256 of included block: 37b877ce5f0d4bc22e450ad70a6192411c61bd76f3ca8c720374e5b345a3972d
########################################################################################

rm(list = ls())
## -------------------- Step 0.1 Loading Packages. --------------------
library(tidyverse)
library(survival)
library(ggsci)
library(table1)
library(survminer)
library(MatchIt)
library(gtsummary)
library(jstable)
library(forestplot)
library(lubridate)
library(ggplot2)
library(cols4all)
options(max.print = 999999)

## -------------------- Step 1 Data Construction and Cleaning --------------------
load("./data/ukb_data_dietary_org.RData")
colnames(ukb_data)

## -------------------- Step 2 数据质控与排除 --------------------
### -------------------- Step 2.1 with completion period less than the 5th percentile --------------------
# we invalidated questionnaires with completion period less than the 5th percentile (Field ID 20082) 
head(ukb_data %>% select(starts_with('p20082')))
table(is.na(ukb_data$p20082_i0))

# 生成回答问题的平均时间，因后续需要求平均值，此处也求平均值p20082_mean
# 删除p20082_mean为NA的行
# 删除less than the 5th percentile的人群
ukb_data <- ukb_data %>%
  # 选出以 p20082_ 开头的列
  mutate(p20082_mean = rowMeans(select(., starts_with("p20082_")), na.rm = TRUE)) %>%
  
  # 删除均值为 NA 的行
  filter(!is.na(p20082_mean)) %>% # 210839，数量对的上
  
  # 计算第5百分位数并过滤
  filter(p20082_mean > quantile(p20082_mean, probs = 0.05)) #197640


### -------------------- Step 2.2 who self-reported weekly changes in daily diet --------------------
# who self-reported weekly changes in daily diet (Field ID 1548)
# Data-Field 1548
# Description:	Variation in diet
# 1	Never/rarely
# 2	Sometimes
# 3	Often
# -1	Do not know
# -3	Prefer not to answer

head(ukb_data %>% select(starts_with('p1548')))

ukb_data %>% count(p1548_i0)

# 保留p1548_i0为1或2的行
ukb_data <- ukb_data %>%
  filter(p1548_i0 %in% c(1, 2)) # 排除15752，数量对的上

### -------------------- Step 2.3 who had doubtful estimates of total energy intake --------------------
# who had doubtful estimates of total energy intake (<500 or >3,500 kcal per day for women and <800 or >4,200 kcal per day for men, Field ID 26002)
head(ukb_data %>% select(starts_with('p26002')))
head(ukb_data %>% select(starts_with('p31')))
# 0	Female
# 1	Male

# 先生成均值
ukb_data <- ukb_data %>%
  # 选出以 p20082_ 开头的列
  mutate(energy_intake_mean = rowMeans(select(., starts_with("p26002_")), na.rm = TRUE))

# 再根据性别和均值进行过滤
ukb_data <- ukb_data %>%
  filter(
    (p31 == 0 & (energy_intake_mean < 500 | energy_intake_mean > 3500)) |
      (p31 == 1 & (energy_intake_mean < 800 | energy_intake_mean > 4200))
  ) # 181888 - 180484 排除1404，数量对的上

# 数据质控完成，剩下需要排除基线患病的人群

## -------------------- Step 3 回答问题的基线时间数据 --------------------

# 110001	Invitation to complete online 24-hour recall dietary questionnaire, acceptance
# 0	No response
# 1	Partial
# 2	Completed

# 110002	Invitation to complete online 24-hour recall dietary questionnaire, date sent

# i0的时间就是基线时间
head(ukb_data %>% select(starts_with('p53')))

head(ukb_data %>% select(starts_with(c('p110002', 'p110001'))))

# 时间为p53_i0和starts_with("p110002_")的列
ukb_data_date <- ukb_data %>%
  # 1. 先把 p53_i0 和所有 p110002_ 列转换成 Date 格式
  mutate(across(c(p53_i0, starts_with("p110002_")), ~ as.Date(.))) %>%
  
  # 2. 宽表转长表
  pivot_longer(
    cols = c(p53_i0, starts_with("p110002_")),
    names_to = "var",
    values_to = "date"
  ) %>%
  
  # 3. 按 eid 分组取最大日期
  group_by(eid) %>%
  summarise(
    date_food_baseline = if (all(is.na(date))) NA_Date_ else max(date, na.rm = TRUE),
    .groups = "drop"
  )


# 合并结果回原始数据框
ukb_data <- ukb_data %>%
  left_join(ukb_data_date, by = "eid")

# ukb_data <- ukb_data %>%
#   mutate(p110001_max = do.call(pmax, c(select(., starts_with("p110001")), na.rm = TRUE)))
# 
# table(ukb_data$p110001_max)
# 
# 
# ukb_data_test <- ukb_data %>%
#   filter(p110001_max == 0)
# 
# head(ukb_data_test %>% select(starts_with(c('p110002', 'p110001'))))

## -------------------- Step 4 疾病诊断时间，用于排除基线患病人群 --------------------
# load('./data/ukb_data_outcomes.RData')

# 该步骤放后面做


## -------------------- Step 5 Food Groups 前期数据整理与函数设定 --------------------
### -------------------- Step 5.1 分组处理每组的变量，连续型不用替换，而80和其他的需要替换 --------------------
# Group数据读取
dietary_id <- read.csv("./data/Step_01_1_file_id_and_description_dietary.csv")
dietary_portion_id <- read.csv("./data/Step_01_1_file_id_and_description_dietary_portion_data.csv")

# 试剂上是202种
length(unique(dietary_id$field_id))
# [1] 202

# 去重，用于后面的合并
dietary_portion_id <- dietary_portion_id %>% 
  distinct(field_id, .keep_all = TRUE)

# 查看哪些重复
dietary_portion_id %>%
  count(field_id) %>%
  filter(n >= 2)

# field_id n
# 1   100950 2
# 2   101020 2
# 3   101090 2
# 4   101160 2

### 数据清洗，合并与整理
dietary_id <- dietary_id %>%
  filter(food_group != "") %>%
  filter(food_group != "FFQ") %>%
  filter(food_group != "date")

colnames(dietary_id)
dietary_id <- dietary_id %>%
  select(field_id, description, category, food_group)

# 用于判断替换的编码，还有是否为连续型
field_data <- read.csv("./data/field.csv")
field_data <- field_data %>%
  select(field_id, encoding_id)

# 合并数据
dietary_id <- dietary_id %>%
  left_join(field_data)


unique(dietary_id$encoding_id)

# 第一种情况的替换
dietary_id_1 <- dietary_id %>%
  filter(!encoding_id %in% c(0, 80))

field_ids <- dietary_id_1$field_id

# 生成 pXXXXXX 前缀向量
p_cols <- paste0("p", field_ids)

# 第一步：值替换
ukb_data <- ukb_data %>%
  mutate(
    across(
      .cols = starts_with(p_cols),
      .fns = ~ case_when(
        . == 444 ~ 0.25,
        . == 555 ~ 0.5,
        . == 200 ~ 2,
        . == 300 ~ 3,
        . == 400 ~ 4,
        . == 500 ~ 5,
        . == 600 ~ 6,
        TRUE ~ .
      )
    )
  )


# 第二情况的替换，category == 80
ukb_data %>% select(starts_with(c('p20090')))
# 这部分进用来判断是否食用，适用于单个食用油的种类
ukb_data <- ukb_data %>%
  mutate(across(
    starts_with("p20090"),
    ~ ifelse(
      is.na(.),
      NA_real_,
      str_count(., ",") + 1
    )
  ))

# 第三种情况的替换，category == 0的连续型不用替换
dietary_id_3 <- dietary_id %>%
  filter(encoding_id == 0)

field_ids <- dietary_id_3$field_id

# 生成 pXXXXXX 前缀向量
p_cols <- paste0("p", field_ids)

ukb_data %>% select(starts_with(c('p26106')))

dietary_portion_id %>%
  filter(field_id %in% field_ids)

# 连续型变量数据转化，除以相应的蛋白质含量，构建统一的片为单位的数据
summary(ukb_data$p26106_i1)
summary(ukb_data$p26138_i1)

# 遍历 p_cols 中的每个前缀
for (p in p_cols) {
  # 从 p 去掉前缀 p，得到 field_id
  field_id <- as.integer(str_remove(p, "^p"))
  
  # 找到对应的 portion_grams
  portion <- dietary_portion_id %>%
    filter(field_id == !!field_id) %>%
    pull(portion_grams)
  
  # 如果找到了 portion_grams
  if (length(portion) == 1) {
    # 找出所有以该 p 开头的列
    cols_to_update <- grep(paste0("^", p), names(ukb_data), value = TRUE)
    
    # 对这些列做除法
    ukb_data <- ukb_data %>%
      mutate(across(all_of(cols_to_update), ~ .x / portion))
  }
}

# 结果验证
summary(ukb_data$p26106_i1)
summary(ukb_data$p26138_i1)

# 第四种情况white食物转化
# field_id n
# 1   100950 2
# 2   101020 2
# 3   101090 2
# 4   101160 2

# 20091	Type of sliced bread eaten

# Coding	Meaning
# 1	white
# 2	mixed
# 3	wholemeal
# 4	seeded
# 5	other

# 20092	Type of baguette eaten
# 20093	Type of large bap eaten
# 20094	Type of bread roll eaten

ukb_data %>% select(starts_with(c('p20092', 'p100950')))

table(ukb_data$p20092_i0)

ukb_data <- ukb_data %>%
  mutate(
    pw100950_i0 = if_else(p20091_i0 == "[1]", p100950_i0, NA_real_),
    pw100950_i1 = if_else(p20091_i1 == "[1]", p100950_i1, NA_real_),
    pw100950_i2 = if_else(p20091_i2 == "[1]", p100950_i2, NA_real_),
    pw100950_i3 = if_else(p20091_i3 == "[1]", p100950_i3, NA_real_),
    pw100950_i4 = if_else(p20091_i4 == "[1]", p100950_i4, NA_real_)
  ) %>%
  mutate(
    pw100950_i0 = if_else(p20092_i0 == "[1]", p101020_i0, NA_real_),
    pw100950_i1 = if_else(p20092_i1 == "[1]", p101020_i1, NA_real_),
    pw100950_i2 = if_else(p20092_i2 == "[1]", p101020_i2, NA_real_),
    pw100950_i3 = if_else(p20092_i3 == "[1]", p101020_i3, NA_real_),
    pw100950_i4 = if_else(p20092_i4 == "[1]", p101020_i4, NA_real_)
  ) %>%
  mutate(
    pw100950_i0 = if_else(p20093_i0 == "[1]", p101090_i0, NA_real_),
    pw100950_i1 = if_else(p20093_i1 == "[1]", p101090_i1, NA_real_),
    pw100950_i2 = if_else(p20093_i2 == "[1]", p101090_i2, NA_real_),
    pw100950_i3 = if_else(p20093_i3 == "[1]", p101090_i3, NA_real_),
    pw100950_i4 = if_else(p20093_i4 == "[1]", p101090_i4, NA_real_)
  ) %>%
  mutate(
    pw100950_i0 = if_else(p20094_i0 == "[1]", p101160_i0, NA_real_),
    pw100950_i1 = if_else(p20094_i1 == "[1]", p101160_i1, NA_real_),
    pw100950_i2 = if_else(p20094_i2 == "[1]", p101160_i2, NA_real_),
    pw100950_i3 = if_else(p20094_i3 == "[1]", p101160_i3, NA_real_),
    pw100950_i4 = if_else(p20094_i4 == "[1]", p101160_i4, NA_real_)
  )


### -------------------- Step 5.2 定义饮食数据的处理函数 --------------------
# 定义饮食数据的处理函数
# 通用函数：处理任意食物组
process_food_group <- function(data, prefix, field_ids) {
  library(dplyr)
  library(tidyr)
  
  # Step 1: 生成 pXXXXXX 前缀向量
  p_cols <- paste0("p", field_ids)
  
  # Step 2: 生成均值列
  for (fid in field_ids) {
    colname <- paste0(prefix, "_p", fid)
    data[[colname]] <- rowMeans(select(data, starts_with(paste0("p", fid))), na.rm = TRUE)
    data[[colname]] <- replace_na(data[[colname]], 0)
  }
  
  # Step 3: 生成总和列
  sum_col <- paste0(prefix)
  data[[sum_col]] <- rowSums(select(data, starts_with(paste0(prefix, "_p"))), na.rm = TRUE)
  
  # Step 5: 五分位分组
  quintile_col <- paste0(prefix, "_quintile")
  data[[quintile_col]] <- paste0("Q", ntile(data[[sum_col]], 5))
  
  # Step 6: 中位数分组（按中位数升序编号）
  quintile_medians <- data %>%
    group_by(.data[[quintile_col]]) %>%
    summarise(median_val = median(.data[[sum_col]], na.rm = TRUE), .groups = "drop") %>%
    arrange(median_val) %>%
    mutate(group = paste0("G", row_number()))
  
  group_col <- paste0(prefix, "_group")
  data <- data %>%
    left_join(quintile_medians %>% select(!!quintile_col, group),
              by = setNames(quintile_col, quintile_col)) %>%
    rename(!!group_col := group)
  
  # Step 6: 生成范围列（按中位数排序）
  group_stats <- data %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      median_val = median(.data[[sum_col]], na.rm = TRUE),
      min_val = min(.data[[sum_col]], na.rm = TRUE),
      max_val = max(.data[[sum_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_val) %>%
    mutate(range = paste0(round(median_val, 2), " (", round(min_val, 2), "-", round(max_val, 2), ")"))
  
  range_col <- paste0(prefix, "_group_range")
  data <- data %>%
    left_join(group_stats %>% select(!!group_col, range),
              by = setNames(group_col, group_col)) %>%
    rename(!!range_col := range)
  
  # Step 7: 生成数值型范围分组列（去重后按中位数顺序）
  range_num_col <- paste0(prefix, "_group_range_num")
  unique_levels <- unique(group_stats$range)  # 按中位数顺序且去重
  data[[range_num_col]] <- match(data[[range_col]], unique_levels)
  
  return(data)
}


## -------------------- Step 6 Food Groups 衍生 -------------------
### -------------------- Step 6.1 Berries -------------------
data <- dietary_id %>% filter(food_group == "Berries")
head(ukb_data %>% select(starts_with(c('p104470', 'p104480'))))

ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "berries",
  field_ids = data$field_id
)

ukb_data %>% count(berries_group_range)
ukb_data %>% count(berries_group_range, berries_group_range_num)
ukb_data %>% count(berries_quintile)
ukb_data %>% count(berries_group)



### -------------------- Step 6.2 Citrus fruits -------------------
data <- dietary_id %>% filter(food_group == "Citrus fruits")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "citrus_fruits",
  field_ids = data$field_id
)

### -------------------- Step 6.3 Other fruits -------------------
data <- dietary_id %>% filter(food_group == "Other fruits")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "other_fruits",
  field_ids = data$field_id
)

### -------------------- Step 6.4 Green leafy vegetables -------------------
data <- dietary_id %>% filter(food_group == "Green leafy vegetables")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "green_leafy_veg",
  field_ids = data$field_id
)

### -------------------- Step 6.5 Red/Orange vegetables -------------------
data <- dietary_id %>% filter(food_group == "Red/Orange vegetables")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "red_orange_veg",
  field_ids = data$field_id
)

### -------------------- Step 6.6 Other vegetables -------------------
data <- dietary_id %>% filter(food_group == "Other vegetables")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "other_veg",
  field_ids = data$field_id
)

### -------------------- Step 6.7 Potatoes -------------------
data <- dietary_id %>% filter(food_group == "Potatoes")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "potatoes",
  field_ids = data$field_id
)

### -------------------- Step 6.8 Refined grains -------------------
data <- dietary_id %>% filter(food_group == "Refined grains")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "refined_grains",
  field_ids = c(data$field_id, 'w100950', 'w101020', 'w101090', 'w101160') # 包含白面包的变量
)

### -------------------- Step 6.9 Whole grains -------------------
data <- dietary_id %>% filter(food_group == "Whole grains")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "whole_grains",
  field_ids = data$field_id
)

### -------------------- Step 6.10 Legumes -------------------
data <- dietary_id %>% filter(food_group == "Legumes")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "legumes",
  field_ids = data$field_id
)

### -------------------- Step 6.11 Nuts -------------------
data <- dietary_id %>% filter(food_group == "Nuts")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "nuts",
  field_ids = data$field_id
)

### -------------------- Step 6.12 Eggs -------------------
data <- dietary_id %>% filter(food_group == "Eggs")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "eggs",
  field_ids = data$field_id
)

### -------------------- Step 6.13 High-fat dairy products -------------------
data <- dietary_id %>% filter(food_group == "High-fat dairy products")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "high_fat_dairy",
  field_ids = data$field_id
)

### -------------------- Step 6.14 Low-fat dairy products -------------------
data <- dietary_id %>% filter(food_group == "Low-fat dairy products")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "low_fat_dairy",
  field_ids = data$field_id
)

### -------------------- Step 6.15 Red meats -------------------
data <- dietary_id %>% filter(food_group == "Red meats")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "red_meats",
  field_ids = data$field_id
)

### -------------------- Step 6.16 Processed meats -------------------
data <- dietary_id %>% filter(food_group == "Processed meats")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "processed_meats",
  field_ids = data$field_id
)

### -------------------- Step 6.17 Poultry -------------------
data <- dietary_id %>% filter(food_group == "Poultry")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "poultry",
  field_ids = data$field_id
)

### -------------------- Step 6.18 Fish and other seafood -------------------
data <- dietary_id %>% filter(food_group == "Fish and other seafood")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "fish_seafood",
  field_ids = data$field_id
)

### -------------------- Step 6.19 Organ meats -------------------
data <- dietary_id %>% filter(food_group == "Organ meats")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "organ_meats",
  field_ids = data$field_id
)

ukb_data %>% select(starts_with(c("organ_meats"))) %>% head(20)

max(ukb_data$organ_meats)
test_data <- ukb_data %>% select(eid, starts_with(c("p103090", "p103100", "organ_meats"))) %>% filter(organ_meats > 10)


### -------------------- Step 6.20 Butter -------------------
data <- dietary_id %>% filter(food_group == "Butter")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "butter",
  field_ids = data$field_id
)

### -------------------- Step 6.21 Margarine -------------------
data <- dietary_id %>% filter(food_group == "Margarine")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "margarine",
  field_ids = data$field_id
)

### -------------------- Step 6.22 Olive oil -------------------
data <- dietary_id %>% filter(food_group == "Olive oil")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "olive_oil",
  field_ids = data$field_id
)

### -------------------- Step 6.23 Other oil -------------------
data <- dietary_id %>% filter(food_group == "Other oil")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "other_oil",
  field_ids = data$field_id
)

ukb_data %>% count(other_oil_group_range, other_oil_group_range_num)

### -------------------- Step 6.24 Beer -------------------
data <- dietary_id %>% filter(food_group == "Beer")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "beer",
  field_ids = data$field_id
)

### -------------------- Step 6.25 Coffee -------------------
data <- dietary_id %>% filter(food_group == "Coffee")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "coffee",
  field_ids = data$field_id
)

### -------------------- Step 6.26 Fruit juices -------------------
data <- dietary_id %>% filter(food_group == "Fruit juices")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "fruit_juices",
  field_ids = data$field_id
)

### -------------------- Step 6.27 Liquor -------------------
data <- dietary_id %>% filter(food_group == "Liquor")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "liquor",
  field_ids = data$field_id
)

### -------------------- Step 6.28 Sweetened beverages -------------------
data <- dietary_id %>% filter(food_group == "Sweetened beverages")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "sweetened_bev",
  field_ids = data$field_id
)

### -------------------- Step 6.29 Tea -------------------
data <- dietary_id %>% filter(food_group == "Tea")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "tea",
  field_ids = data$field_id
)

### -------------------- Step 6.30 Water -------------------
data <- dietary_id %>% filter(food_group == "Water")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "water",
  field_ids = data$field_id
)

### -------------------- Step 6.31 Wine -------------------
data <- dietary_id %>% filter(food_group == "Wine")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "wine",
  field_ids = data$field_id
)

### -------------------- Step 6.32 Pizza -------------------
data <- dietary_id %>% filter(food_group == "Pizza")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "pizza",
  field_ids = data$field_id
)

### -------------------- Step 6.33 Snacks and pastries -------------------
data <- dietary_id %>% filter(food_group == "Snacks and pastries")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "snacks_pastries",
  field_ids = data$field_id
)

### -------------------- Step 6.34 Fried and fast foods -------------------
data <- dietary_id %>% filter(food_group == "Fried and fast foods")
ukb_data <- process_food_group(
  data = ukb_data,
  prefix = "fried_fast_foods",
  field_ids = data$field_id
)


## -------------------- Step 7 结果保存 -------------------
colnames(ukb_data)

### -------------------- Step 7.1 饮食结果保存 -------------------
col_str <- paste(names(ukb_data), collapse = ", ")
col_str

unique(dietary_id$food_group)

ukb_data_food_group <- ukb_data %>%
  select(eid, date_food_baseline,
         berries, berries_quintile, berries_group, berries_group_range, berries_group_range_num,
         citrus_fruits, citrus_fruits_quintile, citrus_fruits_group, citrus_fruits_group_range, citrus_fruits_group_range_num,
         other_fruits, other_fruits_quintile, other_fruits_group, other_fruits_group_range, other_fruits_group_range_num,
         green_leafy_veg, green_leafy_veg_quintile, green_leafy_veg_group, green_leafy_veg_group_range, green_leafy_veg_group_range_num,
         red_orange_veg, red_orange_veg_quintile, red_orange_veg_group, red_orange_veg_group_range, red_orange_veg_group_range_num,
         other_veg, other_veg_quintile, other_veg_group, other_veg_group_range, other_veg_group_range_num,
         potatoes, potatoes_quintile, potatoes_group, potatoes_group_range, potatoes_group_range_num,
         refined_grains, refined_grains_quintile, refined_grains_group, refined_grains_group_range, refined_grains_group_range_num,
         whole_grains, whole_grains_quintile, whole_grains_group, whole_grains_group_range, whole_grains_group_range_num,
         legumes, legumes_quintile, legumes_group, legumes_group_range, legumes_group_range_num,
         nuts, nuts_quintile, nuts_group, nuts_group_range, nuts_group_range_num,
         eggs, eggs_quintile, eggs_group, eggs_group_range, eggs_group_range_num,
         high_fat_dairy, high_fat_dairy_quintile, high_fat_dairy_group, high_fat_dairy_group_range, high_fat_dairy_group_range_num,
         low_fat_dairy, low_fat_dairy_quintile, low_fat_dairy_group, low_fat_dairy_group_range, low_fat_dairy_group_range_num,
         red_meats, red_meats_quintile, red_meats_group, red_meats_group_range, red_meats_group_range_num,
         processed_meats, processed_meats_quintile, processed_meats_group, processed_meats_group_range, processed_meats_group_range_num,
         poultry, poultry_quintile, poultry_group, poultry_group_range, poultry_group_range_num,
         fish_seafood, fish_seafood_quintile, fish_seafood_group, fish_seafood_group_range, fish_seafood_group_range_num,
         organ_meats, organ_meats_quintile, organ_meats_group, organ_meats_group_range, organ_meats_group_range_num,
         butter, butter_quintile, butter_group, butter_group_range, butter_group_range_num,
         margarine, margarine_quintile, margarine_group, margarine_group_range, margarine_group_range_num,
         olive_oil, olive_oil_quintile, olive_oil_group, olive_oil_group_range, olive_oil_group_range_num,
         other_oil, other_oil_quintile, other_oil_group, other_oil_group_range, other_oil_group_range_num,
         beer, beer_quintile, beer_group, beer_group_range, beer_group_range_num,
         coffee, coffee_quintile, coffee_group, coffee_group_range, coffee_group_range_num,
         fruit_juices, fruit_juices_quintile, fruit_juices_group, fruit_juices_group_range, fruit_juices_group_range_num,
         liquor, liquor_quintile, liquor_group, liquor_group_range, liquor_group_range_num,
         sweetened_bev, sweetened_bev_quintile, sweetened_bev_group, sweetened_bev_group_range, sweetened_bev_group_range_num,
         tea, tea_quintile, tea_group, tea_group_range, tea_group_range_num,
         water, water_quintile, water_group, water_group_range, water_group_range_num,
         wine, wine_quintile, wine_group, wine_group_range, wine_group_range_num,
         pizza, pizza_quintile, pizza_group, pizza_group_range, pizza_group_range_num,
         snacks_pastries, snacks_pastries_quintile, snacks_pastries_group, snacks_pastries_group_range, snacks_pastries_group_range_num,
         fried_fast_foods, fried_fast_foods_quintile, fried_fast_foods_group, fried_fast_foods_group_range, fried_fast_foods_group_range_num)


# 结果保存
save(ukb_data_food_group, file = './data/ukb_data_food_group.RData')

### -------------------- Step 7.2 随访时间结果保存 -------------------

ukb_data_food_date_group <- ukb_data %>%
  select(eid, date_food_baseline)

# 结果保存
save(ukb_data_food_date_group, file = './data/ukb_data_food_date_group.RData')

########################################################################################
# END ANALYSIS BLOCK 01/11: Step_02_1_1
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 02/11: Step_03_1_0_SPLIT
# Source: Step_03_1_0_ukb_70_30_split_manifest_and_cohort_generation(1).R
# Provenance: Library source dated 2026-11-23; SPLIT_DIR normalized to the shared path used by all downstream Step_03/Step_04 scripts
# SHA256 of included block: 5a965f9e598a80396a5233831814a5606262adc9a745bca3c830eea4d9ece436
########################################################################################

# =============================================================================
# STEP_03_1_0
# UK Biobank 70/30 split manifest and cohort generation
#
# Purpose:
#   1. Reconstruct the eligible depression cohort using the original
#      Step_03_1_2 REF2 cohort-definition rules.
#   2. Create the 70/30 split manifest once, or reuse the existing manifest
#      without replacing it.
#   3. Save all split-related outputs directly under:
#        ./results/Step_03/ukb_70_30_split/
#
# Important:
#   - This script does NOT run food-group association models or generate plots.
#   - An existing ukb_70_30_split_manifest.csv is NEVER overwritten.
#   - If an existing manifest does not match the current eligible cohort or
#     split settings, the script stops; it does not silently create a new split.
#   - No files or folders outside SPLIT_DIR are created, modified, or deleted.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# 0. Global settings
# -----------------------------------------------------------------------------

SPLIT_SEED <- 20261103L
TRAIN_PROP <- 0.70

SPLIT_DIR <- "./results/Step_03/ukb_70_30_split"
SPLIT_MANIFEST_FILE <- file.path(
  SPLIT_DIR,
  "ukb_70_30_split_manifest.csv"
)

dir.create(
  SPLIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(SPLIT_DIR)) {
  stop("STOP: Could not create output directory: ", SPLIT_DIR)
}

# -----------------------------------------------------------------------------
# 1. Load the four original UK Biobank source datasets
# -----------------------------------------------------------------------------

load_required_object <- function(path, object_name) {

  if (!file.exists(path)) {
    stop("STOP: Required input file was not found: ", path)
  }

  loaded_objects <- load(path, envir = .GlobalEnv)

  if (!object_name %in% loaded_objects ||
      !exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
    stop(
      "STOP: Input file does not contain the expected object '",
      object_name, "': ", path
    )
  }

  message("Loaded: ", path, " -> ", object_name)
  invisible(TRUE)
}

load_required_object(
  "./data/ukb_data_food_group.RData",
  "ukb_data_food_group"
)

load_required_object(
  "./data/ukb_data_outcomes.RData",
  "ukb_data_outcomes"
)

load_required_object(
  "./data/ukb_data_cancer_outcomes.RData",
  "ukb_data_cancer_outcomes"
)

load_required_object(
  "./data/ukb_data_covariates.RData",
  "ukb_data_covariates"
)

# -----------------------------------------------------------------------------
# 2. Reconstruct the original analysis dataset
# -----------------------------------------------------------------------------
# Keep the original left_join() behavior and source order.
# No alternative data source or replacement mapping is introduced.

ukb_data <- ukb_data_food_group %>%
  left_join(ukb_data_outcomes) %>%
  left_join(ukb_data_cancer_outcomes) %>%
  left_join(ukb_data_covariates)

# -----------------------------------------------------------------------------
# 3. Validate columns required by the original cohort construction
# -----------------------------------------------------------------------------

required_columns <- c(
  "depression_status",
  "twnsend",
  "career",
  "smoking_status",
  "drink_status",
  "bmi",
  "education_attainment",
  "physical_activity_proup",
  "sedentary_group",
  "central_obesity",
  "triglycerides",
  "hdl_cholesterol",
  "household_income",
  "energy_intake",
  "t2diabetes_baseline_status",
  "t2diabetes_baseline_dx",
  "all_stroke_baseline_status",
  "all_stroke_baseline_dx",
  "hypertension_baseline_status",
  "hypertension_baseline_dx",
  "dyslipidemia_baseline_status",
  "dyslipidemia_baseline_dx",
  "AF_baseline_status",
  "AF_baseline_dx",
  "HF_baseline_status",
  "HF_baseline_dx",
  "MI_baseline_status",
  "MI_baseline_dx",
  "all_cancer_baseline_status",
  "all_cancer_baseline_dx",
  "lost_follow_up_status",
  "depression_time",
  "depression_baseline_dx"
)

missing_required_columns <- setdiff(
  required_columns,
  names(ukb_data)
)

if (length(missing_required_columns) > 0L) {
  stop(
    "STOP: Required columns are missing after data combination: ",
    paste(missing_required_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# 4. Original covariate complete-case exclusions
# -----------------------------------------------------------------------------
# These exclusions are retained exactly from the source script.

n_before_covariate_filter <- nrow(ukb_data)

ukb_data <- ukb_data %>%
  filter(!is.na(twnsend)) %>%
  filter(!is.na(career)) %>%
  filter(!is.na(smoking_status)) %>%
  filter(!is.na(drink_status)) %>%
  filter(!is.na(bmi)) %>%
  filter(!is.na(education_attainment)) %>%
  filter(!is.na(physical_activity_proup)) %>%
  filter(!is.na(sedentary_group)) %>%
  filter(!is.na(central_obesity)) %>%
  filter(!is.na(triglycerides)) %>%
  filter(!is.na(hdl_cholesterol)) %>%
  filter(!is.na(household_income)) %>%
  filter(!is.na(energy_intake))

n_after_covariate_filter <- nrow(ukb_data)

# -----------------------------------------------------------------------------
# 5. Original derived baseline disease-status variables
# -----------------------------------------------------------------------------
# Retained so the saved cohort RDS preserves the original data structure.

ukb_data <- ukb_data %>%
  mutate(
    t2diabetes_status_final = pmax(
      t2diabetes_baseline_status,
      t2diabetes_baseline_dx,
      na.rm = TRUE
    ),
    all_stroke_status_final = pmax(
      all_stroke_baseline_status,
      all_stroke_baseline_dx,
      na.rm = TRUE
    ),
    hypertension_status_final = pmax(
      hypertension_baseline_status,
      hypertension_baseline_dx,
      na.rm = TRUE
    ),
    dyslipidemia_status_final = pmax(
      dyslipidemia_baseline_status,
      dyslipidemia_baseline_dx,
      na.rm = TRUE
    ),
    AF_status_final = pmax(
      AF_baseline_status,
      AF_baseline_dx,
      na.rm = TRUE
    ),
    HF_status_final = pmax(
      HF_baseline_status,
      HF_baseline_dx,
      na.rm = TRUE
    ),
    MI_status_final = pmax(
      MI_baseline_status,
      MI_baseline_dx,
      na.rm = TRUE
    ),
    all_cancer_status_final = pmax(
      all_cancer_baseline_status,
      all_cancer_baseline_dx,
      na.rm = TRUE
    ),
    cvd_status_final = pmax(
      all_stroke_status_final,
      AF_status_final,
      HF_status_final,
      MI_status_final,
      na.rm = TRUE
    )
  )

# -----------------------------------------------------------------------------
# 6. Original follow-up and depression baseline eligibility filters
# -----------------------------------------------------------------------------

n_before_followup_filter <- nrow(ukb_data)

ukb_data <- ukb_data %>%
  filter(lost_follow_up_status == 0)

n_after_followup_filter <- nrow(ukb_data)

# The original source script defines the eligible depression cohort using:
#   depression_time > 0
#   depression_baseline_dx != 1
#
# filter() also removes rows for which either condition evaluates to NA,
# matching the original dplyr filtering behavior.

ukb_data_depression <- ukb_data %>%
  filter(depression_time > 0) %>%
  filter(depression_baseline_dx != 1)

if (nrow(ukb_data_depression) == 0L) {
  stop("STOP: The eligible depression cohort contains zero participants.")
}

message(
  "Cohort construction completed. ",
  "N before covariate exclusions = ", n_before_covariate_filter,
  "; after covariate exclusions = ", n_after_covariate_filter,
  "; before follow-up exclusion = ", n_before_followup_filter,
  "; after follow-up exclusion = ", n_after_followup_filter,
  "; eligible depression cohort N = ", nrow(ukb_data_depression)
)

# -----------------------------------------------------------------------------
# 7. Original food-group factor coding
# -----------------------------------------------------------------------------
# Preserve the 34 food-group factor conversions used before the original
# training/test RDS files were saved. For each variable, the reference level
# remains the second sorted observed value, followed by the other sorted values.
#
# This step changes factor coding only; it does not change cohort membership.

food_group_range_variables <- c(
  "berries_group_range",
  "citrus_fruits_group_range",
  "other_fruits_group_range",
  "green_leafy_veg_group_range",
  "red_orange_veg_group_range",
  "other_veg_group_range",
  "potatoes_group_range",
  "refined_grains_group_range",
  "whole_grains_group_range",
  "legumes_group_range",
  "nuts_group_range",
  "eggs_group_range",
  "high_fat_dairy_group_range",
  "low_fat_dairy_group_range",
  "red_meats_group_range",
  "processed_meats_group_range",
  "poultry_group_range",
  "fish_seafood_group_range",
  "organ_meats_group_range",
  "butter_group_range",
  "margarine_group_range",
  "olive_oil_group_range",
  "other_oil_group_range",
  "beer_group_range",
  "coffee_group_range",
  "fruit_juices_group_range",
  "liquor_group_range",
  "sweetened_bev_group_range",
  "tea_group_range",
  "water_group_range",
  "wine_group_range",
  "pizza_group_range",
  "snacks_pastries_group_range",
  "fried_fast_foods_group_range"
)

missing_food_group_columns <- setdiff(
  food_group_range_variables,
  names(ukb_data_depression)
)

if (length(missing_food_group_columns) > 0L) {
  stop(
    "STOP: Food-group range columns required for original factor coding ",
    "are missing: ",
    paste(missing_food_group_columns, collapse = ", ")
  )
}

for (variable_name in food_group_range_variables) {

  observed_values <- sort(
    unique(ukb_data_depression[[variable_name]])
  )

  if (length(observed_values) < 2L) {
    stop(
      "STOP: Fewer than two observed categories for ",
      variable_name,
      ". Original reference-level coding cannot be reproduced."
    )
  }

  factor_levels <- c(
    observed_values[2],
    observed_values[-2]
  )

  ukb_data_depression[[variable_name]] <- factor(
    ukb_data_depression[[variable_name]],
    levels = factor_levels,
    ordered = FALSE
  )
}

# -----------------------------------------------------------------------------
# 8. Identify and validate the participant ID column
# -----------------------------------------------------------------------------

ID_CANDIDATES <- c(
  "eid",
  "f.eid",
  "participant_id",
  "participant.id",
  "participantid",
  "ukb_id",
  "ukb_eid",
  "id"
)

id_match <- match(
  ID_CANDIDATES,
  tolower(names(ukb_data_depression))
)

id_match <- id_match[!is.na(id_match)]

if (length(id_match) > 0L) {

  ID_COL <- names(ukb_data_depression)[id_match[1]]

} else {

  ID_COL <- names(ukb_data_depression)[1]

  warning(
    "No standard participant ID column name was found. ",
    "The first column will be used as ID: ", ID_COL,
    ". Please verify that this is the unique UK Biobank participant ID."
  )
}

id_values <- trimws(
  as.character(ukb_data_depression[[ID_COL]])
)

if (anyNA(id_values) || any(id_values == "")) {
  stop("STOP: Missing or blank participant IDs detected.")
}

if (anyDuplicated(id_values) > 0L) {
  stop(
    "STOP: Participant IDs are not unique in the eligible depression cohort. ",
    "Resolve duplicated participant rows before splitting."
  )
}

eligible_ids <- sort(
  id_values,
  method = "radix"
)

n_eligible <- length(eligible_ids)
n_train_expected <- as.integer(round(TRAIN_PROP * n_eligible))
n_test_expected <- n_eligible - n_train_expected

if (n_train_expected < 1L || n_test_expected < 1L) {
  stop("STOP: The eligible cohort is too small for a 70/30 split.")
}

# -----------------------------------------------------------------------------
# 9. Create or reuse the frozen shared manifest
# -----------------------------------------------------------------------------
# Critical rule:
#   If the manifest already exists, read it and NEVER overwrite it.
#   If it does not exist, create it once using the frozen seed.

manifest_preexisting <- file.exists(SPLIT_MANIFEST_FILE)

if (manifest_preexisting) {

  message(
    "Existing manifest detected. It will be reused and will NOT be replaced: ",
    normalizePath(SPLIT_MANIFEST_FILE, mustWork = TRUE)
  )

  split_manifest <- read.csv(
    SPLIT_MANIFEST_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(
      participant_id = "character",
      split = "character",
      id_column = "character",
      split_seed = "integer",
      training_fraction = "numeric"
    )
  )

} else {

  message(
    "No existing manifest found. Creating the frozen 70/30 manifest."
  )

  set.seed(SPLIT_SEED)

  training_ids_created <- sample(
    eligible_ids,
    size = n_train_expected,
    replace = FALSE
  )

  split_manifest <- data.frame(
    participant_id = eligible_ids,
    split = ifelse(
      eligible_ids %in% training_ids_created,
      "training",
      "test"
    ),
    id_column = ID_COL,
    split_seed = SPLIT_SEED,
    training_fraction = TRAIN_PROP,
    stringsAsFactors = FALSE
  )

  # This is the only code path that writes the manifest.
  # The existence check above prevents replacement of an existing manifest.
  write.csv(
    split_manifest,
    file = SPLIT_MANIFEST_FILE,
    row.names = FALSE,
    na = ""
  )
}

# -----------------------------------------------------------------------------
# 10. Strict manifest validation
# -----------------------------------------------------------------------------

required_manifest_columns <- c(
  "participant_id",
  "split",
  "id_column",
  "split_seed",
  "training_fraction"
)

if (!all(required_manifest_columns %in% names(split_manifest))) {
  stop(
    "STOP: The shared manifest has an unexpected structure. ",
    "The existing manifest has not been modified."
  )
}

if (anyNA(split_manifest$participant_id) ||
    any(trimws(split_manifest$participant_id) == "")) {
  stop("STOP: Missing or blank participant IDs in the shared manifest.")
}

if (anyDuplicated(split_manifest$participant_id) > 0L) {
  stop("STOP: Duplicate participant IDs in the shared manifest.")
}

if (!setequal(split_manifest$participant_id, eligible_ids)) {

  missing_from_manifest <- length(
    setdiff(eligible_ids, split_manifest$participant_id)
  )

  extra_in_manifest <- length(
    setdiff(split_manifest$participant_id, eligible_ids)
  )

  stop(
    "STOP: Current eligible cohort does not exactly match the shared manifest. ",
    "IDs missing from manifest = ", missing_from_manifest,
    "; extra IDs in manifest = ", extra_in_manifest,
    ". The existing manifest was NOT replaced. ",
    "Check the source data and cohort inclusion rules."
  )
}

if (anyNA(split_manifest$id_column) ||
    any(split_manifest$id_column != ID_COL)) {
  stop(
    "STOP: Participant ID column differs from the frozen manifest. ",
    "The existing manifest was NOT replaced."
  )
}

if (anyNA(split_manifest$split_seed) ||
    any(split_manifest$split_seed != SPLIT_SEED)) {
  stop(
    "STOP: Split seed differs from the frozen manifest. ",
    "The existing manifest was NOT replaced."
  )
}

if (anyNA(split_manifest$training_fraction) ||
    any(abs(split_manifest$training_fraction - TRAIN_PROP) > 1e-12)) {
  stop(
    "STOP: Training fraction differs from the frozen manifest. ",
    "The existing manifest was NOT replaced."
  )
}

if (!all(split_manifest$split %in% c("training", "test"))) {
  stop("STOP: The manifest contains invalid split labels.")
}

if (sum(split_manifest$split == "training") != n_train_expected ||
    sum(split_manifest$split == "test") != n_test_expected) {
  stop(
    "STOP: Training/test counts in the manifest do not match the 70/30 rule."
  )
}

# Align the manifest to the sorted eligible ID order.
split_manifest <- split_manifest[
  match(eligible_ids, split_manifest$participant_id),
  ,
  drop = FALSE
]

if (anyNA(split_manifest$split)) {
  stop("STOP: Some eligible participants have no split assignment.")
}

training_ids <- split_manifest$participant_id[
  split_manifest$split == "training"
]

test_ids <- split_manifest$participant_id[
  split_manifest$split == "test"
]

# -----------------------------------------------------------------------------
# 11. Save the shared participant ID lists
# -----------------------------------------------------------------------------

write.csv(
  data.frame(participant_id = training_ids),
  file = file.path(
    SPLIT_DIR,
    "ukb_70_30_training_ids.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(participant_id = test_ids),
  file = file.path(
    SPLIT_DIR,
    "ukb_70_30_test_ids.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 12. Apply the frozen split to the factorized eligible cohort
# -----------------------------------------------------------------------------

ukb_data_depression_full <- ukb_data_depression

id_values_factorized <- trimws(
  as.character(ukb_data_depression_full[[ID_COL]])
)

if (anyNA(id_values_factorized) ||
    any(id_values_factorized == "") ||
    anyDuplicated(id_values_factorized) > 0L) {
  stop("STOP: Participant ID QC failed after factor conversion.")
}

split_lookup <- setNames(
  split_manifest$split,
  split_manifest$participant_id
)

ukb_data_depression_full$analysis_split <- unname(
  split_lookup[id_values_factorized]
)

if (anyNA(ukb_data_depression_full$analysis_split)) {
  stop(
    "STOP: Some eligible participants could not be assigned ",
    "to training or test."
  )
}

ukb_data_depression_train <- ukb_data_depression_full %>%
  filter(analysis_split == "training")

ukb_data_depression_test <- ukb_data_depression_full %>%
  filter(analysis_split == "test")

# -----------------------------------------------------------------------------
# 13. Exact ID identity, overlap, and reconstruction checks
# -----------------------------------------------------------------------------

train_ids_observed <- sort(
  as.character(ukb_data_depression_train[[ID_COL]]),
  method = "radix"
)

test_ids_observed <- sort(
  as.character(ukb_data_depression_test[[ID_COL]]),
  method = "radix"
)

if (!identical(
  train_ids_observed,
  sort(training_ids, method = "radix")
)) {
  stop("STOP: Training participant IDs do not match the shared manifest.")
}

if (!identical(
  test_ids_observed,
  sort(test_ids, method = "radix")
)) {
  stop("STOP: Test participant IDs do not match the shared manifest.")
}

train_test_overlap_n <- length(
  intersect(train_ids_observed, test_ids_observed)
)

if (train_test_overlap_n != 0L) {
  stop("STOP: Participant overlap detected between training and test sets.")
}

if (!identical(
  sort(c(train_ids_observed, test_ids_observed), method = "radix"),
  eligible_ids
)) {
  stop(
    "STOP: Training/test union does not reconstruct the eligible cohort."
  )
}

if (nrow(ukb_data_depression_train) != n_train_expected ||
    nrow(ukb_data_depression_test) != n_test_expected) {
  stop("STOP: Saved cohort row counts do not match the expected 70/30 split.")
}

# -----------------------------------------------------------------------------
# 14. Cohort summary
# -----------------------------------------------------------------------------

summarise_depression_split <- function(dat, split_label) {

  tibble(
    split = split_label,
    participants = nrow(dat),
    depression_events = sum(
      dat$depression_status == 1,
      na.rm = TRUE
    ),
    non_events = sum(
      dat$depression_status == 0,
      na.rm = TRUE
    ),
    missing_depression_status = sum(
      is.na(dat$depression_status)
    )
  )
}

split_cohort_summary <- bind_rows(
  summarise_depression_split(
    ukb_data_depression_full,
    "full_cohort"
  ),
  summarise_depression_split(
    ukb_data_depression_train,
    "training"
  ),
  summarise_depression_split(
    ukb_data_depression_test,
    "test"
  )
) %>%
  mutate(
    event_percent = ifelse(
      participants > 0,
      100 * depression_events / participants,
      NA_real_
    )
  )

write.csv(
  split_cohort_summary,
  file = file.path(
    SPLIT_DIR,
    "ukb_70_30_split_cohort_summary.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 15. Save training and test cohort RDS files
# -----------------------------------------------------------------------------
# Both files are saved directly in SPLIT_DIR; no main/ref0 subfolder is used.

saveRDS(
  ukb_data_depression_train,
  file = file.path(
    SPLIT_DIR,
    "ukb_depression_training_cohort.rds"
  )
)

saveRDS(
  ukb_data_depression_test,
  file = file.path(
    SPLIT_DIR,
    "ukb_depression_test_cohort.rds"
  )
)

# -----------------------------------------------------------------------------
# 16. Split run QC
# -----------------------------------------------------------------------------

manifest_md5 <- unname(
  tools::md5sum(SPLIT_MANIFEST_FILE)
)

split_run_qc <- data.frame(
  script_tag = "Step_03_1_0",
  id_column = ID_COL,
  split_seed = SPLIT_SEED,
  training_fraction = TRAIN_PROP,
  eligible_n = nrow(ukb_data_depression_full),
  training_n = nrow(ukb_data_depression_train),
  test_n = nrow(ukb_data_depression_test),
  eligible_events = sum(
    ukb_data_depression_full$depression_status == 1,
    na.rm = TRUE
  ),
  training_events = sum(
    ukb_data_depression_train$depression_status == 1,
    na.rm = TRUE
  ),
  test_events = sum(
    ukb_data_depression_test$depression_status == 1,
    na.rm = TRUE
  ),
  train_ids_match_manifest = TRUE,
  test_ids_match_manifest = TRUE,
  train_test_overlap_n = train_test_overlap_n,
  manifest_md5 = manifest_md5,
  manifest_reused = manifest_preexisting,
  stringsAsFactors = FALSE
)

write.csv(
  split_run_qc,
  file = file.path(
    SPLIT_DIR,
    "split_run_qc.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 17. Final QC report and output inventory
# -----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP_03_1_0: UKB 70/30 SPLIT GENERATION COMPLETED\n")
cat("============================================================\n")
cat("Output directory: ", normalizePath(SPLIT_DIR), "\n", sep = "")
cat("Manifest reused (not replaced): ", manifest_preexisting, "\n", sep = "")
cat("Eligible cohort N: ", nrow(ukb_data_depression_full), "\n", sep = "")
cat("Training N: ", nrow(ukb_data_depression_train), "\n", sep = "")
cat("Test N: ", nrow(ukb_data_depression_test), "\n", sep = "")
cat("Training/test overlap: ", train_test_overlap_n, "\n", sep = "")
cat("Manifest MD5: ", manifest_md5, "\n", sep = "")
cat("------------------------------------------------------------\n")
print(split_cohort_summary)
cat("------------------------------------------------------------\n")
print(split_run_qc)
cat("------------------------------------------------------------\n")
cat("Files in the split output directory:\n")
print(list.files(SPLIT_DIR, full.names = FALSE))
cat("============================================================\n")

# Final hard check: all required outputs must exist in the one target folder.
required_output_files <- c(
  "ukb_70_30_split_manifest.csv",
  "ukb_70_30_training_ids.csv",
  "ukb_70_30_test_ids.csv",
  "ukb_70_30_split_cohort_summary.csv",
  "ukb_depression_training_cohort.rds",
  "ukb_depression_test_cohort.rds",
  "split_run_qc.csv"
)

missing_output_files <- required_output_files[
  !file.exists(file.path(SPLIT_DIR, required_output_files))
]

if (length(missing_output_files) > 0L) {
  stop(
    "STOP: One or more expected split outputs were not created: ",
    paste(missing_output_files, collapse = ", ")
  )
}

message(
  "PASS: All Step_03_1_0 outputs are saved under ",
  normalizePath(SPLIT_DIR),
  ". No analysis models or plots were generated."
)

########################################################################################
# END ANALYSIS BLOCK 02/11: Step_03_1_0_SPLIT
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 03/11: Step_03_1_0_MODEL2_REFERENCE
# Source: Step_03_1_0_association_between_food_group_and_depression_70_30_MAIN_ROTATING_LOWEST_RISK_REFERENCE_V1.R
# Provenance: Library source dated 2026-11-05
# SHA256 of included block: 2680d715a0f5bf94627c65a8c732c2d72178b968dc8e1baa3f8d7865d3a352b6
########################################################################################

# Step_03_1_1 MAIN — Lowest-risk reference rotation (full script)
# Based on the supplied Step_03_1_1_association_between_food_group_and_depression_70_30_MAIN(6).R
# The original frozen 70/30 manifest (seed 20261103) is reused without changes.
# Original food categorization, Cox adjustments, trend formula and forestplot
# display settings are preserved. Outputs go into a new directory to avoid
# overwriting the previously frozen MAIN figures and CSV files.
# Refitting all category references is a categorical Cox analysis only; it does
# not automatically resolve RCS reference-choice or optimum-threshold concerns.

rm(list = ls())
## -------------------- Step 0 Loading Packages. --------------------
library(tidyverse)
library(survival)
# library(survivalAnalysis)
library(ggsci)
library(cols4all) # 新的颜色包
library(table1)
library(survminer)
library(gtsummary)
library(jstable)
library(forestplot)
library(glue)
library(Hmisc)
# table 1
library(gtsummary)
library(forestplot)
library(EValue)

## -------------------- Step 1 Data prepare --------------------
### -------------------- Step 1.1 loading data --------------------
# food
load('./data/ukb_data_food_group.RData')
# prs
# load('./data/ukb_data_prs.RData')
# outcomes: diseases
load('./data/ukb_data_outcomes.RData')
# outcomes: cancers
load('./data/ukb_data_cancer_outcomes.RData')
# covariates
load("./data/ukb_data_covariates.RData")

### -------------------- Step 1.2 data combination --------------------
ukb_data <- ukb_data_food_group %>%
  # left_join(ukb_data_prs) %>%
  left_join(ukb_data_outcomes) %>%
  left_join(ukb_data_cancer_outcomes) %>%
  left_join(ukb_data_covariates)

colnames(ukb_data)

table(ukb_data$depression_status)

### -------------------- Step 1.3 缺失值处理和疾病状态生成 --------------------
# 查看缺失值情况
data.frame(colSums(is.na(ukb_data[, -1])))

#### -------------------- Step 1.3.1 删除协变量缺失值 --------------------
ukb_data <- ukb_data %>%
  # filter(!is.na(depression_prs)) %>%
  filter(!is.na(twnsend)) %>%
  filter(!is.na(career)) %>%
  filter(!is.na(smoking_status)) %>%
  filter(!is.na(drink_status)) %>%
  filter(!is.na(bmi)) %>%
  filter(!is.na(education_attainment)) %>%
  filter(!is.na(physical_activity_proup)) %>%
  filter(!is.na(sedentary_group)) %>%
  filter(!is.na(central_obesity)) %>%
  filter(!is.na(triglycerides)) %>%
  filter(!is.na(hdl_cholesterol)) %>%
  filter(!is.na(household_income)) %>%
  filter(!is.na(energy_intake))

#### -------------------- Step 1.3.2 协变量诊断衍生 --------------------
colnames(ukb_data)

# 整个随访时间出现的疾病状态
ukb_data %>% select(HF_status)
# accelerometry 90011	End time of wear时间前出现的疾病状态
ukb_data %>% select(HF_baseline_status)

# 所有的诊断时间都减去基线时间，即90011	End time of wear

# 协变量中的疾病状态生成
ukb_data <- ukb_data %>%
  mutate(
    t2diabetes_status_final = pmax(t2diabetes_baseline_status, t2diabetes_baseline_dx, na.rm = TRUE),
    all_stroke_status_final = pmax(all_stroke_baseline_status, all_stroke_baseline_dx, na.rm = TRUE),
    hypertension_status_final = pmax(hypertension_baseline_status, hypertension_baseline_dx, na.rm = TRUE),
    dyslipidemia_status_final = pmax(dyslipidemia_baseline_status, dyslipidemia_baseline_dx, na.rm = TRUE),
    AF_status_final = pmax(AF_baseline_status, AF_baseline_dx, na.rm = TRUE),
    HF_status_final = pmax(HF_baseline_status, HF_baseline_dx, na.rm = TRUE),
    MI_status_final = pmax(MI_baseline_status, MI_baseline_dx, na.rm = TRUE),
    all_cancer_status_final = pmax(all_cancer_baseline_status, all_cancer_baseline_dx, na.rm = TRUE),
    cvd_status_final = pmax(all_stroke_status_final, AF_status_final, HF_status_final, MI_status_final, na.rm = TRUE)
  )

# # 肿瘤状态的NA值替换为0
# ukb_data <- ukb_data %>%
#   mutate(
#     all_cancer_status_final = ifelse(is.na(all_cancer_status_final), 0, all_cancer_status_final)
#   )

ukb_data %>% count(all_cancer_status_final)
ukb_data %>% count(HF_status_final)

## -------------------- Step 2 Data filter --------------------
### -------------------- Step 2.1 lost to follow up --------------------
ukb_data %>% count(lost_follow_up_status)
ukb_data <- ukb_data %>% 
  filter(lost_follow_up_status == 0)

### -------------------- Step 2.2 基线状态过滤 --------------------
colnames(ukb_data %>% select(ends_with("_time")))
colnames(ukb_data %>% select(ends_with("_dx")))
colnames(ukb_data %>% select(ends_with("_meds")))


# list of diseases
events <- c(
  Cs(
    depression
  )
)


table(ukb_data$depression_status)

### -------------------- Step 2.3 depression基线ICD和自我报告状态过滤 --------------------

for (event in events) {
  assign(paste0("ukb_data_", event), ukb_data %>%
           filter(eval(parse(text = paste0(event, "_time > 0")))) %>%
           filter(eval(parse(text = paste0(event, "_baseline_dx != 1")))))
}



## -------------------- Step 2.4 70/30 split manifest --------------------
# This block is intentionally identical in the main and reference-0 scripts,
# except for the script tag and whether the script is allowed to create the
# shared manifest. The eligible cohort has already been defined above.
#
# The split is created once from participant IDs, without using food-outcome
# associations, P values, HRs, or model performance. The same manifest is
# then reused and strictly checked by both scripts.

SPLIT_SEED <- 20261103L
TRAIN_PROP <- 0.70
SPLIT_DIR <- "./results/Step_03/ukb_70_30_split"
SPLIT_MANIFEST_FILE <- file.path(SPLIT_DIR, "ukb_70_30_split_manifest.csv")
SPLIT_CREATE_IF_MISSING <- FALSE  # Frozen original manifest is required
SPLIT_SCRIPT_TAG <- "main"

dir.create(SPLIT_DIR, recursive = TRUE, showWarnings = FALSE)

# Identify the participant ID column. The original script excludes the first
# column when checking missingness; if no standard ID name is found, that
# first column is used as a fallback and is subjected to strict uniqueness QC.
ID_CANDIDATES <- c(
  "eid", "f.eid", "participant_id", "participant.id",
  "participantid", "ukb_id", "ukb_eid", "id"
)
id_match <- match(ID_CANDIDATES, tolower(names(ukb_data_depression)))
id_match <- id_match[!is.na(id_match)]

if (length(id_match) > 0L) {
  ID_COL <- names(ukb_data_depression)[id_match[1]]
} else {
  ID_COL <- names(ukb_data_depression)[1]
  message(
    "WARNING: No standard participant ID column name was found. ",
    "Using the first column as ID: ", ID_COL,
    ". Please verify this is the unique UK Biobank participant identifier."
  )
}

id_values <- trimws(as.character(ukb_data_depression[[ID_COL]]))

if (anyNA(id_values) || any(id_values == "")) {
  stop("STOP: Missing/blank participant IDs detected before the 70/30 split.")
}
if (anyDuplicated(id_values) > 0L) {
  stop(
    "STOP: Participant IDs are not unique in the eligible depression cohort. ",
    "Resolve duplicated participant rows before splitting."
  )
}

eligible_ids <- sort(id_values, method = "radix")
n_eligible <- length(eligible_ids)
n_train_expected <- as.integer(round(TRAIN_PROP * n_eligible))
n_test_expected <- n_eligible - n_train_expected

if (n_train_expected < 1L || n_test_expected < 1L) {
  stop("STOP: The eligible cohort is too small for a 70/30 split.")
}

# Main script creates the manifest if absent. The ref0 script must read the
# manifest created by the main script; it must never create a second split.
if (file.exists(SPLIT_MANIFEST_FILE)) {

  split_manifest <- read.csv(
    SPLIT_MANIFEST_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(
      participant_id = "character",
      split = "character",
      id_column = "character",
      split_seed = "integer",
      training_fraction = "numeric"
    )
  )

} else {

  if (!SPLIT_CREATE_IF_MISSING) {
    stop(
      "STOP: Shared 70/30 split manifest was not found: ",
      SPLIT_MANIFEST_FILE,
      ". Run Step_03_1_1 (main script) first to create it."
    )
  }

  set.seed(SPLIT_SEED)
  training_ids_created <- sample(
    eligible_ids,
    size = n_train_expected,
    replace = FALSE
  )

  split_manifest <- data.frame(
    participant_id = eligible_ids,
    split = ifelse(
      eligible_ids %in% training_ids_created,
      "training",
      "test"
    ),
    id_column = ID_COL,
    split_seed = SPLIT_SEED,
    training_fraction = TRAIN_PROP,
    stringsAsFactors = FALSE
  )

  write.csv(
    split_manifest,
    file = SPLIT_MANIFEST_FILE,
    row.names = FALSE,
    na = ""
  )
}

# Strict manifest validation: any cohort mismatch stops the script.
required_manifest_cols <- c(
  "participant_id", "split", "id_column",
  "split_seed", "training_fraction"
)
if (!all(required_manifest_cols %in% names(split_manifest))) {
  stop("STOP: The shared split manifest has an unexpected structure.")
}
if (anyNA(split_manifest$participant_id) ||
    any(trimws(split_manifest$participant_id) == "")) {
  stop("STOP: Missing participant IDs in the shared split manifest.")
}
if (anyDuplicated(split_manifest$participant_id) > 0L) {
  stop("STOP: Duplicate participant IDs in the shared split manifest.")
}
if (!setequal(split_manifest$participant_id, eligible_ids)) {
  missing_from_manifest <- length(setdiff(eligible_ids, split_manifest$participant_id))
  extra_in_manifest <- length(setdiff(split_manifest$participant_id, eligible_ids))
  stop(
    "STOP: Current eligible cohort does not exactly match the shared manifest. ",
    "IDs missing from manifest = ", missing_from_manifest,
    "; extra IDs in manifest = ", extra_in_manifest,
    ". Do not proceed with inconsistent cohorts."
  )
}
if (anyNA(split_manifest$id_column) ||
    any(split_manifest$id_column != ID_COL)) {
  stop("STOP: Participant ID column differs between the two scripts.")
}
if (anyNA(split_manifest$split_seed) ||
    any(split_manifest$split_seed != SPLIT_SEED)) {
  stop("STOP: Split seed differs from the frozen shared manifest.")
}
if (anyNA(split_manifest$training_fraction) ||
    any(abs(split_manifest$training_fraction - TRAIN_PROP) > 1e-12)) {
  stop("STOP: Training fraction differs from the frozen shared manifest.")
}
if (!all(split_manifest$split %in% c("training", "test"))) {
  stop("STOP: The manifest contains invalid split labels.")
}
if (sum(split_manifest$split == "training") != n_train_expected ||
    sum(split_manifest$split == "test") != n_test_expected) {
  stop("STOP: Training/test counts in the manifest do not match the 70/30 rule.")
}

# Align the manifest to the current eligible cohort's sorted ID order.
split_manifest <- split_manifest[
  match(eligible_ids, split_manifest$participant_id),
  ,
  drop = FALSE
]

training_ids <- split_manifest$participant_id[
  split_manifest$split == "training"
]
test_ids <- split_manifest$participant_id[
  split_manifest$split == "test"
]

# Shared ID lists: these are identical for the main and ref0 scripts.
write.csv(
  data.frame(participant_id = training_ids),
  file = file.path(SPLIT_DIR, "ukb_70_30_training_ids.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(participant_id = test_ids),
  file = file.path(SPLIT_DIR, "ukb_70_30_test_ids.csv"),
  row.names = FALSE
)

message(
  "70/30 manifest validated. Eligible N = ", n_eligible,
  "; training N = ", length(training_ids),
  "; test N = ", length(test_ids),
  "; ID column = ", ID_COL,
  "; manifest = ", normalizePath(SPLIT_MANIFEST_FILE, mustWork = FALSE)
)

## -------------------- Step 3 factor转化 --------------------

### -------------------- Step 3.1 berries --------------------
ukb_data_depression %>% count(berries_group_range)
levels_berries <- {
  vals <- sort(unique(ukb_data_depression$berries_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(berries_group_range = factor(
    berries_group_range,
    levels = levels_berries,
    ordered = FALSE
  ))

### -------------------- Step 3.2 citrus_fruits --------------------
ukb_data_depression %>% count(citrus_fruits_group_range)
levels_citrus_fruits <- {
  vals <- sort(unique(ukb_data_depression$citrus_fruits_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(citrus_fruits_group_range = factor(
    citrus_fruits_group_range,
    levels = levels_citrus_fruits,
    ordered = FALSE
  ))

### -------------------- Step 3.3 other_fruits --------------------
ukb_data_depression %>% count(other_fruits_group_range)
levels_other_fruits <- {
  vals <- sort(unique(ukb_data_depression$other_fruits_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(other_fruits_group_range = factor(
    other_fruits_group_range,
    levels = levels_other_fruits,
    ordered = FALSE
  ))

### -------------------- Step 3.4 green_leafy_veg --------------------
ukb_data_depression %>% count(green_leafy_veg_group_range)
levels_green_leafy_veg <- {
  vals <- sort(unique(ukb_data_depression$green_leafy_veg_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(green_leafy_veg_group_range = factor(
    green_leafy_veg_group_range,
    levels = levels_green_leafy_veg,
    ordered = FALSE
  ))

### -------------------- Step 3.5 red_orange_veg --------------------
ukb_data_depression %>% count(red_orange_veg_group_range)
levels_red_orange_veg <- {
  vals <- sort(unique(ukb_data_depression$red_orange_veg_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(red_orange_veg_group_range = factor(
    red_orange_veg_group_range,
    levels = levels_red_orange_veg,
    ordered = FALSE
  ))

### -------------------- Step 3.6 other_veg --------------------
ukb_data_depression %>% count(other_veg_group_range)
levels_other_veg <- {
  vals <- sort(unique(ukb_data_depression$other_veg_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(other_veg_group_range = factor(
    other_veg_group_range,
    levels = levels_other_veg,
    ordered = FALSE
  ))

### -------------------- Step 3.7 potatoes --------------------
ukb_data_depression %>% count(potatoes_group_range)
levels_potatoes <- {
  vals <- sort(unique(ukb_data_depression$potatoes_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(potatoes_group_range = factor(
    potatoes_group_range,
    levels = levels_potatoes,
    ordered = FALSE
  ))

### -------------------- Step 3.8 refined_grains --------------------
ukb_data_depression %>% count(refined_grains_group_range)
levels_refined_grains <- {
  vals <- sort(unique(ukb_data_depression$refined_grains_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(refined_grains_group_range = factor(
    refined_grains_group_range,
    levels = levels_refined_grains,
    ordered = FALSE
  ))

### -------------------- Step 3.9 whole_grains --------------------
ukb_data_depression %>% count(whole_grains_group_range)
levels_whole_grains <- {
  vals <- sort(unique(ukb_data_depression$whole_grains_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(whole_grains_group_range = factor(
    whole_grains_group_range,
    levels = levels_whole_grains,
    ordered = FALSE
  ))

### -------------------- Step 3.10 legumes --------------------
ukb_data_depression %>% count(legumes_group_range)
levels_legumes <- {
  vals <- sort(unique(ukb_data_depression$legumes_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(legumes_group_range = factor(
    legumes_group_range,
    levels = levels_legumes,
    ordered = FALSE
  ))

### -------------------- Step 3.11 nuts --------------------
ukb_data_depression %>% count(nuts_group_range)
levels_nuts <- {
  vals <- sort(unique(ukb_data_depression$nuts_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(nuts_group_range = factor(
    nuts_group_range,
    levels = levels_nuts,
    ordered = FALSE
  ))

### -------------------- Step 3.12 eggs --------------------
ukb_data_depression %>% count(eggs_group_range)
levels_eggs <- {
  vals <- sort(unique(ukb_data_depression$eggs_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(eggs_group_range = factor(
    eggs_group_range,
    levels = levels_eggs,
    ordered = FALSE
  ))

### -------------------- Step 3.13 high_fat_dairy --------------------
ukb_data_depression %>% count(high_fat_dairy_group_range)
levels_high_fat_dairy <- {
  vals <- sort(unique(ukb_data_depression$high_fat_dairy_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(high_fat_dairy_group_range = factor(
    high_fat_dairy_group_range,
    levels = levels_high_fat_dairy,
    ordered = FALSE
  ))

### -------------------- Step 3.14 low_fat_dairy --------------------
ukb_data_depression %>% count(low_fat_dairy_group_range)
levels_low_fat_dairy <- {
  vals <- sort(unique(ukb_data_depression$low_fat_dairy_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(low_fat_dairy_group_range = factor(
    low_fat_dairy_group_range,
    levels = levels_low_fat_dairy,
    ordered = FALSE
  ))

### -------------------- Step 3.15 red_meats --------------------
ukb_data_depression %>% count(red_meats_group_range)
levels_red_meats <- {
  vals <- sort(unique(ukb_data_depression$red_meats_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(red_meats_group_range = factor(
    red_meats_group_range,
    levels = levels_red_meats,
    ordered = FALSE
  ))

### -------------------- Step 3.16 processed_meats --------------------
ukb_data_depression %>% count(processed_meats_group_range)
levels_processed_meats <- {
  vals <- sort(unique(ukb_data_depression$processed_meats_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(processed_meats_group_range = factor(
    processed_meats_group_range,
    levels = levels_processed_meats,
    ordered = FALSE
  ))

### -------------------- Step 3.17 poultry --------------------
ukb_data_depression %>% count(poultry_group_range)
levels_poultry <- {
  vals <- sort(unique(ukb_data_depression$poultry_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(poultry_group_range = factor(
    poultry_group_range,
    levels = levels_poultry,
    ordered = FALSE
  ))

### -------------------- Step 3.18 fish_seafood --------------------
ukb_data_depression %>% count(fish_seafood_group_range)
levels_fish_seafood <- {
  vals <- sort(unique(ukb_data_depression$fish_seafood_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(fish_seafood_group_range = factor(
    fish_seafood_group_range,
    levels = levels_fish_seafood,
    ordered = FALSE
  ))

### -------------------- Step 3.19 organ_meats --------------------
ukb_data_depression %>% count(organ_meats_group_range)
levels_organ_meats <- {
  vals <- sort(unique(ukb_data_depression$organ_meats_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(organ_meats_group_range = factor(
    organ_meats_group_range,
    levels = levels_organ_meats,
    ordered = FALSE
  ))

### -------------------- Step 3.20 butter --------------------
ukb_data_depression %>% count(butter_group_range)
levels_butter <- {
  vals <- sort(unique(ukb_data_depression$butter_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(butter_group_range = factor(
    butter_group_range,
    levels = levels_butter,
    ordered = FALSE
  ))

### -------------------- Step 3.21 margarine -------------------
ukb_data_depression %>% count(margarine_group_range)
levels_margarine <- {
  vals <- sort(unique(ukb_data_depression$margarine_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(margarine_group_range = factor(
    margarine_group_range,
    levels = levels_margarine,
    ordered = FALSE
  ))

### -------------------- Step 3.22 olive_oil -------------------
ukb_data_depression %>% count(olive_oil_group_range)
levels_olive_oil <- {
  vals <- sort(unique(ukb_data_depression$olive_oil_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(olive_oil_group_range = factor(
    olive_oil_group_range,
    levels = levels_olive_oil,
    ordered = FALSE
  ))

### -------------------- Step 3.23 other_oil -------------------
ukb_data_depression %>% count(other_oil_group_range)
levels_other_oil <- {
  vals <- sort(unique(ukb_data_depression$other_oil_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(other_oil_group_range = factor(
    other_oil_group_range,
    levels = levels_other_oil,
    ordered = FALSE
  ))

### -------------------- Step 3.24 beer -------------------
ukb_data_depression %>% count(beer_group_range)
levels_beer <- {
  vals <- sort(unique(ukb_data_depression$beer_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(beer_group_range = factor(
    beer_group_range,
    levels = levels_beer,
    ordered = FALSE
  ))


### -------------------- Step 3.25 coffee -------------------
ukb_data_depression %>% count(coffee_group_range)

levels_coffee <- {
  vals <- sort(unique(ukb_data_depression$coffee_group_range))
  c(vals[2], vals[-2])   # 这里是你的特殊顺序调整
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(coffee_group_range = factor(
    coffee_group_range,
    levels = levels_coffee,
    ordered = FALSE
  ))

### -------------------- Step 3.26 fruit_juices -------------------
ukb_data_depression %>% count(fruit_juices_group_range)
levels_fruit_juices <- {
  vals <- sort(unique(ukb_data_depression$fruit_juices_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(fruit_juices_group_range = factor(
    fruit_juices_group_range,
    levels = levels_fruit_juices,
    ordered = FALSE
  ))

### -------------------- Step 3.27 liquor -------------------
ukb_data_depression %>% count(liquor_group_range)
levels_liquor <- {
  vals <- sort(unique(ukb_data_depression$liquor_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(liquor_group_range = factor(
    liquor_group_range,
    levels = levels_liquor,
    ordered = FALSE
  ))

### -------------------- Step 3.28 sweetened_bev -------------------
ukb_data_depression %>% count(sweetened_bev_group_range)
levels_sweetened_bev <- {
  vals <- sort(unique(ukb_data_depression$sweetened_bev_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(sweetened_bev_group_range = factor(
    sweetened_bev_group_range,
    levels = levels_sweetened_bev,
    ordered = FALSE
  ))

### -------------------- Step 3.29 tea -------------------
ukb_data_depression %>% count(tea_group_range)
levels_tea <- {
  vals <- sort(unique(ukb_data_depression$tea_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(tea_group_range = factor(
    tea_group_range,
    levels = levels_tea,
    ordered = FALSE
  ))

### -------------------- Step 3.30 water -------------------
ukb_data_depression %>% count(water_group_range)
levels_water <- {
  vals <- sort(unique(ukb_data_depression$water_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(water_group_range = factor(
    water_group_range,
    levels = levels_water,
    ordered = FALSE
  ))

### -------------------- Step 3.31 wine -------------------
ukb_data_depression %>% count(wine_group_range)
levels_wine <- {
  vals <- sort(unique(ukb_data_depression$wine_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(wine_group_range = factor(
    wine_group_range,
    levels = levels_wine,
    ordered = FALSE
  ))

### -------------------- Step 3.32 pizza -------------------
ukb_data_depression %>% count(pizza_group_range)
levels_pizza <- {
  vals <- sort(unique(ukb_data_depression$pizza_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(pizza_group_range = factor(
    pizza_group_range,
    levels = levels_pizza,
    ordered = FALSE
  ))

### -------------------- Step 3.33 snacks_pastries -------------------
ukb_data_depression %>% count(snacks_pastries_group_range)
levels_snacks_pastries <- {
  vals <- sort(unique(ukb_data_depression$snacks_pastries_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(snacks_pastries_group_range = factor(
    snacks_pastries_group_range,
    levels = levels_snacks_pastries,
    ordered = FALSE
  ))

### -------------------- Step 3.34 fried_fast_foods -------------------
ukb_data_depression %>% count(fried_fast_foods_group_range)
levels_fried_fast_foods <- {
  vals <- sort(unique(ukb_data_depression$fried_fast_foods_group_range))
  c(vals[2], vals[-2])
}

ukb_data_depression <- ukb_data_depression %>%
  mutate(fried_fast_foods_group_range = factor(
    fried_fast_foods_group_range,
    levels = levels_fried_fast_foods,
    ordered = FALSE
  ))



## -------------------- Step 3.35 Apply the shared 70/30 split --------------------
# Step 3 factor conversion above only formats exposure categories and does not
# use depression associations or outcome-derived selection. It is retained on
# the full eligible cohort to preserve the original category/reference coding.
# All outcome-associated food-group models begin in Step 4 and will use TRAINING
# only, because the existing model calls below continue to use
# ukb_data_depression as their data object.

SCRIPT_SPLIT_DIR <- file.path(SPLIT_DIR, SPLIT_SCRIPT_TAG)
dir.create(SCRIPT_SPLIT_DIR, recursive = TRUE, showWarnings = FALSE)

ukb_data_depression_full <- ukb_data_depression

id_values_factorized <- trimws(
  as.character(ukb_data_depression_full[[ID_COL]])
)

if (anyNA(id_values_factorized) ||
    any(id_values_factorized == "") ||
    anyDuplicated(id_values_factorized) > 0L) {
  stop("STOP: Participant ID QC failed after Step 3 factor conversion.")
}

split_lookup <- setNames(
  split_manifest$split,
  split_manifest$participant_id
)

ukb_data_depression_full$analysis_split <- unname(
  split_lookup[id_values_factorized]
)

if (anyNA(ukb_data_depression_full$analysis_split)) {
  stop("STOP: Some eligible participants could not be assigned to training/test.")
}

ukb_data_depression_train <- ukb_data_depression_full %>%
  filter(analysis_split == "training")

ukb_data_depression_test <- ukb_data_depression_full %>%
  filter(analysis_split == "test")

# Verify exact ID identity against the shared manifest.
train_ids_observed <- sort(
  as.character(ukb_data_depression_train[[ID_COL]]),
  method = "radix"
)
test_ids_observed <- sort(
  as.character(ukb_data_depression_test[[ID_COL]]),
  method = "radix"
)

if (!identical(train_ids_observed, sort(training_ids, method = "radix"))) {
  stop("STOP: Training participant IDs do not match the shared manifest.")
}
if (!identical(test_ids_observed, sort(test_ids, method = "radix"))) {
  stop("STOP: Test participant IDs do not match the shared manifest.")
}
if (length(intersect(train_ids_observed, test_ids_observed)) != 0L) {
  stop("STOP: Participant overlap detected between training and test sets.")
}
if (!identical(
  sort(c(train_ids_observed, test_ids_observed), method = "radix"),
  eligible_ids
)) {
  stop("STOP: Training/test union does not reconstruct the eligible cohort.")
}

summarise_depression_split <- function(dat, split_label) {
  tibble(
    split = split_label,
    participants = nrow(dat),
    depression_events = sum(dat$depression_status == 1, na.rm = TRUE),
    non_events = sum(dat$depression_status == 0, na.rm = TRUE),
    missing_depression_status = sum(is.na(dat$depression_status))
  )
}

split_cohort_summary <- bind_rows(
  summarise_depression_split(ukb_data_depression_full, "full_cohort"),
  summarise_depression_split(ukb_data_depression_train, "training"),
  summarise_depression_split(ukb_data_depression_test, "test")
) %>%
  mutate(
    event_percent = ifelse(
      participants > 0,
      100 * depression_events / participants,
      NA_real_
    )
  )

write.csv(
  split_cohort_summary,
  file = file.path(SPLIT_DIR, "ukb_70_30_split_cohort_summary.csv"),
  row.names = FALSE
)

# Save factorized datasets separately for the main and ref0 coding variants.
# The participant allocation is shared; factor reference coding remains as in
# each original script.
saveRDS(
  ukb_data_depression_train,
  file = file.path(SCRIPT_SPLIT_DIR, "ukb_depression_training_cohort.rds")
)
saveRDS(
  ukb_data_depression_test,
  file = file.path(SCRIPT_SPLIT_DIR, "ukb_depression_test_cohort.rds")
)

manifest_md5 <- unname(tools::md5sum(SPLIT_MANIFEST_FILE))

split_run_qc <- data.frame(
  script_tag = SPLIT_SCRIPT_TAG,
  id_column = ID_COL,
  split_seed = SPLIT_SEED,
  training_fraction = TRAIN_PROP,
  eligible_n = nrow(ukb_data_depression_full),
  training_n = nrow(ukb_data_depression_train),
  test_n = nrow(ukb_data_depression_test),
  eligible_events = sum(
    ukb_data_depression_full$depression_status == 1,
    na.rm = TRUE
  ),
  training_events = sum(
    ukb_data_depression_train$depression_status == 1,
    na.rm = TRUE
  ),
  test_events = sum(
    ukb_data_depression_test$depression_status == 1,
    na.rm = TRUE
  ),
  train_ids_match_manifest = TRUE,
  test_ids_match_manifest = TRUE,
  train_test_overlap_n = 0L,
  manifest_md5 = manifest_md5,
  stringsAsFactors = FALSE
)

write.csv(
  split_run_qc,
  file = file.path(SCRIPT_SPLIT_DIR, "split_run_qc.csv"),
  row.names = FALSE
)

cat("\n===== 70/30 SPLIT QC =====\n")
print(split_cohort_summary)
print(split_run_qc)

# CRITICAL: preserve all original Step 4/5 model calls and route them to the
# training set simply by replacing this working object. The test set remains
# available as ukb_data_depression_test and is saved for later frozen-score
# internal validation; it is NOT used for food-group screening in this script.
ukb_data_depression <- ukb_data_depression_train

## -------------------- Step 4 food group Model 1 --------------------
paste(colnames(ukb_data), collapse = " + ")

# =============================================================================
# REFERENCE ROTATION EXTENSION (training cohort only)
# Keep the original Cox covariates, categorical intake groups, trend definitions,
# forest-plot row order, column layout, colors, ticks and PDF dimensions.
# Each observed intake category is treated as the reference once. The final
# display uses the lowest adjusted relative-risk category as reference.
# Lower category relative risk is a DATA-DERIVED estimate, not proof of an
# optimal intake or a causal protection. This does NOT change the RCS model.
# =============================================================================
ROTATE_RESULTS_DIR <- "./results/Step_03/reference_rotation_lowest_risk"
dir.create(ROTATE_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

ROTATION_RECORDS <- new.env(parent = emptyenv())
ROTATION_RECORDS$selected <- list()
ROTATION_RECORDS$all_reference_comparisons <- list()
ROTATION_RECORDS$diagnostics <- list()
ROTATION_RECORDS$warnings <- list()

rotation_format_p <- function(x) {
  if (length(x) != 1L || is.na(x)) return(NA_character_)
  if (x < 0.001) return("< 0.001")
  # Avoid displaying a truly significant P = 0.0496 as a misleading "0.050".
  if (x >= 0.0495 && x < 0.0505) return(formatC(x, format = "f", digits = 4))
  formatC(x, format = "f", digits = 3)
}

rotation_parse_p <- function(x) {
  ans <- suppressWarnings(as.numeric(as.character(x)))
  lower_label <- grepl("<", as.character(x), fixed = TRUE)
  ans[lower_label & !is.na(lower_label)] <- 0.0005
  ans
}

rotation_cox_fit <- function(data, formula, fit_tag = "unspecified") {
  # The original calls use survival::coxph with the default Efron ties method.
  warnings_seen <- character(0)
  fit <- withCallingHandlers(
    survival::coxph(formula, data = data, ties = "efron"),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      # Keep the warning visible in the run log; do not silently suppress it.
    }
  )
  if (length(warnings_seen) > 0L) {
    ROTATION_RECORDS$warnings[[length(ROTATION_RECORDS$warnings) + 1L]] <-
      data.frame(
        fit_tag = fit_tag, warning = warnings_seen,
        stringsAsFactors = FALSE
      )
  }
  fit
}

rotation_log_risks <- function(fit, exposure_var, exposure_levels, reference) {
  b <- stats::coef(fit)
  other <- exposure_levels[exposure_levels != reference]
  term <- paste0(exposure_var, other)
  if (!all(term %in% names(b))) {
    stop("STOP: Exposure coefficients missing for ", exposure_var,
         " (reference ", reference, "). Expected: ",
         paste(setdiff(term, names(b)), collapse = "; "))
  }
  vals <- setNames(rep(0, length(exposure_levels)), exposure_levels)
  vals[other] <- unname(b[term])
  if (any(!is.finite(vals))) {
    stop("STOP: Non-finite categorical log-HR values: ", exposure_var)
  }
  vals
}

rotation_fit_one_food <- function(data, var, var_num, group_name,
                                  ref_level, levels_order, variable_name,
                                  model_name, covariate_text) {
  # Do not use 30% test data. The calling script has already replaced
  # ukb_data_depression with the frozen training-cohort object.
  if (!is.data.frame(data) || !is.factor(data[[var]])) {
    stop("STOP: Expected factor exposure for ", var)
  }
  if (!var_num %in% names(data)) {
    stop("STOP: Missing original trend variable: ", var_num)
  }
  cov_names <- trimws(strsplit(covariate_text, "+", fixed = TRUE)[[1]])
  need <- unique(c("depression_time", "depression_status", var, cov_names))
  missing <- setdiff(need, names(data))
  if (length(missing) > 0L) {
    stop("STOP: Model variables missing for ", var, ": ",
         paste(missing, collapse = ", "))
  }

  # Fixed complete-case analytical sample ensures all rotated comparisons
  # use the same participants, covariates and number of depression events.
  keep <- stats::complete.cases(data[, need, drop = FALSE])
  d <- data[keep, , drop = FALSE]
  cat_vals <- as.character(d[[var]])
  observed_levels <- as.character(levels_order)
  observed_levels <- observed_levels[observed_levels %in% cat_vals]
  observed_levels <- unique(observed_levels)
  if (length(observed_levels) < 2L ||
      any(!cat_vals %in% observed_levels)) {
    stop("STOP: Insufficient or unmapped intake groups for ", var)
  }
  d[[var]] <- factor(cat_vals, levels = observed_levels, ordered = FALSE)
  formula_categorical <- stats::as.formula(paste0(
    "survival::Surv(depression_time, depression_status) ~ ",
    var, " + ", covariate_text
  ))
  formula_trend <- stats::as.formula(paste0(
    "survival::Surv(depression_time, depression_status) ~ as.numeric(",
    var_num, ") + ", covariate_text
  ))

  # First category is only a numerical anchor for finding the category with
  # minimum estimated adjusted risk. All references are then actually refitted.
  original_reference <- observed_levels[1]
  anchor_fit <- rotation_cox_fit(d, formula_categorical,
                                 paste(model_name, var, "initial", sep = "__"))
  anchor_log_hr <- rotation_log_risks(
    anchor_fit, var, observed_levels, original_reference
  )
  selected_reference <- observed_levels[which.min(anchor_log_hr)]

  fit_by_reference <- setNames(vector("list", length(observed_levels)), observed_levels)
  comparison_rows <- list()
  diagnostics_rows <- list()
  k <- 0L
  for (current_ref in observed_levels) {
    d[[var]] <- stats::relevel(d[[var]], ref = current_ref)
    fit <- rotation_cox_fit(d, formula_categorical,
                            paste(model_name, var, current_ref, sep = "__"))
    cf <- as.matrix(summary(fit)$coefficients)
    ci <- as.matrix(summary(fit)$conf.int)
    this_log_hr <- rotation_log_risks(fit, var, observed_levels, current_ref)
    if (!identical(as.integer(fit$n), as.integer(anchor_fit$n)) ||
        !identical(as.integer(fit$nevent), as.integer(anchor_fit$nevent))) {
      stop("STOP: Rotating reference changed model N/events for ", var)
    }

    # Reference invariance: the same categorical model must yield identical
    # pairwise fitted-risk contrasts irrespective of the baseline category.
    expected <- anchor_log_hr - unname(anchor_log_hr[current_ref])
    max_delta <- max(abs(this_log_hr - expected))
    if (!is.finite(max_delta) || max_delta > 1e-4) {
      stop("STOP: Reference rotation changed adjusted log-HR contrasts for ",
           var, " at ", current_ref, "; max delta=", signif(max_delta, 6))
    }
    fit_by_reference[[current_ref]] <- fit
    diagnostics_rows[[current_ref]] <- data.frame(
      model = model_name, group = group_name, food_variable = var,
      reference_category = current_ref, model_n = fit$n,
      model_events = fit$nevent, max_contrast_delta = max_delta,
      loglik_full = unname(fit$loglik[2]),
      stringsAsFactors = FALSE
    )
    for (lv in observed_levels) {
      k <- k + 1L
      if (lv == current_ref) {
        hr <- 1; low <- 1; high <- 1; p <- NA_real_; beta <- 0
      } else {
        term <- paste0(var, lv)
        if (!term %in% rownames(cf) || !term %in% rownames(ci)) {
          stop("STOP: Cannot map rotated coefficient for ", var,
               ", reference=", current_ref, ", level=", lv)
        }
        beta <- as.numeric(cf[term, "coef"])
        hr <- exp(beta)
        low <- as.numeric(ci[term, "lower .95"])
        high <- as.numeric(ci[term, "upper .95"])
        p <- as.numeric(cf[term, "Pr(>|z|)"])
        if (any(!is.finite(c(beta, hr, low, high, p)))) {
          stop("STOP: Invalid Cox coefficient or CI: ", var, " / ", current_ref,
               " vs ", lv)
        }
      }
      comparison_rows[[k]] <- data.frame(
        model = model_name, group = group_name, food_variable = var,
        reference_category = current_ref, comparison_category = lv,
        is_selected_lowest_risk_reference = current_ref == selected_reference,
        is_reference_row = lv == current_ref,
        beta = beta, hazard_ratio = hr, ci_low = low, ci_high = high,
        p_value_raw = p, hr_ci = if (lv == current_ref) "Ref." else
          sprintf("%.2f (%.2f, %.2f)", hr, low, high),
        model_n = fit$n, model_events = fit$nevent,
        stringsAsFactors = FALSE
      )
    }
  }
  fit_selected <- fit_by_reference[[selected_reference]]
  selection_log_hr <- anchor_log_hr - unname(anchor_log_hr[selected_reference])
  if (!all(selection_log_hr >= -1e-6)) {
    stop("STOP: Lowest-risk category selection mismatch for ", var)
  }

  # Keep the original trend model exactly as defined in the user's script.
  # The trend predictor and all adjustment covariates remain unchanged.
  fit_trend <- rotation_cox_fit(data, formula_trend,
                                  paste(model_name, var, "trend", sep = "__"))
  trend_coef <- as.matrix(summary(fit_trend)$coefficients)
  trend_term <- paste0("as.numeric(", var_num, ")")
  if (!trend_term %in% rownames(trend_coef)) {
    stop("STOP: Cannot identify original P-for-trend coefficient for ", var)
  }
  trend_p <- as.numeric(trend_coef[trend_term, "Pr(>|z|)"])

  # Original figure row ordering and column schema: unchanged.
  display_rows <- lapply(observed_levels, function(lv) {
    if (lv == selected_reference) {
      hr <- 1; low <- 1; high <- 1; p_raw <- NA_real_
      ci_label <- "Ref."; p_label <- "Ref."
    } else {
      key <- paste0(var, lv)
      tab <- as.matrix(summary(fit_selected)$coefficients)
      ci <- as.matrix(summary(fit_selected)$conf.int)
      hr <- exp(unname(tab[key, "coef"]))
      low <- unname(ci[key, "lower .95"])
      high <- unname(ci[key, "upper .95"])
      p_raw <- unname(tab[key, "Pr(>|z|)"])
      ci_label <- sprintf("%.2f (%.2f, %.2f)", hr, low, high)
      p_label <- rotation_format_p(p_raw)
    }
    data.frame(
      group = group_name, model = model_name, variable = lv,
      hr_ci = ci_label, p_value = p_label,
      hazard_ratio = hr, ci_low = low, ci_high = high,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, display_rows)

  # Match the original forestplot Events/Total (entire frozen training cohort)
  # while recording the complete-case Cox model n and events separately.
  counts <- data.frame(
    variable = observed_levels,
    Events_Total = vapply(observed_levels, function(lv) {
      included <- !is.na(data[[var]]) & as.character(data[[var]]) == lv
      sprintf("%s/%s",
        formatC(sum(data$depression_status[included] == 1, na.rm = TRUE),
                format = "d", big.mark = ","),
        formatC(sum(included), format = "d", big.mark = ","))
    }, character(1)),
    stringsAsFactors = FALSE
  )
  result <- merge(result, counts, by = "variable", sort = FALSE, all.x = TRUE)
  result <- result[match(observed_levels, result$variable), , drop = FALSE]
  trend_row <- data.frame(
    group = group_name, model = model_name, variable = "P for trend",
    hr_ci = NA_character_, p_value = rotation_format_p(trend_p),
    hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    Events_Total = NA_character_, stringsAsFactors = FALSE
  )
  result <- rbind(result, trend_row)
  result <- result[, c("group", "model", "variable", "Events_Total",
                       "hr_ci", "p_value", "hazard_ratio", "ci_low", "ci_high")]

  key <- paste(model_name, var, sep = "__")
  ROTATION_RECORDS$selected[[key]] <- data.frame(
    model = model_name, group = group_name, food_variable = var,
    original_reference_category = as.character(ref_level),
    lowest_risk_reference_category = selected_reference,
    reference_changed = selected_reference != as.character(ref_level),
    n_categories = length(observed_levels),
    model_n = fit_selected$n, model_events = fit_selected$nevent,
    original_trend_p_raw = trend_p,
    chosen_reference_is_boundary_group =
      selected_reference %in% c(observed_levels[1], tail(observed_levels, 1)),
    stringsAsFactors = FALSE
  )
  ROTATION_RECORDS$all_reference_comparisons[[key]] <-
    do.call(rbind, comparison_rows)
  ROTATION_RECORDS$diagnostics[[key]] <- do.call(rbind, diagnostics_rows)
  result
}
# =============================================================================

### -------------------- Step 4.0 模型定义与函数封装 --------------------

# Model 1 adjusted for age_food + sex + ethnicity + education_attainment + twnsend + career + smoking_status + drink_status + bmi_group + physical_activity_proup + sedentary_group + household_income + energy_intake

run_food_analysis_model1 <- function(data, var, var_num, group_name, ref_level,
                                     levels_order, variable_name = NULL,
                                     model_name = "model1") {
  rotation_fit_one_food(
    data = data, var = var, var_num = var_num, group_name = group_name,
    ref_level = ref_level, levels_order = levels_order,
    variable_name = variable_name, model_name = model_name,
    covariate_text = "age_food + sex + ethnicity + education_attainment + twnsend + career + smoking_status + drink_status + bmi_group + physical_activity_proup + sedentary_group + household_income + energy_intake"
  )
}

run_food_analysis_two_level_model1 <- function(data, var, var_num, group_name, ref_level,
                                     levels_order, variable_name = NULL,
                                     model_name = "model1") {
  rotation_fit_one_food(
    data = data, var = var, var_num = var_num, group_name = group_name,
    ref_level = ref_level, levels_order = levels_order,
    variable_name = variable_name, model_name = model_name,
    covariate_text = "age_food + sex + ethnicity + education_attainment + twnsend + career + smoking_status + drink_status + bmi_group + physical_activity_proup + sedentary_group + household_income + energy_intake"
  )
}

### -------------------- Step 4.1 berries --------------------
levels_berries_ordered <- {
  vals <- levels_berries
  c(vals[2], vals[-2])
}


berries_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "berries_group_range",
  var_num = "berries_group_range_num",
  group_name = "Berries",
  ref_level = levels_berries_ordered[2],
  levels_order = levels_berries_ordered,
  model_name = "model1"
)

berries_model1_result

### -------------------- Step 4.2 citrus_fruits --------------------
levels_citrus_fruits_ordered <- {
  vals <- levels_citrus_fruits
  c(vals[2], vals[-2])
}

citrus_fruits_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "citrus_fruits_group_range",
  var_num = "citrus_fruits_group_range_num",
  group_name = "Citrus fruits",
  ref_level = levels_citrus_fruits_ordered[2],
  levels_order = levels_citrus_fruits_ordered,
  model_name = "model1"
)

### -------------------- Step 4.3 other_fruits --------------------
levels_other_fruits_ordered <- {
  vals <- levels_other_fruits
  c(vals[2], vals[-2])
}

other_fruits_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "other_fruits_group_range",
  var_num = "other_fruits_group_range_num",
  group_name = "Other fruits",
  ref_level = levels_other_fruits_ordered[2],
  levels_order = levels_other_fruits_ordered,
  model_name = "model1"
)

### -------------------- Step 4.4 green_leafy_veg --------------------
levels_green_leafy_veg_ordered <- {
  vals <- levels_green_leafy_veg
  c(vals[2], vals[-2])
}

green_leafy_veg_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "green_leafy_veg_group_range",
  var_num = "green_leafy_veg_group_range_num",
  group_name = "Green leafy vegetables",
  ref_level = levels_green_leafy_veg_ordered[2],
  levels_order = levels_green_leafy_veg_ordered,
  model_name = "model1"
)

### -------------------- Step 4.5 red_orange_veg --------------------
levels_red_orange_veg_ordered <- {
  vals <- levels_red_orange_veg
  c(vals[2], vals[-2])
}

red_orange_veg_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "red_orange_veg_group_range",
  var_num = "red_orange_veg_group_range_num",
  group_name = "Red/orange vegetables",
  ref_level = levels_red_orange_veg_ordered[2],
  levels_order = levels_red_orange_veg_ordered,
  model_name = "model1"
)

### -------------------- Step 4.6 other_veg --------------------
levels_other_veg_ordered <- {
  vals <- levels_other_veg
  c(vals[2], vals[-2])
}

other_veg_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "other_veg_group_range",
  var_num = "other_veg_group_range_num",
  group_name = "Other vegetables",
  ref_level = levels_other_veg_ordered[2],
  levels_order = levels_other_veg_ordered,
  model_name = "model1"
)

### -------------------- Step 4.7 potatoes --------------------
levels_potatoes_ordered <- {
  vals <- levels_potatoes
  c(vals[2], vals[-2])
}

potatoes_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "potatoes_group_range",
  var_num = "potatoes_group_range_num",
  group_name = "Potatoes",
  ref_level = levels_potatoes_ordered[2],
  levels_order = levels_potatoes_ordered,
  model_name = "model1"
)

### -------------------- Step 4.8 refined_grains --------------------
levels_refined_grains_ordered <- {
  vals <- levels_refined_grains
  c(vals[2], vals[-2])
}

refined_grains_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "refined_grains_group_range",
  var_num = "refined_grains_group_range_num",
  group_name = "Refined grains",
  ref_level = levels_refined_grains_ordered[2],
  levels_order = levels_refined_grains_ordered,
  model_name = "model1"
)

### -------------------- Step 4.9 whole_grains --------------------
levels_whole_grains_ordered <- {
  vals <- levels_whole_grains
  c(vals[2], vals[-2])
}

whole_grains_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "whole_grains_group_range",
  var_num = "whole_grains_group_range_num",
  group_name = "Whole grains",
  ref_level = levels_whole_grains_ordered[2],
  levels_order = levels_whole_grains_ordered,
  model_name = "model1"
)

### -------------------- Step 4.10 legumes --------------------
levels_legumes_ordered <- {
  vals <- levels_legumes
  c(vals[2], vals[-2])
}

legumes_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "legumes_group_range",
  var_num = "legumes_group_range_num",
  group_name = "Legumes",
  ref_level = levels_legumes_ordered[2],
  levels_order = levels_legumes_ordered,
  model_name = "model1"
)


### -------------------- Step 4.11 nuts --------------------
levels_nuts_ordered <- {
  vals <- levels_nuts
  c(vals[2], vals[-2])
}

nuts_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "nuts_group_range",
  var_num = "nuts_group_range_num",
  group_name = "Nuts",
  ref_level = levels_nuts_ordered[2],
  levels_order = levels_nuts_ordered,
  model_name = "model1"
)

### -------------------- Step 4.12 eggs --------------------
levels_eggs_ordered <- {
  vals <- levels_eggs
  c(vals[2], vals[-2])
}

eggs_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "eggs_group_range",
  var_num = "eggs_group_range_num",
  group_name = "Eggs",
  ref_level = levels_eggs_ordered[2],
  levels_order = levels_eggs_ordered,
  model_name = "model1"
)

### -------------------- Step 4.13 high_fat_dairy --------------------
levels_high_fat_dairy_ordered <- {
  vals <- levels_high_fat_dairy
  c(vals[2], vals[-2])
}

high_fat_dairy_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "high_fat_dairy_group_range",
  var_num = "high_fat_dairy_group_range_num",
  group_name = "High fat dairy",
  ref_level = levels_high_fat_dairy_ordered[2],
  levels_order = levels_high_fat_dairy_ordered,
  model_name = "model1"
)

### -------------------- Step 4.14 low_fat_dairy --------------------
levels_low_fat_dairy_ordered <- {
  vals <- levels_low_fat_dairy
  c(vals[2], vals[-2])
}

low_fat_dairy_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "low_fat_dairy_group_range",
  var_num = "low_fat_dairy_group_range_num",
  group_name = "Low fat dairy",
  ref_level = levels_low_fat_dairy_ordered[2],
  levels_order = levels_low_fat_dairy_ordered,
  model_name = "model1"
)

### -------------------- Step 4.15 red_meats --------------------
levels_red_meats_ordered <- {
  vals <- levels_red_meats
  c(vals[2], vals[-2])
}

red_meats_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "red_meats_group_range",
  var_num = "red_meats_group_range_num",
  group_name = "Red meats",
  ref_level = levels_red_meats_ordered[2],
  levels_order = levels_red_meats_ordered,
  model_name = "model1"
)

### -------------------- Step 4.16 processed_meats --------------------
levels_processed_meats_ordered <- {
  vals <- levels_processed_meats
  c(vals[2], vals[-2])
}

processed_meats_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "processed_meats_group_range",
  var_num = "processed_meats_group_range_num",
  group_name = "Processed meats",
  ref_level = levels_processed_meats_ordered[2],
  levels_order = levels_processed_meats_ordered,
  model_name = "model1"
)

### -------------------- Step 4.17 poultry --------------------
levels_poultry_ordered <- {
  vals <- levels_poultry
  c(vals[2], vals[-2])
}

poultry_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "poultry_group_range",
  var_num = "poultry_group_range_num",
  group_name = "Poultry",
  ref_level = levels_poultry_ordered[2],
  levels_order = levels_poultry_ordered,
  model_name = "model1"
)

### -------------------- Step 4.18 fish_seafood --------------------
levels_fish_seafood_ordered <- {
  vals <- levels_fish_seafood
  c(vals[2], vals[-2])
}

fish_seafood_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "fish_seafood_group_range",
  var_num = "fish_seafood_group_range_num",
  group_name = "Fish and seafood",
  ref_level = levels_fish_seafood_ordered[2],
  levels_order = levels_fish_seafood_ordered,
  model_name = "model1"
)

### -------------------- Step 4.19 organ_meats --------------------
levels_organ_meats_ordered <- {
  vals <- levels_organ_meats
  c(vals[2], vals[-2])
}

organ_meats_model1_result <- run_food_analysis_two_level_model1(
  data = ukb_data_depression,
  var = "organ_meats_group_range",
  var_num = "organ_meats_group_range_num",
  group_name = "Organ meats",
  ref_level = levels_organ_meats_ordered[2],
  variable_name = levels_organ_meats_ordered[1],
  levels_order = levels_organ_meats_ordered,
  model_name = "model1"
)

### -------------------- Step 4.20 butter --------------------
levels_butter_ordered <- {
  vals <- levels_butter
  c(vals[2], vals[-2])
}

butter_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "butter_group_range",
  var_num = "butter_group_range_num",
  group_name = "Butter",
  ref_level = levels_butter_ordered[2],
  levels_order = levels_butter_ordered,
  model_name = "model1"
)


### -------------------- Step 4.21 margarine --------------------
levels_margarine_ordered <- {
  vals <- levels_margarine
  c(vals[2], vals[-2])
}

margarine_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "margarine_group_range",
  var_num = "margarine_group_range_num",
  group_name = "Margarine",
  ref_level = levels_margarine_ordered[2],
  levels_order = levels_margarine_ordered,
  model_name = "model1"
)

### -------------------- Step 4.22 olive_oil --------------------
levels_olive_oil_ordered <- {
  vals <- levels_olive_oil
  c(vals[2], vals[-2])
}

olive_oil_model1_result <- run_food_analysis_two_level_model1(
  data = ukb_data_depression,
  var = "olive_oil_group_range",
  var_num = "olive_oil_group_range_num",
  group_name = "Olive oil",
  ref_level = levels_olive_oil_ordered[2],
  variable_name = levels_olive_oil_ordered[1],
  levels_order = levels_olive_oil_ordered,
  model_name = "model1"
)

### -------------------- Step 4.23 other_oil --------------------
levels_other_oil_ordered <- {
  vals <- levels_other_oil
  c(vals[2], vals[-2])
}

other_oil_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "other_oil_group_range",
  var_num = "other_oil_group_range_num",
  group_name = "Other oil",
  ref_level = levels_other_oil_ordered[2],
  levels_order = levels_other_oil_ordered,
  model_name = "model1"
)

### -------------------- Step 4.24 beer --------------------
levels_beer_ordered <- {
  vals <- levels_beer
  c(vals[2], vals[-2])
}

beer_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "beer_group_range",
  var_num = "beer_group_range_num",
  group_name = "Beer",
  ref_level = levels_beer_ordered[2],
  levels_order = levels_beer_ordered,
  model_name = "model1"
)

### -------------------- Step 4.25 coffee --------------------
levels_coffee_ordered <- {
  vals <- levels_coffee
  c(vals[2], vals[-2])
}

coffee_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "coffee_group_range",
  var_num = "coffee_group_range_num",
  group_name = "Coffee",
  ref_level = levels_coffee_ordered[2],
  levels_order = levels_coffee_ordered,
  model_name = "model1"
)

### -------------------- Step 4.26 fruit_juices --------------------
levels_fruit_juices_ordered <- {
  vals <- levels_fruit_juices
  c(vals[2], vals[-2])
}

fruit_juices_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "fruit_juices_group_range",
  var_num = "fruit_juices_group_range_num",
  group_name = "Fruit juices",
  ref_level = levels_fruit_juices_ordered[2],
  levels_order = levels_fruit_juices_ordered,
  model_name = "model1"
)

### -------------------- Step 4.27 liquor --------------------
levels_liquor_ordered <- {
  vals <- levels_liquor
  c(vals[2], vals[-2])
}

liquor_model1_result <- run_food_analysis_two_level_model1(
  data = ukb_data_depression,
  var = "liquor_group_range",
  var_num = "liquor_group_range_num",
  group_name = "Liquor",
  ref_level = levels_liquor_ordered[2],
  variable_name = levels_liquor_ordered[1],
  levels_order = levels_liquor_ordered,
  model_name = "model1"
)

### -------------------- Step 4.28 sweetened_bev --------------------
levels_sweetened_bev_ordered <- {
  vals <- levels_sweetened_bev
  c(vals[2], vals[-2])
}

sweetened_bev_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "sweetened_bev_group_range",
  var_num = "sweetened_bev_group_range_num",
  group_name = "Sweetened beverages",
  ref_level = levels_sweetened_bev_ordered[2],
  levels_order = levels_sweetened_bev_ordered,
  model_name = "model1"
)

### -------------------- Step 4.29 tea --------------------
levels_tea_ordered <- {
  vals <- levels_tea
  c(vals[2], vals[-2])
}

tea_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "tea_group_range",
  var_num = "tea_group_range_num",
  group_name = "Tea",
  ref_level = levels_tea_ordered[2],
  levels_order = levels_tea_ordered,
  model_name = "model1"
)

### -------------------- Step 4.30 water --------------------
levels_water_ordered <- {
  vals <- levels_water
  c(vals[2], vals[-2])
}

water_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "water_group_range",
  var_num = "water_group_range_num",
  group_name = "Water",
  ref_level = levels_water_ordered[2],
  levels_order = levels_water_ordered,
  model_name = "model1"
)

### -------------------- Step 4.31 wine --------------------
levels_wine_ordered <- {
  vals <- levels_wine
  c(vals[2], vals[-2])
}

wine_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "wine_group_range",
  var_num = "wine_group_range_num",
  group_name = "Wine",
  ref_level = levels_wine_ordered[2],
  levels_order = levels_wine_ordered,
  model_name = "model1"
)

### -------------------- Step 4.32 pizza --------------------
levels_pizza_ordered <- {
  vals <- levels_pizza
  c(vals[2], vals[-2])
}

pizza_model1_result <- run_food_analysis_two_level_model1(
  data = ukb_data_depression,
  var = "pizza_group_range",
  var_num = "pizza_group_range_num",
  group_name = "Pizza",
  ref_level = levels_pizza_ordered[2],
  variable_name = levels_pizza_ordered[1],
  levels_order = levels_pizza_ordered,
  model_name = "model1"
)


### -------------------- Step 4.33 snacks_pastries --------------------
levels_snacks_pastries_ordered <- {
  vals <- levels_snacks_pastries
  c(vals[2], vals[-2])
}

snacks_pastries_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "snacks_pastries_group_range",
  var_num = "snacks_pastries_group_range_num",
  group_name = "Snacks and pastries",
  ref_level = levels_snacks_pastries_ordered[2],
  levels_order = levels_snacks_pastries_ordered,
  model_name = "model1"
)

### -------------------- Step 4.34 fried_fast_foods --------------------
levels_fried_fast_foods_ordered <- {
  vals <- levels_fried_fast_foods
  c(vals[2], vals[-2])
}

fried_fast_foods_model1_result <- run_food_analysis_model1(
  data = ukb_data_depression,
  var = "fried_fast_foods_group_range",
  var_num = "fried_fast_foods_group_range_num",
  group_name = "Fried and fast foods",
  ref_level = levels_fried_fast_foods_ordered[2],
  levels_order = levels_fried_fast_foods_ordered,
  model_name = "model1"
)

### -------------------- Step 4.35 results combination --------------------
model1_final_result <- bind_rows(
  berries_model1_result,
  citrus_fruits_model1_result,
  other_fruits_model1_result,
  green_leafy_veg_model1_result,
  red_orange_veg_model1_result,
  other_veg_model1_result,
  potatoes_model1_result,
  refined_grains_model1_result,
  whole_grains_model1_result,
  legumes_model1_result,
  nuts_model1_result,
  eggs_model1_result,
  high_fat_dairy_model1_result,
  low_fat_dairy_model1_result,
  red_meats_model1_result,
  processed_meats_model1_result,
  poultry_model1_result,
  fish_seafood_model1_result,
  organ_meats_model1_result,
  butter_model1_result,
  margarine_model1_result,
  olive_oil_model1_result,
  other_oil_model1_result,
  beer_model1_result,
  coffee_model1_result,
  fruit_juices_model1_result,
  liquor_model1_result,
  sweetened_bev_model1_result,
  tea_model1_result,
  water_model1_result,
  wine_model1_result,
  pizza_model1_result,
  snacks_pastries_model1_result,
  fried_fast_foods_model1_result
)  %>%
  # 仅保留食物的第一行
  mutate(group = if_else(duplicated(group), NA_character_, group))


## -------------------- Step 5 food group Model 2 --------------------
paste(colnames(ukb_data), collapse = " + ")

### -------------------- Step 5.0 模型定义与函数封装 --------------------
# Model 2 adjusted for age_food + sex + ethnicity + education_attainment + twnsend + career + smoking_status + drink_status + bmi_group + physical_activity_proup + sedentary_group + household_income + energy_intake + central_obesity + triglycerides + hdl_cholesterol + t2diabetes_status_final + hypertension_status_final + dyslipidemia_status_final + cvd_status_final + all_cancer_status_final



run_food_analysis_model2 <- function(data, var, var_num, group_name, ref_level,
                                     levels_order, variable_name = NULL,
                                     model_name = "model2") {
  rotation_fit_one_food(
    data = data, var = var, var_num = var_num, group_name = group_name,
    ref_level = ref_level, levels_order = levels_order,
    variable_name = variable_name, model_name = model_name,
    covariate_text = "age_food + sex + ethnicity + education_attainment + twnsend + career + smoking_status + drink_status + bmi_group + physical_activity_proup + sedentary_group + household_income + energy_intake + central_obesity + triglycerides + hdl_cholesterol + t2diabetes_status_final + hypertension_status_final + dyslipidemia_status_final + cvd_status_final + all_cancer_status_final"
  )
}

run_food_analysis_two_level_model2 <- function(data, var, var_num, group_name, ref_level,
                                     levels_order, variable_name = NULL,
                                     model_name = "model2") {
  rotation_fit_one_food(
    data = data, var = var, var_num = var_num, group_name = group_name,
    ref_level = ref_level, levels_order = levels_order,
    variable_name = variable_name, model_name = model_name,
    covariate_text = "age_food + sex + ethnicity + education_attainment + twnsend + career + smoking_status + drink_status + bmi_group + physical_activity_proup + sedentary_group + household_income + energy_intake + central_obesity + triglycerides + hdl_cholesterol + t2diabetes_status_final + hypertension_status_final + dyslipidemia_status_final + cvd_status_final + all_cancer_status_final"
  )
}

### -------------------- Step 5.1 berries --------------------
berries_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "berries_group_range",
  var_num = "berries_group_range_num",
  group_name = "Berries",
  ref_level = levels_berries_ordered[2],
  levels_order = levels_berries_ordered,
  model_name = "model2"
)

berries_model2_result

### -------------------- Step 5.2 citrus_fruits --------------------
citrus_fruits_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "citrus_fruits_group_range",
  var_num = "citrus_fruits_group_range_num",
  group_name = "Citrus fruits",
  ref_level = levels_citrus_fruits_ordered[2],
  levels_order = levels_citrus_fruits_ordered,
  model_name = "model2"
)

### -------------------- Step 5.3 other_fruits --------------------
other_fruits_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "other_fruits_group_range",
  var_num = "other_fruits_group_range_num",
  group_name = "Other fruits",
  ref_level = levels_other_fruits_ordered[2],
  levels_order = levels_other_fruits_ordered,
  model_name = "model2"
)

### -------------------- Step 5.4 green_leafy_veg --------------------
green_leafy_veg_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "green_leafy_veg_group_range",
  var_num = "green_leafy_veg_group_range_num",
  group_name = "Green leafy vegetables",
  ref_level = levels_green_leafy_veg_ordered[2],
  levels_order = levels_green_leafy_veg_ordered,
  model_name = "model2"
)

### -------------------- Step 5.5 red_orange_veg --------------------
red_orange_veg_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "red_orange_veg_group_range",
  var_num = "red_orange_veg_group_range_num",
  group_name = "Red/orange vegetables",
  ref_level = levels_red_orange_veg_ordered[2],
  levels_order = levels_red_orange_veg_ordered,
  model_name = "model2"
)

### -------------------- Step 5.6 other_veg --------------------
other_veg_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "other_veg_group_range",
  var_num = "other_veg_group_range_num",
  group_name = "Other vegetables",
  ref_level = levels_other_veg_ordered[2],
  levels_order = levels_other_veg_ordered,
  model_name = "model2"
)

### -------------------- Step 5.7 potatoes --------------------
potatoes_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "potatoes_group_range",
  var_num = "potatoes_group_range_num",
  group_name = "Potatoes",
  ref_level = levels_potatoes_ordered[2],
  levels_order = levels_potatoes_ordered,
  model_name = "model2"
)

### -------------------- Step 5.8 refined_grains --------------------
refined_grains_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "refined_grains_group_range",
  var_num = "refined_grains_group_range_num",
  group_name = "Refined grains",
  ref_level = levels_refined_grains_ordered[2],
  levels_order = levels_refined_grains_ordered,
  model_name = "model2"
)

### -------------------- Step 5.9 whole_grains --------------------
whole_grains_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "whole_grains_group_range",
  var_num = "whole_grains_group_range_num",
  group_name = "Whole grains",
  ref_level = levels_whole_grains_ordered[2],
  levels_order = levels_whole_grains_ordered,
  model_name = "model2"
)

### -------------------- Step 5.10 legumes --------------------
legumes_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "legumes_group_range",
  var_num = "legumes_group_range_num",
  group_name = "Legumes",
  ref_level = levels_legumes_ordered[2],
  levels_order = levels_legumes_ordered,
  model_name = "model2"
)


### -------------------- Step 5.11 nuts --------------------
nuts_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "nuts_group_range",
  var_num = "nuts_group_range_num",
  group_name = "Nuts",
  ref_level = levels_nuts_ordered[2],
  levels_order = levels_nuts_ordered,
  model_name = "model2"
)

### -------------------- Step 5.12 eggs --------------------
eggs_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "eggs_group_range",
  var_num = "eggs_group_range_num",
  group_name = "Eggs",
  ref_level = levels_eggs_ordered[2],
  levels_order = levels_eggs_ordered,
  model_name = "model2"
)

### -------------------- Step 5.13 high_fat_dairy --------------------
high_fat_dairy_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "high_fat_dairy_group_range",
  var_num = "high_fat_dairy_group_range_num",
  group_name = "High fat dairy",
  ref_level = levels_high_fat_dairy_ordered[2],
  levels_order = levels_high_fat_dairy_ordered,
  model_name = "model2"
)

### -------------------- Step 5.14 low_fat_dairy --------------------
low_fat_dairy_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "low_fat_dairy_group_range",
  var_num = "low_fat_dairy_group_range_num",
  group_name = "Low fat dairy",
  ref_level = levels_low_fat_dairy_ordered[2],
  levels_order = levels_low_fat_dairy_ordered,
  model_name = "model2"
)

### -------------------- Step 5.15 red_meats --------------------
red_meats_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "red_meats_group_range",
  var_num = "red_meats_group_range_num",
  group_name = "Red meats",
  ref_level = levels_red_meats_ordered[2],
  levels_order = levels_red_meats_ordered,
  model_name = "model2"
)

### -------------------- Step 5.16 processed_meats --------------------
processed_meats_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "processed_meats_group_range",
  var_num = "processed_meats_group_range_num",
  group_name = "Processed meats",
  ref_level = levels_processed_meats_ordered[2],
  levels_order = levels_processed_meats_ordered,
  model_name = "model2"
)

### -------------------- Step 5.17 poultry --------------------
poultry_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "poultry_group_range",
  var_num = "poultry_group_range_num",
  group_name = "Poultry",
  ref_level = levels_poultry_ordered[2],
  levels_order = levels_poultry_ordered,
  model_name = "model2"
)

### -------------------- Step 5.18 fish_seafood --------------------
fish_seafood_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "fish_seafood_group_range",
  var_num = "fish_seafood_group_range_num",
  group_name = "Fish and seafood",
  ref_level = levels_fish_seafood_ordered[2],
  levels_order = levels_fish_seafood_ordered,
  model_name = "model2"
)

### -------------------- Step 5.19 organ_meats --------------------
organ_meats_model2_result <- run_food_analysis_two_level_model2(
  data = ukb_data_depression,
  var = "organ_meats_group_range",
  var_num = "organ_meats_group_range_num",
  group_name = "Organ meats",
  ref_level = levels_organ_meats_ordered[2],
  variable_name = levels_organ_meats_ordered[1],
  levels_order = levels_organ_meats_ordered,
  model_name = "model2"
)

### -------------------- Step 5.20 butter --------------------
butter_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "butter_group_range",
  var_num = "butter_group_range_num",
  group_name = "Butter",
  ref_level = levels_butter_ordered[2],
  levels_order = levels_butter_ordered,
  model_name = "model2"
)


### -------------------- Step 5.21 margarine --------------------
margarine_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "margarine_group_range",
  var_num = "margarine_group_range_num",
  group_name = "Margarine",
  ref_level = levels_margarine_ordered[2],
  levels_order = levels_margarine_ordered,
  model_name = "model2"
)

### -------------------- Step 5.22 olive_oil --------------------
olive_oil_model2_result <- run_food_analysis_two_level_model2(
  data = ukb_data_depression,
  var = "olive_oil_group_range",
  var_num = "olive_oil_group_range_num",
  group_name = "Olive oil",
  ref_level = levels_olive_oil_ordered[2],
  variable_name = levels_olive_oil_ordered[1],
  levels_order = levels_olive_oil_ordered,
  model_name = "model2"
)

### -------------------- Step 5.23 other_oil --------------------
other_oil_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "other_oil_group_range",
  var_num = "other_oil_group_range_num",
  group_name = "Other oil",
  ref_level = levels_other_oil_ordered[2],
  levels_order = levels_other_oil_ordered,
  model_name = "model2"
)

### -------------------- Step 5.24 beer --------------------
beer_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "beer_group_range",
  var_num = "beer_group_range_num",
  group_name = "Beer",
  ref_level = levels_beer_ordered[2],
  levels_order = levels_beer_ordered,
  model_name = "model2"
)

### -------------------- Step 5.25 coffee --------------------
coffee_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "coffee_group_range",
  var_num = "coffee_group_range_num",
  group_name = "Coffee",
  ref_level = levels_coffee_ordered[2],
  levels_order = levels_coffee_ordered,
  model_name = "model2"
)

### -------------------- Step 5.26 fruit_juices --------------------
fruit_juices_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "fruit_juices_group_range",
  var_num = "fruit_juices_group_range_num",
  group_name = "Fruit juices",
  ref_level = levels_fruit_juices_ordered[2],
  levels_order = levels_fruit_juices_ordered,
  model_name = "model2"
)

### -------------------- Step 5.27 liquor --------------------
liquor_model2_result <- run_food_analysis_two_level_model2(
  data = ukb_data_depression,
  var = "liquor_group_range",
  var_num = "liquor_group_range_num",
  group_name = "Liquor",
  ref_level = levels_liquor_ordered[2],
  variable_name = levels_liquor_ordered[1],
  levels_order = levels_liquor_ordered,
  model_name = "model2"
)

### -------------------- Step 5.28 sweetened_bev --------------------
sweetened_bev_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "sweetened_bev_group_range",
  var_num = "sweetened_bev_group_range_num",
  group_name = "Sweetened beverages",
  ref_level = levels_sweetened_bev_ordered[2],
  levels_order = levels_sweetened_bev_ordered,
  model_name = "model2"
)

### -------------------- Step 5.29 tea --------------------
tea_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "tea_group_range",
  var_num = "tea_group_range_num",
  group_name = "Tea",
  ref_level = levels_tea_ordered[2],
  levels_order = levels_tea_ordered,
  model_name = "model2"
)

### -------------------- Step 5.30 water --------------------
water_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "water_group_range",
  var_num = "water_group_range_num",
  group_name = "Water",
  ref_level = levels_water_ordered[2],
  levels_order = levels_water_ordered,
  model_name = "model2"
)

### -------------------- Step 5.31 wine --------------------
wine_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "wine_group_range",
  var_num = "wine_group_range_num",
  group_name = "Wine",
  ref_level = levels_wine_ordered[2],
  levels_order = levels_wine_ordered,
  model_name = "model2"
)

### -------------------- Step 5.32 pizza --------------------
pizza_model2_result <- run_food_analysis_two_level_model2(
  data = ukb_data_depression,
  var = "pizza_group_range",
  var_num = "pizza_group_range_num",
  group_name = "Pizza",
  ref_level = levels_pizza_ordered[2],
  variable_name = levels_pizza_ordered[1],
  levels_order = levels_pizza_ordered,
  model_name = "model2"
)


### -------------------- Step 5.33 snacks_pastries --------------------
snacks_pastries_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "snacks_pastries_group_range",
  var_num = "snacks_pastries_group_range_num",
  group_name = "Snacks and pastries",
  ref_level = levels_snacks_pastries_ordered[2],
  levels_order = levels_snacks_pastries_ordered,
  model_name = "model2"
)

### -------------------- Step 5.34 fried_fast_foods --------------------
fried_fast_foods_model2_result <- run_food_analysis_model2(
  data = ukb_data_depression,
  var = "fried_fast_foods_group_range",
  var_num = "fried_fast_foods_group_range_num",
  group_name = "Fried and fast foods",
  ref_level = levels_fried_fast_foods_ordered[2],
  levels_order = levels_fried_fast_foods_ordered,
  model_name = "model2"
)

### -------------------- Step 5.35 results combination --------------------
model2_final_result <- bind_rows(
  berries_model2_result,
  citrus_fruits_model2_result,
  other_fruits_model2_result,
  green_leafy_veg_model2_result,
  red_orange_veg_model2_result,
  other_veg_model2_result,
  potatoes_model2_result,
  refined_grains_model2_result,
  whole_grains_model2_result,
  legumes_model2_result,
  nuts_model2_result,
  eggs_model2_result,
  high_fat_dairy_model2_result,
  low_fat_dairy_model2_result,
  red_meats_model2_result,
  processed_meats_model2_result,
  poultry_model2_result,
  fish_seafood_model2_result,
  organ_meats_model2_result,
  butter_model2_result,
  margarine_model2_result,
  olive_oil_model2_result,
  other_oil_model2_result,
  beer_model2_result,
  coffee_model2_result,
  fruit_juices_model2_result,
  liquor_model2_result,
  sweetened_bev_model2_result,
  tea_model2_result,
  water_model2_result,
  wine_model2_result,
  pizza_model2_result,
  snacks_pastries_model2_result,
  fried_fast_foods_model2_result
) %>%
  # 仅保留食物的第一行
  mutate(group = if_else(duplicated(group), NA_character_, group))


### -------------------- Step 5.36 All-reference audit and QC --------------------
# The original 34 food groups are fit with every available category as reference
# in Model 1 and Model 2. Final figure results use the selected minimum-HR ref.
rotation_selected_reference <- dplyr::bind_rows(ROTATION_RECORDS$selected)
rotation_all_reference_comparisons <- dplyr::bind_rows(
  ROTATION_RECORDS$all_reference_comparisons
)
rotation_fit_diagnostics <- dplyr::bind_rows(ROTATION_RECORDS$diagnostics)
rotation_fit_warnings <- if (length(ROTATION_RECORDS$warnings) > 0L) {
  dplyr::bind_rows(ROTATION_RECORDS$warnings)
} else {
  data.frame(fit_tag = character(0), warning = character(0),
             stringsAsFactors = FALSE)
}
if (nrow(ukb_data_depression_train) != 76578L ||
    nrow(ukb_data_depression_test) != 32819L ||
    sum(ukb_data_depression_train$depression_status == 1, na.rm = TRUE) != 2122L ||
    sum(ukb_data_depression_test$depression_status == 1, na.rm = TRUE) != 875L) {
  stop("STOP: Frozen 70/30 participant counts or events changed unexpectedly.")
}
if (nrow(rotation_selected_reference) != 68L ||
    dplyr::n_distinct(rotation_selected_reference$food_variable) != 34L ||
    !setequal(rotation_selected_reference$model, c("model1", "model2"))) {
  stop("STOP: Expected 34 food groups with both Model 1 and Model 2 rotation results.")
}
if (any(!is.finite(rotation_fit_diagnostics$max_contrast_delta)) ||
    any(rotation_fit_diagnostics$max_contrast_delta > 1e-4)) {
  stop("STOP: Cross-reference Cox contrast invariance QC failed.")
}
if (!all(rotation_selected_reference$model_n > 0) ||
    anyNA(rotation_selected_reference$lowest_risk_reference_category)) {
  stop("STOP: Invalid category-reference selection QC.")
}
rotation_counts <- rotation_selected_reference %>%
  dplyr::count(model, name = "food_groups")
rotation_qc <- data.frame(
  selected_rows = nrow(rotation_selected_reference),
  food_groups = dplyr::n_distinct(rotation_selected_reference$food_variable),
  categorical_fits = nrow(rotation_fit_diagnostics),
  pairwise_comparison_rows = nrow(rotation_all_reference_comparisons),
  reference_invariance_max_delta = max(rotation_fit_diagnostics$max_contrast_delta),
  warning_rows = nrow(rotation_fit_warnings),
  training_n = nrow(ukb_data_depression_train),
  training_events = sum(ukb_data_depression_train$depression_status == 1, na.rm = TRUE),
  test_n_used_for_fitting = 0L,
  manifest_md5 = manifest_md5,
  overall_status = if (nrow(rotation_fit_warnings) == 0L) "PASS" else
    "PASS_WITH_WARNINGS",
  stringsAsFactors = FALSE
)
utils::write.csv(rotation_selected_reference,
  file.path(ROTATE_RESULTS_DIR, "food_group_reference_rotation_selected_training.csv"),
  row.names = FALSE)
utils::write.csv(rotation_all_reference_comparisons,
  file.path(ROTATE_RESULTS_DIR, "food_group_reference_rotation_all_comparisons_training.csv"),
  row.names = FALSE)
utils::write.csv(rotation_fit_diagnostics,
  file.path(ROTATE_RESULTS_DIR, "food_group_reference_rotation_fit_diagnostics_training.csv"),
  row.names = FALSE)
utils::write.csv(rotation_fit_warnings,
  file.path(ROTATE_RESULTS_DIR, "food_group_reference_rotation_fit_warnings_training.csv"),
  row.names = FALSE)
utils::write.csv(rotation_qc,
  file.path(ROTATE_RESULTS_DIR, "food_group_reference_rotation_QC_training.csv"),
  row.names = FALSE)
cat("\n===== Reference rotation audit =====\n")
print(rotation_counts)
print(rotation_qc)

## -------------------- Step 6 数据导出及森林图绘制 --------------------
### -------------------- Step 6.1 数据导出 --------------------
### -------------------- Step 6.1.1 Model 1 --------------------
write.csv(model1_final_result, file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model1_results.csv", row.names = FALSE)

### -------------------- Step 6.1.2 Model 1 sig food group --------------------
model1_final_result_sig <- model1_final_result %>%
  # Only the first row prints each food name; restore it for group-wise QC.
  tidyr::fill(group, .direction = "down") %>%
  mutate(p_value_num = rotation_parse_p(p_value)) %>%
  group_by(group) %>%
  # Selection uses unrounded 95% CIs, never displayed/rounded P labels.
  filter(any((ci_low > 1 | ci_high < 1) & !is.na(hazard_ratio),
             na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(group = if_else(duplicated(group), NA_character_, group))

length(unique(model1_final_result_sig$group))

# 保存结果
write.csv(model1_final_result_sig, file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model1_results_sig.csv", row.names = FALSE)

### -------------------- Step 6.1.3 Model 2 --------------------
write.csv(model2_final_result, file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model2_results.csv", row.names = FALSE)

### -------------------- Step 6.1.4 Model 2 sig food group --------------------
model2_final_result_sig <- model2_final_result %>%
  # Only the first row prints each food name; restore it for group-wise QC.
  tidyr::fill(group, .direction = "down") %>%
  mutate(p_value_num = rotation_parse_p(p_value)) %>%
  group_by(group) %>%
  # Selection uses unrounded 95% CIs, never displayed/rounded P labels.
  filter(any((ci_low > 1 | ci_high < 1) & !is.na(hazard_ratio),
             na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(group = if_else(duplicated(group), NA_character_, group))

length(unique(model2_final_result_sig$group))

# 保存结果
write.csv(model2_final_result_sig, file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model2_results_sig.csv", row.names = FALSE)


### -------------------- Step 6.2 forestplot --------------------
mycolors <- c4a("friendly13")

### -------------------- Step 6.2.1 Model 1 --------------------
forest_model1_final_result <- model1_final_result

# 根据forestplot包生成颜色的要求，先生成备用数据，并删除缺失hazard_ratio
model_col_number <- forest_model1_final_result %>%
  filter(!is.na(hazard_ratio))

# 把 p_value 转成数值
model_col_number$p_value_num <- rotation_parse_p(model_col_number$p_value)

# 找出显著行的行号（p<0.05 且 hazard_ratio 非 NA）
sig_rows <- which(
  !is.na(model_col_number$hazard_ratio) &
    (model_col_number$ci_low > 1 | model_col_number$ci_high < 1)
)

sig_rows

# 定义 fn
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i %in% sig_rows){
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1])
    } else {
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3])
    }
  }
})


# 查看数据格式
colnames(forest_model1_final_result)
head(forest_model1_final_result)


# 数据标签准备
tabletext_HR <- cbind(c("Food Groups", forest_model1_final_result$group),
                      c("Consumption Level", forest_model1_final_result$variable),
                      c("Events/Total", forest_model1_final_result$Events_Total),
                      c("HR (95% CI)", forest_model1_final_result$hr_ci),
                      c("P value", forest_model1_final_result$p_value))

# 绘图
model1_forestplot <- forestplot(
  labeltext = tabletext_HR,
  graph.pos = 5,
  fn.ci_norm = fn, # g更改颜色用fn方案
  boxsize = 0.2,
  mean = c(NA, forest_model1_final_result$hazard_ratio),
  lower = c(NA, forest_model1_final_result$ci_low),
  upper = c(NA, forest_model1_final_result$ci_high),
  xlog = T,
  col = fpColors(zero = "black"),
  
  #-------------森林图细节
  zero = 1, #设置无效线参照HR值1
  graphwidth = unit(0.10,"npc"),
  lty.zero = 3, # 无效线线条类型
  lwd.zero = 1.5, # 无效线宽度
  
  #----------坐标轴设置
  xlab="Hazard ratio (95% CI)",#X轴标签
  lwd.xaxis = 1.5, #X轴线宽
  # clip=c(0.4,8),#X轴范围
  xticks = c(exp(-log(1.5)), 1, 1.5)
  # ci.vertices = TRUE,
  # ci.vertices.height = 0.05
)

# model1_forestplot

# 导出图片
pdf(file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model1_forestplot.pdf", width = 15, height = 35)
print(model1_forestplot)
dev.off()


### -------------------- Step 6.2.2 Model 1 sig food group --------------------
forest_model1_final_result <- model1_final_result_sig

# 根据forestplot包生成颜色的要求，先生成备用数据，并删除缺失hazard_ratio
model_col_number <- forest_model1_final_result %>%
  filter(!is.na(hazard_ratio))

# 把 p_value 转成数值
model_col_number$p_value_num <- rotation_parse_p(model_col_number$p_value)

# 找出显著行的行号（p<0.05 且 hazard_ratio 非 NA）
sig_rows <- which(
  !is.na(model_col_number$hazard_ratio) &
    (model_col_number$ci_low > 1 | model_col_number$ci_high < 1)
)

sig_rows

# 定义 fn
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i %in% sig_rows){
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1])
    } else {
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3])
    }
  }
})


# 查看数据格式
colnames(forest_model1_final_result)
head(forest_model1_final_result)


# 数据标签准备
tabletext_HR <- cbind(c("Food Groups", forest_model1_final_result$group),
                      c("Consumption Level", forest_model1_final_result$variable),
                      c("Events/Total", forest_model1_final_result$Events_Total),
                      c("HR (95% CI)", forest_model1_final_result$hr_ci),
                      c("P value", forest_model1_final_result$p_value))

# 绘图
model1_forestplot <- forestplot(
  labeltext = tabletext_HR,
  graph.pos = 5,
  fn.ci_norm = fn, # g更改颜色用fn方案
  boxsize = 0.2,
  mean = c(NA, forest_model1_final_result$hazard_ratio),
  lower = c(NA, forest_model1_final_result$ci_low),
  upper = c(NA, forest_model1_final_result$ci_high),
  xlog = T,
  col = fpColors(zero = "black"),
  
  #-------------森林图细节
  zero = 1, #设置无效线参照HR值1
  graphwidth = unit(0.10,"npc"),
  lty.zero = 3, # 无效线线条类型
  lwd.zero = 1.5, # 无效线宽度
  
  #----------坐标轴设置
  xlab="Hazard ratio (95% CI)",#X轴标签
  lwd.xaxis = 1.5, #X轴线宽
  # clip=c(0.4,8),#X轴范围
  xticks = c(exp(-log(1.5)), 1, 1.5)
  # ci.vertices = TRUE,
  # ci.vertices.height = 0.05
)

# model1_forestplot

# 导出图片
pdf(file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model1_forestplot_sig.pdf", width = 15, height = 20)
print(model1_forestplot)
dev.off()

### -------------------- Step 6.2.3 Model 2 --------------------
forest_model2_final_result <- model2_final_result

# 根据forestplot包生成颜色的要求，先生成备用数据，并删除缺失hazard_ratio
model_col_number <- forest_model2_final_result %>%
  filter(!is.na(hazard_ratio))

# 把 p_value 转成数值
model_col_number$p_value_num <- rotation_parse_p(model_col_number$p_value)

# 找出显著行的行号（p<0.05 且 hazard_ratio 非 NA）
sig_rows <- which(
  !is.na(model_col_number$hazard_ratio) &
    (model_col_number$ci_low > 1 | model_col_number$ci_high < 1)
)

sig_rows

# 定义 fn
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i %in% sig_rows){
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1])
    } else {
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3])
    }
  }
})


# 查看数据格式
colnames(forest_model2_final_result)
head(forest_model2_final_result)


# 数据标签准备
tabletext_HR <- cbind(c("Food Groups", forest_model2_final_result$group),
                      c("Consumption Level", forest_model2_final_result$variable),
                      c("Events/Total", forest_model2_final_result$Events_Total),
                      c("HR (95% CI)", forest_model2_final_result$hr_ci),
                      c("P value", forest_model2_final_result$p_value))

# 绘图
model2_forestplot <- forestplot(
  labeltext = tabletext_HR,
  graph.pos = 5,
  fn.ci_norm = fn, # g更改颜色用fn方案
  boxsize = 0.2,
  mean = c(NA, forest_model2_final_result$hazard_ratio),
  lower = c(NA, forest_model2_final_result$ci_low),
  upper = c(NA, forest_model2_final_result$ci_high),
  xlog = T,
  col = fpColors(zero = "black"),
  
  #-------------森林图细节
  zero = 1, #设置无效线参照HR值1
  graphwidth = unit(0.10,"npc"),
  lty.zero = 3, # 无效线线条类型
  lwd.zero = 1.5, # 无效线宽度
  
  #----------坐标轴设置
  xlab="Hazard ratio (95% CI)",#X轴标签
  lwd.xaxis = 1.5, #X轴线宽
  # clip=c(0.4,8),#X轴范围
  xticks = c(exp(-log(1.5)), 1, 1.5)
  # ci.vertices = TRUE,
  # ci.vertices.height = 0.05
)

# model2_forestplot

# 导出图片
pdf(file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model2_forestplot.pdf", width = 15, height = 16)
print(model2_forestplot)
dev.off()


### -------------------- Step 6.2.4 Model 2 sig food group --------------------
forest_model2_final_result <- model2_final_result_sig

# 根据forestplot包生成颜色的要求，先生成备用数据，并删除缺失hazard_ratio
model_col_number <- forest_model2_final_result %>%
  filter(!is.na(hazard_ratio))

# 把 p_value 转成数值
model_col_number$p_value_num <- rotation_parse_p(model_col_number$p_value)

# 找出显著行的行号（p<0.05 且 hazard_ratio 非 NA）
sig_rows <- which(
  !is.na(model_col_number$hazard_ratio) &
    (model_col_number$ci_low > 1 | model_col_number$ci_high < 1)
)

sig_rows

# 定义 fn
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i %in% sig_rows){
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1])
    } else {
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3])
    }
  }
})


# 查看数据格式
colnames(forest_model2_final_result)
head(forest_model2_final_result)


# 数据标签准备
tabletext_HR <- cbind(c("Food Groups", forest_model2_final_result$group),
                      c("Consumption Level", forest_model2_final_result$variable),
                      c("Events/Total", forest_model2_final_result$Events_Total),
                      c("HR (95% CI)", forest_model2_final_result$hr_ci),
                      c("P value", forest_model2_final_result$p_value))

# 绘图
model2_forestplot <- forestplot(
  labeltext = tabletext_HR,
  graph.pos = 5,
  fn.ci_norm = fn, # g更改颜色用fn方案
  boxsize = 0.2,
  mean = c(NA, forest_model2_final_result$hazard_ratio),
  lower = c(NA, forest_model2_final_result$ci_low),
  upper = c(NA, forest_model2_final_result$ci_high),
  xlog = T,
  col = fpColors(zero = "black"),
  
  #-------------森林图细节
  zero = 1, #设置无效线参照HR值1
  graphwidth = unit(0.10,"npc"),
  lty.zero = 3, # 无效线线条类型
  lwd.zero = 1.5, # 无效线宽度
  
  #----------坐标轴设置
  xlab="Hazard ratio (95% CI)",#X轴标签
  lwd.xaxis = 1.5, #X轴线宽
  # clip=c(0.4,8),#X轴范围
  xticks = c(exp(-log(1.5)), 1, 1.5)
  # ci.vertices = TRUE,
  # ci.vertices.height = 0.05
)

# model2_forestplot

# 导出图片
pdf(file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model2_forestplot_sig.pdf", width = 15, height = 16)
print(model2_forestplot)
dev.off()


### -------------------- Step 6.2.5 Model 2 上半部分 --------------------
model2_final_result_84 <- model2_final_result %>%
  slice(1:84)

forest_model2_final_result <- model2_final_result_84

# 根据forestplot包生成颜色的要求，先生成备用数据，并删除缺失hazard_ratio
model_col_number <- forest_model2_final_result %>%
  filter(!is.na(hazard_ratio))

# 把 p_value 转成数值
model_col_number$p_value_num <- rotation_parse_p(model_col_number$p_value)

# 找出显著行的行号（p<0.05 且 hazard_ratio 非 NA）
sig_rows <- which(
  !is.na(model_col_number$hazard_ratio) &
    (model_col_number$ci_low > 1 | model_col_number$ci_high < 1)
)

sig_rows

# 定义 fn
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i %in% sig_rows){
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1])
    } else {
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3])
    }
  }
})


# 查看数据格式
colnames(forest_model2_final_result)
head(forest_model2_final_result)


# 数据标签准备
tabletext_HR <- cbind(c("Food Groups", forest_model2_final_result$group),
                      c("Consumption Level", forest_model2_final_result$variable),
                      c("Events/Total", forest_model2_final_result$Events_Total),
                      c("HR (95% CI)", forest_model2_final_result$hr_ci),
                      c("P value", forest_model2_final_result$p_value))

# 绘图
model2_forestplot <- forestplot(
  labeltext = tabletext_HR,
  graph.pos = 5,
  fn.ci_norm = fn, # g更改颜色用fn方案
  boxsize = 0.2,
  mean = c(NA, forest_model2_final_result$hazard_ratio),
  lower = c(NA, forest_model2_final_result$ci_low),
  upper = c(NA, forest_model2_final_result$ci_high),
  xlog = T,
  col = fpColors(zero = "black"),
  
  #-------------森林图细节
  zero = 1, #设置无效线参照HR值1
  graphwidth = unit(0.10,"npc"),
  lty.zero = 3, # 无效线线条类型
  lwd.zero = 1.5, # 无效线宽度
  
  #----------坐标轴设置
  xlab="Hazard ratio (95% CI)",#X轴标签
  lwd.xaxis = 1.5, #X轴线宽
  # clip=c(0.4,8),#X轴范围
  xticks = c(exp(-log(1.5)), 1, 1.5),
  
  # ci.vertices = TRUE,
  # ci.vertices.height = 0.05
  # -------- 字体设置
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.1),   # 左侧标签 fontface = "bold",
    ticks = gpar(cex = 0.8),   # 坐标轴刻度
    xlab  = gpar(cex = 0.9),   # X轴标题 
    # title = gpar(fontface = "bold", cex = 1.4)  # 图标题（如果有）
  )
)
# model2_forestplot

# 导出图片
pdf(file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model2_forestplot_84.pdf", width = 15, height = 20)
print(model2_forestplot)
dev.off()


### -------------------- Step 6.2.6 Model 2 下半部分 --------------------
model2_final_result_170 <- model2_final_result %>%
  slice(85:170)

forest_model2_final_result <- model2_final_result_170

# 根据forestplot包生成颜色的要求，先生成备用数据，并删除缺失hazard_ratio
model_col_number <- forest_model2_final_result %>%
  filter(!is.na(hazard_ratio))

# 把 p_value 转成数值
model_col_number$p_value_num <- rotation_parse_p(model_col_number$p_value)

# 找出显著行的行号（p<0.05 且 hazard_ratio 非 NA）
sig_rows <- which(
  !is.na(model_col_number$hazard_ratio) &
    (model_col_number$ci_low > 1 | model_col_number$ci_high < 1)
)

sig_rows

# 定义 fn
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i %in% sig_rows){
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1])
    } else {
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3])
    }
  }
})


# 查看数据格式
colnames(forest_model2_final_result)
head(forest_model2_final_result)


# 数据标签准备
tabletext_HR <- cbind(c("Food Groups", forest_model2_final_result$group),
                      c("Consumption Level", forest_model2_final_result$variable),
                      c("Events/Total", forest_model2_final_result$Events_Total),
                      c("HR (95% CI)", forest_model2_final_result$hr_ci),
                      c("P value", forest_model2_final_result$p_value))

# 绘图
model2_forestplot <- forestplot(
  labeltext = tabletext_HR,
  graph.pos = 5,
  fn.ci_norm = fn, # g更改颜色用fn方案
  boxsize = 0.2,
  mean = c(NA, forest_model2_final_result$hazard_ratio),
  lower = c(NA, forest_model2_final_result$ci_low),
  upper = c(NA, forest_model2_final_result$ci_high),
  xlog = T,
  col = fpColors(zero = "black"),
  
  #-------------森林图细节
  zero = 1, #设置无效线参照HR值1
  graphwidth = unit(0.10,"npc"),
  lty.zero = 3, # 无效线线条类型
  lwd.zero = 1.5, # 无效线宽度
  
  #----------坐标轴设置
  xlab="Hazard ratio (95% CI)",#X轴标签
  lwd.xaxis = 1.5, #X轴线宽
  # clip=c(0.4,8),#X轴范围
  xticks = c(exp(-log(1.5)), 1, 1.5),
  
  # ci.vertices = TRUE,
  # ci.vertices.height = 0.05
  # -------- 字体设置
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.1),   # 左侧标签 fontface = "bold",
    ticks = gpar(cex = 0.8),   # 坐标轴刻度
    xlab  = gpar(cex = 0.9),   # X轴标题 
    # title = gpar(fontface = "bold", cex = 1.4)  # 图标题（如果有）
  )
)
# model2_forestplot

# 导出图片
pdf(file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model2_forestplot_170.pdf", width = 15, height = 20)
print(model2_forestplot)
dev.off()

### -------------------- Step 6.2.7 Model 1 上半部分 --------------------
model1_final_result_84 <- model1_final_result %>%
  slice(1:84)

forest_model1_final_result <- model1_final_result_84

# 根据forestplot包生成颜色的要求，先生成备用数据，并删除缺失hazard_ratio
model_col_number <- forest_model1_final_result %>%
  filter(!is.na(hazard_ratio))

# 把 p_value 转成数值
model_col_number$p_value_num <- rotation_parse_p(model_col_number$p_value)

# 找出显著行的行号（p<0.05 且 hazard_ratio 非 NA）
sig_rows <- which(
  !is.na(model_col_number$hazard_ratio) &
    (model_col_number$ci_low > 1 | model_col_number$ci_high < 1)
)

sig_rows

# 定义 fn
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i %in% sig_rows){
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1])
    } else {
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3])
    }
  }
})


# 查看数据格式
colnames(forest_model1_final_result)
head(forest_model1_final_result)


# 数据标签准备
tabletext_HR <- cbind(c("Food Groups", forest_model1_final_result$group),
                      c("Consumption Level", forest_model1_final_result$variable),
                      c("Events/Total", forest_model1_final_result$Events_Total),
                      c("HR (95% CI)", forest_model1_final_result$hr_ci),
                      c("P value", forest_model1_final_result$p_value))

# 绘图
model1_forestplot <- forestplot(
  labeltext = tabletext_HR,
  graph.pos = 5,
  fn.ci_norm = fn, # g更改颜色用fn方案
  boxsize = 0.2,
  mean = c(NA, forest_model1_final_result$hazard_ratio),
  lower = c(NA, forest_model1_final_result$ci_low),
  upper = c(NA, forest_model1_final_result$ci_high),
  xlog = T,
  col = fpColors(zero = "black"),
  
  #-------------森林图细节
  zero = 1, #设置无效线参照HR值1
  graphwidth = unit(0.10,"npc"),
  lty.zero = 3, # 无效线线条类型
  lwd.zero = 1.5, # 无效线宽度
  
  #----------坐标轴设置
  xlab="Hazard ratio (95% CI)",#X轴标签
  lwd.xaxis = 1.5, #X轴线宽
  # clip=c(0.4,8),#X轴范围
  xticks = c(exp(-log(1.5)), 1, 1.5),
  
  # ci.vertices = TRUE,
  # ci.vertices.height = 0.05
  # -------- 字体设置
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.1),   # 左侧标签 fontface = "bold",
    ticks = gpar(cex = 0.8),   # 坐标轴刻度
    xlab  = gpar(cex = 0.9),   # X轴标题 
    # title = gpar(fontface = "bold", cex = 1.4)  # 图标题（如果有）
  )
)
# model1_forestplot

# 导出图片
pdf(file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model1_forestplot_84.pdf", width = 15, height = 20)
print(model1_forestplot)
dev.off()


### -------------------- Step 6.2.8 Model 1 下半部分 --------------------


model1_final_result_170 <- model1_final_result %>%
  slice(85:170)

forest_model1_final_result <- model1_final_result_170

# 根据forestplot包生成颜色的要求，先生成备用数据，并删除缺失hazard_ratio
model_col_number <- forest_model1_final_result %>%
  filter(!is.na(hazard_ratio))

# 把 p_value 转成数值
model_col_number$p_value_num <- rotation_parse_p(model_col_number$p_value)

# 找出显著行的行号（p<0.05 且 hazard_ratio 非 NA）
sig_rows <- which(
  !is.na(model_col_number$hazard_ratio) &
    (model_col_number$ci_low > 1 | model_col_number$ci_high < 1)
)

sig_rows

# 定义 fn
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i %in% sig_rows){
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1])
    } else {
      fpDrawCircleCI(..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3])
    }
  }
})


# 查看数据格式
colnames(forest_model1_final_result)
head(forest_model1_final_result)


# 数据标签准备
tabletext_HR <- cbind(c("Food Groups", forest_model1_final_result$group),
                      c("Consumption Level", forest_model1_final_result$variable),
                      c("Events/Total", forest_model1_final_result$Events_Total),
                      c("HR (95% CI)", forest_model1_final_result$hr_ci),
                      c("P value", forest_model1_final_result$p_value))

# 绘图
model1_forestplot <- forestplot(
  labeltext = tabletext_HR,
  graph.pos = 5,
  fn.ci_norm = fn, # g更改颜色用fn方案
  boxsize = 0.2,
  mean = c(NA, forest_model1_final_result$hazard_ratio),
  lower = c(NA, forest_model1_final_result$ci_low),
  upper = c(NA, forest_model1_final_result$ci_high),
  xlog = T,
  col = fpColors(zero = "black"),
  
  #-------------森林图细节
  zero = 1, #设置无效线参照HR值1
  graphwidth = unit(0.10,"npc"),
  lty.zero = 3, # 无效线线条类型
  lwd.zero = 1.5, # 无效线宽度
  
  #----------坐标轴设置
  xlab="Hazard ratio (95% CI)",#X轴标签
  lwd.xaxis = 1.5, #X轴线宽
  # clip=c(0.4,8),#X轴范围
  xticks = c(exp(-log(1.5)), 1, 1.5),
  
  # ci.vertices = TRUE,
  # ci.vertices.height = 0.05
  # -------- 字体设置
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.1),   # 左侧标签 fontface = "bold",
    ticks = gpar(cex = 0.8),   # 坐标轴刻度
    xlab  = gpar(cex = 0.9),   # X轴标题 
    # title = gpar(fontface = "bold", cex = 1.4)  # 图标题（如果有）
  )
)
# model1_forestplot

# 导出图片
pdf(file = "./results/Step_03/reference_rotation_lowest_risk/food_group_model1_forestplot_170.pdf", width = 15, height = 20)
print(model1_forestplot)
dev.off()


### -------------------- Step 6.3 Verify all original-style figures and CSVs ----
rotation_expected_pdf <- c(
  "food_group_model1_forestplot.pdf",
  "food_group_model1_forestplot_sig.pdf",
  "food_group_model1_forestplot_84.pdf",
  "food_group_model1_forestplot_170.pdf",
  "food_group_model2_forestplot.pdf",
  "food_group_model2_forestplot_sig.pdf",
  "food_group_model2_forestplot_84.pdf",
  "food_group_model2_forestplot_170.pdf"
)
rotation_expected_csv <- c(
  "food_group_model1_results.csv",
  "food_group_model1_results_sig.csv",
  "food_group_model2_results.csv",
  "food_group_model2_results_sig.csv"
)
rotation_visual_qc <- data.frame(
  type = c(rep("PDF", length(rotation_expected_pdf)),
           rep("CSV", length(rotation_expected_csv))),
  file = c(rotation_expected_pdf, rotation_expected_csv),
  stringsAsFactors = FALSE
)
rotation_visual_qc$path <- file.path(ROTATE_RESULTS_DIR, rotation_visual_qc$file)
rotation_visual_qc$exists <- file.exists(rotation_visual_qc$path)
rotation_visual_qc$bytes <- as.numeric(file.info(rotation_visual_qc$path)$size)
rotation_visual_qc$status <- ifelse(
  rotation_visual_qc$exists & is.finite(rotation_visual_qc$bytes) &
    rotation_visual_qc$bytes > 0,
  "PASS", "FAIL"
)
utils::write.csv(rotation_visual_qc,
  file.path(ROTATE_RESULTS_DIR, "food_group_reference_rotation_figure_and_csv_QC.csv"),
  row.names = FALSE)
print(rotation_visual_qc)
if (any(rotation_visual_qc$status != "PASS")) {
  stop("STOP: One or more of the original-style forest plots/CSVs were not generated.")
}
cat("\n============================================================\n")
cat("REFERENCE ROTATION MAIN COMPLETED — 34 food groups, Model 1/2\n")
cat("All original-style forestplot PDFs and result CSVs have been verified.\n")
cat("Output directory: ", normalizePath(ROTATE_RESULTS_DIR, mustWork = FALSE), "\n", sep = "")
cat("Training-only analysis; split manifest unchanged.\n")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 03/11: Step_03_1_0_MODEL2_REFERENCE
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 04/11: Step_03_2_7_FINAL15_AUDIT
# Source: Step_03_2_7_final15_original_table_item_level_portion.R
# Provenance: local current source
# SHA256 of included block: 2c6d0fbdc966dfa350846aed02c076a89e0135f9a36fdd9972a880f8fc9fb17f
########################################################################################

# =============================================================================
# Step_03_2_7 final15 original-table item-level and portion audit
# V1.0
#
# PURPOSE
#   1. Freeze the 85 rows belonging to the final 15 RWDID components from the
#      user's original 206-row food classification / portion table.
#   2. Validate every listed item, Field ID, food-group assignment, encoding,
#      and portion value against the project metadata files.
#   3. Audit item-level raw responses across the original five instances i0-i4.
#   4. Reproduce the actual Step_02_1_1 transformations and group totals,
#      including the final surviving Bread roll pseudo-component only.
#   5. Compare the original NA-to-zero processing with a diagnostic scenario
#      that preserves an all-five-instance-missing item as NA.
#   6. Report participant/event counts, quintile/category shifts, and the
#      joint 15-component evaluable cohort in Full / Training / Test.
#
# STRICT SCOPE
#   - Source mapping: only the user's attached 206-row original food table;
#     this script freezes its 85 rows belonging to the final 15 groups.
#   - No new food groups, no external dietary indicators, no questionnaire
#     Gate fields, no p20082, and no alternative dietary indices.
#   - p20091-p20094 are used only as the existing Step_02_1_1 white-bread
#     type selectors needed to reproduce its pseudo-component. They are not
#     treated as additional food exposures or food groups.
#   - The final Step_02_1_1 script overwrites pw100950_i0-i4 sequentially.
#     The surviving assignment is p101160 (Bread roll) when p20094 == '[1]'.
#     The other three white-bread pseudo rows are audited but are NOT added
#     to the final Refined grains score.
#   - Portion grams are used exactly as in Step_02_1_1: for encoding_id == 0,
#     divide the raw amount by the first distinct portion_grams value when
#     exactly one value is available. Do NOT multiply frequency by grams.
#   - Non-continuous code recoding is reproduced exactly; rowMeans over the
#     available i0-i4 cells is retained. This script does not change the
#     original main score, quintiles, categories, beta weights, cohort split,
#     or Cox models.
#   - The NA-preserved score is a diagnostic sensitivity only, not a new
#     primary exposure and not a claim that all missing cells mean non-intake.
#
# INPUTS (run from the UKB project root)
#   ./data/ukb_data_dietary_org.RData
#   ./data/ukb_data_food_group.RData
#   ./data/Step_01_1_file_id_and_description_dietary.csv
#   ./data/Step_01_1_file_id_and_description_dietary_portion_data.csv
#   ./data/field.csv
#   ./results/Step_03/ukb_70_30_split/main/ukb_depression_training_cohort.rds
#   ./results/Step_03/ukb_70_30_split/main/ukb_depression_test_cohort.rds
#
# OUTPUT
#   ./results/Step_03/Step_03_2_7_final15_original_table_item_level_portion_audit_v1_0/
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
})

options(stringsAsFactors = FALSE, scipen = 999)

# -----------------------------------------------------------------------------
# 0. Paths and frozen cohort constants
# -----------------------------------------------------------------------------

RAW_RDATA <- "./data/ukb_data_dietary_org.RData"
FOOD_GROUP_RDATA <- "./data/ukb_data_food_group.RData"
DIETARY_METADATA <- "./data/Step_01_1_file_id_and_description_dietary.csv"
PORTION_METADATA <- "./data/Step_01_1_file_id_and_description_dietary_portion_data.csv"
FIELD_METADATA <- "./data/field.csv"

SPLIT_MAIN_DIR <- "./results/Step_03/ukb_70_30_split/main"
TRAIN_RDS <- file.path(SPLIT_MAIN_DIR, "ukb_depression_training_cohort.rds")
TEST_RDS <- file.path(SPLIT_MAIN_DIR, "ukb_depression_test_cohort.rds")

OUT_DIR <- paste0(
  "./results/Step_03/",
  "Step_03_2_7_final15_original_table_item_level_portion_audit_v1_0"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

INSTANCES <- 0:4
TOL <- 1e-8
ATTACHMENT_TOTAL_ROWS <- 206L
EXPECTED_FINAL15_ITEM_ROWS <- 85L
EXPECTED_EFFECTIVE_COMPONENTS <- 82L

# -----------------------------------------------------------------------------
# 1. Frozen final 15 RWDID components (original manuscript order)
# -----------------------------------------------------------------------------

food_specs <- data.frame(
  order = seq_len(15L),
  food_group = c(
    "Citrus fruits", "Red/Orange vegetables", "Other vegetables",
    "Potatoes", "Refined grains", "Eggs", "Low-fat dairy products",
    "Red meats", "Processed meats", "Fish and other seafood", "Coffee",
    "Fruit juices", "Sweetened beverages", "Water", "Wine"
  ),
  prefix = c(
    "citrus_fruits", "red_orange_veg", "other_veg", "potatoes",
    "refined_grains", "eggs", "low_fat_dairy", "red_meats",
    "processed_meats", "fish_seafood", "coffee", "fruit_juices",
    "sweetened_bev", "water", "wine"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    intake_var = prefix,
    quintile_var = paste0(prefix, "_quintile"),
    group_var = paste0(prefix, "_group"),
    category_var = paste0(prefix, "_group_range"),
    group_range_num_var = paste0(prefix, "_group_range_num")
  )

EXPECTED_FINAL15 <- c(
  "citrus_fruits", "red_orange_veg", "other_veg", "potatoes",
  "refined_grains", "eggs", "low_fat_dairy", "red_meats",
  "processed_meats", "fish_seafood", "coffee", "fruit_juices",
  "sweetened_bev", "water", "wine"
)

if (nrow(food_specs) != 15L || !identical(food_specs$prefix, EXPECTED_FINAL15)) {
  stop("STOP: Frozen final 15 component order is invalid.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 2. Frozen 85-row subset of the user's original 206-row table
#    Portion values are copied from that attachment, not re-estimated.
# -----------------------------------------------------------------------------

item_map <- data.frame(
  item_order = seq_len(85),
  category = c(
    "Fruit",
    "Fruit",
    "Fruit",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Vegetable",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Cereal & Bread",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Dairy & Egg",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Meat & Fish",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage",
    "Water & beverage"
  ),
  food_group = c(
    "Citrus fruits",
    "Citrus fruits",
    "Citrus fruits",
    "Red/Orange vegetables",
    "Red/Orange vegetables",
    "Red/Orange vegetables",
    "Red/Orange vegetables",
    "Red/Orange vegetables",
    "Red/Orange vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Other vegetables",
    "Potatoes",
    "Potatoes",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Refined grains",
    "Eggs",
    "Eggs",
    "Eggs",
    "Eggs",
    "Eggs",
    "Low-fat dairy products",
    "Low-fat dairy products",
    "Low-fat dairy products",
    "Low-fat dairy products",
    "Low-fat dairy products",
    "Low-fat dairy products",
    "Red meats",
    "Red meats",
    "Red meats",
    "Processed meats",
    "Processed meats",
    "Processed meats",
    "Fish and other seafood",
    "Fish and other seafood",
    "Fish and other seafood",
    "Fish and other seafood",
    "Fish and other seafood",
    "Fish and other seafood",
    "Fish and other seafood",
    "Coffee",
    "Coffee",
    "Coffee",
    "Coffee",
    "Coffee",
    "Coffee",
    "Fruit juices",
    "Fruit juices",
    "Fruit juices",
    "Sweetened beverages",
    "Sweetened beverages",
    "Sweetened beverages",
    "Sweetened beverages",
    "Sweetened beverages",
    "Sweetened beverages",
    "Sweetened beverages",
    "Water",
    "Wine",
    "Wine",
    "Wine"
  ),
  food_item = c(
    "Grapefruit intake",
    "Orange intake",
    "Satsuma intake",
    "Butternut squash intake",
    "Carrot intake",
    "Sweet pepper intake",
    "Sweet potato intake",
    "Fresh tomato intake",
    "Tinned tomato intake",
    "Mixed vegetable intake",
    "Vegetable pieces intake",
    "Coleslaw intake",
    "Side salad intake",
    "Avocado intake",
    "Beetroot intake",
    "Cauliflower intake",
    "Celery intake",
    "Courgette intake",
    "Cucumber intake",
    "Garlic intake",
    "Leek intake",
    "Onion intake",
    "Parsnip intake",
    "Turnip/swede intake",
    "Other vegetables intake",
    "Boiled/baked potatoes intake",
    "Mashed potato intake",
    "Naan bread intake",
    "Garlic bread intake",
    "Crispbread intake",
    "Oatcakes intake",
    "Other bread intake",
    "White pasta intake",
    "White rice intake",
    "Sushi intake",
    "Snackpot intake",
    "Couscous intake",
    "Sliced bread intake (white)",
    "Baguette intake (white)",
    "Bap intake (white)",
    "Bread roll intake (white)",
    "Whole egg intake",
    "Omelette intake",
    "Eggs in sandwiches intake",
    "Scotch egg intake",
    "Other egg intake",
    "Low fat yogurt",
    "Semi skimmed milk",
    "Skimmed milk and cholesterol-lowering milk",
    "Low fat hard cheese intake",
    "Low fat cheese spread intake",
    "Cottage cheese intake",
    "Beef intake",
    "Pork intake",
    "Lamb intake",
    "Sausage intake",
    "Bacon intake",
    "Ham intake",
    "Tinned tuna intake",
    "Oily fish intake",
    "White fish intake",
    "Other fish intake",
    "Prawns intake",
    "Lobster/crab intake",
    "Shellfish intake",
    "Instant coffee intake",
    "Filtered coffee intake",
    "Cappuccino intake",
    "Latte intake",
    "Espresso intake",
    "Other coffee type",
    "Orange juice intake",
    "Grapefruit juice intake",
    "Pure fruit/vegetable juice intake",
    "Fizzy drink intake",
    "Squash intake",
    "Fruit smoothie intake",
    "Dairy smoothie intake",
    "Hot chocolate intake",
    "Low calorie drink intake",
    "Low calorie hot chocolate intake",
    "Drinking water intake",
    "Fortified wine intake",
    "Red wine intake",
    "White wine intake"
  ),
  attachment_field_id = c(
    "104490",
    "104530",
    "104540",
    "104150",
    "104170",
    "104290",
    "104330",
    "104340",
    "104350",
    "104060",
    "104070",
    "104080",
    "104090",
    "104100",
    "104130",
    "104180",
    "104190",
    "104200",
    "104210",
    "104220",
    "104230",
    "104260",
    "104270",
    "104360",
    "104380",
    "104030",
    "104050",
    "101230",
    "101240",
    "101250",
    "101260",
    "101270",
    "102710",
    "102730",
    "102750",
    "102760",
    "102770",
    "w-100950",
    "w-101020",
    "w-101090",
    "w-101160",
    "102940",
    "102950",
    "102960",
    "102970",
    "102980",
    "26102",
    "26131",
    "26133",
    "102810",
    "102850",
    "102870",
    "103020",
    "103030",
    "103040",
    "103010",
    "103070",
    "103080",
    "103150",
    "103160",
    "103190",
    "103230",
    "103200",
    "103210",
    "103220",
    "100250",
    "100270",
    "100290",
    "100300",
    "100310",
    "100330",
    "100190",
    "100200",
    "100210",
    "100170",
    "100180",
    "100220",
    "100230",
    "100550",
    "100160",
    "100540",
    "100150",
    "26151",
    "26152",
    "26153"
  ),
  portion_grams_attachment = c(148, 120, 70, 130, 60, 133, 130, 85, 135, 90, 90, 90, 66, 136, 48, 90, 30, 90, 60, 5, 80, 60, 65, 60, 60, 175, 180, 160, 20, 10, 13, 45, 230, 150, 278, 280, 150, 36, 95, 90, 60, 50, 120, 120, 120, 75, 125, 259, 259, 40, 15, 60, 120, 120, 120, 30, 46, 23, 92, 100, 120, 100, 60, 85, 40, 190, 190, 30, 30, 30, 190, 250, 250, 250, 330, 250, 260, 260, 143, 330, 143, 250, 50, 175, 175),
  stringsAsFactors = FALSE
) %>%
  mutate(
    group_order = match(food_group, food_specs$food_group),
    item_type = ifelse(startsWith(attachment_field_id, "w-"), "white_bread_pseudo", "original_field"),
    source_field_id = ifelse(item_type == "white_bread_pseudo", sub("^w-", "", attachment_field_id), attachment_field_id),
    pseudo_selector_field_id = dplyr::case_when(
      source_field_id == "100950" & item_type == "white_bread_pseudo" ~ "20091",
      source_field_id == "101020" & item_type == "white_bread_pseudo" ~ "20092",
      source_field_id == "101090" & item_type == "white_bread_pseudo" ~ "20093",
      source_field_id == "101160" & item_type == "white_bread_pseudo" ~ "20094",
      TRUE ~ NA_character_
    ),
    step02_pseudo_id = ifelse(item_type == "white_bread_pseudo", paste0("w", source_field_id), NA_character_),
    effective_in_final_step02_score = item_type == "original_field" | source_field_id == "101160",
    effective_step02_slot = dplyr::case_when(
      item_type == "original_field" ~ paste0("p", source_field_id),
      source_field_id == "101160" ~ "pw100950 (overwritten slot; final content is Bread roll)",
      TRUE ~ "Not retained in final Step_02_1_1 score"
    ),
    final_step02_mapping_note = dplyr::case_when(
      item_type == "original_field" ~ "Direct original Field ID component",
      source_field_id == "101160" ~ "Effective Bread roll source: p101160 selected by p20094; stored in overwritten pw100950 slot",
      source_field_id == "100950" ~ "Not independently retained: pw100950 was overwritten by later assignments; final slot carries Bread roll",
      source_field_id %in% c("101020", "101090") ~ "No corresponding pw101020/pw101090 columns survive the final sequential overwrite; zero contribution",
      TRUE ~ "Review"
    )
  )

if (nrow(item_map) != EXPECTED_FINAL15_ITEM_ROWS) {
  stop("STOP: Expected exactly 85 final-15 rows from the original 206-row table.", call. = FALSE)
}
if (sum(item_map$item_type == "original_field") != 81L ||
    sum(item_map$item_type == "white_bread_pseudo") != 4L ||
    sum(item_map$effective_in_final_step02_score) != 82L) {
  stop("STOP: Expected 81 direct Field IDs, four white-bread pseudo rows, and 82 effective Step_02_1_1 components.", call. = FALSE)
}
if (anyNA(item_map$portion_grams_attachment) ||
    any(!is.finite(item_map$portion_grams_attachment)) ||
    any(item_map$portion_grams_attachment <= 0)) {
  stop("STOP: One or more attached portion values are missing or invalid.", call. = FALSE)
}
if (sum(item_map$item_type == "original_field") != 81L ||
    sum(item_map$item_type == "white_bread_pseudo") != 4L ||
    sum(item_map$effective_in_final_step02_score) != 82L) {
  stop("STOP: Expected 81 direct Field IDs, four white-bread pseudo rows, and 82 effective Step_02_1_1 components.", call. = FALSE)
}
if (anyDuplicated(item_map$attachment_field_id) > 0L) {
  stop("STOP: Duplicate attachment Field IDs/pseudo IDs in the final-15 mapping.", call. = FALSE)
}

expected_item_counts <- data.frame(
  food_group = food_specs$food_group,
  expected_items_n = c(3L, 6L, 16L, 2L, 14L, 5L, 6L, 3L, 3L, 7L, 6L, 3L, 7L, 1L, 3L),
  stringsAsFactors = FALSE
)
observed_item_counts <- item_map %>% count(food_group, name = "observed_items_n")
item_count_check <- expected_item_counts %>%
  left_join(observed_item_counts, by = "food_group") %>%
  mutate(observed_items_n = ifelse(is.na(observed_items_n), 0L, observed_items_n),
         count_matches = expected_items_n == observed_items_n)
write.csv(item_count_check,
          file.path(OUT_DIR, "00a_final15_attachment_item_count_check.csv"),
          row.names = FALSE, na = "")
if (any(!item_count_check$count_matches)) {
  stop("STOP: Original attachment item counts do not match the frozen 15-group specification.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 3. Utility functions
# -----------------------------------------------------------------------------

stop_now <- function(...) stop(paste0(...), call. = FALSE)

normalize_id <- function(x) {
  if (is.numeric(x)) {
    z <- format(x, scientific = FALSE, trim = TRUE, digits = 22)
  } else {
    z <- as.character(x)
  }
  z <- trimws(z)
  z[z == ""] <- NA_character_
  z <- sub("\\.0+$", "", z)
  z
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))

collapse_unique <- function(x, sep = " | ") {
  z <- unique(as.character(x[!is.na(x) & as.character(x) != ""]))
  if (length(z) == 0L) return(NA_character_)
  paste(z, collapse = sep)
}

load_named_rdata <- function(path, object_name) {
  if (!file.exists(path)) stop_now("Required file not found: ", path)
  e <- new.env(parent = globalenv())
  loaded_names <- load(path, envir = e)
  if (!object_name %in% loaded_names) {
    stop_now("Object '", object_name, "' not found in ", path,
             ". Objects found: ", paste(loaded_names, collapse = ", "))
  }
  get(object_name, envir = e, inherits = FALSE)
}

get_instance_columns <- function(data, field_id) {
  prefix <- paste0("p", field_id, "_i")
  cols <- names(data)[startsWith(names(data), prefix)]
  expected <- paste0("p", field_id, "_i", INSTANCES)
  cols <- expected[expected %in% names(data)]
  unexpected <- setdiff(names(data)[startsWith(names(data), paste0("p", field_id))], expected)
  if (length(unexpected) > 0L) {
    stop_now("Unexpected raw columns for Field ID ", field_id, ": ",
             paste(unexpected, collapse = " | "))
  }
  cols
}

get_selector_columns <- function(data, field_id) {
  expected <- paste0("p", field_id, "_i", INSTANCES)
  missing <- setdiff(expected, names(data))
  if (length(missing) > 0L) {
    stop_now("Missing original Step_02_1_1 bread-type selector column(s): ",
             paste(missing, collapse = " | "))
  }
  expected
}

to_numeric_checked <- function(x, column_name) {
  if (is.numeric(x)) return(as.numeric(x))
  original_nonmissing <- !is.na(x)
  converted <- suppressWarnings(as.numeric(as.character(x)))
  bad <- original_nonmissing & is.na(converted)
  if (any(bad)) {
    examples <- unique(as.character(x[bad]))
    stop_now("Non-numeric values in raw dietary column ", column_name,
             ". Examples: ", paste(utils::head(examples, 10L), collapse = " | "))
  }
  converted
}

transform_field_block <- function(raw_data, field_id, encoding_id, portion_grams) {
  cols <- get_instance_columns(raw_data, field_id)
  if (length(cols) == 0L) stop_now("No i0-i4 raw columns found for Field ID ", field_id, ".")

  raw_block <- raw_data[, cols, drop = FALSE]
  transformed <- vector("list", length(cols))
  names(transformed) <- cols

  for (j in seq_along(cols)) {
    col <- cols[j]
    x <- to_numeric_checked(raw_block[[col]], col)

    if (!(encoding_id %in% c(0, 80))) {
      x[x == 444] <- 0.25
      x[x == 555] <- 0.5
      x[x == 200] <- 2
      x[x == 300] <- 3
      x[x == 400] <- 4
      x[x == 500] <- 5
      x[x == 600] <- 6
    } else if (encoding_id == 0) {
      if (length(portion_grams) == 1L && !is.na(portion_grams) &&
          is.finite(portion_grams) && portion_grams > 0) {
        x <- x / portion_grams
      }
    }
    # encoding_id == 80 remains unchanged, exactly as Step_02_1_1.
    transformed[[j]] <- x
  }

  transformed_block <- as.data.frame(transformed, check.names = FALSE)
  component_mean <- rowMeans(transformed_block, na.rm = TRUE)
  component_original <- component_mean
  component_original[is.na(component_original)] <- 0
  component_na_preserved <- component_mean

  raw_matrix <- as.matrix(raw_block)
  transformed_matrix <- as.matrix(transformed_block)
  list(
    raw_block = raw_block,
    transformed_block = transformed_block,
    component_original = component_original,
    component_na_preserved = component_na_preserved,
    raw_observed_n = rowSums(!is.na(raw_matrix)),
    transformed_observed_n = rowSums(!is.na(transformed_matrix)),
    transformed_positive_n = rowSums(transformed_matrix > 0, na.rm = TRUE),
    transformed_zero_n = rowSums(transformed_matrix == 0, na.rm = TRUE),
    transformed_negative_n = rowSums(transformed_matrix < 0, na.rm = TRUE),
    columns = cols
  )
}

calculate_event_rate <- function(events, non_events) {
  denominator <- events + non_events
  ifelse(denominator > 0, 100 * events / denominator, NA_real_)
}

write_tsv_console <- function(x) {
  utils::write.table(as.data.frame(x), file = "", sep = "\t",
                     row.names = FALSE, col.names = TRUE,
                     quote = FALSE, na = "NA")
}

compare_numeric <- function(a, b, tolerance = TOL) {
  a <- as_numeric_safe(a)
  b <- as_numeric_safe(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tolerance)
}

# -----------------------------------------------------------------------------
# 4. Load original data, derived groups, and frozen MAIN cohorts
# -----------------------------------------------------------------------------

message("[1/10] Loading raw dietary data, derived food groups, and frozen cohorts ...")
raw_data <- load_named_rdata(RAW_RDATA, "ukb_data")
food_data <- load_named_rdata(FOOD_GROUP_RDATA, "ukb_data_food_group")
if (!file.exists(TRAIN_RDS)) stop_now("Training cohort RDS not found: ", TRAIN_RDS)
if (!file.exists(TEST_RDS)) stop_now("Test cohort RDS not found: ", TEST_RDS)
train_data <- readRDS(TRAIN_RDS)
test_data <- readRDS(TEST_RDS)

required_raw <- "eid"
if (!required_raw %in% names(raw_data) || !required_raw %in% names(food_data) ||
    !all(c("eid", "depression_status") %in% names(train_data)) ||
    !all(c("eid", "depression_status") %in% names(test_data))) {
  stop_now("One or more source objects lack required EID/depression_status fields.")
}

raw_data$eid_audit <- normalize_id(raw_data$eid)
food_data$eid_audit <- normalize_id(food_data$eid)
train_data$eid_audit <- normalize_id(train_data$eid)
test_data$eid_audit <- normalize_id(test_data$eid)
train_data$depression_status <- as_numeric_safe(train_data$depression_status)
test_data$depression_status <- as_numeric_safe(test_data$depression_status)

for (object_name in c("raw_data", "food_data", "train_data", "test_data")) {
  object_value <- get(object_name)
  if (anyNA(object_value$eid_audit) || anyDuplicated(object_value$eid_audit) > 0L) {
    stop_now(object_name, " contains missing or duplicated EIDs.")
  }
}
if (length(intersect(train_data$eid_audit, test_data$eid_audit)) > 0L) {
  stop_now("Training and test EIDs overlap.")
}
full_data <- bind_rows(train_data, test_data)
if (anyDuplicated(full_data$eid_audit) > 0L) stop_now("Full cohort contains duplicated EIDs.")
if (length(setdiff(food_data$eid_audit, raw_data$eid_audit)) > 0L) {
  stop_now("Some derived food-group EIDs are absent from raw dietary data.")
}
if (length(setdiff(full_data$eid_audit, food_data$eid_audit)) > 0L) {
  stop_now("Some frozen cohort EIDs are absent from derived food-group data.")
}
if (any(!is.na(full_data$depression_status) & !(full_data$depression_status %in% c(0, 1)))) {
  stop_now("depression_status contains values other than 0, 1, or NA.")
}

raw_index <- match(food_data$eid_audit, raw_data$eid_audit)
if (anyNA(raw_index)) stop_now("Raw/derived EID matching failed.")
raw_cohort <- raw_data[raw_index, , drop = FALSE]
if (!identical(raw_cohort$eid_audit, food_data$eid_audit)) {
  stop_now("Raw/derived EID ordering validation failed.")
}

cohort_event_summary <- bind_rows(
  full_data %>% mutate(sample = "Full"),
  train_data %>% mutate(sample = "Training"),
  test_data %>% mutate(sample = "Test")
) %>%
  group_by(sample) %>%
  summarise(
    participants_n = n(),
    depression_events_n = sum(depression_status == 1, na.rm = TRUE),
    non_events_n = sum(depression_status == 0, na.rm = TRUE),
    outcome_missing_n = sum(is.na(depression_status)),
    event_rate_percent = calculate_event_rate(depression_events_n, non_events_n),
    .groups = "drop"
  )

expected_cohorts <- data.frame(
  sample = c("Full", "Training", "Test"),
  participants_n = c(109397L, 76578L, 32819L),
  depression_events_n = c(2997L, 2122L, 875L),
  stringsAsFactors = FALSE
)
cohort_check <- cohort_event_summary %>%
  left_join(expected_cohorts, by = "sample", suffix = c("", "_expected")) %>%
  mutate(counts_match = participants_n == participants_n_expected &
           depression_events_n == depression_events_n_expected & outcome_missing_n == 0L)
write.csv(cohort_check, file.path(OUT_DIR, "00_frozen_cohort_count_check.csv"),
          row.names = FALSE, na = "")
if (any(!cohort_check$counts_match)) {
  stop_now("Frozen MAIN cohort/event counts do not match the established 70/30 split.")
}

# -----------------------------------------------------------------------------
# 5. Read metadata and validate all 85 attachment rows
# -----------------------------------------------------------------------------

message("[2/10] Validating attachment Field IDs, descriptions, encodings, and portions ...")
required_metadata_files <- c(DIETARY_METADATA, PORTION_METADATA, FIELD_METADATA)
missing_metadata_files <- required_metadata_files[!file.exists(required_metadata_files)]
if (length(missing_metadata_files) > 0L) {
  stop_now("Missing metadata file(s): ", paste(missing_metadata_files, collapse = "; "))
}

dietary_dict_raw <- read.csv(DIETARY_METADATA, stringsAsFactors = FALSE, check.names = FALSE)
portion_raw <- read.csv(PORTION_METADATA, stringsAsFactors = FALSE, check.names = FALSE)
field_raw <- read.csv(FIELD_METADATA, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("field_id", "description", "category", "food_group") %in% names(dietary_dict_raw))) {
  stop_now("Dietary dictionary must contain field_id, description, category, food_group.")
}
if (!all(c("field_id", "portion_grams") %in% names(portion_raw))) {
  stop_now("Portion metadata must contain field_id and portion_grams.")
}
if (!all(c("field_id", "encoding_id") %in% names(field_raw))) {
  stop_now("field.csv must contain field_id and encoding_id.")
}

dietary_dict <- dietary_dict_raw %>%
  mutate(
    field_id_chr = normalize_id(field_id),
    description_norm = tolower(trimws(as.character(description))),
    category_norm = tolower(trimws(as.character(category))),
    food_group_norm = tolower(trimws(as.character(food_group)))
  )

portion_raw <- portion_raw %>%
  mutate(
    field_id_chr = normalize_id(field_id),
    portion_grams_num = suppressWarnings(as.numeric(as.character(portion_grams)))
  )
# This is the exact first-row de-duplication used in Step_02_1_1.
portion_first <- portion_raw %>% distinct(field_id_chr, .keep_all = TRUE)
portion_duplicate_summary <- portion_raw %>%
  group_by(field_id_chr) %>%
  summarise(
    portion_rows_n = n(),
    portion_values_n = n_distinct(portion_grams_num[!is.na(portion_grams_num)]),
    portion_values = collapse_unique(portion_grams_num),
    .groups = "drop"
  )

field_meta <- field_raw %>%
  mutate(
    field_id_chr = normalize_id(field_id),
    encoding_id_num = suppressWarnings(as.numeric(as.character(encoding_id)))
  ) %>%
  distinct(field_id_chr, .keep_all = TRUE) %>%
  select(field_id_chr, encoding_id = encoding_id_num)

item_map <- item_map %>%
  mutate(
    source_portion_grams = portion_first$portion_grams_num[
      match(source_field_id, portion_first$field_id_chr)
    ],
    encoding_id = field_meta$encoding_id[
      match(source_field_id, field_meta$field_id_chr)
    ],
    raw_columns = vapply(source_field_id,
                         function(fid) paste(get_instance_columns(raw_cohort, fid), collapse = " | "),
                         character(1)),
    raw_columns_n = vapply(source_field_id,
                           function(fid) length(get_instance_columns(raw_cohort, fid)),
                           integer(1)),
    portion_metadata_match = !is.na(source_portion_grams) &
      abs(source_portion_grams - portion_grams_attachment) <= TOL,
    portion_duplicate_values_consistent = vapply(source_field_id, function(fid) {
      z <- portion_duplicate_summary$portion_values_n[
        match(fid, portion_duplicate_summary$field_id_chr)
      ]
      !is.na(z) && z <= 1L
    }, logical(1)),
    field_encoding_present = !is.na(encoding_id)
  )

# Exact dictionary matching for the 81 ordinary Field ID rows.
item_map$dictionary_group_description_match <- NA
item_map$dictionary_category_match <- NA
for (i in seq_len(nrow(item_map))) {
  if (item_map$item_type[i] == "original_field") {
    hit <- dietary_dict$field_id_chr == item_map$source_field_id[i] &
      dietary_dict$food_group_norm == tolower(trimws(item_map$food_group[i])) &
      dietary_dict$description_norm == tolower(trimws(item_map$food_item[i]))
    item_map$dictionary_group_description_match[i] <- any(hit, na.rm = TRUE)
    if (any(hit, na.rm = TRUE)) {
      item_map$dictionary_category_match[i] <- any(
        dietary_dict$category_norm[hit] == tolower(trimws(item_map$category[i])),
        na.rm = TRUE
      )
    } else {
      item_map$dictionary_category_match[i] <- FALSE
    }
  }
}

# Pseudo white-bread rows are checked against their underlying amount Field IDs
# and portion/encoding metadata; they are not separate numeric UKB Field IDs.
item_map$source_field_exists_in_dictionary <- vapply(
  item_map$source_field_id,
  function(fid) any(dietary_dict$field_id_chr == fid, na.rm = TRUE),
  logical(1)
)

# All four selectors belong to the original Step_02_1_1 white-bread mapping.
selector_presence_rows <- list()
for (i in which(item_map$item_type == "white_bread_pseudo")) {
  selector_cols <- get_selector_columns(raw_cohort, item_map$pseudo_selector_field_id[i])
  selector_presence_rows[[length(selector_presence_rows) + 1L]] <- data.frame(
    attachment_field_id = item_map$attachment_field_id[i],
    food_item = item_map$food_item[i],
    selector_field_id = item_map$pseudo_selector_field_id[i],
    selector_columns_n = length(selector_cols),
    selector_columns = paste(selector_cols, collapse = " | "),
    stringsAsFactors = FALSE
  )
}
selector_presence_summary <- bind_rows(selector_presence_rows)

mapping_validation <- item_map %>%
  select(item_order, group_order, category, food_group, food_item,
         attachment_field_id, source_field_id, item_type,
         pseudo_selector_field_id, step02_pseudo_id,
         portion_grams_attachment, source_portion_grams, portion_metadata_match,
         portion_duplicate_values_consistent, encoding_id, field_encoding_present,
         raw_columns_n, raw_columns, dictionary_group_description_match,
         dictionary_category_match, source_field_exists_in_dictionary,
         effective_in_final_step02_score, effective_step02_slot,
         final_step02_mapping_note)

write.csv(mapping_validation,
          file.path(OUT_DIR, "02_attachment_85_item_mapping_and_metadata_validation.csv"),
          row.names = FALSE, na = "")
write.csv(selector_presence_summary,
          file.path(OUT_DIR, "02a_original_step02_bread_selector_dependency_audit.csv"),
          row.names = FALSE, na = "")
write.csv(portion_duplicate_summary,
          file.path(OUT_DIR, "02b_portion_metadata_duplicate_audit.csv"),
          row.names = FALSE, na = "")

bad_direct_dictionary <- item_map %>%
  filter(item_type == "original_field",
         is.na(dictionary_group_description_match) | !dictionary_group_description_match |
           is.na(dictionary_category_match) | !dictionary_category_match)
bad_metadata <- item_map %>%
  filter(!portion_metadata_match | !portion_duplicate_values_consistent |
           !field_encoding_present | raw_columns_n == 0L | !source_field_exists_in_dictionary)
if (nrow(bad_direct_dictionary) > 0L || nrow(bad_metadata) > 0L) {
  write.csv(bad_direct_dictionary,
            file.path(OUT_DIR, "02c_original_table_dictionary_mismatch_review.csv"),
            row.names = FALSE, na = "")
  write.csv(bad_metadata,
            file.path(OUT_DIR, "02d_portion_encoding_or_raw_column_mismatch_review.csv"),
            row.names = FALSE, na = "")
  stop_now("Attachment-to-project metadata validation failed. Review files 02c and 02d before proceeding.")
}

# All four pseudo rows must have the expected raw amount columns and selectors.
if (nrow(selector_presence_summary) != 4L || any(selector_presence_summary$selector_columns_n != 5L)) {
  stop_now("The four original white-bread selector mappings are incomplete.")
}

# -----------------------------------------------------------------------------
# 6. Transform all direct fields and white-bread pseudo-item diagnostics
# -----------------------------------------------------------------------------

message("[3/10] Reproducing original Step_02_1_1 item transformations ...")

component_original <- list()
component_na_preserved <- list()
item_status_vectors <- list()
item_level_rows <- list()
item_level_k <- 0L

# 6A. 81 ordinary original Field ID items.
direct_map <- item_map %>% filter(item_type == "original_field")
for (i in seq_len(nrow(direct_map))) {
  row <- direct_map[i, ]
  fid <- row$source_field_id
  transformed <- transform_field_block(
    raw_data = raw_cohort,
    field_id = fid,
    encoding_id = row$encoding_id,
    portion_grams = row$source_portion_grams
  )
  key <- paste0("direct_", fid)
  component_original[[key]] <- transformed$component_original
  component_na_preserved[[key]] <- transformed$component_na_preserved

  status <- ifelse(
    transformed$raw_observed_n == 0L,
    "ALL_FIVE_INSTANCES_MISSING",
    ifelse(transformed$raw_observed_n == length(INSTANCES),
           "ALL_FIVE_INSTANCES_OBSERVED", "PARTIAL_INSTANCE_RESPONSE")
  )
  item_status_vectors[[as.character(row$item_order)]] <- status

  item_level_k <- item_level_k + 1L
  item_level_rows[[item_level_k]] <- data.frame(
    item_order = row$item_order,
    group_order = row$group_order,
    food_group = row$food_group,
    food_item = row$food_item,
    attachment_field_id = row$attachment_field_id,
    source_field_id = fid,
    item_type = "original_field",
    effective_in_final_step02_score = TRUE,
    encoding_id = row$encoding_id,
    portion_grams_attachment = row$portion_grams_attachment,
    raw_instance_columns_n = length(transformed$columns),
    raw_instance_columns = paste(transformed$columns, collapse = " | "),
    participants_n = nrow(raw_cohort),
    participants_any_raw_response_n = sum(transformed$raw_observed_n > 0L),
    participants_all_five_instances_missing_n = sum(transformed$raw_observed_n == 0L),
    participants_partial_instance_response_n = sum(
      transformed$raw_observed_n > 0L & transformed$raw_observed_n < length(INSTANCES)
    ),
    raw_cells_observed_n = sum(transformed$raw_observed_n),
    raw_cells_missing_n = nrow(raw_cohort) * length(transformed$columns) - sum(transformed$raw_observed_n),
    transformed_positive_cells_n = sum(transformed$transformed_positive_n),
    transformed_zero_cells_n = sum(transformed$transformed_zero_n),
    transformed_negative_cells_n = sum(transformed$transformed_negative_n),
    participants_component_zero_original_n = sum(transformed$component_original == 0, na.rm = TRUE),
    participants_component_positive_original_n = sum(transformed$component_original > 0, na.rm = TRUE),
    participants_component_missing_na_preserved_n = sum(is.na(transformed$component_na_preserved)),
    participant_response_status_levels = collapse_unique(status),
    stringsAsFactors = FALSE
  )
}

# 6B. Four white-bread rows in the attachment. Their source amount fields and
#     selectors are inspected as item diagnostics. Only the Bread roll row
#     (w-101160) is used in the final Step_02_1_1 Refined grains score.
pseudo_map <- item_map %>% filter(item_type == "white_bread_pseudo")
pseudo_item_component_original <- list()
pseudo_item_component_na <- list()

for (i in seq_len(nrow(pseudo_map))) {
  row <- pseudo_map[i, ]
  fid <- row$source_field_id
  selector_fid <- row$pseudo_selector_field_id
  amount <- transform_field_block(
    raw_data = raw_cohort,
    field_id = fid,
    encoding_id = row$encoding_id,
    portion_grams = row$source_portion_grams
  )
  selector_cols <- get_selector_columns(raw_cohort, selector_fid)
  selector_block <- raw_cohort[, selector_cols, drop = FALSE]
  selected_matrix <- matrix(NA_real_, nrow = nrow(raw_cohort), ncol = length(INSTANCES))
  selector_yes_matrix <- matrix(FALSE, nrow = nrow(raw_cohort), ncol = length(INSTANCES))
  selector_observed_matrix <- matrix(FALSE, nrow = nrow(raw_cohort), ncol = length(INSTANCES))
  amount_missing_when_selected_matrix <- matrix(FALSE, nrow = nrow(raw_cohort), ncol = length(INSTANCES))

  for (j in seq_along(INSTANCES)) {
    selector_value <- as.character(selector_block[[j]])
    selector_observed_matrix[, j] <- !is.na(selector_value)
    selector_yes_matrix[, j] <- !is.na(selector_value) & selector_value == "[1]"
    selected_matrix[, j] <- ifelse(
      selector_yes_matrix[, j],
      amount$transformed_block[[j]],
      NA_real_
    )
    amount_missing_when_selected_matrix[, j] <-
      selector_yes_matrix[, j] & is.na(amount$transformed_block[[j]])
  }

  selected_block <- as.data.frame(selected_matrix)
  selected_mean <- rowMeans(selected_block, na.rm = TRUE)
  selected_original <- selected_mean
  selected_original[is.na(selected_original)] <- 0

  selector_observed_n <- rowSums(selector_observed_matrix)
  selector_yes_n <- rowSums(selector_yes_matrix)
  selected_amount_observed_n <- rowSums(!is.na(selected_matrix))
  selected_amount_missing_n <- rowSums(amount_missing_when_selected_matrix)

  # NA-preserved diagnostic: if every selector is observed and none selects
  # this item, the item contribution is a supported zero. If selection or
  # selector status is unresolved and no selected amount is observed, retain NA.
  selected_na_preserved <- selected_mean
  known_not_selected <- selector_observed_n == length(INSTANCES) & selector_yes_n == 0L
  selected_na_preserved[known_not_selected] <- 0
  unresolved_no_amount <- selected_amount_observed_n == 0L & !known_not_selected
  selected_na_preserved[unresolved_no_amount] <- NA_real_

  pseudo_item_component_original[[fid]] <- selected_original
  pseudo_item_component_na[[fid]] <- selected_na_preserved

  used_in_score <- isTRUE(row$effective_in_final_step02_score)
  if (used_in_score) {
    # Final Step_02_1_1 assigns this Bread roll vector to pw100950_i0-i4.
    component_original[["refined_grains_final_bread_roll_pw100950"]] <- selected_original
    component_na_preserved[["refined_grains_final_bread_roll_pw100950"]] <- selected_na_preserved
  }

  status <- ifelse(
    selector_observed_n == 0L,
    "SELECTOR_ALL_FIVE_INSTANCES_MISSING",
    ifelse(
      selected_amount_observed_n > 0L,
      "SELECTED_AMOUNT_OBSERVED",
      ifelse(
        known_not_selected,
        "SELECTOR_OBSERVED_ITEM_NOT_SELECTED",
        ifelse(selected_amount_missing_n > 0L,
               "SELECTED_BUT_AMOUNT_MISSING", "SELECTOR_PARTIAL_OR_UNRESOLVED")
      )
    )
  )
  item_status_vectors[[as.character(row$item_order)]] <- status

  item_level_k <- item_level_k + 1L
  item_level_rows[[item_level_k]] <- data.frame(
    item_order = row$item_order,
    group_order = row$group_order,
    food_group = row$food_group,
    food_item = row$food_item,
    attachment_field_id = row$attachment_field_id,
    source_field_id = fid,
    item_type = "white_bread_pseudo",
    effective_in_final_step02_score = used_in_score,
    encoding_id = row$encoding_id,
    portion_grams_attachment = row$portion_grams_attachment,
    raw_instance_columns_n = length(amount$columns),
    raw_instance_columns = paste(amount$columns, collapse = " | "),
    selector_field_id = selector_fid,
    selector_columns_n = length(selector_cols),
    selector_columns = paste(selector_cols, collapse = " | "),
    participants_n = nrow(raw_cohort),
    participants_any_raw_response_n = sum(amount$raw_observed_n > 0L),
    participants_all_five_instances_missing_n = sum(amount$raw_observed_n == 0L),
    participants_partial_instance_response_n = sum(
      amount$raw_observed_n > 0L & amount$raw_observed_n < length(INSTANCES)
    ),
    raw_cells_observed_n = sum(amount$raw_observed_n),
    raw_cells_missing_n = nrow(raw_cohort) * length(amount$columns) - sum(amount$raw_observed_n),
    transformed_positive_cells_n = sum(amount$transformed_positive_n),
    transformed_zero_cells_n = sum(amount$transformed_zero_n),
    transformed_negative_cells_n = sum(amount$transformed_negative_n),
    participants_selector_any_observed_n = sum(selector_observed_n > 0L),
    participants_selector_all_five_missing_n = sum(selector_observed_n == 0L),
    participants_selector_selected_at_least_once_n = sum(selector_yes_n > 0L),
    participants_selected_amount_observed_n = sum(selected_amount_observed_n > 0L),
    selected_amount_observed_cells_n = sum(selected_amount_observed_n),
    selected_but_amount_missing_cells_n = sum(selected_amount_missing_n),
    participants_component_zero_original_n = sum(selected_original == 0, na.rm = TRUE),
    participants_component_positive_original_n = sum(selected_original > 0, na.rm = TRUE),
    participants_component_missing_na_preserved_n = sum(is.na(selected_na_preserved)),
    participant_response_status_levels = collapse_unique(status),
    final_step02_mapping_note = row$final_step02_mapping_note,
    stringsAsFactors = FALSE
  )
}

item_level_summary <- bind_rows(item_level_rows) %>% arrange(item_order)
if (nrow(item_level_summary) != EXPECTED_FINAL15_ITEM_ROWS) {
  stop_now("Item-level response audit did not produce exactly 85 rows.")
}

# -----------------------------------------------------------------------------
# 7. Reconstruct each final 15 group exactly and build NA-preserved diagnostic
# -----------------------------------------------------------------------------

message("[4/10] Reconstructing the 15 original group totals ...")

# Group-specific effective item keys: 81 direct rows plus only the final
# Bread roll pseudo-component for Refined grains (82 effective components).
effective_item_map <- item_map %>% filter(effective_in_final_step02_score)
if (nrow(effective_item_map) != EXPECTED_EFFECTIVE_COMPONENTS) {
  stop_now("Expected 82 effective Step_02_1_1 components (81 direct + final Bread roll).")
}

# The final Step_02_1_1 Refined grains routine passes four pseudo IDs. Only
# pw100950 has source columns, and its final contents are the Bread roll
# assignment. The other three pseudo slots are zero after the original
# replace_na(..., 0) behavior; they are recorded, not reclassified.
pseudo_slot_audit <- data.frame(
  pseudo_id = c("w100950", "w101020", "w101090", "w101160"),
  matching_pw_columns = c("pw100950_i0-pw100950_i4", NA, NA, NA),
  final_slot_status = c(
    "PW100950_EXISTS_BUT_FINAL_CONTENT_IS_BREAD_ROLL_FROM_P101160",
    "NO_PW101020_COLUMNS_AFTER_OVERWRITE_ZERO_CONTRIBUTION",
    "NO_PW101090_COLUMNS_AFTER_OVERWRITE_ZERO_CONTRIBUTION",
    "NO_PW101160_COLUMNS_AFTER_OVERWRITE_ZERO_CONTRIBUTION"
  ),
  effective_attachment_item = c("w-101160", NA, NA, NA),
  stringsAsFactors = FALSE
)

# Build a stable vector of component keys for each food group.
group_component_keys <- vector("list", nrow(food_specs))
group_reconstruction_rows <- vector("list", nrow(food_specs))
group_score_original <- data.frame(eid_audit = food_data$eid_audit, stringsAsFactors = FALSE)
group_score_na <- data.frame(eid_audit = food_data$eid_audit, stringsAsFactors = FALSE)
group_missing_component_n <- data.frame(eid_audit = food_data$eid_audit, stringsAsFactors = FALSE)
group_status_vectors <- list()

for (g in seq_len(nrow(food_specs))) {
  sp <- food_specs[g, ]
  rows_g <- effective_item_map %>% filter(food_group == sp$food_group)
  keys <- character(0)
  for (j in seq_len(nrow(rows_g))) {
    row <- rows_g[j, ]
    if (row$item_type == "original_field") {
      keys <- c(keys, paste0("direct_", row$source_field_id))
    } else if (row$attachment_field_id == "w-101160") {
      keys <- c(keys, "refined_grains_final_bread_roll_pw100950")
    }
  }
  if (length(keys) == 0L) stop_now("No effective source components found for ", sp$food_group)
  group_component_keys[[g]] <- keys

  original_matrix <- do.call(cbind, lapply(keys, function(k) component_original[[k]]))
  strict_matrix <- do.call(cbind, lapply(keys, function(k) component_na_preserved[[k]]))
  if (is.null(dim(original_matrix))) original_matrix <- matrix(original_matrix, ncol = 1L)
  if (is.null(dim(strict_matrix))) strict_matrix <- matrix(strict_matrix, ncol = 1L)

  reconstructed_original <- rowSums(original_matrix, na.rm = TRUE)
  reconstructed_strict <- rowSums(strict_matrix, na.rm = FALSE)
  missing_component_n <- rowSums(is.na(strict_matrix))
  strict_evaluable <- missing_component_n == 0L

  saved_total <- as_numeric_safe(food_data[[sp$intake_var]])
  difference <- reconstructed_original - saved_total
  match_flag <- compare_numeric(reconstructed_original, saved_total)

  group_score_original[[sp$prefix]] <- reconstructed_original
  group_score_na[[sp$prefix]] <- reconstructed_strict
  group_missing_component_n[[paste0(sp$prefix, "_missing_item_n")]] <- missing_component_n
  group_status_vectors[[sp$prefix]] <- ifelse(
    strict_evaluable,
    "ALL_EFFECTIVE_ITEMS_HAVE_EVALUABLE_SOURCE_STATUS",
    "ONE_OR_MORE_EFFECTIVE_ITEMS_ALL_INSTANCES_MISSING_OR_UNRESOLVED"
  )

  group_reconstruction_rows[[g]] <- data.frame(
    order = sp$order,
    food_group = sp$food_group,
    prefix = sp$prefix,
    attachment_table_rows_n = sum(item_map$food_group == sp$food_group),
    effective_step02_components_n = length(keys),
    effective_component_keys = paste(keys, collapse = ";"),
    participants_n = nrow(food_data),
    saved_original_zero_n = sum(saved_total == 0, na.rm = TRUE),
    saved_original_positive_n = sum(saved_total > 0, na.rm = TRUE),
    reconstructed_original_zero_n = sum(reconstructed_original == 0, na.rm = TRUE),
    reconstructed_original_positive_n = sum(reconstructed_original > 0, na.rm = TRUE),
    na_preserved_diagnostic_evaluable_n = sum(strict_evaluable),
    na_preserved_diagnostic_missing_n = sum(!strict_evaluable),
    participants_any_effective_item_all_five_missing_n = sum(missing_component_n > 0L),
    reconstruction_match_n = sum(match_flag),
    reconstruction_mismatch_n = sum(!match_flag),
    max_absolute_difference = if (any(is.finite(difference))) max(abs(difference[is.finite(difference)])) else NA_real_,
    reconstruction_status = if (all(match_flag)) "PASS" else "REVIEW_REQUIRED",
    stringsAsFactors = FALSE
  )
}

group_reconstruction_summary <- bind_rows(group_reconstruction_rows)
write.csv(group_reconstruction_summary,
          file.path(OUT_DIR, "06_final15_group_reconstruction_and_na_preserved_summary.csv"),
          row.names = FALSE, na = "")
write.csv(pseudo_slot_audit,
          file.path(OUT_DIR, "04_refined_grains_final_step02_pseudo_slot_audit.csv"),
          row.names = FALSE, na = "")

if (any(group_reconstruction_summary$reconstruction_status != "PASS")) {
  stop_now("One or more final-15 group totals do not reproduce the saved Step_02_1_1 totals.")
}

# -----------------------------------------------------------------------------
# 8. Create diagnostic quintiles/categories without changing frozen originals
# -----------------------------------------------------------------------------

message("[5/10] Calculating diagnostic-only quintiles and category shifts ...")

diagnostic_coding <- data.frame(eid_audit = food_data$eid_audit, stringsAsFactors = FALSE)
quintile_shift_rows <- vector("list", nrow(food_specs))

for (g in seq_len(nrow(food_specs))) {
  sp <- food_specs[g, ]
  strict_score <- group_score_na[[sp$prefix]]
  original_quintile <- as.character(food_data[[sp$quintile_var]])
  original_range_num <- as.character(food_data[[sp$group_range_num_var]])
  diagnostic_quintile <- rep(NA_character_, length(strict_score))
  evaluable <- !is.na(strict_score)
  if (sum(evaluable) > 0L) {
    diagnostic_quintile[evaluable] <- paste0("Q", dplyr::ntile(strict_score[evaluable], 5L))
  }

  diagnostic_range <- rep(NA_character_, length(strict_score))
  diagnostic_range_num <- rep(NA_integer_, length(strict_score))
  if (sum(evaluable) > 0L) {
    temp <- data.frame(
      q = diagnostic_quintile[evaluable],
      score = strict_score[evaluable],
      stringsAsFactors = FALSE
    )
    q_medians <- temp %>%
      group_by(q) %>%
      summarise(median_val = median(score, na.rm = TRUE), .groups = "drop") %>%
      arrange(median_val) %>%
      mutate(group = paste0("G", row_number()))
    temp <- temp %>% left_join(q_medians %>% select(q, group), by = "q")
    group_stats <- temp %>%
      group_by(group) %>%
      summarise(
        median_val = median(score, na.rm = TRUE),
        min_val = min(score, na.rm = TRUE),
        max_val = max(score, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(median_val) %>%
      mutate(range = paste0(round(median_val, 2), " (", round(min_val, 2), "-", round(max_val, 2), ")"))
    temp <- temp %>% left_join(group_stats %>% select(group, range), by = "group")
    unique_ranges <- unique(group_stats$range)
    diagnostic_range[evaluable] <- temp$range
    diagnostic_range_num[evaluable] <- match(temp$range, unique_ranges)
  }

  diagnostic_coding[[paste0(sp$prefix, "_diagnostic_quintile")]] <- diagnostic_quintile
  diagnostic_coding[[paste0(sp$prefix, "_diagnostic_range")]] <- diagnostic_range
  diagnostic_coding[[paste0(sp$prefix, "_diagnostic_range_num")]] <- diagnostic_range_num

  q_comparable <- evaluable & !is.na(original_quintile) & !is.na(diagnostic_quintile)
  c_comparable <- evaluable & !is.na(original_range_num) & !is.na(diagnostic_range_num)
  quintile_shift_rows[[g]] <- data.frame(
    food_group = sp$food_group,
    prefix = sp$prefix,
    participants_full_food_data_n = nrow(food_data),
    na_preserved_diagnostic_evaluable_n = sum(evaluable),
    na_preserved_diagnostic_missing_n = sum(!evaluable),
    original_and_diagnostic_quintile_comparable_n = sum(q_comparable),
    diagnostic_quintile_changed_n = sum(original_quintile[q_comparable] != diagnostic_quintile[q_comparable]),
    diagnostic_quintile_changed_percent = if (sum(q_comparable) > 0L) {
      100 * mean(original_quintile[q_comparable] != diagnostic_quintile[q_comparable])
    } else NA_real_,
    original_and_diagnostic_category_comparable_n = sum(c_comparable),
    diagnostic_category_changed_n = sum(original_range_num[c_comparable] != as.character(diagnostic_range_num[c_comparable])),
    diagnostic_category_changed_percent = if (sum(c_comparable) > 0L) {
      100 * mean(original_range_num[c_comparable] != as.character(diagnostic_range_num[c_comparable]))
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}

quintile_shift_summary <- bind_rows(quintile_shift_rows)
write.csv(quintile_shift_summary,
          file.path(OUT_DIR, "09_final15_diagnostic_quintile_category_shift_summary.csv"),
          row.names = FALSE, na = "")

# -----------------------------------------------------------------------------
# 9. Item-level and group-level event summaries in frozen Full/Training/Test
# -----------------------------------------------------------------------------

message("[6/10] Summarising item and group status by depression events ...")

cohort_list <- list(Full = full_data, Training = train_data, Test = test_data)
item_event_rows <- list()
item_event_k <- 0L
group_event_rows <- list()
group_event_k <- 0L
group_sensitivity_rows <- list()
group_sensitivity_k <- 0L
quintile_sample_rows <- list()
quintile_sample_k <- 0L

for (sample_name in names(cohort_list)) {
  cohort <- cohort_list[[sample_name]]
  idx <- match(cohort$eid_audit, food_data$eid_audit)
  if (anyNA(idx)) stop_now("Cohort-to-food-data EID matching failed for ", sample_name)

  # Item-level events: the 85 attachment rows remain visible, including the
  # three white-bread rows that are not used in the final Step_02 score.
  for (i in seq_len(nrow(item_map))) {
    row <- item_map[i, ]
    status <- item_status_vectors[[as.character(row$item_order)]][idx]
    tab <- data.frame(
      status = status,
      depression_status = cohort$depression_status,
      stringsAsFactors = FALSE
    ) %>%
      group_by(status) %>%
      summarise(
        participants_n = n(),
        depression_events_n = sum(depression_status == 1, na.rm = TRUE),
        non_events_n = sum(depression_status == 0, na.rm = TRUE),
        outcome_missing_n = sum(is.na(depression_status)),
        event_rate_percent = calculate_event_rate(depression_events_n, non_events_n),
        .groups = "drop"
      )
    if (nrow(tab) > 0L) {
      tab$sample <- sample_name
      tab$item_order <- row$item_order
      tab$food_group <- row$food_group
      tab$food_item <- row$food_item
      tab$attachment_field_id <- row$attachment_field_id
      tab$effective_in_final_step02_score <- row$effective_in_final_step02_score
      item_event_k <- item_event_k + 1L
      item_event_rows[[item_event_k]] <- tab
    }
  }

  # Group-level score and missingness status.
  for (g in seq_len(nrow(food_specs))) {
    sp <- food_specs[g, ]
    original_score <- as_numeric_safe(food_data[[sp$intake_var]][idx])
    strict_score <- as_numeric_safe(group_score_na[[sp$prefix]][idx])
    missing_n <- group_missing_component_n[[paste0(sp$prefix, "_missing_item_n")]][idx]
    strict_status <- ifelse(
      is.na(strict_score),
      "NA_PRESERVED_DIAGNOSTIC_MISSING",
      "NA_PRESERVED_DIAGNOSTIC_EVALUABLE"
    )
    original_status <- ifelse(
      is.na(original_score), "ORIGINAL_SCORE_MISSING",
      ifelse(original_score == 0, "ORIGINAL_DERIVED_ZERO", "ORIGINAL_DERIVED_POSITIVE")
    )

    status_tab <- data.frame(
      strict_status = strict_status,
      depression_status = cohort$depression_status,
      stringsAsFactors = FALSE
    ) %>%
      group_by(strict_status) %>%
      summarise(
        participants_n = n(),
        depression_events_n = sum(depression_status == 1, na.rm = TRUE),
        non_events_n = sum(depression_status == 0, na.rm = TRUE),
        outcome_missing_n = sum(is.na(depression_status)),
        event_rate_percent = calculate_event_rate(depression_events_n, non_events_n),
        .groups = "drop"
      )
    status_tab$sample <- sample_name
    status_tab$food_group <- sp$food_group
    group_event_k <- group_event_k + 1L
    group_event_rows[[group_event_k]] <- status_tab

    evaluable <- !is.na(strict_score)
    group_sensitivity_k <- group_sensitivity_k + 1L
    group_sensitivity_rows[[group_sensitivity_k]] <- data.frame(
      sample = sample_name,
      food_group = sp$food_group,
      participants_n = length(idx),
      original_score_zero_n = sum(original_score == 0, na.rm = TRUE),
      original_score_zero_events_n = sum(cohort$depression_status[original_score == 0] == 1, na.rm = TRUE),
      original_score_zero_non_events_n = sum(cohort$depression_status[original_score == 0] == 0, na.rm = TRUE),
      original_score_zero_event_rate_percent = calculate_event_rate(
        sum(cohort$depression_status[original_score == 0] == 1, na.rm = TRUE),
        sum(cohort$depression_status[original_score == 0] == 0, na.rm = TRUE)
      ),
      original_score_positive_n = sum(original_score > 0, na.rm = TRUE),
      original_score_positive_events_n = sum(cohort$depression_status[original_score > 0] == 1, na.rm = TRUE),
      original_score_positive_non_events_n = sum(cohort$depression_status[original_score > 0] == 0, na.rm = TRUE),
      original_score_positive_event_rate_percent = calculate_event_rate(
        sum(cohort$depression_status[original_score > 0] == 1, na.rm = TRUE),
        sum(cohort$depression_status[original_score > 0] == 0, na.rm = TRUE)
      ),
      na_preserved_evaluable_n = sum(evaluable),
      na_preserved_missing_n = sum(!evaluable),
      na_preserved_evaluable_percent = 100 * mean(evaluable),
      na_preserved_events_n = sum(cohort$depression_status[evaluable] == 1, na.rm = TRUE),
      na_preserved_non_events_n = sum(cohort$depression_status[evaluable] == 0, na.rm = TRUE),
      na_preserved_event_rate_percent = calculate_event_rate(
        sum(cohort$depression_status[evaluable] == 1, na.rm = TRUE),
        sum(cohort$depression_status[evaluable] == 0, na.rm = TRUE)
      ),
      original_median = median(original_score, na.rm = TRUE),
      na_preserved_median = if (any(evaluable)) median(strict_score[evaluable], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )

    # Quintile/category shifts among the specified cohort participants.
    diag_q <- diagnostic_coding[[paste0(sp$prefix, "_diagnostic_quintile")]][idx]
    diag_cat <- diagnostic_coding[[paste0(sp$prefix, "_diagnostic_range_num")]][idx]
    orig_q <- as.character(food_data[[sp$quintile_var]][idx])
    orig_cat <- as.character(food_data[[sp$group_range_num_var]][idx])
    q_ok <- !is.na(diag_q) & !is.na(orig_q)
    c_ok <- !is.na(diag_cat) & !is.na(orig_cat)
    quintile_sample_k <- quintile_sample_k + 1L
    quintile_sample_rows[[quintile_sample_k]] <- data.frame(
      sample = sample_name,
      food_group = sp$food_group,
      participants_n = length(idx),
      diagnostic_evaluable_n = sum(!is.na(strict_score)),
      diagnostic_missing_n = sum(is.na(strict_score)),
      quintile_comparable_n = sum(q_ok),
      quintile_changed_n = sum(orig_q[q_ok] != diag_q[q_ok]),
      quintile_changed_percent = if (sum(q_ok) > 0L) 100 * mean(orig_q[q_ok] != diag_q[q_ok]) else NA_real_,
      category_comparable_n = sum(c_ok),
      category_changed_n = sum(orig_cat[c_ok] != as.character(diag_cat[c_ok])),
      category_changed_percent = if (sum(c_ok) > 0L) 100 * mean(orig_cat[c_ok] != as.character(diag_cat[c_ok])) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

item_event_summary <- bind_rows(item_event_rows) %>%
  arrange(match(sample, c("Full", "Training", "Test")), item_order, status)
group_event_summary <- bind_rows(group_event_rows) %>%
  arrange(match(sample, c("Full", "Training", "Test")), match(food_group, food_specs$food_group), strict_status)
group_sensitivity_summary <- bind_rows(group_sensitivity_rows) %>%
  arrange(match(sample, c("Full", "Training", "Test")), match(food_group, food_specs$food_group))
quintile_sample_summary <- bind_rows(quintile_sample_rows) %>%
  arrange(match(sample, c("Full", "Training", "Test")), match(food_group, food_specs$food_group))

# -----------------------------------------------------------------------------
# 10. Joint 15-component complete diagnostic cohort and final QC
# -----------------------------------------------------------------------------

message("[7/10] Summarising joint 15-component completeness and QC ...")
strict_matrix_all <- as.matrix(group_score_na[, EXPECTED_FINAL15, drop = FALSE])
joint_complete <- rowSums(is.na(strict_matrix_all)) == 0L
participant_joint <- data.frame(
  eid_audit = food_data$eid_audit,
  strict_15_component_complete = joint_complete,
  groups_with_na_preserved_missing_n = rowSums(is.na(strict_matrix_all)),
  stringsAsFactors = FALSE
)

joint_rows <- list()
for (sample_name in names(cohort_list)) {
  cohort <- cohort_list[[sample_name]]
  idx <- match(cohort$eid_audit, participant_joint$eid_audit)
  ok <- participant_joint$strict_15_component_complete[idx]
  joint_rows[[length(joint_rows) + 1L]] <- data.frame(
    sample = sample_name,
    participants_n = nrow(cohort),
    strict_15_component_complete_n = sum(ok),
    strict_15_component_incomplete_n = sum(!ok),
    strict_15_component_complete_events_n = sum(cohort$depression_status[ok] == 1, na.rm = TRUE),
    strict_15_component_complete_non_events_n = sum(cohort$depression_status[ok] == 0, na.rm = TRUE),
    strict_15_component_complete_event_rate_percent = calculate_event_rate(
      sum(cohort$depression_status[ok] == 1, na.rm = TRUE),
      sum(cohort$depression_status[ok] == 0, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}
joint_complete_summary <- bind_rows(joint_rows)

participant_joint_full <- participant_joint %>%
  inner_join(full_data %>% select(eid_audit, depression_status), by = "eid_audit")

# Full eligible participant-level count only; no source EID list is exported.
participant_joint_count_summary <- participant_joint_full %>%
  count(groups_with_na_preserved_missing_n, strict_15_component_complete,
        name = "participants_n") %>%
  arrange(groups_with_na_preserved_missing_n)

qc_summary <- data.frame(
  check = c(
    "Original attachment total row count",
    "Final 15 attachment item rows",
    "Final 15 food-group count and order",
    "Effective Step_02_1_1 components (81 direct + final Bread roll)",
    "All attachment portion values match first distinct project portion metadata",
    "All direct item descriptions/groups match original dietary dictionary",
    "All source field encoding IDs are available",
    "All raw i0-i4 amount columns are available",
    "All original white-bread selector columns are available",
    "All 15 reconstructed group totals match saved Step_02_1_1 totals",
    "Frozen MAIN cohort/event counts match"
  ),
  result = c(
    ifelse(ATTACHMENT_TOTAL_ROWS == 206L, "PASS", "FAIL"),
    ifelse(nrow(item_map) == 85L, "PASS", "FAIL"),
    ifelse(identical(food_specs$prefix, EXPECTED_FINAL15), "PASS", "FAIL"),
    ifelse(nrow(effective_item_map) == 82L, "PASS", "FAIL"),
    ifelse(all(item_map$portion_metadata_match & item_map$portion_duplicate_values_consistent), "PASS", "FAIL"),
    ifelse(nrow(bad_direct_dictionary) == 0L, "PASS", "FAIL"),
    ifelse(all(item_map$field_encoding_present), "PASS", "FAIL"),
    ifelse(all(item_map$raw_columns_n == 5L), "PASS", "FAIL"),
    ifelse(nrow(selector_presence_summary) == 4L && all(selector_presence_summary$selector_columns_n == 5L), "PASS", "FAIL"),
    ifelse(all(group_reconstruction_summary$reconstruction_status == "PASS"), "PASS", "FAIL"),
    ifelse(all(cohort_check$counts_match), "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 11. Export audit tables and a concise run report
# -----------------------------------------------------------------------------

message("[8/10] Exporting Step_03_2_7 outputs ...")

write.csv(item_map,
          file.path(OUT_DIR, "01_final15_attachment_85_item_mapping_frozen.csv"),
          row.names = FALSE, na = "")
write.csv(mapping_validation,
          file.path(OUT_DIR, "02_attachment_85_item_mapping_and_metadata_validation.csv"),
          row.names = FALSE, na = "")
write.csv(item_level_summary,
          file.path(OUT_DIR, "03_final15_item_level_raw_response_summary.csv"),
          row.names = FALSE, na = "")
write.csv(selector_presence_summary,
          file.path(OUT_DIR, "03a_original_step02_bread_selector_dependency_audit.csv"),
          row.names = FALSE, na = "")
write.csv(pseudo_slot_audit,
          file.path(OUT_DIR, "04_refined_grains_final_step02_pseudo_slot_audit.csv"),
          row.names = FALSE, na = "")
write.csv(item_event_summary,
          file.path(OUT_DIR, "05_final15_item_level_response_status_event_summary.csv"),
          row.names = FALSE, na = "")
write.csv(group_reconstruction_summary,
          file.path(OUT_DIR, "06_final15_group_reconstruction_and_na_preserved_summary.csv"),
          row.names = FALSE, na = "")
write.csv(group_event_summary,
          file.path(OUT_DIR, "07_final15_group_na_preserved_status_event_summary.csv"),
          row.names = FALSE, na = "")
write.csv(group_sensitivity_summary,
          file.path(OUT_DIR, "08_final15_group_original_vs_na_preserved_sensitivity.csv"),
          row.names = FALSE, na = "")
write.csv(quintile_shift_summary,
          file.path(OUT_DIR, "09_final15_diagnostic_quintile_category_shift_full_food_data.csv"),
          row.names = FALSE, na = "")
write.csv(quintile_sample_summary,
          file.path(OUT_DIR, "10_final15_diagnostic_quintile_category_shift_by_cohort.csv"),
          row.names = FALSE, na = "")
write.csv(joint_complete_summary,
          file.path(OUT_DIR, "11_joint_15_component_na_preserved_cohort_event_summary.csv"),
          row.names = FALSE, na = "")
write.csv(participant_joint_count_summary,
          file.path(OUT_DIR, "12_joint_15_component_missing_count_distribution.csv"),
          row.names = FALSE, na = "")
write.csv(cohort_event_summary,
          file.path(OUT_DIR, "13_frozen_cohort_event_count_reconciliation.csv"),
          row.names = FALSE, na = "")
write.csv(qc_summary,
          file.path(OUT_DIR, "14_Step_03_2_7_qc_summary.csv"),
          row.names = FALSE, na = "")

report_lines <- c(
  "Step_03_2_7 final15 original-table item-level and portion audit V1.0",
  "====================================================================",
  paste0("Original attachment rows: ", ATTACHMENT_TOTAL_ROWS),
  paste0("Final-15 mapping rows: ", nrow(item_map)),
  paste0("Direct original Field ID rows: ", sum(item_map$item_type == "original_field")),
  paste0("White-bread pseudo rows in attachment: ", sum(item_map$item_type == "white_bread_pseudo")),
  paste0("Effective Step_02_1_1 score components: ", nrow(effective_item_map)),
  paste0("Raw dietary participants: ", nrow(raw_data)),
  paste0("Derived food-group participants: ", nrow(food_data)),
  paste0("Training N/events: ", nrow(train_data), "/", sum(train_data$depression_status == 1, na.rm = TRUE)),
  paste0("Test N/events: ", nrow(test_data), "/", sum(test_data$depression_status == 1, na.rm = TRUE)),
  paste0("Full eligible N/events: ", nrow(full_data), "/", sum(full_data$depression_status == 1, na.rm = TRUE)),
  paste0("15-group exact reconstruction mismatches: ", sum(group_reconstruction_summary$reconstruction_mismatch_n)),
  paste0("Joint 15-component NA-preserved complete N in Full: ", joint_complete_summary$strict_15_component_complete_n[joint_complete_summary$sample == "Full"]),
  "",
  "Interpretation constraints:",
  "- The original Step_02_1_1 score is reconstructed; no source score is modified.",
  "- Portion grams are used only for the original encoding_id==0 division rule; no frequency*portion multiplication is performed.",
  "- The NA-preserved scenario is diagnostic only and does not assert that all missing cells represent non-intake.",
  "- The final Step_02_1_1 refined-grains pseudo mapping retains only the Bread roll source p101160 selected by p20094 in the overwritten pw100950 slot.",
  "- Three other white-bread pseudo rows are retained in the attachment mapping audit but do not contribute to the final score after the sequential overwrite.",
  "- No Gate fields, p20082, extra food groups, beta-weight recalculation, or Cox refitting are used."
)
writeLines(report_lines, con = file.path(OUT_DIR, "15_Step_03_2_7_run_report.txt"), useBytes = TRUE)

# Console summaries use write.table to avoid print.default/na.print issues.
cat("\n============================================================\n")
cat("Step_03_2_7 final15 original-table item-level and portion audit completed\n")
cat("============================================================\n")
cat("Attachment rows: ", ATTACHMENT_TOTAL_ROWS, "\n", sep = "")
cat("Final-15 table rows: ", nrow(item_map), "\n", sep = "")
cat("Effective Step_02_1_1 components: ", nrow(effective_item_map), "\n", sep = "")
cat("Raw dietary N: ", nrow(raw_data), "\n", sep = "")
cat("Full eligible N/events: ", nrow(full_data), "/",
    sum(full_data$depression_status == 1, na.rm = TRUE), "\n", sep = "")
cat("Reconstruction mismatches: ",
    sum(group_reconstruction_summary$reconstruction_mismatch_n), "\n", sep = "")
cat("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = "")
cat("\nQC summary:\n")
write_tsv_console(qc_summary)
cat("\nGroup reconstruction summary:\n")
write_tsv_console(group_reconstruction_summary %>%
                    select(food_group, effective_step02_components_n,
                           reconstruction_mismatch_n,
                           na_preserved_diagnostic_evaluable_n,
                           na_preserved_diagnostic_missing_n,
                           reconstruction_status))
cat("\nJoint 15-component diagnostic cohort:\n")
write_tsv_console(joint_complete_summary)
cat("\nREAD-ONLY AUDIT: original scores, quintiles, categories, cohorts, beta weights, and Cox models were not modified.\n")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 04/11: Step_03_2_7_FINAL15_AUDIT
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 05/11: Step_04_1_1_1_RCS_BASE
# Source: Step_04_1_1_1_rcs_legacy_anchored.R
# Provenance: project source archive data_analysis_v2.zip
# SHA256 of included block: 74a058d5864f061827353bee3daf4166f6de0dce0979bb48c074f734011b49e6
########################################################################################

# =============================================================================
# Step_04_1_1_RCS_LEGACY_ANCHORED_LOWEST_RISK_REFERENCE_USER_RANGES.R
#
# Purpose
#   applies the user's supplied 15-component scoring-interval table as
#   explicit primary scoring anchors. The V5 Cox/RCS workflow is retained to
#   provide the same model diagnostics, lowest-risk reference, P10/P50/P90,
#   knot information, legacy-reference comparison, and component-level audits.
#   Importantly, the supplied anchor boundaries are NOT silently moved to
#   obtain statistical significance. Cox/RCS findings are reported as
#   validation/context; they do not redefine the prespecified boundaries.
#
# Key principles
#   1. Use the exact frozen Step_03_1_0 70% training cohort; never re-split.
#   2. Keep the exact 15 components and existing Model 3 covariates.
#   3. Main Figure S4 RCS reference = median of the Step_03_1_0 lowest-risk
#      category. The HR at this reference is exactly 1.
#   4. Fit the original minimum-intake-reference model only for legacy
#      comparison; it does not determine the primary scoring boundaries.
#   5. Apply the exact value1/value2 anchors supplied by the user. Finite
#      value2 means an inclusive middle interval; NA uses the explicit
#      direction map documented below (low or high one-sided interval).
#   6. Export the previous V5 data-derived interval alongside the anchor
#      so boundary differences and direction changes can be audited.
#   7. RCS 5%/10%/20% minimum-risk plateaus and CI-defined legacy intervals are
#      diagnostics; they are not used to move the anchor boundaries.
#
# Frozen cohort and model framework
#   - Exact Step_03_1_0 frozen 70% training RDS and participant-ID manifest;
#   - 15 current food components and the established Model 3 covariates;
#   - training-only RCS and scoring; the 30% test cohort is not accessed;
#   - P10/P50/P90 annotations and three-knot RCS diagnostics retained.
#
# output
#   All results are written to:
#     ./results/Step_04/Step_04_1_1_1
# =============================================================================

# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(rms)
  library(cols4all)
  library(cowplot)
})

if (!requireNamespace("officer", quietly = TRUE) ||
    !requireNamespace("flextable", quietly = TRUE)) {
  stop(
    "STOP: V5 Word-table export requires the R packages 'officer' and ",
    "'flextable'. Install them with install.packages(c('officer','flextable')) ",
    "and rerun this script."
  )
}

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 0. Frozen paths and analysis settings
# -----------------------------------------------------------------------------

SPLIT_DIR <- "./results/Step_03/ukb_70_30_split"
MAIN_SPLIT_DIR <- file.path(SPLIT_DIR, "main")
SPLIT_MANIFEST_FILE <- file.path(
  SPLIT_DIR, "ukb_70_30_split_manifest.csv"
)
TRAIN_RDS <- file.path(
  MAIN_SPLIT_DIR, "ukb_depression_training_cohort.rds"
)
REF_CSV <- paste0(
  "./results/Step_03/reference_rotation_lowest_risk/",
  "food_group_model2_results_sig.csv"
)

OUT_DIR <- "./results/Step_04/Step_04_1_1_1"
PRED_DIR <- file.path(OUT_DIR, "predictions")
PLOT_DIR <- file.path(OUT_DIR, "plots")
AUDIT_DIR <- file.path(OUT_DIR, "audit")

dir.create(PRED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPECTED_SPLIT_SEED <- 20261103L
EXPECTED_TRAIN_PROP <- 0.70

PREDICTION_GRID_N <- 501L
PREDICTION_SUPPORT_PROBS <- c(0.10, 0.90)
RANGE_TOLERANCES <- c(0.05, 0.10, 0.20)

# The original Model 3 adjustment set is retained without modification.
MODEL3_COVARIATES <- c(
  "age_food",
  "sex",
  "ethnicity",
  "education_attainment",
  "twnsend",
  "career",
  "smoking_status",
  "drink_status",
  "bmi_group",
  "physical_activity_proup",
  "sedentary_group",
  "household_income",
  "energy_intake",
  "central_obesity",
  "triglycerides",
  "hdl_cholesterol",
  "t2diabetes_status_final",
  "hypertension_status_final",
  "dyslipidemia_status_final",
  "cvd_status_final",
  "all_cancer_status_final"
)

# -----------------------------------------------------------------------------
# 1. Read the exact Step_03_1_0 training cohort and validate its manifest
# -----------------------------------------------------------------------------

if (!file.exists(SPLIT_MANIFEST_FILE)) {
  stop(
    "STOP: Frozen Step_03_1_0 split manifest is missing: ",
    SPLIT_MANIFEST_FILE,
    ". Run Step_03_1_0 first. This script will not create a new split."
  )
}
if (!file.exists(TRAIN_RDS)) {
  stop(
    "STOP: Step_03_1_0 training RDS is missing: ", TRAIN_RDS,
    ". Run Step_03_1_0 first. This script will not reconstruct or resample it."
  )
}
if (!file.exists(REF_CSV)) {
  stop(
    "STOP: Step_03_1_0 Model 2 significant-results CSV is missing: ",
    REF_CSV
  )
}

split_manifest <- read.csv(
  SPLIT_MANIFEST_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(
    participant_id = "character",
    split = "character",
    id_column = "character",
    split_seed = "integer",
    training_fraction = "numeric"
  )
)

required_manifest_cols <- c(
  "participant_id", "split", "id_column",
  "split_seed", "training_fraction"
)
if (!all(required_manifest_cols %in% names(split_manifest))) {
  stop("STOP: Unexpected Step_03_1_0 split-manifest structure.")
}
if (anyNA(split_manifest$participant_id) ||
    any(trimws(split_manifest$participant_id) == "") ||
    anyDuplicated(split_manifest$participant_id) > 0L) {
  stop("STOP: Invalid or duplicated participant IDs in the split manifest.")
}
if (!all(split_manifest$split %in% c("training", "test"))) {
  stop("STOP: Invalid split labels in the frozen manifest.")
}
if (length(unique(split_manifest$id_column)) != 1L ||
    anyNA(split_manifest$id_column)) {
  stop("STOP: The manifest must contain exactly one non-missing ID column.")
}
if (length(unique(split_manifest$split_seed)) != 1L ||
    anyNA(split_manifest$split_seed) ||
    unique(split_manifest$split_seed) != EXPECTED_SPLIT_SEED) {
  stop("STOP: Split seed does not match Step_03_1_0 (20261103).")
}
if (length(unique(split_manifest$training_fraction)) != 1L ||
    anyNA(split_manifest$training_fraction) ||
    abs(unique(split_manifest$training_fraction) - EXPECTED_TRAIN_PROP) > 1e-12) {
  stop("STOP: Training fraction does not match Step_03_1_0 (70%).")
}

ID_COL <- unique(split_manifest$id_column)
training_ids <- sort(
  split_manifest$participant_id[split_manifest$split == "training"],
  method = "radix"
)
test_ids <- sort(
  split_manifest$participant_id[split_manifest$split == "test"],
  method = "radix"
)

ukb_data_depression_rcs <- readRDS(TRAIN_RDS)

if (!is.data.frame(ukb_data_depression_rcs)) {
  stop("STOP: The Step_03_1_0 training RDS is not a data.frame.")
}
if (!ID_COL %in% names(ukb_data_depression_rcs)) {
  stop("STOP: ID column from manifest is absent in the training RDS: ", ID_COL)
}

training_ids_observed <- sort(
  trimws(as.character(ukb_data_depression_rcs[[ID_COL]])),
  method = "radix"
)
if (anyNA(training_ids_observed) ||
    any(training_ids_observed == "") ||
    anyDuplicated(training_ids_observed) > 0L) {
  stop("STOP: Invalid or duplicated participant IDs in the training RDS.")
}
if (!identical(training_ids_observed, training_ids)) {
  stop(
    "STOP: Training RDS IDs do not exactly match the Step_03_1_0 manifest. ",
    "No RCS models have been fitted."
  )
}
if (length(intersect(training_ids_observed, test_ids)) != 0L) {
  stop("STOP: Test participants are present in the RCS training RDS.")
}
if (nrow(ukb_data_depression_rcs) != length(training_ids)) {
  stop("STOP: Training RDS row count does not match manifest training IDs.")
}
if ("analysis_split" %in% names(ukb_data_depression_rcs) &&
    any(as.character(ukb_data_depression_rcs$analysis_split) != "training")) {
  stop("STOP: Non-training rows were found in the Step_03_1_0 training RDS.")
}
if (!all(c("depression_time", "depression_status") %in%
         names(ukb_data_depression_rcs))) {
  stop("STOP: Depression survival variables are missing from training RDS.")
}

TRAINING_QC <- tibble(
  split_source = "Step_03_1_0 frozen training RDS",
  manifest = SPLIT_MANIFEST_FILE,
  training_rds = TRAIN_RDS,
  id_column = ID_COL,
  split_seed = EXPECTED_SPLIT_SEED,
  training_fraction = EXPECTED_TRAIN_PROP,
  manifest_eligible_n = nrow(split_manifest),
  manifest_training_n = length(training_ids),
  manifest_test_n = length(test_ids),
  training_rds_n = nrow(ukb_data_depression_rcs),
  training_depression_events = sum(
    ukb_data_depression_rcs$depression_status == 1, na.rm = TRUE
  ),
  training_ids_exact_match = TRUE,
  test_ids_excluded = TRUE
)

write.csv(
  TRAINING_QC,
  file.path(AUDIT_DIR, "RCS_training_cohort_QC.csv"),
  row.names = FALSE
)

cat("\n===== FROZEN TRAINING COHORT QC =====\n")
print(TRAINING_QC)

# -----------------------------------------------------------------------------
# 2. Exact 15 components from Step_03_1_0 Model 2 significant results
# -----------------------------------------------------------------------------

FOOD_MAP <- tibble(
  group = c(
    "Citrus fruits",
    "Red/orange vegetables",
    "Other vegetables",
    "Potatoes",
    "Refined grains",
    "Eggs",
    "Low fat dairy",
    "Red meats",
    "Processed meats",
    "Fish and seafood",
    "Coffee",
    "Fruit juices",
    "Sweetened beverages",
    "Water",
    "Wine"
  ),
  food = c(
    "citrus_fruits",
    "red_orange_veg",
    "other_veg",
    "potatoes",
    "refined_grains",
    "eggs",
    "low_fat_dairy",
    "red_meats",
    "processed_meats",
    "fish_seafood",
    "coffee",
    "fruit_juices",
    "sweetened_bev",
    "water",
    "wine"
  ),
  category_var = c(
    "citrus_fruits_group_range",
    "red_orange_veg_group_range",
    "other_veg_group_range",
    "potatoes_group_range",
    "refined_grains_group_range",
    "eggs_group_range",
    "low_fat_dairy_group_range",
    "red_meats_group_range",
    "processed_meats_group_range",
    "fish_seafood_group_range",
    "coffee_group_range",
    "fruit_juices_group_range",
    "sweetened_bev_group_range",
    "water_group_range",
    "wine_group_range"
  ),
  display_name = c(
    "Citrus fruits",
    "Red/Orange vegetables",
    "Other vegetables",
    "Potatoes",
    "Refined grains",
    "Eggs",
    "Low-fat dairy",
    "Red meats",
    "Processed meats",
    "Fish and seafood",
    "Coffee",
    "Fruit juices",
    "Sweetened beverages",
    "Water",
    "Wine"
  ),
  x_break = c(1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  x_limit_fixed = c(7, 9, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  annotation_x = c(1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

if (nrow(FOOD_MAP) != 15L ||
    anyDuplicated(FOOD_MAP$food) > 0L ||
    anyDuplicated(FOOD_MAP$group) > 0L ||
    anyDuplicated(FOOD_MAP$category_var) > 0L) {
  stop("STOP: Internal 15-component map is invalid.")
}

if (any(c("green_leafy_veg", "high_fat_dairy") %in% FOOD_MAP$food)) {
  stop("STOP: Excluded components entered the 15-component RCS map.")
}
if (!all(c("fruit_juices", "water") %in% FOOD_MAP$food)) {
  stop("STOP: Fruit juices and water must be included in the corrected 15.")
}


# -----------------------------------------------------------------------------
# User-specified primary scoring anchors (servings/day).
#
# Direction is explicit because value2 = NA alone cannot distinguish a low-
# intake cutoff (x <= value1) from a high-intake open tail (x >= value1).
# The directions below preserve the prior high-intake direction where the
# previous range audit supported it. Zero cutoffs for potatoes/red meat are
# interpreted as low-intake/no-intake rules to avoid an all-one score.
# Low-fat dairy retains the original one-sided low-intake rule. Processed meats
# uses the one-sided low-intake cutoff x <= 2.7 as specified in this revision.
# If the researcher's intended direction differs, change only this table.
# -----------------------------------------------------------------------------
V51_RANGE_ANCHORS <- tibble(
  food = c(
    "citrus_fruits", "fruit_juices", "red_orange_veg", "other_veg",
    "potatoes", "refined_grains", "eggs", "water", "low_fat_dairy",
    "red_meats", "processed_meats", "fish_seafood", "coffee",
    "sweetened_bev", "wine"
  ),
  value1 = c(
    0.034338, 0, 0.01, 0.01, 0, 3, 0.030151, 0.01, 3.071574,
    0, 2.7, 0, 0.01, 0, 0
  ),
  value2 = c(
    1.682583, 1.37, NA_real_, 4.62, NA_real_, NA_real_, 1.869347,
    3.08, NA_real_, NA_real_, NA_real_, 1.98, 4.5, 2.08, 2
  ),
  range_direction = c(
    "middle", "middle", "high", "middle", "low", "high", "middle",
    "middle", "low", "low", "v5_derived", "middle", "middle", "middle",
    "middle"
  ),
  direction_basis = c(
    "User-specified bounded interval",
    "User-specified bounded interval",
    "Preserved high-intake direction from prior component-range specification",
    "User-specified bounded interval",
    "Zero-intake cutoff interpreted as low-side rule (x <= 0); high-side would score all nonnegative intakes as 1",
    "Preserved high-intake direction from prior component-range specification",
    "User-specified bounded interval",
    "User-specified bounded interval",
    "Original one-sided low-intake scoring direction retained",
    "Zero-intake cutoff interpreted as low-side rule (x <= 0); high-side would score all nonnegative intakes as 1",
    "Use the V5-derived one-sided direction (low/high) for the 2.7 cutoff; do not infer direction from NA alone",
    "User-specified bounded interval",
    "User-specified bounded interval",
    "User-specified bounded interval (replaces prior one-sided sweetened-beverage cutoff)",
    "User-specified bounded interval"
  )
)

if (nrow(V51_RANGE_ANCHORS) != 15L ||
    anyDuplicated(V51_RANGE_ANCHORS$food) > 0L ||
    !setequal(V51_RANGE_ANCHORS$food, FOOD_MAP$food)) {
  stop("STOP: anchor table must contain exactly the frozen 15 foods.")
}
V51_RANGE_ANCHORS <- V51_RANGE_ANCHORS[
  match(FOOD_MAP$food, V51_RANGE_ANCHORS$food), , drop = FALSE
]
if (!identical(V51_RANGE_ANCHORS$food, FOOD_MAP$food)) {
  stop("STOP: anchor table could not be aligned to FOOD_MAP.")
}
if (any(!is.finite(V51_RANGE_ANCHORS$value1)) ||
    any(!V51_RANGE_ANCHORS$range_direction %in% c("middle", "low", "high", "v5_derived"))) {
  stop("STOP: anchor values/directions are invalid.")
}
if (any(V51_RANGE_ANCHORS$range_direction == "middle" &
        (!is.finite(V51_RANGE_ANCHORS$value2) |
         V51_RANGE_ANCHORS$value2 < V51_RANGE_ANCHORS$value1))) {
  stop("STOP: Every middle interval must have finite ordered bounds.")
}
if (any(V51_RANGE_ANCHORS$range_direction %in% c("low", "high", "v5_derived") &
        !is.na(V51_RANGE_ANCHORS$value2))) {
  stop("STOP: one-sided intervals must have value2 = NA.")
}
write.csv(
  V51_RANGE_ANCHORS,
  file.path(AUDIT_DIR, "user_prespecified_range_anchors.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

model2_sig <- read.csv(
  REF_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_result_cols <- c("group", "model", "variable", "hr_ci")
if (!all(required_result_cols %in% names(model2_sig))) {
  stop(
    "STOP: Model 2 results CSV is missing columns: ",
    paste(setdiff(required_result_cols, names(model2_sig)), collapse = ", ")
  )
}

model2_sig <- model2_sig %>%
  tidyr::fill(group, .direction = "down") %>%
  filter(model == "model2")

csv_groups <- unique(model2_sig$group)
if (!setequal(csv_groups, FOOD_MAP$group)) {
  stop(
    "STOP: The Model 2 significant-results CSV does not contain the exact ",
    "15-component set. Missing: ",
    paste(setdiff(FOOD_MAP$group, csv_groups), collapse = "; "),
    ". Extra: ",
    paste(setdiff(csv_groups, FOOD_MAP$group), collapse = "; ")
  )
}

reference_rows <- model2_sig %>%
  filter(group %in% FOOD_MAP$group, hr_ci == "Ref.") %>%
  select(group, variable, hr_ci)

reference_counts <- reference_rows %>%
  count(group, name = "n_reference_rows")

if (nrow(reference_counts) != 15L ||
    any(reference_counts$n_reference_rows != 1L)) {
  stop(
    "STOP: Each of the 15 components must have exactly one Model 2 Ref. row."
  )
}

REFERENCE_SOURCE <- FOOD_MAP %>%
  left_join(reference_rows, by = "group") %>%
  rename(reference_category_label = variable)

if (nrow(REFERENCE_SOURCE) != 15L ||
    anyNA(REFERENCE_SOURCE$reference_category_label)) {
  stop("STOP: Could not map the 15 lowest-risk reference categories.")
}

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

format_p_value <- function(p) {
  if (length(p) != 1L || is.na(p) || !is.finite(p)) return("= NA")
  if (p < 0.001) return("< 0.001")
  paste0("= ", sprintf("%.3f", p))
}

fmt_num <- function(x) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
  if (abs(x) < 1e-10) return("0")
  formatC(x, format = "f", digits = 2)
}

parse_leading_numeric <- function(x) {
  pattern <- "^[[:space:]]*[-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?"
  out <- regmatches(
    as.character(x),
    regexpr(pattern, as.character(x), perl = TRUE)
  )
  out[out == ""] <- NA_character_
  suppressWarnings(as.numeric(out))
}

# Choose three distinct knots. The primary rule is full-sample P10/P50/P90.
# If ties caused by a large zero-intake mass collapse those quantiles, use
# zero plus the positive-consumer P50 and P90. If this is still impossible,
# use three distinct observed support values and record that fallback.
choose_three_knots <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 100L || length(unique(x)) < 3L) {
    stop("STOP: Fewer than three distinct exposure values.")
  }

  q_full <- as.numeric(stats::quantile(
    x, probs = c(0.10, 0.50, 0.90),
    na.rm = TRUE, names = FALSE, type = 7
  ))
  names(q_full) <- c("P10", "P50", "P90")

  if (length(unique(q_full)) == 3L && all(diff(q_full) > 0)) {
    return(list(
      knots = q_full,
      method = "full_training_sample_P10_P50_P90",
      positive_quantiles = c(NA_real_, NA_real_, NA_real_)
    ))
  }

  zero_prop <- mean(x <= 1e-12)
  x_positive <- x[x > 1e-12]

  if (zero_prop > 0 && length(unique(x_positive)) >= 2L) {
    q_positive <- as.numeric(stats::quantile(
      x_positive, probs = c(0.10, 0.50, 0.90),
      na.rm = TRUE, names = FALSE, type = 7
    ))
    names(q_positive) <- c("P10_positive", "P50_positive", "P90_positive")

    zero_anchored <- c(0, q_positive["P50_positive"], q_positive["P90_positive"])
    if (length(unique(zero_anchored)) == 3L &&
        all(is.finite(zero_anchored)) &&
        all(diff(zero_anchored) > 0)) {
      return(list(
        knots = unname(zero_anchored),
        method = "zero_anchored_plus_positive_consumer_P50_P90",
        positive_quantiles = q_positive
      ))
    }
  }

  # Last-resort fallback for heavily tied/discrete non-zero distributions.
  # This is not silently described as P10/P50/P90.
  ux <- sort(unique(x))
  if (length(ux) >= 3L) {
    idx <- unique(pmax(
      1L,
      pmin(length(ux), as.integer(round(
        c(0.10, 0.50, 0.90) * (length(ux) - 1L)
      ) + 1L))
    ))
    if (length(idx) == 3L) {
      fallback_knots <- ux[idx]
      if (all(diff(fallback_knots) > 0)) {
        return(list(
          knots = fallback_knots,
          method = "distinct_observed_support_fallback",
          positive_quantiles = c(NA_real_, NA_real_, NA_real_)
        ))
      }
    }
  }

  stop(
    "STOP: Three distinct RCS knots cannot be defined without arbitrary ",
    "extrapolation. Inspect the exposure distribution."
  )
}

linear_crossing_x <- function(x1, y1, x2, y2, y_target) {
  if (any(!is.finite(c(x1, y1, x2, y2, y_target)))) return(NA_real_)
  if (abs(y2 - y1) < .Machine$double.eps^0.5) return(mean(c(x1, x2)))
  x1 + (y_target - y1) * (x2 - x1) / (y2 - y1)
}

get_minimum_risk_plateau <- function(x, y, tolerance) {
  ok <- is.finite(x) & is.finite(y) & y > 0
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 5L) {
    return(c(lower = NA_real_, upper = NA_real_, threshold = NA_real_))
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  i_min <- which.min(y)
  y_min <- y[i_min]
  threshold <- y_min * (1 + tolerance)

  i_left <- i_min
  while (i_left > 1L && y[i_left - 1L] <= threshold) {
    i_left <- i_left - 1L
  }

  i_right <- i_min
  while (i_right < length(y) && y[i_right + 1L] <= threshold) {
    i_right <- i_right + 1L
  }

  lower <- if (i_left == 1L) {
    x[1]
  } else {
    linear_crossing_x(
      x[i_left - 1L], y[i_left - 1L],
      x[i_left], y[i_left], threshold
    )
  }

  upper <- if (i_right == length(y)) {
    x[length(x)]
  } else {
    linear_crossing_x(
      x[i_right], y[i_right],
      x[i_right + 1L], y[i_right + 1L],
      threshold
    )
  }

  c(lower = lower, upper = upper, threshold = threshold)
}

classify_rcs_shape <- function(pred, food_var, p_overall, p_non_linear) {
  pred <- pred %>%
    filter(
      is.finite(.data[[food_var]]),
      is.finite(yhat),
      yhat > 0
    ) %>%
    arrange(.data[[food_var]])

  x <- pred[[food_var]]
  y <- pred$yhat

  if (length(y) < 5L) {
    return(list(
      curve_shape = "insufficient_prediction_points",
      shape_evidence = "not_assessable",
      x_min = NA_real_,
      hr_min = NA_real_,
      left_tail_ratio = NA_real_,
      right_tail_ratio = NA_real_,
      range_eligible = FALSE
    ))
  }

  i_min <- which.min(y)
  i_max <- which.max(y)
  x_min <- x[i_min]
  hr_min <- y[i_min]

  left_ratio <- y[1] / hr_min
  right_ratio <- y[length(y)] / hr_min
  rise_cut <- 1.02

  if (i_min > 1L && i_min < length(y) &&
      left_ratio > rise_cut && right_ratio > rise_cut) {
    curve_shape <- "U-shaped"
  } else if (i_min == 1L && right_ratio > rise_cut) {
    curve_shape <- "J-shaped_increasing-risk"
  } else if (i_min == length(y) && left_ratio > rise_cut) {
    curve_shape <- "reverse-J-shaped_decreasing-risk"
  } else if (i_max > 1L && i_max < length(y) &&
             y[1] < y[i_max] / rise_cut &&
             y[length(y)] < y[i_max] / rise_cut) {
    curve_shape <- "inverted-U-shaped"
  } else if (y[length(y)] > y[1] * rise_cut) {
    curve_shape <- "monotonic-increasing-risk"
  } else if (y[1] > y[length(y)] * rise_cut) {
    curve_shape <- "monotonic-decreasing-risk"
  } else {
    curve_shape <- "flat-or-complex"
  }

  nonlinear_supported <- is.finite(p_non_linear) && p_non_linear < 0.05
  overall_supported <- is.finite(p_overall) && p_overall < 0.05

  shape_supported <- if (curve_shape == "U-shaped") {
    nonlinear_supported
  } else if (curve_shape %in% c(
    "J-shaped_increasing-risk",
    "reverse-J-shaped_decreasing-risk"
  )) {
    nonlinear_supported || overall_supported
  } else if (curve_shape %in% c(
    "monotonic-increasing-risk",
    "monotonic-decreasing-risk"
  )) {
    overall_supported
  } else {
    FALSE
  }

  shape_evidence <- if (shape_supported) {
    paste0(curve_shape, "_with_statistical_support")
  } else if (curve_shape %in% c(
    "U-shaped",
    "J-shaped_increasing-risk",
    "reverse-J-shaped_decreasing-risk",
    "monotonic-increasing-risk",
    "monotonic-decreasing-risk"
  )) {
    paste0(curve_shape, "_appearance_without_sufficient_statistical_support")
  } else {
    "no_clear_single_optimum_shape"
  }

  list(
    curve_shape = curve_shape,
    shape_evidence = shape_evidence,
    x_min = x_min,
    hr_min = hr_min,
    left_tail_ratio = left_ratio,
    right_tail_ratio = right_ratio,
    range_eligible = shape_supported
  )
}

# Category direction is determined from the observed ordering of the intake
# categories in the Step_03_1_0 training data, not from alphabetical labels.
get_category_interval <- function(data, food_var, category_var, ref_label) {
  x <- data[[food_var]]
  cat <- as.character(data[[category_var]])

  group_stats <- tibble(
    category = cat,
    intake = x
  ) %>%
    filter(!is.na(category), is.finite(intake)) %>%
    group_by(category) %>%
    summarise(
      category_n = n(),
      category_median = median(intake),
      category_min = min(intake),
      category_max = max(intake),
      .groups = "drop"
    )

  # Calculate category-specific event counts using explicit row masks.
  group_stats$category_events <- vapply(
    group_stats$category,
    function(g) {
      sum(data$depression_status[
        !is.na(cat) & cat == g
      ] == 1, na.rm = TRUE)
    },
    numeric(1)
  )

  group_stats <- group_stats %>%
    arrange(category_median, category_min, category_max)

  if (!ref_label %in% group_stats$category) {
    stop(
      "STOP: Lowest-risk category label from Step_03_1_0 is not found in ",
      category_var, ": ", ref_label,
      ". Available labels: ",
      paste(group_stats$category, collapse = " | ")
    )
  }

  ref_row <- group_stats %>%
    filter(category == ref_label)

  if (nrow(ref_row) != 1L || ref_row$category_n[1] < 1L) {
    stop("STOP: Invalid lowest-risk category data for ", food_var)
  }

  overall_min <- min(x, na.rm = TRUE)
  overall_max <- max(x, na.rm = TRUE)
  span <- overall_max - overall_min
  tol <- max(1e-10, abs(span) * 1e-8)

  ref_min <- ref_row$category_min[1]
  ref_max <- ref_row$category_max[1]
  ref_median <- median(
    x[!is.na(cat) & cat == ref_label],
    na.rm = TRUE
  )

  ref_rank <- match(ref_label, group_stats$category)
  n_categories <- nrow(group_stats)

  if (abs(ref_max - ref_min) <= tol) {
    direction <- "point"
  } else if (ref_rank == 1L) {
    direction <- "low"
  } else if (ref_rank == n_categories) {
    direction <- "high"
  } else {
    direction <- "middle"
  }

  if (direction == "low") {
    value1 <- ref_max
    value2 <- NA_real_
  } else if (direction == "high") {
    value1 <- ref_min
    value2 <- NA_real_
  } else {
    value1 <- ref_min
    value2 <- ref_max
  }

  list(
    category_stats = group_stats,
    reference_category = ref_label,
    reference_category_rank = ref_rank,
    n_categories = n_categories,
    reference_category_n = ref_row$category_n[1],
    reference_category_events = ref_row$category_events[1],
    reference_category_median = ref_median,
    reference_category_min = ref_min,
    reference_category_max = ref_max,
    overall_min = overall_min,
    overall_max = overall_max,
    range_direction = direction,
    value1 = value1,
    value2 = value2
  )
}

# -----------------------------------------------------------------------------
# 4. Fit the 15 RCS models and create the original-style plots
# -----------------------------------------------------------------------------

col <- cols4all::c4a("friendly13")
col_line <- cols4all::c4a("friendly13")

# Exact colours requested for the three empirical intake-percentile markers.
PERCENTILE_COLORS <- c(
  P10 = "#77AADD",
  P50 = "#009988",
  P90 = "#9467BD"
)
if (!identical(
  unname(as.character(col[3:5])),
  unname(PERCENTILE_COLORS)
)) {
  stop(
    "STOP: The friendly13 palette indices col[3:5] no longer match the ",
    "requested P10/P50/P90 colours. Check the palette before plotting."
  )
}

range_results <- vector("list", nrow(FOOD_MAP))
model_results <- vector("list", nrow(FOOD_MAP))
plot_list <- vector("list", nrow(FOOD_MAP))
names(plot_list) <- FOOD_MAP$food

reference_audit_list <- vector("list", nrow(FOOD_MAP))
zero_knot_audit_list <- vector("list", nrow(FOOD_MAP))
percentile_plot_audit_list <- vector("list", nrow(FOOD_MAP))
legacy_prediction_list <- vector("list", nrow(FOOD_MAP))
names(legacy_prediction_list) <- FOOD_MAP$food

for (i in seq_len(nrow(FOOD_MAP))) {

  food_var <- FOOD_MAP$food[i]
  category_var <- FOOD_MAP$category_var[i]
  display_name <- FOOD_MAP$display_name[i]
  ref_label <- REFERENCE_SOURCE$reference_category_label[
    match(food_var, REFERENCE_SOURCE$food)
  ]

  message(
    "\n[RCS ", i, "/15] ", food_var,
    " | lowest-risk category = ", ref_label
  )

  missing_vars <- setdiff(
    c(
      "depression_time", "depression_status",
      food_var, category_var, MODEL3_COVARIATES
    ),
    names(ukb_data_depression_rcs)
  )
  if (length(missing_vars) > 0L) {
    stop(
      "STOP: Missing variables for ", food_var, ": ",
      paste(missing_vars, collapse = ", ")
    )
  }

  required_vars <- unique(c(
    "depression_time", "depression_status",
    food_var, category_var, MODEL3_COVARIATES
  ))

  keep <- stats::complete.cases(
    ukb_data_depression_rcs[, required_vars, drop = FALSE]
  )
  fit_data <- ukb_data_depression_rcs[keep, , drop = FALSE]

  if (nrow(fit_data) < 100L ||
      sum(fit_data$depression_status == 1, na.rm = TRUE) < 20L) {
    stop("STOP: Insufficient complete-case training data/events for ", food_var)
  }

  x <- as.numeric(fit_data[[food_var]])
  if (any(!is.finite(x))) {
    stop("STOP: Non-finite food intake values remain for ", food_var)
  }
  if (any(x < -1e-8)) {
    stop("STOP: Negative food intake values detected for ", food_var)
  }
  if (length(unique(x)) < 3L) {
    stop("STOP: Fewer than three unique intake values for ", food_var)
  }

  # The selected reference category is taken directly from Step_03_1_0.
  # Its median is recomputed in the exact RCS complete-case training sample.
  category_info <- get_category_interval(
    data = fit_data,
    food_var = food_var,
    category_var = category_var,
    ref_label = ref_label
  )

  ref_value <- category_info$reference_category_median
  # The original script used min(food intake) as the datadist reference.
  # Re-estimate that legacy reference on the full frozen training cohort,
  # matching the original script's datadist source.
  training_min_all <- min(
    as.numeric(ukb_data_depression_rcs[[food_var]]), na.rm = TRUE
  )
  training_max_all <- max(
    as.numeric(ukb_data_depression_rcs[[food_var]]), na.rm = TRUE
  )
  legacy_reference_value <- training_min_all
  if (!is.finite(legacy_reference_value) ||
      !is.finite(training_max_all) ||
      training_max_all < training_min_all ||
      training_min_all < -1e-8) {
    stop("STOP: Invalid full-training intake support for ", food_var)
  }
  legacy_reference_type <- if (abs(legacy_reference_value) <= 1e-10) {
    "zero_intake"
  } else {
    "minimum_observed_intake"
  }
  ref_label_median <- parse_leading_numeric(ref_label)

  if (!is.finite(ref_value) ||
      ref_value < min(x) - 1e-8 ||
      ref_value > max(x) + 1e-8) {
    stop("STOP: Invalid lowest-risk category median for ", food_var)
  }

  if (length(ref_label_median) == 1L &&
      is.finite(ref_label_median) &&
      abs(ref_label_median - ref_value) >
        max(0.03, abs(ref_value) * 0.03)) {
    warning(
      "Rounded median in the Step_03_1_0 category label differs from the ",
      "complete-case training median for ", food_var,
      ". The exact training median is used as the RCS reference."
    )
  }

  # Full training-sample percentiles are reported for every component.
  p_values <- as.numeric(stats::quantile(
    x, probs = c(0.10, 0.50, 0.90),
    na.rm = TRUE, names = FALSE, type = 7
  ))
  names(p_values) <- c("P10", "P50", "P90")

  zero_prop <- mean(x <= 1e-12)
  knot_info <- choose_three_knots(x)
  knots <- as.numeric(knot_info$knots)

  if (length(knots) != 3L ||
      any(!is.finite(knots)) ||
      any(diff(knots) <= 0)) {
    stop("STOP: Invalid three-knot specification for ", food_var)
  }

  # P10/P90 define the central prediction support, matching the original
  # Predict-based display range; the fitted model still uses all training data.
  support_values <- as.numeric(stats::quantile(
    x,
    probs = PREDICTION_SUPPORT_PROBS,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
  support_low <- support_values[1]
  support_high <- support_values[2]
  support_method <- "P10_P90"

  # If a very large zero mass makes P10 and P90 identical, retain the same
  # zero-aware analysis but expand the prediction support to P1-P99.
  if (!is.finite(support_low) ||
      !is.finite(support_high) ||
      support_high <= support_low) {
    support_values <- as.numeric(stats::quantile(
      x,
      probs = c(0.01, 0.99),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    ))
    support_low <- support_values[1]
    support_high <- support_values[2]
    support_method <- "P1_P99_zero_mass_fallback"
  }

  if (!is.finite(support_low) ||
      !is.finite(support_high) ||
      support_high <= support_low) {
    stop(
      "STOP: Exposure has no usable central prediction support for ",
      food_var, ". A reliable continuous RCS cannot be displayed."
    )
  }

  # Use a food-specific datadist so the reference point is the median of the
  # Step_03_1_0 lowest-risk category, not automatically zero/minimum intake.
  ddist <- rms::datadist(
    fit_data[, unique(c(food_var, MODEL3_COVARIATES)), drop = FALSE]
  )
  ddist$limits["Adjust to", food_var] <- ref_value
  options(datadist = "ddist")

  knot_text <- paste(sprintf("%.15g", knots), collapse = ", ")

  rcs_formula <- stats::as.formula(paste0(
    "survival::Surv(depression_time, depression_status == 1) ~ ",
    "rcs(", food_var, ", c(", knot_text, ")) + ",
    paste(MODEL3_COVARIATES, collapse = " + ")
  ))

  fit_rcs3 <- rms::cph(
    rcs_formula,
    data = fit_data,
    x = TRUE,
    y = TRUE,
    surv = TRUE,
    method = "efron"
  )

  # Refit the identical linear predictor with the original minimum-intake
  # reference. Coefficients must be identical; only the displayed HR origin
  # changes. This model is used to reconstruct the legacy interval boundaries,
  # not as the primary Figure S4 reference or the sole basis for final scoring.
  # Match the original script's datadist source: the full frozen training
  # analytic cohort, not only the model-specific complete-case subset.
  ddist_legacy <- rms::datadist(
    ukb_data_depression_rcs[, unique(c(food_var, MODEL3_COVARIATES)), drop = FALSE]
  )
  ddist_legacy$limits["Adjust to", food_var] <- legacy_reference_value
  options(datadist = "ddist_legacy")
  fit_rcs3_legacy <- rms::cph(
    rcs_formula,
    data = fit_data,
    x = TRUE,
    y = TRUE,
    surv = TRUE,
    method = "efron"
  )
  options(datadist = "ddist")

  coef_main <- stats::coef(fit_rcs3)
  coef_legacy <- stats::coef(fit_rcs3_legacy)
  if (!identical(names(coef_main), names(coef_legacy)) ||
      max(abs(coef_main - coef_legacy), na.rm = TRUE) > 1e-8) {
    stop("STOP: Legacy-reference refit changed model coefficients for ", food_var)
  }

  # Do not use fit_rcs3$n or fit_rcs3$nevent here:
  # rms::cph may store n as a two-element vector (observations, events),
  # while $nevent is not a reliable scalar field. Derive both directly from
  # the exact complete-case data passed to cph().
  model_n_actual <- as.integer(nrow(fit_data))
  model_events_actual <- as.integer(sum(
    fit_data$depression_status == 1,
    na.rm = TRUE
  ))

  if (length(model_n_actual) != 1L ||
      length(model_events_actual) != 1L ||
      model_n_actual < 1L ||
      model_events_actual < 0L ||
      model_events_actual > model_n_actual) {
    stop("STOP: Invalid model N/events derived from fit_data for ", food_var)
  }

  anova_rcs3 <- as.matrix(stats::anova(fit_rcs3))
  if (nrow(anova_rcs3) < 2L || ncol(anova_rcs3) < 3L) {
    stop("STOP: Unexpected rms::anova output for ", food_var)
  }

  p_overall <- as.numeric(anova_rcs3[1, 3])
  p_non_linear <- as.numeric(anova_rcs3[2, 3])

  # Prediction grid: dense P10-P90 support plus exact P10/P50/P90, knots and
  # the frozen lowest-risk-category median. The graph styling itself is kept.
  grid_low <- min(
    support_low, ref_value, legacy_reference_value,
    category_info$value1, na.rm = TRUE
  )
  grid_high <- max(
    support_high,
    ref_value,
    legacy_reference_value,
    category_info$value1,
    category_info$value2,
    na.rm = TRUE
  )

  pred_grid <- sort(unique(c(
    seq(support_low, support_high, length.out = PREDICTION_GRID_N),
    seq(grid_low, grid_high, length.out = PREDICTION_GRID_N),
    p_values,
    knots,
    ref_value,
    legacy_reference_value,
    category_info$value1,
    category_info$value2
  )))
  pred_grid <- pred_grid[
    is.finite(pred_grid) &
      pred_grid >= min(min(x), legacy_reference_value) - 1e-8 &
      pred_grid <= max(x) + 1e-8
  ]

  pred_args <- c(
    list(fit_rcs3),
    stats::setNames(list(pred_grid), food_var),
    list(
      fun = exp,
      type = "predictions",
      ref.zero = TRUE,
      conf.int = 0.95
    )
  )
  fit_rcs3_pre <- as.data.frame(do.call(rms::Predict, pred_args))

  if (!food_var %in% names(fit_rcs3_pre) ||
      !all(c("yhat", "lower", "upper") %in% names(fit_rcs3_pre))) {
    stop("STOP: Unexpected rms::Predict output for ", food_var)
  }

  fit_rcs3_pre <- fit_rcs3_pre %>%
    arrange(.data[[food_var]]) %>%
    distinct(.data[[food_var]], .keep_all = TRUE)

  # Switch the datadist option back to the legacy reference while creating
  # legacy predictions, then restore the primary lowest-risk reference.
  options(datadist = "ddist_legacy")
  legacy_pred_args <- c(
    list(fit_rcs3_legacy),
    stats::setNames(list(pred_grid), food_var),
    list(
      fun = exp,
      type = "predictions",
      ref.zero = TRUE,
      conf.int = 0.95
    )
  )
  legacy_pred <- as.data.frame(do.call(rms::Predict, legacy_pred_args)) %>%
    arrange(.data[[food_var]]) %>%
    distinct(.data[[food_var]], .keep_all = TRUE)
  if (!all(c(food_var, "yhat", "lower", "upper") %in% names(legacy_pred)) ||
      any(!is.finite(legacy_pred$yhat)) ||
      any(!is.finite(legacy_pred$lower)) ||
      any(!is.finite(legacy_pred$upper))) {
    stop("STOP: Invalid legacy-reference predictions for ", food_var)
  }
  legacy_ref_args <- c(
    list(fit_rcs3_legacy),
    stats::setNames(list(legacy_reference_value), food_var),
    list(fun = exp, type = "predictions", ref.zero = TRUE, conf.int = 0.95)
  )
  legacy_ref_prediction <- as.data.frame(
    do.call(rms::Predict, legacy_ref_args)
  )
  if (nrow(legacy_ref_prediction) != 1L ||
      !is.finite(legacy_ref_prediction$yhat[1]) ||
      abs(legacy_ref_prediction$yhat[1] - 1) > 1e-6) {
    stop("STOP: Legacy-reference HR is not 1 at its reference for ", food_var)
  }
  options(datadist = "ddist")
  legacy_pred <- legacy_pred %>%
    mutate(
      food = food_var,
      display_name = display_name,
      legacy_reference_value = legacy_reference_value,
      legacy_reference_type = legacy_reference_type
    )
  legacy_prediction_list[[food_var]] <- legacy_pred
  utils::write.csv(
    legacy_pred,
    file.path(PRED_DIR, paste0(food_var, "_rcs3_legacy_min_reference.csv")),
    row.names = FALSE,
    na = ""
  )

  # Hard QC: HR at the exact lowest-risk-category median must be 1.
  ref_args <- c(
    list(fit_rcs3),
    stats::setNames(list(ref_value), food_var),
    list(
      fun = exp,
      type = "predictions",
      ref.zero = TRUE,
      conf.int = 0.95
    )
  )
  ref_prediction <- as.data.frame(do.call(rms::Predict, ref_args))

  if (nrow(ref_prediction) != 1L ||
      !is.finite(ref_prediction$yhat[1]) ||
      abs(ref_prediction$yhat[1] - 1) > 1e-6) {
    stop(
      "STOP: RCS HR at the Step_03_1_0 lowest-risk-category median is not 1 ",
      "for ", food_var, "."
    )
  }

  if (!any(abs(fit_rcs3_pre[[food_var]] - ref_value) < 1e-10)) {
    fit_rcs3_pre <- bind_rows(fit_rcs3_pre, ref_prediction) %>%
      arrange(.data[[food_var]])
  }

  # RCS shape and plateau are used to validate a legacy-anchored interval;
  # they do not automatically narrow it to a literal 5% intersection.
  pred_support <- fit_rcs3_pre %>%
    filter(
      .data[[food_var]] >= support_low - 1e-10,
      .data[[food_var]] <= support_high + 1e-10
    ) %>%
    arrange(.data[[food_var]])

  shape_info <- classify_rcs_shape(
    pred = pred_support,
    food_var = food_var,
    p_overall = p_overall,
    p_non_linear = p_non_linear
  )

  plateau_list <- lapply(RANGE_TOLERANCES, function(tol) {
    get_minimum_risk_plateau(
      x = pred_support[[food_var]],
      y = pred_support$yhat,
      tolerance = tol
    )
  })
  names(plateau_list) <- paste0(
    "tol_", sprintf("%02d", round(RANGE_TOLERANCES * 100))
  )

  plateau_05 <- plateau_list[["tol_05"]]
  plateau_10 <- plateau_list[["tol_10"]]
  plateau_20 <- plateau_list[["tol_20"]]

  category_low <- category_info$reference_category_min
  category_high <- category_info$reference_category_max
  category_direction <- category_info$range_direction

  plateau_overlap <- is.finite(plateau_05["lower"]) &&
    is.finite(plateau_05["upper"]) &&
    plateau_05["upper"] >= category_low &&
    plateau_05["lower"] <= category_high

  plateau_overlap_low <- if (plateau_overlap) {
    max(category_low, unname(plateau_05["lower"]))
  } else {
    NA_real_
  }
  plateau_overlap_high <- if (plateau_overlap) {
    min(category_high, unname(plateau_05["upper"]))
  } else {
    NA_real_
  }

  # Check whether the categorical lowest-risk group and RCS curve have
  # directionally compatible patterns. This is an audit flag, not a filter
  # that silently removes a food component.
  direction_consistency <- case_when(
    # A single-value Cox category is retained as a point rule and does not
    # require spline refinement.
    category_direction == "point" &&
      category_low == category_high ~ "zero_or_single_value_category",
    # Statistical support must be checked BEFORE directional compatibility.
    # In V3.2 this condition appeared after the direction-specific branches,
    # allowing an unsupported reverse-J/U-shaped appearance to pass the gate.
    !isTRUE(shape_info$range_eligible) ~ "RCS_not_confirmatory",
    category_direction == "low" &&
      shape_info$curve_shape %in% c(
        "J-shaped_increasing-risk", "monotonic-increasing-risk"
      ) ~ "consistent",
    category_direction == "high" &&
      shape_info$curve_shape %in% c(
        "reverse-J-shaped_decreasing-risk", "monotonic-decreasing-risk"
      ) ~ "consistent",
    category_direction == "middle" &&
      shape_info$curve_shape == "U-shaped" ~ "consistent",
    TRUE ~ "review_category_RCS_discordance"
  )

  # ---------------------------------------------------------------------------
  # Legacy interval reconstruction (exact original decision rules).
  # The original manuscript script used min(intake) as the reference and
  # selected the first/last grid points for which the entire pointwise 95% CI
  # was below HR=1. Low-fat dairy used the first point with the entire CI above
  # HR=1 and scored intake <= value1. Sweetened beverages also retained only
  # value1 in the original food_rcs_range object and were scored intake <= value1.
  # ---------------------------------------------------------------------------
  get_legacy_grid_boundary <- function(dat, direction = c("below1", "above1"),
                                       which_end = c("first", "last")) {
    direction <- match.arg(direction)
    which_end <- match.arg(which_end)
    idx <- if (direction == "below1") {
      which(dat$yhat < 1 & dat$lower < 1 & dat$upper < 1)
    } else {
      which(dat$yhat > 1 & dat$lower > 1 & dat$upper > 1)
    }
    if (length(idx) == 0L) return(NA_real_)
    if (which_end == "first") return(dat[[food_var]][min(idx)])
    dat[[food_var]][max(idx)]
  }

  legacy_beneficial_idx <- which(
    legacy_pred$yhat < 1 & legacy_pred$lower < 1 & legacy_pred$upper < 1
  )
  legacy_harmful_idx <- which(
    legacy_pred$yhat > 1 & legacy_pred$lower > 1 & legacy_pred$upper > 1
  )
  count_runs <- function(idx) {
    if (length(idx) == 0L) return(0L)
    as.integer(sum(c(TRUE, diff(idx) > 1L)))
  }

  if (food_var %in% c("low_fat_dairy", "sweetened_bev")) {
    legacy_value1 <- get_legacy_grid_boundary(
      legacy_pred, direction = "above1", which_end = "first"
    )
    legacy_value2 <- NA_real_
    legacy_range_direction <- "low"
    legacy_range_rule <- if (food_var == "low_fat_dairy") {
      "Original rule: first point with yhat/lower/upper > 1; score x <= value1"
    } else {
      paste0(
        "Original stored rule: value1 is first point with yhat/lower/upper > 1; ",
        "value2 was not retained in the original range object; score x <= value1"
      )
    }
    legacy_region_runs <- count_runs(legacy_harmful_idx)
    legacy_boundary_status <- if (is.finite(legacy_value1)) {
      "EXTRACTED_ORIGINAL_ONE_SIDED_LOW_CUTOFF"
    } else {
      "NO_LEGACY_CUTOFF_FOUND"
    }
  } else {
    legacy_value1 <- get_legacy_grid_boundary(
      legacy_pred, direction = "below1", which_end = "first"
    )
    legacy_value2 <- get_legacy_grid_boundary(
      legacy_pred, direction = "below1", which_end = "last"
    )
    legacy_range_direction <- "middle"
    legacy_range_rule <- paste(
      "Original rule: first/last points with yhat/lower/upper < 1;",
      "the entire pointwise 95% CI is below 1."
    )
    legacy_region_runs <- count_runs(legacy_beneficial_idx)
    legacy_boundary_status <- if (!is.finite(legacy_value1) ||
                                  !is.finite(legacy_value2)) {
      "NO_LEGACY_BENEFICIAL_INTERVAL_FOUND"
    } else if (legacy_value2 < legacy_value1) {
      "INVALID_LEGACY_BOUNDARY_ORDER"
    } else if (legacy_region_runs > 1L) {
      "DISJOINT_LEGACY_REGIONS_OUTER_ENVELOPE_NOT_USED"
    } else {
      "EXTRACTED_SINGLE_CONNECTED_LEGACY_INTERVAL"
    }
  }

  legacy_score_lower <- if (legacy_range_direction == "low") {
    training_min_all
  } else {
    legacy_value1
  }
  legacy_score_upper <- if (legacy_range_direction == "low") {
    legacy_value1
  } else {
    legacy_value2
  }
  legacy_interval_valid <- is.finite(legacy_value1) && (
    if (legacy_range_direction == "low") {
      TRUE
    } else {
      is.finite(legacy_value2) &&
        legacy_value2 >= legacy_value1 &&
        legacy_region_runs == 1L
    }
  )

  category_score_lower <- if (category_direction == "low") {
    training_min_all
  } else {
    category_low
  }
  category_score_upper <- if (category_direction == "high") {
    training_max_all
  } else if (category_direction == "low") {
    category_high
  } else {
    category_high
  }

  interval_tol <- max(
    1e-8,
    abs(category_info$overall_max - category_info$overall_min) * 1e-8
  )
  legacy_category_overlap <- legacy_interval_valid &&
    max(legacy_score_lower, category_score_lower) <=
      min(legacy_score_upper, category_score_upper) + interval_tol
  legacy_plateau_overlap <- legacy_interval_valid &&
    is.finite(plateau_05["lower"]) &&
    is.finite(plateau_05["upper"]) &&
    max(legacy_score_lower, unname(plateau_05["lower"])) <=
      min(legacy_score_upper, unname(plateau_05["upper"])) + interval_tol

  # Legacy interval retention is deliberately conditional. A supported
  # high-intake protective pattern is handled separately: the legacy lower
  # boundary is retained, but the finite legacy upper boundary is opened.
  legacy_shape_compatible <- direction_consistency %in% c(
    "consistent", "RCS_not_confirmatory", "zero_or_single_value_category"
  )
  legacy_candidate_accepted <- legacy_interval_valid &&
    legacy_category_overlap && legacy_shape_compatible &&
    (legacy_plateau_overlap || !isTRUE(shape_info$range_eligible) ||
       category_direction == "high")
  # If the corrected RCS is significantly and monotonically protective toward
  # higher intake, it can bridge a gap between the old finite legacy interval
  # and the highest-intake Cox category. This is the specific fish/seafood
  # safeguard against retaining an artificial legacy upper cap.
  legacy_high_tail_bridge <- category_direction == "high" &&
    direction_consistency == "consistent" &&
    legacy_interval_valid &&
    legacy_range_direction == "middle" &&
    shape_info$curve_shape %in% c(
      "reverse-J-shaped_decreasing-risk", "monotonic-decreasing-risk"
    ) &&
    is.finite(category_low) &&
    legacy_value1 <= category_low + interval_tol &&
    legacy_value2 < category_low - interval_tol

  category_value1 <- category_info$value1
  category_value2 <- category_info$value2
  value1 <- category_value1
  value2 <- category_value2
  range_direction <- category_direction
  range_source <- "cox_category_fallback_no_compatible_legacy_interval"
  legacy_interval_retained <- FALSE
  refined_by_RCS_05 <- FALSE
  legacy_validation_status <- if (!legacy_interval_valid) {
    legacy_boundary_status
  } else if (!legacy_category_overlap) {
    "LEGACY_INTERVAL_DOES_NOT_OVERLAP_COX_LOWEST_RISK_CATEGORY"
  } else if (!legacy_shape_compatible) {
    "LEGACY_INTERVAL_DISCORDANT_WITH_COX_RCS_DIRECTION"
  } else if (!legacy_plateau_overlap && isTRUE(shape_info$range_eligible) &&
             category_direction != "high") {
    "LEGACY_INTERVAL_DOES_NOT_OVERLAP_CORRECTED_RCS_5PCT_PLATEAU"
  } else {
    "LEGACY_INTERVAL_ELIGIBLE_FOR_RETENTION"
  }

  if (category_direction == "high" &&
      direction_consistency == "consistent" &&
      legacy_range_direction == "middle" &&
      (legacy_candidate_accepted || legacy_high_tail_bridge)) {
    # Preserve the old lower boundary when it is compatible with the selected
    # high-intake Cox category, but do not retain a finite upper cutoff when
    # the corrected RCS supports decreasing risk toward higher intake.
    value1 <- min(legacy_value1, category_low)
    value2 <- NA_real_
    range_direction <- "high"
    range_source <- "legacy_lower_bound_retained_open_upper_tail_for_high_intake_pattern"
    legacy_interval_retained <- TRUE
    legacy_validation_status <- if (legacy_high_tail_bridge) {
      "LEGACY_LOWER_BOUND_RETAINED_WITH_MONOTONIC_HIGH_TAIL_BRIDGE"
    } else {
      "LEGACY_LOWER_BOUND_RETAINED_UPPER_TAIL_OPENED"
    }
  } else if (legacy_candidate_accepted &&
             legacy_range_direction == "low" &&
             category_direction %in% c("low", "point")) {
    # Preserve the original low-intake score direction. If the Cox lowest-risk
    # category extends farther than the old cutoff, expand only enough to keep
    # that category included; never reverse the old score direction.
    value1 <- max(legacy_value1, category_high)
    value2 <- NA_real_
    range_direction <- "low"
    range_source <- "legacy_low_cutoff_retained_and_category_anchored"
    legacy_interval_retained <- TRUE
    legacy_validation_status <- "LEGACY_LOW_CUTOFF_RETAINED"
  } else if (legacy_candidate_accepted &&
             legacy_range_direction == "middle" &&
             category_direction %in% c("low", "middle", "point")) {
    # Preserve the original first/last CI-defined interval rather than forcing
    # a narrower Cox-category/RCS-plateau intersection.
    value1 <- legacy_value1
    value2 <- legacy_value2
    range_direction <- "middle"
    range_source <- "legacy_CI_interval_retained_after_category_RCS_validation"
    legacy_interval_retained <- TRUE
    legacy_validation_status <- "LEGACY_CI_INTERVAL_RETAINED"
  } else if (category_direction == "point") {
    range_source <- "cox_category_point_rule_no_compatible_legacy_interval"
  } else if (legacy_interval_valid && !legacy_category_overlap) {
    range_source <- "cox_category_fallback_legacy_category_discordance"
  } else if (legacy_interval_valid && !legacy_plateau_overlap &&
             isTRUE(shape_info$range_eligible) &&
             category_direction != "high") {
    range_source <- "cox_category_fallback_legacy_RCS_plateau_discordance"
  }

  # Preserve the unmodified V5 algorithm result for side-by-side auditing.
  v5_derived_value1 <- value1
  v5_derived_value2 <- value2
  v5_derived_range_direction <- range_direction
  v5_derived_range_source <- range_source
  v5_derived_score_rule <- if (range_direction == "low") {
    paste0("x <= ", fmt_num(value1))
  } else if (range_direction == "high") {
    paste0("x >= ", fmt_num(value1))
  } else if (range_direction == "middle") {
    paste0(fmt_num(value1), " <= x <= ", fmt_num(value2))
  } else {
    paste0("x = ", fmt_num(value1))
  }

  anchor_row <- V51_RANGE_ANCHORS[
    match(food_var, V51_RANGE_ANCHORS$food), , drop = FALSE
  ]
  if (nrow(anchor_row) != 1L || anyNA(anchor_row$range_direction)) {
    stop("STOP: Missing anchor definition for ", food_var)
  }
  anchor_value1 <- as.numeric(anchor_row$value1[1])
  anchor_value2 <- as.numeric(anchor_row$value2[1])
  anchor_direction <- as.character(anchor_row$range_direction[1])
  anchor_direction_basis <- as.character(anchor_row$direction_basis[1])
  if (identical(anchor_direction, "v5_derived")) {
    anchor_direction <- as.character(v5_derived_range_direction)
    if (!anchor_direction %in% c("low", "high")) {
      stop(
        "STOP: The one-sided processed_meats anchor requires a V5-derived ",
        "low/high direction, but V5 derived: ", anchor_direction,
        ". Review the direction before scoring."
      )
    }
  }

  # The primary boundaries are fixed to the supplied anchor table.
  # RCS/Cox results remain available as diagnostics and do not move these values.
  value1 <- anchor_value1
  value2 <- anchor_value2
  range_direction <- anchor_direction
  range_source <- "user_prespecified_anchor_interval"
  legacy_interval_retained <- FALSE
  refined_by_RCS_05 <- FALSE
  # Keep the original V5 legacy-validation status intact as a diagnostic.
  # The final range source below separately records that applies the
  # prespecified anchor table.

  anchor_tol <- max(1e-8, abs(training_max_all - training_min_all) * 1e-8)
  if (value1 < training_min_all - anchor_tol ||
      value1 > training_max_all + anchor_tol ||
      (range_direction == "middle" &&
       (value2 < training_min_all - anchor_tol ||
        value2 > training_max_all + anchor_tol))) {
    stop(
      "STOP: anchor boundary is outside the observed training support for ",
      food_var, ". Anchor=", value1, ", ", value2,
      "; support=", training_min_all, " to ", training_max_all
    )
  }

  score_lower <- case_when(
    range_direction == "low" ~ training_min_all,
    range_direction == "high" ~ value1,
    range_direction %in% c("middle", "point") ~ value1,
    TRUE ~ NA_real_
  )
  score_upper <- case_when(
    range_direction == "low" ~ value1,
    range_direction == "high" ~ Inf,
    range_direction %in% c("middle", "point") ~ value2,
    TRUE ~ NA_real_
  )
  score_rule <- case_when(
    range_direction == "low" ~ paste0("x <= ", fmt_num(value1)),
    range_direction == "high" ~ paste0("x >= ", fmt_num(value1)),
    range_direction == "middle" ~ paste0(
      fmt_num(value1), " <= x <= ", fmt_num(value2)
    ),
    range_direction == "point" ~ paste0("x = ", fmt_num(value1)),
    TRUE ~ "not_defined"
  )

  legacy_rule_label <- if (legacy_range_direction == "low" &&
                           is.finite(legacy_value1)) {
    paste0("x <= ", fmt_num(legacy_value1))
  } else if (legacy_range_direction == "middle" &&
             is.finite(legacy_value1) && is.finite(legacy_value2)) {
    paste0(fmt_num(legacy_value1), " <= x <= ", fmt_num(legacy_value2))
  } else {
    "not_defined"
  }

  food_rcs_range <- tibble(
    food = food_var,
    value1 = value1,
    value2 = value2,
    range_direction = range_direction,
    range_source = range_source,
    score_lower = score_lower,
    score_upper = score_upper,
    score_rule = score_rule,
    v5_derived_value1 = v5_derived_value1,
    v5_derived_value2 = v5_derived_value2,
    v5_derived_range_direction = v5_derived_range_direction,
    v5_derived_range_source = v5_derived_range_source,
    v5_derived_score_rule = v5_derived_score_rule,
    anchor_value1 = anchor_value1,
    anchor_value2 = anchor_value2,
    anchor_range_direction = anchor_direction,
    anchor_direction_basis = anchor_direction_basis,
    anchor_override_applied = TRUE,
    category_value1 = category_value1,
    category_value2 = category_value2,
    refined_by_RCS_5pct = refined_by_RCS_05,
    category_RCS_direction_consistency = direction_consistency,
    reference_category = ref_label,
    reference_median = ref_value,
    category_range_low = category_low,
    category_range_high = category_high,
    legacy_reference_value = legacy_reference_value,
    legacy_reference_type = legacy_reference_type,
    legacy_value1 = legacy_value1,
    legacy_value2 = legacy_value2,
    legacy_range_direction = legacy_range_direction,
    legacy_range_rule = legacy_range_rule,
    legacy_rule_label = legacy_rule_label,
    legacy_boundary_status = legacy_boundary_status,
    legacy_region_runs = legacy_region_runs,
    legacy_score_lower = legacy_score_lower,
    legacy_score_upper = legacy_score_upper,
    legacy_interval_valid = legacy_interval_valid,
    legacy_overlaps_category = legacy_category_overlap,
    legacy_overlaps_5pct_plateau = legacy_plateau_overlap,
    legacy_high_tail_bridge = legacy_high_tail_bridge,
    legacy_interval_retained = legacy_interval_retained,
    legacy_validation_status = legacy_validation_status
  )

  # Export the complete prediction grid with exact percentile/knot/reference
  # flags for reproducibility and for downstream histogram overlays if needed.
  prediction_export <- fit_rcs3_pre %>%
    mutate(
      food = food_var,
      display_name = display_name,
      is_reference_value = abs(.data[[food_var]] - ref_value) < 1e-10,
      is_p10 = abs(.data[[food_var]] - p_values["P10"]) < 1e-10,
      is_p50 = abs(.data[[food_var]] - p_values["P50"]) < 1e-10,
      is_p90 = abs(.data[[food_var]] - p_values["P90"]) < 1e-10,
      is_rcs_knot_1 = abs(.data[[food_var]] - knots[1]) < 1e-10,
      is_rcs_knot_2 = abs(.data[[food_var]] - knots[2]) < 1e-10,
      is_rcs_knot_3 = abs(.data[[food_var]] - knots[3]) < 1e-10
    )

  write.csv(
    prediction_export,
    file.path(PRED_DIR, paste0(food_var, "_rcs3.csv")),
    row.names = FALSE
  )

  range_results[[i]] <- food_rcs_range

  model_results[[i]] <- tibble(
    group = FOOD_MAP$group[i],
    food = food_var,
    display_name = display_name,
    category_variable = category_var,
    reference_category_label = ref_label,
    reference_median_label = ref_label_median,
    reference_median_used = ref_value,
    reference_HR = as.numeric(ref_prediction$yhat[1]),
    reference_category_rank = category_info$reference_category_rank,
    number_of_categories = category_info$n_categories,
    reference_category_n = category_info$reference_category_n,
    reference_category_events = category_info$reference_category_events,
    category_range_low = category_low,
    category_range_high = category_high,
    category_range_direction = category_direction,
    overall_training_min = training_min_all,
    overall_training_max = training_max_all,
    model_n = model_n_actual,
    model_events = model_events_actual,
    zero_intake_n = sum(x <= 1e-12),
    zero_intake_percent = 100 * zero_prop,
    p10 = unname(p_values["P10"]),
    p50 = unname(p_values["P50"]),
    p90 = unname(p_values["P90"]),
    rcs_knot_1 = knots[1],
    rcs_knot_2 = knots[2],
    rcs_knot_3 = knots[3],
    rcs_knot_method = knot_info$method,
    positive_p10 = unname(knot_info$positive_quantiles[1]),
    positive_p50 = unname(knot_info$positive_quantiles[2]),
    positive_p90 = unname(knot_info$positive_quantiles[3]),
    prediction_support_low = support_low,
    prediction_support_high = support_high,
    prediction_support_method = support_method,
    p_overall = p_overall,
    p_non_linear = p_non_linear,
    curve_shape = shape_info$curve_shape,
    shape_evidence = shape_info$shape_evidence,
    x_min = shape_info$x_min,
    hr_min = shape_info$hr_min,
    left_tail_hr_ratio_to_min = shape_info$left_tail_ratio,
    right_tail_hr_ratio_to_min = shape_info$right_tail_ratio,
    rcs_range_eligible = shape_info$range_eligible,
    plateau_05_low = unname(plateau_05["lower"]),
    plateau_05_high = unname(plateau_05["upper"]),
    plateau_10_low = unname(plateau_10["lower"]),
    plateau_10_high = unname(plateau_10["upper"]),
    plateau_20_low = unname(plateau_20["lower"]),
    plateau_20_high = unname(plateau_20["upper"]),
    plateau_05_overlaps_reference_category = plateau_overlap,
    plateau_05_overlap_low = plateau_overlap_low,
    plateau_05_overlap_high = plateau_overlap_high,
    category_RCS_direction_consistency = direction_consistency,
    primary_range_source = range_source,
    primary_range_direction = range_direction,
    primary_range_value1 = value1,
    primary_range_value2 = value2,
    v5_derived_value1 = v5_derived_value1,
    v5_derived_value2 = v5_derived_value2,
    v5_derived_range_direction = v5_derived_range_direction,
    v5_derived_range_source = v5_derived_range_source,
    v5_derived_score_rule = v5_derived_score_rule,
    anchor_value1 = anchor_value1,
    anchor_value2 = anchor_value2,
    anchor_range_direction = anchor_direction,
    anchor_direction_basis = anchor_direction_basis,
    anchor_override_applied = TRUE,
    category_range_value1 = category_value1,
    category_range_value2 = category_value2,
    refined_by_RCS_5pct = refined_by_RCS_05,
    score_lower = score_lower,
    score_upper = score_upper,
    score_rule = score_rule,
    legacy_reference_value = legacy_reference_value,
    legacy_reference_type = legacy_reference_type,
    legacy_value1 = legacy_value1,
    legacy_value2 = legacy_value2,
    legacy_range_direction = legacy_range_direction,
    legacy_range_rule = legacy_range_rule,
    legacy_boundary_status = legacy_boundary_status,
    legacy_region_runs = legacy_region_runs,
    legacy_overlaps_category = legacy_category_overlap,
    legacy_overlaps_5pct_plateau = legacy_plateau_overlap,
    legacy_high_tail_bridge = legacy_high_tail_bridge,
    legacy_interval_retained = legacy_interval_retained,
    legacy_validation_status = legacy_validation_status
  )

  reference_audit_list[[i]] <- tibble(
    group = FOOD_MAP$group[i],
    food = food_var,
    category_variable = category_var,
    reference_category_label = ref_label,
    reference_median_label = ref_label_median,
    reference_median_used = ref_value,
    reference_category_min = category_low,
    reference_category_max = category_high,
    reference_category_n = category_info$reference_category_n,
    reference_category_events = category_info$reference_category_events,
    reference_category_rank = category_info$reference_category_rank,
    number_of_categories = category_info$n_categories,
    category_range_direction = category_direction,
    range_value1 = value1,
    range_value2 = value2,
    range_direction = range_direction,
    anchor_direction_basis = anchor_direction_basis
  )

  zero_knot_audit_list[[i]] <- tibble(
    food = food_var,
    n_complete_case = length(x),
    zero_intake_n = sum(x <= 1e-12),
    zero_intake_percent = 100 * zero_prop,
    positive_intake_n = sum(x > 1e-12),
    full_sample_p10 = unname(p_values["P10"]),
    full_sample_p50 = unname(p_values["P50"]),
    full_sample_p90 = unname(p_values["P90"]),
    actual_knot_1 = knots[1],
    actual_knot_2 = knots[2],
    actual_knot_3 = knots[3],
    knot_method = knot_info$method,
    prediction_support_low = support_low,
    prediction_support_high = support_high,
    prediction_support_method = support_method
  )

  # ---------------------------------------------------------------------------
  # V5 RCS plots: retain the fitted curve and interval boundaries, while
  # explicitly displaying empirical P10/P50/P90 values above the plotting area.
  # The percentile lines use the exact colours requested by the user.
  # ---------------------------------------------------------------------------

  label_overall_rcs3 <- paste0(
    "P for overall ",
    format_p_value(p_overall)
  )
  label_non_linear_rcs3 <- paste0(
    "P for non-linear ",
    format_p_value(p_non_linear)
  )

  plot_x_max <- if (is.finite(FOOD_MAP$x_limit_fixed[i])) {
    FOOD_MAP$x_limit_fixed[i]
  } else {
    max(
      fit_rcs3_pre[[food_var]],
      ref_value,
      value1,
      ifelse(is.finite(value2), value2, -Inf),
      na.rm = TRUE
    )
  }
  plot_x_max <- max(plot_x_max, 0.01)

  if (any(p_values < -1e-8 | p_values > plot_x_max + 1e-8)) {
    stop(
      "STOP: At least one P10/P50/P90 marker would fall outside the displayed ",
      "x-axis for ", food_var, ". Review the fixed x-axis limit."
    )
  }

  annotation_x <- FOOD_MAP$annotation_x[i]
  plot_support_upper <- fit_rcs3_pre$upper[
    fit_rcs3_pre[[food_var]] >= support_low - 1e-10 &
      fit_rcs3_pre[[food_var]] <= support_high + 1e-10
  ]
  plot_text_y <- max(plot_support_upper, na.rm = TRUE)

  plot_y_span <- diff(range(
    c(fit_rcs3_pre$lower, fit_rcs3_pre$upper, 1),
    finite = TRUE
  ))
  if (!is.finite(plot_y_span) || plot_y_span <= 0) {
    plot_y_span <- max(0.10, abs(plot_text_y) * 0.10)
  }
  percentile_label_y <- plot_text_y + 0.18 * plot_y_span

  percentile_plot_data <- tibble(
    percentile = factor(
      names(p_values),
      levels = c("P10", "P50", "P90")
    ),
    percentile_value = as.numeric(p_values)
  )

  # Labels are arranged in a dedicated top row in P10/P50/P90 order. This
  # prevents overlap when two or more empirical percentiles are tied.
  percentile_label_data <- percentile_plot_data %>%
    mutate(
      label = paste0(
        as.character(percentile), " = ",
        vapply(percentile_value, fmt_num, character(1))
      ),
      label_x = plot_x_max * c(0.17, 0.50, 0.83),
      label_y = percentile_label_y
    )

  p <- ggplot() +
    # The confidence ribbon is deliberately drawn before the fitted curve.
    geom_ribbon(
      data = fit_rcs3_pre,
      aes(x = .data[[food_var]], ymin = lower, ymax = upper),
      alpha = 0.10,
      fill = col_line[1]
    ) +
    geom_vline(
      data = percentile_plot_data,
      aes(xintercept = percentile_value, colour = percentile),
      linetype = "dashed",
      linewidth = 0.65,
      alpha = 0.95,
      show.legend = TRUE
    ) +
    geom_vline(
      xintercept = ref_value,
      linetype = "dotdash",
      linewidth = 0.65,
      colour = "black",
      alpha = 0.85
    ) +
    geom_line(
      data = fit_rcs3_pre,
      aes(x = .data[[food_var]], y = yhat),
      linetype = "solid",
      linewidth = 1.0,
      alpha = 0.85,
      colour = col_line[1]
    ) +
    geom_text(
      data = percentile_label_data,
      aes(
        x = label_x,
        y = label_y,
        label = label,
        colour = percentile
      ),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 4.0,
      fontface = "bold",
      vjust = 0
    ) +
    geom_text(
      data = NULL,
      aes(
        x = annotation_x,
        y = plot_text_y,
        label = paste(
          label_overall_rcs3,
          label_non_linear_rcs3,
          paste0("RCS reference = ", fmt_num(ref_value)),
          sep = "\n"
        )
      ),
      size = 5.0,
      fontface = "bold"
    ) +
    scale_colour_manual(
      name = "Training-set percentile",
      values = PERCENTILE_COLORS,
      breaks = c("P10", "P50", "P90"),
      drop = FALSE
    ) +
    guides(
      colour = guide_legend(
        direction = "horizontal",
        title.position = "left",
        override.aes = list(
          linetype = "dashed",
          linewidth = 1.1,
          alpha = 1
        )
      )
    ) +
    scale_y_continuous(
      name = "Hazard ratio (95% CI)",
      expand = expansion(mult = c(0.02, 0.32))
    ) +
    scale_x_continuous(
      name = display_name,
      breaks = seq(0, plot_x_max, by = FOOD_MAP$x_break[i]),
      limits = c(0, plot_x_max)
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(color = "black", size = 18, face = "bold"),
      axis.title = element_text(color = "black", size = 18, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.width = grid::unit(0.42, "inches"),
      plot.margin = margin(9, 12, 8, 12)
    ) +
    geom_hline(
      yintercept = 1,
      linetype = 2,
      linewidth = 0.5,
      colour = col_line[2],
      alpha = 0.4
    )

  get_curve_y <- function(x0) {
    stats::approx(
      x = fit_rcs3_pre[[food_var]],
      y = fit_rcs3_pre$yhat,
      xout = x0,
      rule = 2,
      ties = "ordered"
    )$y
  }

  if (range_direction == "middle") {
    boundary_x <- c(value1, value2)
    boundary_y <- vapply(boundary_x, get_curve_y, numeric(1))
    boundary_data <- tibble(
      x = boundary_x,
      y = boundary_y,
      label = sprintf("%.2f", boundary_x)
    )

    p <- p +
      geom_vline(
        xintercept = boundary_x,
        linetype = 2,
        size = 0.5,
        colour = col_line[2],
        alpha = 0.6
      ) +
      geom_point(
        data = boundary_data,
        aes(x = x, y = y),
        colour = col_line[3],
        size = 3
      ) +
      geom_text(
        data = boundary_data,
        aes(x = x, y = y, label = label),
        vjust = -1,
        colour = "black",
        size = 6,
        fontface = "bold"
      )
  } else {
    boundary_y <- get_curve_y(value1)
    boundary_data <- tibble(
      x = value1,
      y = boundary_y,
      label = sprintf("%.2f", value1)
    )

    p <- p +
      geom_vline(
        xintercept = value1,
        linetype = 2,
        size = 0.5,
        colour = col_line[2],
        alpha = 0.6
      ) +
      geom_point(
        data = boundary_data,
        aes(x = x, y = y),
        colour = col_line[3],
        size = 3
      ) +
      geom_text(
        data = boundary_data,
        aes(x = x, y = y, label = label),
        vjust = -1,
        colour = "black",
        size = 6,
        fontface = "bold"
      )
  }

  plot_list[[food_var]] <- p

  percentile_plot_audit_list[[i]] <- tibble(
    food = food_var,
    display_name = display_name,
    p10 = unname(p_values["P10"]),
    p50 = unname(p_values["P50"]),
    p90 = unname(p_values["P90"]),
    p10_colour = unname(PERCENTILE_COLORS["P10"]),
    p50_colour = unname(PERCENTILE_COLORS["P50"]),
    p90_colour = unname(PERCENTILE_COLORS["P90"]),
    percentile_lines_dashed = TRUE,
    numeric_labels_at_top = TRUE,
    rcs_reference_value = ref_value,
    reference_vertical_line_shown = TRUE,
    reference_numeric_label_shown = TRUE,
    legend_included = TRUE,
    all_percentiles_within_displayed_x_axis = all(
      p_values >= -1e-8 & p_values <= plot_x_max + 1e-8
    )
  )

  # Preserve original individual PDF dimensions and file naming.
  ggsave(
    filename = file.path(PLOT_DIR, paste0(food_var, "_rcs3.pdf")),
    plot = p,
    width = 10,
    height = 6,
    device = "pdf"
  )

  message(
    "[DONE] ", food_var,
    " | N=", model_n_actual,
    " | events=", model_events_actual,
    " | reference=", fmt_num(ref_value),
    " | primary interval=",
    if (range_direction == "middle") {
      paste0(fmt_num(value1), "–", fmt_num(value2))
    } else if (range_direction == "low") {
      paste0("≤ ", fmt_num(value1))
    } else if (range_direction == "high") {
      paste0("≥ ", fmt_num(value1))
    } else {
      fmt_num(value1)
    },
    " | RCS shape=", shape_info$curve_shape
  )
}

# -----------------------------------------------------------------------------
# 5. Export range/model/reference/knot audits and combine the 15 original plots
# -----------------------------------------------------------------------------

food_rcs_range_all <- bind_rows(range_results)
rcs_model_summary <- bind_rows(model_results)
rcs_reference_audit <- bind_rows(reference_audit_list)
rcs_zero_knot_audit <- bind_rows(zero_knot_audit_list)
rcs_percentile_plot_audit <- bind_rows(percentile_plot_audit_list)
legacy_predictions_all <- bind_rows(legacy_prediction_list, .id = "food_list_name")

if (nrow(legacy_predictions_all) == 0L ||
    length(unique(legacy_predictions_all$food)) != 15L) {
  stop("STOP: Legacy-reference prediction export does not cover all 15 foods.")
}
utils::write.csv(
  legacy_predictions_all,
  file.path(AUDIT_DIR, "RCS_legacy_minimum_reference_predictions_all_15.csv"),
  row.names = FALSE,
  na = ""
)

if (nrow(rcs_percentile_plot_audit) != 15L ||
    any(!rcs_percentile_plot_audit$all_percentiles_within_displayed_x_axis)) {
  stop("STOP: P10/P50/P90 plot-annotation audit failed.")
}

write.csv(
  rcs_percentile_plot_audit,
  file.path(AUDIT_DIR, "RCS_percentile_plot_annotation_QC.csv"),
  row.names = FALSE,
  na = ""
)

# Exact component reconciliation is performed on normalized character IDs.
# This avoids false failures caused only by vector attributes/names while still
# rejecting missing, extra, duplicated, or unexpected food components.
expected_food_ids <- trimws(as.character(unname(FOOD_MAP$food)))

reconcile_component_order <- function(x, object_name) {
  if (!is.data.frame(x) || !"food" %in% names(x)) {
    stop("STOP: ", object_name, " must be a data.frame with a food column.")
  }

  observed_food_ids <- trimws(as.character(unname(x$food)))

  if (length(observed_food_ids) != length(expected_food_ids) ||
      anyNA(observed_food_ids) ||
      any(observed_food_ids == "") ||
      anyDuplicated(observed_food_ids) > 0L) {
    stop(
      "STOP: ", object_name,
      " has an invalid row count, missing IDs, or duplicated components. ",
      "Observed rows=", nrow(x),
      "; expected rows=", length(expected_food_ids),
      "; observed IDs=",
      paste(observed_food_ids, collapse = ", ")
    )
  }

  missing_ids <- setdiff(expected_food_ids, observed_food_ids)
  extra_ids <- setdiff(observed_food_ids, expected_food_ids)

  if (length(missing_ids) > 0L || length(extra_ids) > 0L) {
    stop(
      "STOP: Component set mismatch in ", object_name,
      ". Missing: ", paste(missing_ids, collapse = ", "),
      ". Extra: ", paste(extra_ids, collapse = ", ")
    )
  }

  x <- x[match(expected_food_ids, observed_food_ids), , drop = FALSE]
  rownames(x) <- NULL
  x
}

food_rcs_range_all <- reconcile_component_order(
  food_rcs_range_all, "food_rcs_range_all"
)
rcs_model_summary <- reconcile_component_order(
  rcs_model_summary, "rcs_model_summary"
)

if (nrow(food_rcs_range_all) != 15L ||
    nrow(rcs_model_summary) != 15L) {
  stop("STOP: Final RCS outputs do not contain exactly 15 unique components.")
}

v51_anchor_comparison <- food_rcs_range_all %>%
  transmute(
    food,
    V5_derived_value1 = v5_derived_value1,
    V5_derived_value2 = v5_derived_value2,
    V5_derived_direction = v5_derived_range_direction,
    V5_derived_rule = v5_derived_score_rule,
    anchor_value1 = anchor_value1,
    anchor_value2 = anchor_value2,
    anchor_direction = anchor_range_direction,
    direction_basis = anchor_direction_basis,
    final_rule = score_rule,
    value1_absolute_difference = abs(anchor_value1 - v5_derived_value1),
    value2_absolute_difference = ifelse(
      is.finite(anchor_value2) & is.finite(v5_derived_value2),
      abs(anchor_value2 - v5_derived_value2), NA_real_
    ),
    direction_matches_V5 = anchor_range_direction == v5_derived_range_direction,
    one_sided_status_matches_V5 = is.na(anchor_value2) == is.na(v5_derived_value2),
    V5_derived_source = v5_derived_range_source,
    training_support_min = rcs_model_summary$overall_training_min,
    training_support_max = rcs_model_summary$overall_training_max,
    anchor_override_applied = anchor_override_applied
  )
write.csv(
  v51_anchor_comparison,
  file.path(AUDIT_DIR, "user_anchor_vs_V5_derived_interval_comparison.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

if (any(!is.finite(rcs_model_summary$reference_median_used)) ||
    any(abs(rcs_model_summary$reference_HR - 1) > 1e-6)) {
  stop("STOP: One or more RCS reference-point QC checks failed.")
}

write.csv(
  food_rcs_range_all,
  file.path(OUT_DIR, "food_rcs_range_all.csv"),
  row.names = FALSE
)
write.csv(
  food_rcs_range_all %>%
    select(
      food, value1, value2, range_direction, range_source, score_rule,
      v5_derived_value1, v5_derived_value2, v5_derived_range_direction,
      v5_derived_range_source, anchor_direction_basis, anchor_override_applied,
      legacy_reference_value, legacy_reference_type,
      legacy_value1, legacy_value2, legacy_range_direction,
      legacy_range_rule, legacy_rule_label, legacy_boundary_status,
      legacy_region_runs, legacy_overlaps_category,
      legacy_overlaps_5pct_plateau, legacy_high_tail_bridge,
      legacy_interval_retained,
      legacy_validation_status
    ),
  file.path(AUDIT_DIR, "food_rcs_range_legacy_comparison.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  rcs_model_summary,
  file.path(AUDIT_DIR, "RCS_model_summary_and_range_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  rcs_reference_audit,
  file.path(AUDIT_DIR, "RCS_lowest_risk_reference_category_audit.csv"),
  row.names = FALSE
)
write.csv(
  rcs_zero_knot_audit,
  file.path(AUDIT_DIR, "RCS_zero_intake_P10_P50_P90_knot_audit.csv"),
  row.names = FALSE
)

# Preserve the original 5-column x 3-row combined figure layout and size.
# Individual plots retain their percentile legends; the combined figure uses
# one shared legend to avoid repeating the same key 15 times.
rcs_legend <- cowplot::get_legend(
  plot_list[[FOOD_MAP$food[1]]] +
    theme(legend.position = "bottom")
)
if (is.null(rcs_legend)) {
  stop("STOP: Could not extract the shared P10/P50/P90 legend.")
}

rcs_plots <- lapply(
  unname(plot_list[FOOD_MAP$food]),
  function(z) z + theme(legend.position = "none")
)

rcs_combined_body <- cowplot::plot_grid(
  plotlist = rcs_plots,
  ncol = 5,
  nrow = 3,
  labels = LETTERS[1:15],
  label_size = 36,
  label_fontface = "bold",
  label_y = 1.1
)

rcs_combined <- cowplot::plot_grid(
  rcs_combined_body,
  rcs_legend,
  ncol = 1,
  rel_heights = c(1, 0.045)
)

ggsave(
  filename = file.path(PLOT_DIR, "rcs3_all_food_combined.pdf"),
  plot = rcs_combined,
  width = 38,
  height = 15,
  device = "pdf"
)

# -----------------------------------------------------------------------------
# 6. Original interval-range visualization, now supporting low/middle/high
# -----------------------------------------------------------------------------

food_rcs_range_display <- food_rcs_range_all %>%
  left_join(
    FOOD_MAP %>% select(food, display_name),
    by = "food"
  ) %>%
  mutate(
    min_value = 0,
    max_value = vapply(
      food,
      function(f) {
        max(
          ukb_data_depression_rcs[[f]],
          na.rm = TRUE
        )
      },
      numeric(1)
    ),
    display_name = factor(display_name, levels = FOOD_MAP$display_name)
  )

intervals <- food_rcs_range_display %>%
  rowwise() %>%
  mutate(
    y_start = case_when(
      range_direction == "middle" ~ value1,
      range_direction == "low" ~ min_value,
      range_direction == "high" ~ value1,
      range_direction == "point" ~ value1,
      TRUE ~ NA_real_
    ),
    y_end = case_when(
      range_direction == "middle" ~ value2,
      range_direction == "low" ~ value1,
      range_direction == "high" ~ max_value,
      range_direction == "point" ~ value2,
      TRUE ~ NA_real_
    ),
    label = case_when(
      range_direction == "middle" ~ paste0(fmt_num(value1), "–", fmt_num(value2)),
      range_direction == "low" ~ paste0(fmt_num(min_value), "–", fmt_num(value1)),
      range_direction == "high" ~ paste0("≥ ", fmt_num(value1)),
      range_direction == "point" ~ fmt_num(value1),
      TRUE ~ NA_character_
    ),
    y_mid = (y_start + y_end) / 2
  ) %>%
  ungroup() %>%
  filter(is.finite(y_start), is.finite(y_end))

# Original interval plot styling is retained.
p_food_group <- ggplot(food_rcs_range_display) +
  geom_segment(
    aes(
      x = forcats::fct_rev(factor(display_name)),
      xend = display_name,
      y = min_value,
      yend = max_value
    ),
    size = 2.5,
    color = col[1]
  ) +
  geom_segment(
    data = filter(food_rcs_range_display, range_direction == "middle"),
    aes(x = display_name, xend = display_name, y = value1, yend = value2),
    size = 2.5,
    color = col[3]
  ) +
  geom_segment(
    data = filter(food_rcs_range_display, range_direction == "low"),
    aes(x = display_name, xend = display_name, y = min_value, yend = value1),
    size = 2.5,
    color = col[3]
  ) +
  geom_segment(
    data = filter(food_rcs_range_display, range_direction == "high"),
    aes(x = display_name, xend = display_name, y = value1, yend = max_value),
    size = 2.5,
    color = col[3]
  ) +
  geom_segment(
    data = filter(food_rcs_range_display, range_direction == "point"),
    aes(x = display_name, xend = display_name, y = value1, yend = value2),
    size = 2.5,
    color = col[3]
  ) +
  geom_point(
    data = filter(food_rcs_range_display, range_direction == "point"),
    aes(x = display_name, y = value1),
    size = 2.5,
    color = col[3]
  ) +
  geom_text(
    data = intervals,
    aes(x = display_name, y = y_mid, label = label),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  coord_flip() +
  labs(y = "Consumption level (serving/day)", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "bold", colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black"),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11, face = "bold")
  )

ggsave(
  filename = file.path(PLOT_DIR, "p_food_group.pdf"),
  plot = p_food_group,
  width = 8,
  height = 4.2,
  device = "pdf"
)

# Original tile-plot style retained; range direction is now explicit.
step <- 0.25

tile_list <- lapply(seq_len(nrow(food_rcs_range_display)), function(i) {
  z <- food_rcs_range_display[i, , drop = FALSE]
  x_low <- z$min_value[1]
  x_high <- z$max_value[1]

  starts <- seq(x_low, x_high, by = step)
  if (length(starts) == 0L) starts <- x_low
  ends <- pmin(starts + step, x_high)
  mids <- (starts + ends) / 2

  v1 <- z$value1[1]
  v2 <- z$value2[1]
  direction <- z$range_direction[1]

  point <- if (direction == "middle") {
    as.numeric(mids >= v1 & mids <= v2)
  } else if (direction == "low") {
    as.numeric(mids <= v1)
  } else if (direction == "high") {
    as.numeric(mids >= v1)
  } else if (direction == "point") {
    # A single-value interval must not be expanded to an entire 0.25-serving
    # tile. It is drawn separately as an exact point marker below.
    rep(NA_real_, length(mids))
  } else {
    rep(NA_real_, length(mids))
  }

  tibble(
    food = as.character(z$display_name[1]),
    x_start = starts,
    x_end = ends,
    x_mid = mids,
    point = point
  )
})

food_rcs_range_all_tiles <- bind_rows(tile_list) %>%
  mutate(
    food = factor(food, levels = rev(FOOD_MAP$display_name))
  )

p_food_group_tiles <- ggplot(
  food_rcs_range_all_tiles,
  aes(x = x_mid, y = food, fill = factor(point))
) +
  geom_tile(
    data = filter(food_rcs_range_all_tiles, !is.na(point)),
    aes(width = pmax(x_end - x_start, 0.001)),
    color = "white",
    height = 0.9
  ) +
  geom_point(
    data = filter(food_rcs_range_display, range_direction == "point"),
    aes(x = value1, y = display_name),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.5,
    stroke = 0.8,
    fill = col[2],
    colour = "black"
  ) +
  scale_fill_manual(
    values = c("0" = col[1], "1" = col[2]),
    name = "Point",
    na.value = "grey85"
  ) +
  scale_x_continuous(name = "Servings/day", expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

ggsave(
  filename = file.path(PLOT_DIR, "p_food_group_tiles.pdf"),
  plot = p_food_group_tiles,
  width = 8,
  height = 4.2,
  device = "pdf"
)


# -----------------------------------------------------------------------------
# 6B. Final interval derivation table, frozen binary scores, and figures
# -----------------------------------------------------------------------------
# This block exports the user-prespecified anchor intervals. The fitted
# Cox/RCS models and V5.0-derived intervals are retained as diagnostics; they
# do not alter the fixed primary boundaries.
#
# Binary score:
#   1 = intake lies within the frozen primary interval
#   0 = intake lies outside the frozen primary interval
#   NA = intake is missing/non-finite; no imputation is performed.
# -----------------------------------------------------------------------------

message("\n[6B/7] Exporting interval derivation, binary scores, and new figures ...")

# 6B.1 Validate the exact 15-component interval objects before deriving scores.
required_range_columns <- c(
  "food", "value1", "value2", "range_direction", "range_source",
  "score_lower", "score_upper", "score_rule",
  "category_RCS_direction_consistency", "refined_by_RCS_5pct",
  "legacy_value1", "legacy_value2", "legacy_boundary_status",
  "legacy_high_tail_bridge", "legacy_interval_retained",
  "legacy_validation_status", "v5_derived_value1", "v5_derived_value2",
  "v5_derived_range_direction", "anchor_value1", "anchor_value2",
  "anchor_range_direction", "anchor_direction_basis", "anchor_override_applied"
)
required_model_columns <- c(
  "food", "group", "display_name", "reference_category_label",
  "reference_median_used", "category_range_low", "category_range_high",
  "category_range_direction", "overall_training_min", "overall_training_max",
  "p10", "p50", "p90", "p_overall", "p_non_linear", "curve_shape",
  "shape_evidence", "x_min", "rcs_range_eligible", "plateau_05_low", "plateau_05_high",
  "plateau_10_low", "plateau_10_high", "plateau_20_low", "plateau_20_high",
  "plateau_05_overlaps_reference_category",
  "category_RCS_direction_consistency", "primary_range_source",
  "primary_range_direction", "primary_range_value1", "primary_range_value2",
  "v5_derived_value1", "v5_derived_value2", "v5_derived_range_direction",
  "anchor_value1", "anchor_value2", "anchor_range_direction",
  "anchor_direction_basis", "anchor_override_applied",
  "score_lower", "score_upper", "score_rule", "model_n", "model_events"
)

missing_range_columns <- setdiff(required_range_columns, names(food_rcs_range_all))
missing_model_columns <- setdiff(required_model_columns, names(rcs_model_summary))

if (length(missing_range_columns) > 0L) {
  stop(
    "STOP: food_rcs_range_all is missing columns required for final scoring: ",
    paste(missing_range_columns, collapse = ", ")
  )
}
if (length(missing_model_columns) > 0L) {
  stop(
    "STOP: rcs_model_summary is missing columns required for the derivation table: ",
    paste(missing_model_columns, collapse = ", ")
  )
}

food_rcs_range_all <- reconcile_component_order(
  food_rcs_range_all, "food_rcs_range_all before final scoring"
)
rcs_model_summary <- reconcile_component_order(
  rcs_model_summary, "rcs_model_summary before final scoring"
)

if (nrow(food_rcs_range_all) != 15L ||
    nrow(rcs_model_summary) != 15L ||
    !identical(as.character(food_rcs_range_all$food), as.character(FOOD_MAP$food)) ||
    !identical(as.character(rcs_model_summary$food), as.character(FOOD_MAP$food))) {
  stop("STOP: Final interval/scoring derivation requires the exact frozen 15 foods.")
}

if (!all(food_rcs_range_all$anchor_override_applied) ||
    !identical(as.numeric(food_rcs_range_all$value1),
               as.numeric(V51_RANGE_ANCHORS$value1)) ||
    !identical(as.character(food_rcs_range_all$range_direction),
               as.character(food_rcs_range_all$anchor_range_direction)) ||
    !all((is.na(food_rcs_range_all$value2) & is.na(V51_RANGE_ANCHORS$value2)) |
         (!is.na(food_rcs_range_all$value2) & !is.na(V51_RANGE_ANCHORS$value2) &
          abs(food_rcs_range_all$value2 - V51_RANGE_ANCHORS$value2) <= 1e-12))) {
  stop("STOP: Final interval values do not exactly match the user anchor table.")
}

if (anyNA(food_rcs_range_all$range_direction) ||
    any(food_rcs_range_all$range_direction == "not_defined") ||
    anyNA(food_rcs_range_all$score_lower) ||
    anyNA(food_rcs_range_all$score_upper) ||
    anyNA(food_rcs_range_all$score_rule) ||
    any(food_rcs_range_all$score_rule == "not_defined")) {
  stop("STOP: At least one food component has an undefined binary scoring interval.")
}

# Human-readable interval formatters. Numeric fields remain available in the
# same table; these labels are intended for manuscript/supplementary review.
format_category_interval <- function(direction, category_low, category_high,
                                      overall_min, overall_max) {
  if (direction == "low") {
    return(paste0("[", fmt_num(overall_min), ", ", fmt_num(category_high), "]"))
  }
  if (direction == "high") {
    return(paste0("[", fmt_num(category_low), ", ", fmt_num(overall_max), "]"))
  }
  if (direction == "middle") {
    return(paste0("[", fmt_num(category_low), ", ", fmt_num(category_high), "]"))
  }
  if (direction == "point") {
    return(paste0("{", fmt_num(category_low), "}"))
  }
  "Not defined"
}

format_final_interval <- function(direction, lower, upper, overall_min, overall_max) {
  if (direction == "low") {
    return(paste0("[", fmt_num(overall_min), ", ", fmt_num(upper), "]"))
  }
  if (direction == "high") {
    return(paste0("[", fmt_num(lower), ", ", fmt_num(overall_max), "] (upper tail)"))
  }
  if (direction == "middle") {
    return(paste0("[", fmt_num(lower), ", ", fmt_num(upper), "]"))
  }
  if (direction == "point") {
    return(paste0("{", fmt_num(lower), "} (single value)"))
  }
  "Not defined"
}

format_rcs_plateau <- function(lower, upper) {
  if (length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper)) {
    return("Not estimable")
  }
  paste0("[", fmt_num(lower), ", ", fmt_num(upper), "]")
}

# 6B.2 Build the one-row-per-food derivation table. It explicitly shows:
#   (i) Cox lowest-risk category and its observed intake range;
#   (ii) the lowest-risk-category median used as the primary RCS reference;
#   (iii) the legacy minimum-intake reference and original CI-based interval;
#   (iv) category overlap, corrected-RCS shape and plateau validation;
#   (v) whether the legacy interval was retained or the Cox category was used;
#   (vi) the final binary scoring rule.
interval_derivation_table <- FOOD_MAP %>%
  select(group, food, display_name, category_var) %>%
  left_join(
    rcs_model_summary %>%
      select(
        food,
        reference_category_label,
        reference_median_used,
        reference_category_rank,
        number_of_categories,
        reference_category_n,
        reference_category_events,
        category_range_low,
        category_range_high,
        category_range_direction,
        overall_training_min,
        overall_training_max,
        model_n,
        model_events,
        p10, p50, p90,
        rcs_knot_1, rcs_knot_2, rcs_knot_3,
        rcs_knot_method,
        p_overall,
        p_non_linear,
        curve_shape,
        shape_evidence,
        x_min,
        rcs_range_eligible,
        plateau_05_low,
        plateau_05_high,
        plateau_10_low,
        plateau_10_high,
        plateau_20_low,
        plateau_20_high,
        plateau_05_overlaps_reference_category,
        category_RCS_direction_consistency,
        primary_range_source,
        primary_range_direction,
        primary_range_value1,
        primary_range_value2,
        v5_derived_value1,
        v5_derived_value2,
        v5_derived_range_direction,
        v5_derived_range_source,
        v5_derived_score_rule,
        anchor_value1,
        anchor_value2,
        anchor_range_direction,
        anchor_direction_basis,
        anchor_override_applied,
        refined_by_RCS_5pct,
        score_lower,
        score_upper,
        score_rule
      ),
    by = "food"
  ) %>%
  left_join(
    food_rcs_range_all %>%
      select(
        food,
        cox_category_value1 = category_value1,
        cox_category_value2 = category_value2,
        final_value1 = value1,
        final_value2 = value2,
        final_range_direction = range_direction,
        final_range_source = range_source,
        final_score_lower = score_lower,
        final_score_upper = score_upper,
        final_score_rule = score_rule,
        final_direction_consistency = category_RCS_direction_consistency,
        final_refined_by_RCS_5pct = refined_by_RCS_5pct,
        legacy_reference_value,
        legacy_reference_type,
        legacy_value1,
        legacy_value2,
        legacy_range_direction,
        legacy_range_rule,
        legacy_rule_label,
        legacy_boundary_status,
        legacy_region_runs,
        legacy_score_lower,
        legacy_score_upper,
        legacy_interval_valid,
        legacy_overlaps_category,
        legacy_overlaps_5pct_plateau,
        legacy_high_tail_bridge,
        legacy_interval_retained,
        legacy_validation_status
      ),
    by = "food",
    suffix = c("", "_range")
  ) %>%
  mutate(
    cox_lowest_risk_category = reference_category_label,
    cox_reference_category_observed_range = paste0(
      "[", vapply(category_range_low, fmt_num, character(1)),
      ", ", vapply(category_range_high, fmt_num, character(1)), "]"
    ),
    cox_candidate_interval = mapply(
      format_category_interval,
      category_range_direction,
      category_range_low,
      category_range_high,
      overall_training_min,
      overall_training_max,
      USE.NAMES = FALSE
    ),
    rcs_reference_value = reference_median_used,
    rcs_5pct_plateau = mapply(
      format_rcs_plateau,
      plateau_05_low,
      plateau_05_high,
      USE.NAMES = FALSE
    ),
    rcs_10pct_plateau = mapply(
      format_rcs_plateau,
      plateau_10_low,
      plateau_10_high,
      USE.NAMES = FALSE
    ),
    rcs_20pct_plateau = mapply(
      format_rcs_plateau,
      plateau_20_low,
      plateau_20_high,
      USE.NAMES = FALSE
    ),
    legacy_candidate_interval = ifelse(
      legacy_interval_valid,
      mapply(
        format_final_interval,
        legacy_range_direction,
        legacy_score_lower,
        legacy_score_upper,
        overall_training_min,
        overall_training_max,
        USE.NAMES = FALSE
      ),
      "Not available"
    ),
    v5_derived_candidate_interval = mapply(
      format_final_interval,
      v5_derived_range_direction,
      v5_derived_value1,
      v5_derived_value2,
      overall_training_min,
      overall_training_max,
      USE.NAMES = FALSE
    ),
    anchor_comparison = paste0(
      "fixed anchor; direction=", anchor_range_direction,
      "; V5.0 direction=", v5_derived_range_direction,
      "; direction_match=", anchor_range_direction == v5_derived_range_direction
    ),
    interval_derivation = "USER-ANCHORED: supplied 15-component scoring interval fixed before score construction",
    interval_derivation_reason = paste0(
      "The primary binary scoring interval was set to the prespecified anchor table supplied for this revision. ",
      "The interval boundaries were not moved in response to the observed Cox P values, RCS P values, confidence intervals, or score-association results. ",
      "The previous V5 data-derived interval, lowest-risk Cox category, corrected-reference RCS curve, minimum-risk plateaus, and legacy CI-based interval are retained as diagnostics for transparency. ",
      "An operational score-1 interval must not be described as a statistically significant protective range or clinical recommendation unless the corresponding association results independently support that interpretation. ",
      "Direction basis: ", anchor_direction_basis, "."
    ),
    final_primary_interval = mapply(
      format_final_interval,
      final_range_direction,
      final_score_lower,
      final_score_upper,
      overall_training_min,
      overall_training_max,
      USE.NAMES = FALSE
    ),
    binary_score_1_rule = paste0("1 if intake is within: ", final_score_rule),
    binary_score_0_rule = paste0("0 if intake is outside: ", final_score_rule),
    missing_value_rule = "NA; no imputation"
  ) %>%
  arrange(match(food, FOOD_MAP$food)) %>%
  select(
    group,
    food,
    display_name,
    category_variable = category_var,
    cox_lowest_risk_category,
    cox_reference_category_rank = reference_category_rank,
    cox_reference_category_n = reference_category_n,
    cox_reference_category_events = reference_category_events,
    cox_reference_category_observed_range,
    cox_candidate_interval,
    cox_category_min = category_range_low,
    cox_category_max = category_range_high,
    cox_category_direction = category_range_direction,
    legacy_reference_value,
    legacy_reference_type,
    legacy_candidate_interval,
    legacy_range_rule,
    legacy_boundary_status,
    legacy_region_runs,
    legacy_overlaps_category,
    legacy_overlaps_5pct_plateau,
    legacy_high_tail_bridge,
    legacy_interval_retained,
    legacy_validation_status,
    rcs_reference_value,
    p10, p50, p90,
    rcs_knot_1, rcs_knot_2, rcs_knot_3,
    rcs_knot_method,
    p_overall,
    p_non_linear,
    rcs_curve_shape = curve_shape,
    rcs_shape_evidence = shape_evidence,
    rcs_minimum_risk_point = x_min,
    rcs_5pct_plateau,
    rcs_10pct_plateau,
    rcs_20pct_plateau,
    rcs_5pct_overlaps_cox_interval = plateau_05_overlaps_reference_category,
    cox_rcs_direction_consistency = final_direction_consistency,
    v5_derived_candidate_interval,
    anchor_comparison,
    anchor_direction_basis,
    interval_derivation,
    interval_derivation_reason,
    interval_source_code = final_range_source,
    final_interval_direction = final_range_direction,
    final_interval_lower = final_score_lower,
    final_interval_upper = final_score_upper,
    final_primary_interval,
    binary_score_1_rule,
    binary_score_0_rule,
    missing_value_rule,
    model_n,
    model_events
  )

if (nrow(interval_derivation_table) != 15L ||
    anyDuplicated(interval_derivation_table$food) > 0L ||
    !identical(as.character(interval_derivation_table$food),
               as.character(FOOD_MAP$food))) {
  stop("STOP: Interval derivation table failed exact 15-component reconciliation.")
}

write.csv(
  interval_derivation_table,
  file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

# V5 manuscript-ready Word table. The CSV above remains the complete audit
# table; the Word version is a compact, publication-ready presentation table.
format_p_word <- function(p) {
  if (length(p) != 1L || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

format_final_interval_word <- function(direction, lower, upper) {
  if (direction == "high") {
    return(paste0("≥ ", fmt_num(lower), " (open upper tail)"))
  }
  if (direction == "low") {
    return(paste0("≤ ", fmt_num(upper)))
  }
  if (direction == "middle") {
    return(paste0(fmt_num(lower), "–", fmt_num(upper)))
  }
  if (direction == "point") {
    return(paste0("= ", fmt_num(lower)))
  }
  "Not defined"
}

word_table_data <- interval_derivation_table %>%
  mutate(
    Component = display_name,
    `Cox Ref. category` = cox_lowest_risk_category,
    `Cox candidate interval` = cox_candidate_interval,
    `V5.0-derived interval` = v5_derived_candidate_interval,
    `V5.0 legacy diagnostic` = legacy_validation_status,
    `RCS reference` = vapply(rcs_reference_value, fmt_num, character(1)),
    `P10 / P50 / P90` = paste(
      vapply(p10, fmt_num, character(1)),
      vapply(p50, fmt_num, character(1)),
      vapply(p90, fmt_num, character(1)),
      sep = " / "
    ),
    `RCS 5% plateau` = rcs_5pct_plateau,
    `RCS evidence / derivation` = paste0(
      rcs_curve_shape,
      "\nP overall = ", vapply(p_overall, format_p_word, character(1)),
      "; P non-linear = ", vapply(p_non_linear, format_p_word, character(1)),
      "\n",
      "supplied anchor applied; RCS/Cox results are diagnostics, not boundary selectors"
    ),
    `Final interval (score = 1)` = mapply(
      format_final_interval_word,
      final_interval_direction,
      final_interval_lower,
      final_interval_upper,
      USE.NAMES = FALSE
    )
  ) %>%
  select(
    Component,
    `Cox Ref. category`,
    `Cox candidate interval`,
    `V5.0-derived interval`,
    `V5.0 legacy diagnostic`,
    `RCS reference`,
    `P10 / P50 / P90`,
    `RCS 5% plateau`,
    `RCS evidence / derivation`,
    `Final interval (score = 1)`
  )

if (nrow(word_table_data) != 15L ||
    anyNA(word_table_data$Component) ||
    any(word_table_data$Component == "")) {
  stop("STOP: The Word-ready interval table failed 15-component validation.")
}

word_table_title <- paste(
  "Supplementary Table Sx.",
  "Prespecified intake anchors and RCS diagnostic evidence"
)

word_footnotes <- c(
  paste0(
    "Note. All component-specific analyses used the exact frozen Step_03_1_0 ",
    "70% training cohort (N = ", format(nrow(ukb_data_depression_rcs), big.mark = ","),
    "; depression events = ",
    format(sum(ukb_data_depression_rcs$depression_status == 1, na.rm = TRUE),
           big.mark = ","), "). The 30% test set was not used for RCS fitting or ",
    "interval specification. Intake is expressed as servings/day unless otherwise specified."
  ),
  paste0(
    "Abbreviations: Cox, Cox proportional hazards regression; HR, hazard ratio; ",
    "RCS, restricted cubic spline; Ref., reference category; P10/P50/P90, the ",
    "10th/50th/90th percentiles of component intake."
  ),
  paste0(
    "1. Primary RCS reference. The primary Figure S4 curve is referenced to ",
    "the median intake within the Step_03_1_0 lowest-risk Model 2 category, ",
    "recalculated in the component-specific RCS complete-case training sample. ",
    "The HR at this reference value is 1. The empirical P10/P50/P90 values are ",
    "shown separately from the actual spline knots."
  ),
  paste0(
    "2. primary scoring intervals. The 15 primary binary scoring intervals ",
    "are the prespecified component-specific anchors supplied for this revision ",
    "and are listed in the anchor CSV. Finite value1/value2 pairs define ",
    "inclusive middle intervals. For one-sided intervals, the direction is ",
    "explicitly recorded: low means intake <= value1; high means intake >= value1. ",
    "The direction is not inferred from NA alone. For processed meats, the ",
    "2.7 one-sided cutoff preserves the one-sided low/high direction derived ",
    "by the preceding V5 algorithm; this direction is exported in the final ",
    "food_rcs_range_all.csv and the anchor-comparison audit."
  ),
  paste0(
    "3. Role of Cox/RCS analyses. The lowest-risk Cox category, corrected-reference ",
    "RCS shape, pointwise confidence intervals, and 5%/10%/20% minimum-risk ",
    "plateaus are reported as diagnostic evidence. They do not automatically ",
    "move the prespecified anchor boundaries. The prior V5 data-derived ",
    "interval and the legacy minimum-intake-reference CI interval are retained ",
    "in the audit outputs for comparison."
  ),
  paste0(
    "4. Interpretation. An operational score-1 interval is not necessarily a ",
    "statistically significant protective range and is not a clinical intake ",
    "recommendation. A protective association is claimed only where the fitted ",
    "association estimate and its confidence interval provide the relevant ",
    "statistical support. In particular, absence of a statistically significant ",
    "beneficial region is not converted into evidence of benefit by the scoring rule."
  ),
  paste0(
    "5. Binary scoring. A component receives 1 when observed intake satisfies ",
    "its inclusive anchor rule and 0 otherwise. Missing component intake ",
    "remains NA; no imputation is performed. The complete 0-15 RWDID total is ",
    "calculated only for participants with observed values for all 15 components."
  )
)

ft <- flextable::flextable(word_table_data)
ft <- flextable::theme_booktabs(ft)
ft <- flextable::font(ft, fontname = "Arial", part = "all")
ft <- flextable::fontsize(ft, size = 7, part = "all")
ft <- flextable::fontsize(ft, size = 7.5, part = "header")
ft <- flextable::bold(ft, part = "header", bold = TRUE)
ft <- flextable::bg(ft, part = "header", bg = "#EAF2F8")
ft <- flextable::color(ft, part = "header", color = "#1F2937")
ft <- flextable::align(ft, align = "center", part = "header")
ft <- flextable::align(ft, j = c(1, 2, 3, 4, 6, 7, 8, 10), align = "center", part = "body")
ft <- flextable::align(ft, j = c(5, 9), align = "left", part = "body")
ft <- flextable::valign(ft, valign = "center", part = "all")
ft <- flextable::padding(ft, padding = 3, part = "all")
# Word tables wrap cell text automatically when fixed column widths are used.
# flextable does not export a word_wrap() function; do not call it here.
ft <- flextable::width(
  ft,
  j = 1:10,
  width = c(1.05, 1.00, 1.10, 1.00, 1.35, 1.10, 0.72, 0.85, 1.55, 1.65)
)
ft <- flextable::set_table_properties(ft, layout = "fixed", width = 1)
ft <- flextable::fit_to_width(ft, max_width = 10.35)

# Use style names actually present in the default officer template:
# "table title", "Normal", and "heading 2" (case-sensitive).
doc <- officer::read_docx()
doc <- officer::body_set_default_section(
  doc,
  value = officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(
      top = 0.42, bottom = 0.42, left = 0.48, right = 0.48
    )
  )
)
doc <- officer::body_add_par(
  doc,
  word_table_title,
  style = "table title"
)
doc <- officer::body_add_par(
  doc,
  paste0(
    "Fifteen components; frozen 70% training cohort; ",
    "user-prespecified component-specific scoring anchors with Cox/RCS diagnostics."
  ),
  style = "Normal"
)
doc <- flextable::body_add_flextable(doc, value = ft)
doc <- officer::body_add_par(doc, "Footnotes and interval-derivation rules", style = "heading 2")
for (note in word_footnotes) {
  doc <- officer::body_add_par(doc, note, style = "Normal")
}
WORD_TABLE_FILE <- file.path(
  OUT_DIR, "Table_Sx_interval_derivation_summary.docx"
)
print(doc, target = WORD_TABLE_FILE)

if (!file.exists(WORD_TABLE_FILE) || file.info(WORD_TABLE_FILE)$size < 5000) {
  stop("STOP: The formatted Word interval-derivation table was not created.")
}

message("[6B/7] Word supplementary table exported: ", WORD_TABLE_FILE)


# A compact scoring-rule file is also exported for direct use in the next
# RWDID construction step. The primary derivation table remains the full audit.
binary_scoring_rules <- food_rcs_range_all %>%
  left_join(
    FOOD_MAP %>% select(food, group, display_name),
    by = "food"
  ) %>%
  left_join(
    rcs_model_summary %>%
      select(food, overall_training_min, overall_training_max,
             reference_category_label, reference_median_used,
             plateau_05_low, plateau_05_high),
    by = "food"
  ) %>%
  mutate(
    interval_source = case_when(
      range_source == "user_prespecified_anchor_interval" ~
        "USER-PRESPECIFIED ANCHOR",
      range_source == "legacy_CI_interval_retained_after_category_RCS_validation" ~
        "LEGACY CI INTERVAL RETAINED",
      range_source == "legacy_low_cutoff_retained_and_category_anchored" ~
        "LEGACY LOW CUTOFF RETAINED",
      range_source == "legacy_lower_bound_retained_open_upper_tail_for_high_intake_pattern" ~
        "LEGACY LOWER BOUND; OPEN UPPER TAIL",
      TRUE ~ "COX CATEGORY FALLBACK"
    ),
    score_1 = paste0("1 if ", score_rule),
    score_0 = paste0("0 otherwise"),
    missing_rule = "NA; no imputation"
  ) %>%
  arrange(match(food, FOOD_MAP$food)) %>%
  select(
    group, food, display_name,
    reference_category_label, reference_median_used,
    category_range_low, category_range_high,
    plateau_05_low,
    plateau_05_high,
    range_source, interval_source, range_direction,
    value1, value2, score_lower, score_upper,
    legacy_reference_value, legacy_reference_type,
    legacy_value1, legacy_value2, legacy_range_direction,
    legacy_boundary_status, legacy_validation_status,
    legacy_high_tail_bridge, legacy_interval_retained,
    score_rule, score_1, score_0, missing_rule,
    overall_training_min, overall_training_max
  )

# 6B.3 Apply the frozen binary rules to the exact Step_03_1_0 training RDS.
# Missing exposure values remain NA. No new split, refitting, imputation, or
# test-set access occurs in this step.
score_one_component <- function(x, direction, value1, value2, food_name) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  score <- rep(NA_integer_, length(x))
  observed <- !is.na(x)

  if (any(!is.finite(x[observed]))) {
    stop("STOP: Non-finite non-missing exposure values found for ", food_name)
  }
  if (any(x[observed] < -1e-8)) {
    stop("STOP: Negative intake values found while scoring ", food_name)
  }

  if (!any(observed)) {
    return(score)
  }

  scale_value <- max(1, abs(value1), ifelse(is.finite(value2), abs(value2), 0))
  tol <- max(1e-8, scale_value * 1e-8)

  if (direction == "middle") {
    if (!is.finite(value1) || !is.finite(value2) || value2 < value1) {
      stop("STOP: Invalid middle-interval boundaries for ", food_name)
    }
    score[observed] <- as.integer(
      x[observed] >= value1 - tol & x[observed] <= value2 + tol
    )
  } else if (direction == "low") {
    if (!is.finite(value1)) {
      stop("STOP: Invalid low-interval upper boundary for ", food_name)
    }
    score[observed] <- as.integer(x[observed] <= value1 + tol)
  } else if (direction == "high") {
    if (!is.finite(value1)) {
      stop("STOP: Invalid high-interval lower boundary for ", food_name)
    }
    score[observed] <- as.integer(x[observed] >= value1 - tol)
  } else if (direction == "point") {
    if (!is.finite(value1)) {
      stop("STOP: Invalid point interval for ", food_name)
    }
    score[observed] <- as.integer(abs(x[observed] - value1) <= tol)
  } else {
    stop("STOP: Unsupported interval direction for ", food_name, ": ", direction)
  }

  score
}

score_data <- data.frame(
  participant_id = trimws(as.character(ukb_data_depression_rcs[[ID_COL]])),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
names(score_data)[1] <- ID_COL

score_columns <- paste0("score_", FOOD_MAP$food)

for (i in seq_len(nrow(FOOD_MAP))) {
  food_var <- FOOD_MAP$food[i]
  range_row <- food_rcs_range_all[
    match(food_var, food_rcs_range_all$food), ,
    drop = FALSE
  ]

  score_data[[score_columns[i]]] <- score_one_component(
    x = ukb_data_depression_rcs[[food_var]],
    direction = as.character(range_row$range_direction[1]),
    value1 = as.numeric(range_row$value1[1]),
    value2 = as.numeric(range_row$value2[1]),
    food_name = food_var
  )
}

score_matrix <- as.matrix(score_data[, score_columns, drop = FALSE])
storage.mode(score_matrix) <- "numeric"
score_data$n_components_scored <- as.integer(rowSums(!is.na(score_matrix)))

complete_15_components <- score_data$n_components_scored == 15L
score_data$RWDID_total_score_complete15 <- NA_integer_

if (any(complete_15_components)) {
  score_data$RWDID_total_score_complete15[complete_15_components] <-
    as.integer(rowSums(
      score_matrix[complete_15_components, , drop = FALSE],
      na.rm = FALSE
    ))
}

score_data$analysis_split <- "training"

if (nrow(score_data) != nrow(ukb_data_depression_rcs) ||
    anyNA(score_data[[ID_COL]]) ||
    any(score_data[[ID_COL]] == "") ||
    anyDuplicated(score_data[[ID_COL]]) > 0L) {
  stop("STOP: Participant-level binary scoring output failed ID/row reconciliation.")
}

write.csv(
  score_data,
  file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

# Component-level 0/1 distribution for audit and manuscript reporting.
binary_score_distribution <- bind_rows(lapply(seq_len(nrow(FOOD_MAP)), function(i) {
  food_var <- FOOD_MAP$food[i]
  score_vec <- score_data[[score_columns[i]]]
  observed_scores <- score_vec[!is.na(score_vec)]
  n_scored <- length(observed_scores)
  n_zero <- sum(observed_scores == 0L)
  n_one <- sum(observed_scores == 1L)

  tibble(
    group = FOOD_MAP$group[i],
    food = food_var,
    display_name = FOOD_MAP$display_name[i],
    n_training = length(score_vec),
    n_missing = sum(is.na(score_vec)),
    n_scored = n_scored,
    n_score_0 = n_zero,
    n_score_1 = n_one,
    pct_score_0_among_scored = if (n_scored > 0L) 100 * n_zero / n_scored else NA_real_,
    pct_score_1_among_scored = if (n_scored > 0L) 100 * n_one / n_scored else NA_real_
  )
}))

if (nrow(binary_score_distribution) != 15L ||
    any(binary_score_distribution$n_score_0 +
          binary_score_distribution$n_score_1 !=
          binary_score_distribution$n_scored)) {
  stop("STOP: Component-level binary score distribution failed reconciliation.")
}

if (any(binary_score_distribution$n_score_0 == 0L) ||
    any(binary_score_distribution$n_score_1 == 0L)) {
  print(binary_score_distribution)
  stop(
    "STOP: At least one anchor produces no 0/1 variation in the training set. ",
    "Do not estimate a beta coefficient for a non-varying component; review its ",
    "prespecified direction/cutpoint before proceeding."
  )
}

write.csv(
  binary_score_distribution,
  file.path(AUDIT_DIR, "RWDID_binary_score_distribution_training.csv"),
  row.names = FALSE,
  na = ""
)

total_score_observed <- score_data$RWDID_total_score_complete15[
  !is.na(score_data$RWDID_total_score_complete15)
]
total_score_counts <- tabulate(total_score_observed + 1L, nbins = 16L)

total_score_distribution <- tibble(
  RWDID_total_score = 0:15,
  n = as.integer(total_score_counts),
  percent_among_complete15 = if (length(total_score_observed) > 0L) {
    100 * as.integer(total_score_counts) / length(total_score_observed)
  } else {
    NA_real_
  },
  complete_15_component_n = length(total_score_observed),
  incomplete_component_n = sum(!complete_15_components)
)

write.csv(
  total_score_distribution,
  file.path(AUDIT_DIR, "RWDID_total_score_distribution_training.csv"),
  row.names = FALSE,
  na = ""
)

# 6B.4 Construct the body-style 0/1 scoring scheme figure.
# The three columns are Low / Medium / High. For one-sided rules, adjacent
# columns with the same score are merged, matching the supplied example.
make_scoring_scheme_cells <- function(i) {
  food_var <- FOOD_MAP$food[i]
  r <- food_rcs_range_all[
    match(food_var, food_rcs_range_all$food), ,
    drop = FALSE
  ]
  m <- rcs_model_summary[
    match(food_var, rcs_model_summary$food), ,
    drop = FALSE
  ]

  direction <- as.character(r$range_direction[1])
  lower <- as.numeric(r$score_lower[1])
  upper <- as.numeric(r$score_upper[1])
  x_min <- as.numeric(m$overall_training_min[1])
  x_max <- as.numeric(m$overall_training_max[1])
  y <- nrow(FOOD_MAP) - i + 1L

  tol <- max(1e-8, abs(x_max - x_min) * 1e-8, abs(lower) * 1e-8)
  cells <- list()

  add_cell <- function(xmin, xmax, score, label) {
    cells[[length(cells) + 1L]] <<- tibble(
      food = food_var,
      display_name = FOOD_MAP$display_name[i],
      row_y = y,
      xmin = xmin,
      xmax = xmax,
      ymin = y - 0.44,
      ymax = y + 0.44,
      score = as.integer(score),
      label = label
    )
  }

  if (direction == "point") {
    point_value <- lower
    # A single observed value is not a consumption band. Show the whole
    # background as score 0 and add a separate exact-value score-1 marker.
    add_cell(
      0, 3, 0L,
      paste0("Only x = ", fmt_num(point_value), " scores 1")
    )
  } else if (is.infinite(upper)) {
    if (lower <= x_min + tol) {
      add_cell(0, 3, 1L, "All observed intake")
    } else {
      add_cell(0, 1, 0L, paste0("< ", fmt_num(lower)))
      add_cell(1, 3, 1L, paste0("≥ ", fmt_num(lower)))
    }
  } else if (lower <= x_min + tol && upper >= x_max - tol) {
    add_cell(0, 3, 1L, "All observed intake")
  } else if (lower <= x_min + tol) {
    add_cell(0, 2, 1L, paste0("≤ ", fmt_num(upper)))
    add_cell(2, 3, 0L, paste0("> ", fmt_num(upper)))
  } else if (upper >= x_max - tol) {
    add_cell(0, 1, 0L, paste0("< ", fmt_num(lower)))
    add_cell(1, 3, 1L, paste0("≥ ", fmt_num(lower)))
  } else {
    add_cell(0, 1, 0L, paste0("< ", fmt_num(lower)))
    add_cell(1, 2, 1L, paste0(fmt_num(lower), "–", fmt_num(upper)))
    add_cell(2, 3, 0L, paste0("> ", fmt_num(upper)))
  }

  bind_rows(cells)
}

scoring_scheme_cells <- bind_rows(lapply(seq_len(nrow(FOOD_MAP)), make_scoring_scheme_cells)) %>%
  mutate(score_factor = factor(score, levels = c(0, 1), labels = c("0 point", "1 point")))

if (nrow(scoring_scheme_cells) < 15L ||
    any(!scoring_scheme_cells$score %in% c(0L, 1L)) ||
    any(!is.finite(scoring_scheme_cells$xmin)) ||
    any(!is.finite(scoring_scheme_cells$xmax))) {
  stop("STOP: The 15-component binary scoring scheme could not be constructed.")
}

# Apply a 20% tint of each original score colour (80% white mixed with
# the original hue) to create a restrained, print-friendly pastel palette.
mix_with_white <- function(hex, colour_fraction = 0.20) {
  if (!is.finite(colour_fraction) ||
      colour_fraction < 0 || colour_fraction > 1) {
    stop("STOP: colour_fraction must be between 0 and 1.")
  }
  rgb <- grDevices::col2rgb(hex) / 255
  mixed <- 1 - colour_fraction * (1 - rgb)
  grDevices::rgb(mixed[1], mixed[2], mixed[3])
}

SCORE_FILL_COLORS <- c(
  "0" = mix_with_white(col[1], 0.20),
  "1" = mix_with_white(col[2], 0.20)
)

point_marker_data <- food_rcs_range_all %>%
  filter(range_direction == "point") %>%
  left_join(
    rcs_model_summary %>%
      select(food, overall_training_min, overall_training_max),
    by = "food"
  ) %>%
  left_join(
    FOOD_MAP %>% select(food, display_name),
    by = "food"
  ) %>%
  mutate(
    row_y = nrow(FOOD_MAP) - match(food, FOOD_MAP$food) + 1L,
    marker_x = if_else(
      is.finite(overall_training_min) &
        is.finite(overall_training_max) &
        overall_training_max > overall_training_min,
      0.5 + 2 * (value1 - overall_training_min) /
        (overall_training_max - overall_training_min),
      1.5
    )
  )

p_binary_scoring_scheme <- ggplot(scoring_scheme_cells) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = factor(score)
    ),
    colour = "white",
    linewidth = 0.9
  ) +
  geom_point(
    data = point_marker_data,
    aes(x = marker_x, y = row_y),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.6,
    stroke = 0.9,
    fill = SCORE_FILL_COLORS[["1"]],
    colour = "#243B53"
  ) +
  geom_text(
    aes(
      x = (xmin + xmax) / 2,
      y = row_y,
      label = label
    ),
    size = 6.0,
    colour = "black",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = SCORE_FILL_COLORS,
    breaks = c("0", "1"),
    labels = c("0 point", "1 point"),
    name = "Score:"
  ) +
  scale_x_continuous(
    breaks = c(0.5, 1.5, 2.5),
    labels = c("Low", "Medium", "High"),
    limits = c(0, 3),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(nrow(FOOD_MAP)),
    labels = rev(FOOD_MAP$display_name),
    expand = expansion(add = c(0.55, 0.75))
  ) +
  labs(
    title = "Consumption level (servings per day)",
    x = "Grade of consumption level",
    y = NULL,
    tag = "b"
  ) +
  guides(
    fill = guide_legend(
      direction = "horizontal",
      title.position = "left",
      title.hjust = 0,
      override.aes = list(colour = NA)
    )
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(colour = "black", size = 20, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 19, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 22, face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
    plot.tag = element_text(face = "bold", size = 32),
    plot.tag.position = c(0.005, 1.02),
    plot.margin = margin(18, 18, 14, 26)
  )

ggsave(
  filename = file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.pdf"),
  plot = p_binary_scoring_scheme,
  width = 8.33,
  height = 9.2,
  device = "pdf"
)
ggsave(
  filename = file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.png"),
  plot = p_binary_scoring_scheme,
  width = 8.33,
  height = 9.2,
  dpi = 600,
  bg = "white"
)

# 6B.5 Journal-style anchor-and-diagnostic workflow.
# The primary scoring intervals are fixed from the supplied anchor table;
# Cox/RCS results are diagnostic and do not move the anchor boundaries.
workflow_nodes <- tibble(
  node = c("cox", "rcs", "gate", "confirmed", "not_confirmed", "score"),
  xmin = c(0.35, 6.25, 0.75, 0.35, 6.25, 1.65),
  xmax = c(5.75, 11.65, 11.25, 5.75, 11.65, 10.35),
  ymin = c(7.45, 7.45, 5.15, 2.55, 2.55, 0.25),
  ymax = c(9.35, 9.35, 6.85, 4.35, 4.35, 1.75),
  fill = c(
    "#EEF5FB", "#EDF7F5", "#F5F7FA",
    "#EAF5F2", "#F4F5F7", "#EEF3F8"
  ),
  border = c(
    "#77AADD", "#009988", "#526579",
    "#009988", "#64748B", "#334E68"
  )
)

workflow_connectors <- tibble(
  x = c(3.05, 8.95, 3.05, 8.95, 3.05, 8.95),
  y = c(7.45, 7.45, 5.15, 5.15, 2.55, 2.55),
  xend = c(4.20, 7.80, 3.05, 8.95, 4.70, 7.30),
  yend = c(6.85, 6.85, 4.35, 4.35, 1.75, 1.75)
)

workflow_criteria <- tibble(
  xmin = c(1.05, 4.45, 7.85),
  xmax = c(4.15, 7.55, 10.95),
  ymin = 5.43,
  ymax = 6.12,
  label = c(
    "Anchor boundaries match\nthe supplied 15-food table",
    "Anchor direction is explicit\nfor every one-sided interval",
    "Anchor limits are within\nobserved training support"
  )
)

workflow_headers <- tibble(
  x = c(3.05, 8.95, 6.00, 3.05, 8.95, 6.00),
  y = c(8.95, 8.95, 6.57, 4.00, 4.00, 1.38),
  label = c(
    "01  FROZEN TRAINING COHORT",
    "02  PRESPECIFIED 15-COMPONENT ANCHORS",
    "03  ANCHOR QC — ALL THREE REQUIRED",
    "QC PASSED",
    "QC FAILED",
    "04  APPLY THE FIXED BINARY SCORING RULES"
  ),
  colour = c(
    "#24415B", "#176B63", "#263746",
    "#176B63", "#475569", "#24415B"
  )
)

workflow_bodies <- tibble(
  x = c(3.05, 8.95, 3.05, 8.95, 6.00),
  y = c(8.23, 8.20, 3.30, 3.30, 0.78),
  label = c(
    paste0(
      "Use the exact Step_03_1_0 70% training cohort\n",
      "No re-splitting; no use of the held-out test set"
    ),
    paste0(
      "Apply the supplied value1/value2 for all 15 foods\n",
      "Use explicit low/high directions when value2 is NA\n",
      "Do not move boundaries to chase statistical significance"
    ),
    paste0(
      "Fit the corrected-reference RCS and legacy-reference RCS\n",
      "Report reference, P10/P50/P90, knots, shape and plateaus\n",
      "Compare the anchors with the previous V5 ranges"
    ),
    paste0(
      "Stop if a boundary is outside training support\n",
      "Stop if a direction/bound pair is invalid\n",
      "Stop if any component has no 0/1 variation"
    ),
    paste0(
      "Score 1 = intake satisfies the fixed anchor   |   ",
      "Score 0 = intake outside   |   NA = missing\n",
      "Cox/RCS are diagnostic; anchor boundaries are not data-tuned"
    )
  )
)

p_interval_algorithm <- ggplot() +
  geom_segment(
    data = workflow_connectors,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "#526170",
    linewidth = 0.75,
    arrow = grid::arrow(
      type = "closed",
      length = grid::unit(0.16, "inches")
    )
  ) +
  geom_rect(
    data = workflow_nodes,
    aes(
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax,
      fill = fill, colour = border
    ),
    linewidth = 0.75
  ) +
  geom_rect(
    data = workflow_criteria,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "white",
    colour = "#CBD5E1",
    linewidth = 0.45
  ) +
  geom_text(
    data = workflow_criteria,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3.0,
    lineheight = 1.05,
    colour = "#263746",
    fontface = "bold"
  ) +
  geom_text(
    data = workflow_headers,
    aes(x = x, y = y, label = label),
    colour = workflow_headers$colour,
    size = 3.35,
    fontface = "bold"
  ) +
  geom_text(
    data = workflow_bodies,
    aes(x = x, y = y, label = label),
    colour = "#263746",
    size = 2.85,
    lineheight = 1.12
  ) +
  annotate(
    "text", x = 3.05, y = 4.72,
    label = "YES", size = 3.0, fontface = "bold",
    colour = "#176B63"
  ) +
  annotate(
    "text", x = 8.95, y = 4.72,
    label = "NO", size = 3.0, fontface = "bold",
    colour = "#475569"
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(
    xlim = c(0, 12),
    ylim = c(0, 10),
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    title = "prespecified component anchors and binary scoring",
    caption = paste0(
      "All RCS diagnostics use the exact Step_03_1_0 frozen 70% training cohort. ",
      "Primary scoring anchors are supplied in advance; the 30% test set is untouched."
    )
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 16,
      colour = "#1F2937", margin = margin(b = 12)
    ),
    plot.caption = element_text(
      hjust = 0.5, size = 9, colour = "#475569",
      margin = margin(t = 10)
    ),
    plot.margin = margin(16, 16, 12, 16)
  )

ggsave(
  filename = file.path(PLOT_DIR, "Figure_Sx_interval_derivation_schema.pdf"),
  plot = p_interval_algorithm,
  width = 11.5,
  height = 7.2,
  device = "pdf"
)
ggsave(
  filename = file.path(PLOT_DIR, "Figure_Sx_interval_derivation_schema.png"),
  plot = p_interval_algorithm,
  width = 11.5,
  height = 7.2,
  dpi = 600,
  bg = "white"
)

# Export the compact rules after the full interval table has been created.
write.csv(
  binary_scoring_rules,
  file.path(OUT_DIR, "RWDID_binary_scoring_rules_15_components.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

message("[6B/7] Interval derivation table, binary scoring files, and figures exported.")

# -----------------------------------------------------------------------------
# 7. Final QC
# -----------------------------------------------------------------------------

expected_prediction_files <- file.path(
  PRED_DIR,
  paste0(FOOD_MAP$food, "_rcs3.csv")
)
expected_legacy_prediction_files <- file.path(
  PRED_DIR,
  paste0(FOOD_MAP$food, "_rcs3_legacy_min_reference.csv")
)
expected_plot_files <- file.path(
  PLOT_DIR, paste0(FOOD_MAP$food, "_rcs3.pdf")
)

final_qc <- tibble(
  check = c(
    "training_RDS_matches_frozen_Step03_manifest",
    "test_participants_excluded",
    "exactly_15_Model2_selected_components",
    "green_leafy_veg_excluded",
    "high_fat_dairy_excluded",
    "fruit_juices_included",
    "water_included",
    "one_rotating_lowest_risk_reference_per_component",
    "all_15_RCS_models_fitted",
    "all_15_reference_HR_equal_1",
    "all_15_P10_P50_P90_values_exported",
    "all_15_three_distinct_knots",
    "all_15_primary_category_ranges_exported",
    "all_15_prediction_CSVs_written",
    "all_15_legacy_reference_prediction_CSVs_written",
    "combined_legacy_reference_predictions_cover_all_15_foods",
    "legacy_interval_validation_status_complete_for_all_15_foods",
    "all_15_individual_PDFs_written",
    "no_RCS_refinement_without_statistical_support",
    "all_15_RCS_plots_show_P10_P50_P90_reference_and_legend",
    "percentile_plot_annotation_audit_written",
    "combined_15_component_PDF_written_with_shared_percentile_legend",
    "interval_range_PDF_written",
    "tile_range_PDF_written",
    "15_food_interval_derivation_table_written",
    "15_food_binary_scoring_rules_written",
    "participant_level_training_binary_scores_written",
    "component_0_1_distribution_audit_written",
    "complete_15_component_total_score_distribution_written",
    "interval_derivation_algorithm_PDF_written",
    "interval_derivation_algorithm_PNG_written",
    "binary_scoring_scheme_PDF_written",
    "binary_scoring_scheme_PNG_written",
    "single_value_intervals_shown_as_exact_markers",
    "Word_interval_derivation_table_written",
    "all_15_binary_scoring_rules_defined",
    "all_15_interval_derivation_rows_reconciled",
    "all_15_final_ranges_exactly_match_anchor_table",
    "all_15_components_have_both_binary_score_values",
    "anchor_vs_V5_derived_comparison_written"
  ),
  status = c(
    TRUE,
    TRUE,
    nrow(FOOD_MAP) == 15L,
    !("green_leafy_veg" %in% FOOD_MAP$food),
    !("high_fat_dairy" %in% FOOD_MAP$food),
    "fruit_juices" %in% FOOD_MAP$food,
    "water" %in% FOOD_MAP$food,
    nrow(rcs_reference_audit) == 15L,
    nrow(rcs_model_summary) == 15L,
    all(abs(rcs_model_summary$reference_HR - 1) <= 1e-6),
    all(is.finite(rcs_model_summary$p10) &
          is.finite(rcs_model_summary$p50) &
          is.finite(rcs_model_summary$p90)),
    all(is.finite(rcs_model_summary$rcs_knot_1) &
          is.finite(rcs_model_summary$rcs_knot_2) &
          is.finite(rcs_model_summary$rcs_knot_3) &
          rcs_model_summary$rcs_knot_1 < rcs_model_summary$rcs_knot_2 &
          rcs_model_summary$rcs_knot_2 < rcs_model_summary$rcs_knot_3),
    nrow(food_rcs_range_all) == 15L &&
      all(is.finite(food_rcs_range_all$value1)),
    all(file.exists(expected_prediction_files)),
    all(file.exists(expected_legacy_prediction_files)),
    nrow(legacy_predictions_all) > 0L &&
      length(unique(legacy_predictions_all$food)) == 15L,
    nrow(food_rcs_range_all) == 15L &&
      all(!is.na(food_rcs_range_all$legacy_validation_status)) &&
      all(!is.na(food_rcs_range_all$legacy_interval_retained)),
    all(file.exists(expected_plot_files)),
    all(
      !food_rcs_range_all$refined_by_RCS_5pct |
        rcs_model_summary$rcs_range_eligible
    ),
    nrow(rcs_percentile_plot_audit) == 15L &&
      all(rcs_percentile_plot_audit$percentile_lines_dashed) &&
      all(rcs_percentile_plot_audit$numeric_labels_at_top) &&
      all(rcs_percentile_plot_audit$reference_vertical_line_shown) &&
      all(rcs_percentile_plot_audit$reference_numeric_label_shown) &&
      all(rcs_percentile_plot_audit$legend_included) &&
      all(rcs_percentile_plot_audit$all_percentiles_within_displayed_x_axis),
    file.exists(file.path(AUDIT_DIR, "RCS_percentile_plot_annotation_QC.csv")),
    file.exists(file.path(PLOT_DIR, "rcs3_all_food_combined.pdf")),
    file.exists(file.path(PLOT_DIR, "p_food_group.pdf")),
    file.exists(file.path(PLOT_DIR, "p_food_group_tiles.pdf")),
    file.exists(file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.csv")) &&
      nrow(interval_derivation_table) == 15L,
    file.exists(file.path(OUT_DIR, "RWDID_binary_scoring_rules_15_components.csv")) &&
      nrow(binary_scoring_rules) == 15L,
    file.exists(file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv")) &&
      nrow(score_data) == nrow(ukb_data_depression_rcs),
    file.exists(file.path(AUDIT_DIR, "RWDID_binary_score_distribution_training.csv")) &&
      nrow(binary_score_distribution) == 15L,
    file.exists(file.path(AUDIT_DIR, "RWDID_total_score_distribution_training.csv")) &&
      nrow(total_score_distribution) == 16L,
    file.exists(file.path(PLOT_DIR, "Figure_Sx_interval_derivation_schema.pdf")),
    file.exists(file.path(PLOT_DIR, "Figure_Sx_interval_derivation_schema.png")),
    file.exists(file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.pdf")),
    file.exists(file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.png")),
    nrow(point_marker_data) ==
      sum(food_rcs_range_all$range_direction == "point"),
    file.exists(file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.docx")),
    all(!is.na(food_rcs_range_all$score_rule) &
          food_rcs_range_all$score_rule != "not_defined"),
    nrow(interval_derivation_table) == 15L &&
      identical(as.character(interval_derivation_table$food),
                as.character(FOOD_MAP$food)),
    all(food_rcs_range_all$anchor_override_applied) &&
      identical(as.numeric(food_rcs_range_all$value1), as.numeric(V51_RANGE_ANCHORS$value1)) &&
      identical(as.character(food_rcs_range_all$range_direction), as.character(food_rcs_range_all$anchor_range_direction)) &&
      all((is.na(food_rcs_range_all$value2) & is.na(V51_RANGE_ANCHORS$value2)) |
          (!is.na(food_rcs_range_all$value2) & !is.na(V51_RANGE_ANCHORS$value2) &
           abs(food_rcs_range_all$value2 - V51_RANGE_ANCHORS$value2) <= 1e-12)),
    nrow(binary_score_distribution) == 15L &&
      all(binary_score_distribution$n_score_0 > 0L) &&
      all(binary_score_distribution$n_score_1 > 0L),
    file.exists(file.path(AUDIT_DIR, "user_anchor_vs_V5_derived_interval_comparison.csv"))
  )
)

write.csv(
  final_qc,
  file.path(AUDIT_DIR, "RCS_final_QC.csv"),
  row.names = FALSE
)

if (any(!final_qc$status)) {
  print(final_qc)
  stop("STOP: One or more final RCS QC checks failed.")
}

cat("\n===== RCS MODEL AND INTERVAL SUMMARY =====\n")
print(
  rcs_model_summary %>%
    select(
      group,
      reference_category_label,
      reference_median_used,
      category_range_low,
      category_range_high,
      category_range_direction,
      p_overall,
      p_non_linear,
      curve_shape,
      plateau_05_low,
      plateau_05_high,
      plateau_05_overlaps_reference_category,
      refined_by_RCS_5pct,
      primary_range_value1,
      primary_range_value2,
      primary_range_direction,
      v5_derived_value1,
      v5_derived_value2,
      v5_derived_range_direction,
      legacy_value1,
      legacy_value2,
      legacy_range_direction,
      legacy_overlaps_category,
      legacy_overlaps_5pct_plateau,
      legacy_high_tail_bridge,
      legacy_interval_retained,
      legacy_validation_status,
      score_rule
    )
)

cat("\n============================================================\n")
cat("Step_04_1_1 RCS USER-ANCHORED RANGES COMPLETED\n")
cat("Training cohort: exact Step_03_1_0 frozen 70% RDS; no re-split.\n")
cat("Frozen 15-component set retained, including Water; no re-split or test-set fitting.\n")
cat("P10/P50/P90: coloured dashed markers and top-row numeric labels on all RCS plots.\n")
cat("Percentile colours: P10 #77AADD; P50 #009988; P90 #9467BD.\n")
cat("Combined RCS figure uses one shared percentile legend.\n")
cat("Binary scoring figure: 20% colour tints, enlarged bold text, near-square canvas.\n")
cat("Supplementary table: CSV audit plus formatted Word version with detailed footnotes.\n")
cat("Primary intervals: fixed to the supplied 15-component anchor table; Cox/RCS results are diagnostics, not boundary selectors.\n")
cat("One-sided intervals use explicit low/high directions recorded in user_prespecified_range_anchors.csv.\n")
cat("RCS 5%/10%/20% plateaus and V5.0-derived intervals are comparison diagnostics only.\n")
cat("Legacy-reference curves are audit-only; Figure S4 uses the lowest-risk-category median as HR=1.\n")
cat("Missing component intake remains NA; no imputation.\n")
cat("Continuous 0-1 scoring sensitivity analysis is not run in this script.\n")
cat("Word table: ", file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.docx"), "\n", sep = "")
cat("Interval CSV: ", file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.csv"), "\n", sep = "")
cat("Algorithm figure: ", file.path(PLOT_DIR, "Figure_Sx_interval_derivation_schema.pdf"), "\n", sep = "")
cat("Binary scoring figure: ", file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.pdf"), "\n", sep = "")
cat("Training scores: ", file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv"), "\n", sep = "")
cat("Anchor comparison: ", file.path(AUDIT_DIR, "user_anchor_vs_V5_derived_interval_comparison.csv"), "\n", sep = "")
cat("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 05/11: Step_04_1_1_1_RCS_BASE
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 06/11: Step_04_1_1_2_INTERVAL_RECONCILIATION
# Source: Step_04_1_1_2_data_derived_interval_reconciliation.R
# Provenance: project source archive data_analysis_v2.zip
# SHA256 of included block: 67d3206aedbaea68d885afb8aaa0ccd3dd25eb4764850e8e03f6d93a8946c569
########################################################################################

# =============================================================================
# Step_04_1_1_2_DATA_DERIVED_INTERVAL_RECONCILIATION.R
#
# Purpose:
#   Reconcile the six reviewer-critical component intervals using the FROZEN
#   RCS predictions and the frozen lowest-risk reference-rotation results.
#   This script does NOT refit the RCS models or alter the frozen 70/30 split.
#
#   It promotes the current-data RCS 5% minimum-risk plateau to the primary
#   interval evidence, then cross-checks its boundaries against Model 2
#   lowest-risk Cox categories and their pairwise contrasts.
#
#   The other nine components retain their frozen scoring intervals.
#   The six revised component scores are recalculated from the frozen training
#   RDS. The resulting files are designed for Beta V5.2 downstream validation.
#
# Key rules:
#   - U-shaped: start from the corrected-reference 5% minimum-risk plateau.
#     If a significantly higher-risk adjacent Cox category overlaps the
#     plateau, truncate the plateau at that category boundary.
#   - Processed meats (reverse-J / decreasing-risk): use the lower edge of the
#     5% plateau as a high-intake scoring threshold. Report the RCS plateau
#     support and the wider highest-risk Cox category separately.
#   - Sweetened beverages (J-shaped / increasing-risk): start at the lowest-
#     risk zero-intake category and extend only through contiguous adjacent
#     categories whose pairwise 95% CI includes HR=1, capped by the 5% plateau.
#   - These are operational scoring intervals / minimum-risk regions, NOT
#     automatically statistically significant protective ranges or advice.
#
# Frozen inputs:
#   ./results/Step_04/Step_04_1_1_1
#   ./results/Step_03/reference_rotation_lowest_risk
#   ./results/Step_03/ukb_70_30_split/main/ukb_depression_training_cohort.rds
#
# Outputs:
#   ./results/Step_04/Step_04_1_1_2
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

options(stringsAsFactors = FALSE, scipen = 999)

# -----------------------------------------------------------------------------
# 0. Paths and frozen analysis settings
# -----------------------------------------------------------------------------

V51_DIR <- "./results/Step_04/Step_04_1_1_1"
ROTATION_DIR <- "./results/Step_03/reference_rotation_lowest_risk"
TRAIN_RDS <- "./results/Step_03/ukb_70_30_split/main/ukb_depression_training_cohort.rds"

OUT_DIR <- "./results/Step_04/Step_04_1_1_2"
AUDIT_DIR <- file.path(OUT_DIR, "audit")
PRED_DIR <- file.path(OUT_DIR, "predictions")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PRED_DIR, recursive = TRUE, showWarnings = FALSE)

EXPECTED_TRAIN_N <- 76578L
EXPECTED_TRAIN_EVENTS <- 2122L
ID_COL <- "eid"

TARGET_FOODS <- c(
  "citrus_fruits",
  "fruit_juices",
  "eggs",
  "low_fat_dairy",
  "processed_meats",
  "sweetened_bev"
)

# -----------------------------------------------------------------------------
# 1. Input file checks
# -----------------------------------------------------------------------------

V51_RANGE_FILE <- file.path(V51_DIR, "food_rcs_range_all.csv")
V51_RULE_FILE <- file.path(
  V51_DIR, "RWDID_binary_scoring_rules_15_components.csv"
)
V51_SCORE_FILE <- file.path(
  V51_DIR, "RWDID_binary_component_scores_training.csv"
)
V51_QC_FILE <- file.path(V51_DIR, "audit", "RCS_final_QC.csv")
V51_MODEL_FILE <- file.path(
  V51_DIR, "audit", "RCS_model_summary_and_range_diagnostics.csv"
)
V51_REFERENCE_AUDIT_FILE <- file.path(
  V51_DIR, "audit", "RCS_lowest_risk_reference_category_audit.csv"
)
V51_ANCHOR_FILE <- file.path(
  V51_DIR, "audit", "user_prespecified_range_anchors.csv"
)

ROTATION_SELECTED_FILE <- file.path(
  ROTATION_DIR, "food_group_reference_rotation_selected_training.csv"
)
ROTATION_PAIRWISE_FILE <- file.path(
  ROTATION_DIR, "food_group_reference_rotation_all_comparisons_training.csv"
)

required_files <- c(
  V51_RANGE_FILE, V51_RULE_FILE, V51_SCORE_FILE, V51_QC_FILE,
  V51_MODEL_FILE, V51_REFERENCE_AUDIT_FILE, V51_ANCHOR_FILE,
  ROTATION_SELECTED_FILE, ROTATION_PAIRWISE_FILE, TRAIN_RDS
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "STOP: Required frozen / reference-rotation / training files are missing:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 2. Read and validate frozen outputs
# -----------------------------------------------------------------------------

v51_ranges <- read.csv(
  V51_RANGE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
v51_rules <- read.csv(
  V51_RULE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
v51_scores <- read.csv(
  V51_SCORE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
v51_qc <- read.csv(
  V51_QC_FILE, stringsAsFactors = FALSE, check.names = FALSE
)
rcs_diag <- read.csv(
  V51_MODEL_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
reference_audit <- read.csv(
  V51_REFERENCE_AUDIT_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
target_anchors <- read.csv(
  V51_ANCHOR_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
rotation_selected <- read.csv(
  ROTATION_SELECTED_FILE, stringsAsFactors = FALSE, check.names = FALSE
)
rotation_pairwise <- read.csv(
  ROTATION_PAIRWISE_FILE, stringsAsFactors = FALSE, check.names = FALSE
)

if (!all(as.logical(v51_qc$status))) {
  stop("STOP: Frozen RCS QC is not entirely PASS.", call. = FALSE)
}
if (nrow(v51_ranges) != 15L ||
    nrow(v51_rules) != 15L ||
    nrow(v51_scores) != EXPECTED_TRAIN_N ||
    nrow(rcs_diag) != 15L ||
    nrow(reference_audit) != 15L) {
  stop("STOP: Frozen files do not have the expected 15-component / training dimensions.")
}
if (!all(TARGET_FOODS %in% v51_ranges$food) ||
    !all(TARGET_FOODS %in% rcs_diag$food) ||
    !all(TARGET_FOODS %in% reference_audit$food)) {
  stop("STOP: One or more target foods are absent from the frozen outputs.")
}
if (!all(c("model", "food_variable", "lowest_risk_reference_category") %in%
         names(rotation_selected))) {
  stop("STOP: Reference-rotation selected file has an unexpected schema.")
}
if (!all(c(
  "model", "food_variable", "reference_category", "comparison_category",
  "is_selected_lowest_risk_reference", "is_reference_row",
  "hazard_ratio", "ci_low", "ci_high", "p_value_raw",
  "model_n", "model_events"
) %in% names(rotation_pairwise))) {
  stop("STOP: Reference-rotation pairwise file has an unexpected schema.")
}

# Check that the frozen RCS fit and the independent rotation analysis used the
# same training sample/event count for the six target foods.
if (any(rcs_diag$model_n[match(TARGET_FOODS, rcs_diag$food)] != EXPECTED_TRAIN_N) ||
    any(rcs_diag$model_events[match(TARGET_FOODS, rcs_diag$food)] != EXPECTED_TRAIN_EVENTS)) {
  stop("STOP: RCS target-food sample size/events differ from the frozen training cohort.")
}

# -----------------------------------------------------------------------------
# 3. Helpers for category labels, intervals, and score rules
# -----------------------------------------------------------------------------

parse_category_bounds <- function(label) {
  pattern <- "\\(([0-9]+\\.?[0-9]*|\\.[0-9]+)[[:space:]]*-[[:space:]]*([0-9]+\\.?[0-9]*|\\.[0-9]+)\\)"
  m <- regexec(pattern, as.character(label), perl = TRUE)
  z <- regmatches(as.character(label), m)[[1]]
  if (length(z) != 3L) return(c(low = NA_real_, high = NA_real_))
  c(low = as.numeric(z[2]), high = as.numeric(z[3]))
}

parse_category_median <- function(label) {
  pattern <- "^[[:space:]]*([0-9]+\\.?[0-9]*|\\.[0-9]+)"
  m <- regmatches(
    as.character(label),
    regexec(pattern, as.character(label), perl = TRUE)
  )[[1]]
  if (length(m) < 2L) return(NA_real_)
  as.numeric(m[2])
}

fmt <- function(x, digits = 3L) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  if (abs(x) < 1e-10) return("0")
  formatC(x, format = "f", digits = digits)
}

range_label <- function(direction, value1, value2) {
  if (direction == "middle") {
    return(paste0(fmt(value1), "–", fmt(value2)))
  }
  if (direction == "low") {
    return(paste0("≤", fmt(value1)))
  }
  if (direction == "high") {
    return(paste0("≥", fmt(value1)))
  }
  if (direction == "point") {
    return(paste0("=", fmt(value1)))
  }
  "Undefined"
}

score_one <- function(x, direction, value1, value2, food) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  out <- rep(NA_integer_, length(x))
  observed <- !is.na(x)

  if (any(!is.finite(x[observed]))) {
    stop("STOP: Non-finite non-missing intake values for ", food)
  }
  if (any(x[observed] < -1e-8)) {
    stop("STOP: Negative intake values for ", food)
  }
  if (!any(observed)) return(out)

  scale_value <- max(
    1,
    abs(value1),
    ifelse(is.finite(value2), abs(value2), 0)
  )
  tol <- max(1e-8, scale_value * 1e-8)

  if (direction == "middle") {
    if (!is.finite(value1) || !is.finite(value2) || value2 < value1) {
      stop("STOP: Invalid middle interval for ", food)
    }
    out[observed] <- as.integer(
      x[observed] >= value1 - tol & x[observed] <= value2 + tol
    )
  } else if (direction == "low") {
    if (!is.finite(value1)) stop("STOP: Invalid low interval for ", food)
    out[observed] <- as.integer(x[observed] <= value1 + tol)
  } else if (direction == "high") {
    if (!is.finite(value1)) stop("STOP: Invalid high interval for ", food)
    out[observed] <- as.integer(x[observed] >= value1 - tol)
  } else if (direction == "point") {
    if (!is.finite(value1)) stop("STOP: Invalid point interval for ", food)
    out[observed] <- as.integer(abs(x[observed] - value1) <= tol)
  } else {
    stop("STOP: Unsupported score direction for ", food)
  }
  out
}

get_contiguous_intervals <- function(pred, food, condition, label) {
  pred <- pred[order(pred[[food]]), , drop = FALSE]
  keep <- which(condition(pred))
  if (length(keep) == 0L) {
    return(tibble(evidence = label, interval = "None", lower = NA_real_,
                  upper = NA_real_, n_grid_points = 0L))
  }
  split_id <- cumsum(c(1L, as.integer(diff(keep) > 1L)))
  groups <- split(keep, split_id)
  bind_rows(lapply(groups, function(ii) {
    tibble(
      evidence = label,
      interval = paste0(fmt(pred[[food]][min(ii)]), "–",
                        fmt(pred[[food]][max(ii)])),
      lower = pred[[food]][min(ii)],
      upper = pred[[food]][max(ii)],
      n_grid_points = length(ii)
    )
  }))
}

# -----------------------------------------------------------------------------
# 4. Derive data-based intervals for the six target foods
# -----------------------------------------------------------------------------

derive_interval <- function(food) {
  d <- rcs_diag[rcs_diag$food == food, , drop = FALSE]
  a <- reference_audit[reference_audit$food == food, , drop = FALSE]

  if (nrow(d) != 1L || nrow(a) != 1L) {
    stop("STOP: Expected exactly one diagnostic row for ", food)
  }

  selected <- rotation_selected[
    rotation_selected$food_variable == paste0(food, "_group_range") &
      rotation_selected$model %in% c("model1", "model2"),
    , drop = FALSE
  ]
  if (nrow(selected) != 2L ||
      length(unique(selected$lowest_risk_reference_category)) != 1L) {
    stop(
      "STOP: Model 1 and Model 2 do not agree on the selected lowest-risk ",
      "reference category for ", food,
      ". Review reference rotation before deriving the interval."
    )
  }

  ref_label <- as.character(selected$lowest_risk_reference_category[1])
  if (!identical(ref_label, as.character(d$reference_category_label[1]))) {
    stop(
      "STOP: Frozen RCS reference category does not match the reference-rotation ",
      "result for ", food, ". RCS=", d$reference_category_label[1],
      "; rotation=", ref_label
    )
  }

  ref_bounds <- parse_category_bounds(ref_label)
  ref_median <- parse_category_median(ref_label)
  if (any(!is.finite(ref_bounds)) || !is.finite(ref_median)) {
    stop("STOP: Could not parse the selected reference category for ", food)
  }

  pair <- rotation_pairwise[
    rotation_pairwise$model == "model2" &
      rotation_pairwise$food_variable == paste0(food, "_group_range") &
      rotation_pairwise$reference_category == ref_label &
      as.logical(rotation_pairwise$is_selected_lowest_risk_reference) &
      !as.logical(rotation_pairwise$is_reference_row),
    , drop = FALSE
  ]
  if (nrow(pair) < 1L) {
    stop("STOP: No Model 2 pairwise comparisons found for ", food)
  }

  pair_bounds <- t(vapply(
    pair$comparison_category,
    parse_category_bounds,
    numeric(2)
  ))
  pair$category_low <- pair_bounds[, "low"]
  pair$category_high <- pair_bounds[, "high"]
  pair$category_median <- vapply(
    pair$comparison_category, parse_category_median, numeric(1)
  )
  if (any(!is.finite(pair$category_low)) ||
      any(!is.finite(pair$category_high)) ||
      any(!is.finite(pair$category_median))) {
    stop("STOP: Could not parse comparison category boundaries for ", food)
  }

  plateau_low <- as.numeric(d$plateau_05_low[1])
  plateau_high <- as.numeric(d$plateau_05_high[1])
  support_low <- as.numeric(d$prediction_support_low[1])
  support_high <- as.numeric(d$prediction_support_high[1])
  if (!all(is.finite(c(plateau_low, plateau_high, support_low, support_high))) ||
      plateau_high < plateau_low ||
      plateau_low < support_low - 1e-8 ||
      plateau_high > support_high + 1e-8) {
    stop("STOP: Invalid or unsupported 5% RCS plateau for ", food)
  }

  direction <- "middle"
  value1 <- plateau_low
  value2 <- plateau_high
  rule_basis <- "Corrected-reference RCS 5% minimum-risk plateau"
  cox_boundary_note <- "Cox lowest-risk category used as cross-check"

  if (food %in% c(
    "citrus_fruits", "fruit_juices", "eggs", "low_fat_dairy"
  )) {
    # U-shaped components: start with the RCS 5% plateau. If a significantly
    # higher-risk Cox category overlaps a plateau edge, trim only that edge.
    tol <- 1e-8

    lower_sig <- pair[
      pair$category_median < ref_median &
        pair$ci_low > 1 &
        pair$category_high >= plateau_low - tol &
        pair$category_high <= ref_bounds["low"] + tol,
      , drop = FALSE
    ]
    if (nrow(lower_sig) > 0L) {
      value1 <- max(value1, ref_bounds["low"])
      cox_boundary_note <- paste0(
        cox_boundary_note,
        "; lower boundary checked against significantly higher-risk lower category"
      )
    }

    upper_sig <- pair[
      pair$category_median > ref_median &
        pair$ci_low > 1 &
        pair$category_low <= plateau_high + tol &
        pair$category_low >= ref_bounds["high"] - tol,
      , drop = FALSE
    ]
    if (nrow(upper_sig) > 0L) {
      first_sig_upper <- upper_sig[which.min(upper_sig$category_low), , drop = FALSE]
      value2 <- min(value2, first_sig_upper$category_low[1])
      cox_boundary_note <- paste0(
        cox_boundary_note,
        "; upper boundary capped at first significantly higher-risk Cox category"
      )
    }

    if (value2 <= value1) {
      stop("STOP: Cox boundary checks collapsed the interval for ", food)
    }
    direction <- "middle"
    rule_basis <- paste(
      "5% RCS minimum-risk plateau;",
      "Cox category boundary used only when an adjacent higher-risk category",
      "significantly overlapped the plateau edge"
    )

  } else if (food == "processed_meats") {
    # The selected Cox reference is the highest intake category. The lower edge
    # is the larger of the 5% plateau lower boundary and Cox category lower edge.
    # Scoring remains high-sided; the RCS plateau high edge is retained as an
    # evidence-support field and is NOT hidden as if RCS supported infinity.
    value1 <- max(plateau_low, ref_bounds["low"])
    value2 <- NA_real_
    direction <- "high"
    rule_basis <- paste(
      "High-intake direction supported by the lowest-risk Cox category;",
      "threshold set at the lower edge of the corrected-reference RCS 5% plateau"
    )
    cox_boundary_note <- paste0(
      "Lowest-risk Cox category = ", ref_label,
      "; categorical support extends to ", fmt(ref_bounds["high"]),
      ", while the continuous RCS 5% plateau ends at ", fmt(plateau_high)
    )

  } else if (food == "sweetened_bev") {
    # J-shaped/increasing-risk component: start at the zero-intake reference
    # category and extend only through contiguous categories whose pairwise
    # 95% CI includes HR=1, with the extension capped by the 5% RCS plateau.
    current_upper <- min(plateau_high, ref_bounds["high"])
    if (ref_bounds["low"] > support_low + 1e-8) {
      stop("STOP: Sweetened-beverage reference category is not at the low boundary.")
    }

    ordered_pair <- pair[order(pair$category_low, pair$category_median), , drop = FALSE]
    for (j in seq_len(nrow(ordered_pair))) {
      z <- ordered_pair[j, , drop = FALSE]
      contiguous <- z$category_low[1] <= current_upper + 1e-8 &&
        z$category_high[1] > current_upper + 1e-8
      if (!contiguous) next

      ci_includes_one <- z$ci_low[1] <= 1 && z$ci_high[1] >= 1
      if (!ci_includes_one) {
        # Stop at the first adjacent category whose 95% CI excludes HR=1.
        break
      }
      current_upper <- min(plateau_high, z$category_high[1])
      if (current_upper >= plateau_high - 1e-8) break
    }

    value1 <- ref_bounds["low"]
    value2 <- current_upper
    direction <- "middle"
    if (!is.finite(value2) || value2 <= value1) {
      stop("STOP: Could not derive a non-zero sweetened-beverage interval.")
    }
    rule_basis <- paste(
      "Zero-intake lowest-risk Cox category extended through contiguous",
      "non-significantly different category/categories, capped by the 5% RCS plateau"
    )
    cox_boundary_note <- paste0(
      "0–1 category CI includes HR=1; the next 1–2 category has lower 95% CI >1"
    )
  }

  # The resulting data-derived interval must remain inside the frozen exposure
  # support used by the corrected-reference RCS predictions.
  if (value1 < support_low - 1e-8 ||
      (direction == "middle" && value2 > support_high + 1e-8)) {
    stop("STOP: Derived interval falls outside RCS prediction support for ", food)
  }

  tibble(
    food = food,
    display_name = as.character(d$display_name[1]),
    reference_category = ref_label,
    reference_median = ref_median,
    reference_category_low = ref_bounds["low"],
    reference_category_high = ref_bounds["high"],
    curve_shape = as.character(d$curve_shape[1]),
    p_overall = as.numeric(d$p_overall[1]),
    p_non_linear = as.numeric(d$p_non_linear[1]),
    plateau_05_low = plateau_low,
    plateau_05_high = plateau_high,
    prediction_support_low = support_low,
    prediction_support_high = support_high,
    derived_value1 = as.numeric(value1),
    derived_value2 = as.numeric(value2),
    derived_direction = direction,
    derived_range = range_label(direction, value1, value2),
    derivation_basis = rule_basis,
    cox_crosscheck = cox_boundary_note
  )
}

derived_six <- bind_rows(lapply(TARGET_FOODS, derive_interval))

# -----------------------------------------------------------------------------
# 5. Extract pointwise significant-benefit / increased-risk regions
# -----------------------------------------------------------------------------

pointwise_evidence <- bind_rows(lapply(TARGET_FOODS, function(food) {
  pred_file <- file.path(V51_DIR, "predictions", paste0(food, "_rcs3.csv"))
  if (!file.exists(pred_file)) {
    stop("STOP: Missing frozen corrected-reference prediction file: ", pred_file)
  }
  pred <- read.csv(
    pred_file, stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c(food, "yhat", "lower", "upper")
  if (!all(required %in% names(pred))) {
    stop("STOP: Prediction file lacks required columns for ", food)
  }
  support_row <- rcs_diag[rcs_diag$food == food, , drop = FALSE]
  support_low <- as.numeric(support_row$prediction_support_low[1])
  support_high <- as.numeric(support_row$prediction_support_high[1])
  pred <- pred[
    pred[[food]] >= support_low - 1e-10 &
      pred[[food]] <= support_high + 1e-10,
    , drop = FALSE
  ]
  pred <- pred[order(pred[[food]]), , drop = FALSE]

  beneficial <- get_contiguous_intervals(
    pred, food, function(z) z$upper < 1,
    "Pointwise statistically significant beneficial region (upper 95% CI < 1)"
  )
  increased <- get_contiguous_intervals(
    pred, food, function(z) z$lower > 1,
    "Pointwise statistically significant increased-risk region (lower 95% CI > 1)"
  )
  bind_rows(
    beneficial %>% mutate(food = food),
    increased %>% mutate(food = food)
  ) %>% select(food, everything())
}))

# -----------------------------------------------------------------------------
# 6. Compare target intervals with current-data-derived intervals
# -----------------------------------------------------------------------------

target_direction_map <- target_anchors$range_direction
names(target_direction_map) <- target_anchors$food
# The processed-meat target was explicitly high-sided in the frozen
# derived output; resolve the anchor's "v5_derived" label from that output.
for (food in TARGET_FOODS) {
  if (identical(target_direction_map[[food]], "v5_derived")) {
    target_direction_map[[food]] <- as.character(
      v51_ranges$range_direction[match(food, v51_ranges$food)]
    )
  }
}

target_comparison <- derived_six %>%
  left_join(
    target_anchors %>%
      filter(food %in% TARGET_FOODS) %>%
      transmute(
        food,
        target_value1 = as.numeric(value1),
        target_value2 = as.numeric(value2),
        target_direction_raw = as.character(range_direction)
      ),
    by = "food"
  ) %>%
  mutate(
    target_direction = vapply(
      food,
      function(f) as.character(target_direction_map[[f]]),
      character(1)
    ),
    target_range = mapply(
      range_label, target_direction, target_value1, target_value2,
      USE.NAMES = FALSE
    ),
    final_judgement = case_when(
      food == "citrus_fruits" ~
        "Partial compatibility: target lower bound is below the current 5% plateau; the data-derived range is 0.558–2.000.",
      food == "fruit_juices" ~
        "Partial compatibility: the target includes zero/very-low intake below the 5% plateau; the selected Cox category lies within the plateau.",
      food == "eggs" ~
        "Target is too broad at the upper end: the 1–10 intake Cox category is significantly higher-risk than the lowest-risk category; upper boundary capped at 1.000.",
      food == "low_fat_dairy" ~
        "Direction mismatch: the target x<=3.071574 includes high intake, whereas the current RCS is U-shaped and the middle Cox category is lowest-risk; use a bounded middle interval.",
      food == "processed_meats" ~
        "Direction concordant: target threshold 2.700 lies inside the 2.271–3.000 RCS plateau; open upper scoring direction is supported categorically to the highest Cox category, but continuous RCS support ends at 3.000.",
      food == "sweetened_bev" ~
        "Target is too broad: the 1–2 and 2–30 categories have significantly higher risk than zero intake; the adjacent 0–1 category is not significantly different and lies within the 5% plateau.",
      TRUE ~ "Review component-specific diagnostics."
    )
  ) %>%
  select(
    food, display_name, target_range, derived_range,
    curve_shape, p_overall, p_non_linear,
    reference_category, plateau_05_low, plateau_05_high,
    derived_value1, derived_value2, derived_direction,
    derivation_basis, cox_crosscheck, final_judgement
  )

# -----------------------------------------------------------------------------
# 7. Load frozen training cohort and recalculate six participant-level scores
# -----------------------------------------------------------------------------

training_data <- readRDS(TRAIN_RDS)
if (!is.data.frame(training_data)) {
  stop("STOP: Frozen training RDS did not contain a data.frame/tibble.")
}
if (!all(c(ID_COL, TARGET_FOODS) %in% names(training_data))) {
  stop(
    "STOP: Training RDS lacks ID or target food intake columns: ",
    paste(setdiff(c(ID_COL, TARGET_FOODS), names(training_data)), collapse = ", ")
  )
}
if (nrow(training_data) != EXPECTED_TRAIN_N ||
    anyDuplicated(as.character(training_data[[ID_COL]])) > 0L) {
  stop("STOP: Frozen training RDS dimensions or participant IDs are invalid.")
}

normalize_id <- function(x) {
  z <- trimws(as.character(x))
  z[z == ""] <- NA_character_
  sub("\\.0+$", "", z)
}

training_ids <- normalize_id(training_data[[ID_COL]])
score_ids <- normalize_id(v51_scores[[ID_COL]])

if (anyNA(training_ids) || anyNA(score_ids) ||
    anyDuplicated(training_ids) > 0L ||
    anyDuplicated(score_ids) > 0L ||
    !setequal(training_ids, score_ids)) {
  stop("STOP: participant score IDs do not exactly match the frozen training RDS.")
}

score_data <- v51_scores
score_data[[ID_COL]] <- score_ids
training_match <- match(score_ids, training_ids)
if (anyNA(training_match)) {
  stop("STOP: Could not align frozen food intake values to participant scores.")
}

score_col_map <- c(
  citrus_fruits = "score_citrus_fruits",
  fruit_juices = "score_fruit_juices",
  eggs = "score_eggs",
  low_fat_dairy = "score_low_fat_dairy",
  processed_meats = "score_processed_meats",
  sweetened_bev = "score_sweetened_bev"
)

for (food in TARGET_FOODS) {
  r <- derived_six[derived_six$food == food, , drop = FALSE]
  x <- suppressWarnings(as.numeric(training_data[[food]][training_match]))
  score_data[[score_col_map[[food]]]] <- score_one(
    x = x,
    direction = as.character(r$derived_direction[1]),
    value1 = as.numeric(r$derived_value1[1]),
    value2 = as.numeric(r$derived_value2[1]),
    food = food
  )
}

score_cols <- c(
  "score_citrus_fruits",
  "score_red_orange_veg",
  "score_other_veg",
  "score_potatoes",
  "score_refined_grains",
  "score_eggs",
  "score_low_fat_dairy",
  "score_red_meats",
  "score_processed_meats",
  "score_fish_seafood",
  "score_coffee",
  "score_fruit_juices",
  "score_sweetened_bev",
  "score_water",
  "score_wine"
)
if (!all(score_cols %in% names(score_data))) {
  stop("STOP: One or more of the frozen 15 component score columns are missing.")
}

score_matrix <- as.matrix(score_data[, score_cols, drop = FALSE])
storage.mode(score_matrix) <- "numeric"
score_data$n_components_scored <- as.integer(rowSums(!is.na(score_matrix)))
score_data$RWDID_total_score_complete15 <- ifelse(
  score_data$n_components_scored == 15L,
  as.integer(rowSums(score_matrix, na.rm = TRUE)),
  NA_integer_
)
score_data$analysis_split <- "training"

# -----------------------------------------------------------------------------
# 8. Update full 15-component range and scoring-rule exports
# -----------------------------------------------------------------------------

new_ranges <- v51_ranges
new_rules <- v51_rules

for (food in TARGET_FOODS) {
  r <- derived_six[derived_six$food == food, , drop = FALSE]
  idx <- match(food, new_ranges$food)
  ridx <- match(food, new_rules$food)

  direction <- as.character(r$derived_direction[1])
  value1 <- as.numeric(r$derived_value1[1])
  value2 <- as.numeric(r$derived_value2[1])

  score_lower <- if (direction == "low") {
    min(as.numeric(training_data[[food]]), na.rm = TRUE)
  } else {
    value1
  }
  score_upper <- if (direction == "high") Inf else value2
  rule <- if (direction == "middle") {
    paste0(fmt(value1), " <= x <= ", fmt(value2))
  } else if (direction == "high") {
    paste0("x >= ", fmt(value1))
  } else if (direction == "low") {
    paste0("x <= ", fmt(value1))
  } else {
    paste0("x = ", fmt(value1))
  }

  # Update the range file.
  new_ranges$value1[idx] <- value1
  new_ranges$value2[idx] <- value2
  new_ranges$range_direction[idx] <- direction
  new_ranges$range_source[idx] <- "V5.2_data_derived_RCS_5pct_plus_Model2_Cox_crosscheck"
  new_ranges$score_lower[idx] <- score_lower
  new_ranges$score_upper[idx] <- score_upper
  new_ranges$score_rule[idx] <- rule
  if ("anchor_override_applied" %in% names(new_ranges)) {
    new_ranges$anchor_override_applied[idx] <- FALSE
  }

  # Retain the original user target as a comparison field; it is not used
  # to determine the V5.2 primary interval.
  if (!"data_derived_value1" %in% names(new_ranges)) {
    new_ranges$data_derived_value1 <- NA_real_
    new_ranges$data_derived_value2 <- NA_real_
    new_ranges$data_derived_range_direction <- NA_character_
    new_ranges$data_derived_range_source <- NA_character_
  }
  new_ranges$data_derived_value1[idx] <- value1
  new_ranges$data_derived_value2[idx] <- value2
  new_ranges$data_derived_range_direction[idx] <- direction
  new_ranges$data_derived_range_source[idx] <- "RCS_5pct_plateau_with_Model2_Cox_crosscheck"

  # Update the downstream compact scoring-rule file.
  new_rules$value1[ridx] <- value1
  new_rules$value2[ridx] <- value2
  new_rules$range_direction[ridx] <- direction
  new_rules$score_lower[ridx] <- score_lower
  new_rules$score_upper[ridx] <- score_upper
  new_rules$score_rule[ridx] <- rule
  new_rules$score_1[ridx] <- paste0("1 if ", rule)
  new_rules$score_0[ridx] <- "0 otherwise"
  new_rules$missing_rule[ridx] <- "NA; no imputation"

  if ("interval_source" %in% names(new_rules)) {
    new_rules$interval_source[ridx] <- "V5.2 data-derived RCS/Cox interval"
  }
  if ("range_source" %in% names(new_rules)) {
    new_rules$range_source[ridx] <- "V5.2_data_derived_RCS_5pct_plus_Model2_Cox_crosscheck"
  }
  if (!"data_derived_value1" %in% names(new_rules)) {
    new_rules$data_derived_value1 <- NA_real_
    new_rules$data_derived_value2 <- NA_real_
    new_rules$data_derived_range_direction <- NA_character_
  }
  new_rules$data_derived_value1[ridx] <- value1
  new_rules$data_derived_value2[ridx] <- value2
  new_rules$data_derived_range_direction[ridx] <- direction
}

# -----------------------------------------------------------------------------
# 9. Distribution and reconciliation audits
# -----------------------------------------------------------------------------

binary_distribution <- bind_rows(lapply(score_cols, function(sc) {
  z <- score_data[[sc]]
  tibble(
    score_column = sc,
    n_training = length(z),
    n_score_0 = sum(z == 0, na.rm = TRUE),
    n_score_1 = sum(z == 1, na.rm = TRUE),
    n_missing = sum(is.na(z)),
    pct_score_1 = 100 * mean(z == 1, na.rm = TRUE)
  )
}))

total_distribution <- as.data.frame(table(
  factor(score_data$RWDID_total_score_complete15, levels = 0:15),
  useNA = "no"
))
names(total_distribution) <- c("total_score", "n")
total_distribution$total_score <- as.integer(as.character(total_distribution$total_score))
total_distribution$pct <- 100 * total_distribution$n / sum(total_distribution$n)

# Keep a focused pairwise table for the six target foods.
pairwise_evidence <- rotation_pairwise %>%
  filter(
    model == "model2",
    food_variable %in% paste0(TARGET_FOODS, "_group_range"),
    as.logical(is_selected_lowest_risk_reference)
  ) %>%
  mutate(
    comparison_category_low = vapply(
      comparison_category, function(z) parse_category_bounds(z)["low"], numeric(1)
    ),
    comparison_category_high = vapply(
      comparison_category, function(z) parse_category_bounds(z)["high"], numeric(1)
    ),
    significantly_higher_than_reference = ci_low > 1,
    significantly_lower_than_reference = ci_high < 1,
    CI_includes_1 = ci_low <= 1 & ci_high >= 1
  )

# Target vs derived comparison includes the exact user target values and
# the evidence-derived final values.
target_comparison$derived_range <- mapply(
  range_label,
  derived_six$derived_direction[match(target_comparison$food, derived_six$food)],
  derived_six$derived_value1[match(target_comparison$food, derived_six$food)],
  derived_six$derived_value2[match(target_comparison$food, derived_six$food)],
  USE.NAMES = FALSE
)

# -----------------------------------------------------------------------------
# 10. Final QC gates
# -----------------------------------------------------------------------------

qc <- tibble(
  check = c(
    "frozen_QC_all_pass",
    "exact_training_ID_match",
    "training_N_matches_76578",
    "all_15_components_retained",
    "six_target_foods_derived",
    "Model1_Model2_lowest_reference_agreement_for_targets",
    "RCS_reference_matches_reference_rotation_for_targets",
    "all_six_derived_ranges_valid",
    "all_six_derived_ranges_within_RCS_support",
    "processed_meats_high_threshold_inside_RCS_plateau",
    "sweetened_beverages_cutoff_uses_adjacent_Cox_comparison",
    "all_15_scores_binary_or_missing",
    "all_15_components_have_both_0_and_1",
    "participant_score_file_has_expected_training_rows",
    "range_and_rule_files_reconcile"
  ),
  status = c(
    all(as.logical(v51_qc$status)),
    setequal(training_ids, score_ids),
    nrow(score_data) == EXPECTED_TRAIN_N,
    length(score_cols) == 15L,
    nrow(derived_six) == 6L && setequal(derived_six$food, TARGET_FOODS),
    all(vapply(TARGET_FOODS, function(food) {
      z <- rotation_selected[
        rotation_selected$food_variable == paste0(food, "_group_range") &
          rotation_selected$model %in% c("model1", "model2"),
        "lowest_risk_reference_category"
      ]
      length(z) == 2L && length(unique(z)) == 1L
    }, logical(1))),
    all(vapply(TARGET_FOODS, function(food) {
      z <- rotation_selected[
        rotation_selected$food_variable == paste0(food, "_group_range") &
          rotation_selected$model == "model2",
        "lowest_risk_reference_category"
      ]
      length(z) == 1L &&
        identical(as.character(z), as.character(
          rcs_diag$reference_category_label[match(food, rcs_diag$food)]
        ))
    }, logical(1))),
    all(derived_six$derived_direction %in% c("middle", "high")) &&
      all(is.finite(derived_six$derived_value1)) &&
      all(derived_six$derived_direction != "middle" |
            (is.finite(derived_six$derived_value2) &
               derived_six$derived_value2 >= derived_six$derived_value1)),
    all(derived_six$derived_value1 >= derived_six$prediction_support_low - 1e-8) &&
      all(derived_six$derived_direction != "middle" |
            derived_six$derived_value2 <= derived_six$prediction_support_high + 1e-8),
    with(derived_six[derived_six$food == "processed_meats", ],
         derived_value1 >= plateau_05_low - 1e-8 &
           derived_value1 <= plateau_05_high + 1e-8),
    with(derived_six[derived_six$food == "sweetened_bev", ],
         abs(derived_value1) < 1e-8 && abs(derived_value2 - 1.0) < 1e-6),
    all(vapply(score_cols, function(sc) {
      z <- score_data[[sc]]
      all(is.na(z) | z %in% c(0, 1))
    }, logical(1))),
    all(binary_distribution$n_score_0 > 0L &
          binary_distribution$n_score_1 > 0L),
    nrow(score_data) == EXPECTED_TRAIN_N,
    identical(as.character(new_ranges$food), as.character(new_rules$food)) &&
      all(new_ranges$range_direction == new_rules$range_direction) &&
      all(abs(new_ranges$value1 - new_rules$value1) < 1e-10) &&
      all((is.na(new_ranges$value2) & is.na(new_rules$value2)) |
            (!is.na(new_ranges$value2) & !is.na(new_rules$value2) &
               abs(new_ranges$value2 - new_rules$value2) < 1e-10))
  )
)

# -----------------------------------------------------------------------------
# 11. Write all V5.2 outputs
# -----------------------------------------------------------------------------

write.csv(
  new_ranges, file.path(OUT_DIR, "food_rcs_range_all.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  new_rules, file.path(OUT_DIR, "RWDID_binary_scoring_rules_15_components.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  score_data, file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  derived_six,
  file.path(AUDIT_DIR, "Table_Sx_six_component_data_derived_intervals.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  target_comparison,
  file.path(AUDIT_DIR, "Table_Sx_target_vs_data_derived_ranges.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  pairwise_evidence,
  file.path(AUDIT_DIR, "Table_Sx_Model2_lowest_risk_pairwise_evidence.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  pointwise_evidence,
  file.path(AUDIT_DIR, "Table_Sx_RCS_pointwise_CI_regions.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  binary_distribution,
  file.path(AUDIT_DIR, "RWDID_binary_score_distribution_training.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  total_distribution,
  file.path(AUDIT_DIR, "RWDID_total_score_distribution_training.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  qc, file.path(AUDIT_DIR, "V5_2_final_QC.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# Preserve frozen source diagnostics in the new output folder for traceability.
file.copy(V51_MODEL_FILE, file.path(AUDIT_DIR, "frozen_RCS_model_summary.csv"),
          overwrite = TRUE)
file.copy(V51_REFERENCE_AUDIT_FILE,
          file.path(AUDIT_DIR, "frozen_lowest_risk_reference_audit.csv"),
          overwrite = TRUE)
file.copy(V51_QC_FILE, file.path(AUDIT_DIR, "source_RCS_final_QC.csv"),
          overwrite = TRUE)

if (any(!qc$status)) {
  print(qc)
  stop("STOP: V5.2 final QC failed. Outputs were retained for inspection.")
}

cat("\n============================================================\n")
cat("STEP_04_1_1 DATA-DERIVED INTERVAL RECONCILIATION V5.2 COMPLETED\n")
cat("Frozen training N: ", nrow(score_data), "\n", sep = "")
cat("Training events: ", EXPECTED_TRAIN_EVENTS, "\n", sep = "")
cat("RCS models were NOT refitted; frozen prediction grids were used.\n")
cat("Six target component intervals were derived from corrected-reference RCS\n")
cat("5% minimum-risk plateaus with Model 2 lowest-risk Cox cross-checks.\n")
cat("Other nine component intervals remain unchanged from frozen.\n")
cat("Score missingness remains NA; no imputation was performed.\n")
cat("QC: ", sum(qc$status), "/", nrow(qc), " PASS\n", sep = "")
cat("New interval CSV: ", file.path(OUT_DIR, "food_rcs_range_all.csv"), "\n", sep = "")
cat("New training scores: ",
    file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv"), "\n", sep = "")
cat("Target comparison: ",
    file.path(AUDIT_DIR, "Table_Sx_target_vs_data_derived_ranges.csv"), "\n", sep = "")
cat("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 06/11: Step_04_1_1_2_INTERVAL_RECONCILIATION
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 07/11: Step_04_1_1_3_RCS20_OVERLAP
# Source: Step_04_1_1_3__rcs_plateau_overlap.R
# Provenance: project source archive data_analysis_v2.zip
# SHA256 of included block: 4550ac4991ea72ff67a85a796bf58fb8c2f7241e2469dd2013843130f7efd968
########################################################################################

# =============================================================================
# STEP_04_1_1 RWDID V5.3 - MINIMAL-CHANGE RCS PLATEAU OVERLAP SENSITIVITY
#
# Primary analysis:
#   Frozen is preserved without modification.
#
# Sensitivity analyses:
#   RCS05 / RCS10 / RCS20:
#     For six reviewer-critical components, intersect the frozen scoring
#     interval with the corresponding corrected-reference RCS 5%, 10%, or 20%
#     minimum-risk plateau.
#   The other nine components retain their exact scoring rules.
#
# This is NOT an HRS-outcome-driven optimization. HRS results are not read by
# this script and cannot influence the interval definitions.
#
# Frozen inputs:
#   ./results/Step_04/Step_04_1_1_1/
#   ./results/Step_03/ukb_70_30_split/main/ukb_depression_training_cohort.rds
#
# Output:
#   ./results/Step_04/Step_04_1_1_3/Step_04_1_1_3_RCS05_OVERLAP/
#   ./results/Step_04/Step_04_1_1_3/Step_04_1_1_3_RCS10_OVERLAP/
#   ./results/Step_04/Step_04_1_1_3/Step_04_1_1_3_RCS20_OVERLAP/
#   ./results/Step_04/Step_04_1_1_3/audit/
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

options(stringsAsFactors = FALSE, scipen = 999)

# -----------------------------------------------------------------------------
# 0. Paths and frozen design
# -----------------------------------------------------------------------------

V51_DIR <- "./results/Step_04/Step_04_1_1_1"
TRAIN_RDS <- "./results/Step_03/ukb_70_30_split/main/ukb_depression_training_cohort.rds"

MASTER_OUT <- "./results/Step_04/Step_04_1_1_3"
MASTER_AUDIT <- file.path(MASTER_OUT, "audit")
dir.create(MASTER_AUDIT, recursive = TRUE, showWarnings = FALSE)

EXPECTED_TRAIN_N <- 76578L
EXPECTED_TRAIN_EVENTS <- 2122L
ID_COL <- "eid"

# Six components specifically raised during the reviewer-response discussion.
TARGET_FOODS <- c(
  "citrus_fruits",
  "fruit_juices",
  "eggs",
  "low_fat_dairy",
  "processed_meats",
  "sweetened_bev"
)

VARIANTS <- c("RCS05", "RCS10", "RCS20")
PLATEAU_SUFFIX <- c(
  RCS05 = "05",
  RCS10 = "10",
  RCS20 = "20"
)

V51_RANGE_FILE <- file.path(V51_DIR, "food_rcs_range_all.csv")
V51_RULE_FILE <- file.path(
  V51_DIR, "RWDID_binary_scoring_rules_15_components.csv"
)
V51_SCORE_FILE <- file.path(
  V51_DIR, "RWDID_binary_component_scores_training.csv"
)
V51_QC_FILE <- file.path(V51_DIR, "audit", "RCS_final_QC.csv")
V51_DIAG_FILE <- file.path(
  V51_DIR, "audit", "RCS_model_summary_and_range_diagnostics.csv"
)
ROTATION_DIR <- "./results/Step_03/reference_rotation_lowest_risk"
ROTATION_SELECTED_FILE <- file.path(
  ROTATION_DIR, "food_group_reference_rotation_selected_training.csv"
)
ROTATION_PAIRWISE_FILE <- file.path(
  ROTATION_DIR, "food_group_reference_rotation_all_comparisons_training.csv"
)

required_files <- c(
  V51_RANGE_FILE, V51_RULE_FILE, V51_SCORE_FILE,
  V51_QC_FILE, V51_DIAG_FILE, TRAIN_RDS,
  ROTATION_SELECTED_FILE, ROTATION_PAIRWISE_FILE
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "STOP: Missing frozen inputs:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1. Read and validate the frozen outputs
# -----------------------------------------------------------------------------

v51_ranges <- read.csv(
  V51_RANGE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
v51_rules <- read.csv(
  V51_RULE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
v51_scores <- read.csv(
  V51_SCORE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
v51_qc <- read.csv(
  V51_QC_FILE, stringsAsFactors = FALSE, check.names = FALSE
)
rcs_diag <- read.csv(
  V51_DIAG_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
rotation_selected <- read.csv(
  ROTATION_SELECTED_FILE, stringsAsFactors = FALSE, check.names = FALSE
)
rotation_pairwise <- read.csv(
  ROTATION_PAIRWISE_FILE, stringsAsFactors = FALSE, check.names = FALSE
)

if (!all(c(
  "model", "food_variable", "lowest_risk_reference_category",
  "model_n", "model_events"
) %in% names(rotation_selected)) ||
    !all(c(
      "model", "food_variable", "reference_category", "comparison_category",
      "is_selected_lowest_risk_reference", "is_reference_row",
      "hazard_ratio", "ci_low", "ci_high", "p_value_raw"
    ) %in% names(rotation_pairwise))) {
  stop("STOP: Reference-rotation files have an unexpected schema.")
}

if (nrow(v51_ranges) != 15L || nrow(v51_rules) != 15L ||
    nrow(rcs_diag) != 15L || nrow(v51_scores) != EXPECTED_TRAIN_N) {
  stop("STOP: Frozen dimensions are not 15 components / 76,578 training rows.")
}
if (!all(as.logical(v51_qc$status))) {
  stop("STOP: Frozen final QC is not entirely PASS.")
}
if (!all(TARGET_FOODS %in% v51_ranges$food) ||
    !all(TARGET_FOODS %in% rcs_diag$food)) {
  stop("STOP: One or more target foods are missing from frozen files.")
}
if (!all(c(
  "plateau_05_low", "plateau_05_high",
  "plateau_10_low", "plateau_10_high",
  "plateau_20_low", "plateau_20_high",
  "prediction_support_low", "prediction_support_high",
  "curve_shape", "p_overall", "p_non_linear"
) %in% names(rcs_diag))) {
  stop("STOP: RCS diagnostic file lacks the required plateau/evidence fields.")
}

# -----------------------------------------------------------------------------
# 2. Validate exact training cohort and participant IDs
# -----------------------------------------------------------------------------

training_data <- readRDS(TRAIN_RDS)
if (!is.data.frame(training_data)) {
  stop("STOP: Frozen training RDS is not a data.frame.")
}
required_training_cols <- c(ID_COL, "depression_status", TARGET_FOODS)
if (!all(required_training_cols %in% names(training_data))) {
  stop(
    "STOP: Frozen training RDS lacks required columns: ",
    paste(setdiff(required_training_cols, names(training_data)), collapse = ", ")
  )
}
if (nrow(training_data) != EXPECTED_TRAIN_N ||
    anyDuplicated(as.character(training_data[[ID_COL]])) > 0L) {
  stop("STOP: Frozen training RDS N or participant IDs are invalid.")
}
if (anyNA(training_data$depression_status) ||
    !all(training_data$depression_status %in% c(0, 1)) ||
    sum(training_data$depression_status == 1L) != EXPECTED_TRAIN_EVENTS) {
  stop("STOP: Frozen training event coding/count does not match 2,122 events.")
}
if (!all(c(ID_COL, "analysis_split", "n_components_scored",
           "RWDID_total_score_complete15") %in% names(v51_scores))) {
  stop("STOP: participant-level score file has an unexpected schema.")
}

normalize_id <- function(x) {
  z <- trimws(as.character(x))
  z[z == ""] <- NA_character_
  sub("\\.0+$", "", z)
}

train_ids <- normalize_id(training_data[[ID_COL]])
v51_ids <- normalize_id(v51_scores[[ID_COL]])

if (anyNA(train_ids) || anyNA(v51_ids) ||
    anyDuplicated(train_ids) > 0L || anyDuplicated(v51_ids) > 0L ||
    !setequal(train_ids, v51_ids)) {
  stop("STOP: Frozen training RDS IDs do not exactly match score IDs.")
}
if (any(as.character(v51_scores$analysis_split) != "training")) {
  stop("STOP: score file contains non-training rows.")
}

# Align raw intake data to the participant order in the frozen score file.
train_match <- match(v51_ids, train_ids)
if (anyNA(train_match)) {
  stop("STOP: Could not align raw training intake values to frozen scores.")
}
training_aligned <- training_data[train_match, , drop = FALSE]

# -----------------------------------------------------------------------------
# 3. Scoring helpers: exactly the direction-aware numeric logic
# -----------------------------------------------------------------------------

score_one_component <- function(x, lower, upper, food_name) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  out <- rep(NA_integer_, length(x))
  observed <- !is.na(x)

  if (any(!is.finite(x[observed]))) {
    stop("STOP: Non-finite non-missing intake values for ", food_name)
  }
  if (any(x[observed] < -1e-8)) {
    stop("STOP: Negative intake values for ", food_name)
  }
  if (!any(observed)) return(out)

  scale_value <- max(1, abs(lower), ifelse(is.finite(upper), abs(upper), 0))
  tol <- max(1e-8, scale_value * 1e-8)

  if (is.infinite(upper) && upper > 0) {
    out[observed] <- as.integer(x[observed] >= lower - tol)
  } else {
    out[observed] <- as.integer(
      x[observed] >= lower - tol & x[observed] <= upper + tol
    )
  }
  out
}

format_num <- function(x, digits = 6L) {
  if (length(x) != 1L || !is.finite(x)) return(NA_character_)
  if (abs(x) < 1e-12) return("0")
  formatC(x, format = "f", digits = digits)
}

make_rule <- function(lower, upper) {
  paste0(format_num(lower), " <= x <= ", format_num(upper))
}

# Convert the frozen rule to numeric closed bounds for intersection.
get_v51_bounds <- function(row) {
  direction <- as.character(row$range_direction[1])
  value1 <- as.numeric(row$value1[1])
  value2 <- as.numeric(row$value2[1])
  score_lower <- as.numeric(row$score_lower[1])
  score_upper <- as.numeric(row$score_upper[1])

  if (direction == "middle") {
    c(lower = value1, upper = value2)
  } else if (direction == "low") {
    c(lower = score_lower, upper = value1)
  } else if (direction == "high") {
    c(lower = value1, upper = score_upper)
  } else if (direction == "point") {
    c(lower = value1, upper = value1)
  } else {
    stop("STOP: Unsupported range direction for ", row$food[1])
  }
}

# -----------------------------------------------------------------------------
# 4. Build overlap audit and three sensitivity range sets
# -----------------------------------------------------------------------------

all_overlap_rows <- list()
variant_results <- list()

for (variant in VARIANTS) {
  suffix <- PLATEAU_SUFFIX[[variant]]
  low_col <- paste0("plateau_", suffix, "_low")
  high_col <- paste0("plateau_", suffix, "_high")

  v_ranges <- v51_ranges
  v_rules <- v51_rules
  derived_rows <- list()

  for (food in TARGET_FOODS) {
    idx_range <- match(food, v_ranges$food)
    idx_rule <- match(food, v_rules$food)
    idx_diag <- match(food, rcs_diag$food)

    if (anyNA(c(idx_range, idx_rule, idx_diag))) {
      stop("STOP: Cannot match all frozen records for ", food)
    }

    range_row <- v_ranges[idx_range, , drop = FALSE]
    diag_row <- rcs_diag[idx_diag, , drop = FALSE]

    selected_ref <- rotation_selected[
      rotation_selected$model == "model2" &
        rotation_selected$food_variable == paste0(food, "_group_range"),
      , drop = FALSE
    ]
    if (nrow(selected_ref) != 1L ||
        as.character(selected_ref$lowest_risk_reference_category[1]) !=
          as.character(diag_row$reference_category_label[1]) ||
        as.integer(selected_ref$model_n[1]) != EXPECTED_TRAIN_N ||
        as.integer(selected_ref$model_events[1]) != EXPECTED_TRAIN_EVENTS) {
      stop(
        "STOP: Model 2 lowest-risk reference rotation does not match frozen ",
        "RCS reference/sample for ", food
      )
    }

    cox_pair <- rotation_pairwise[
      rotation_pairwise$model == "model2" &
        rotation_pairwise$food_variable == paste0(food, "_group_range") &
        as.logical(rotation_pairwise$is_selected_lowest_risk_reference) &
        !as.logical(rotation_pairwise$is_reference_row),
      , drop = FALSE
    ]
    if (nrow(cox_pair) < 1L) {
      stop("STOP: No Model 2 pairwise reference comparisons for ", food)
    }
    higher_risk_categories <- as.character(
      cox_pair$comparison_category[cox_pair$ci_low > 1]
    )
    lower_risk_categories <- as.character(
      cox_pair$comparison_category[cox_pair$ci_high < 1]
    )
    cox_ref_label <- as.character(selected_ref$lowest_risk_reference_category[1])

    anchor <- get_v51_bounds(range_row)
    plateau_low <- as.numeric(diag_row[[low_col]][1])
    plateau_high <- as.numeric(diag_row[[high_col]][1])
    support_low <- as.numeric(diag_row$prediction_support_low[1])
    support_high <- as.numeric(diag_row$prediction_support_high[1])

    if (!is.finite(anchor["lower"]) || is.na(anchor["upper"]) ||
        is.nan(anchor["upper"]) ||
        !all(is.finite(c(plateau_low, plateau_high,
                         support_low, support_high))) ||
        anchor["upper"] < anchor["lower"] ||
        plateau_high < plateau_low ||
        plateau_low < support_low - 1e-8 ||
        plateau_high > support_high + 1e-8) {
      stop("STOP: Invalid frozen anchor or RCS plateau for ", food, " / ", variant)
    }

    overlap_low <- max(anchor["lower"], plateau_low)
    overlap_high <- min(anchor["upper"], plateau_high)

    if (overlap_high < overlap_low - 1e-8) {
      stop(
        "STOP: interval and ", variant,
        " plateau do not overlap for ", food,
        ". No nearest-boundary substitution will be made."
      )
    }

    overlap_low <- max(0, overlap_low)
    overlap_high <- max(overlap_low, overlap_high)

    # Default scoring rule is the finite intersection.
    scoring_direction <- "middle"
    scoring_value1 <- overlap_low
    scoring_value2 <- overlap_high
    scoring_lower <- overlap_low
    scoring_upper <- overlap_high
    new_rule <- make_rule(overlap_low, overlap_high)
    tail_note <- "No open-tail extension."

    # Processed meats is a pre-specified decreasing-risk / high-intake pattern.
    # Its corrected-reference RCS plateau reaches the upper prediction-support
    # boundary (3.0), while the frozen Model 2 lowest-risk reference is the
    # highest intake category (2-13). Retain the open upper scoring tail
    # from the RCS-supported threshold, but explicitly state that continuous
    # RCS evidence does not extend beyond 3.0; the tail continuation is
    # categorical Cox support only.
    if (food == "processed_meats") {
      ref_label <- as.character(diag_row$reference_category_label[1])
      if (!grepl("2[[:space:]]*-[[:space:]]*13", ref_label) ||
          !identical(ref_label, cox_ref_label) ||
          abs(plateau_high - support_high) > 1e-8 ||
          as.character(diag_row$curve_shape[1]) !=
            "reverse-J-shaped_decreasing-risk") {
        stop(
          "STOP: Processed-meat open-tail retention conditions are not met. ",
          "Reference=", ref_label,
          "; plateau upper=", plateau_high,
          "; support upper=", support_high
        )
      }
      scoring_direction <- "high"
      scoring_value1 <- overlap_low
      scoring_value2 <- NA_real_
      scoring_lower <- overlap_low
      scoring_upper <- Inf
      new_rule <- paste0("x >= ", format_num(overlap_low))
      tail_note <- paste0(
        "Open upper tail retained because Model 2 lowest-risk category is ",
        ref_label,
        "; continuous RCS plateau/support ends at ",
        format_num(plateau_high),
        ", so values above that point are supported categorically, not by ",
        "continuous RCS."
      )
    }

    derived_rows[[food]] <- tibble(
      variant = variant,
      food = food,
      display_name = as.character(diag_row$display_name[1]),
      curve_shape = as.character(diag_row$curve_shape[1]),
      p_overall = as.numeric(diag_row$p_overall[1]),
      p_non_linear = as.numeric(diag_row$p_non_linear[1]),
      Model2_lowest_risk_reference_category = cox_ref_label,
      Model2_significantly_higher_risk_categories = ifelse(
        length(higher_risk_categories) > 0L,
        paste(higher_risk_categories, collapse = "; "),
        "None"
      ),
      Model2_significantly_lower_risk_categories = ifelse(
        length(lower_risk_categories) > 0L,
        paste(lower_risk_categories, collapse = "; "),
        "None"
      ),
      Model2_pairwise_discordance_flag = length(higher_risk_categories) > 0L,
      direction = as.character(range_row$range_direction[1]),
      value1 = as.numeric(range_row$value1[1]),
      value2 = as.numeric(range_row$value2[1]),
      anchor_lower = unname(anchor["lower"]),
      anchor_upper = unname(anchor["upper"]),
      RCS_plateau_low = plateau_low,
      RCS_plateau_high = plateau_high,
      RCS_prediction_support_low = support_low,
      RCS_prediction_support_high = support_high,
      overlap_low = overlap_low,
      overlap_high = overlap_high,
      overlap_range = paste0(format_num(overlap_low), "–", format_num(overlap_high)),
      scoring_direction = scoring_direction,
      scoring_value1 = scoring_value1,
      scoring_value2 = scoring_value2,
      scoring_rule = new_rule,
      open_upper_tail_categorical_only = food == "processed_meats",
      tail_evidence_note = tail_note,
      interval_width = unname(anchor["upper"] - anchor["lower"]),
      overlap_width = overlap_high - overlap_low,
      percent_of_width_retained = ifelse(
        is.finite(anchor["upper"] - anchor["lower"]) &&
          (anchor["upper"] - anchor["lower"]) > 0,
        100 * (overlap_high - overlap_low) /
          (anchor["upper"] - anchor["lower"]),
        NA_real_
      ),
      interval_source = paste0(
        "V5.3_intersection_of_frozen_anchor_and_RCS_",
        suffix, "pct_minimum_risk_plateau"
      )
    )

    # Update only the six target components in this sensitivity copy.
    v_ranges$value1[idx_range] <- scoring_value1
    v_ranges$value2[idx_range] <- scoring_value2
    v_ranges$range_direction[idx_range] <- scoring_direction
    v_ranges$score_lower[idx_range] <- scoring_lower
    v_ranges$score_upper[idx_range] <- scoring_upper
    v_ranges$score_rule[idx_range] <- new_rule
    v_ranges$range_source[idx_range] <- paste0(
      "V5.3_RCS_", suffix, "pct_plateau_intersected_with_anchor"
    )
    if ("interval_source" %in% names(v_ranges)) {
      v_ranges$interval_source[idx_range] <- paste0(
        "V5.3 RCS ", suffix, "% plateau intersected with frozen interval"
      )
    }
    if (!"original_value1" %in% names(v_ranges)) {
      v_ranges$original_value1 <- NA_real_
      v_ranges$original_value2 <- NA_real_
      v_ranges$original_direction <- NA_character_
    }
    v_ranges$original_value1[idx_range] <- range_row$value1
    v_ranges$original_value2[idx_range] <- range_row$value2
    v_ranges$original_direction[idx_range] <- range_row$range_direction
    v_ranges$RCS_plateau_low[idx_range] <- plateau_low
    v_ranges$RCS_plateau_high[idx_range] <- plateau_high

    # Update matching scoring-rule record.
    v_rules$value1[idx_rule] <- scoring_value1
    v_rules$value2[idx_rule] <- scoring_value2
    v_rules$range_direction[idx_rule] <- scoring_direction
    v_rules$score_lower[idx_rule] <- scoring_lower
    v_rules$score_upper[idx_rule] <- scoring_upper
    v_rules$score_rule[idx_rule] <- new_rule
    v_rules$score_1[idx_rule] <- paste0("1 if ", new_rule)
    v_rules$score_0[idx_rule] <- "0 otherwise"
    v_rules$missing_rule[idx_rule] <- "NA; no imputation"
    v_rules$range_source[idx_rule] <- paste0(
      "V5.3_RCS_", suffix, "pct_plateau_intersected_with_anchor"
    )
    v_rules$interval_source[idx_rule] <- paste0(
      "V5.3 RCS ", suffix, "% plateau intersected with frozen interval"
    )
    if (!"original_value1" %in% names(v_rules)) {
      v_rules$original_value1 <- NA_real_
      v_rules$original_value2 <- NA_real_
      v_rules$original_direction <- NA_character_
    }
    v_rules$original_value1[idx_rule] <- range_row$value1
    v_rules$original_value2[idx_rule] <- range_row$value2
    v_rules$original_direction[idx_rule] <- range_row$range_direction
    v_rules$RCS_plateau_low[idx_rule] <- plateau_low
    v_rules$RCS_plateau_high[idx_rule] <- plateau_high
  }

  overlap_table <- bind_rows(derived_rows)
  all_overlap_rows[[variant]] <- overlap_table

  # Recalculate the six revised binary components from the exact frozen
  # training intake values. The other nine scores are copied unchanged from
  # the frozen participant-level score file.
  score_data <- v51_scores
  for (food in TARGET_FOODS) {
    z <- overlap_table[overlap_table$food == food, , drop = FALSE]
    score_col <- paste0("score_", food)
    score_data[[score_col]] <- score_one_component(
      x = training_aligned[[food]],
      lower = as.numeric(z$scoring_value1[1]),
      upper = if (is.na(z$scoring_value2[1])) Inf else as.numeric(z$scoring_value2[1]),
      food_name = food
    )
  }

  score_cols <- paste0(
    "score_",
    c(
      "citrus_fruits", "red_orange_veg", "other_veg", "potatoes",
      "refined_grains", "eggs", "low_fat_dairy", "red_meats",
      "processed_meats", "fish_seafood", "coffee", "fruit_juices",
      "sweetened_bev", "water", "wine"
    )
  )
  if (!all(score_cols %in% names(score_data))) {
    stop("STOP: One or more of the 15 score columns are missing.")
  }

  score_matrix <- as.matrix(score_data[, score_cols, drop = FALSE])
  storage.mode(score_matrix) <- "numeric"
  score_data$n_components_scored <- as.integer(rowSums(!is.na(score_matrix)))
  score_data$RWDID_total_score_complete15 <- ifelse(
    score_data$n_components_scored == 15L,
    as.integer(rowSums(score_matrix, na.rm = TRUE)),
    NA_integer_
  )
  score_data$analysis_split <- "training"

  # Component-level distribution and changed-score audit.
  component_distribution <- bind_rows(lapply(score_cols, function(sc) {
    z <- score_data[[sc]]
    z0 <- v51_scores[[sc]]
    tibble(
      variant = variant,
      score_column = sc,
      n_training = length(z),
      n_score_0 = sum(z == 0, na.rm = TRUE),
      n_score_1 = sum(z == 1, na.rm = TRUE),
      n_missing = sum(is.na(z)),
      pct_score_1 = 100 * mean(z == 1, na.rm = TRUE),
      n_changed_vs_V5_1 = sum(
        (is.na(z) != is.na(z0)) |
          (!is.na(z) & !is.na(z0) & z != z0),
        na.rm = TRUE
      )
    )
  }))

  complete <- score_data$n_components_scored == 15L &
    !is.na(score_data$RWDID_total_score_complete15)
  primary_complete <- v51_scores$n_components_scored == 15L &
    !is.na(v51_scores$RWDID_total_score_complete15)

  common <- complete & primary_complete
  score_cor <- if (sum(common) >= 3L) {
    suppressWarnings(cor(
      score_data$RWDID_total_score_complete15[common],
      v51_scores$RWDID_total_score_complete15[common],
      method = "pearson"
    ))
  } else NA_real_
  score_spearman <- if (sum(common) >= 3L) {
    suppressWarnings(cor(
      score_data$RWDID_total_score_complete15[common],
      v51_scores$RWDID_total_score_complete15[common],
      method = "spearman"
    ))
  } else NA_real_

  total_score_comparison <- tibble(
    variant = variant,
    n_common_complete15 = sum(common),
    pearson_r_vs_V5_1 = score_cor,
    spearman_rho_vs_V5_1 = score_spearman,
    mean_total_score = mean(
      v51_scores$RWDID_total_score_complete15[common], na.rm = TRUE
    ),
    sensitivity_mean_total_score = mean(
      score_data$RWDID_total_score_complete15[common], na.rm = TRUE
    ),
    sd_total_score = sd(
      v51_scores$RWDID_total_score_complete15[common], na.rm = TRUE
    ),
    sensitivity_sd_total_score = sd(
      score_data$RWDID_total_score_complete15[common], na.rm = TRUE
    )
  )

  # Variant-specific output directory follows the exact file contract expected
  # by the Beta training script.
  variant_dir <- file.path(
    "./results/Step_04/Step_04_1_1_3",
    paste0("Step_04_1_1_3_", variant, "_OVERLAP")
  )
  variant_audit <- file.path(variant_dir, "audit")
  dir.create(variant_audit, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    v_ranges, file.path(variant_dir, "food_rcs_range_all.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
  write.csv(
    v_rules,
    file.path(variant_dir, "RWDID_binary_scoring_rules_15_components.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
  write.csv(
    score_data,
    file.path(variant_dir, "RWDID_binary_component_scores_training.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
  write.csv(
    overlap_table,
    file.path(variant_audit, "component_interval_overlap.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
  cox_pairwise_export <- rotation_pairwise %>%
    filter(
      model == "model2",
      food_variable %in% paste0(TARGET_FOODS, "_group_range"),
      as.logical(is_selected_lowest_risk_reference)
    ) %>%
    mutate(
      CI_includes_1 = ci_low <= 1 & ci_high >= 1,
      significantly_higher_risk_than_reference = ci_low > 1,
      significantly_lower_risk_than_reference = ci_high < 1
    )
  write.csv(
    cox_pairwise_export,
    file.path(variant_audit, "Model2_lowest_risk_pairwise_crosscheck.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
  write.csv(
    component_distribution,
    file.path(variant_audit, "binary_component_score_distribution.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
  write.csv(
    total_score_comparison,
    file.path(variant_audit, "total_score_comparison_vs_V5_1.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )

  # Preserve the exact frozen diagnostics used to define this sensitivity.
  file.copy(V51_DIAG_FILE,
            file.path(variant_audit, "frozen_RCS_model_summary_and_range_diagnostics.csv"),
            overwrite = TRUE)
  file.copy(V51_QC_FILE,
            file.path(variant_audit, "source_RCS_final_QC.csv"),
            overwrite = TRUE)

  # Variant QC.
  unchanged_foods <- setdiff(v51_ranges$food, TARGET_FOODS)
  unchanged_score_cols <- paste0("score_", unchanged_foods)
  qc <- tibble(
    check = c(
      "source_QC_all_pass",
      "Model2_reference_rotation_matches_RCS_for_all_six",
      "exact_training_ID_match",
      "expected_training_N",
      "exactly_15_components",
      "six_target_foods_reconciled",
      "all_target_intersections_nonempty",
      "all_target_intersections_inside_RCS_plateau",
      "all_target_intersections_inside_anchor",
      "non_target_nine_ranges_unchanged",
      "non_target_nine_scores_unchanged",
      "all_scores_are_0_1_or_NA",
      "all_15_components_have_both_0_and_1",
      "participant_score_rows_preserved",
      "total_score_only_for_complete_15"
    ),
    status = c(
      all(as.logical(v51_qc$status)),
      all(vapply(TARGET_FOODS, function(food) {
        z <- rotation_selected[
          rotation_selected$model == "model2" &
            rotation_selected$food_variable == paste0(food, "_group_range"),
          , drop = FALSE
        ]
        nrow(z) == 1L &&
          as.character(z$lowest_risk_reference_category[1]) ==
            as.character(rcs_diag$reference_category_label[match(food, rcs_diag$food)])
      }, logical(1))),
      setequal(train_ids, v51_ids),
      nrow(score_data) == EXPECTED_TRAIN_N,
      length(score_cols) == 15L,
      setequal(overlap_table$food, TARGET_FOODS) && nrow(overlap_table) == 6L,
      all(overlap_table$overlap_high >= overlap_table$overlap_low),
      all(overlap_table$overlap_low >= overlap_table$RCS_plateau_low - 1e-8) &&
        all(overlap_table$overlap_high <= overlap_table$RCS_plateau_high + 1e-8),
      all(overlap_table$overlap_low >= overlap_table$anchor_lower - 1e-8) &&
        all(overlap_table$overlap_high <= overlap_table$anchor_upper + 1e-8),
      all(v_ranges$value1[match(unchanged_foods, v_ranges$food)] ==
            v51_ranges$value1[match(unchanged_foods, v51_ranges$food)]) &&
        all(as.character(v_ranges$range_direction[match(unchanged_foods, v_ranges$food)]) ==
              as.character(v51_ranges$range_direction[match(unchanged_foods, v51_ranges$food)])) &&
        all(as.character(v_ranges$score_rule[match(unchanged_foods, v_ranges$food)]) ==
              as.character(v51_ranges$score_rule[match(unchanged_foods, v51_ranges$food)])) &&
        all((is.na(v_ranges$value2[match(unchanged_foods, v_ranges$food)]) &
               is.na(v51_ranges$value2[match(unchanged_foods, v51_ranges$food)])) |
              (!is.na(v_ranges$value2[match(unchanged_foods, v_ranges$food)]) &
                 !is.na(v51_ranges$value2[match(unchanged_foods, v51_ranges$food)]) &
                 v_ranges$value2[match(unchanged_foods, v_ranges$food)] ==
                   v51_ranges$value2[match(unchanged_foods, v51_ranges$food)])),
      all(vapply(unchanged_score_cols, function(sc) {
        identical(score_data[[sc]], v51_scores[[sc]])
      }, logical(1))),
      all(vapply(score_cols, function(sc) {
        z <- score_data[[sc]]
        all(is.na(z) | z %in% c(0, 1))
      }, logical(1))),
      all(component_distribution$n_score_0 > 0L &
            component_distribution$n_score_1 > 0L),
      nrow(score_data) == nrow(v51_scores),
      all(is.na(score_data$RWDID_total_score_complete15) |
            score_data$n_components_scored == 15L)
    )
  )

  write.csv(
    qc, file.path(variant_audit, "final_QC.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )

  if (any(!qc$status)) {
    print(qc)
    stop("STOP: V5.3 QC failed for ", variant, ". Outputs retained for review.")
  }

  variant_results[[variant]] <- list(
    overlap = overlap_table,
    distribution = component_distribution,
    total = total_score_comparison,
    qc = qc
  )
}

# -----------------------------------------------------------------------------
# 5. Master evidence tables
# -----------------------------------------------------------------------------

all_overlap <- bind_rows(all_overlap_rows)
all_distribution <- bind_rows(lapply(variant_results, `[[`, "distribution"))
all_total_comparison <- bind_rows(lapply(variant_results, `[[`, "total"))
all_qc <- bind_rows(lapply(names(variant_results), function(v) {
  x <- variant_results[[v]]$qc
  x$variant <- v
  x
}))

write.csv(
  all_overlap,
  file.path(MASTER_AUDIT, "RCS05_RCS10_RCS20_overlap_comparison.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  all_distribution,
  file.path(MASTER_AUDIT, "all_variant_component_score_distributions.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  all_total_comparison,
  file.path(MASTER_AUDIT, "total_score_comparison_vs_V5_1.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  all_qc,
  file.path(MASTER_AUDIT, "all_variant_QC.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
file.copy(
  ROTATION_SELECTED_FILE,
  file.path(MASTER_AUDIT, "source_Model2_lowest_risk_reference_rotation.csv"),
  overwrite = TRUE
)
file.copy(
  ROTATION_PAIRWISE_FILE,
  file.path(MASTER_AUDIT, "source_Model2_pairwise_comparisons.csv"),
  overwrite = TRUE
)

if (any(!all_qc$status)) {
  stop("STOP: At least one V5.3 variant failed QC.")
}

cat("\n============================================================\n")
cat("RWDID V5.3 MINIMAL-CHANGE RCS-PLATEAU SENSITIVITY COMPLETED\n")
cat("Frozen primary files were not modified.\n")
cat("Training N: ", EXPECTED_TRAIN_N, "\n", sep = "")
cat("Training events: ", EXPECTED_TRAIN_EVENTS, "\n", sep = "")
cat("Variants generated: RCS05, RCS10, RCS20\n")
cat("Only six reviewer-critical components were intersected with RCS plateaus.\n")
cat("The other nine components were retained exactly from V5.1.\n")
cat("HRS outcomes and HRS P values were not read by this script.\n")
cat("Master audit: ", normalizePath(MASTER_AUDIT, mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 07/11: Step_04_1_1_3_RCS20_OVERLAP
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 08/11: Step_04_1_1_4_FINAL_RCS_BINARY
# Source: Step_04_1_1_4_RCS_analysis_between_food_and_depression_FINAL_B500_PARALLEL_BOOTSTRAP_5_10_20_FONT_REVISION.R
# Provenance: Library source dated 2026-11-23
# SHA256 of included block: 5bdcedc3279558f02ebfdcb6814bf75db35c2f10643ccc701f2533d73517e54d
########################################################################################

# =============================================================================
# Step_04_1_1_RCS_analysis_between_food_and_depression_FINAL.R
#
# FINAL: RCS analysis and binary component scoring ONLY.
# Beta-weight estimation and weighted-score association remain in the separate
# Step_04_2_1_beta_weighted_diet_score_training script.
# The frozen V5.3 15-component intervals are imported; this script
# re-fits Model 3 RCS curves and exports the full reviewer-facing Table Sx.
#
# Run from the UK Biobank project root AFTER frozen Step_03_1_0 is available.
# Inputs: frozen Step_03_1_0 training RDS/manifest, Model 2 rotating-reference
# results CSV, and the existing 15-component food_rcs_range_all.csv anchor.
#
# 1. Verify the exact frozen 70/30 split and corrected 15 components.
# 2. Retain the 15 anchor intervals unchanged; fit adjusted 3-knot RCS models
#    to report P values, P10/P50/P90, curve shape, CI status and compatibility.
# 3. Run parallel 5%/10%/20% bootstrap plateau-boundary stability audits and
#    audit P90-P99 tail support.
# 4. Apply unchanged anchor intervals to training participants and generate
#    RCS, interval-comparison and binary-scoring figures.
# 5. Export a supplementary table with anchor, validation, bootstrap and tail
#    audits, followed by unified final QC. No additional report-export script
#    is required.
#
# Fixed model covariates, statistical settings and plotting aesthetics are
# retained from the validated analysis, including Fruit juices and Water,
# the 5x3 RCS P-value correction, and light PiYG score-map colors.
# Word section labels are unnumbered, tables may flow across pages, and no
# forced page breaks introduce empty half-pages. Word-table options are kept
# compatible with flextable releases that do not accept repeat_headers.
#
# RANGE_ANCHOR_CSV and FROZEN_RCS_DIAG_FILE are frozen V5.3 RCS20 inputs.
# This script derives the two documented interval intersections, exports
# the final ranges and binary component scores, and does NOT estimate betas.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(rms)
  library(cols4all)
  library(cowplot)
})

if (!requireNamespace("officer", quietly = TRUE) ||
    !requireNamespace("flextable", quietly = TRUE)) {
  stop(
    "STOP: Integrated Word-table export requires the R packages 'officer' and ",
    "'flextable'. Install them with install.packages(c('officer','flextable')) ",
    "and rerun this script."
  )
}

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 0. Frozen paths and analysis settings
# -----------------------------------------------------------------------------

SPLIT_DIR <- "./results/Step_03/ukb_70_30_split"
MAIN_SPLIT_DIR <- file.path(SPLIT_DIR, "main")
SPLIT_MANIFEST_FILE <- file.path(
  SPLIT_DIR, "ukb_70_30_split_manifest.csv"
)
TRAIN_RDS <- file.path(
  MAIN_SPLIT_DIR, "ukb_depression_training_cohort.rds"
)
REF_CSV <- paste0(
  "./results/Step_03/reference_rotation_lowest_risk/",
  "food_group_model2_results_sig.csv"
)
RANGE_ANCHOR_CSV <- paste0(
  "./results/Step_04/Step_04_1_1_3/Step_04_1_1_3_RCS20_OVERLAP/",
  "food_rcs_range_all.csv"
)
FROZEN_RCS_DIAG_FILE <- paste0(
  "./results/Step_04/Step_04_1_1_3/Step_04_1_1_3_RCS20_OVERLAP/",
  "audit/frozen_RCS_model_summary_and_range_diagnostics.csv"
)

OUT_DIR <- "./results/Step_04/Step_04_1_1_4"
PRED_DIR <- file.path(OUT_DIR, "predictions")
PLOT_DIR <- file.path(OUT_DIR, "plots")
AUDIT_DIR <- file.path(OUT_DIR, "audit")

dir.create(PRED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPECTED_SPLIT_SEED <- 20261103L
EXPECTED_TRAIN_PROP <- 0.70

PREDICTION_GRID_N <- 501L
PREDICTION_SUPPORT_PROBS <- c(0.10, 0.90)
# Add only a small right-side margin to the data-adaptive RCS x-axis.
# The displayed range is based on the actual prediction grid, knots,
# reference value, and final interval boundaries—not arbitrary fixed limits.
RCS_X_AXIS_PADDING_PROP <- 0.04
RANGE_TOLERANCES <- c(0.05, 0.10, 0.20)

# Reproducible nonparametric case bootstrap for parallel plateau-boundary
# stability. Each component is resampled/fitted B=500 times; every bootstrap
# curve is evaluated at all three tolerances (5%, 10%, 20%), so no additional
# model refits are needed for the sensitivity thresholds. The original 70/30
# split, 15 components, Cox reference, RCS knots and frozen scoring rules remain
# unchanged. The 20% tolerance is the primary near-minimum-risk threshold;
# 5% and 10% are stricter sensitivity thresholds.
BOOTSTRAP_B <- 500L
BOOTSTRAP_SEED <- 20261105L
BOOTSTRAP_MIN_SUCCESS_RATE <- 0.95
BOOTSTRAP_RANGE_TOLERANCES <- c(0.05, 0.10, 0.20)
BOOTSTRAP_BOUNDARY_CI_LEVEL <- 0.95
BOOTSTRAP_BOUNDARY_CI_METHOD <- "percentile bootstrap (2.5th and 97.5th percentiles)"

# P90+ tail sample/event support: descriptive checks, not proof of benefit.
TAIL_MIN_N_P90_P99 <- 100L
TAIL_MIN_EVENTS_P90_P99 <- 10L
TAIL_MIN_N_P99_PLUS <- 50L
TAIL_MIN_EVENTS_P99_PLUS <- 5L


# The original Model 3 adjustment set is retained without modification.
MODEL3_COVARIATES <- c(
  "age_food",
  "sex",
  "ethnicity",
  "education_attainment",
  "twnsend",
  "career",
  "smoking_status",
  "drink_status",
  "bmi_group",
  "physical_activity_proup",
  "sedentary_group",
  "household_income",
  "energy_intake",
  "central_obesity",
  "triglycerides",
  "hdl_cholesterol",
  "t2diabetes_status_final",
  "hypertension_status_final",
  "dyslipidemia_status_final",
  "cvd_status_final",
  "all_cancer_status_final"
)

# -----------------------------------------------------------------------------
# 1. Read the exact Step_03_1_0 training cohort and validate its manifest
# -----------------------------------------------------------------------------

if (!file.exists(SPLIT_MANIFEST_FILE)) {
  stop(
    "STOP: Frozen Step_03_1_0 split manifest is missing: ",
    SPLIT_MANIFEST_FILE,
    ". Run Step_03_1_0 first. This script will not create a new split."
  )
}
if (!file.exists(TRAIN_RDS)) {
  stop(
    "STOP: Step_03_1_0 training RDS is missing: ", TRAIN_RDS,
    ". Run Step_03_1_0 first. This script will not reconstruct or resample it."
  )
}
if (!file.exists(REF_CSV)) {
  stop(
    "STOP: Step_03_1_0 Model 2 significant-results CSV is missing: ",
    REF_CSV
  )
}
if (!file.exists(RANGE_ANCHOR_CSV)) {
  stop(
    "STOP: The frozen 15-component range-anchor CSV is missing: ",
    RANGE_ANCHOR_CSV,
    ". This script will not reconstruct or invent replacement cut-points."
  )
}

split_manifest <- read.csv(
  SPLIT_MANIFEST_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(
    participant_id = "character",
    split = "character",
    id_column = "character",
    split_seed = "integer",
    training_fraction = "numeric"
  )
)

required_manifest_cols <- c(
  "participant_id", "split", "id_column",
  "split_seed", "training_fraction"
)
if (!all(required_manifest_cols %in% names(split_manifest))) {
  stop("STOP: Unexpected Step_03_1_0 split-manifest structure.")
}
if (anyNA(split_manifest$participant_id) ||
    any(trimws(split_manifest$participant_id) == "") ||
    anyDuplicated(split_manifest$participant_id) > 0L) {
  stop("STOP: Invalid or duplicated participant IDs in the split manifest.")
}
if (!all(split_manifest$split %in% c("training", "test"))) {
  stop("STOP: Invalid split labels in the frozen manifest.")
}
if (length(unique(split_manifest$id_column)) != 1L ||
    anyNA(split_manifest$id_column)) {
  stop("STOP: The manifest must contain exactly one non-missing ID column.")
}
if (length(unique(split_manifest$split_seed)) != 1L ||
    anyNA(split_manifest$split_seed) ||
    unique(split_manifest$split_seed) != EXPECTED_SPLIT_SEED) {
  stop("STOP: Split seed does not match Step_03_1_0 (20261103).")
}
if (length(unique(split_manifest$training_fraction)) != 1L ||
    anyNA(split_manifest$training_fraction) ||
    abs(unique(split_manifest$training_fraction) - EXPECTED_TRAIN_PROP) > 1e-12) {
  stop("STOP: Training fraction does not match Step_03_1_0 (70%).")
}

ID_COL <- unique(split_manifest$id_column)
training_ids <- sort(
  split_manifest$participant_id[split_manifest$split == "training"],
  method = "radix"
)
test_ids <- sort(
  split_manifest$participant_id[split_manifest$split == "test"],
  method = "radix"
)

ukb_data_depression_rcs <- readRDS(TRAIN_RDS)

if (!is.data.frame(ukb_data_depression_rcs)) {
  stop("STOP: The Step_03_1_0 training RDS is not a data.frame.")
}
if (!ID_COL %in% names(ukb_data_depression_rcs)) {
  stop("STOP: ID column from manifest is absent in the training RDS: ", ID_COL)
}

training_ids_observed <- sort(
  trimws(as.character(ukb_data_depression_rcs[[ID_COL]])),
  method = "radix"
)
if (anyNA(training_ids_observed) ||
    any(training_ids_observed == "") ||
    anyDuplicated(training_ids_observed) > 0L) {
  stop("STOP: Invalid or duplicated participant IDs in the training RDS.")
}
if (!identical(training_ids_observed, training_ids)) {
  stop(
    "STOP: Training RDS IDs do not exactly match the Step_03_1_0 manifest. ",
    "No RCS models have been fitted."
  )
}
if (length(intersect(training_ids_observed, test_ids)) != 0L) {
  stop("STOP: Test participants are present in the RCS training RDS.")
}
if (nrow(ukb_data_depression_rcs) != length(training_ids)) {
  stop("STOP: Training RDS row count does not match manifest training IDs.")
}
if ("analysis_split" %in% names(ukb_data_depression_rcs) &&
    any(as.character(ukb_data_depression_rcs$analysis_split) != "training")) {
  stop("STOP: Non-training rows were found in the Step_03_1_0 training RDS.")
}
if (!all(c("depression_time", "depression_status") %in%
         names(ukb_data_depression_rcs))) {
  stop("STOP: Depression survival variables are missing from training RDS.")
}

TRAINING_QC <- tibble(
  split_source = "Step_03_1_0 frozen training RDS",
  manifest = SPLIT_MANIFEST_FILE,
  training_rds = TRAIN_RDS,
  id_column = ID_COL,
  split_seed = EXPECTED_SPLIT_SEED,
  training_fraction = EXPECTED_TRAIN_PROP,
  manifest_eligible_n = nrow(split_manifest),
  manifest_training_n = length(training_ids),
  manifest_test_n = length(test_ids),
  training_rds_n = nrow(ukb_data_depression_rcs),
  training_depression_events = sum(
    ukb_data_depression_rcs$depression_status == 1, na.rm = TRUE
  ),
  training_ids_exact_match = TRUE,
  test_ids_excluded = TRUE
)

write.csv(
  TRAINING_QC,
  file.path(AUDIT_DIR, "RCS_training_cohort_QC.csv"),
  row.names = FALSE
)

cat("\n===== FROZEN TRAINING COHORT QC =====\n")
print(TRAINING_QC)

# -----------------------------------------------------------------------------
# 2. Exact 15 components from Step_03_1_0 Model 2 significant results
# -----------------------------------------------------------------------------

FOOD_MAP <- tibble(
  group = c(
    "Citrus fruits",
    "Red/orange vegetables",
    "Other vegetables",
    "Potatoes",
    "Refined grains",
    "Eggs",
    "Low fat dairy",
    "Red meats",
    "Processed meats",
    "Fish and seafood",
    "Coffee",
    "Fruit juices",
    "Sweetened beverages",
    "Water",
    "Wine"
  ),
  food = c(
    "citrus_fruits",
    "red_orange_veg",
    "other_veg",
    "potatoes",
    "refined_grains",
    "eggs",
    "low_fat_dairy",
    "red_meats",
    "processed_meats",
    "fish_seafood",
    "coffee",
    "fruit_juices",
    "sweetened_bev",
    "water",
    "wine"
  ),
  category_var = c(
    "citrus_fruits_group_range",
    "red_orange_veg_group_range",
    "other_veg_group_range",
    "potatoes_group_range",
    "refined_grains_group_range",
    "eggs_group_range",
    "low_fat_dairy_group_range",
    "red_meats_group_range",
    "processed_meats_group_range",
    "fish_seafood_group_range",
    "coffee_group_range",
    "fruit_juices_group_range",
    "sweetened_bev_group_range",
    "water_group_range",
    "wine_group_range"
  ),
  display_name = c(
    "Citrus fruits",
    "Red/Orange vegetables",
    "Other vegetables",
    "Potatoes",
    "Refined grains",
    "Eggs",
    "Low-fat dairy",
    "Red meats",
    "Processed meats",
    "Fish and seafood",
    "Coffee",
    "Fruit juices",
    "Sweetened beverages",
    "Water",
    "Wine"
  ),
  x_break = c(1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  annotation_x = c(1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

if (nrow(FOOD_MAP) != 15L ||
    anyDuplicated(FOOD_MAP$food) > 0L ||
    anyDuplicated(FOOD_MAP$group) > 0L ||
    anyDuplicated(FOOD_MAP$category_var) > 0L) {
  stop("STOP: Internal 15-component map is invalid.")
}

if (any(c("green_leafy_veg", "high_fat_dairy") %in% FOOD_MAP$food)) {
  stop("STOP: Excluded components entered the 15-component RCS map.")
}
if (!all(c("fruit_juices", "water") %in% FOOD_MAP$food)) {
  stop("STOP: Fruit juices and water must be included in the corrected 15.")
}

# -----------------------------------------------------------------------------
# 2B. Read and freeze the previously established 15-component scoring ranges
# -----------------------------------------------------------------------------
# The V5.3 RCS20 overlap file is the frozen source of truth and carries an
# explicit range_direction column. Preserve that direction exactly:
#   middle = value1 <= intake <= value2
#   low    = intake <= value1
#   high   = intake >= value1
#   point  = intake equals value1
# Do NOT infer low versus high from a missing value2: both are one-sided rules.
# The current Cox/RCS fits below are validation diagnostics only; they do not
# replace or move the imported V5.3 RCS20 scoring intervals.
range_anchor_raw <- read.csv(
  RANGE_ANCHOR_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

required_anchor_cols <- c("food", "value1", "value2", "range_direction")
if (!all(required_anchor_cols %in% names(range_anchor_raw))) {
  stop(
    "STOP: Frozen V5.3 RCS20 range-anchor CSV must contain food, value1, ",
    "value2 and explicit range_direction. Missing: ",
    paste(setdiff(required_anchor_cols, names(range_anchor_raw)), collapse = ", "),
    ". Direction will not be guessed from value2 because both low- and ",
    "high-intake one-sided rules can have a missing upper bound."
  )
}

anchor_direction_import <- trimws(as.character(range_anchor_raw$range_direction))
if (anyNA(anchor_direction_import) ||
    any(!anchor_direction_import %in% c("low", "middle", "high", "point"))) {
  stop(
    "STOP: Frozen V5.3 RCS20 anchor contains missing/unsupported ",
    "range_direction values."
  )
}

range_anchor <- tibble(
  food = trimws(as.character(range_anchor_raw$food)),
  anchor_value1 = suppressWarnings(as.numeric(range_anchor_raw$value1)),
  anchor_value2 = suppressWarnings(as.numeric(range_anchor_raw$value2)),
  anchor_direction = anchor_direction_import,
  anchor_range_source = if ("range_source" %in% names(range_anchor_raw)) {
    as.character(range_anchor_raw$range_source)
  } else {
    rep("frozen_interval_anchor_source_not_provided", nrow(range_anchor_raw))
  }
)

if (nrow(range_anchor) != 15L ||
    anyNA(range_anchor$food) ||
    any(range_anchor$food == "") ||
    anyDuplicated(range_anchor$food) > 0L ||
    anyNA(range_anchor$anchor_value1) ||
    any(!is.finite(range_anchor$anchor_value1)) ||
    any(!is.na(range_anchor$anchor_value2) &
          !is.finite(range_anchor$anchor_value2) &
          !(range_anchor$anchor_direction == "high" &
              is.infinite(range_anchor$anchor_value2) &
              range_anchor$anchor_value2 > 0)) ||
    any(!range_anchor$anchor_direction %in% c("low", "middle", "high", "point")) ||
    anyNA(range_anchor$anchor_range_source) ||
    any(trimws(range_anchor$anchor_range_source) == "") ||
    any(range_anchor$anchor_direction == "middle" &
          (is.na(range_anchor$anchor_value2) |
           range_anchor$anchor_value2 < range_anchor$anchor_value1)) ||
    any(range_anchor$anchor_direction == "low" &
          !is.na(range_anchor$anchor_value2)) ||
    any(range_anchor$anchor_direction == "high" &
          !is.na(range_anchor$anchor_value2) &
          !(is.infinite(range_anchor$anchor_value2) &
              range_anchor$anchor_value2 > 0)) ||
    any(range_anchor$anchor_direction == "point" &
          !is.na(range_anchor$anchor_value2) &
          abs(range_anchor$anchor_value2 - range_anchor$anchor_value1) > 1e-8)) {
  stop("STOP: Invalid, duplicated or non-15-component range-anchor CSV.")
}

if (!setequal(range_anchor$food, FOOD_MAP$food)) {
  stop(
    "STOP: Range-anchor component set differs from the frozen 15. Missing: ",
    paste(setdiff(FOOD_MAP$food, range_anchor$food), collapse = ", "),
    ". Extra: ",
    paste(setdiff(range_anchor$food, FOOD_MAP$food), collapse = ", ")
  )
}

range_anchor <- range_anchor[
  match(FOOD_MAP$food, range_anchor$food), , drop = FALSE
]
rownames(range_anchor) <- NULL

# interval repair is performed HERE, in Step_04_1_1, before binary
# scoring and before the separate Step_04_2_1 beta-weight workflow.
# Rule: frozen V5.3 RCS20 20% plateau intersected with Model 2 lowest-risk
# category for Potatoes and Red meats only. No Cox P-value search is performed.
if (!file.exists(FROZEN_RCS_DIAG_FILE)) {
  stop("STOP: Frozen V5.3 RCS20 diagnostic file is missing: ", FROZEN_RCS_DIAG_FILE)
}
frozen_rcs_diag <- read.csv(
  FROZEN_RCS_DIAG_FILE, stringsAsFactors = FALSE,
  check.names = FALSE, na.strings = c("", "NA")
)
required_frozen_diag <- c(
  "food", "plateau_20_low", "plateau_20_high",
  "category_range_low", "category_range_high"
)
missing_frozen_diag <- setdiff(required_frozen_diag, names(frozen_rcs_diag))
if (length(missing_frozen_diag) > 0L) {
  stop("STOP: Frozen RCS20 diagnostics missing: ",
       paste(missing_frozen_diag, collapse = ", "))
}
repair_foods <- c("potatoes", "red_meats")
repair_diag <- frozen_rcs_diag[
  match(repair_foods, as.character(frozen_rcs_diag$food)), , drop = FALSE
]
if (nrow(repair_diag) != 2L || anyNA(repair_diag$food)) {
  stop("STOP: Could not identify Potatoes and Red meats in frozen RCS20 diagnostics.")
}
repair_numeric <- c(
  "plateau_20_low", "plateau_20_high",
  "category_range_low", "category_range_high"
)
for (nm in repair_numeric) repair_diag[[nm]] <- as.numeric(repair_diag[[nm]])
if (any(!is.finite(as.matrix(repair_diag[, repair_numeric, drop = FALSE])))) {
  stop("STOP: Non-finite frozen RCS20/Model 2 boundaries for the two repair foods.")
}
repair_intervals <- repair_diag %>%
  transmute(
    food = as.character(food),
    RCS20_plateau_low = plateau_20_low,
    RCS20_plateau_high = plateau_20_high,
    Model2_lowest_risk_category_low = category_range_low,
    Model2_lowest_risk_category_high = category_range_high,
    final_interval_low = pmax(plateau_20_low, category_range_low),
    final_interval_high = pmin(plateau_20_high, category_range_high),
    derivation = "Intersection of frozen V5.3 RCS20 20% plateau and Model 2 lowest-risk category"
  )
if (any(repair_intervals$final_interval_high <= repair_intervals$final_interval_low)) {
  stop("STOP: interval intersection is empty for Potatoes or Red meats.")
}
expected_repair <- tibble(
  food = c("potatoes", "red_meats"),
  expected_low = c(1.00, 1.00),
  expected_high = c(1.50, 1.75)
) %>%
  left_join(repair_intervals, by = "food")
if (any(abs(expected_repair$final_interval_low - expected_repair$expected_low) > 1e-6) ||
    any(abs(expected_repair$final_interval_high - expected_repair$expected_high) > 1e-6)) {
  print(expected_repair)
  stop("STOP: repair intervals differ from frozen expected values.")
}
for (j in seq_len(nrow(repair_intervals))) {
  k <- match(repair_intervals$food[j], range_anchor$food)
  range_anchor$anchor_value1[k] <- repair_intervals$final_interval_low[j]
  range_anchor$anchor_value2[k] <- repair_intervals$final_interval_high[j]
  range_anchor$anchor_direction[k] <- "middle"
  range_anchor$anchor_range_source[k] <-
    "RCS20_20pct_plateau_intersected_with_Model2_lowest_risk_category"
}
range_anchor <- range_anchor %>%
  mutate(
    anchor_source_file = normalizePath(
      RANGE_ANCHOR_CSV, winslash = "/", mustWork = TRUE
    )
  )

write.csv(
  repair_intervals,
  file.path(AUDIT_DIR, "two_component_interval_derivation.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

write.csv(
  range_anchor,
  file.path(AUDIT_DIR, "frozen_15_component_range_anchor_QC.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

model2_sig <- read.csv(
  REF_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_result_cols <- c("group", "model", "variable", "hr_ci")
if (!all(required_result_cols %in% names(model2_sig))) {
  stop(
    "STOP: Model 2 results CSV is missing columns: ",
    paste(setdiff(required_result_cols, names(model2_sig)), collapse = ", ")
  )
}

model2_sig <- model2_sig %>%
  tidyr::fill(group, .direction = "down") %>%
  filter(model == "model2")

csv_groups <- unique(model2_sig$group)
if (!setequal(csv_groups, FOOD_MAP$group)) {
  stop(
    "STOP: The Model 2 significant-results CSV does not contain the exact ",
    "15-component set. Missing: ",
    paste(setdiff(FOOD_MAP$group, csv_groups), collapse = "; "),
    ". Extra: ",
    paste(setdiff(csv_groups, FOOD_MAP$group), collapse = "; ")
  )
}

reference_rows <- model2_sig %>%
  filter(group %in% FOOD_MAP$group, hr_ci == "Ref.") %>%
  select(group, variable, hr_ci)

reference_counts <- reference_rows %>%
  count(group, name = "n_reference_rows")

if (nrow(reference_counts) != 15L ||
    any(reference_counts$n_reference_rows != 1L)) {
  stop(
    "STOP: Each of the 15 components must have exactly one Model 2 Ref. row."
  )
}

REFERENCE_SOURCE <- FOOD_MAP %>%
  left_join(reference_rows, by = "group") %>%
  rename(reference_category_label = variable)

if (nrow(REFERENCE_SOURCE) != 15L ||
    anyNA(REFERENCE_SOURCE$reference_category_label)) {
  stop("STOP: Could not map the 15 lowest-risk reference categories.")
}

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

format_p_value <- function(p) {
  if (length(p) != 1L || is.na(p) || !is.finite(p)) return("= NA")

  # Display all P values to three decimal places. Add an asterisk when the
  # unrounded P value is < 0.05, so a value such as 0.049934 is displayed as
  # "= 0.050*" without losing its nominal significance status.
  if (p < 0.001) return("< 0.001*")

  paste0(
    "= ", sprintf("%.3f", p),
    if (p < 0.05) "*" else ""
  )
}

fmt_num <- function(x) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
  if (abs(x) < 1e-10) return("0")
  formatC(x, format = "f", digits = 2)
}

# Percentiles are printed to three decimals so the plotted P10/P50/P90 labels
# preserve more of the empirical quantile precision than ordinary range labels.
fmt_quantile <- function(x) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
  if (abs(x) < 1e-10) return("0")
  formatC(x, format = "f", digits = 3)
}

parse_leading_numeric <- function(x) {
  pattern <- "^[[:space:]]*[-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?"
  out <- regmatches(
    as.character(x),
    regexpr(pattern, as.character(x), perl = TRUE)
  )
  out[out == ""] <- NA_character_
  suppressWarnings(as.numeric(out))
}

# Choose three distinct knots. The primary rule is full-sample P10/P50/P90.
# If ties caused by a large zero-intake mass collapse those quantiles, use
# zero plus the positive-consumer P50 and P90. If this is still impossible,
# use three distinct observed support values and record that fallback.
choose_three_knots <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 100L || length(unique(x)) < 3L) {
    stop("STOP: Fewer than three distinct exposure values.")
  }

  q_full <- as.numeric(stats::quantile(
    x, probs = c(0.10, 0.50, 0.90),
    na.rm = TRUE, names = FALSE, type = 7
  ))
  names(q_full) <- c("P10", "P50", "P90")

  if (length(unique(q_full)) == 3L && all(diff(q_full) > 0)) {
    return(list(
      knots = q_full,
      method = "full_training_sample_P10_P50_P90",
      positive_quantiles = c(NA_real_, NA_real_, NA_real_)
    ))
  }

  zero_prop <- mean(x <= 1e-12)
  x_positive <- x[x > 1e-12]

  if (zero_prop > 0 && length(unique(x_positive)) >= 2L) {
    q_positive <- as.numeric(stats::quantile(
      x_positive, probs = c(0.10, 0.50, 0.90),
      na.rm = TRUE, names = FALSE, type = 7
    ))
    names(q_positive) <- c("P10_positive", "P50_positive", "P90_positive")

    zero_anchored <- c(0, q_positive["P50_positive"], q_positive["P90_positive"])
    if (length(unique(zero_anchored)) == 3L &&
        all(is.finite(zero_anchored)) &&
        all(diff(zero_anchored) > 0)) {
      return(list(
        knots = unname(zero_anchored),
        method = "zero_anchored_plus_positive_consumer_P50_P90",
        positive_quantiles = q_positive
      ))
    }
  }

  # Last-resort fallback for heavily tied/discrete non-zero distributions.
  # This is not silently described as P10/P50/P90.
  ux <- sort(unique(x))
  if (length(ux) >= 3L) {
    idx <- unique(pmax(
      1L,
      pmin(length(ux), as.integer(round(
        c(0.10, 0.50, 0.90) * (length(ux) - 1L)
      ) + 1L))
    ))
    if (length(idx) == 3L) {
      fallback_knots <- ux[idx]
      if (all(diff(fallback_knots) > 0)) {
        return(list(
          knots = fallback_knots,
          method = "distinct_observed_support_fallback",
          positive_quantiles = c(NA_real_, NA_real_, NA_real_)
        ))
      }
    }
  }

  stop(
    "STOP: Three distinct RCS knots cannot be defined without arbitrary ",
    "extrapolation. Inspect the exposure distribution."
  )
}

linear_crossing_x <- function(x1, y1, x2, y2, y_target) {
  if (any(!is.finite(c(x1, y1, x2, y2, y_target)))) return(NA_real_)
  if (abs(y2 - y1) < .Machine$double.eps^0.5) return(mean(c(x1, x2)))
  x1 + (y_target - y1) * (x2 - x1) / (y2 - y1)
}

get_minimum_risk_plateau <- function(x, y, tolerance) {
  ok <- is.finite(x) & is.finite(y) & y > 0
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 5L) {
    return(c(lower = NA_real_, upper = NA_real_, threshold = NA_real_))
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  i_min <- which.min(y)
  y_min <- y[i_min]
  threshold <- y_min * (1 + tolerance)

  i_left <- i_min
  while (i_left > 1L && y[i_left - 1L] <= threshold) {
    i_left <- i_left - 1L
  }

  i_right <- i_min
  while (i_right < length(y) && y[i_right + 1L] <= threshold) {
    i_right <- i_right + 1L
  }

  lower <- if (i_left == 1L) {
    x[1]
  } else {
    linear_crossing_x(
      x[i_left - 1L], y[i_left - 1L],
      x[i_left], y[i_left], threshold
    )
  }

  upper <- if (i_right == length(y)) {
    x[length(x)]
  } else {
    linear_crossing_x(
      x[i_right], y[i_right],
      x[i_right + 1L], y[i_right + 1L],
      threshold
    )
  }

  c(lower = lower, upper = upper, threshold = threshold)
}

classify_rcs_shape <- function(pred, food_var, p_overall, p_non_linear) {
  pred <- pred %>%
    filter(is.finite(.data[[food_var]]), is.finite(yhat), yhat > 0) %>%
    arrange(.data[[food_var]])
  x <- pred[[food_var]]
  y <- pred$yhat
  if (length(y) < 5L) {
    return(list(
      curve_shape = "insufficient_prediction_points",
      raw_curve_shape = "insufficient_prediction_points",
      shape_evidence = "not_assessable",
      x_min = NA_real_, hr_min = NA_real_,
      left_tail_ratio = NA_real_, right_tail_ratio = NA_real_,
      nonlinear_supported = FALSE,
      overall_supported = FALSE,
      range_eligible = FALSE
    ))
  }

  i_min <- which.min(y)
  i_max <- which.max(y)
  x_min <- x[i_min]
  hr_min <- y[i_min]
  left_ratio <- y[1] / hr_min
  right_ratio <- y[length(y)] / hr_min
  rise_cut <- 1.02

  raw_shape <- if (i_min > 1L && i_min < length(y) &&
                   left_ratio > rise_cut && right_ratio > rise_cut) {
    "U-shaped"
  } else if (i_min == 1L && right_ratio > rise_cut) {
    "J-shaped_increasing-risk"
  } else if (i_min == length(y) && left_ratio > rise_cut) {
    "reverse-J-shaped_decreasing-risk"
  } else if (i_max > 1L && i_max < length(y) &&
             y[1] < y[i_max] / rise_cut &&
             y[length(y)] < y[i_max] / rise_cut) {
    "inverted-U-shaped"
  } else if (y[length(y)] > y[1] * rise_cut) {
    "monotonic-increasing-risk"
  } else if (y[1] > y[length(y)] * rise_cut) {
    "monotonic-decreasing-risk"
  } else {
    "flat-or-complex"
  }

  nonlinear_supported <- is.finite(p_non_linear) && p_non_linear < 0.05
  overall_supported <- is.finite(p_overall) && p_overall < 0.05
  nonlinear_shapes <- c(
    "U-shaped", "J-shaped_increasing-risk",
    "reverse-J-shaped_decreasing-risk", "inverted-U-shaped"
  )

  # A significant overall test does not establish a J/reverse-J shape.
  # When nonlinear evidence is absent, label the supported overall trend by
  # its observed endpoint direction; retain the raw geometry separately.
  curve_shape <- raw_shape
  if (raw_shape %in% nonlinear_shapes && !nonlinear_supported) {
    if (overall_supported && y[length(y)] > y[1] * (1 + 1e-8)) {
      curve_shape <- "monotonic-increasing-risk"
    } else if (overall_supported && y[1] > y[length(y)] * (1 + 1e-8)) {
      curve_shape <- "monotonic-decreasing-risk"
    } else {
      curve_shape <- "flat-or-complex"
    }
  }

  shape_supported <- if (curve_shape %in% nonlinear_shapes) {
    nonlinear_supported
  } else if (curve_shape %in% c(
    "monotonic-increasing-risk", "monotonic-decreasing-risk"
  )) {
    overall_supported
  } else {
    FALSE
  }

  evidence <- if (curve_shape %in% nonlinear_shapes && shape_supported) {
    paste0(curve_shape, "_with_non_linear_evidence")
  } else if (curve_shape %in% c(
    "monotonic-increasing-risk", "monotonic-decreasing-risk"
  ) && shape_supported) {
    paste0(curve_shape, "_with_overall_association_no_claim_of_nonlinearity")
  } else {
    paste0("not_confirmed;_raw_geometry=", raw_shape)
  }

  list(
    curve_shape = curve_shape,
    raw_curve_shape = raw_shape,
    shape_evidence = evidence,
    x_min = x_min,
    hr_min = hr_min,
    left_tail_ratio = left_ratio,
    right_tail_ratio = right_ratio,
    nonlinear_supported = nonlinear_supported,
    overall_supported = overall_supported,
    range_eligible = shape_supported
  )
}

# Category direction is determined from the observed ordering of the intake
# categories in the Step_03_1_0 training data, not from alphabetical labels.
get_category_interval <- function(data, food_var, category_var, ref_label) {
  x <- data[[food_var]]
  cat <- as.character(data[[category_var]])

  group_stats <- tibble(
    category = cat,
    intake = x
  ) %>%
    filter(!is.na(category), is.finite(intake)) %>%
    group_by(category) %>%
    summarise(
      category_n = n(),
      category_median = median(intake),
      category_min = min(intake),
      category_max = max(intake),
      .groups = "drop"
    )

  # Calculate category-specific event counts using explicit row masks.
  group_stats$category_events <- vapply(
    group_stats$category,
    function(g) {
      sum(data$depression_status[
        !is.na(cat) & cat == g
      ] == 1, na.rm = TRUE)
    },
    numeric(1)
  )

  group_stats <- group_stats %>%
    arrange(category_median, category_min, category_max)

  if (!ref_label %in% group_stats$category) {
    stop(
      "STOP: Lowest-risk category label from Step_03_1_0 is not found in ",
      category_var, ": ", ref_label,
      ". Available labels: ",
      paste(group_stats$category, collapse = " | ")
    )
  }

  ref_row <- group_stats %>%
    filter(category == ref_label)

  if (nrow(ref_row) != 1L || ref_row$category_n[1] < 1L) {
    stop("STOP: Invalid lowest-risk category data for ", food_var)
  }

  overall_min <- min(x, na.rm = TRUE)
  overall_max <- max(x, na.rm = TRUE)
  span <- overall_max - overall_min
  tol <- max(1e-10, abs(span) * 1e-8)

  ref_min <- ref_row$category_min[1]
  ref_max <- ref_row$category_max[1]
  ref_median <- median(
    x[!is.na(cat) & cat == ref_label],
    na.rm = TRUE
  )

  ref_rank <- match(ref_label, group_stats$category)
  n_categories <- nrow(group_stats)

  if (abs(ref_max - ref_min) <= tol) {
    direction <- "point"
  } else if (ref_rank == 1L) {
    direction <- "low"
  } else if (ref_rank == n_categories) {
    direction <- "high"
  } else {
    direction <- "middle"
  }

  if (direction == "low") {
    value1 <- ref_max
    value2 <- NA_real_
  } else if (direction == "high") {
    value1 <- ref_min
    value2 <- NA_real_
  } else {
    value1 <- ref_min
    value2 <- ref_max
  }

  list(
    category_stats = group_stats,
    reference_category = ref_label,
    reference_category_rank = ref_rank,
    n_categories = n_categories,
    reference_category_n = ref_row$category_n[1],
    reference_category_events = ref_row$category_events[1],
    reference_category_median = ref_median,
    reference_category_min = ref_min,
    reference_category_max = ref_max,
    overall_min = overall_min,
    overall_max = overall_max,
    range_direction = direction,
    value1 = value1,
    value2 = value2
  )
}

# -----------------------------------------------------------------------------
# 3B. Conditional bootstrap of spline plateau boundaries, and tail-support QC
# -----------------------------------------------------------------------------

summarise_bootstrap_cutoffs <- function(fit, food_var, pred_support,
                                         reference_value, category_low,
                                         category_high, category_direction,
                                         main_rule, model_index) {
  grid <- as.numeric(pred_support[[food_var]])
  if (length(grid) < 40L || any(!is.finite(grid)) ||
      any(diff(grid) <= 0)) {
    stop("STOP: Invalid bootstrap dose grid: ", food_var)
  }

  tolerance_tags <- sprintf(
    "%02d", round(BOOTSTRAP_RANGE_TOLERANCES * 100)
  )
  if (!identical(tolerance_tags, c("05", "10", "20"))) {
    stop("STOP: Parallel bootstrap requires exactly 5%, 10% and 20% tolerances.")
  }

  # contrast.rms returns model-matrix contrasts, with all covariates held
  # fixed so they cancel. The exact spline basis/knots from the fitted rms
  # model is reused. Never reconstruct the RCS basis by hand.
  a <- stats::setNames(list(grid), food_var)
  b <- stats::setNames(list(reference_value), food_var)
  contrast_grid <- rms::contrast(fit, a = a, b = b, expand = FALSE)
  X <- as.matrix(contrast_grid$X)
  if (nrow(X) != length(grid)) {
    stop("STOP: Contrast grid size mismatch for ", food_var)
  }

  boot <- rms::bootcov(
    fit, B = BOOTSTRAP_B, coef.reps = TRUE,
    loglik = FALSE, pr = FALSE,
    seed = as.integer(BOOTSTRAP_SEED + model_index)
  )
  bc <- as.matrix(boot$boot.Coef)
  fit_coef <- stats::coef(fit)
  boot_coef <- stats::coef(boot)

  # ---------------------------------------------------------------------------
  # Robust coefficient alignment for rms::contrast() + rms::bootcov().
  # ---------------------------------------------------------------------------
  model_coef_names <- names(fit_coef)
  boot_model_coef_names <- names(boot_coef)
  X_names_original <- colnames(X)
  bc_names_original <- colnames(bc)

  if (is.null(model_coef_names) || anyDuplicated(model_coef_names) ||
      length(fit_coef) != ncol(X) ||
      length(fit_coef) != ncol(bc) ||
      length(boot_model_coef_names) != length(fit_coef)) {
    stop(
      "STOP: Bootstrap coefficient dimensions/names are invalid for ", food_var,
      ". ncol(X)=", ncol(X),
      "; ncol(boot.Coef)=", ncol(bc),
      "; length(coef(fit))=", length(fit_coef),
      "; length(coef(bootcov fit))=", length(boot_model_coef_names), "."
    )
  }

  if (!setequal(boot_model_coef_names, model_coef_names)) {
    stop(
      "STOP: rms::bootcov model coefficient set differs from original fit for ",
      food_var,
      ". Original={", paste(model_coef_names, collapse = " | "),
      "}; bootcov={", paste(boot_model_coef_names, collapse = " | "), "}."
    )
  }

  if (!is.null(bc_names_original) &&
      setequal(bc_names_original, boot_model_coef_names)) {
    bc <- bc[, boot_model_coef_names, drop = FALSE]
  } else if (is.null(bc_names_original)) {
    colnames(bc) <- boot_model_coef_names
  } else {
    stop(
      "STOP: boot.Coef columns cannot be mapped to rms::bootcov coefficients for ",
      food_var,
      ". boot.Coef columns={", paste(bc_names_original, collapse = " | "),
      "}; bootcov coefficient names={",
      paste(boot_model_coef_names, collapse = " | "), "}."
    )
  }

  if (!is.null(X_names_original) &&
      setequal(X_names_original, model_coef_names)) {
    X <- X[, model_coef_names, drop = FALSE]
  }

  contrast_target <- as.numeric(contrast_grid$Contrast)
  contrast_from_X <- as.numeric(X %*% fit_coef)
  if (length(contrast_target) != nrow(X) ||
      length(contrast_from_X) != length(contrast_target) ||
      !isTRUE(all.equal(
        contrast_from_X, contrast_target,
        tolerance = 1e-8, check.attributes = FALSE
      ))) {
    stop(
      "STOP: Contrast design-matrix order could not be validated for ", food_var,
      ". X columns={", paste(X_names_original, collapse = " | "),
      "}; model coefficients={", paste(model_coef_names, collapse = " | "),
      "}; boot.Coef columns={", paste(colnames(bc), collapse = " | "),
      "}; maximum absolute contrast reconstruction error=",
      if (length(contrast_from_X) == length(contrast_target)) {
        format(max(abs(contrast_from_X - contrast_target)), scientific = TRUE)
      } else {
        "not-computable"
      }, "."
    )
  }

  # rms::bootcov() stores coefficient replicates in model coefficient order.
  colnames(X) <- model_coef_names
  if (!identical(colnames(bc), model_coef_names)) {
    bc <- bc[, model_coef_names, drop = FALSE]
  }

  coefficient_alignment_method <- if (
    !is.null(X_names_original) &&
      identical(X_names_original, model_coef_names) &&
      !is.null(bc_names_original) &&
      identical(bc_names_original, model_coef_names)
  ) {
    "exact_name_match"
  } else {
    "validated_by_rms_contrast_estimate_and_bootcov_model_order"
  }
  contrast_reconstruction_max_abs_error <- max(
    abs(contrast_from_X - contrast_target)
  )

  if (nrow(bc) < ceiling(BOOTSTRAP_B * BOOTSTRAP_MIN_SUCCESS_RATE)) {
    stop("STOP: Insufficient successful RCS bootstrap refits for ", food_var,
         ": ", nrow(bc), " of ", BOOTSTRAP_B)
  }

  # Each row is one participant-resampled Model 3 Cox/RCS curve. The SAME
  # bootstrap curve is evaluated at 5%, 10% and 20%; thresholds are not refit
  # separately and do not alter the original fitted model or scoring anchor.
  log_hr <- bc %*% t(X)
  replicate_results <- vector("list", nrow(log_hr))
  bounded_category <- identical(category_direction, "middle") &&
    is.finite(category_low) && is.finite(category_high) &&
    category_high > category_low

  for (j in seq_len(nrow(log_hr))) {
    lp <- as.numeric(log_hr[j, ])
    if (any(!is.finite(lp))) {
      stop("STOP: Non-finite bootstrap log-HR for ", food_var, ", rep ", j)
    }
    y <- exp(pmax(-500, pmin(500, lp - min(lp))))

    plateau_list <- stats::setNames(
      lapply(BOOTSTRAP_RANGE_TOLERANCES, function(tol) {
        get_minimum_risk_plateau(x = grid, y = y, tolerance = tol)
      }),
      tolerance_tags
    )

    overlap_by_tag <- stats::setNames(vector("list", length(tolerance_tags)),
                                      tolerance_tags)
    intersection_by_tag <- stats::setNames(vector("list", length(tolerance_tags)),
                                            tolerance_tags)

    for (tag in tolerance_tags) {
      p <- plateau_list[[tag]]
      overlap <- is.finite(p["lower"]) && is.finite(p["upper"]) &&
        p["lower"] <= category_high && p["upper"] >= category_low
      overlap_by_tag[[tag]] <- isTRUE(overlap)

      int_lower <- NA_real_
      int_upper <- NA_real_
      if (bounded_category && isTRUE(overlap)) {
        int_lower <- max(category_low, unname(p["lower"]))
        int_upper <- min(category_high, unname(p["upper"]))
        if (!is.finite(int_lower) || !is.finite(int_upper) ||
            int_upper <= int_lower) {
          int_lower <- int_upper <- NA_real_
        }
      }
      intersection_by_tag[[tag]] <- c(lower = int_lower, upper = int_upper)
    }

    # Preserve the legacy current-model candidate-boundary bootstrap, which is
    # specifically based on the 5% rule. It remains an audit and never changes
    # the frozen primary scoring interval.
    plateau_05 <- plateau_list[["05"]]
    overlap_05 <- overlap_by_tag[["05"]]
    candidate_lower <- NA_real_
    candidate_upper <- NA_real_
    if (main_rule == "bounded_intersection" && isTRUE(overlap_05)) {
      candidate_lower <- max(category_low, unname(plateau_05["lower"]))
      candidate_upper <- min(category_high, unname(plateau_05["upper"]))
      if (candidate_upper <= candidate_lower) {
        candidate_lower <- candidate_upper <- NA_real_
      }
    }

    endpoint_risk <- c(y[1], y[length(y)])
    geometry_u <- which.min(y) > 1L && which.min(y) < length(y) &&
      all(endpoint_risk / min(y) > 1.02)

    row_values <- list(food = food_var, bootstrap_rep = j)
    for (tag in tolerance_tags) {
      p <- plateau_list[[tag]]
      int <- intersection_by_tag[[tag]]
      row_values[[paste0("plateau_", tag, "_lower")]] <- unname(p["lower"])
      row_values[[paste0("plateau_", tag, "_upper")]] <- unname(p["upper"])
      row_values[[paste0("plateau_", tag, "_threshold")]] <- unname(p["threshold"])
      row_values[[paste0("plateau_", tag, "_overlaps_cox_category")]] <-
        overlap_by_tag[[tag]]
      row_values[[paste0("plateau_", tag, "_intersection_lower")]] <-
        unname(int["lower"])
      row_values[[paste0("plateau_", tag, "_intersection_upper")]] <-
        unname(int["upper"])
    }

    row_values$x_min <- grid[which.min(y)]
    row_values$u_geometry <- geometry_u
    # Backward-compatible alias for the former 5% overlap field.
    row_values$plateau_overlaps_cox_category <- overlap_05
    row_values$conditional_candidate_lower <- unname(candidate_lower)
    row_values$conditional_candidate_upper <- unname(candidate_upper)
    row_values$main_rule <- main_rule
    replicate_results[[j]] <- tibble::as_tibble(row_values)
  }

  d <- dplyr::bind_rows(replicate_results)
  alpha <- (1 - BOOTSTRAP_BOUNDARY_CI_LEVEL) / 2
  ci <- function(v) {
    v <- v[is.finite(v)]
    min_n <- max(30L, floor(0.8 * nrow(d)))
    if (length(v) < min_n) return(c(NA_real_, NA_real_))
    unname(stats::quantile(
      v, probs = c(alpha, 1 - alpha), names = FALSE, type = 7
    ))
  }

  summary_values <- list(
    food = food_var,
    bootstrap_requested = as.integer(BOOTSTRAP_B),
    bootstrap_successful = nrow(d),
    bootstrap_conditional_on_frozen_cox_category_and_rcs_knots = TRUE,
    bootstrap_reference_value_frozen = TRUE,
    bootstrap_rcs_knots_frozen = TRUE,
    bootstrap_prediction_grid_frozen = TRUE,
    bootstrap_thresholds_evaluated = paste(tolerance_tags, collapse = ";"),
    bootstrap_boundary_ci_level = BOOTSTRAP_BOUNDARY_CI_LEVEL,
    bootstrap_boundary_ci_method = BOOTSTRAP_BOUNDARY_CI_METHOD,
    bootstrap_main_rule = main_rule,
    bootstrap_coefficient_alignment_method = coefficient_alignment_method,
    bootstrap_contrast_matrix_columns = as.integer(ncol(X)),
    bootstrap_coefficient_columns = as.integer(ncol(bc)),
    bootstrap_contrast_reconstruction_max_abs_error =
      contrast_reconstruction_max_abs_error,
    bootstrap_x_min_ci_lower = ci(d$x_min)[1],
    bootstrap_x_min_ci_upper = ci(d$x_min)[2]
  )

  for (tag in tolerance_tags) {
    lower_name <- paste0("plateau_", tag, "_lower")
    upper_name <- paste0("plateau_", tag, "_upper")
    overlap_name <- paste0("plateau_", tag, "_overlaps_cox_category")
    int_lower_name <- paste0("plateau_", tag, "_intersection_lower")
    int_upper_name <- paste0("plateau_", tag, "_intersection_upper")

    lower_ci <- ci(d[[lower_name]])
    upper_ci <- ci(d[[upper_name]])
    width <- d[[upper_name]] - d[[lower_name]]
    width_ci <- ci(width)
    int_lower_ci <- ci(d[[int_lower_name]])
    int_upper_ci <- ci(d[[int_upper_name]])
    int_width <- d[[int_upper_name]] - d[[int_lower_name]]
    int_width_ci <- ci(int_width)

    summary_values[[paste0("bootstrap_plateau_", tag, "_lower_ci_lower")]] <- lower_ci[1]
    summary_values[[paste0("bootstrap_plateau_", tag, "_lower_ci_upper")]] <- lower_ci[2]
    summary_values[[paste0("bootstrap_plateau_", tag, "_upper_ci_lower")]] <- upper_ci[1]
    summary_values[[paste0("bootstrap_plateau_", tag, "_upper_ci_upper")]] <- upper_ci[2]
    summary_values[[paste0("bootstrap_plateau_", tag, "_lower_estimable_n")]] <-
      sum(is.finite(d[[lower_name]]))
    summary_values[[paste0("bootstrap_plateau_", tag, "_upper_estimable_n")]] <-
      sum(is.finite(d[[upper_name]]))
    summary_values[[paste0("bootstrap_plateau_", tag, "_lower_estimable_fraction")]] <-
      mean(is.finite(d[[lower_name]]))
    summary_values[[paste0("bootstrap_plateau_", tag, "_upper_estimable_fraction")]] <-
      mean(is.finite(d[[upper_name]]))
    summary_values[[paste0("bootstrap_plateau_", tag, "_width_median")]] <-
      if (any(is.finite(width))) stats::median(width, na.rm = TRUE) else NA_real_
    summary_values[[paste0("bootstrap_plateau_", tag, "_width_ci_lower")]] <- width_ci[1]
    summary_values[[paste0("bootstrap_plateau_", tag, "_width_ci_upper")]] <- width_ci[2]
    summary_values[[paste0("bootstrap_plateau_", tag, "_cox_overlap_fraction")]] <-
      mean(as.logical(d[[overlap_name]]), na.rm = TRUE)
    summary_values[[paste0("bootstrap_plateau_", tag, "_intersection_lower_ci_lower")]] <-
      int_lower_ci[1]
    summary_values[[paste0("bootstrap_plateau_", tag, "_intersection_lower_ci_upper")]] <-
      int_lower_ci[2]
    summary_values[[paste0("bootstrap_plateau_", tag, "_intersection_upper_ci_lower")]] <-
      int_upper_ci[1]
    summary_values[[paste0("bootstrap_plateau_", tag, "_intersection_upper_ci_upper")]] <-
      int_upper_ci[2]
    summary_values[[paste0("bootstrap_plateau_", tag, "_intersection_estimable_fraction")]] <-
      mean(is.finite(d[[int_lower_name]]) & is.finite(d[[int_upper_name]]))
    summary_values[[paste0("bootstrap_plateau_", tag, "_intersection_width_median")]] <-
      if (any(is.finite(int_width))) stats::median(int_width, na.rm = TRUE) else NA_real_
    summary_values[[paste0("bootstrap_plateau_", tag, "_intersection_width_ci_lower")]] <-
      int_width_ci[1]
    summary_values[[paste0("bootstrap_plateau_", tag, "_intersection_width_ci_upper")]] <-
      int_width_ci[2]
  }

  # Legacy 5% candidate-cutpoint fields retained for downstream compatibility.
  summary_values$bootstrap_conditional_cutpoint_lower_ci_lower <-
    ci(d$conditional_candidate_lower)[1]
  summary_values$bootstrap_conditional_cutpoint_lower_ci_upper <-
    ci(d$conditional_candidate_lower)[2]
  summary_values$bootstrap_conditional_cutpoint_upper_ci_lower <-
    ci(d$conditional_candidate_upper)[1]
  summary_values$bootstrap_conditional_cutpoint_upper_ci_upper <-
    ci(d$conditional_candidate_upper)[2]
  summary_values$bootstrap_cox_overlap_fraction <-
    summary_values$bootstrap_plateau_05_cox_overlap_fraction
  summary_values$bootstrap_u_geometry_fraction <- mean(d$u_geometry)
  summary_values$bootstrap_conditional_cutpoint_estimable_fraction <-
    mean(is.finite(d$conditional_candidate_lower))

  list(
    summary = tibble::as_tibble(summary_values),
    replicates = d
  )
}

quantile_tail_support <- function(fit, food_var, x, event) {
  q <- as.numeric(stats::quantile(
    x, probs = c(0.90, 0.95, 0.99),
    na.rm = TRUE, names = FALSE, type = 7
  ))
  names(q) <- c("P90", "P95", "P99")
  bins <- tibble(
    segment = c("below_P90", "P90_to_P95", "P95_to_P99", "P99_and_above"),
    n = c(
      sum(x < q[1]), sum(x >= q[1] & x < q[2]),
      sum(x >= q[2] & x < q[3]), sum(x >= q[3])
    ),
    events = c(
      sum(event[x < q[1]] == 1),
      sum(event[x >= q[1] & x < q[2]] == 1),
      sum(event[x >= q[2] & x < q[3]] == 1),
      sum(event[x >= q[3]] == 1)
    )
  )
  if (sum(bins$n) != length(x) || sum(bins$events) != sum(event == 1)) {
    stop("STOP: P90/P95/P99 tail bin counts do not reconcile: ", food_var)
  }
  n_90_99 <- sum(x >= q[1] & x < q[3])
  events_90_99 <- sum(event[x >= q[1] & x < q[3]] == 1)
  n_99 <- sum(x >= q[3]); events_99 <- sum(event[x >= q[3]] == 1)
  eligible <- n_90_99 >= TAIL_MIN_N_P90_P99 &&
    events_90_99 >= TAIL_MIN_EVENTS_P90_P99 &&
    n_99 >= TAIL_MIN_N_P99_PLUS && events_99 >= TAIL_MIN_EVENTS_P99_PLUS

  # Independent of the plot x-axis: tail contrasts are evaluated at P90/P99
  # using the same fitted model and its original Cox lowest-risk reference.
  tail_contrast <- rms::contrast(
    fit,
    a = stats::setNames(list(q[3]), food_var),
    b = stats::setNames(list(q[1]), food_var)
  )
  contrast <- as.numeric(tail_contrast$Contrast[1])
  se <- as.numeric(tail_contrast$SE[1])
  z <- stats::qnorm(0.975)
  if (!is.finite(contrast) || !is.finite(se)) {
    stop("STOP: Nonfinite P99/P90 RCS contrast: ", food_var)
  }
  tail <- tibble(
    food = food_var,
    empirical_p90 = q[1], empirical_p95 = q[2], empirical_p99 = q[3],
    n_p90_to_p99 = n_90_99, events_p90_to_p99 = events_90_99,
    n_p99_plus = n_99, events_p99_plus = events_99,
    tail_support_meets_descriptive_thresholds = eligible,
    hr_p99_vs_p90 = exp(contrast),
    hr_p99_vs_p90_lower95 = exp(contrast - z * se),
    hr_p99_vs_p90_upper95 = exp(contrast + z * se),
    p99_vs_p90_ci_excludes_1 = (contrast - z * se > 0) ||
      (contrast + z * se < 0),
    tail_inference = if (eligible) {
      "P90-P99 descriptive support met; interpret observational association cautiously"
    } else {
      "Sparse tail or tied quantiles; P99 contrast considered exploratory only"
    }
  )
  bins$food <- food_var
  bins$p90 <- q[1]; bins$p95 <- q[2]; bins$p99 <- q[3]
  list(summary = tail, bins = bins)
}

# -----------------------------------------------------------------------------
# 4. Fit the 15 RCS models and create the original-style plots
# -----------------------------------------------------------------------------

col <- cols4all::c4a("friendly13")
col_line <- cols4all::c4a("friendly13")

# Exact colours requested for the three empirical intake-percentile markers.
PERCENTILE_COLORS <- c(
  P10 = "#77AADD",
  P50 = "#009988",
  P90 = "#9467BD"
)
if (!identical(
  unname(as.character(col[3:5])),
  unname(PERCENTILE_COLORS)
)) {
  stop(
    "STOP: The friendly13 palette indices col[3:5] no longer match the ",
    "requested P10/P50/P90 colours. Check the palette before plotting."
  )
}

range_results <- vector("list", nrow(FOOD_MAP))
model_results <- vector("list", nrow(FOOD_MAP))
plot_list <- vector("list", nrow(FOOD_MAP))
names(plot_list) <- FOOD_MAP$food

reference_audit_list <- vector("list", nrow(FOOD_MAP))
zero_knot_audit_list <- vector("list", nrow(FOOD_MAP))
percentile_plot_audit_list <- vector("list", nrow(FOOD_MAP))
pvalue_annotation_list <- vector("list", nrow(FOOD_MAP))
bootstrap_summary_list <- vector("list", nrow(FOOD_MAP))
bootstrap_replicates_list <- vector("list", nrow(FOOD_MAP))
tail_summary_list <- vector("list", nrow(FOOD_MAP))
tail_bins_list <- vector("list", nrow(FOOD_MAP))

for (i in seq_len(nrow(FOOD_MAP))) {

  food_var <- FOOD_MAP$food[i]
  category_var <- FOOD_MAP$category_var[i]
  display_name <- FOOD_MAP$display_name[i]
  ref_label <- REFERENCE_SOURCE$reference_category_label[
    match(food_var, REFERENCE_SOURCE$food)
  ]
  anchor_row <- range_anchor[
    match(food_var, range_anchor$food), , drop = FALSE
  ]
  if (nrow(anchor_row) != 1L ||
      !is.finite(anchor_row$anchor_value1[1])) {
    stop("STOP: Missing/invalid frozen range anchor for ", food_var)
  }
  anchor_value1 <- as.numeric(anchor_row$anchor_value1[1])
  anchor_value2 <- as.numeric(anchor_row$anchor_value2[1])

  message(
    "\n[RCS ", i, "/15] ", food_var,
    " | lowest-risk category = ", ref_label
  )

  missing_vars <- setdiff(
    c(
      "depression_time", "depression_status",
      food_var, category_var, MODEL3_COVARIATES
    ),
    names(ukb_data_depression_rcs)
  )
  if (length(missing_vars) > 0L) {
    stop(
      "STOP: Missing variables for ", food_var, ": ",
      paste(missing_vars, collapse = ", ")
    )
  }

  required_vars <- unique(c(
    "depression_time", "depression_status",
    food_var, category_var, MODEL3_COVARIATES
  ))

  keep <- stats::complete.cases(
    ukb_data_depression_rcs[, required_vars, drop = FALSE]
  )
  fit_data <- ukb_data_depression_rcs[keep, , drop = FALSE]

  if (nrow(fit_data) < 100L ||
      sum(fit_data$depression_status == 1, na.rm = TRUE) < 20L) {
    stop("STOP: Insufficient complete-case training data/events for ", food_var)
  }

  x <- as.numeric(fit_data[[food_var]])
  if (any(!is.finite(x))) {
    stop("STOP: Non-finite food intake values remain for ", food_var)
  }
  if (any(x < -1e-8)) {
    stop("STOP: Negative food intake values detected for ", food_var)
  }
  if (length(unique(x)) < 3L) {
    stop("STOP: Fewer than three unique intake values for ", food_var)
  }

  # The selected reference category is taken directly from Step_03_1_0.
  # Its median is recomputed in the exact RCS complete-case training sample.
  category_info <- get_category_interval(
    data = fit_data,
    food_var = food_var,
    category_var = category_var,
    ref_label = ref_label
  )

  ref_value <- category_info$reference_category_median
  ref_label_median <- parse_leading_numeric(ref_label)

  if (!is.finite(ref_value) ||
      ref_value < min(x) - 1e-8 ||
      ref_value > max(x) + 1e-8) {
    stop("STOP: Invalid lowest-risk category median for ", food_var)
  }

  if (length(ref_label_median) == 1L &&
      is.finite(ref_label_median) &&
      abs(ref_label_median - ref_value) >
        max(0.03, abs(ref_value) * 0.03)) {
    warning(
      "Rounded median in the Step_03_1_0 category label differs from the ",
      "complete-case training median for ", food_var,
      ". The exact training median is used as the RCS reference."
    )
  }

  # Full training-sample percentiles are reported for every component.
  p_values <- as.numeric(stats::quantile(
    x, probs = c(0.10, 0.50, 0.90),
    na.rm = TRUE, names = FALSE, type = 7
  ))
  names(p_values) <- c("P10", "P50", "P90")

  zero_prop <- mean(x <= 1e-12)
  knot_info <- choose_three_knots(x)
  knots <- as.numeric(knot_info$knots)

  if (length(knots) != 3L ||
      any(!is.finite(knots)) ||
      any(diff(knots) <= 0)) {
    stop("STOP: Invalid three-knot specification for ", food_var)
  }

  # P10/P90 define the central prediction support, matching the original
  # Predict-based display range; the fitted model still uses all training data.
  support_values <- as.numeric(stats::quantile(
    x,
    probs = PREDICTION_SUPPORT_PROBS,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
  support_low <- support_values[1]
  support_high <- support_values[2]
  support_method <- "P10_P90"

  # If a very large zero mass makes P10 and P90 identical, retain the same
  # zero-aware analysis but expand the prediction support to P1-P99.
  if (!is.finite(support_low) ||
      !is.finite(support_high) ||
      support_high <= support_low) {
    support_values <- as.numeric(stats::quantile(
      x,
      probs = c(0.01, 0.99),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    ))
    support_low <- support_values[1]
    support_high <- support_values[2]
    support_method <- "P1_P99_zero_mass_fallback"
  }

  if (!is.finite(support_low) ||
      !is.finite(support_high) ||
      support_high <= support_low) {
    stop(
      "STOP: Exposure has no usable central prediction support for ",
      food_var, ". A reliable continuous RCS cannot be displayed."
    )
  }

  # Use a food-specific datadist so the reference point is the median of the
  # Step_03_1_0 lowest-risk category, not automatically zero/minimum intake.
  ddist <- rms::datadist(
    fit_data[, unique(c(food_var, MODEL3_COVARIATES)), drop = FALSE]
  )
  ddist$limits["Adjust to", food_var] <- ref_value
  options(datadist = "ddist")

  knot_text <- paste(sprintf("%.15g", knots), collapse = ", ")

  rcs_formula <- stats::as.formula(paste0(
    "survival::Surv(depression_time, depression_status == 1) ~ ",
    "rcs(", food_var, ", c(", knot_text, ")) + ",
    paste(MODEL3_COVARIATES, collapse = " + ")
  ))

  fit_rcs3 <- rms::cph(
    rcs_formula,
    data = fit_data,
    x = TRUE,
    y = TRUE,
    surv = TRUE,
    method = "efron"
  )

  # Do not use fit_rcs3$n or fit_rcs3$nevent here:
  # rms::cph may store n as a two-element vector (observations, events),
  # while $nevent is not a reliable scalar field. Derive both directly from
  # the exact complete-case data passed to cph().
  model_n_actual <- as.integer(nrow(fit_data))
  model_events_actual <- as.integer(sum(
    fit_data$depression_status == 1,
    na.rm = TRUE
  ))

  if (length(model_n_actual) != 1L ||
      length(model_events_actual) != 1L ||
      model_n_actual < 1L ||
      model_events_actual < 0L ||
      model_events_actual > model_n_actual) {
    stop("STOP: Invalid model N/events derived from fit_data for ", food_var)
  }

  anova_rcs3 <- as.matrix(stats::anova(fit_rcs3))
  if (nrow(anova_rcs3) < 2L || ncol(anova_rcs3) < 3L) {
    stop("STOP: Unexpected rms::anova output for ", food_var)
  }

  p_overall <- as.numeric(anova_rcs3[1, 3])
  p_non_linear <- as.numeric(anova_rcs3[2, 3])

  # Prediction grid: dense P10-P90 support plus exact P10/P50/P90, knots and
  # the frozen lowest-risk-category median. The graph styling itself is kept.
  grid_low <- min(
    support_low, ref_value, category_info$value1,
    anchor_value1, anchor_value2, na.rm = TRUE
  )
  grid_high <- max(
    support_high,
    ref_value,
    category_info$value1,
    category_info$value2,
    anchor_value1,
    anchor_value2,
    na.rm = TRUE
  )

  pred_grid <- sort(unique(c(
    seq(support_low, support_high, length.out = PREDICTION_GRID_N),
    seq(grid_low, grid_high, length.out = PREDICTION_GRID_N),
    p_values,
    knots,
    ref_value,
    category_info$value1,
    category_info$value2,
    anchor_value1,
    anchor_value2
  )))
  pred_grid <- pred_grid[
    is.finite(pred_grid) &
      pred_grid >= min(x) - 1e-8 &
      pred_grid <= max(x) + 1e-8
  ]

  pred_args <- c(
    list(fit_rcs3),
    stats::setNames(list(pred_grid), food_var),
    list(
      fun = exp,
      type = "predictions",
      ref.zero = TRUE,
      conf.int = 0.95
    )
  )
  fit_rcs3_pre <- as.data.frame(do.call(rms::Predict, pred_args))

  if (!food_var %in% names(fit_rcs3_pre) ||
      !all(c("yhat", "lower", "upper") %in% names(fit_rcs3_pre))) {
    stop("STOP: Unexpected rms::Predict output for ", food_var)
  }

  fit_rcs3_pre <- fit_rcs3_pre %>%
    arrange(.data[[food_var]]) %>%
    distinct(.data[[food_var]], .keep_all = TRUE)

  # Hard QC: HR at the exact lowest-risk-category median must be 1.
  ref_args <- c(
    list(fit_rcs3),
    stats::setNames(list(ref_value), food_var),
    list(
      fun = exp,
      type = "predictions",
      ref.zero = TRUE,
      conf.int = 0.95
    )
  )
  ref_prediction <- as.data.frame(do.call(rms::Predict, ref_args))

  if (nrow(ref_prediction) != 1L ||
      !is.finite(ref_prediction$yhat[1]) ||
      abs(ref_prediction$yhat[1] - 1) > 1e-6) {
    stop(
      "STOP: RCS HR at the Step_03_1_0 lowest-risk-category median is not 1 ",
      "for ", food_var, "."
    )
  }

  if (!any(abs(fit_rcs3_pre[[food_var]] - ref_value) < 1e-10)) {
    fit_rcs3_pre <- bind_rows(fit_rcs3_pre, ref_prediction) %>%
      arrange(.data[[food_var]])
  }

  # RCS shape and plateau are estimated for current-model validation; the
  # frozen scoring anchor is not altered by current P values or CI patterns.
  pred_support <- fit_rcs3_pre %>%
    filter(
      .data[[food_var]] >= support_low - 1e-10,
      .data[[food_var]] <= support_high + 1e-10
    ) %>%
    arrange(.data[[food_var]])

  shape_info <- classify_rcs_shape(
    pred = pred_support,
    food_var = food_var,
    p_overall = p_overall,
    p_non_linear = p_non_linear
  )

  plateau_list <- lapply(RANGE_TOLERANCES, function(tol) {
    get_minimum_risk_plateau(
      x = pred_support[[food_var]],
      y = pred_support$yhat,
      tolerance = tol
    )
  })
  names(plateau_list) <- paste0(
    "tol_", sprintf("%02d", round(RANGE_TOLERANCES * 100))
  )

  plateau_05 <- plateau_list[["tol_05"]]
  plateau_10 <- plateau_list[["tol_10"]]
  plateau_20 <- plateau_list[["tol_20"]]

  category_low <- category_info$reference_category_min
  category_high <- category_info$reference_category_max
  category_direction <- category_info$range_direction

  plateau_overlap <- is.finite(plateau_05["lower"]) &&
    is.finite(plateau_05["upper"]) &&
    plateau_05["upper"] >= category_low &&
    plateau_05["lower"] <= category_high

  plateau_overlap_low <- if (plateau_overlap) {
    max(category_low, unname(plateau_05["lower"]))
  } else {
    NA_real_
  }
  plateau_overlap_high <- if (plateau_overlap) {
    min(category_high, unname(plateau_05["upper"]))
  } else {
    NA_real_
  }

  # Check whether the categorical lowest-risk group and RCS curve have
  # directionally compatible patterns. This is an audit flag, not a filter
  # that silently removes a food component.
  direction_consistency <- case_when(
    # A single-value Cox category is retained as a point rule and does not
    # require spline refinement.
    category_direction == "point" &&
      category_low == category_high ~ "zero_or_single_value_category",
    # Statistical support must be checked BEFORE directional compatibility.
    # In V3.2 this condition appeared after the direction-specific branches,
    # allowing an unsupported reverse-J/U-shaped appearance to pass the gate.
    !isTRUE(shape_info$range_eligible) ~ "RCS_not_confirmatory",
    category_direction == "low" &&
      shape_info$curve_shape %in% c(
        "J-shaped_increasing-risk", "monotonic-increasing-risk"
      ) ~ "consistent",
    category_direction == "high" &&
      shape_info$curve_shape %in% c(
        "reverse-J-shaped_decreasing-risk", "monotonic-decreasing-risk"
      ) ~ "consistent",
    category_direction == "middle" &&
      shape_info$curve_shape == "U-shaped" ~ "consistent",
    TRUE ~ "review_category_RCS_discordance"
  )

  # ---------------------------------------------------------------------------
  # Candidate interval: reproduce the current Cox-category + RCS-plateau rule
  # for transparent comparison. This is an audit candidate, not the final
  # scoring interval in V1.3.
  # ---------------------------------------------------------------------------
  category_value1 <- category_info$value1
  category_value2 <- category_info$value2
  candidate_value1 <- category_value1
  candidate_value2 <- category_value2
  candidate_range_direction <- category_direction
  candidate_range_source <- "Step_03_1_0_lowest_risk_category_only"
  candidate_refined_by_RCS_05 <- FALSE

  eligible_for_plateau_refinement <-
    isTRUE(shape_info$nonlinear_supported) &&
    shape_info$curve_shape %in% c(
      "U-shaped", "J-shaped_increasing-risk",
      "reverse-J-shaped_decreasing-risk"
    )

  if (isTRUE(eligible_for_plateau_refinement) &&
      direction_consistency == "consistent" &&
      isTRUE(plateau_overlap)) {
    candidate_lower <- max(
      category_low,
      unname(plateau_05["lower"])
    )
    candidate_upper <- min(
      category_high,
      unname(plateau_05["upper"])
    )

    if (candidate_range_direction == "low") {
      if (is.finite(candidate_lower) &&
          is.finite(candidate_upper) &&
          candidate_upper > candidate_lower) {
        lower_tol <- max(
          1e-10,
          abs(category_info$overall_max - category_info$overall_min) * 1e-8
        )
        if (candidate_lower <= category_info$overall_min + lower_tol) {
          candidate_value1 <- candidate_upper
          candidate_value2 <- NA_real_
          candidate_range_direction <- "low"
        } else {
          candidate_value1 <- candidate_lower
          candidate_value2 <- candidate_upper
          candidate_range_direction <- "middle"
        }
        candidate_refined_by_RCS_05 <- TRUE
      }
    } else if (candidate_range_direction == "high") {
      if (is.finite(candidate_lower) &&
          candidate_lower <= category_high) {
        candidate_value1 <- candidate_lower
        candidate_value2 <- NA_real_
        candidate_range_direction <- "high"
        candidate_refined_by_RCS_05 <- TRUE
      }
    } else if (candidate_range_direction == "middle") {
      if (is.finite(candidate_lower) &&
          is.finite(candidate_upper) &&
          candidate_upper > candidate_lower) {
        candidate_value1 <- candidate_lower
        candidate_value2 <- candidate_upper
        candidate_range_direction <- "middle"
        candidate_refined_by_RCS_05 <- TRUE
      }
    }
  }

  if (candidate_refined_by_RCS_05) {
    if (candidate_range_direction == "high") {
      candidate_range_source <-
        "lowest_risk_category_RCS_lower_bound_refined_open_upper_tail"
    } else {
      candidate_range_source <-
        "lowest_risk_category_intersect_RCS_5pct_plateau"
    }
  }

  candidate_score_lower <- case_when(
    candidate_range_direction == "low" ~ category_info$overall_min,
    candidate_range_direction == "high" ~ candidate_value1,
    candidate_range_direction %in% c("middle", "point") ~ candidate_value1,
    TRUE ~ NA_real_
  )
  candidate_score_upper <- case_when(
    candidate_range_direction == "low" ~ candidate_value1,
    candidate_range_direction == "high" ~ Inf,
    candidate_range_direction %in% c("middle", "point") ~ candidate_value2,
    TRUE ~ NA_real_
  )
  candidate_score_rule <- case_when(
    candidate_range_direction == "low" ~
      paste0("x <= ", fmt_num(candidate_value1)),
    candidate_range_direction == "high" ~
      paste0("x >= ", fmt_num(candidate_value1)),
    candidate_range_direction == "middle" ~ paste0(
      fmt_num(candidate_value1), " <= x <= ", fmt_num(candidate_value2)
    ),
    candidate_range_direction == "point" ~
      paste0("x = ", fmt_num(candidate_value1)),
    TRUE ~ "not_defined"
  )

  # ---------------------------------------------------------------------------
  # Primary interval: retain the range_direction/value1/value2 derived above.
  # Low/high one-sided intervals and point intervals are preserved explicitly.
  # Open tails are represented by observed training support only for overlap
  # and plotting audits; the scoring rule itself remains open-ended.
  # Current Cox/RCS results do not move cut-points or optimize Cox P values.
  # ---------------------------------------------------------------------------
  anchor_range_direction <- as.character(anchor_row$anchor_direction[1])

  anchor_effective_lower <- if (anchor_range_direction == "low") {
    category_info$overall_min
  } else {
    anchor_value1
  }
  anchor_effective_upper <- if (anchor_range_direction == "low") {
    anchor_value1
  } else if (anchor_range_direction == "high") {
    category_info$overall_max
  } else if (anchor_range_direction == "point") {
    anchor_value1
  } else {
    anchor_value2
  }

  anchor_support_overlap <- is.finite(anchor_effective_lower) &&
    is.finite(anchor_effective_upper) &&
    anchor_effective_upper >= category_info$overall_min - 1e-8 &&
    anchor_effective_lower <= category_info$overall_max + 1e-8
  if (!anchor_support_overlap) {
    stop(
      "STOP: Frozen anchor interval does not overlap observed training support for ",
      food_var, ". Anchor=[", fmt_num(anchor_effective_lower), ", ",
      fmt_num(anchor_effective_upper), "]; observed=[",
      fmt_num(category_info$overall_min), ", ",
      fmt_num(category_info$overall_max), "]."
    )
  }

  anchor_category_overlap <- anchor_effective_upper >= category_low - 1e-8 &&
    anchor_effective_lower <= category_high + 1e-8

  anchor_plateau_overlap <- is.finite(plateau_05["lower"]) &&
    is.finite(plateau_05["upper"]) &&
    anchor_effective_upper >= unname(plateau_05["lower"]) - 1e-8 &&
    anchor_effective_lower <= unname(plateau_05["upper"]) + 1e-8

  anchor_prediction_mask <- fit_rcs3_pre[[food_var]] >=
    anchor_effective_lower - 1e-8 &
    fit_rcs3_pre[[food_var]] <= anchor_effective_upper + 1e-8

  anchor_significant_protective_anywhere <- any(
    fit_rcs3_pre$upper < 1,
    na.rm = TRUE
  )
  anchor_significant_protective_in_range <- any(
    fit_rcs3_pre$upper[anchor_prediction_mask] < 1,
    na.rm = TRUE
  )

  anchor_min_predicted_hr <- if (any(anchor_prediction_mask)) {
    min(fit_rcs3_pre$yhat[anchor_prediction_mask], na.rm = TRUE)
  } else {
    NA_real_
  }
  anchor_max_predicted_hr <- if (any(anchor_prediction_mask)) {
    max(fit_rcs3_pre$yhat[anchor_prediction_mask], na.rm = TRUE)
  } else {
    NA_real_
  }

  anchor_validation_status <- case_when(
    !isTRUE(shape_info$range_eligible) ~
      "Current RCS association pattern not statistically confirmatory",
    anchor_category_overlap && anchor_plateau_overlap ~
      "Anchor overlaps both current Cox lowest-risk category and RCS 5% plateau",
    anchor_category_overlap && !anchor_plateau_overlap ~
      "Anchor overlaps Cox category; RCS 5% plateau does not overlap",
    !anchor_category_overlap && anchor_plateau_overlap ~
      "Anchor overlaps RCS 5% plateau; Cox category does not overlap",
    TRUE ~
      "Anchor does not overlap current Cox category or RCS 5% plateau"
  )

  anchor_protective_status <- if (anchor_significant_protective_in_range) {
    "At least one point has nominal pointwise 95% CI upper bound < 1 within the anchor"
  } else {
    "No point has nominal pointwise 95% CI upper bound < 1 within the anchor"
  }

  # Final score boundaries are exactly those in the frozen anchor CSV.
  value1 <- anchor_value1
  value2 <- anchor_value2
  range_direction <- anchor_range_direction
  range_source <- as.character(anchor_row$anchor_range_source[1])
  refined_by_RCS_05 <- FALSE

  score_lower <- anchor_effective_lower
  score_upper <- if (range_direction == "high") Inf else anchor_effective_upper
  score_rule <- dplyr::case_when(
    range_direction == "low" ~ paste0("x <= ", fmt_num(value1)),
    range_direction == "high" ~ paste0("x >= ", fmt_num(value1)),
    range_direction == "point" ~ paste0("x = ", fmt_num(value1)),
    range_direction == "middle" ~ paste0(
      fmt_num(value1), " <= x <= ", fmt_num(value2)
    ),
    TRUE ~ "not_defined"
  )

  food_rcs_range <- tibble(
    food = food_var,
    value1 = value1,
    value2 = value2,
    range_direction = range_direction,
    range_source = range_source,
    score_lower = score_lower,
    score_upper = score_upper,
    score_rule = score_rule,
    anchor_value1 = anchor_value1,
    anchor_value2 = anchor_value2,
    anchor_direction = anchor_range_direction,
    anchor_range_source = range_source,
    anchor_effective_lower = anchor_effective_lower,
    anchor_effective_upper = anchor_effective_upper,
    anchor_category_overlap = anchor_category_overlap,
    anchor_plateau_05_overlap = anchor_plateau_overlap,
    anchor_support_overlap = anchor_support_overlap,
    anchor_significant_protective_anywhere =
      anchor_significant_protective_anywhere,
    anchor_significant_protective_in_range =
      anchor_significant_protective_in_range,
    anchor_min_predicted_hr = anchor_min_predicted_hr,
    anchor_max_predicted_hr = anchor_max_predicted_hr,
    anchor_protective_status = anchor_protective_status,
    anchor_validation_status = anchor_validation_status,
    candidate_value1 = candidate_value1,
    candidate_value2 = candidate_value2,
    candidate_range_direction = candidate_range_direction,
    candidate_range_source = candidate_range_source,
    candidate_score_lower = candidate_score_lower,
    candidate_score_upper = candidate_score_upper,
    candidate_score_rule = candidate_score_rule,
    candidate_refined_by_RCS_5pct = candidate_refined_by_RCS_05,
    category_value1 = category_value1,
    category_value2 = category_value2,
    refined_by_RCS_5pct = refined_by_RCS_05,
    category_RCS_direction_consistency = direction_consistency,
    reference_category = ref_label,
    reference_median = ref_value,
    category_range_low = category_low,
    category_range_high = category_high
  )

  # Export the complete prediction grid with exact percentile/knot/reference
  # flags for reproducibility and for downstream histogram overlays if needed.
  prediction_export <- fit_rcs3_pre %>%
    mutate(
      food = food_var,
      display_name = display_name,
      is_reference_value = abs(.data[[food_var]] - ref_value) < 1e-10,
      is_p10 = abs(.data[[food_var]] - p_values["P10"]) < 1e-10,
      is_p50 = abs(.data[[food_var]] - p_values["P50"]) < 1e-10,
      is_p90 = abs(.data[[food_var]] - p_values["P90"]) < 1e-10,
      is_rcs_knot_1 = abs(.data[[food_var]] - knots[1]) < 1e-10,
      is_rcs_knot_2 = abs(.data[[food_var]] - knots[2]) < 1e-10,
      is_rcs_knot_3 = abs(.data[[food_var]] - knots[3]) < 1e-10
    )

  write.csv(
    prediction_export,
    file.path(PRED_DIR, paste0(food_var, "_rcs3.csv")),
    row.names = FALSE
  )

  range_results[[i]] <- food_rcs_range

  model_results[[i]] <- tibble(
    group = FOOD_MAP$group[i],
    food = food_var,
    display_name = display_name,
    category_variable = category_var,
    reference_category_label = ref_label,
    reference_median_label = ref_label_median,
    reference_median_used = ref_value,
    reference_HR = as.numeric(ref_prediction$yhat[1]),
    reference_category_rank = category_info$reference_category_rank,
    number_of_categories = category_info$n_categories,
    reference_category_n = category_info$reference_category_n,
    reference_category_events = category_info$reference_category_events,
    category_range_low = category_low,
    category_range_high = category_high,
    category_range_direction = category_direction,
    overall_training_min = category_info$overall_min,
    overall_training_max = category_info$overall_max,
    model_n = model_n_actual,
    model_events = model_events_actual,
    zero_intake_n = sum(x <= 1e-12),
    zero_intake_percent = 100 * zero_prop,
    p10 = unname(p_values["P10"]),
    p50 = unname(p_values["P50"]),
    p90 = unname(p_values["P90"]),
    rcs_knot_1 = knots[1],
    rcs_knot_2 = knots[2],
    rcs_knot_3 = knots[3],
    rcs_knot_method = knot_info$method,
    positive_p10 = unname(knot_info$positive_quantiles[1]),
    positive_p50 = unname(knot_info$positive_quantiles[2]),
    positive_p90 = unname(knot_info$positive_quantiles[3]),
    prediction_support_low = support_low,
    prediction_support_high = support_high,
    prediction_support_method = support_method,
    p_overall = p_overall,
    p_non_linear = p_non_linear,
    curve_shape = shape_info$curve_shape,
    raw_curve_shape = shape_info$raw_curve_shape,
    shape_evidence = shape_info$shape_evidence,
    nonlinear_evidence_supported = shape_info$nonlinear_supported,
    overall_association_supported = shape_info$overall_supported,
    eligible_for_plateau_refinement = eligible_for_plateau_refinement,
    x_min = shape_info$x_min,
    hr_min = shape_info$hr_min,
    left_tail_hr_ratio_to_min = shape_info$left_tail_ratio,
    right_tail_hr_ratio_to_min = shape_info$right_tail_ratio,
    rcs_range_eligible = shape_info$range_eligible,
    plateau_05_low = unname(plateau_05["lower"]),
    plateau_05_high = unname(plateau_05["upper"]),
    plateau_05_threshold_hr = unname(plateau_05["threshold"]),
    plateau_10_low = unname(plateau_10["lower"]),
    plateau_10_high = unname(plateau_10["upper"]),
    plateau_10_threshold_hr = unname(plateau_10["threshold"]),
    plateau_20_low = unname(plateau_20["lower"]),
    plateau_20_high = unname(plateau_20["upper"]),
    plateau_20_threshold_hr = unname(plateau_20["threshold"]),
    plateau_05_overlaps_reference_category = plateau_overlap,
    plateau_05_overlap_low = plateau_overlap_low,
    plateau_05_overlap_high = plateau_overlap_high,
    category_RCS_direction_consistency = direction_consistency,
    primary_range_source = range_source,
    primary_range_direction = range_direction,
    primary_range_value1 = value1,
    primary_range_value2 = value2,
    anchor_effective_lower = anchor_effective_lower,
    anchor_effective_upper = anchor_effective_upper,
    anchor_category_overlap = anchor_category_overlap,
    anchor_plateau_05_overlap = anchor_plateau_overlap,
    anchor_significant_protective_anywhere =
      anchor_significant_protective_anywhere,
    anchor_significant_protective_in_range =
      anchor_significant_protective_in_range,
    anchor_protective_status = anchor_protective_status,
    anchor_validation_status = anchor_validation_status,
    candidate_range_source = candidate_range_source,
    candidate_range_direction = candidate_range_direction,
    candidate_range_value1 = candidate_value1,
    candidate_range_value2 = candidate_value2,
    candidate_refined_by_RCS_5pct = candidate_refined_by_RCS_05,
    candidate_score_rule = candidate_score_rule,
    category_range_value1 = category_value1,
    category_range_value2 = category_value2,
    refined_by_RCS_5pct = refined_by_RCS_05,
    score_lower = score_lower,
    score_upper = score_upper,
    score_rule = score_rule
  )

  reference_audit_list[[i]] <- tibble(
    group = FOOD_MAP$group[i],
    food = food_var,
    category_variable = category_var,
    reference_category_label = ref_label,
    reference_median_label = ref_label_median,
    reference_median_used = ref_value,
    reference_category_min = category_low,
    reference_category_max = category_high,
    reference_category_n = category_info$reference_category_n,
    reference_category_events = category_info$reference_category_events,
    reference_category_rank = category_info$reference_category_rank,
    number_of_categories = category_info$n_categories,
    category_range_direction = category_direction,
    range_value1 = value1,
    range_value2 = value2
  )

  zero_knot_audit_list[[i]] <- tibble(
    food = food_var,
    n_complete_case = length(x),
    zero_intake_n = sum(x <= 1e-12),
    zero_intake_percent = 100 * zero_prop,
    positive_intake_n = sum(x > 1e-12),
    full_sample_p10 = unname(p_values["P10"]),
    full_sample_p50 = unname(p_values["P50"]),
    full_sample_p90 = unname(p_values["P90"]),
    actual_knot_1 = knots[1],
    actual_knot_2 = knots[2],
    actual_knot_3 = knots[3],
    knot_method = knot_info$method,
    prediction_support_low = support_low,
    prediction_support_high = support_high,
    prediction_support_method = support_method
  )

  # ---------------------------------------------------------------------------
  # Integrated robustness audits: parallel 5%/10%/20% participant bootstrap
  # plus P90-to-P99 tail support. The same B=500 curves feed all three plateau
  # tolerances. No resampling changes the frozen main-fit coefficients or rules.
  # ---------------------------------------------------------------------------
  # Bootstrap the current-model plateau/candidate boundaries, not the frozen
  # scoring anchors. The anchor itself is not re-selected in bootstrap samples.
  interval_bootstrap_mode <- if (
    isTRUE(candidate_refined_by_RCS_05) &&
      candidate_range_direction == "middle"
  ) {
    "bounded_intersection"
  } else if (candidate_range_direction == "high") {
    "high_open_tail_cox_only"
  } else {
    "low_or_point_cox_only"
  }

  boot_audit <- summarise_bootstrap_cutoffs(
    fit = fit_rcs3,
    food_var = food_var,
    pred_support = pred_support,
    reference_value = ref_value,
    category_low = category_low,
    category_high = category_high,
    category_direction = category_direction,
    main_rule = interval_bootstrap_mode,
    model_index = i
  )
  bootstrap_summary_list[[i]] <- boot_audit$summary
  bootstrap_replicates_list[[i]] <- boot_audit$replicates

  tail_audit <- quantile_tail_support(
    fit = fit_rcs3, food_var = food_var,
    x = x, event = as.integer(fit_data$depression_status == 1)
  )
  tail_summary_list[[i]] <- tail_audit$summary
  tail_bins_list[[i]] <- tail_audit$bins

  # ---------------------------------------------------------------------------
  # V4 RCS plots: retain the fitted curve and interval boundaries, while
  # explicitly displaying empirical P10/P50/P90 values above the plotting area.
  # The percentile lines use the exact colours requested by the user.
  # ---------------------------------------------------------------------------

  label_overall_rcs3 <- paste0(
    "P for overall ",
    format_p_value(p_overall)
  )
  label_non_linear_rcs3 <- paste0(
    "P for non-linear ",
    format_p_value(p_non_linear)
  )

  # Data-adaptive display range:
  # include every point in the actual prediction grid, the three empirical
  # percentiles, the frozen reference value, and the final score boundaries.
  # This prevents the citrus/red-orange panels from being stretched to arbitrary
  # x-axis maxima (7 and 9 in FIX2), while retaining all relevant plotted values.
  plot_x_data_max <- max(
    c(
      fit_rcs3_pre[[food_var]],
      p_values,
      ref_value,
      value1,
      value2
    ),
    na.rm = TRUE
  )

  if (!is.finite(plot_x_data_max) || plot_x_data_max <= 0) {
    stop("STOP: Could not establish a positive data-adaptive x-axis maximum for ",
         food_var)
  }

  plot_x_max <- max(
    plot_x_data_max * (1 + RCS_X_AXIS_PADDING_PROP),
    0.01
  )

  if (any(p_values < -1e-8 | p_values > plot_x_max + 1e-8) ||
      any(fit_rcs3_pre[[food_var]] < -1e-8 |
            fit_rcs3_pre[[food_var]] > plot_x_max + 1e-8)) {
    stop(
      "STOP: Percentile or prediction values fall outside the data-adaptive ",
      "x-axis for ", food_var, "."
    )
  }

  annotation_x <- FOOD_MAP$annotation_x[i]
  plot_support_upper <- fit_rcs3_pre$upper[
    fit_rcs3_pre[[food_var]] >= support_low - 1e-10 &
      fit_rcs3_pre[[food_var]] <= support_high + 1e-10
  ]
  plot_text_y <- max(plot_support_upper, na.rm = TRUE)

  plot_y_span <- diff(range(
    c(fit_rcs3_pre$lower, fit_rcs3_pre$upper, 1),
    finite = TRUE
  ))
  if (!is.finite(plot_y_span) || plot_y_span <= 0) {
    plot_y_span <- max(0.10, abs(plot_text_y) * 0.10)
  }
  percentile_label_y <- plot_text_y + 0.18 * plot_y_span

  percentile_plot_data <- tibble(
    percentile = factor(
      names(p_values),
      levels = c("P10", "P50", "P90")
    ),
    percentile_value = as.numeric(p_values)
  )

  # Labels are arranged in a dedicated top row in P10/P50/P90 order. This
  # prevents overlap when two or more empirical percentiles are tied.
  percentile_label_data <- percentile_plot_data %>%
    mutate(
      label = paste0(
        as.character(percentile), " = ",
        vapply(percentile_value, fmt_quantile, character(1))
      ),
      label_x = plot_x_max * c(0.17, 0.50, 0.83),
      label_y = percentile_label_y
    )

  # Put the exact empirical P10/P50/P90 values ON THE X AXIS as well as
  # retaining the coloured dashed vertical lines and top-of-panel labels.
  # Tied percentiles (e.g. a large zero-intake mass) share one tick label.
  x_axis_base_breaks <- seq(
    0, plot_x_max, by = FOOD_MAP$x_break[i]
  )
  x_axis_breaks <- sort(unique(c(x_axis_base_breaks, as.numeric(p_values))))
  x_axis_labels <- vapply(x_axis_breaks, function(b) {
    pct_names <- names(p_values)[
      abs(as.numeric(p_values) - b) <= max(1e-8, abs(plot_x_max) * 1e-8)
    ]
    if (length(pct_names) > 0L) {
      paste0(fmt_quantile(b), "\n", paste(pct_names, collapse = "/"))
    } else {
      fmt_num(b)
    }
  }, character(1))

  # Freeze per-panel values BEFORE constructing ggplot, so combined plot
  # redrawing cannot capture Wine's P-values through loop-scoped variables.
  pvalue_annotation_data <- tibble(
    x = as.numeric(annotation_x), y = as.numeric(plot_text_y),
    label = paste(label_overall_rcs3, label_non_linear_rcs3, sep = "\n"),
    p_overall = p_overall, p_non_linear = p_non_linear,
    food = food_var
  )
  pvalue_annotation_list[[i]] <- pvalue_annotation_data

  p <- ggplot() +
    # The confidence ribbon is deliberately drawn before the fitted curve.
    geom_ribbon(
      data = fit_rcs3_pre,
      aes(x = .data[[food_var]], ymin = lower, ymax = upper),
      alpha = 0.10,
      fill = col_line[1]
    ) +
    geom_vline(
      data = percentile_plot_data,
      aes(xintercept = percentile_value, colour = percentile),
      linetype = "dashed",
      linewidth = 0.65,
      alpha = 0.95,
      show.legend = TRUE
    ) +
    geom_line(
      data = fit_rcs3_pre,
      aes(x = .data[[food_var]], y = yhat),
      linetype = "solid",
      linewidth = 1.0,
      alpha = 0.85,
      colour = col_line[1]
    ) +
    geom_text(
      data = percentile_label_data,
      aes(
        x = label_x,
        y = label_y,
        label = label,
        colour = percentile
      ),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 5.33,
      fontface = "bold",
      vjust = 0
    ) +
    geom_text(
      data = pvalue_annotation_data,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 6.67,
      fontface = "bold"
    ) +
    scale_colour_manual(
      name = "Training-set percentile",
      values = PERCENTILE_COLORS,
      breaks = c("P10", "P50", "P90"),
      drop = FALSE
    ) +
    guides(
      colour = guide_legend(
        direction = "horizontal",
        title.position = "left",
        override.aes = list(
          linetype = "dashed",
          linewidth = 1.1,
          alpha = 1
        )
      )
    ) +
    scale_y_continuous(
      name = "Hazard ratio (95% CI)",
      expand = expansion(mult = c(0.02, 0.32))
    ) +
    scale_x_continuous(
      name = display_name,
      breaks = x_axis_breaks,
      labels = x_axis_labels,
      limits = c(0, plot_x_max),
      guide = ggplot2::guide_axis(n.dodge = 2, check.overlap = FALSE)
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(color = "black", size = 18, face = "bold"),
      axis.text.x = element_text(color = "black", size = 18, face = "bold"),
      axis.title = element_text(color = "black", size = 18, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 13.3, face = "bold"),
      legend.text = element_text(size = 13.3),
      legend.key.width = grid::unit(0.56, "inches"),
      plot.margin = margin(9, 12, 8, 12)
    ) +
    geom_hline(
      yintercept = 1,
      linetype = 2,
      linewidth = 0.5,
      colour = col_line[2],
      alpha = 0.4
    )

  get_curve_y <- function(x0) {
    stats::approx(
      x = fit_rcs3_pre[[food_var]],
      y = fit_rcs3_pre$yhat,
      xout = x0,
      rule = 2,
      ties = "ordered"
    )$y
  }

  if (range_direction == "middle") {
    boundary_x <- c(value1, value2)
    boundary_y <- vapply(boundary_x, get_curve_y, numeric(1))
    boundary_data <- tibble(
      x = boundary_x,
      y = boundary_y,
      label = sprintf("%.2f", boundary_x)
    )

    p <- p +
      geom_vline(
        xintercept = boundary_x,
        linetype = 2,
        size = 0.5,
        colour = col_line[2],
        alpha = 0.6
      ) +
      geom_point(
        data = boundary_data,
        aes(x = x, y = y),
        colour = col_line[3],
        size = 3
      ) +
      geom_text(
        data = boundary_data,
        aes(x = x, y = y, label = label),
        vjust = -1,
        colour = "black",
        size = 8,
        fontface = "bold"
      )
  } else {
    boundary_y <- get_curve_y(value1)
    boundary_data <- tibble(
      x = value1,
      y = boundary_y,
      label = sprintf("%.2f", value1)
    )

    p <- p +
      geom_vline(
        xintercept = value1,
        linetype = 2,
        size = 0.5,
        colour = col_line[2],
        alpha = 0.6
      ) +
      geom_point(
        data = boundary_data,
        aes(x = x, y = y),
        colour = col_line[3],
        size = 3
      ) +
      geom_text(
        data = boundary_data,
        aes(x = x, y = y, label = label),
        vjust = -1,
        colour = "black",
        size = 8,
        fontface = "bold"
      )
  }

  # The 5th geometry is the explicitly frozen P-value annotation text layer.
  # Store its position for exact ggplot-layer vs summary reconciliation.
  attr(p, "frozen_pvalue_layer_index") <- 5L
  plot_list[[food_var]] <- p

  percentile_plot_audit_list[[i]] <- tibble(
    food = food_var,
    display_name = display_name,
    p10 = unname(p_values["P10"]),
    p50 = unname(p_values["P50"]),
    p90 = unname(p_values["P90"]),
    p10_colour = unname(PERCENTILE_COLORS["P10"]),
    p50_colour = unname(PERCENTILE_COLORS["P50"]),
    p90_colour = unname(PERCENTILE_COLORS["P90"]),
    percentile_lines_dashed = TRUE,
    numeric_labels_at_top = TRUE,
    percentile_values_on_x_axis = TRUE,
    x_axis_percentile_labels = paste(x_axis_labels[
      vapply(x_axis_breaks, function(b) {
        any(abs(as.numeric(p_values) - b) <=
              max(1e-8, abs(plot_x_max) * 1e-8))
      }, logical(1))
    ], collapse = " | "),
    legend_included = TRUE,
    prediction_grid_min = min(fit_rcs3_pre[[food_var]], na.rm = TRUE),
    prediction_grid_max = max(fit_rcs3_pre[[food_var]], na.rm = TRUE),
    displayed_x_min = 0,
    displayed_x_max = plot_x_max,
    x_axis_padding_percent = 100 * (
      plot_x_max / plot_x_data_max - 1
    ),
    x_axis_method = paste0(
      "data_adaptive_prediction_grid_percentiles_reference_and_interval_values_with_P10_P50_P90_axis_ticks_",
      sprintf("%d", round(100 * RCS_X_AXIS_PADDING_PROP)),
      "pct_right_padding"
    ),
    all_percentiles_within_displayed_x_axis = all(
      p_values >= -1e-8 & p_values <= plot_x_max + 1e-8
    ),
    all_prediction_values_within_displayed_x_axis = all(
      fit_rcs3_pre[[food_var]] >= -1e-8 &
        fit_rcs3_pre[[food_var]] <= plot_x_max + 1e-8
    )
  )

  # Preserve original individual PDF dimensions and file naming.
  ggsave(
    filename = file.path(PLOT_DIR, paste0(food_var, "_rcs3.pdf")),
    plot = p,
    width = 10,
    height = 6,
    device = "pdf"
  )

  message(
    "[DONE] ", food_var,
    " | N=", model_n_actual,
    " | events=", model_events_actual,
    " | reference=", fmt_num(ref_value),
    " | primary interval=",
    if (range_direction == "middle") {
      paste0(fmt_num(value1), "–", fmt_num(value2))
    } else if (range_direction == "low") {
      paste0("≤ ", fmt_num(value1))
    } else if (range_direction == "high") {
      paste0("≥ ", fmt_num(value1))
    } else {
      fmt_num(value1)
    },
    " | RCS shape=", shape_info$curve_shape
  )
}

# -----------------------------------------------------------------------------
# 5. Export range/model/reference/knot audits and combine the 15 original plots
# -----------------------------------------------------------------------------

food_rcs_range_all <- bind_rows(range_results)
rcs_model_summary <- bind_rows(model_results)
rcs_reference_audit <- bind_rows(reference_audit_list)
rcs_zero_knot_audit <- bind_rows(zero_knot_audit_list)
rcs_percentile_plot_audit <- bind_rows(percentile_plot_audit_list)

if (nrow(rcs_percentile_plot_audit) != 15L ||
    any(!rcs_percentile_plot_audit$all_percentiles_within_displayed_x_axis)) {
  stop("STOP: P10/P50/P90 plot-annotation audit failed.")
}

write.csv(
  rcs_percentile_plot_audit,
  file.path(AUDIT_DIR, "RCS_percentile_plot_annotation_QC.csv"),
  row.names = FALSE,
  na = ""
)

# Exact component reconciliation is performed on normalized character IDs.
# This avoids false failures caused only by vector attributes/names while still
# rejecting missing, extra, duplicated, or unexpected food components.
expected_food_ids <- trimws(as.character(unname(FOOD_MAP$food)))

reconcile_component_order <- function(x, object_name) {
  if (!is.data.frame(x) || !"food" %in% names(x)) {
    stop("STOP: ", object_name, " must be a data.frame with a food column.")
  }

  observed_food_ids <- trimws(as.character(unname(x$food)))

  if (length(observed_food_ids) != length(expected_food_ids) ||
      anyNA(observed_food_ids) ||
      any(observed_food_ids == "") ||
      anyDuplicated(observed_food_ids) > 0L) {
    stop(
      "STOP: ", object_name,
      " has an invalid row count, missing IDs, or duplicated components. ",
      "Observed rows=", nrow(x),
      "; expected rows=", length(expected_food_ids),
      "; observed IDs=",
      paste(observed_food_ids, collapse = ", ")
    )
  }

  missing_ids <- setdiff(expected_food_ids, observed_food_ids)
  extra_ids <- setdiff(observed_food_ids, expected_food_ids)

  if (length(missing_ids) > 0L || length(extra_ids) > 0L) {
    stop(
      "STOP: Component set mismatch in ", object_name,
      ". Missing: ", paste(missing_ids, collapse = ", "),
      ". Extra: ", paste(extra_ids, collapse = ", ")
    )
  }

  x <- x[match(expected_food_ids, observed_food_ids), , drop = FALSE]
  rownames(x) <- NULL
  x
}

food_rcs_range_all <- reconcile_component_order(
  food_rcs_range_all, "food_rcs_range_all"
)
rcs_model_summary <- reconcile_component_order(
  rcs_model_summary, "rcs_model_summary"
)

if (nrow(food_rcs_range_all) != 15L ||
    nrow(rcs_model_summary) != 15L) {
  stop("STOP: Final RCS outputs do not contain exactly 15 unique components.")
}

if (any(!is.finite(rcs_model_summary$reference_median_used)) ||
    any(abs(rcs_model_summary$reference_HR - 1) > 1e-6)) {
  stop("STOP: One or more RCS reference-point QC checks failed.")
}

# Bootstrap and tail audits are part of this self-contained script and are not
# read from a previous run.
bootstrap_boundary_audit <- bind_rows(bootstrap_summary_list)
bootstrap_boundary_replicates <- bind_rows(bootstrap_replicates_list)

# Long-format participant-bootstrap output: one row per component x replicate
# x tolerance (15 x 500 x 3 = 22,500 rows). Lower/upper boundaries and the
# threshold HR are retained together for transparent downstream review.
bootstrap_plateau_values_long <- bootstrap_boundary_replicates %>%
  dplyr::select(
    food, bootstrap_rep, x_min, u_geometry, main_rule,
    tidyselect::matches("^plateau_(05|10|20)_(lower|upper|threshold)$")
  ) %>%
  tidyr::pivot_longer(
    cols = tidyselect::matches("^plateau_(05|10|20)_(lower|upper|threshold)$"),
    names_to = c("tolerance_percent", ".value"),
    names_pattern = "plateau_(05|10|20)_(lower|upper|threshold)"
  )
bootstrap_plateau_overlap_long <- bootstrap_boundary_replicates %>%
  dplyr::select(
    food, bootstrap_rep,
    tidyselect::matches("^plateau_(05|10|20)_overlaps_cox_category$")
  ) %>%
  tidyr::pivot_longer(
    cols = tidyselect::matches("^plateau_(05|10|20)_overlaps_cox_category$"),
    names_to = "tolerance_percent",
    names_pattern = "plateau_(05|10|20)_overlaps_cox_category",
    values_to = "plateau_overlaps_cox_category"
  )
bootstrap_plateau_intersection_long <- bootstrap_boundary_replicates %>%
  dplyr::select(
    food, bootstrap_rep,
    tidyselect::matches("^plateau_(05|10|20)_intersection_(lower|upper)$")
  ) %>%
  tidyr::pivot_longer(
    cols = tidyselect::matches("^plateau_(05|10|20)_intersection_(lower|upper)$"),
    names_to = c("tolerance_percent", ".value"),
    names_pattern = "plateau_(05|10|20)_intersection_(lower|upper)"
  ) %>%
  dplyr::rename(intersection_lower = lower, intersection_upper = upper)

bootstrap_boundary_replicates_long <- bootstrap_plateau_values_long %>%
  dplyr::left_join(
    bootstrap_plateau_overlap_long,
    by = c("food", "bootstrap_rep", "tolerance_percent")
  ) %>%
  dplyr::left_join(
    bootstrap_plateau_intersection_long,
    by = c("food", "bootstrap_rep", "tolerance_percent")
  ) %>%
  dplyr::mutate(
    tolerance_percent = as.integer(tolerance_percent),
    tolerance = tolerance_percent / 100,
    plateau_width = upper - lower,
    intersection_width = intersection_upper - intersection_lower
  ) %>%
  dplyr::arrange(food, bootstrap_rep, tolerance_percent)

tail_support_audit <- bind_rows(tail_summary_list)
tail_support_bins <- bind_rows(tail_bins_list)
rcs_pvalue_annotation_audit <- bind_rows(pvalue_annotation_list)

if (nrow(bootstrap_boundary_audit) != 15L ||
    nrow(bootstrap_boundary_replicates_long) !=
      15L * BOOTSTRAP_B * length(BOOTSTRAP_RANGE_TOLERANCES) ||
    nrow(tail_support_audit) != 15L ||
    nrow(rcs_pvalue_annotation_audit) != 15L ||
    anyDuplicated(bootstrap_boundary_audit$food) > 0L ||
    anyDuplicated(tail_support_audit$food) > 0L) {
  stop("STOP: In-run robustness audits do not cover exactly 15 components.")
}

write.csv(
  bootstrap_boundary_audit,
  file.path(AUDIT_DIR, "RCS_bootstrap_boundary_stability_summary.csv"),
  row.names = FALSE
)
write.csv(
  bootstrap_boundary_replicates,
  file.path(AUDIT_DIR, "RCS_bootstrap_boundary_stability_replicates.csv"),
  row.names = FALSE
)
write.csv(
  bootstrap_boundary_replicates_long,
  file.path(AUDIT_DIR, "RCS_bootstrap_boundary_stability_replicates_long.csv"),
  row.names = FALSE
)
write.csv(
  tail_support_audit,
  file.path(AUDIT_DIR, "RCS_P90_P99_tail_support_summary.csv"),
  row.names = FALSE
)
write.csv(
  tail_support_bins,
  file.path(AUDIT_DIR, "RCS_P90_P99_tail_bins.csv"),
  row.names = FALSE
)
write.csv(
  rcs_pvalue_annotation_audit,
  file.path(AUDIT_DIR, "RCS_panel_Pvalue_text_audit.csv"),
  row.names = FALSE
)

# Publication-ready long sensitivity table: one row per food component and
# tolerance, with full-model boundary estimates and percentile-bootstrap CIs.
bootstrap_boundary_sensitivity_summary <- dplyr::bind_rows(
  lapply(c("05", "10", "20"), function(tag) {
    model_part <- rcs_model_summary %>%
      dplyr::transmute(
        food,
        display_name = FOOD_MAP$display_name[match(food, FOOD_MAP$food)],
        tolerance_percent = as.integer(tag),
        tolerance = as.integer(tag) / 100,
        model_plateau_lower = .data[[paste0("plateau_", tag, "_low")]],
        model_plateau_upper = .data[[paste0("plateau_", tag, "_high")]],
        model_plateau_threshold_hr = .data[[paste0("plateau_", tag, "_threshold_hr")]]
      )
    boot_part <- bootstrap_boundary_audit %>%
      dplyr::transmute(
        food,
        bootstrap_requested,
        bootstrap_successful,
        lower_ci_lower = .data[[paste0("bootstrap_plateau_", tag, "_lower_ci_lower")]],
        lower_ci_upper = .data[[paste0("bootstrap_plateau_", tag, "_lower_ci_upper")]],
        upper_ci_lower = .data[[paste0("bootstrap_plateau_", tag, "_upper_ci_lower")]],
        upper_ci_upper = .data[[paste0("bootstrap_plateau_", tag, "_upper_ci_upper")]],
        lower_estimable_fraction = .data[[paste0("bootstrap_plateau_", tag, "_lower_estimable_fraction")]],
        upper_estimable_fraction = .data[[paste0("bootstrap_plateau_", tag, "_upper_estimable_fraction")]],
        plateau_width_median = .data[[paste0("bootstrap_plateau_", tag, "_width_median")]],
        plateau_width_ci_lower = .data[[paste0("bootstrap_plateau_", tag, "_width_ci_lower")]],
        plateau_width_ci_upper = .data[[paste0("bootstrap_plateau_", tag, "_width_ci_upper")]],
        cox_category_overlap_fraction = .data[[paste0("bootstrap_plateau_", tag, "_cox_overlap_fraction")]],
        intersection_lower_ci_lower = .data[[paste0("bootstrap_plateau_", tag, "_intersection_lower_ci_lower")]],
        intersection_lower_ci_upper = .data[[paste0("bootstrap_plateau_", tag, "_intersection_lower_ci_upper")]],
        intersection_upper_ci_lower = .data[[paste0("bootstrap_plateau_", tag, "_intersection_upper_ci_lower")]],
        intersection_upper_ci_upper = .data[[paste0("bootstrap_plateau_", tag, "_intersection_upper_ci_upper")]],
        intersection_estimable_fraction = .data[[paste0("bootstrap_plateau_", tag, "_intersection_estimable_fraction")]],
        bootstrap_boundary_ci_method
      )
    dplyr::left_join(model_part, boot_part, by = "food")
  })
) %>%
  dplyr::arrange(match(food, FOOD_MAP$food), tolerance_percent)

if (nrow(bootstrap_boundary_sensitivity_summary) != 45L ||
    anyDuplicated(paste(bootstrap_boundary_sensitivity_summary$food,
                        bootstrap_boundary_sensitivity_summary$tolerance_percent)) > 0L) {
  stop("STOP: Parallel 5%/10%/20% bootstrap summary must contain 45 unique rows.")
}

write.csv(
  bootstrap_boundary_sensitivity_summary,
  file.path(AUDIT_DIR, "RCS_bootstrap_boundary_stability_sensitivity_long.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
write.csv(
  bootstrap_boundary_sensitivity_summary,
  file.path(OUT_DIR, "Table_Sx_bootstrap_stability_sensitivity.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# Figure Sx: original fitted plateau boundaries (points) and percentile
# bootstrap 95% CIs (horizontal segments), faceted by food group. These are
# boundary-estimation intervals, not HR confidence bands or noninferiority CIs.
boundary_plot_data <- dplyr::bind_rows(lapply(c("05", "10", "20"), function(tag) {
  dplyr::bind_rows(lapply(c("lower", "upper"), function(side) {
    estimate_col <- if (side == "lower") {
      paste0("plateau_", tag, "_low")
    } else {
      paste0("plateau_", tag, "_high")
    }
    tibble::tibble(
      food = rcs_model_summary$food,
      display_name = FOOD_MAP$display_name[match(rcs_model_summary$food, FOOD_MAP$food)],
      tolerance_percent = as.integer(tag),
      boundary = side,
      estimate = rcs_model_summary[[estimate_col]],
      ci_lower = bootstrap_boundary_audit[[paste0(
        "bootstrap_plateau_", tag, "_", side, "_ci_lower"
      )]],
      ci_upper = bootstrap_boundary_audit[[paste0(
        "bootstrap_plateau_", tag, "_", side, "_ci_upper"
      )]]
    )
  }))
})) %>%
  dplyr::mutate(
    tolerance_label = factor(
      paste0(tolerance_percent, "%"), levels = c("5%", "10%", "20%")
    ),
    plot_label = factor(
      paste0(tolerance_percent, "% ", boundary),
      levels = rev(c("5% lower", "5% upper", "10% lower", "10% upper",
                     "20% lower", "20% upper"))
    )
  )

p_bootstrap_boundary_stability <- ggplot2::ggplot(
  boundary_plot_data,
  ggplot2::aes(x = estimate, y = plot_label, colour = tolerance_label)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = ci_lower, xend = ci_upper,
                 y = plot_label, yend = plot_label),
    linewidth = 0.55, na.rm = TRUE
  ) +
  ggplot2::geom_point(size = 1.35, na.rm = TRUE) +
  ggplot2::facet_wrap(~ display_name, scales = "free_x", ncol = 3) +
  ggplot2::scale_colour_manual(
    values = c("5%" = "#4477AA", "10%" = "#EE7733", "20%" = "#228833"),
    name = "Relative tolerance"
  ) +
  ggplot2::labs(
    title = "Bootstrap stability of RCS near-minimum-risk plateau boundaries",
    subtitle = paste0(
      "Points: original-model boundaries; segments: ",
      round(100 * BOOTSTRAP_BOUNDARY_CI_LEVEL),
      "% percentile bootstrap CIs; B = ", BOOTSTRAP_B,
      " per component"
    ),
    x = "Dietary intake (servings/day)", y = NULL,
    caption = paste(
      "The 5%, 10% and 20% thresholds are applied to the same bootstrap curves.",
      "Intervals quantify boundary uncertainty conditional on the original RCS knots, reference and prediction grid.",
      "They are not HR confidence bands or formal noninferiority limits."
    )
  ) +
  ggplot2::theme_classic(base_size = 8.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 8, hjust = 0.5),
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    axis.text.x = ggplot2::element_text(size = 7),
    axis.text.y = ggplot2::element_text(size = 7),
    axis.title.x = ggplot2::element_text(face = "bold", size = 9),
    legend.position = "bottom",
    plot.caption = ggplot2::element_text(size = 6.8, hjust = 0.5),
    panel.spacing = grid::unit(0.8, "lines")
  )

ggplot2::ggsave(
  file.path(PLOT_DIR, "Figure_Sx_bootstrap_boundary_stability_5_10_20.pdf"),
  p_bootstrap_boundary_stability,
  width = 12, height = 13, units = "in", device = "pdf", limitsize = FALSE
)
ggplot2::ggsave(
  file.path(PLOT_DIR, "Figure_Sx_bootstrap_boundary_stability_5_10_20.png"),
  p_bootstrap_boundary_stability,
  width = 12, height = 13, units = "in", dpi = 600, limitsize = FALSE
)

# Verify the label stored in EACH ACTUAL ggplot panel layer, not merely a
# separate audit object. A mismatch must stop export of the combined PDF.
for (i in seq_len(nrow(FOOD_MAP))) {
  food_var <- FOOD_MAP$food[i]
  p <- plot_list[[food_var]]
  k <- attr(p, "frozen_pvalue_layer_index")
  if (length(k) != 1L || k < 1L || k > length(p$layers)) {
    stop("STOP: Missing P-value geom layer index in panel: ", food_var)
  }
  label_data <- p$layers[[k]]$data
  expected <- rcs_pvalue_annotation_audit[
    rcs_pvalue_annotation_audit$food == food_var, , drop = FALSE
  ]
  model_row <- rcs_model_summary[
    rcs_model_summary$food == food_var, , drop = FALSE
  ]
  if (nrow(label_data) != 1L || nrow(expected) != 1L ||
      nrow(model_row) != 1L ||
      !identical(as.character(label_data$label[1]),
                 as.character(expected$label[1])) ||
      abs(as.numeric(label_data$p_overall[1]) - model_row$p_overall[1]) > 1e-12 ||
      abs(as.numeric(label_data$p_non_linear[1]) - model_row$p_non_linear[1]) > 1e-12 ||
      !identical(as.character(label_data$food[1]), food_var)) {
    stop("STOP: Combined-panel P-value annotation differs from model: ", food_var)
  }
}

write.csv(
  food_rcs_range_all,
  file.path(OUT_DIR, "food_rcs_range_all_final.csv"),
  row.names = FALSE
)
write.csv(
  rcs_model_summary,
  file.path(AUDIT_DIR, "RCS_model_summary_and_range_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  rcs_reference_audit,
  file.path(AUDIT_DIR, "RCS_lowest_risk_reference_category_audit.csv"),
  row.names = FALSE
)
write.csv(
  rcs_zero_knot_audit,
  file.path(AUDIT_DIR, "RCS_zero_intake_P10_P50_P90_knot_audit.csv"),
  row.names = FALSE
)

# Preserve the original 5-column x 3-row combined figure layout and size.
# Individual plots retain their percentile legends; the combined figure uses
# one shared legend to avoid repeating the same key 15 times.
rcs_legend <- cowplot::get_legend(
  plot_list[[FOOD_MAP$food[1]]] +
    theme(legend.position = "bottom")
)
if (is.null(rcs_legend)) {
  stop("STOP: Could not extract the shared P10/P50/P90 legend.")
}

rcs_plots <- lapply(
  unname(plot_list[FOOD_MAP$food]),
  function(z) z + theme(
    legend.position = "none",
    plot.margin = margin(t = 18, r = 8, b = 8, l = 8)
  )
)

rcs_combined_body <- cowplot::plot_grid(
  plotlist = rcs_plots,
  ncol = 5,
  nrow = 3,
  labels = LETTERS[1:15],
  label_size = 36,
  label_fontface = "bold",
  label_y = 1.03,
  label_x = 0.02,
  hjust = 0,
  vjust = 1
)

# Explicitly define the significance marker used in the P-value annotations.
# A single asterisk denotes nominal P < 0.05 based on the unrounded value.
rcs_significance_note <- cowplot::ggdraw() +
  cowplot::draw_label(
    "* P < 0.05 (based on unrounded P values); P values are displayed to three decimal places.",
    x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
    size = 12, colour = "#334155"
  )

rcs_combined <- cowplot::plot_grid(
  rcs_combined_body,
  rcs_legend,
  rcs_significance_note,
  ncol = 1,
  rel_heights = c(1, 0.045, 0.025)
)

ggsave(
  filename = file.path(PLOT_DIR, "rcs3_all_food_combined.pdf"),
  plot = rcs_combined,
  width = 39,
  height = 15.8,
  device = "pdf"
)

# -----------------------------------------------------------------------------
# 6. Original interval-range visualization, now supporting low/middle/high
# -----------------------------------------------------------------------------

food_rcs_range_display <- food_rcs_range_all %>%
  left_join(
    FOOD_MAP %>% select(food, display_name),
    by = "food"
  ) %>%
  mutate(
    min_value = 0,
    max_value = vapply(
      food,
      function(f) {
        max(
          ukb_data_depression_rcs[[f]],
          na.rm = TRUE
        )
      },
      numeric(1)
    ),
    display_name = factor(display_name, levels = FOOD_MAP$display_name)
  )

intervals <- food_rcs_range_display %>%
  rowwise() %>%
  mutate(
    y_start = case_when(
      range_direction == "middle" ~ value1,
      range_direction == "low" ~ min_value,
      range_direction == "high" ~ value1,
      range_direction == "point" ~ value1,
      TRUE ~ NA_real_
    ),
    y_end = case_when(
      range_direction == "middle" ~ value2,
      range_direction == "low" ~ value1,
      range_direction == "high" ~ max_value,
      range_direction == "point" ~ value2,
      TRUE ~ NA_real_
    ),
    label = case_when(
      range_direction == "middle" ~ paste0(fmt_num(value1), "–", fmt_num(value2)),
      range_direction == "low" ~ paste0(fmt_num(min_value), "–", fmt_num(value1)),
      range_direction == "high" ~ paste0("≥ ", fmt_num(value1)),
      range_direction == "point" ~ fmt_num(value1),
      TRUE ~ NA_character_
    ),
    y_mid = (y_start + y_end) / 2
  ) %>%
  ungroup() %>%
  filter(is.finite(y_start), is.finite(y_end))

# Original interval plot styling is retained.
p_food_group <- ggplot(food_rcs_range_display) +
  geom_segment(
    aes(
      x = forcats::fct_rev(factor(display_name)),
      xend = display_name,
      y = min_value,
      yend = max_value
    ),
    size = 2.5,
    color = col[1]
  ) +
  geom_segment(
    data = filter(food_rcs_range_display, range_direction == "middle"),
    aes(x = display_name, xend = display_name, y = value1, yend = value2),
    size = 2.5,
    color = col[3]
  ) +
  geom_segment(
    data = filter(food_rcs_range_display, range_direction == "low"),
    aes(x = display_name, xend = display_name, y = min_value, yend = value1),
    size = 2.5,
    color = col[3]
  ) +
  geom_segment(
    data = filter(food_rcs_range_display, range_direction == "high"),
    aes(x = display_name, xend = display_name, y = value1, yend = max_value),
    size = 2.5,
    color = col[3]
  ) +
  geom_segment(
    data = filter(food_rcs_range_display, range_direction == "point"),
    aes(x = display_name, xend = display_name, y = value1, yend = value2),
    size = 2.5,
    color = col[3]
  ) +
  geom_point(
    data = filter(food_rcs_range_display, range_direction == "point"),
    aes(x = display_name, y = value1),
    size = 2.5,
    color = col[3]
  ) +
  geom_text(
    data = intervals,
    aes(x = display_name, y = y_mid, label = label),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  coord_flip() +
  labs(y = "Consumption level (serving/day)", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "bold", colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black"),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11, face = "bold")
  )

ggsave(
  filename = file.path(PLOT_DIR, "p_food_group.pdf"),
  plot = p_food_group,
  width = 8,
  height = 4.2,
  device = "pdf"
)

# Original tile-plot style retained; range direction is now explicit.
step <- 0.25

tile_list <- lapply(seq_len(nrow(food_rcs_range_display)), function(i) {
  z <- food_rcs_range_display[i, , drop = FALSE]
  x_low <- z$min_value[1]
  x_high <- z$max_value[1]

  starts <- seq(x_low, x_high, by = step)
  if (length(starts) == 0L) starts <- x_low
  ends <- pmin(starts + step, x_high)
  mids <- (starts + ends) / 2

  v1 <- z$value1[1]
  v2 <- z$value2[1]
  direction <- z$range_direction[1]

  point <- if (direction == "middle") {
    as.numeric(mids >= v1 & mids <= v2)
  } else if (direction == "low") {
    as.numeric(mids <= v1)
  } else if (direction == "high") {
    as.numeric(mids >= v1)
  } else if (direction == "point") {
    # A single-value interval must not be expanded to an entire 0.25-serving
    # tile. It is drawn separately as an exact point marker below.
    rep(NA_real_, length(mids))
  } else {
    rep(NA_real_, length(mids))
  }

  tibble(
    food = as.character(z$display_name[1]),
    x_start = starts,
    x_end = ends,
    x_mid = mids,
    point = point
  )
})

food_rcs_range_all_tiles <- bind_rows(tile_list) %>%
  mutate(
    food = factor(food, levels = rev(FOOD_MAP$display_name))
  )

p_food_group_tiles <- ggplot(
  food_rcs_range_all_tiles,
  aes(x = x_mid, y = food, fill = factor(point))
) +
  geom_tile(
    data = filter(food_rcs_range_all_tiles, !is.na(point)),
    aes(width = pmax(x_end - x_start, 0.001)),
    color = "white",
    height = 0.9
  ) +
  geom_point(
    data = filter(food_rcs_range_display, range_direction == "point"),
    aes(x = value1, y = display_name),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.5,
    stroke = 0.8,
    fill = col[2],
    colour = "black"
  ) +
  scale_fill_manual(
    values = c("0" = col[1], "1" = col[2]),
    name = "Point",
    na.value = "grey85"
  ) +
  scale_x_continuous(name = "Servings/day", expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

ggsave(
  filename = file.path(PLOT_DIR, "p_food_group_tiles.pdf"),
  plot = p_food_group_tiles,
  width = 8,
  height = 4.2,
  device = "pdf"
)


# -----------------------------------------------------------------------------
# 6B. Final interval derivation table, frozen binary scores, and figures
# -----------------------------------------------------------------------------
# This block does NOT refit or alter any RCS model and does NOT modify the
# existing individual/combined RCS plots or the existing interval/tile plots.
#
# It uses the already-frozen primary interval imported from
# food_rcs_range_all.csv. Current Cox-category/RCS candidate intervals, plateau
# overlap and pointwise CI status are reported as diagnostics; they do not
# replace or move the frozen scoring boundaries in this run.
#
# Binary score:
#   1 = intake lies within the frozen primary interval
#   0 = intake lies outside the frozen primary interval
#   NA = intake is missing/non-finite; no imputation is performed.
# -----------------------------------------------------------------------------

message("\n[6B/7] Exporting interval derivation, binary scores, and new figures ...")

# 6B.1 Validate the exact 15-component interval objects before deriving scores.
required_range_columns <- c(
  "food", "value1", "value2", "range_direction", "range_source",
  "score_lower", "score_upper", "score_rule",
  "anchor_value1", "anchor_value2", "anchor_direction",
  "anchor_category_overlap", "anchor_plateau_05_overlap",
  "anchor_significant_protective_anywhere",
  "anchor_significant_protective_in_range",
  "anchor_validation_status", "candidate_value1", "candidate_value2",
  "candidate_range_direction", "candidate_range_source",
  "candidate_score_rule", "candidate_refined_by_RCS_5pct",
  "category_RCS_direction_consistency", "refined_by_RCS_5pct"
)
required_model_columns <- c(
  "food", "group", "display_name", "reference_category_label",
  "reference_median_used", "category_range_low", "category_range_high",
  "category_range_direction", "overall_training_min", "overall_training_max",
  "p10", "p50", "p90", "p_overall", "p_non_linear", "curve_shape",
  "shape_evidence", "x_min", "rcs_range_eligible", "plateau_05_low", "plateau_05_high",
  "plateau_10_low", "plateau_10_high", "plateau_20_low", "plateau_20_high",
  "plateau_05_overlaps_reference_category",
  "category_RCS_direction_consistency", "primary_range_source",
  "primary_range_direction", "primary_range_value1", "primary_range_value2",
  "anchor_effective_lower", "anchor_effective_upper",
  "anchor_category_overlap", "anchor_plateau_05_overlap",
  "anchor_significant_protective_anywhere",
  "anchor_significant_protective_in_range",
  "anchor_protective_status", "anchor_validation_status",
  "candidate_range_source", "candidate_range_direction",
  "candidate_range_value1", "candidate_range_value2",
  "candidate_refined_by_RCS_5pct", "candidate_score_rule",
  "score_lower", "score_upper", "score_rule", "model_n", "model_events"
)

missing_range_columns <- setdiff(required_range_columns, names(food_rcs_range_all))
missing_model_columns <- setdiff(required_model_columns, names(rcs_model_summary))

if (length(missing_range_columns) > 0L) {
  stop(
    "STOP: food_rcs_range_all is missing columns required for final scoring: ",
    paste(missing_range_columns, collapse = ", ")
  )
}
if (length(missing_model_columns) > 0L) {
  stop(
    "STOP: rcs_model_summary is missing columns required for the derivation table: ",
    paste(missing_model_columns, collapse = ", ")
  )
}

food_rcs_range_all <- reconcile_component_order(
  food_rcs_range_all, "food_rcs_range_all before final scoring"
)
rcs_model_summary <- reconcile_component_order(
  rcs_model_summary, "rcs_model_summary before final scoring"
)

if (nrow(food_rcs_range_all) != 15L ||
    nrow(rcs_model_summary) != 15L ||
    !identical(as.character(food_rcs_range_all$food), as.character(FOOD_MAP$food)) ||
    !identical(as.character(rcs_model_summary$food), as.character(FOOD_MAP$food))) {
  stop("STOP: Final interval/scoring derivation requires the exact frozen 15 foods.")
}

if (!isTRUE(all.equal(
  as.numeric(food_rcs_range_all$value1),
  as.numeric(range_anchor$anchor_value1),
  tolerance = 1e-12,
  check.attributes = FALSE
)) ||
    !isTRUE(all.equal(
      as.numeric(food_rcs_range_all$value2),
      as.numeric(range_anchor$anchor_value2),
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
  stop("STOP: Final score boundaries no longer exactly match the frozen range-anchor CSV.")
}

range_anchor_comparison <- food_rcs_range_all %>%
  select(
    food, value1, value2, range_direction, score_rule,
    candidate_value1, candidate_value2, candidate_range_direction,
    candidate_range_source, candidate_score_rule,
    anchor_category_overlap, anchor_plateau_05_overlap,
    anchor_significant_protective_anywhere,
    anchor_significant_protective_in_range,
    anchor_protective_status, anchor_validation_status
  ) %>%
  left_join(
    rcs_model_summary %>%
      select(
        food, p10, p50, p90, p_overall, p_non_linear, curve_shape,
        plateau_05_low, plateau_05_high
      ),
    by = "food"
  )

write.csv(
  range_anchor_comparison,
  file.path(AUDIT_DIR, "frozen_anchor_vs_current_Cox_RCS_comparison.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
write.csv(
  food_rcs_range_all %>% select(food, value1, value2),
  file.path(OUT_DIR, "food_rcs_range_all_15comp_downstream.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

if (anyNA(food_rcs_range_all$range_direction) ||
    any(food_rcs_range_all$range_direction == "not_defined") ||
    anyNA(food_rcs_range_all$score_lower) ||
    anyNA(food_rcs_range_all$score_upper) ||
    anyNA(food_rcs_range_all$score_rule) ||
    any(food_rcs_range_all$score_rule == "not_defined")) {
  stop("STOP: At least one food component has an undefined binary scoring interval.")
}

# Human-readable interval formatters. Numeric fields remain available in the
# same table; these labels are intended for manuscript/supplementary review.
format_category_interval <- function(direction, category_low, category_high,
                                      overall_min, overall_max) {
  if (direction == "low") {
    return(paste0("[", fmt_num(overall_min), ", ", fmt_num(category_high), "]"))
  }
  if (direction == "high") {
    return(paste0("[", fmt_num(category_low), ", ", fmt_num(overall_max), "]"))
  }
  if (direction == "middle") {
    return(paste0("[", fmt_num(category_low), ", ", fmt_num(category_high), "]"))
  }
  if (direction == "point") {
    return(paste0("{", fmt_num(category_low), "}"))
  }
  "Not defined"
}

format_final_interval <- function(direction, lower, upper, overall_min, overall_max) {
  if (direction == "low") {
    return(paste0("[", fmt_num(overall_min), ", ", fmt_num(upper), "]"))
  }
  if (direction == "high") {
    return(paste0("[", fmt_num(lower), ", ", fmt_num(overall_max), "] (upper tail)"))
  }
  if (direction == "middle") {
    return(paste0("[", fmt_num(lower), ", ", fmt_num(upper), "]"))
  }
  if (direction == "point") {
    return(paste0("{", fmt_num(lower), "} (single value)"))
  }
  "Not defined"
}

format_rcs_plateau <- function(lower, upper) {
  if (length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper)) {
    return("Not estimable")
  }
  paste0("[", fmt_num(lower), ", ", fmt_num(upper), "]")
}

# 6B.2 Build the one-row-per-food derivation table. It explicitly shows:
#   (i) Cox lowest-risk category and its observed intake range;
#   (ii) the category median used as the RCS reference;
#   (iii) P10/P50/P90 and RCS 5/10/20% plateaus;
#   (iv) overlap/direction/support checks;
#   (v) exact frozen anchor interval, current-model candidate and compatibility;
#   (vi) whether a nominal pointwise protective segment is present;
#   (vii) the final binary scoring rule.
interval_derivation_table <- FOOD_MAP %>%
  select(group, food, display_name, category_var) %>%
  left_join(
    rcs_model_summary %>%
      select(
        food,
        reference_category_label,
        reference_median_used,
        reference_category_rank,
        number_of_categories,
        reference_category_n,
        reference_category_events,
        category_range_low,
        category_range_high,
        category_range_direction,
        overall_training_min,
        overall_training_max,
        model_n,
        model_events,
        p10, p50, p90,
        rcs_knot_1, rcs_knot_2, rcs_knot_3,
        rcs_knot_method,
        p_overall,
        p_non_linear,
        curve_shape,
        shape_evidence,
        x_min,
        rcs_range_eligible,
        plateau_05_low,
        plateau_05_high,
        plateau_05_threshold_hr,
        plateau_10_low,
        plateau_10_high,
        plateau_10_threshold_hr,
        plateau_20_low,
        plateau_20_high,
        plateau_20_threshold_hr,
        plateau_05_overlaps_reference_category
      ),
    by = "food"
  ) %>%
  left_join(
    food_rcs_range_all %>%
      select(
        food,
        cox_category_value1 = category_value1,
        cox_category_value2 = category_value2,
        final_value1 = value1,
        final_value2 = value2,
        final_range_direction = range_direction,
        final_range_source = range_source,
        final_score_lower = score_lower,
        final_score_upper = score_upper,
        final_score_rule = score_rule,
        final_direction_consistency = category_RCS_direction_consistency,
        final_refined_by_RCS_5pct = refined_by_RCS_5pct,
        anchor_effective_lower,
        anchor_effective_upper,
        anchor_category_overlap,
        anchor_plateau_05_overlap,
        anchor_significant_protective_anywhere,
        anchor_significant_protective_in_range,
        anchor_protective_status,
        anchor_validation_status,
        candidate_value1,
        candidate_value2,
        candidate_range_direction,
        candidate_range_source,
        candidate_score_lower,
        candidate_score_upper,
        candidate_score_rule,
        candidate_refined_by_RCS_5pct
      ),
    by = "food",
    suffix = c("", "_range")
  ) %>%
  mutate(
    cox_lowest_risk_category = reference_category_label,
    cox_reference_category_observed_range = paste0(
      "[", vapply(category_range_low, fmt_num, character(1)),
      ", ", vapply(category_range_high, fmt_num, character(1)), "]"
    ),
    cox_candidate_interval = mapply(
      format_category_interval,
      category_range_direction,
      category_range_low,
      category_range_high,
      overall_training_min,
      overall_training_max,
      USE.NAMES = FALSE
    ),
    rcs_reference_value = reference_median_used,
    rcs_5pct_plateau = mapply(
      format_rcs_plateau,
      plateau_05_low,
      plateau_05_high,
      USE.NAMES = FALSE
    ),
    rcs_10pct_plateau = mapply(
      format_rcs_plateau,
      plateau_10_low,
      plateau_10_high,
      USE.NAMES = FALSE
    ),
    rcs_20pct_plateau = mapply(
      format_rcs_plateau,
      plateau_20_low,
      plateau_20_high,
      USE.NAMES = FALSE
    ),
    candidate_interval = mapply(
      format_final_interval,
      candidate_range_direction,
      candidate_score_lower,
      candidate_score_upper,
      overall_training_min,
      overall_training_max,
      USE.NAMES = FALSE
    ),
    interval_derivation = if_else(
      food %in% c("potatoes", "red_meats"),
      "V5.3.1: frozen RCS20 20% plateau intersected with Model 2 lowest-risk category",
      "Retained unchanged from the frozen V5.3 RCS20 interval framework"
    ),
    interval_derivation_reason = case_when(
      food %in% c("potatoes", "red_meats") &
        !anchor_significant_protective_in_range ~
        paste0(
          "The final interval is the intersection of the frozen V5.3 RCS20 20% ",
          "near-minimum-risk plateau and the Model 2 lowest-risk category. ",
          "No point within the final interval met the pointwise 95% CI upper-bound <1 ",
          "criterion; it is not described as statistically significantly protective. ",
          "Current RCS re-fit is reported as validation: ", anchor_validation_status, "."
        ),
      food %in% c("potatoes", "red_meats") ~
        paste0(
          "The final interval is the intersection of the frozen V5.3 RCS20 20% ",
          "near-minimum-risk plateau and the Model 2 lowest-risk category. ",
          "At least one point had a pointwise 95% CI upper bound <1, but the ",
          "interval was not selected by searching for significant P values. ",
          "Current RCS re-fit: ", anchor_validation_status, "."
        ),
      !anchor_significant_protective_anywhere ~
        paste0(
          "The frozen V5.3 interval was retained unchanged. No intake value in the ",
          "current RCS prediction grid had a nominal pointwise 95% CI upper bound <1; ",
          "the interval is an operational scoring rule, not a statistically ",
          "significant protective range. Current validation: ", anchor_validation_status, "."
        ),
      !anchor_significant_protective_in_range ~
        paste0(
          "The frozen V5.3 interval was retained unchanged. No nominally significant ",
          "protective segment was observed within this interval; it is not described ",
          "as a significantly beneficial interval. Current validation: ",
          anchor_validation_status, "."
        ),
      TRUE ~
        paste0(
          "The frozen V5.3 interval was retained unchanged. At least one point within ",
          "the interval had a nominal pointwise 95% CI upper bound below 1; the ",
          "interval itself was not selected by searching for significant P values. ",
          "Current validation: ", anchor_validation_status, "."
        )
    ),
    final_primary_interval = mapply(
      format_final_interval,
      final_range_direction,
      final_score_lower,
      final_score_upper,
      overall_training_min,
      overall_training_max,
      USE.NAMES = FALSE
    ),
    binary_score_1_rule = paste0("1 if intake is within: ", final_score_rule),
    binary_score_0_rule = paste0("0 if intake is outside: ", final_score_rule),
    missing_value_rule = "NA; no imputation"
  ) %>%
  arrange(match(food, FOOD_MAP$food)) %>%
  select(
    group,
    food,
    display_name,
    category_variable = category_var,
    cox_lowest_risk_category,
    cox_reference_category_rank = reference_category_rank,
    cox_reference_category_n = reference_category_n,
    cox_reference_category_events = reference_category_events,
    cox_reference_category_observed_range,
    cox_candidate_interval,
    cox_category_min = category_range_low,
    cox_category_max = category_range_high,
    cox_category_direction = category_range_direction,
    rcs_reference_value,
    p10, p50, p90,
    rcs_knot_1, rcs_knot_2, rcs_knot_3,
    rcs_knot_method,
    p_overall,
    p_non_linear,
    rcs_curve_shape = curve_shape,
    rcs_shape_evidence = shape_evidence,
    rcs_minimum_risk_point = x_min,
    rcs_5pct_plateau,
    rcs_10pct_plateau,
    rcs_20pct_plateau,
    # Keep the model-derived HR threshold required by the integrated Word report.
    plateau_20_threshold_hr,
    rcs_5pct_overlaps_cox_interval = plateau_05_overlaps_reference_category,
    cox_rcs_direction_consistency = final_direction_consistency,
    candidate_interval,
    candidate_interval_source = candidate_range_source,
    candidate_interval_rule = candidate_score_rule,
    candidate_refined_by_RCS_5pct,
    anchor_overlaps_cox_category = anchor_category_overlap,
    anchor_overlaps_rcs_5pct_plateau = anchor_plateau_05_overlap,
    anchor_validation_status,
    nominal_95ci_protective_anywhere = anchor_significant_protective_anywhere,
    nominal_95ci_protective_within_anchor = anchor_significant_protective_in_range,
    anchor_protective_status,
    interval_derivation,
    interval_derivation_reason,
    interval_source_code = final_range_source,
    final_interval_direction = final_range_direction,
    final_interval_lower = final_score_lower,
    final_interval_upper = final_score_upper,
    final_primary_interval,
    binary_score_1_rule,
    binary_score_0_rule,
    missing_value_rule,
    model_n,
    model_events
  )

# Append conditional bootstrap CI and empirical P90-P99 tail support.
interval_derivation_table <- interval_derivation_table %>%
  left_join(bootstrap_boundary_audit, by = "food") %>%
  left_join(tail_support_audit, by = "food") %>%
  mutate(
    plateau_20_rule = "Predicted HR(x) <= 1.20 x minimum predicted HR",
    interval_interpretation = case_when(
      nominal_95ci_protective_within_anchor ~
        "At least one intake value has pointwise 95% CI upper bound < 1; not all values are necessarily significant",
      TRUE ~
        "Scoring interval is not a statistically significant protective interval; interpret as a scoring/near-minimum-risk range only"
    )
  )

if (nrow(interval_derivation_table) != 15L ||
    anyDuplicated(interval_derivation_table$food) > 0L ||
    !identical(as.character(interval_derivation_table$food),
               as.character(FOOD_MAP$food))) {
  stop("STOP: Interval derivation table failed exact 15-component reconciliation.")
}

write.csv(
  interval_derivation_table,
  file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

# ----------------------------------------------------------------------------
# 6B.3  INTEGRATED SUPPLEMENTARY TABLE EXPORT (SINGLE-SCRIPT FINAL)
#
# Source of truth: interval_derivation_table calculated IN MEMORY above.
# The report never reads a prior Step_04_1_1 CSV or sources another R script.
# Three journal-ready sections are created in ONE Word document:
#   A. Frozen component intervals and current Cox/RCS validation
#   B. Conditional bootstrap stability (B resampled model fits)
#   C. High-intake distribution support and model-based P99/P90 contrasts
# The original complete numeric audit CSV is preserved unchanged.
# Word headings use formatted Normal paragraphs (not numbered heading styles).
# Word tables flow automatically and repeat headers; there are NO forced page
# breaks and NO separate report-run dependencies.
# ----------------------------------------------------------------------------

fmt_report <- function(x, digits = 2L) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sprintf(paste0("%.", digits, "f"), as.numeric(x))
}
fmt_report_p <- function(p) {
  if (length(p) != 1L || is.na(p) || !is.finite(p)) return("NA")

  # Keep a consistent three-decimal display in the Word table. The asterisk
  # reflects the unrounded P value and therefore preserves P < 0.05 when the
  # displayed value rounds to 0.050 (e.g. 0.049934 -> 0.050*).
  if (p < 0.001) return("<0.001*")

  paste0(
    sprintf("%.3f", as.numeric(p)),
    if (p < 0.05) "*" else ""
  )
}
fmt_report_ci <- function(lower, upper, digits = 2L) {
  if (length(lower) != 1L || length(upper) != 1L ||
      is.na(lower) || is.na(upper) ||
      !is.finite(lower) || !is.finite(upper)) return("NA")
  paste0(fmt_report(lower, digits), "–", fmt_report(upper, digits))
}
fmt_report_pct <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  paste0(fmt_report(100 * as.numeric(x), digits = 1L), "%")
}
fmt_report_interval <- function(direction, lower, upper) {
  if (is.na(direction)) return("Not defined")
  if (direction == "middle") return(fmt_report_ci(lower, upper))
  if (direction == "high") return(paste0(">= ", fmt_report(lower), " (open upper tail)"))
  if (direction == "low") return(paste0("<= ", fmt_report(upper)))
  if (direction == "point") return(paste0("= ", fmt_report(lower)))
  "Not defined"
}
fmt_report_source <- function(source, direction) {
  # Match the exact interval_source_code values exported by the current
  # V5.3.1/B500 analysis. Do not infer the derivation from interval direction:
  # the source code is the authoritative record of how each final interval
  # was constructed.
  if (identical(
    source,
    "V5.3_RCS_20pct_plateau_intersected_with_anchor"
  )) {
    return(
      "RCS 20% near-minimum-risk plateau intersected with the predefined component-specific scoring anchor"
    )
  }

  if (identical(
    source,
    "user_prespecified_anchor_interval"
  )) {
    return(
      "Predefined component-specific scoring interval retained; current Cox/RCS results used for validation"
    )
  }

  if (identical(
    source,
    "RCS20_20pct_plateau_intersected_with_Model2_lowest_risk_category"
  )) {
    return(
      "RCS 20% near-minimum-risk plateau intersected with the Model 2 lowest-risk category"
    )
  }

  # Fail visibly rather than silently assigning an incorrect Cox-only label
  # if a future analysis introduces a new interval-source code.
  paste0("Unmapped interval source: ", as.character(source))
}

# Check all fields before making a report, so output cannot silently drop an
# uncertainty statistic or use an incomplete stale table from another run.
report_required_cols <- c(
  "food", "display_name", "cox_lowest_risk_category", "cox_candidate_interval",
  "cox_category_direction", "rcs_reference_value", "p10", "p50", "p90",
  "rcs_knot_1", "rcs_knot_2", "rcs_knot_3", "rcs_knot_method",
  "rcs_5pct_plateau", "rcs_10pct_plateau", "rcs_20pct_plateau",
  "plateau_20_threshold_hr",
  "p_overall", "p_non_linear", "rcs_curve_shape", "interval_source_code",
  "anchor_validation_status", "anchor_protective_status",
  "nominal_95ci_protective_anywhere", "nominal_95ci_protective_within_anchor",
  "final_interval_direction", "final_interval_lower", "final_interval_upper",
  "binary_score_1_rule", "model_n", "model_events",
  "bootstrap_requested", "bootstrap_successful",
  "bootstrap_plateau_05_lower_ci_lower", "bootstrap_plateau_05_lower_ci_upper",
  "bootstrap_plateau_05_upper_ci_lower", "bootstrap_plateau_05_upper_ci_upper",
  "bootstrap_conditional_cutpoint_lower_ci_lower",
  "bootstrap_conditional_cutpoint_lower_ci_upper",
  "bootstrap_conditional_cutpoint_upper_ci_lower",
  "bootstrap_conditional_cutpoint_upper_ci_upper",
  "bootstrap_conditional_cutpoint_estimable_fraction",
  "bootstrap_cox_overlap_fraction", "bootstrap_u_geometry_fraction",
  "bootstrap_coefficient_columns", "bootstrap_contrast_matrix_columns",
  "bootstrap_contrast_reconstruction_max_abs_error",
  "empirical_p90", "empirical_p95", "empirical_p99", "n_p90_to_p99",
  "events_p90_to_p99", "n_p99_plus", "events_p99_plus",
  "hr_p99_vs_p90", "hr_p99_vs_p90_lower95", "hr_p99_vs_p90_upper95",
  "p99_vs_p90_ci_excludes_1", "tail_support_meets_descriptive_thresholds"
)
parallel_bootstrap_required_cols <- unlist(lapply(c("05", "10", "20"), function(tag) {
  c(
    paste0("bootstrap_plateau_", tag, "_lower_ci_lower"),
    paste0("bootstrap_plateau_", tag, "_lower_ci_upper"),
    paste0("bootstrap_plateau_", tag, "_upper_ci_lower"),
    paste0("bootstrap_plateau_", tag, "_upper_ci_upper"),
    paste0("bootstrap_plateau_", tag, "_lower_estimable_fraction"),
    paste0("bootstrap_plateau_", tag, "_upper_estimable_fraction"),
    paste0("bootstrap_plateau_", tag, "_cox_overlap_fraction"),
    paste0("bootstrap_plateau_", tag, "_intersection_lower_ci_lower"),
    paste0("bootstrap_plateau_", tag, "_intersection_lower_ci_upper"),
    paste0("bootstrap_plateau_", tag, "_intersection_upper_ci_lower"),
    paste0("bootstrap_plateau_", tag, "_intersection_upper_ci_upper"),
    paste0("bootstrap_plateau_", tag, "_intersection_estimable_fraction")
  )
}))
report_required_cols <- unique(c(report_required_cols, parallel_bootstrap_required_cols))
missing_report_cols <- setdiff(report_required_cols, names(interval_derivation_table))
if (length(missing_report_cols) > 0L) {
  stop("STOP: Integrated report is missing model-derived columns: ",
       paste(missing_report_cols, collapse = ", "))
}
report_data <- interval_derivation_table[
  match(FOOD_MAP$food, interval_derivation_table$food), , drop = FALSE
]
if (nrow(report_data) != 15L ||
    !identical(as.character(report_data$food), as.character(FOOD_MAP$food)) ||
    anyDuplicated(report_data$food) > 0L) {
  stop("STOP: Report component mapping differs from frozen 15-component order.")
}

# Section A: exactly one row per component, with non-rounded P values retained
# in the audit; the publication display uses 4 significant digits for P >= .001.
word_table_data <- tibble::tibble(
  Component = as.character(report_data$display_name),
  `Cox Ref. category` = as.character(report_data$cox_lowest_risk_category),
  `Cox candidate interval` = as.character(report_data$cox_candidate_interval),
  `RCS reference` = vapply(report_data$rcs_reference_value,
                           fmt_report, character(1)),
  `P10 / P50 / P90` = paste(
    vapply(report_data$p10, function(z) fmt_report(z, digits = 3L), character(1)),
    vapply(report_data$p50, function(z) fmt_report(z, digits = 3L), character(1)),
    vapply(report_data$p90, function(z) fmt_report(z, digits = 3L), character(1)), sep = " / "
  ) %>% paste0(
    "\nRCS knots = ",
    paste(
      vapply(report_data$rcs_knot_1, fmt_report, character(1)),
      vapply(report_data$rcs_knot_2, fmt_report, character(1)),
      vapply(report_data$rcs_knot_3, fmt_report, character(1)),
      sep = " / "
    ),
    "\nKnot method = ", report_data$rcs_knot_method
  ),
  `RCS 5% plateau` = as.character(report_data$rcs_5pct_plateau),
  `RCS 10% plateau` = as.character(report_data$rcs_10pct_plateau),
  `RCS 20% near-minimum-risk plateau` = paste0(
    as.character(report_data$rcs_20pct_plateau),
    "\nHR threshold = ",
    vapply(report_data$plateau_20_threshold_hr, fmt_report, character(1))
  ),
  `RCS evidence / derivation` = paste0(
    as.character(report_data$rcs_curve_shape),
    "\nP overall = ", vapply(report_data$p_overall, fmt_report_p, character(1)),
    "; P non-linear = ",
    vapply(report_data$p_non_linear, fmt_report_p, character(1)),
    "\n95% CI upper <1 in evaluated grid / within anchor = ",
    ifelse(report_data$nominal_95ci_protective_anywhere, "Yes", "No"),
    " / ",
    ifelse(report_data$nominal_95ci_protective_within_anchor, "Yes", "No"),
    "\n", as.character(report_data$anchor_validation_status),
    "\n",
    mapply(fmt_report_source,
           as.character(report_data$interval_source_code),
           as.character(report_data$final_interval_direction),
           USE.NAMES = FALSE)
  ),
  `Final interval (score = 1)` = mapply(
    fmt_report_interval,
    as.character(report_data$final_interval_direction),
    report_data$final_interval_lower,
    report_data$final_interval_upper,
    USE.NAMES = FALSE
  )
)

# Section B: 'overlap' and 'U-shape geometry retention' measure DIFFERENT
# properties. Bootstrap cutpoint intervals refer to current-model candidate
# boundaries only; the frozen primary anchor is not re-selected in bootstrap.
supported_u <- as.character(report_data$rcs_curve_shape) == "U-shaped" &
  is.finite(report_data$p_non_linear) & report_data$p_non_linear < 0.05

word_bootstrap_data <- tibble::tibble(
  Component = as.character(report_data$display_name),
  `Successful refits` = paste0(report_data$bootstrap_successful,
                               "/", report_data$bootstrap_requested),
  `U-shape geometry retention` = vapply(seq_len(nrow(report_data)), function(j) {
    if (!supported_u[j]) return("Not a U-shape target")
    fmt_report_pct(report_data$bootstrap_u_geometry_fraction[j])
  }, character(1)),
  `Cox/RCS plateau overlap` = vapply(report_data$bootstrap_cox_overlap_fraction,
                                     fmt_report_pct, character(1)),
  `RCS 5% plateau boundary 95% CI (lower; upper)` = paste0(
    "Lower: ", mapply(fmt_report_ci,
                       report_data$bootstrap_plateau_05_lower_ci_lower,
                       report_data$bootstrap_plateau_05_lower_ci_upper,
                       USE.NAMES = FALSE),
    "\nUpper: ", mapply(fmt_report_ci,
                        report_data$bootstrap_plateau_05_upper_ci_lower,
                        report_data$bootstrap_plateau_05_upper_ci_upper,
                        USE.NAMES = FALSE)
  ),
  `Current-model candidate cutpoint 95% CI (lower; upper)` = vapply(seq_len(nrow(report_data)),
    function(j) {
      if (!is.finite(report_data$bootstrap_conditional_cutpoint_estimable_fraction[j]) ||
          report_data$bootstrap_conditional_cutpoint_estimable_fraction[j] <= 0) {
        return("NA (candidate boundary not conditionally estimable)")
      }
      paste0(
        "Lower: ", fmt_report_ci(
          report_data$bootstrap_conditional_cutpoint_lower_ci_lower[j],
          report_data$bootstrap_conditional_cutpoint_lower_ci_upper[j]),
        "\nUpper: ", fmt_report_ci(
          report_data$bootstrap_conditional_cutpoint_upper_ci_lower[j],
          report_data$bootstrap_conditional_cutpoint_upper_ci_upper[j])
      )
    }, character(1)),
  `Cutpoint estimable` = vapply(
    report_data$bootstrap_conditional_cutpoint_estimable_fraction,
    fmt_report_pct, character(1))
)

# Section B2: parallel 5%/10%/20% plateau-boundary sensitivity results.
format_parallel_plateau_ci <- function(j, tag) {
  paste0(
    "Lower: ", fmt_report_ci(
      report_data[[paste0("bootstrap_plateau_", tag, "_lower_ci_lower")]][j],
      report_data[[paste0("bootstrap_plateau_", tag, "_lower_ci_upper")]][j]
    ),
    "\nUpper: ", fmt_report_ci(
      report_data[[paste0("bootstrap_plateau_", tag, "_upper_ci_lower")]][j],
      report_data[[paste0("bootstrap_plateau_", tag, "_upper_ci_upper")]][j]
    )
  )
}
format_parallel_intersection_ci <- function(j) {
  estimable <- report_data$bootstrap_plateau_20_intersection_estimable_fraction[j]
  category_direction <- as.character(report_data$cox_category_direction[j])
  if (!identical(category_direction, "middle")) {
    return("Not applicable (one-sided or point Cox category)")
  }
  if (!is.finite(estimable) || estimable < 0.80) {
    return(paste0(
      "Not estimable: ", fmt_report_pct(estimable),
      " of bootstrap intersections available (<80% required)"
    ))
  }
  paste0(
    "Lower: ", fmt_report_ci(
      report_data$bootstrap_plateau_20_intersection_lower_ci_lower[j],
      report_data$bootstrap_plateau_20_intersection_lower_ci_upper[j]
    ),
    "\nUpper: ", fmt_report_ci(
      report_data$bootstrap_plateau_20_intersection_upper_ci_lower[j],
      report_data$bootstrap_plateau_20_intersection_upper_ci_upper[j]
    ),
    "\nEstimable: ", fmt_report_pct(estimable)
  )
}
word_parallel_bootstrap_data <- tibble::tibble(
  Component = as.character(report_data$display_name),
  `Successful refits` = paste0(
    report_data$bootstrap_successful, "/", report_data$bootstrap_requested
  ),
  `5% plateau boundary 95% CI` = vapply(
    seq_len(nrow(report_data)), format_parallel_plateau_ci,
    character(1), tag = "05"
  ),
  `10% plateau boundary 95% CI` = vapply(
    seq_len(nrow(report_data)), format_parallel_plateau_ci,
    character(1), tag = "10"
  ),
  `20% plateau boundary 95% CI` = vapply(
    seq_len(nrow(report_data)), format_parallel_plateau_ci,
    character(1), tag = "20"
  ),
  `Cox-category overlap (5% / 10% / 20%)` = vapply(
    seq_len(nrow(report_data)), function(j) paste0(
      "5%: ", fmt_report_pct(report_data$bootstrap_plateau_05_cox_overlap_fraction[j]),
      "\n10%: ", fmt_report_pct(report_data$bootstrap_plateau_10_cox_overlap_fraction[j]),
      "\n20%: ", fmt_report_pct(report_data$bootstrap_plateau_20_cox_overlap_fraction[j])
    ), character(1)
  ),
  `20% bounded intersection boundary 95% CI` = vapply(
    seq_len(nrow(report_data)), format_parallel_intersection_ci, character(1)
  )
)

# Section C: empirical support and exploratory covariate-adjusted model contrast.
word_tail_data <- tibble::tibble(
  Component = as.character(report_data$display_name),
  `P90 / P95 / P99` = paste(
    vapply(report_data$empirical_p90, fmt_report, character(1)),
    vapply(report_data$empirical_p95, fmt_report, character(1)),
    vapply(report_data$empirical_p99, fmt_report, character(1)), sep = " / "
  ),
  `P90-P99 n / events` = paste0(
    format(report_data$n_p90_to_p99, big.mark = ",", trim = TRUE), " / ",
    format(report_data$events_p90_to_p99, big.mark = ",", trim = TRUE)
  ),
  `>=P99 n / events` = paste0(
    format(report_data$n_p99_plus, big.mark = ",", trim = TRUE), " / ",
    format(report_data$events_p99_plus, big.mark = ",", trim = TRUE)
  ),
  `HR P99 vs P90 (95% CI)` = paste0(
    sprintf("%.3f", report_data$hr_p99_vs_p90), " (",
    sprintf("%.3f", report_data$hr_p99_vs_p90_lower95), "–",
    sprintf("%.3f", report_data$hr_p99_vs_p90_upper95), ")",
    ifelse(report_data$p99_vs_p90_ci_excludes_1, "*", "")
  ),
  `95% CI excludes 1` = ifelse(report_data$p99_vs_p90_ci_excludes_1,
                               "Yes*", "No"),
  `Descriptive tail support` = ifelse(
    report_data$tail_support_meets_descriptive_thresholds, "Met", "Not met"
  )
)

# Report-specific quality gates evaluated before any Word file is committed.
fish_i <- match("fish_seafood", report_data$food)
fish_p_value <- report_data$p_non_linear[fish_i]
fish_p_display <- fmt_report_p(fish_p_value)
fish_figure_p_display <- format_p_value(fish_p_value)
fish_p_is_significant <- is.finite(fish_p_value) && fish_p_value < 0.05
report_qc <- tibble::tibble(
  check = c(
    "word_four_report_sections_have_exact_15_components",
    "word_parallel_5_10_20_bootstrap_sensitivity_table_complete",
    "all_15_components_have_estimable_20pct_plateau_boundary_CIs",
    "word_all_components_in_frozen_order",
    "word_bootstrap_B_refits_successful_for_15",
    "word_CI_and_geometry_fields_available_for_15",
    "word_non_estimable_candidate_cutpoint_CI_not_claimed",
    "word_tail_support_fields_complete",
    "fish_non_linear_P_display_three_decimals_with_significance_marker",
    "word_coefficient_contrast_dimensions_reconciled"
  ),
  status = c(
    nrow(word_table_data) == 15L &&
      nrow(word_bootstrap_data) == 15L &&
      nrow(word_parallel_bootstrap_data) == 15L &&
      nrow(word_tail_data) == 15L,
    all(report_data$bootstrap_requested == BOOTSTRAP_B) &&
      all(report_data$bootstrap_successful >=
            ceiling(BOOTSTRAP_B * BOOTSTRAP_MIN_SUCCESS_RATE)) &&
      all(vapply(c("05", "10", "20"), function(tag) {
        all(is.finite(report_data[[paste0("bootstrap_plateau_", tag, "_lower_ci_lower")]])) &&
          all(is.finite(report_data[[paste0("bootstrap_plateau_", tag, "_lower_ci_upper")]])) &&
          all(is.finite(report_data[[paste0("bootstrap_plateau_", tag, "_upper_ci_lower")]])) &&
          all(is.finite(report_data[[paste0("bootstrap_plateau_", tag, "_upper_ci_upper")]]))
      }, logical(1))),
    all(report_data$bootstrap_plateau_20_lower_estimable_fraction >= 0.80) &&
      all(report_data$bootstrap_plateau_20_upper_estimable_fraction >= 0.80),
    identical(as.character(report_data$food), as.character(FOOD_MAP$food)),
    all(report_data$bootstrap_requested == BOOTSTRAP_B) &&
      all(report_data$bootstrap_successful == BOOTSTRAP_B),
    all(is.finite(report_data$bootstrap_cox_overlap_fraction)) &&
      all(is.finite(report_data$bootstrap_u_geometry_fraction[supported_u])),
    all(word_bootstrap_data[["Current-model candidate cutpoint 95% CI (lower; upper)"]][
      report_data$bootstrap_conditional_cutpoint_estimable_fraction <= 0] ==
          "NA (candidate boundary not conditionally estimable)"),
    all(is.finite(report_data$hr_p99_vs_p90)) &&
      all(is.finite(report_data$hr_p99_vs_p90_lower95)) &&
      all(is.finite(report_data$hr_p99_vs_p90_upper95)) &&
      all(!is.na(report_data$tail_support_meets_descriptive_thresholds)),
    fish_p_is_significant &&
      identical(fish_p_display, "0.050*") &&
      identical(fish_figure_p_display, "= 0.050*") &&
      grepl(
        paste0("P non-linear = ", fish_p_display),
        word_table_data[["RCS evidence / derivation"]][fish_i],
        fixed = TRUE
      ) &&
      grepl("\\*$", fish_p_display) &&
      grepl("\\*$", fish_figure_p_display),
    all(report_data$bootstrap_coefficient_columns ==
        report_data$bootstrap_contrast_matrix_columns) &&
      all(report_data$bootstrap_contrast_reconstruction_max_abs_error <= 1e-8)
  )
)
write.csv(report_qc,
          file.path(AUDIT_DIR, "RCS_integrated_Word_report_QC.csv"),
          row.names = FALSE)
if (anyNA(report_qc$status) || any(!report_qc$status)) {
  print(report_qc %>% filter(is.na(status) | !status), n = Inf)
  stop("STOP: Integrated Word report QC failed.")
}

# Export the three exactly matching display tables, for editorial audit.
write.csv(word_table_data,
          file.path(OUT_DIR, "Table_Sx_interval_derivation_display.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(word_bootstrap_data,
          file.path(OUT_DIR, "Table_Sx_bootstrap_stability.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(word_parallel_bootstrap_data,
          file.path(OUT_DIR, "Table_Sx_parallel_bootstrap_sensitivity.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(word_tail_data,
          file.path(OUT_DIR, "Table_Sx_high_intake_tail.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

# Table formatting: no heading-style auto-numbering, no unexpected Word
# section breaks, and no attempts to force a 15-row table onto one page.
# A continuous table may flow over pages; Word repeats the header on each.
make_supplementary_ft <- function(dat, widths, font_size = 7) {
  if (ncol(dat) != length(widths) || sum(widths) > 10.85) {
    stop("STOP: Word table width exceeds A4 landscape usable width.")
  }
  ft <- flextable::flextable(dat)
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::fontsize(ft, size = font_size, part = "all")
  ft <- flextable::fontsize(ft, size = font_size + 0.4, part = "header")
  ft <- flextable::bold(ft, bold = TRUE, part = "header")
  ft <- flextable::bg(ft, bg = "#EAF2F8", part = "header")
  ft <- flextable::color(ft, color = "#17324D", part = "header")
  ft <- flextable::align(ft, align = "center", part = "header")
  ft <- flextable::align(ft, align = "center", part = "body")
  ft <- flextable::align(ft, j = 1L, align = "left", part = "body")
  ft <- flextable::valign(ft, valign = "center", part = "all")
  ft <- flextable::padding(ft, padding = 2, part = "all")
  ft <- flextable::width(ft, j = seq_along(widths), width = widths)
  ft <- flextable::set_table_properties(
    ft, layout = "fixed", align = "center",
    # Compatibility: older flextable releases expose only `split` and
    # `keep_with_next` in opts_ft_word(). Header rows are repeated by the
    # Word table header-row behavior; do not pass `repeat_headers` here,
    # because unsupported named arguments make body_add_flextable() fail.
    opts_word = list(split = FALSE, keep_with_next = FALSE)
  )
  ft
}

# Explicit text formatting bypasses any 'heading 2' style numbering installed
# by Word/officer. Only the built-in Normal paragraph style is necessary.
add_supplementary_heading <- function(doc, text, font_size = 10) {
  officer::body_add_fpar(
    doc, officer::fpar(
      officer::ftext(
        text,
        prop = officer::fp_text(
          font.size = font_size, bold = TRUE,
          font.family = "Arial", color = "#17324D"
        )
      ),
      fp_p = officer::fp_par(text.align = "left")
    )
  )
}
add_supplementary_note <- function(doc, text, font_size = 8) {
  officer::body_add_fpar(
    doc, officer::fpar(
      officer::ftext(
        text,
        prop = officer::fp_text(font.size = font_size, font.family = "Arial")
      ),
      fp_p = officer::fp_par(text.align = "left")
    )
  )
}

# Dynamic examples/uncertainty estimates come directly from this run, not old
# hard-coded figures or from a previous report file.
report_row <- function(food) {
  r <- report_data[report_data$food == food, , drop = FALSE]
  if (nrow(r) != 1L) stop("STOP: Report example component missing: ", food)
  r
}
lowfat_row <- report_row("low_fat_dairy")
sweet_row <- report_row("sweetened_bev")
n_u <- sum(supported_u)

word_footnotes <- c(
  paste0(
    "Study population. All analyses used the prespecified, frozen 70% training ",
    "cohort (N=", format(nrow(ukb_data_depression_rcs), big.mark = ","),
    "; depression events=",
    format(sum(ukb_data_depression_rcs$depression_status == 1, na.rm = TRUE),
           big.mark = ","),
    "). The separate 30% test set did not contribute to food-group selection, ",
    "reference-category estimation, RCS fitting or validation. Unless otherwise ",
    "specified, intake is expressed in servings per day."
  ),
  paste0(
    "Abbreviations. Cox, Cox proportional hazards model; RCS, restricted ",
    "cubic spline; HR, hazard ratio; CI, confidence interval; Ref., the ",
    "lowest-estimated-risk category from the rotating-reference categorical ",
    "Cox analysis; P10/P50/P90/P95/P99, empirical exposure percentiles. ",
    "P values are displayed to three decimal places; an asterisk (*) denotes ",
    "nominal P < 0.05 based on the unrounded value. Values below 0.001 are ",
    "displayed as <0.001*."
  ),
  paste0(
    "1. RCS reference and percentile display. For each component, the RCS ",
    "reference was the observed median intake in the lowest-estimated-risk ",
    "categorical Cox group, rather than automatically zero; HR at this ",
    "reference is 1. P10/P50/P90 are the empirical 10th, 50th and 90th ",
    "percentiles of the complete-case training exposure distribution. They ",
    "are displayed as coloured vertical lines, numeric x-axis tick labels and ",
    "numeric labels above each panel. When percentile ties occur because of ",
    "zero-inflation, tied percentiles are identified at the same x-axis value."
  ),
  paste0(
    "2. interval derivation. Thirteen component intervals were retained ",
    "unchanged from the frozen V5.3 RCS20 interval export. For Potatoes and ",
    "Red meats, the final interval is the intersection of the frozen V5.3 RCS20 ",
    "20% near-minimum-risk plateau and the lowest-risk category range from the ",
    "categorical Model 2 analysis (Potatoes 1.00-1.50; Red meats 1.00-1.75 ",
    "servings/day). Bounded middle intervals are inclusive. Low-direction rules ",
    "use intake <= value1; high-direction rules use intake >= value1 with an ",
    "open upper tail; point rules use the specified single value. These are ",
    "operational component-scoring intervals, not clinical intake recommendations."
  ),
  paste0(
    "3. How association evidence is interpreted. The interval is NOT selected ",
    "solely because the pointwise 95% CI lies below HR=1. The RCS curve is ",
    "referenced to the median of the categorical lowest-risk group, so an ",
    "upper 95% CI below 1 indicates a nominally lower hazard than that reference ",
    "at that intake value; it is not, by itself, a definition of the optimal ",
    "range. For every component, the report separately records the overall and ",
    "non-linearity P values, curve shape, 5%/10%/20% near-minimum-risk plateaus, overlap ",
    "with the current Cox category, and whether any pointwise 95% CI upper bound ",
    "is below 1 in the evaluated prediction grid or within the frozen interval. These pointwise CI ",
    "checks are descriptive and are not simultaneous confidence bands."
  ),
  paste0(
    "4. Components without a significantly protective RCS segment. When no ",
    "intake value in the evaluated prediction grid has a pointwise 95% CI upper ",
    "bound below 1, no statistically ",
    "significant protective range is claimed. The pre-existing scoring interval ",
    "is retained only as an operational/historical component rule and is ",
    "explicitly flagged as not confirmed by current RCS evidence. In particular, ",
    "the low-fat dairy anchor is ",
    fmt_report_interval(
      lowfat_row$final_interval_direction,
      lowfat_row$final_interval_lower,
      lowfat_row$final_interval_upper
    ),
    " servings/day, and the sweetened-beverage anchor is ",
    fmt_report_interval(
      sweet_row$final_interval_direction,
      sweet_row$final_interval_lower,
      sweet_row$final_interval_upper
    ),
    " servings/day. Their within-anchor pointwise 95% CI upper<1 status is ",
    ifelse(lowfat_row$nominal_95ci_protective_within_anchor, "present", "not present"),
    " and ",
    ifelse(sweet_row$nominal_95ci_protective_within_anchor, "present", "not present"),
    ", respectively. These thresholds are not described as significantly ",
    "beneficial when the corresponding confidence intervals include 1."
  ),
  paste0(
    "5. RCS pattern and near-minimum plateau. A U/J/reverse-J pattern is ",
    "described as supported nonlinear only when the non-linearity test has ",
    "P<0.05 and the fitted geometry is compatible. A significant overall ",
    "association alone does not establish non-linearity. For each fitted curve, ",
    "the 5%, 10% and 20% plateaus are the contiguous regions around the minimum ",
    "fitted HR where HR(x) is no more than 1.05, 1.10 and 1.20 times the minimum ",
    "fitted HR, respectively. The 20% tolerance is the primary operational ",
    "near-minimum-risk threshold used by the frozen V5.3 RCS20 framework; 5% and ",
    "10% are stricter sensitivity thresholds. These plateaus are point-estimate ",
    "regions; they are not HR confidence limits, clinical targets or proof of ",
    "equivalence. Current-model diagnostics do not move the frozen anchor."
  ),
  paste0(
    "6. Parallel 5%/10%/20% boundary bootstrap. Each component underwent ",
    BOOTSTRAP_B, " participant-level case-bootstrap model refits. The same ",
    "bootstrap coefficient curves were evaluated at all three tolerances, so ",
    "the sensitivity thresholds did not require additional model refits. For ",
    "each tolerance, the reported boundary 95% CIs are the empirical ",
    "2.5th and 97.5th percentiles of the bootstrap distributions of the lower ",
    "and upper plateau boundaries. This quantifies boundary-estimation ",
    "uncertainty conditional on the original RCS knots, Cox reference and ",
    "prediction grid. It is not a confidence band for HR(x), a formal ",
    "noninferiority test, or a bootstrap re-selection of the frozen 15-component ",
    "scoring anchor. The existing 5% candidate-cutpoint audit is retained ",
    "separately for backward compatibility."
  ),
  paste0(
    "7. High-intake tail support. Participant and depression-event counts were ",
    "tabulated for P90-P99 and >=P99. Descriptive support required at least ",
    TAIL_MIN_N_P90_P99, " people and ", TAIL_MIN_EVENTS_P90_P99,
    " events in P90-P99 and at least ", TAIL_MIN_N_P99_PLUS,
    " people and ", TAIL_MIN_EVENTS_P99_PLUS,
    " events at >=P99. The P99-versus-P90 HR is exploratory, without ",
    "multiplicity correction, and cannot establish benefit from indefinite ",
    "high consumption or a clinical upper/lower target."
  ),
  paste0(
    "8. Binary scoring. A participant receives 1 within the frozen anchor ",
    "interval, 0 outside and NA for missing intake, without imputation. The ",
    "0-15 total is calculated only when all 15 components are observed. A ",
    "score of 0 does not prove harm; a score of 1 does not imply causality. ",
    "The updated beta coefficients are applied downstream and do not change ",
    "these frozen intake boundaries."
  )
)

word_title <- paste(
  "Supplementary Table Sx.",
  "Validation of the frozen primary intake intervals and binary component scores"
)
WORD_TABLE_FILE <- file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.docx")

# Output is built in one pass from this run. No old table is read or replaced
# as input, no Word backup copies are created and no versioned duplicate is
# needed. A failed write does not count as a completed final analysis.
doc <- officer::read_docx()
doc <- officer::body_set_default_section(
  doc,
  value = officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(
      top = 0.42, bottom = 0.42, left = 0.42, right = 0.42
    )
  )
)
doc <- add_supplementary_heading(doc, word_title, font_size = 11)
doc <- add_supplementary_note(
  doc, paste0(
    "Fifteen frozen components; current Cox/RCS diagnostics exclusively in the ",
    "70% training cohort; conditional bootstrap and upper-tail audit."
  ), font_size = 8
)

doc <- add_supplementary_heading(
  doc, "A. Frozen component intervals and current Cox/RCS validation", 10
)
doc <- flextable::body_add_flextable(
  doc, value = make_supplementary_ft(
    word_table_data, c(1.05, 0.85, 0.95, 0.65, 1.35, 0.75, 0.75, 1.05, 2.00, 1.10), 7.0
  )
)
doc <- add_supplementary_note(
  doc, "Note. Intake is servings/day; bounded intervals are inclusive; low-direction rules use x <= value1, high-direction rules use x >= value1 (open upper tail), and point rules use the specified single value. ", 7.7
)

doc <- add_supplementary_heading(
  doc, "B1. Conditional bootstrap stability of the 5% current-model candidate boundaries", 10
)
doc <- flextable::body_add_flextable(
  doc, value = make_supplementary_ft(
    word_bootstrap_data, c(1.43, 0.85, 1.38, 1.18, 2.00, 2.06, 1.25), 7.1
  )
)
doc <- add_supplementary_note(
  doc, paste0(
    "Conditional uncertainty only: ", BOOTSTRAP_B,
    " refits per component, frozen knots and categorical reference. ",
    "Cox/RCS overlap is not evidence of U-shape geometry retention. ",
    "The frozen primary anchor is not bootstrapped; NA means the current-model ",
    "candidate boundary was not conditionally estimable."
  ), 7.7
)

doc <- add_supplementary_heading(
  doc, "B2. Parallel 5%/10%/20% bootstrap plateau-boundary sensitivity", 10
)
doc <- flextable::body_add_flextable(
  doc, value = make_supplementary_ft(
    word_parallel_bootstrap_data, c(1.18, 0.82, 1.42, 1.42, 1.42, 1.35, 2.15), 6.8
  )
)
doc <- add_supplementary_note(
  doc, paste0(
    "Note. Each of the ", BOOTSTRAP_B,
    " participant-resampled model curves per component was evaluated at all three tolerances. ",
    "Boundary CIs are percentile intervals for the lower and upper plateau cutpoints, conditional on frozen knots, reference and prediction grid. ",
    "The 20% category-intersection CIs are shown only for bounded two-sided categories; one-sided/point categories are not applicable. ",
    "These are not HR confidence bands or formal noninferiority limits."
  ), 7.2
)

doc <- add_supplementary_heading(
  doc, "C. High-intake tail support and P99-versus-P90 model contrasts", 10
)
doc <- flextable::body_add_flextable(
  doc, value = make_supplementary_ft(
    word_tail_data, c(1.40, 1.50, 1.43, 1.37, 2.34, 1.22, 1.25), 7.1
  )
)
doc <- add_supplementary_note(
  doc, paste0(
    "* Nominal 95% CI excludes 1; exploratory contrast, no multiplicity ",
    "adjustment. Tail support is descriptive and does not establish causality."
  ), 7.7
)

doc <- add_supplementary_heading(doc, "Footnotes and interval-derivation rules", 10)
for (note in word_footnotes) {
  doc <- add_supplementary_note(doc, note, font_size = 8)
}
print(doc, target = WORD_TABLE_FILE)
if (!file.exists(WORD_TABLE_FILE) ||
    !is.finite(file.info(WORD_TABLE_FILE)$size) ||
    file.info(WORD_TABLE_FILE)$size < 10000) {
  stop("STOP: Integrated Word supplementary table was not created correctly.")
}
message("[6B/7] Integrated A/B/C Word supplementary table exported: ",
        WORD_TABLE_FILE)



# A compact scoring-rule file is also exported for direct use in the next
# RWDID construction step. The primary derivation table remains the full audit.
binary_scoring_rules <- food_rcs_range_all %>%
  left_join(
    FOOD_MAP %>% select(food, group, display_name),
    by = "food"
  ) %>%
  left_join(
    rcs_model_summary %>%
      select(food, overall_training_min, overall_training_max,
             reference_category_label, reference_median_used,
             plateau_05_low, plateau_05_high),
    by = "food"
  ) %>%
  mutate(
    interval_source = case_when(
      food %in% c("potatoes", "red_meats") ~
        "FINAL: frozen RCS20 20% plateau intersected with Model 2 lowest-risk category",
      TRUE ~ "V5.3 FINAL: frozen RCS20 interval retained without modification"
    ),
    score_1 = paste0("1 if ", score_rule),
    score_0 = paste0("0 otherwise"),
    missing_rule = "NA; no imputation"
  ) %>%
  arrange(match(food, FOOD_MAP$food)) %>%
  select(
    group, food, display_name,
    reference_category_label, reference_median_used,
    category_range_low, category_range_high,
    plateau_05_low,
    plateau_05_high,
    range_source, interval_source, range_direction,
    value1, value2, score_lower, score_upper,
    score_rule, score_1, score_0, missing_rule,
    anchor_category_overlap, anchor_plateau_05_overlap,
    anchor_significant_protective_anywhere,
    anchor_significant_protective_in_range,
    anchor_protective_status, anchor_validation_status,
    candidate_value1, candidate_value2, candidate_range_direction,
    candidate_range_source, candidate_score_rule,
    overall_training_min, overall_training_max
  )

# 6B.3 Apply the frozen binary rules to the exact Step_03_1_0 training RDS.
# Missing exposure values remain NA. No new split, refitting, imputation, or
# test-set access occurs in this step.
score_one_component <- function(x, direction, value1, value2, food_name) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  score <- rep(NA_integer_, length(x))
  observed <- !is.na(x)

  if (any(!is.finite(x[observed]))) {
    stop("STOP: Non-finite non-missing exposure values found for ", food_name)
  }
  if (any(x[observed] < -1e-8)) {
    stop("STOP: Negative intake values found while scoring ", food_name)
  }

  if (!any(observed)) {
    return(score)
  }

  scale_value <- max(1, abs(value1), ifelse(is.finite(value2), abs(value2), 0))
  tol <- max(1e-8, scale_value * 1e-8)

  if (direction == "middle") {
    if (!is.finite(value1) || !is.finite(value2) || value2 < value1) {
      stop("STOP: Invalid middle-interval boundaries for ", food_name)
    }
    score[observed] <- as.integer(
      x[observed] >= value1 - tol & x[observed] <= value2 + tol
    )
  } else if (direction == "low") {
    if (!is.finite(value1)) {
      stop("STOP: Invalid low-interval upper boundary for ", food_name)
    }
    score[observed] <- as.integer(x[observed] <= value1 + tol)
  } else if (direction == "high") {
    if (!is.finite(value1)) {
      stop("STOP: Invalid high-interval lower boundary for ", food_name)
    }
    score[observed] <- as.integer(x[observed] >= value1 - tol)
  } else if (direction == "point") {
    if (!is.finite(value1)) {
      stop("STOP: Invalid point interval for ", food_name)
    }
    score[observed] <- as.integer(abs(x[observed] - value1) <= tol)
  } else {
    stop("STOP: Unsupported interval direction for ", food_name, ": ", direction)
  }

  score
}

score_data <- data.frame(
  participant_id = trimws(as.character(ukb_data_depression_rcs[[ID_COL]])),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
names(score_data)[1] <- ID_COL

score_columns <- paste0("score_", FOOD_MAP$food)

for (i in seq_len(nrow(FOOD_MAP))) {
  food_var <- FOOD_MAP$food[i]
  range_row <- food_rcs_range_all[
    match(food_var, food_rcs_range_all$food), ,
    drop = FALSE
  ]

  score_data[[score_columns[i]]] <- score_one_component(
    x = ukb_data_depression_rcs[[food_var]],
    direction = as.character(range_row$range_direction[1]),
    value1 = as.numeric(range_row$value1[1]),
    value2 = as.numeric(range_row$value2[1]),
    food_name = food_var
  )
}

score_matrix <- as.matrix(score_data[, score_columns, drop = FALSE])
storage.mode(score_matrix) <- "numeric"
score_data$n_components_scored <- as.integer(rowSums(!is.na(score_matrix)))

complete_15_components <- score_data$n_components_scored == 15L
score_data$RWDID_total_score_complete15 <- NA_integer_

if (any(complete_15_components)) {
  score_data$RWDID_total_score_complete15[complete_15_components] <-
    as.integer(rowSums(
      score_matrix[complete_15_components, , drop = FALSE],
      na.rm = FALSE
    ))
}

score_data$analysis_split <- "training"

if (nrow(score_data) != nrow(ukb_data_depression_rcs) ||
    anyNA(score_data[[ID_COL]]) ||
    any(score_data[[ID_COL]] == "") ||
    anyDuplicated(score_data[[ID_COL]]) > 0L) {
  stop("STOP: Participant-level binary scoring output failed ID/row reconciliation.")
}

write.csv(
  score_data,
  file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

# Component-level 0/1 distribution for audit and manuscript reporting.
binary_score_distribution <- bind_rows(lapply(seq_len(nrow(FOOD_MAP)), function(i) {
  food_var <- FOOD_MAP$food[i]
  score_vec <- score_data[[score_columns[i]]]
  observed_scores <- score_vec[!is.na(score_vec)]
  n_scored <- length(observed_scores)
  n_zero <- sum(observed_scores == 0L)
  n_one <- sum(observed_scores == 1L)

  tibble(
    group = FOOD_MAP$group[i],
    food = food_var,
    display_name = FOOD_MAP$display_name[i],
    n_training = length(score_vec),
    n_missing = sum(is.na(score_vec)),
    n_scored = n_scored,
    n_score_0 = n_zero,
    n_score_1 = n_one,
    pct_score_0_among_scored = if (n_scored > 0L) 100 * n_zero / n_scored else NA_real_,
    pct_score_1_among_scored = if (n_scored > 0L) 100 * n_one / n_scored else NA_real_
  )
}))

if (nrow(binary_score_distribution) != 15L ||
    any(binary_score_distribution$n_score_0 +
          binary_score_distribution$n_score_1 !=
          binary_score_distribution$n_scored)) {
  stop("STOP: Component-level binary score distribution failed reconciliation.")
}

write.csv(
  binary_score_distribution,
  file.path(AUDIT_DIR, "RWDID_binary_score_distribution_training.csv"),
  row.names = FALSE,
  na = ""
)

total_score_observed <- score_data$RWDID_total_score_complete15[
  !is.na(score_data$RWDID_total_score_complete15)
]
total_score_counts <- tabulate(total_score_observed + 1L, nbins = 16L)

total_score_distribution <- tibble(
  RWDID_total_score = 0:15,
  n = as.integer(total_score_counts),
  percent_among_complete15 = if (length(total_score_observed) > 0L) {
    100 * as.integer(total_score_counts) / length(total_score_observed)
  } else {
    NA_real_
  },
  complete_15_component_n = length(total_score_observed),
  incomplete_component_n = sum(!complete_15_components)
)

write.csv(
  total_score_distribution,
  file.path(AUDIT_DIR, "RWDID_total_score_distribution_training.csv"),
  row.names = FALSE,
  na = ""
)

# 6B.4 Construct the body-style 0/1 scoring scheme figure.
# The three columns are Low / Medium / High. For one-sided rules, adjacent
# columns with the same score are merged, matching the supplied example.
make_scoring_scheme_cells <- function(i) {
  food_var <- FOOD_MAP$food[i]
  r <- food_rcs_range_all[
    match(food_var, food_rcs_range_all$food), ,
    drop = FALSE
  ]
  m <- rcs_model_summary[
    match(food_var, rcs_model_summary$food), ,
    drop = FALSE
  ]

  direction <- as.character(r$range_direction[1])
  lower <- as.numeric(r$score_lower[1])
  upper <- as.numeric(r$score_upper[1])
  x_min <- as.numeric(m$overall_training_min[1])
  x_max <- as.numeric(m$overall_training_max[1])
  y <- nrow(FOOD_MAP) - i + 1L

  tol <- max(1e-8, abs(x_max - x_min) * 1e-8, abs(lower) * 1e-8)
  cells <- list()

  add_cell <- function(xmin, xmax, score, label) {
    cells[[length(cells) + 1L]] <<- tibble(
      food = food_var,
      display_name = FOOD_MAP$display_name[i],
      row_y = y,
      xmin = xmin,
      xmax = xmax,
      ymin = y - 0.44,
      ymax = y + 0.44,
      score = as.integer(score),
      label = label
    )
  }

  if (direction == "point") {
    point_value <- lower
    # For x = 0 rules (e.g., Sweetened beverages), match the publication-style
    # display: a single left cell labelled 0 scores 1, and all positive intake
    # values are merged into one score-0 block.
    if (isTRUE(abs(point_value - x_min) <= tol) && isTRUE(abs(point_value) <= tol)) {
      add_cell(0, 1, 1L, "0")
      add_cell(1, 3, 0L, ">0")
    } else {
      # Generic exact-value rule fallback.
      add_cell(
        0, 3, 0L,
        paste0("Only x = ", fmt_num(point_value), " scores 1")
      )
    }
  } else if (is.infinite(upper)) {
    if (lower <= x_min + tol) {
      add_cell(0, 3, 1L, "All observed intake")
    } else {
      add_cell(0, 1, 0L, paste0("< ", fmt_num(lower)))
      add_cell(1, 3, 1L, paste0(">= ", fmt_num(lower)))
    }
  } else if (lower <= x_min + tol && upper >= x_max - tol) {
    add_cell(0, 3, 1L, "All observed intake")
  } else if (lower <= x_min + tol) {
    add_cell(0, 2, 1L, paste0("≤ ", fmt_num(upper)))
    add_cell(2, 3, 0L, paste0("> ", fmt_num(upper)))
  } else if (upper >= x_max - tol) {
    add_cell(0, 1, 0L, paste0("< ", fmt_num(lower)))
    add_cell(1, 3, 1L, paste0(">= ", fmt_num(lower)))
  } else {
    add_cell(0, 1, 0L, paste0("< ", fmt_num(lower)))
    add_cell(1, 2, 1L, paste0(fmt_num(lower), "–", fmt_num(upper)))
    add_cell(2, 3, 0L, paste0("> ", fmt_num(upper)))
  }

  bind_rows(cells)
}

scoring_scheme_cells <- bind_rows(lapply(seq_len(nrow(FOOD_MAP)), make_scoring_scheme_cells)) %>%
  mutate(score_factor = factor(score, levels = c(0, 1), labels = c("0 point", "1 point")))

if (nrow(scoring_scheme_cells) < 15L ||
    any(!scoring_scheme_cells$score %in% c(0L, 1L)) ||
    any(!is.finite(scoring_scheme_cells$xmin)) ||
    any(!is.finite(scoring_scheme_cells$xmax))) {
  stop("STOP: The 15-component binary scoring scheme could not be constructed.")
}

# Use a light 28% PiYG hue fraction (72% white) so black score labels remain clear.
mix_with_white <- function(hex, colour_fraction = 0.28) {
  if (!is.finite(colour_fraction) ||
      colour_fraction < 0 || colour_fraction > 1) {
    stop("STOP: colour_fraction must be between 0 and 1.")
  }
  rgb <- grDevices::col2rgb(hex) / 255
  mixed <- 1 - colour_fraction * (1 - rgb)
  grDevices::rgb(mixed[1], mixed[2], mixed[3])
}

# The user-selected PiYG palette is available under the brewer series in
# cols4all. Resolve its full name explicitly against the installed version.
pi_candidates <- c("brewer.pi_yg", "brewer.piyg", "brewer.PiYG", "piyg", "PiYG")
pi_available <- cols4all::c4a_palettes(type = "div")
pi_found <- pi_candidates[pi_candidates %in% pi_available]
if (length(pi_found) == 0L) {
  stop("STOP: Installed cols4all lacks PiYG. Update cols4all; do not silently substitute another palette.")
}
pi_name <- pi_found[1]
pi_raw <- cols4all::c4a(pi_name, n = 11)
if (length(pi_raw) != 11L) {
  stop("STOP: Unexpected cols4all PiYG length: ", length(pi_raw))
}
# Distinct pink/green categorical anchors, softened for legible black text.
SCORE_FILL_COLORS <- c(
  "0" = mix_with_white(pi_raw[2], 0.28),
  "1" = mix_with_white(pi_raw[10], 0.28)
)
write.csv(
  tibble(score = c(0L, 1L), source_palette = pi_name,
         source_hex = pi_raw[c(2, 10)], displayed_hex = unname(SCORE_FILL_COLORS)),
  file.path(AUDIT_DIR, "Figure_binary_scoring_PiYG_palette_audit.csv"),
  row.names = FALSE
)

point_marker_data <- food_rcs_range_all %>%
  filter(range_direction == "point") %>%
  left_join(
    rcs_model_summary %>%
      select(food, overall_training_min, overall_training_max),
    by = "food"
  ) %>%
  left_join(
    FOOD_MAP %>% select(food, display_name),
    by = "food"
  ) %>%
  mutate(
    row_y = nrow(FOOD_MAP) - match(food, FOOD_MAP$food) + 1L,
    marker_x = if_else(
      is.finite(overall_training_min) &
        is.finite(overall_training_max) &
        overall_training_max > overall_training_min,
      0.5 + 2 * (value1 - overall_training_min) /
        (overall_training_max - overall_training_min),
      1.5
    )
  )

p_binary_scoring_scheme <- ggplot(scoring_scheme_cells) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = factor(score)
    ),
    colour = "white",
    linewidth = 0.9
  ) +
  geom_point(
    data = point_marker_data,
    aes(x = marker_x, y = row_y),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.6,
    stroke = 0.9,
    fill = SCORE_FILL_COLORS[["1"]],
    colour = "#243B53"
  ) +
  geom_text(
    aes(
      x = (xmin + xmax) / 2,
      y = row_y,
      label = label
    ),
    size = 4.0,
    colour = "black",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = SCORE_FILL_COLORS,
    breaks = c("0", "1"),
    labels = c("0 point", "1 point"),
    name = "Score:"
  ) +
  scale_x_continuous(
    breaks = c(0.5, 1.5, 2.5),
    labels = c("Low", "Medium", "High"),
    limits = c(0, 3),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(nrow(FOOD_MAP)),
    labels = rev(FOOD_MAP$display_name),
    expand = expansion(add = c(0.55, 0.75))
  ) +
  labs(
    title = "Consumption level (servings per day)",
    x = "Grade of consumption level",
    y = NULL,
    tag = "b"
  ) +
  guides(
    fill = guide_legend(
      direction = "horizontal",
      title.position = "left",
      title.hjust = 0,
      override.aes = list(colour = NA)
    )
  ) +
  theme_classic(base_size = 13.5) +
  theme(
    axis.text.x = element_text(colour = "black", size = 13.5, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12.75, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13.5, face = "bold"),
    legend.text = element_text(size = 13.5, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17.5),
    plot.tag = element_text(face = "bold", size = 21.5),
    plot.tag.position = c(0.005, 1.02),
    plot.margin = margin(18, 18, 14, 26)
  )

ggsave(
  filename = file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.pdf"),
  plot = p_binary_scoring_scheme,
  width = 8.33,
  height = 9.2,
  device = "pdf"
)
ggsave(
  filename = file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.png"),
  plot = p_binary_scoring_scheme,
  width = 8.33,
  height = 9.2,
  dpi = 600,
  bg = "white"
)

# 6B.5 Revised journal-style methods workflow: evidence-pattern branches,
# tail/bootstrapping audit, Cox/RCS reconciliation, and frozen binary scoring.
# Colours and text sizes follow accessible, high-contrast print guidelines.
flow_boxes <- tibble(
  id = c("cox", "rcs", "type", "u", "monotone", "unsupported",
         "audit", "freeze"),
  xmin = c(0.3,  6.4,  0.7, 0.3, 4.45, 8.6, 1.0, 2.1),
  xmax = c(5.9, 12.0, 11.6, 4.0, 8.15, 12.3,11.6,10.5),
  ymin = c(8.85,8.85,7.1,4.85,4.85,4.85,2.60,0.35),
  ymax = c(10.45,10.45,8.5,6.40,6.40,6.40,4.15,1.95),
  fill = c("#F1F5FA", "#EFF7F6", "#F7F8FA", "#ECF6F3",
           "#EFF5FB", "#F6F6F7", "#F7F8FA", "#EAF3F7"),
  edge = c("#77AADD", "#009988", "#6B7C8B", "#009988",
           "#77AADD", "#87919A", "#778899", "#3F657A")
)
flow_headings <- tibble(
  x = c(3.10,9.20,6.15,2.15,6.30,10.45,6.30,6.30),
  y = c(10.08,10.08,8.15,6.05,6.05,6.05,3.78,1.57),
  label = c(
    "01  CATEGORICAL COX", "02  RESTRICTED CUBIC SPLINE",
    "03  CLASSIFY ASSOCIATION EVIDENCE",
    "A  SUPPORTED U-SHAPE", "B  MONOTONIC PATTERN",
    "C  NOT CONFIRMATORY", "04  UNCERTAINTY & TAIL AUDIT",
    "05  FREEZE BINARY COMPONENT RULES"
  )
)
flow_text <- tibble(
  x = flow_headings$x,
  y = c(9.44,9.42,7.58,5.46,5.42,5.40,3.10,0.94),
  label = c(
    "Lowest-risk Ref. category and observed range\nKeep category median as RCS reference",
    "Fit 3-knot RCS; show P10 / P50 / P90 on x-axis\nExtract 5%, 10%, 20% near-minimum plateaus",
    "Nonlinear P < 0.05 for supported U / J patterns\nOverall P alone does not establish non-linearity",
    "Compare frozen anchor with Cox category and RCS plateau\nReport overlap and pointwise 95% CI <1 status",
    "Report monotonic direction and tail support\nDo not impose an unsupported finite cutoff",
    "No significant protective segment?\nRetain anchor as operational rule; do not claim benefit",
    paste0(
      "Case-bootstrap (B=", BOOTSTRAP_B,
      ") of current-model candidate boundaries\n",
      "P90-P99 sample/event counts and P99-vs-P90 HR (95% CI)"
    ),
    "1 = inside frozen anchor; 0 = outside; NA = missing\nAnchor unchanged; updated beta weights applied downstream"
  )
)
flow_arrows <- tibble(
  x = c(3.1, 9.2, 6.15, 6.15),
  y = c(8.85, 8.85, 4.85, 2.60),
  xend = c(4.45, 7.85, 6.15, 6.15),
  yend = c(8.5, 8.5, 4.15, 1.95)
)
flow_branches <- tibble(
  x = c(2.15,6.30,10.45),
  y = rep(7.1, 3),
  xend = c(2.15,6.30,10.45),
  yend = rep(6.4, 3)
)
flow_join <- tibble(
  x = c(2.15,6.30,10.45), y = rep(4.85, 3),
  xend = rep(6.30,3), yend = rep(4.15,3)
)
p_interval_algorithm <- ggplot() +
  geom_segment(
    data = bind_rows(flow_arrows, flow_branches, flow_join),
    aes(x=x, y=y, xend=xend, yend=yend),
    colour = "#64748B", linewidth = 0.62,
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  ) +
  geom_rect(
    data = flow_boxes,
    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
        fill=fill, colour=edge),
    linewidth = 0.75
  ) +
  geom_text(
    data = flow_headings,
    aes(x=x, y=y, label=label),
    size = 3.45, fontface="bold", colour="#1F3445"
  ) +
  geom_text(
    data = flow_text,
    aes(x=x, y=y, label=label),
    size = 3.15, lineheight=1.12, colour="#253746"
  ) +
  scale_fill_identity() + scale_colour_identity() +
  coord_cartesian(xlim=c(0,12.6), ylim=c(0,10.8), expand=FALSE, clip="off") +
  labs(
    title = "Validating frozen component-specific 0/1 intake rules",
    subtitle = "Frozen 15-component range anchors + current Cox/RCS validation + uncertainty and tail audits",
    caption = paste0(
      "Only the frozen Step_03_1_0 70% training cohort is used for rule derivation. ",
      "Bootstrap uncertainty and tail support are reported, not used to tune the rule to the test set."
    )
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust=.5, face="bold", size=18,
                              colour="#1F2937", margin=margin(b=6)),
    plot.subtitle = element_text(hjust=.5, size=10.5,
                                 colour="#475569", margin=margin(b=14)),
    plot.caption = element_text(hjust=.5, size=9.1,
                                colour="#526579", margin=margin(t=12)),
    plot.margin=margin(18,18,18,18)
  )

ggsave(
  file.path(PLOT_DIR,"Figure_Sx_interval_derivation_schema.pdf"),
  p_interval_algorithm, width=12, height=8.3, device="pdf"
)
ggsave(
  file.path(PLOT_DIR,"Figure_Sx_interval_derivation_schema.png"),
  p_interval_algorithm, width=12, height=8.3, dpi=600, bg="white"
)

# Export the compact rules after the full interval table has been created.
write.csv(
  binary_scoring_rules,
  file.path(OUT_DIR, "RWDID_binary_scoring_rules_15_components.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

message("[6B/7] Interval derivation table, binary scoring files, and figures exported.")

# -----------------------------------------------------------------------------
# 6B.4A Sweetened-beverage frozen-anchor reconciliation
# -----------------------------------------------------------------------------
# In V5.3, sweetened_bev is one of the six components whose frozen V5.1
# bounded interval is intersected with the selected RCS20 minimum-risk plateau.
# It therefore remains a finite "middle" interval (not a forced one-sided
# low-intake rule). This audit checks that the exact imported V5.3 bounds,
# participant scores, score_rule and scoring-scheme labels all agree.
#
# This is a validation-only block: it does not alter either boundary or score.
# -----------------------------------------------------------------------------

same_numeric_or_na <- function(a, b, tolerance = 1e-12) {
  if (length(a) != 1L || length(b) != 1L) return(FALSE)
  if (is.na(a) && is.na(b)) return(TRUE)
  if (is.na(a) || is.na(b)) return(FALSE)
  if (is.infinite(a) || is.infinite(b)) return(isTRUE(a == b))
  is.finite(a) && is.finite(b) && abs(a - b) <= tolerance
}

sweet_anchor <- range_anchor[
  range_anchor$food == "sweetened_bev", , drop = FALSE
]
sweet_final <- food_rcs_range_all[
  food_rcs_range_all$food == "sweetened_bev", , drop = FALSE
]
sweet_cells <- scoring_scheme_cells[
  scoring_scheme_cells$food == "sweetened_bev", , drop = FALSE
]
sweet_model <- rcs_model_summary[
  rcs_model_summary$food == "sweetened_bev", , drop = FALSE
]

if (nrow(sweet_anchor) != 1L || nrow(sweet_final) != 1L ||
    nrow(sweet_model) != 1L ||
    !"score_sweetened_bev" %in% names(score_data)) {
  stop(
    "STOP: Cannot perform the sweetened-beverage frozen-anchor reconciliation; ",
    "the anchor, final range, model summary or participant score column is missing."
  )
}

sweet_anchor_direction <- as.character(sweet_anchor$anchor_direction[1])
sweet_anchor_value1 <- as.numeric(sweet_anchor$anchor_value1[1])
sweet_anchor_value2 <- as.numeric(sweet_anchor$anchor_value2[1])
sweet_final_value1 <- as.numeric(sweet_final$value1[1])
sweet_final_value2 <- as.numeric(sweet_final$value2[1])
sweet_final_direction <- as.character(sweet_final$range_direction[1])
sweet_model_min <- as.numeric(sweet_model$overall_training_min[1])
sweet_model_max <- as.numeric(sweet_model$overall_training_max[1])

sweet_anchor_structure_pass <-
  identical(sweet_anchor_direction, "middle") &&
  is.finite(sweet_anchor_value1) &&
  is.finite(sweet_anchor_value2) &&
  sweet_anchor_value2 > sweet_anchor_value1

sweet_final_anchor_match <-
  identical(sweet_final_direction, sweet_anchor_direction) &&
  same_numeric_or_na(sweet_final_value1, sweet_anchor_value1) &&
  same_numeric_or_na(sweet_final_value2, sweet_anchor_value2)

sweet_expected_rule <- if (sweet_anchor_structure_pass) {
  paste0(
    fmt_num(sweet_anchor_value1), " <= x <= ",
    fmt_num(sweet_anchor_value2)
  )
} else {
  NA_character_
}
sweet_rule_match <- sweet_anchor_structure_pass &&
  identical(as.character(sweet_final$score_rule[1]), sweet_expected_rule)

sweet_exposure <- suppressWarnings(
  as.numeric(as.character(ukb_data_depression_rcs$sweetened_bev))
)
sweet_expected_score <- score_one_component(
  x = sweet_exposure,
  direction = sweet_anchor_direction,
  value1 = sweet_anchor_value1,
  value2 = sweet_anchor_value2,
  food_name = "sweetened_bev"
)
sweet_observed_score <- as.integer(score_data$score_sweetened_bev)
sweet_score_match <- length(sweet_expected_score) == length(sweet_observed_score) &&
  all(
    (is.na(sweet_expected_score) & is.na(sweet_observed_score)) |
      (!is.na(sweet_expected_score) & !is.na(sweet_observed_score) &
         sweet_expected_score == sweet_observed_score)
  )

# Reconstruct the expected display cells from the frozen middle interval and
# the same complete-case training support used by make_scoring_scheme_cells().
sweet_scheme_tol <- max(
  1e-8,
  abs(sweet_model_max - sweet_model_min) * 1e-8,
  abs(sweet_anchor_value1) * 1e-8
)

if (sweet_anchor_structure_pass &&
    is.finite(sweet_model_min) && is.finite(sweet_model_max) &&
    sweet_model_max >= sweet_model_min) {
  if (sweet_anchor_value1 <= sweet_model_min + sweet_scheme_tol &&
      sweet_anchor_value2 >= sweet_model_max - sweet_scheme_tol) {
    sweet_expected_labels <- "All observed intake"
    sweet_expected_cell_scores <- 1L
  } else if (sweet_anchor_value1 <= sweet_model_min + sweet_scheme_tol) {
    sweet_expected_labels <- c(
      paste0("≤ ", fmt_num(sweet_anchor_value2)),
      paste0("> ", fmt_num(sweet_anchor_value2))
    )
    sweet_expected_cell_scores <- c(1L, 0L)
  } else if (sweet_anchor_value2 >= sweet_model_max - sweet_scheme_tol) {
    sweet_expected_labels <- c(
      paste0("< ", fmt_num(sweet_anchor_value1)),
      paste0(">= ", fmt_num(sweet_anchor_value1))
    )
    sweet_expected_cell_scores <- c(0L, 1L)
  } else {
    sweet_expected_labels <- c(
      paste0("< ", fmt_num(sweet_anchor_value1)),
      paste0(
        fmt_num(sweet_anchor_value1), "–",
        fmt_num(sweet_anchor_value2)
      ),
      paste0("> ", fmt_num(sweet_anchor_value2))
    )
    sweet_expected_cell_scores <- c(0L, 1L, 0L)
  }
} else {
  sweet_expected_labels <- character(0)
  sweet_expected_cell_scores <- integer(0)
}

sweet_scheme_match <-
  length(sweet_expected_labels) > 0L &&
  nrow(sweet_cells) == length(sweet_expected_labels) &&
  identical(as.character(sweet_cells$label), sweet_expected_labels) &&
  identical(as.integer(sweet_cells$score), sweet_expected_cell_scores)

sweetened_bev_qc_pass <- all(c(
  sweet_anchor_structure_pass,
  sweet_final_anchor_match,
  sweet_rule_match,
  sweet_score_match,
  sweet_scheme_match
))

sweetened_bev_anchor_qc <- tibble(
  food = "sweetened_bev",
  anchor_direction = sweet_anchor_direction,
  anchor_value1 = sweet_anchor_value1,
  anchor_value2 = sweet_anchor_value2,
  final_direction = sweet_final_direction,
  final_value1 = sweet_final_value1,
  final_value2 = sweet_final_value2,
  final_score_rule = as.character(sweet_final$score_rule[1]),
  expected_score_rule = sweet_expected_rule,
  anchor_is_V5_3_bounded_middle_interval = sweet_anchor_structure_pass,
  final_range_matches_frozen_anchor = sweet_final_anchor_match,
  score_rule_matches_frozen_anchor = sweet_rule_match,
  participant_scores_match_frozen_anchor = sweet_score_match,
  scoring_scheme_matches_frozen_anchor = sweet_scheme_match,
  scoring_scheme_labels = paste(as.character(sweet_cells$label), collapse = " | "),
  scoring_scheme_scores = paste(as.integer(sweet_cells$score), collapse = " | "),
  final_pass = sweetened_bev_qc_pass
)

write.csv(
  sweetened_bev_anchor_qc,
  file.path(AUDIT_DIR, "sweetened_bev_frozen_anchor_reconciliation_QC.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

# -----------------------------------------------------------------------------
# 7. Final QC
# -----------------------------------------------------------------------------

expected_prediction_files <- file.path(
  PRED_DIR, paste0(FOOD_MAP$food, "_rcs3.csv")
)
expected_plot_files <- file.path(
  PLOT_DIR, paste0(FOOD_MAP$food, "_rcs3.pdf")
)

final_qc <- tibble(
  check = c(
    "training_RDS_matches_frozen_Step03_manifest",
    "test_participants_excluded",
    "exactly_15_Model2_selected_components",
    "green_leafy_veg_excluded",
    "high_fat_dairy_excluded",
    "fruit_juices_included",
    "water_included",
    "one_rotating_lowest_risk_reference_per_component",
    "all_15_RCS_models_fitted",
    "all_15_reference_HR_equal_1",
    "all_15_P10_P50_P90_values_exported",
    "all_15_three_distinct_knots",
    "all_15_frozen_anchor_ranges_exported",
    "all_15_prediction_CSVs_written",
    "all_15_individual_PDFs_written",
    "no_RCS_refinement_without_statistical_support",
    "all_15_RCS_plots_show_P10_P50_P90_and_legend",
    "all_15_x_axes_include_numeric_P10_P50_P90_ticks",
    "all_15_x_axes_use_data_adaptive_prediction_support",
    "percentile_plot_annotation_audit_written",
    "frozen_range_anchor_QC_written",
    "frozen_anchor_comparison_and_downstream_range_written",
    "combined_15_component_PDF_written_with_shared_percentile_legend",
    "interval_range_PDF_written",
    "tile_range_PDF_written",
    "15_food_interval_derivation_table_written",
    "15_food_binary_scoring_rules_written",
    "participant_level_training_binary_scores_written",
    "component_0_1_distribution_audit_written",
    "complete_15_component_total_score_distribution_written",
    "interval_derivation_algorithm_PDF_written",
    "interval_derivation_algorithm_PNG_written",
    "binary_scoring_scheme_PDF_written",
    "binary_scoring_scheme_PNG_written",
    "single_value_intervals_rendered_with_valid_point_markers",
    "sweetened_beverages_use_frozen_V5_3_RCS20_bounded_overlap_anchor",
    "sweetened_bev_frozen_anchor_reconciliation_audit_written",
    "Word_interval_derivation_table_written",
    "word_integrated_report_tables_have_15_components_each",
    "word_integrated_report_CI_and_tail_audit_QC_passed",
    "word_integrated_section_display_CSVs_written",
    "combined_panel_Pvalue_layers_match_model_summary_exactly",
    "PiYG_palette_audit_written_and_valid",
    "open_upper_tail_figure_labels_use_ASCII_greater_equal",
    "15_bootstrap_boundary_audits_complete",
    "parallel_5_10_20_bootstrap_sensitivity_summary_written",
    "parallel_bootstrap_replicates_long_written",
    "20pct_bootstrap_boundary_CIs_estimable_for_all_15",
    "bootstrap_boundary_stability_figure_written",
    "word_parallel_bootstrap_sensitivity_table_written",
    "bootstrap_success_rate_at_least_95pct_for_all_15",
    "15_tail_support_audits_complete",
    "all_non_linear_refinements_have_significant_non_linear_P",
    "all_15_binary_scoring_rules_defined",
    "all_15_interval_derivation_rows_reconciled",
    "final_range_export_written"
  ),
  status = c(
    TRUE,
    TRUE,
    nrow(FOOD_MAP) == 15L,
    !("green_leafy_veg" %in% FOOD_MAP$food),
    !("high_fat_dairy" %in% FOOD_MAP$food),
    "fruit_juices" %in% FOOD_MAP$food,
    "water" %in% FOOD_MAP$food,
    nrow(rcs_reference_audit) == 15L,
    nrow(rcs_model_summary) == 15L,
    all(abs(rcs_model_summary$reference_HR - 1) <= 1e-6),
    all(is.finite(rcs_model_summary$p10) &
          is.finite(rcs_model_summary$p50) &
          is.finite(rcs_model_summary$p90)),
    all(is.finite(rcs_model_summary$rcs_knot_1) &
          is.finite(rcs_model_summary$rcs_knot_2) &
          is.finite(rcs_model_summary$rcs_knot_3) &
          rcs_model_summary$rcs_knot_1 < rcs_model_summary$rcs_knot_2 &
          rcs_model_summary$rcs_knot_2 < rcs_model_summary$rcs_knot_3),
    nrow(food_rcs_range_all) == 15L &&
      all(is.finite(food_rcs_range_all$value1)) &&
      isTRUE(all.equal(
        as.numeric(food_rcs_range_all$value1),
        as.numeric(range_anchor$anchor_value1),
        tolerance = 1e-12, check.attributes = FALSE
      )) &&
      isTRUE(all.equal(
        as.numeric(food_rcs_range_all$value2),
        as.numeric(range_anchor$anchor_value2),
        tolerance = 1e-12, check.attributes = FALSE
      )),
    all(file.exists(expected_prediction_files)),
    all(file.exists(expected_plot_files)),
    all(
      !food_rcs_range_all$candidate_refined_by_RCS_5pct |
        rcs_model_summary$rcs_range_eligible
    ),
    nrow(rcs_percentile_plot_audit) == 15L &&
      all(rcs_percentile_plot_audit$percentile_lines_dashed) &&
      all(rcs_percentile_plot_audit$numeric_labels_at_top) &&
      all(rcs_percentile_plot_audit$percentile_values_on_x_axis) &&
      all(rcs_percentile_plot_audit$legend_included) &&
      all(rcs_percentile_plot_audit$all_percentiles_within_displayed_x_axis),
    nrow(rcs_percentile_plot_audit) == 15L &&
      all(rcs_percentile_plot_audit$percentile_values_on_x_axis) &&
      all(nzchar(rcs_percentile_plot_audit$x_axis_percentile_labels)),
    nrow(rcs_percentile_plot_audit) == 15L &&
      all(rcs_percentile_plot_audit$x_axis_method ==
            paste0(
              "data_adaptive_prediction_grid_percentiles_reference_and_interval_values_with_P10_P50_P90_axis_ticks_",
              sprintf("%d", round(100 * RCS_X_AXIS_PADDING_PROP)),
              "pct_right_padding"
            )) &&
      all(rcs_percentile_plot_audit$displayed_x_max >
            rcs_percentile_plot_audit$prediction_grid_max) &&
      all(rcs_percentile_plot_audit$all_prediction_values_within_displayed_x_axis),
    file.exists(file.path(AUDIT_DIR, "RCS_percentile_plot_annotation_QC.csv")),
    file.exists(file.path(AUDIT_DIR, "frozen_15_component_range_anchor_QC.csv")) &&
      nrow(range_anchor) == 15L,
    file.exists(file.path(AUDIT_DIR, "frozen_anchor_vs_current_Cox_RCS_comparison.csv")) &&
      file.exists(file.path(OUT_DIR, "food_rcs_range_all_15comp_downstream.csv")),
    file.exists(file.path(PLOT_DIR, "rcs3_all_food_combined.pdf")),
    file.exists(file.path(PLOT_DIR, "p_food_group.pdf")),
    file.exists(file.path(PLOT_DIR, "p_food_group_tiles.pdf")),
    file.exists(file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.csv")) &&
      nrow(interval_derivation_table) == 15L,
    file.exists(file.path(OUT_DIR, "RWDID_binary_scoring_rules_15_components.csv")) &&
      nrow(binary_scoring_rules) == 15L,
    file.exists(file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv")) &&
      nrow(score_data) == nrow(ukb_data_depression_rcs),
    file.exists(file.path(AUDIT_DIR, "RWDID_binary_score_distribution_training.csv")) &&
      nrow(binary_score_distribution) == 15L,
    file.exists(file.path(AUDIT_DIR, "RWDID_total_score_distribution_training.csv")) &&
      nrow(total_score_distribution) == 16L,
    file.exists(file.path(PLOT_DIR, "Figure_Sx_interval_derivation_schema.pdf")),
    file.exists(file.path(PLOT_DIR, "Figure_Sx_interval_derivation_schema.png")),
    file.exists(file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.pdf")),
    file.exists(file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.png")),
    all(vapply(
      which(food_rcs_range_all$range_direction == "point"),
      function(j) {
        food_var <- as.character(food_rcs_range_all$food[j])
        sum(point_marker_data$food == food_var) == 1L
      },
      logical(1)
    )),
    sweetened_bev_qc_pass,
    file.exists(file.path(AUDIT_DIR, "sweetened_bev_frozen_anchor_reconciliation_QC.csv")) &&
      isTRUE(sweetened_bev_anchor_qc$final_pass[1]),
    file.exists(file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.docx")) &&
      file.info(file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.docx"))$size > 10000,
    nrow(word_table_data) == 15L && nrow(word_bootstrap_data) == 15L &&
      nrow(word_parallel_bootstrap_data) == 15L && nrow(word_tail_data) == 15L &&
      identical(as.character(report_data$food), as.character(FOOD_MAP$food)),
    nrow(report_qc) == 10L && all(report_qc$status) &&
      file.exists(file.path(AUDIT_DIR, "RCS_integrated_Word_report_QC.csv")),
    all(file.exists(file.path(
      OUT_DIR, c("Table_Sx_interval_derivation_display.csv",
                 "Table_Sx_bootstrap_stability.csv",
                 "Table_Sx_parallel_bootstrap_sensitivity.csv",
                 "Table_Sx_bootstrap_stability_sensitivity.csv",
                 "Table_Sx_high_intake_tail.csv")
    ))),
    all(vapply(seq_len(nrow(FOOD_MAP)), function(i) {
      food_var <- FOOD_MAP$food[i]
      p <- plot_list[[food_var]]
      k <- attr(p, "frozen_pvalue_layer_index")
      z <- p$layers[[k]]$data
      m <- rcs_model_summary[rcs_model_summary$food == food_var, , drop=FALSE]
      nrow(z)==1L && nrow(m)==1L &&
        identical(as.character(z$food[1]), food_var) &&
        abs(z$p_overall[1]-m$p_overall[1])<1e-12 &&
        abs(z$p_non_linear[1]-m$p_non_linear[1])<1e-12
    }, logical(1))),
    file.exists(file.path(AUDIT_DIR, "Figure_binary_scoring_PiYG_palette_audit.csv")) &&
      length(SCORE_FILL_COLORS)==2L &&
      grepl("pi_?yg", pi_name, ignore.case=TRUE),
    all(vapply(seq_len(nrow(FOOD_MAP)), function(i) {
      food_var <- FOOD_MAP$food[i]
      r <- food_rcs_range_all[
        food_rcs_range_all$food==food_var, , drop=FALSE]
      if (r$range_direction[1] != "high") return(TRUE)
      zz <- scoring_scheme_cells[
        scoring_scheme_cells$food == food_var &
        scoring_scheme_cells$score == 1L, , drop=FALSE]
      any(grepl(">=", zz$label, fixed=TRUE)) &&
        !any(grepl("≥", zz$label, fixed=TRUE))
    }, logical(1))),
    nrow(bootstrap_boundary_audit)==15L &&
      all(bootstrap_boundary_audit$bootstrap_requested==BOOTSTRAP_B),
    nrow(bootstrap_boundary_sensitivity_summary) == 45L &&
      all(bootstrap_boundary_sensitivity_summary$bootstrap_successful >=
            ceiling(BOOTSTRAP_B * BOOTSTRAP_MIN_SUCCESS_RATE)),
    nrow(bootstrap_boundary_replicates_long) ==
      15L * BOOTSTRAP_B * length(BOOTSTRAP_RANGE_TOLERANCES) &&
      file.exists(file.path(AUDIT_DIR,
        "RCS_bootstrap_boundary_stability_replicates_long.csv")),
    all(is.finite(bootstrap_boundary_audit$bootstrap_plateau_20_lower_ci_lower)) &&
      all(is.finite(bootstrap_boundary_audit$bootstrap_plateau_20_upper_ci_upper)) &&
      all(bootstrap_boundary_audit$bootstrap_plateau_20_lower_estimable_fraction >= 0.80) &&
      all(bootstrap_boundary_audit$bootstrap_plateau_20_upper_estimable_fraction >= 0.80),
    file.exists(file.path(PLOT_DIR,
      "Figure_Sx_bootstrap_boundary_stability_5_10_20.pdf")) &&
      file.exists(file.path(PLOT_DIR,
      "Figure_Sx_bootstrap_boundary_stability_5_10_20.png")),
    file.exists(file.path(OUT_DIR,
      "Table_Sx_parallel_bootstrap_sensitivity.csv")) &&
      nrow(word_parallel_bootstrap_data) == 15L,
    all(bootstrap_boundary_audit$bootstrap_successful >=
          ceiling(BOOTSTRAP_B*BOOTSTRAP_MIN_SUCCESS_RATE)),
    nrow(tail_support_audit)==15L &&
      nrow(tail_support_bins)==60L,
    all(!food_rcs_range_all$candidate_refined_by_RCS_5pct |
          rcs_model_summary$nonlinear_evidence_supported),
    all(!is.na(food_rcs_range_all$score_rule) &
          food_rcs_range_all$score_rule != "not_defined"),
    nrow(interval_derivation_table) == 15L &&
      identical(as.character(interval_derivation_table$food),
                as.character(FOOD_MAP$food)),
    file.exists(file.path(OUT_DIR, "food_rcs_range_all_final.csv"))
  )
)

write.csv(
  final_qc,
  file.path(AUDIT_DIR, "RCS_final_QC.csv"),
  row.names = FALSE
)

if (any(!final_qc$status)) {
  failed_qc <- dplyr::filter(final_qc, !status)
  cat("\\n===== FAILED FINAL QC CHECKS ONLY =====\\n")
  print(failed_qc, n = Inf)
  stop(
    "STOP: ", nrow(failed_qc),
    " of ", nrow(final_qc),
    " final RCS QC checks failed. See failed checks printed above and RCS_final_QC.csv."
  )
}

cat("\n===== RCS MODEL AND INTERVAL SUMMARY =====\n")
print(
  rcs_model_summary %>%
    select(
      group,
      reference_category_label,
      reference_median_used,
      category_range_low,
      category_range_high,
      category_range_direction,
      p_overall,
      p_non_linear,
      curve_shape,
      plateau_05_low,
      plateau_05_high,
      plateau_05_overlaps_reference_category,
      refined_by_RCS_5pct,
      primary_range_value1,
      primary_range_value2,
      primary_range_direction,
      score_rule
    )
)

cat("\n============================================================\n")
cat("Step_04_1_1 RCS FINAL REANALYSIS COMPLETED\n")
cat("Training cohort: exact Step_03_1_0 frozen 70% RDS; no re-split.\n")
cat("All 15 components retained, including Fruit juices and Water.\n")
cat("Primary boundaries exactly match the frozen food_rcs_range_all.csv anchor.\n")
cat("P10/P50/P90: coloured dashed markers, numeric x-axis ticks and top labels on all RCS plots.\n")
cat("Percentile colours: P10 #77AADD; P50 #009988; P90 #9467BD.\n")
cat("Combined RCS figure uses one shared percentile legend.\n")
cat("Binary scoring figure: lightened cols4all PiYG 0/1 palette; frozen V5.3 RCS20 bounded/one-sided rules displayed.\n")
cat("Interval derivation schema: revised with nonlinear/monotonic evidence, bootstrap and tail audits.\n")
cat("Supplementary table: one integrated Word document with A/B1/B2/C tables and detailed footnotes.\n")
cat("Primary intervals: frozen range anchors retained unchanged; current Cox/RCS results are validation diagnostics.\n")
cat("5%/10%/20% near-minimum plateaus are all assessed by parallel B=500 bootstrap; 20% is the primary operational threshold.\n")
cat("Missing component intake remains NA; no imputation.\n")
cat("Parallel 5%/10%/20% bootstrap boundary-stability audit: ", BOOTSTRAP_B,
    " resamples per component; same bootstrap curves evaluated at all three thresholds.\n", sep="")
cat("Bootstrap boundary sensitivity table: ",
    file.path(OUT_DIR, "Table_Sx_parallel_bootstrap_sensitivity.csv"), "\n", sep="")
cat("Bootstrap boundary sensitivity figure: ",
    file.path(PLOT_DIR, "Figure_Sx_bootstrap_boundary_stability_5_10_20.pdf"), "\n", sep="")
cat("Tail audit: P90-P99 persons/events and P99/P90 contrast for all 15 components.\n")
cat("Continuous 0-1 scoring sensitivity analysis is not run in this script.\n")
cat("Word table: ", file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.docx"), "\n", sep = "")
cat("Interval CSV: ", file.path(OUT_DIR, "Table_Sx_interval_derivation_summary.csv"), "\n", sep = "")
cat("Algorithm figure: ", file.path(PLOT_DIR, "Figure_Sx_interval_derivation_schema.pdf"), "\n", sep = "")
cat("Binary scoring figure: ", file.path(PLOT_DIR, "Figure_Sx_binary_scoring_scheme.pdf"), "\n", sep = "")
cat("Training scores: ", file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv"), "\n", sep = "")
cat("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 08/11: Step_04_1_1_4_FINAL_RCS_BINARY
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 09/11: Step_04_1_2_RWDID_C
# Source: Step_04_1_2_RCS_continuous_proportional_scoring_colizzi(2).R
# Provenance: Library source dated 2026-11-23
# SHA256 of included block: ba8963e231d8001b8ba280fe51406c0689782bd4214e7b3eb2fb3e14751f84ea
########################################################################################

# =============================================================================
# Step_04_1_2_CONTINUOUS_PROPORTIONAL_SCORING_COLIZZI_ADAPTED.R
# =============================================================================
# RWDID | anchored Colizzi-adapted continuous proportional scoring
#
# PURPOSE
#   Rebuild the secondary continuous 0-1-per-component score using the FINAL
#   binary scoring intervals exported by Step_04_1_1.
#   file is the sole source of target boundaries and interval directions.
#
# FROZEN DESIGN
#   1. Exact 15-component order from 3.1; no green_leafy_veg or high_fat_dairy.
#   2. Each component is scored continuously from 0 to 1; equal-weight total 0-15.
#   3. interval/direction is the full-score target:
#        high  : full score at/above value1 (adequacy pattern)
#        low   : full score at/below value1 (moderation pattern)
#        middle: full score between value1 and value2 (optimum pattern)
#        point : full score at value1, with a triangular proportional taper
#                towards the training P1 and P99 support boundaries.
#   4. Outside the full-score target, scores change linearly to 0 at the
#      corresponding training-only P1/P99 support boundary.
#   5. Intake is clamped to training P1-P99 before scoring. P1/P99 are estimated
#      from observed food intake in the frozen training cohort only.
#   6. The same frozen scoring functions and P1/P99 values are applied to test.
#   7. Missing intake remains NA; no imputation. Total is calculated only when
#      all 15 component scores are observed.
#   8. This is an equal-weight continuous sensitivity score, NOT the primary
#      binary score and NOT the beta-weighted RWDID.
#   9. No additional HR-based 5%, 10% or 20% tolerance is applied here. The
#      intervals themselves are imported unchanged as score=1 targets.
#  10. The frozen sweetened-beverage rule is a MIDDLE interval:
#      range_direction="middle", value1=0 and value2=2.08, i.e.
#      0 <= intake <= 2.08 receives the full score. Above 2.08, the continuous
#      score declines linearly to 0 at training P99. This is a scoring rule,
#      not a causal or clinical intake recommendation.
#
# RUN FROM UK BIOBANK PROJECT ROOT:
#   Rscript Step_04_1_2_CONTINUOUS_PROPORTIONAL_SCORING_COLIZZI_ADAPTED_FIX3.R
#
# REQUIRED INPUTS (explicit paths; no file searching):
#   ./results/Step_03/ukb_70_30_split/ukb_70_30_split_manifest.csv
#   ./results/Step_03/ukb_70_30_split/main/ukb_depression_training_cohort.rds
#   ./results/Step_03/ukb_70_30_split/main/ukb_depression_test_cohort.rds
#   ./results/Step_04/Step_04_1_1_RCS_FINAL/
#       food_rcs_range_all_final.csv
#       RWDID_binary_scoring_rules_15_components.csv
#       RWDID_binary_component_scores_training.csv
#       audit/RCS_model_summary_and_range_diagnostics.csv
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(tidyverse)
  library(cols4all)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 0. Explicit paths and immutable cohort/component definitions
# -----------------------------------------------------------------------------

SPLIT_DIR <- "./results/Step_03/ukb_70_30_split"
SPLIT_MANIFEST_FILE <- file.path(SPLIT_DIR, "ukb_70_30_split_manifest.csv")
TRAIN_RDS <- file.path(SPLIT_DIR, "main", "ukb_depression_training_cohort.rds")
TEST_RDS <- file.path(SPLIT_DIR, "main", "ukb_depression_test_cohort.rds")

V531_DIR <- "./results/Step_04/Step_04_1_1_4"
RCS_RANGE_FILE <- file.path(V531_DIR, "food_rcs_range_all_final.csv")
RCS_RULE_FILE <- file.path(
  V531_DIR, "RWDID_binary_scoring_rules_15_components.csv"
)
RCS_BINARY_TRAIN_FILE <- file.path(
  V531_DIR, "RWDID_binary_component_scores_training.csv"
)
RCS_DIAG_FILE <- file.path(
  V531_DIR, "audit", "RCS_model_summary_and_range_diagnostics.csv"
)

OUT_DIR <- paste0(
  "./results/Step_04/",
  "Step_04_1_2"
) 
AUDIT_DIR <- file.path(OUT_DIR, "audit")
PLOT_DIR <- file.path(OUT_DIR, "plots")
CURVE_DIR <- file.path(OUT_DIR, "score_curves")
for (p in c(OUT_DIR, AUDIT_DIR, PLOT_DIR, CURVE_DIR)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

EXPECTED_SEED <- 20261103L
EXPECTED_TRAIN_N <- 76578L
EXPECTED_TEST_N <- 32819L
EXPECTED_TRAIN_EVENTS <- 2122L
EXPECTED_TEST_EVENTS <- 875L
N_COMPONENTS <- 15L
NUMERICAL_TOL <- 1e-8
GRID_N <- 501L

FOOD_MAP <- data.frame(
  component_order = seq_len(15L),
  group = c(
    "Citrus fruits", "Red/orange vegetables", "Other vegetables",
    "Potatoes", "Refined grains", "Eggs", "Low-fat dairy", "Red meats",
    "Processed meats", "Fish and seafood", "Coffee", "Fruit juices",
    "Sweetened beverages", "Water", "Wine"
  ),
  food = c(
    "citrus_fruits", "red_orange_veg", "other_veg", "potatoes",
    "refined_grains", "eggs", "low_fat_dairy", "red_meats",
    "processed_meats", "fish_seafood", "coffee", "fruit_juices",
    "sweetened_bev", "water", "wine"
  ),
  display_name = c(
    "Citrus fruits", "Red/Orange vegetables", "Other vegetables",
    "Potatoes", "Refined grains", "Eggs", "Low-fat dairy", "Red meats",
    "Processed meats", "Fish and seafood", "Coffee", "Fruit juices",
    "Sweetened beverages", "Water", "Wine"
  ),
  stringsAsFactors = FALSE
)
FOOD_MAP$score_col <- paste0("score_", FOOD_MAP$food)

if (nrow(FOOD_MAP) != N_COMPONENTS || anyDuplicated(FOOD_MAP$food) ||
    any(c("green_leafy_veg", "high_fat_dairy") %in% FOOD_MAP$food) ||
    !all(c("fruit_juices", "water") %in% FOOD_MAP$food)) {
  stop("STOP: frozen 15-component map is invalid.")
}

# -----------------------------------------------------------------------------
# 1. Strict input and frozen 70/30 split checks
# -----------------------------------------------------------------------------

required_files <- c(
  SPLIT_MANIFEST_FILE, TRAIN_RDS, TEST_RDS,
  RCS_RANGE_FILE, RCS_RULE_FILE, RCS_BINARY_TRAIN_FILE, RCS_DIAG_FILE
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "STOP: Required frozen input file(s) are missing:\n",
    paste(missing_files, collapse = "\n"),
    "\nRun Step_04_1_1 first."
  )
}

manifest <- utils::read.csv(
  SPLIT_MANIFEST_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(participant_id = "character", split = "character",
                 id_column = "character", split_seed = "integer",
                 training_fraction = "numeric")
)
manifest_cols <- c(
  "participant_id", "split", "id_column", "split_seed", "training_fraction"
)
if (!all(manifest_cols %in% names(manifest)) ||
    anyNA(manifest[, manifest_cols]) ||
    anyDuplicated(manifest$participant_id) > 0L ||
    any(!manifest$split %in% c("training", "test")) ||
    length(unique(manifest$id_column)) != 1L ||
    length(unique(manifest$split_seed)) != 1L ||
    unique(manifest$split_seed) != EXPECTED_SEED ||
    length(unique(manifest$training_fraction)) != 1L ||
    abs(unique(manifest$training_fraction) - 0.70) > 1e-12) {
  stop("STOP: Frozen Step_03_1_0 split manifest failed integrity checks.")
}

ID_COL <- unique(manifest$id_column)
train_ids <- sort(
  trimws(as.character(manifest$participant_id[manifest$split == "training"])),
  method = "radix"
)
test_ids <- sort(
  trimws(as.character(manifest$participant_id[manifest$split == "test"])),
  method = "radix"
)
if (length(train_ids) != EXPECTED_TRAIN_N ||
    length(test_ids) != EXPECTED_TEST_N ||
    length(intersect(train_ids, test_ids)) > 0L) {
  stop("STOP: Frozen 70/30 participant counts or separation changed.")
}

training_data <- readRDS(TRAIN_RDS)
test_data <- readRDS(TEST_RDS)
if (!is.data.frame(training_data) || !is.data.frame(test_data) ||
    !ID_COL %in% names(training_data) || !ID_COL %in% names(test_data)) {
  stop("STOP: Frozen training/test RDS or participant ID column is invalid.")
}

check_ids <- function(d, expected, split_name) {
  ids <- trimws(as.character(d[[ID_COL]]))
  if (anyNA(ids) || any(ids == "") || anyDuplicated(ids) > 0L ||
      !identical(sort(ids, method = "radix"), expected)) {
    stop("STOP: ", split_name, " RDS IDs do not match the frozen manifest.")
  }
  if ("analysis_split" %in% names(d) &&
      any(as.character(d$analysis_split) != split_name)) {
    stop("STOP: Non-", split_name, " rows were found in the corresponding RDS.")
  }
  invisible(TRUE)
}
check_ids(training_data, train_ids, "training")
check_ids(test_data, test_ids, "test")

if (nrow(training_data) != EXPECTED_TRAIN_N ||
    nrow(test_data) != EXPECTED_TEST_N ||
    sum(training_data$depression_status == 1, na.rm = TRUE) != EXPECTED_TRAIN_EVENTS ||
    sum(test_data$depression_status == 1, na.rm = TRUE) != EXPECTED_TEST_EVENTS) {
  stop("STOP: Frozen training/test N or event counts differ from expectations.")
}

missing_food_train <- setdiff(FOOD_MAP$food, names(training_data))
missing_food_test <- setdiff(FOOD_MAP$food, names(test_data))
if (length(missing_food_train) || length(missing_food_test)) {
  stop(
    "STOP: Food intake variables missing. Training: ",
    paste(missing_food_train, collapse = ", "),
    "; test: ", paste(missing_food_test, collapse = ", ")
  )
}
for (f in FOOD_MAP$food) {
  if (!is.numeric(training_data[[f]]) || !is.numeric(test_data[[f]])) {
    stop("STOP: Food intake variable must be numeric in both cohorts: ", f)
  }
}

cohort_qc <- data.frame(
  check = c(
    "training_N", "test_N", "training_events", "test_events",
    "split_seed", "training_test_ID_overlap"
  ),
  observed = c(
    nrow(training_data), nrow(test_data),
    sum(training_data$depression_status == 1, na.rm = TRUE),
    sum(test_data$depression_status == 1, na.rm = TRUE),
    unique(manifest$split_seed), length(intersect(train_ids, test_ids))
  ),
  expected = c(
    EXPECTED_TRAIN_N, EXPECTED_TEST_N, EXPECTED_TRAIN_EVENTS,
    EXPECTED_TEST_EVENTS, EXPECTED_SEED, 0
  ),
  stringsAsFactors = FALSE
)
cohort_qc$status <- ifelse(cohort_qc$observed == cohort_qc$expected, "PASS", "FAIL")
utils::write.csv(
  cohort_qc, file.path(AUDIT_DIR, "frozen_70_30_cohort_QC.csv"),
  row.names = FALSE
)
if (any(cohort_qc$status != "PASS")) {
  stop("STOP: Frozen cohort QC failed.")
}

# -----------------------------------------------------------------------------
# 2. Import and reconcile the exact interval definitions
# -----------------------------------------------------------------------------

rcs_ranges <- utils::read.csv(
  RCS_RANGE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
v5_rules <- utils::read.csv(
  RCS_RULE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
range_required <- c(
  "food", "value1", "value2", "range_direction", "range_source",
  "score_lower", "score_upper", "score_rule"
)
rule_required <- c(
  "food", "range_direction", "value1", "value2", "score_lower",
  "score_upper", "score_rule", "score_1", "score_0", "missing_rule"
)
if (!all(range_required %in% names(rcs_ranges))) {
  stop("STOP: range file lacks: ",
       paste(setdiff(range_required, names(rcs_ranges)), collapse = ", "))
}
if (!all(rule_required %in% names(v5_rules))) {
  stop("STOP: scoring-rule file lacks: ",
       paste(setdiff(rule_required, names(v5_rules)), collapse = ", "))
}

check_food_set <- function(d, label) {
  if (nrow(d) != N_COMPONENTS || anyDuplicated(as.character(d$food)) > 0L ||
      !setequal(as.character(d$food), FOOD_MAP$food)) {
    stop("STOP: ", label, " does not contain exactly the frozen 15 foods.")
  }
  d[match(FOOD_MAP$food, as.character(d$food)), , drop = FALSE]
}
rcs_ranges <- check_food_set(rcs_ranges, "range file")
v5_rules <- check_food_set(v5_rules, "binary scoring-rule file")

same_numeric <- function(a, b, tolerance = 1e-10) {
  a <- suppressWarnings(as.numeric(a))
  b <- suppressWarnings(as.numeric(b))
  if (length(a) != length(b)) return(FALSE)
  same_na <- is.na(a) & is.na(b)
  same_finite <- is.finite(a) & is.finite(b) & abs(a - b) <= tolerance
  same_infinite <- is.infinite(a) & is.infinite(b) & sign(a) == sign(b)
  all(same_na | same_finite | same_infinite)
}

rcs_ranges$value1 <- as.numeric(rcs_ranges$value1)
rcs_ranges$value2 <- as.numeric(rcs_ranges$value2)
rcs_ranges$score_lower <- as.numeric(rcs_ranges$score_lower)
rcs_ranges$score_upper <- as.numeric(rcs_ranges$score_upper)
v5_rules$value1 <- as.numeric(v5_rules$value1)
v5_rules$value2 <- as.numeric(v5_rules$value2)
v5_rules$score_lower <- as.numeric(v5_rules$score_lower)
v5_rules$score_upper <- as.numeric(v5_rules$score_upper)

range_rule_match <- identical(
  as.character(rcs_ranges$range_direction), as.character(v5_rules$range_direction)
) &&
  identical(as.character(rcs_ranges$score_rule), as.character(v5_rules$score_rule)) &&
  same_numeric(rcs_ranges$value1, v5_rules$value1) &&
  same_numeric(rcs_ranges$value2, v5_rules$value2) &&
  same_numeric(rcs_ranges$score_lower, v5_rules$score_lower) &&
  same_numeric(rcs_ranges$score_upper, v5_rules$score_upper)
if (!range_rule_match) {
  stop("STOP: range CSV and binary scoring-rule CSV do not reconcile.")
}

if (anyNA(rcs_ranges$range_direction) ||
    any(!rcs_ranges$range_direction %in% c("middle", "low", "high", "point")) ||
    any(!is.finite(rcs_ranges$value1)) ||
    anyNA(rcs_ranges$range_source) || any(trimws(rcs_ranges$range_source) == "") ||
    anyNA(rcs_ranges$score_rule) || any(rcs_ranges$score_rule == "not_defined")) {
  stop("STOP: One or more interval definitions are invalid or incomplete.")
}

middle_rows <- rcs_ranges$range_direction == "middle"
if (any(!is.finite(rcs_ranges$value2[middle_rows])) ||
    any(rcs_ranges$value2[middle_rows] < rcs_ranges$value1[middle_rows])) {
  stop("STOP: Invalid middle interval.")
}
if (any(rcs_ranges$range_direction == "high" &
        !is.infinite(rcs_ranges$score_upper))) {
  stop("STOP: high-direction intervals must retain an open upper tail.")
}
if (any(rcs_ranges$range_direction == "low" &
        !is.na(rcs_ranges$value2))) {
  stop("STOP: low-direction intervals must have value2=NA.")
}

# specific final repairs and the sweetened-beverage 0-2.08
# middle-interval rule are hard gates.
repair_expectations <- data.frame(
  food = c("potatoes", "red_meats"),
  direction = c("middle", "middle"),
  value1 = c(1.00, 1.00),
  value2 = c(1.50, 1.75),
  stringsAsFactors = FALSE
)
repair_idx <- match(repair_expectations$food, rcs_ranges$food)
repair_pass <- !anyNA(repair_idx) &&
  identical(rcs_ranges$range_direction[repair_idx], repair_expectations$direction) &&
  all(abs(rcs_ranges$value1[repair_idx] - repair_expectations$value1) <= 1e-8) &&
  all(abs(rcs_ranges$value2[repair_idx] - repair_expectations$value2) <= 1e-8)
if (!repair_pass) {
  stop("STOP: Potatoes/Red meats repaired intervals do not match 1.00-1.50 and 1.00-1.75.")
}
sweet_idx <- match("sweetened_bev", rcs_ranges$food)
sweetened_rule_text <- if (!is.na(sweet_idx)) {
  trimws(as.character(rcs_ranges$score_rule[sweet_idx]))
} else {
  NA_character_
}

# Parse the literal interval expression, e.g. "0 <= x <= 2.08".
sweetened_rule_match <- if (!is.na(sweetened_rule_text)) {
  regexec(
    "^([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)\\s*<=\\s*x\\s*<=\\s*([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)$",
    sweetened_rule_text,
    perl = TRUE
  )
} else {
  NULL
}
sweetened_rule_capture <- if (!is.null(sweetened_rule_match)) {
  regmatches(sweetened_rule_text, sweetened_rule_match)[[1]]
} else {
  character(0)
}
sweetened_rule_bounds <- if (length(sweetened_rule_capture) == 3L) {
  suppressWarnings(as.numeric(sweetened_rule_capture[2:3]))
} else {
  c(NA_real_, NA_real_)
}

sweetened_middle_pass <- !is.na(sweet_idx) &&
  identical(trimws(as.character(rcs_ranges$range_direction[sweet_idx])), "middle") &&
  is.finite(rcs_ranges$value1[sweet_idx]) &&
  is.finite(rcs_ranges$value2[sweet_idx]) &&
  abs(rcs_ranges$value1[sweet_idx] - 0) <= NUMERICAL_TOL &&
  abs(rcs_ranges$value2[sweet_idx] - 2.08) <= NUMERICAL_TOL &&
  is.finite(rcs_ranges$score_lower[sweet_idx]) &&
  is.finite(rcs_ranges$score_upper[sweet_idx]) &&
  abs(rcs_ranges$score_lower[sweet_idx] - rcs_ranges$value1[sweet_idx]) <= NUMERICAL_TOL &&
  abs(rcs_ranges$score_upper[sweet_idx] - rcs_ranges$value2[sweet_idx]) <= NUMERICAL_TOL &&
  length(sweetened_rule_bounds) == 2L &&
  all(is.finite(sweetened_rule_bounds)) &&
  abs(sweetened_rule_bounds[1] - rcs_ranges$value1[sweet_idx]) <= NUMERICAL_TOL &&
  abs(sweetened_rule_bounds[2] - rcs_ranges$value2[sweet_idx]) <= NUMERICAL_TOL

if (!sweetened_middle_pass) {
  stop(
    "STOP: Sweetened beverages must reconcile to the frozen MIDDLE rule: ",
    "range_direction='middle', value1=0, value2=2.08 and score_rule='0 <= x <= 2.08'. ",
    "Observed direction=",
    ifelse(is.na(sweet_idx), "MISSING", rcs_ranges$range_direction[sweet_idx]),
    "; value1=",
    ifelse(is.na(sweet_idx), "NA", as.character(rcs_ranges$value1[sweet_idx])),
    "; value2=",
    ifelse(is.na(sweet_idx) || is.na(rcs_ranges$value2[sweet_idx]),
           "NA", as.character(rcs_ranges$value2[sweet_idx])),
    "; score_rule=",
    ifelse(is.na(sweetened_rule_text), "NA", sweetened_rule_text)
  )
}

range_rule_audit <- data.frame(
  component_order = seq_len(N_COMPONENTS),
  food = FOOD_MAP$food,
  range_direction = rcs_ranges$range_direction,
  value1 = rcs_ranges$value1,
  value2 = rcs_ranges$value2,
  score_lower = rcs_ranges$score_lower,
  score_upper = rcs_ranges$score_upper,
  score_rule = rcs_ranges$score_rule,
  range_source = rcs_ranges$range_source,
  range_and_rule_files_match = range_rule_match,
  stringsAsFactors = FALSE
)
utils::write.csv(
  range_rule_audit,
  file.path(AUDIT_DIR, "interval_rule_reconciliation.csv"),
  row.names = FALSE, na = ""
)

# -----------------------------------------------------------------------------
# 3. Import RCS evidence for reporting only; it does NOT alter score boundaries
# -----------------------------------------------------------------------------

rcs_diag <- utils::read.csv(
  RCS_DIAG_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
diag_required <- c(
  "food", "p10", "p50", "p90", "p_overall", "p_non_linear", "curve_shape",
  "plateau_20_low", "plateau_20_high"
)
if (!all(diag_required %in% names(rcs_diag))) {
  stop("STOP: RCS diagnostic file lacks: ",
       paste(setdiff(diag_required, names(rcs_diag)), collapse = ", "))
}
rcs_diag <- check_food_set(rcs_diag, "RCS diagnostic file")

# -----------------------------------------------------------------------------
# 4. Derive proportional score functions from intervals + train P1/P99
# -----------------------------------------------------------------------------

support_rows <- vector("list", N_COMPONENTS)
for (i in seq_len(N_COMPONENTS)) {
  f <- FOOD_MAP$food[i]
  x <- as.numeric(training_data[[f]])
  observed <- !is.na(x)
  if (any(!is.finite(x[observed]))) {
    stop("STOP: Non-finite non-missing training intake values for ", f)
  }
  if (any(x[observed] < -NUMERICAL_TOL)) {
    stop("STOP: Negative training intake values for ", f)
  }
  x_obs <- x[observed]
  if (length(x_obs) < 100L) {
    stop("STOP: Insufficient observed training intake values for ", f)
  }
  p1 <- unname(stats::quantile(x_obs, probs = 0.01, type = 7, names = FALSE))
  p99 <- unname(stats::quantile(x_obs, probs = 0.99, type = 7, names = FALSE))
  if (!is.finite(p1) || !is.finite(p99) || p99 <= p1 + NUMERICAL_TOL) {
    stop("STOP: Training P1-P99 support is degenerate for ", f,
         " (P1=", p1, "; P99=", p99, ").")
  }
  support_rows[[i]] <- data.frame(
    food = f, support_p1 = p1, support_p99 = p99,
    n_training_intake_observed = length(x_obs),
    n_training_intake_missing = sum(!observed),
    training_intake_min = min(x_obs), training_intake_max = max(x_obs),
    stringsAsFactors = FALSE
  )
}
support_table <- do.call(rbind, support_rows)

parameters <- FOOD_MAP
parameters <- merge(
  parameters,
  rcs_ranges[, c("food", "value1", "value2", "range_direction",
                 "range_source", "score_lower", "score_upper", "score_rule")],
  by = "food", all.x = TRUE, sort = FALSE
)
parameters <- parameters[match(FOOD_MAP$food, parameters$food), , drop = FALSE]
parameters <- merge(parameters, support_table, by = "food", all.x = TRUE, sort = FALSE)
parameters <- parameters[match(FOOD_MAP$food, parameters$food), , drop = FALSE]
parameters <- merge(
  parameters,
  rcs_diag[, c("food", "p10", "p50", "p90", "p_overall", "p_non_linear",
               "curve_shape", "plateau_20_low", "plateau_20_high")],
  by = "food", all.x = TRUE, sort = FALSE
)
parameters <- parameters[match(FOOD_MAP$food, parameters$food), , drop = FALSE]

parameters$score_type <- c(
  high = "adequacy", low = "moderation", middle = "optimum", point = "point"
)[as.character(parameters$range_direction)]
parameters$target_low <- ifelse(
  parameters$range_direction %in% c("high", "middle", "point"),
  parameters$value1, NA_real_
)
parameters$target_high <- ifelse(
  parameters$range_direction == "middle", parameters$value2,
  ifelse(parameters$range_direction == "low", parameters$value1,
         ifelse(parameters$range_direction == "point", parameters$value1, NA_real_))
)
parameters$target_intersects_training_support <- mapply(
  function(direction, a, b, L, U) {
    if (direction == "high") return(U >= a)
    if (direction == "low") return(L <= a)
    if (direction == "middle") return(b >= L && a <= U)
    if (direction == "point") return(a >= L - NUMERICAL_TOL && a <= U + NUMERICAL_TOL)
    FALSE
  },
  parameters$range_direction, parameters$value1,
  ifelse(is.na(parameters$value2), parameters$value1, parameters$value2),
  parameters$support_p1, parameters$support_p99,
  USE.NAMES = FALSE
)

if (nrow(parameters) != N_COMPONENTS ||
    !identical(as.character(parameters$food), FOOD_MAP$food) ||
    anyNA(parameters$score_type) ||
    any(!is.finite(parameters$support_p1)) ||
    any(!is.finite(parameters$support_p99)) ||
    any(parameters$support_p99 <= parameters$support_p1)) {
  stop("STOP: proportional scoring parameters failed reconciliation.")
}
if (!all(parameters$target_intersects_training_support)) {
  bad <- parameters$food[!parameters$target_intersects_training_support]
  stop("STOP: target interval does not intersect training P1-P99 support for: ",
       paste(bad, collapse = ", "))
}

# Proportional functions. Full score=1 within the target; score tapers
# linearly to 0 at training P1/P99. The same function is applied to train/test.
score_proportional_v531 <- function(intake, prm) {
  x <- suppressWarnings(as.numeric(intake))
  score <- rep(NA_real_, length(x))
  observed <- !is.na(x)
  if (any(!is.finite(x[observed]))) {
    stop("STOP: Non-finite non-missing intake encountered for ", prm$food)
  }
  if (any(x[observed] < -NUMERICAL_TOL)) {
    stop("STOP: Negative intake encountered for ", prm$food)
  }
  if (!any(observed)) return(score)

  L <- as.numeric(prm$support_p1)
  U <- as.numeric(prm$support_p99)
  a <- as.numeric(prm$value1)
  b <- as.numeric(prm$value2)
  direction <- as.character(prm$range_direction)
  z <- pmin(U, pmax(L, x[observed]))
  s <- rep(NA_real_, length(z))

  if (direction == "high") {
    # Adequacy: proportional rise to the lower cutpoint, then plateau.
    if (a <= L + NUMERICAL_TOL) {
      s <- rep(1, length(z))
    } else {
      s <- (z - L) / (a - L)
    }
  } else if (direction == "low") {
    # Moderation: full score for intake <= the low-direction cutoff
    # (value1), then decline linearly above that cutoff. For sweetened
    # beverages, value1 is the upper cutoff; it is NOT assumed to equal zero.
    if (U <= a + NUMERICAL_TOL) {
      s <- rep(1, length(z))
    } else {
      s <- (U - z) / (U - a)
    }
  } else if (direction == "middle") {
    # Optimum: plateau within [value1,value2], with proportional ramps outside.
    if (!is.finite(b) || b < a - NUMERICAL_TOL) {
      stop("STOP: Invalid middle interval for ", prm$food)
    }
    if (b < L - NUMERICAL_TOL || a > U + NUMERICAL_TOL) {
      stop("STOP: middle interval is outside training support for ", prm$food)
    }
    s <- rep(0, length(z))
    inside <- z >= a - NUMERICAL_TOL & z <= b + NUMERICAL_TOL
    s[inside] <- 1
    below <- z < a - NUMERICAL_TOL
    if (any(below) && a > L + NUMERICAL_TOL) {
      s[below] <- (z[below] - L) / (a - L)
    }
    above <- z > b + NUMERICAL_TOL
    if (any(above) && U > b + NUMERICAL_TOL) {
      s[above] <- (U - z[above]) / (U - b)
    }
  } else if (direction == "point") {
    # Generic point target: triangular taper around value1. Sweetened beverages
    # are not encoded as point; their low-direction cutoff is imported unchanged.
    if (a < L - NUMERICAL_TOL || a > U + NUMERICAL_TOL) {
      stop("STOP: point target is outside training support for ", prm$food)
    }
    s <- rep(0, length(z))
    below <- z < a - NUMERICAL_TOL
    if (any(below) && a > L + NUMERICAL_TOL) {
      s[below] <- (z[below] - L) / (a - L)
    }
    above <- z > a + NUMERICAL_TOL
    if (any(above) && U > a + NUMERICAL_TOL) {
      s[above] <- (U - z[above]) / (U - a)
    }
    at_point <- abs(z - a) <= NUMERICAL_TOL
    s[at_point] <- 1
  } else {
    stop("STOP: Unsupported range direction for ", prm$food, ": ", direction)
  }

  score[observed] <- pmin(1, pmax(0, s))
  score
}

# Verify the actual frozen MIDDLE interval for sweetened beverages:
# score=1 throughout [value1,value2] = [0,2.08], then a linear decline to
# zero at training P99 when P99 is above the frozen upper boundary.
sweetened_parameter_idx <- match("sweetened_bev", parameters$food)
if (is.na(sweetened_parameter_idx)) {
  stop("STOP: Sweetened-beverage parameter row is missing.")
}
sweetened_parameter <- parameters[sweetened_parameter_idx, , drop = FALSE]
sweet_a <- as.numeric(sweetened_parameter$value1)
sweet_b <- as.numeric(sweetened_parameter$value2)
sweet_L <- as.numeric(sweetened_parameter$support_p1)
sweet_U <- as.numeric(sweetened_parameter$support_p99)

if (!is.finite(sweet_a) || !is.finite(sweet_b) ||
    abs(sweet_a) > NUMERICAL_TOL ||
    abs(sweet_b - 2.08) > NUMERICAL_TOL ||
    sweet_L > sweet_b + NUMERICAL_TOL) {
  stop("STOP: Sweetened-beverage frozen 0-2.08 interval is incompatible with training support.")
}

if (sweet_U > sweet_b + NUMERICAL_TOL) {
  sweet_mid_above <- (sweet_b + sweet_U) / 2
  sweetened_test_intakes <- c(0, sweet_a, sweet_b, sweet_mid_above, sweet_U)
  sweetened_expected_values <- c(1, 1, 1, 0.5, 0)
} else {
  # The full observed P1-P99 support lies inside the frozen full-score interval.
  sweetened_test_intakes <- c(0, sweet_a, sweet_b, sweet_U)
  sweetened_expected_values <- rep(1, length(sweetened_test_intakes))
}

sweetened_function_values <- score_proportional_v531(
  sweetened_test_intakes, sweetened_parameter
)
sweetened_continuous_function_pass <-
  identical(as.character(sweetened_parameter$range_direction), "middle") &&
  abs(as.numeric(sweetened_parameter$value1) - 0) <= NUMERICAL_TOL &&
  abs(as.numeric(sweetened_parameter$value2) - 2.08) <= NUMERICAL_TOL &&
  length(sweetened_function_values) == length(sweetened_expected_values) &&
  all(is.finite(sweetened_function_values)) &&
  all(abs(sweetened_function_values - sweetened_expected_values) <= NUMERICAL_TOL)

if (!sweetened_continuous_function_pass) {
  stop(
    "STOP: Sweetened-beverage continuous function does not match the frozen ",
    "MIDDLE interval [0,2.08] and the proportional decline above 2.08."
  )
}

sweetened_function_qc <- data.frame(
  food = "sweetened_bev",
  range_direction = as.character(sweetened_parameter$range_direction),
  frozen_interval_lower_value1 = sweet_a,
  frozen_interval_upper_value2 = sweet_b,
  training_p1 = sweet_L,
  training_p99 = sweet_U,
  tested_intakes = paste(format(sweetened_test_intakes, trim = TRUE), collapse = " | "),
  observed_scores = paste(format(sweetened_function_values, trim = TRUE), collapse = " | "),
  expected_scores = paste(format(sweetened_expected_values, trim = TRUE), collapse = " | "),
  function_matches_frozen_middle_interval = sweetened_continuous_function_pass,
  stringsAsFactors = FALSE
)
utils::write.csv(
  sweetened_function_qc,
  file.path(AUDIT_DIR, "sweetened_bev_continuous_function_QC.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# Save parameters before scoring the held-out test set.
parameters$scoring_function <- ifelse(
  parameters$range_direction == "high",
  "Linear adequacy ramp from training P1 to lower cutpoint; then 1",
  ifelse(
    parameters$range_direction == "low",
    "1 through upper cutpoint; linear moderation decline to training P99",
    ifelse(
      parameters$range_direction == "middle",
      "Linear rise to 1 at lower cutpoint; 1 within interval; linear decline to training P99",
      "Point optimum at value1; triangular linear taper to training P1/P99"
    )
  )
)
parameters$outside_support_rule <- "Clamp intake to training-derived P1-P99 before applying the frozen function"
parameters$missing_rule <- "NA; no imputation; total NA unless all 15 components observed"

utils::write.csv(
  parameters,
  file.path(OUT_DIR, "Table_Sx_continuous_scoring_parameters_15.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
utils::write.csv(
  support_table,
  file.path(AUDIT_DIR, "training_only_P1_P99_support.csv"),
  row.names = FALSE
)

frozen_spec <- list(
  method = "anchored Colizzi-adapted proportional continuous scoring",
  version = "Step_04_1_2_CONTINUOUS_PROPORTIONAL_SCORING_COLIZZI_ADAPTED",
  components = FOOD_MAP$food,
  component_order = FOOD_MAP$food,
  equal_component_weights = rep(1, N_COMPONENTS),
  total_range = c(0, 15),
  source_of_full_score_targets = RCS_RANGE_FILE,
  source_of_rule_reconciliation = RCS_RULE_FILE,
  interval_rule = "range_direction/value1/value2 imported unchanged",
  proportional_support = "training cohort observed food intake P1-P99, type=7, by component",
  outside_support = "clamp to training P1-P99",
  missingness = "no imputation; total requires all 15 components",
  point_interval_rule = "triangular proportional taper around the point target",
  sweetened_beverages = "middle interval [0, 2.08] servings/day; score is 1 within the frozen interval and declines linearly above 2.08 to 0 at training P99 when P99 exceeds 2.08",
  no_additional_hr_tolerance = TRUE,
  no_beta_weighting = TRUE,
  test_outcomes_used_to_define_score = FALSE,
  parameters = parameters,
  input_md5 = as.list(tools::md5sum(required_files))
)
saveRDS(
  frozen_spec,
  file.path(OUT_DIR, "RWDID_continuous_scoring_frozen_spec.rds"),
  version = 2
)

# -----------------------------------------------------------------------------
# 5. Apply identical frozen functions to training and held-out test cohorts
# -----------------------------------------------------------------------------

score_dataset <- function(d, split_name, spec, id_col) {
  prm <- spec$parameters
  result <- data.frame(
    participant_id = trimws(as.character(d[[id_col]])),
    split = split_name,
    stringsAsFactors = FALSE
  )
  for (j in seq_len(nrow(prm))) {
    f <- prm$food[j]
    result[[paste0("intake_", f)]] <- as.numeric(d[[f]])
    result[[paste0("score_", f)]] <- score_proportional_v531(
      d[[f]], prm[j, , drop = FALSE]
    )
  }
  score_cols <- paste0("score_", prm$food)
  score_matrix <- as.matrix(result[, score_cols, drop = FALSE])
  storage.mode(score_matrix) <- "numeric"
  result$n_components_observed <- as.integer(rowSums(is.finite(score_matrix)))
  result$RWDID_continuous_total_15 <- rowSums(score_matrix)
  result$RWDID_continuous_total_15[result$n_components_observed != N_COMPONENTS] <- NA_real_
  result$RWDID_continuous_mean_15 <- result$RWDID_continuous_total_15 / N_COMPONENTS
  result
}

train_scores <- score_dataset(training_data, "training", frozen_spec, ID_COL)
# The frozen scoring specification was saved before applying it to the test set.
test_scores <- score_dataset(test_data, "test", frozen_spec, ID_COL)
SCORE_COLS <- paste0("score_", FOOD_MAP$food)

validate_score_dataset <- function(d, split_name) {
  mat <- as.matrix(d[, SCORE_COLS, drop = FALSE])
  storage.mode(mat) <- "numeric"
  finite <- is.finite(mat)
  if (any(mat[finite] < -NUMERICAL_TOL) ||
      any(mat[finite] > 1 + NUMERICAL_TOL)) {
    stop("STOP: Component score escaped [0,1] in ", split_name)
  }
  total <- rowSums(mat)
  total[rowSums(finite) != N_COMPONENTS] <- NA_real_
  same_na <- isTRUE(all(is.na(total) == is.na(d$RWDID_continuous_total_15)))
  if (!same_na ||
      any(abs(total - d$RWDID_continuous_total_15) > NUMERICAL_TOL, na.rm = TRUE)) {
    stop("STOP: Continuous total is not the exact 15-component equal-weight sum in ", split_name)
  }
  if (any(d$RWDID_continuous_total_15 < -NUMERICAL_TOL |
          d$RWDID_continuous_total_15 > 15 + NUMERICAL_TOL, na.rm = TRUE)) {
    stop("STOP: Total score escaped [0,15] in ", split_name)
  }
  if (any(is.na(d$RWDID_continuous_total_15) !=
          (d$n_components_observed != N_COMPONENTS))) {
    stop("STOP: Missing-total rule failed in ", split_name)
  }
  invisible(TRUE)
}
validate_score_dataset(train_scores, "training")
validate_score_dataset(test_scores, "test")

utils::write.csv(
  train_scores,
  file.path(OUT_DIR, "RWDID_continuous_component_scores_training.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
utils::write.csv(
  test_scores,
  file.path(OUT_DIR, "RWDID_continuous_component_scores_test.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
utils::write.csv(
  train_scores[, c("participant_id", "split", "n_components_observed",
                   "RWDID_continuous_total_15", "RWDID_continuous_mean_15")],
  file.path(OUT_DIR, "RWDID_continuous_total_score_training.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)
utils::write.csv(
  test_scores[, c("participant_id", "split", "n_components_observed",
                  "RWDID_continuous_total_15", "RWDID_continuous_mean_15")],
  file.path(OUT_DIR, "RWDID_continuous_total_score_test.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# -----------------------------------------------------------------------------
# 6. Reconcile against the frozen binary score and summarize agreement
# -----------------------------------------------------------------------------

binary_training <- utils::read.csv(
  RCS_BINARY_TRAIN_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
required_binary_cols <- c(ID_COL, SCORE_COLS, "RWDID_total_score_complete15")
if (!all(required_binary_cols %in% names(binary_training))) {
  stop("STOP: binary training score file lacks: ",
       paste(setdiff(required_binary_cols, names(binary_training)), collapse = ", "))
}
binary_ids <- trimws(as.character(binary_training[[ID_COL]]))
if (anyNA(binary_ids) || any(binary_ids == "") || anyDuplicated(binary_ids)) {
  stop("STOP: Invalid or duplicated IDs in the binary training score file.")
}
binary_match <- match(train_scores$participant_id, binary_ids)
if (anyNA(binary_match) || nrow(binary_training) != nrow(train_scores)) {
  stop("STOP: Binary and continuous training participant IDs do not match exactly.")
}

score_binary_component <- function(x, direction, value1, value2, food_name) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  out <- rep(NA_integer_, length(x))
  observed <- !is.na(x)
  if (any(!is.finite(x[observed])) || any(x[observed] < -1e-8)) {
    stop("STOP: Invalid raw intake while reconstructing binary score for ", food_name)
  }
  if (!any(observed)) return(out)
  scale_value <- max(1, abs(value1), ifelse(is.finite(value2), abs(value2), 0))
  tol <- max(1e-8, scale_value * 1e-8)
  if (direction == "middle") {
    out[observed] <- as.integer(x[observed] >= value1 - tol & x[observed] <= value2 + tol)
  } else if (direction == "low") {
    out[observed] <- as.integer(x[observed] <= value1 + tol)
  } else if (direction == "high") {
    out[observed] <- as.integer(x[observed] >= value1 - tol)
  } else if (direction == "point") {
    out[observed] <- as.integer(abs(x[observed] - value1) <= tol)
  } else {
    stop("STOP: Unsupported binary direction during reconstruction for ", food_name)
  }
  out
}

binary_reconciliation <- vector("list", N_COMPONENTS)
component_agreement <- vector("list", N_COMPONENTS)
for (i in seq_len(N_COMPONENTS)) {
  f <- FOOD_MAP$food[i]
  bcol <- FOOD_MAP$score_col[i]
  imported <- suppressWarnings(as.numeric(binary_training[[bcol]][binary_match]))
  recomputed <- score_binary_component(
    training_data[[f]], rcs_ranges$range_direction[i],
    rcs_ranges$value1[i], rcs_ranges$value2[i], f
  )
  same <- (is.na(imported) & is.na(recomputed)) |
    (!is.na(imported) & !is.na(recomputed) & imported == recomputed)
  binary_reconciliation[[i]] <- data.frame(
    food = f,
    imported_score_1_n = sum(imported == 1, na.rm = TRUE),
    recomputed_score_1_n = sum(recomputed == 1, na.rm = TRUE),
    mismatch_n = sum(!same),
    exact_match_to_binary_score = all(same),
    stringsAsFactors = FALSE
  )
  cont <- train_scores[[bcol]]
  ok <- is.finite(imported) & is.finite(cont)
  component_agreement[[i]] <- data.frame(
    food = f,
    n_matched = sum(ok),
    binary_mean = if (any(ok)) mean(imported[ok]) else NA_real_,
    continuous_mean = if (any(ok)) mean(cont[ok]) else NA_real_,
    pearson = if (sum(ok) > 1L && stats::sd(imported[ok]) > 0 && stats::sd(cont[ok]) > 0)
      stats::cor(imported[ok], cont[ok], method = "pearson") else NA_real_,
    spearman = if (sum(ok) > 1L && stats::sd(imported[ok]) > 0 && stats::sd(cont[ok]) > 0)
      stats::cor(imported[ok], cont[ok], method = "spearman") else NA_real_,
    stringsAsFactors = FALSE
  )
}
binary_reconciliation <- do.call(rbind, binary_reconciliation)
component_agreement <- do.call(rbind, component_agreement)
if (nrow(binary_reconciliation) != N_COMPONENTS ||
    !all(binary_reconciliation$exact_match_to_binary_score)) {
  print(binary_reconciliation)
  stop("STOP: binary scores cannot be reproduced from the imported interval rules.")
}

binary_total_imported <- suppressWarnings(
  as.numeric(binary_training$RWDID_total_score_complete15[binary_match])
)
binary_matrix <- as.matrix(binary_training[binary_match, SCORE_COLS, drop = FALSE])
storage.mode(binary_matrix) <- "numeric"
binary_total_recomputed <- rowSums(binary_matrix)
binary_total_recomputed[rowSums(!is.na(binary_matrix)) != N_COMPONENTS] <- NA_real_
same_binary_total <- (is.na(binary_total_imported) & is.na(binary_total_recomputed)) |
  (!is.na(binary_total_imported) & !is.na(binary_total_recomputed) &
     binary_total_imported == binary_total_recomputed)
if (!all(same_binary_total)) {
  stop("STOP: Imported binary total does not reconcile with its 15 score columns.")
}

continuous_total <- train_scores$RWDID_continuous_total_15
corr_ok <- is.finite(binary_total_imported) & is.finite(continuous_total)
if (sum(corr_ok) < 2L) stop("STOP: Insufficient complete binary/continuous matched scores.")
correlation_audit <- data.frame(
  n_matched = sum(corr_ok),
  binary_mean = mean(binary_total_imported[corr_ok]),
  continuous_mean = mean(continuous_total[corr_ok]),
  pearson = stats::cor(binary_total_imported[corr_ok], continuous_total[corr_ok],
                       method = "pearson"),
  spearman = stats::cor(binary_total_imported[corr_ok], continuous_total[corr_ok],
                        method = "spearman"),
  interpretation = "Training-only descriptive agreement; binary and continuous scores are not interchangeable",
  stringsAsFactors = FALSE
)

utils::write.csv(
  binary_reconciliation,
  file.path(AUDIT_DIR, "binary_score_reconstruction_QC.csv"),
  row.names = FALSE
)
utils::write.csv(
  correlation_audit,
  file.path(AUDIT_DIR, "binary_vs_continuous_training_score_correlation.csv"),
  row.names = FALSE
)
utils::write.csv(
  component_agreement,
  file.path(AUDIT_DIR, "binary_vs_continuous_component_agreement.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7. Component/total distributions and support-boundary audits
# -----------------------------------------------------------------------------

component_distribution <- list()
boundary_audit <- list()
for (split_name in c("training", "test")) {
  d <- if (split_name == "training") train_scores else test_scores
  for (i in seq_len(N_COMPONENTS)) {
    f <- FOOD_MAP$food[i]
    s <- d[[paste0("score_", f)]]
    x <- d[[paste0("intake_", f)]]
    good <- is.finite(s)
    n_valid <- sum(good)
    component_distribution[[length(component_distribution) + 1L]] <- data.frame(
      split = split_name, component_order = i, food = f,
      display_name = FOOD_MAP$display_name[i],
      score_type = parameters$score_type[i],
      n = length(s), n_valid = n_valid, n_missing = sum(!good),
      mean = if (n_valid) mean(s[good]) else NA_real_,
      sd = if (n_valid > 1L) stats::sd(s[good]) else NA_real_,
      median = if (n_valid) stats::median(s[good]) else NA_real_,
      q25 = if (n_valid) unname(stats::quantile(s[good], .25)) else NA_real_,
      q75 = if (n_valid) unname(stats::quantile(s[good], .75)) else NA_real_,
      n_zero = sum(s[good] <= NUMERICAL_TOL),
      n_between = sum(s[good] > NUMERICAL_TOL & s[good] < 1 - NUMERICAL_TOL),
      n_one = sum(s[good] >= 1 - NUMERICAL_TOL),
      pct_zero = if (n_valid) 100 * mean(s[good] <= NUMERICAL_TOL) else NA_real_,
      pct_between = if (n_valid) 100 * mean(s[good] > NUMERICAL_TOL & s[good] < 1 - NUMERICAL_TOL) else NA_real_,
      pct_one = if (n_valid) 100 * mean(s[good] >= 1 - NUMERICAL_TOL) else NA_real_,
      stringsAsFactors = FALSE
    )
    p <- parameters[i, , drop = FALSE]
    boundary_audit[[length(boundary_audit) + 1L]] <- data.frame(
      split = split_name, food = f,
      support_p1_training = p$support_p1,
      support_p99_training = p$support_p99,
      n_below_training_p1 = sum(is.finite(x) & x < p$support_p1),
      n_above_training_p99 = sum(is.finite(x) & x > p$support_p99),
      n_missing_intake = sum(!is.finite(x)),
      n_continuous_score_1 = sum(s >= 1 - NUMERICAL_TOL, na.rm = TRUE),
      support_source = "Frozen training cohort observed intake P1-P99",
      stringsAsFactors = FALSE
    )
  }
}
component_distribution <- do.call(rbind, component_distribution)
boundary_audit <- do.call(rbind, boundary_audit)
if (any(component_distribution$n_zero + component_distribution$n_between +
        component_distribution$n_one != component_distribution$n_valid)) {
  stop("STOP: Continuous component score distribution counts do not reconcile.")
}

summarise_total <- function(d) {
  x <- d$RWDID_continuous_total_15
  good <- x[is.finite(x)]
  data.frame(
    split = unique(d$split), n = nrow(d), n_complete = length(good),
    n_missing_total = sum(!is.finite(x)),
    mean = if (length(good)) mean(good) else NA_real_,
    sd = if (length(good) > 1L) stats::sd(good) else NA_real_,
    median = if (length(good)) stats::median(good) else NA_real_,
    q25 = if (length(good)) unname(stats::quantile(good, .25)) else NA_real_,
    q75 = if (length(good)) unname(stats::quantile(good, .75)) else NA_real_,
    min = if (length(good)) min(good) else NA_real_,
    max = if (length(good)) max(good) else NA_real_,
    stringsAsFactors = FALSE
  )
}
total_distribution <- rbind(
  summarise_total(train_scores), summarise_total(test_scores)
)

utils::write.csv(
  component_distribution,
  file.path(OUT_DIR, "Table_Sx_continuous_component_score_distribution_15.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  total_distribution,
  file.path(OUT_DIR, "Table_Sx_continuous_total_score_distribution.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  boundary_audit,
  file.path(AUDIT_DIR, "training_test_P1_P99_boundary_audit.csv"),
  row.names = FALSE, na = ""
)

# -----------------------------------------------------------------------------
# 8. Frozen scoring curves and journal-ready figures
# -----------------------------------------------------------------------------

curve_list <- vector("list", N_COMPONENTS)
heatmap_list <- vector("list", N_COMPONENTS)
for (i in seq_len(N_COMPONENTS)) {
  p <- parameters[i, , drop = FALSE]
  xs <- seq(p$support_p1, p$support_p99, length.out = GRID_N)
  candidate_boundaries <- c(p$value1, p$value2)
  candidate_boundaries <- candidate_boundaries[
    is.finite(candidate_boundaries) &
      candidate_boundaries >= p$support_p1 &
      candidate_boundaries <= p$support_p99
  ]
  xs <- sort(unique(c(xs, candidate_boundaries)))
  sc <- score_proportional_v531(xs, p)
  if (any(!is.finite(sc)) || any(sc < -NUMERICAL_TOL | sc > 1 + NUMERICAL_TOL)) {
    stop("STOP: Invalid scoring function grid for ", p$food)
  }
  curve_list[[i]] <- data.frame(
    component_order = i, food = p$food, display_name = p$display_name,
    range_direction = p$range_direction, score_type = p$score_type,
    intake = xs, score = sc,
    value1 = p$value1, value2 = p$value2,
    support_p1 = p$support_p1, support_p99 = p$support_p99,
    range_source = p$range_source,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    curve_list[[i]],
    file.path(CURVE_DIR, paste0(p$food, "_proportional_function.csv")),
    row.names = FALSE
  )
  frac <- seq(0, 1, length.out = GRID_N)
  qx <- p$support_p1 + frac * (p$support_p99 - p$support_p1)
  heatmap_list[[i]] <- data.frame(
    component_order = i, food = p$food, display_name = p$display_name,
    support_fraction = frac, intake = qx,
    score = score_proportional_v531(qx, p),
    panel_y = N_COMPONENTS + 1L - i,
    stringsAsFactors = FALSE
  )
}
scoring_curves <- do.call(rbind, curve_list)
heatmap_grid <- do.call(rbind, heatmap_list)
utils::write.csv(
  heatmap_grid,
  file.path(PLOT_DIR, "Figure_Sx_proportional_heatmap_training_grid.csv"),
  row.names = FALSE
)

pal <- cols4all::c4a("friendly13")
if (length(pal) < 8L) stop("STOP: friendly13 palette is unavailable or incomplete.")
curve_colours <- c(
  adequacy = unname(pal[3]),
  moderation = unname(pal[4]),
  optimum = unname(pal[5]),
  point = unname(pal[7])
)
scoring_curves$display_name <- factor(
  scoring_curves$display_name, levels = FOOD_MAP$display_name
)
scoring_curves$score_type <- factor(
  scoring_curves$score_type,
  levels = c("adequacy", "moderation", "optimum", "point")
)
p_func <- ggplot2::ggplot(
  scoring_curves,
  ggplot2::aes(x = intake, y = score, colour = score_type, group = food)
) +
  ggplot2::geom_line(linewidth = 0.9, lineend = "round") +
  ggplot2::facet_wrap(~display_name, scales = "free_x", ncol = 3) +
  ggplot2::scale_colour_manual(
    values = curve_colours,
    breaks = c("adequacy", "moderation", "optimum", "point"),
    labels = c("Adequacy", "Moderation", "Bounded optimum", "Point target")
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  ggplot2::labs(
    title = "anchored continuous proportional component scores",
    subtitle = "Full-score targets imported unchanged from the final intervals",
    x = "Food intake (servings/day)", y = "Continuous component score (0-1)",
    colour = "Scoring function",
    caption = paste(
      "Outside the full-score target, scores taper linearly to zero at",
      "the training-only P1/P99 support boundary. Sweetened beverages receive",
      "a full score from 0 to 2.08 servings/day, then taper proportionally;",
      "this is a scoring rule, not a causal or clinical recommendation."
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10.5),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    axis.title = ggplot2::element_text(size = 9, face = "bold"),
    axis.text = ggplot2::element_text(size = 8, colour = "black"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(size = 8, hjust = 0)
  )
ggplot2::ggsave(
  file.path(PLOT_DIR, "Figure_Sx_continuous_proportional_scoring_functions.pdf"),
  p_func, width = 13, height = 16, units = "in", device = "pdf"
)

# PiYG 15-row heatmap; same physical size as the binary scoring scheme figure.
plot_heatmap <- ggplot2::ggplot(
  heatmap_grid,
  ggplot2::aes(x = 100 * support_fraction, y = panel_y, fill = score)
) +
  ggplot2::geom_tile(width = 100 / (GRID_N - 1L), height = 0.84) +
  ggplot2::scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "PiYG"), limits = c(0, 1),
    breaks = c(0, .25, .5, .75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    name = "Component score",
    guide = ggplot2::guide_colourbar(
      direction = "horizontal", title.position = "top",
      barwidth = grid::unit(3.0, "in"),
      barheight = grid::unit(.14, "in")
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 50, 100),
    labels = c("Low", "Medium", "High"),
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = 15:1, labels = FOOD_MAP$display_name,
    limits = c(.5, 15.5), expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::labs(
    tag = "c", title = "anchored continuous proportional scoring",
    x = "Relative intake along each food's training P1-P99 range", y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 18),
    plot.tag.position = c(.012, .99),
    plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = .5),
    axis.text.y = ggplot2::element_text(face = "bold", size = 15, colour = "black"),
    axis.text.x = ggplot2::element_text(face = "bold", size = 15, colour = "black"),
    axis.title.x = ggplot2::element_text(size = 15, margin = ggplot2::margin(t = 7)),
    axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8.5),
    legend.text = ggplot2::element_text(size = 8),
    plot.caption = ggplot2::element_text(size = 6.9, hjust = .5,
                                         margin = ggplot2::margin(t = 6)),
    plot.margin = ggplot2::margin(18, 12, 10, 10)
  )
ggplot2::ggsave(
  file.path(PLOT_DIR, "Figure_Sx_continuous_proportional_score_heatmap_PiYG.pdf"),
  plot_heatmap, width = 599 / 72, height = 662 / 72, units = "in",
  device = "pdf", limitsize = FALSE
)

# Training-only workflow figure (test cohort is only scored after freezing).
workflow_labels <- c(
  "Frozen Step_03_1_0 training cohort\nN=76,578; depression events=2,122",
  "Import final intervals and scoring directions\n15 components; no interval re-derivation",
  "Estimate training-only exposure support\nP1/P99 for each food component",
  "Define proportional 0-1 score functions\nAdequacy | Moderation | Bounded optimum | Point target",
  "Freeze the 15 functions and parameters\nEqual-weight total score: 0-15",
  "Apply the frozen functions\nTraining and held-out test cohorts",
  "QC and export\nBinary reconciliation | Distributions | Figures"
)
workflow <- data.frame(
  step = 1:7, y = 7:1, label = workflow_labels,
  fill_group = c("training", "interval", "support", "score", "freeze", "apply", "qc"),
  stringsAsFactors = FALSE
)
workflow$xmin <- .08; workflow$xmax <- .92
workflow$ymin <- workflow$y - .34; workflow$ymax <- workflow$y + .34
workflow_fill <- c(
  training = "#F3F7F8", interval = "#E8F5F1", support = "#E9F3FA",
  score = "#F3EDF6", freeze = "#E8F5F1", apply = "#F4F5F7", qc = "#E9F3FA"
)
p_flow <- ggplot2::ggplot(workflow) +
  ggplot2::geom_segment(
    data = workflow[workflow$step < 7L, , drop = FALSE],
    ggplot2::aes(x = .5, y = y - .34, xend = .5, yend = y - .63),
    inherit.aes = FALSE, linewidth = .5, colour = pal[3],
    arrow = grid::arrow(length = grid::unit(.13, "in"), type = "closed")
  ) +
  ggplot2::geom_rect(
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                 fill = fill_group),
    colour = pal[3], linewidth = .6
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = .5, y = y, label = label),
    colour = "#202A32", size = 3.35, lineheight = 1.16
  ) +
  ggplot2::scale_fill_manual(values = workflow_fill, guide = "none") +
  ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ggplot2::scale_y_continuous(limits = c(.4, 7.55), expand = c(0, 0)) +
  ggplot2::labs(
    title = "Development of the anchored continuous proportional RWDID",
    subtitle = "Training-derived parameters are frozen before held-out test scoring"
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = .5),
    plot.subtitle = ggplot2::element_text(size = 10, hjust = .5),
    plot.margin = ggplot2::margin(16, 16, 13, 16)
  )
ggplot2::ggsave(
  file.path(PLOT_DIR, "Figure_Sx_continuous_proportional_scoring_workflow.pdf"),
  p_flow, width = 8.5, height = 11, units = "in", device = "pdf"
)

# -----------------------------------------------------------------------------
# 9. Save loadable objects and final QC
# -----------------------------------------------------------------------------

saveRDS(
  list(
    train_scores = train_scores,
    test_scores = test_scores,
    parameters = parameters,
    frozen_spec = frozen_spec,
    total_distribution = total_distribution,
    binary_continuous_correlation = correlation_audit
  ),
  file.path(OUT_DIR, "RWDID_continuous_scores_and_spec.rds"),
  version = 2
)
save(
  train_scores, test_scores, parameters, frozen_spec,
  total_distribution, correlation_audit,
  file = file.path(OUT_DIR, "RWDID_continuous_scores.RData")
)

files_required_output <- c(
  file.path(OUT_DIR, "Table_Sx_continuous_scoring_parameters_15.csv"),
  file.path(OUT_DIR, "RWDID_continuous_component_scores_training.csv"),
  file.path(OUT_DIR, "RWDID_continuous_component_scores_test.csv"),
  file.path(OUT_DIR, "RWDID_continuous_total_score_training.csv"),
  file.path(OUT_DIR, "RWDID_continuous_total_score_test.csv"),
  file.path(OUT_DIR, "RWDID_continuous_scoring_frozen_spec.rds"),
  file.path(OUT_DIR, "RWDID_continuous_scores_and_spec.rds"),
  file.path(OUT_DIR, "RWDID_continuous_scores.RData"),
  file.path(AUDIT_DIR, "sweetened_bev_continuous_function_QC.csv"),
  file.path(OUT_DIR, "Table_Sx_continuous_component_score_distribution_15.csv"),
  file.path(OUT_DIR, "Table_Sx_continuous_total_score_distribution.csv"),
  file.path(PLOT_DIR, "Figure_Sx_continuous_proportional_scoring_functions.pdf"),
  file.path(PLOT_DIR, "Figure_Sx_continuous_proportional_score_heatmap_PiYG.pdf"),
  file.path(PLOT_DIR, "Figure_Sx_continuous_proportional_scoring_workflow.pdf")
)

qc <- data.frame(
  check = c(
    "frozen_training_N", "frozen_test_N", "frozen_training_events",
    "frozen_test_events", "training_test_ID_disjoint",
    "exact_final_15_component_order", "excluded_components_absent",
    "range_and_rule_files_reconcile", "two_repaired_intervals_match",
    "sweetened_beverage_middle_interval_reconciles_to_frozen_rule",
    "sweetened_beverage_continuous_function_matches_frozen_middle_interval",
    "all_training_P1_P99_supports_finite_and_non_degenerate",
    "all_targets_intersect_training_support", "all_component_scores_in_0_1",
    "training_total_is_equal_weight_sum", "test_total_is_equal_weight_sum",
    "missing_component_implies_missing_total", "binary_scores_reconstructed_exactly",
    "binary_total_reconciles", "continuous_score_is_not_beta_weighted",
    "no_additional_HR_tolerance_applied", "test_uses_frozen_training_parameters",
    "all_expected_outputs_created"
  ),
  pass = c(
    nrow(training_data) == EXPECTED_TRAIN_N,
    nrow(test_data) == EXPECTED_TEST_N,
    sum(training_data$depression_status == 1, na.rm = TRUE) == EXPECTED_TRAIN_EVENTS,
    sum(test_data$depression_status == 1, na.rm = TRUE) == EXPECTED_TEST_EVENTS,
    length(intersect(train_scores$participant_id, test_scores$participant_id)) == 0L,
    identical(as.character(parameters$food), FOOD_MAP$food),
    !any(c("green_leafy_veg", "high_fat_dairy") %in% parameters$food),
    range_rule_match,
    repair_pass,
    sweetened_middle_pass,
    sweetened_continuous_function_pass,
    all(is.finite(parameters$support_p1)) &&
      all(is.finite(parameters$support_p99)) &&
      all(parameters$support_p99 > parameters$support_p1),
    all(parameters$target_intersects_training_support),
    all(vapply(
      list(train_scores, test_scores),
      function(d) {
        m <- as.matrix(d[, SCORE_COLS, drop = FALSE])
        all(m[is.finite(m)] >= -NUMERICAL_TOL & m[is.finite(m)] <= 1 + NUMERICAL_TOL)
      }, logical(1)
    )),
    all(is.na(train_scores$RWDID_continuous_total_15) ==
          (train_scores$n_components_observed != N_COMPONENTS)),
    all(is.na(test_scores$RWDID_continuous_total_15) ==
          (test_scores$n_components_observed != N_COMPONENTS)),
    all(is.na(train_scores$RWDID_continuous_total_15) ==
          (train_scores$n_components_observed != N_COMPONENTS)) &&
      all(is.na(test_scores$RWDID_continuous_total_15) ==
            (test_scores$n_components_observed != N_COMPONENTS)),
    all(binary_reconciliation$exact_match_to_binary_score),
    all(same_binary_total),
    isTRUE(frozen_spec$no_beta_weighting),
    isTRUE(frozen_spec$no_additional_hr_tolerance),
    isTRUE(frozen_spec$test_outcomes_used_to_define_score == FALSE),
    all(file.exists(files_required_output) &
          !is.na(file.info(files_required_output)$size) &
          file.info(files_required_output)$size > 0)
  ),
  stringsAsFactors = FALSE
)
qc$status <- ifelse(qc$pass, "PASS", "FAIL")
utils::write.csv(
  qc, file.path(AUDIT_DIR, "continuous_proportional_scoring_FINAL_QC.csv"),
  row.names = FALSE
)
if (any(qc$status != "PASS")) {
  print(qc)
  stop("STOP: continuous proportional scoring final QC failed.")
}

method_note <- c(
  "RWDID CONTINUOUS PROPORTIONAL SCORING — ANCHORED COLIZZI ADAPTATION",
  "==========================================================================",
  "This is a secondary equal-weight continuous score, separate from the primary binary RWDID and the beta-weighted score.",
  "The 15 full-score target intervals and their directions are imported unchanged from Step_04_1_1.",
  "The range CSV and binary scoring-rule CSV are required to match exactly on direction, values, score limits and score rule.",
  "High direction: proportional adequacy rise from training P1 to value1, then score 1.",
  "Low direction: score 1 through value1, then proportional moderation decline to training P99.",
  "Middle direction: score 1 within [value1,value2], with proportional ramps to zero at training P1/P99.",
  "Point direction: score 1 at value1, with a triangular proportional taper to zero at training P1/P99.",
  "For sweetened beverages, the frozen rule is direction=middle, value1=0 and value2=2.08 (score_rule: 0 <= x <= 2.08). The continuous score is 1 throughout this interval and declines linearly above 2.08 to 0 at training P99 when P99 exceeds 2.08.",
  "Training P1/P99 are calculated from observed food intakes in the frozen training cohort only (type=7), separately for each component.",
  "All test-set scores use the already frozen training P1/P99 and interval functions; no test parameter is estimated.",
  "Intakes outside training P1-P99 are clamped to the nearest training support boundary before scoring.",
  "Each component score is in [0,1]; the 15 components are equally weighted; total range is [0,15].",
  "No imputation is performed. If any component intake is missing, the corresponding total score is NA.",
  "No additional HR-based 5%, 10% or 20% tolerance is used to define continuous score boundaries.",
  "The 20% near-minimum-risk rule is embedded only where it was used in the frozen interval derivation; it is not re-applied here.",
  "This proportional continuous adaptation is not the original Colizzi/EAT-Lancet index and is not a nutritional prescription.",
  "The refined-grains one-sided rule is treated as an adequacy function; it remains exploratory if its RCS evidence is not confirmatory.",
  "The score function is an operational continuous representation of the frozen binary intervals; continuous and binary scores are not interchangeable.",
  "Generated outputs include component-level training/test scores, total scores, frozen parameters, binary-score reconciliation, distributions and figures.",
  "No Cox association model is fitted by this script. Association analysis remains a separate downstream step.",
  "No score results are claimed until the script is run against the frozen project RDS and Step_04_1_1 outputs."
)
writeLines(
  method_note,
  con = file.path(OUT_DIR, "README_COLIZZI_ADAPTED_CONTINUOUS_SCORING.txt"),
  useBytes = TRUE
)

cat("\n============================================================\n")
cat("STEP_04_1_2 | ANCHORED CONTINUOUS PROPORTIONAL SCORE FIX3 COMPLETE\n")
cat("Components: 15 | Equal weight | Component range: 0-1 | Total: 0-15\n")
cat("Training N: ", nrow(train_scores), " | Test N: ", nrow(test_scores), "\n", sep = "")
cat("Training complete scores: ", sum(is.finite(train_scores$RWDID_continuous_total_15)), "\n", sep = "")
cat("Test complete scores: ", sum(is.finite(test_scores$RWDID_continuous_total_15)), "\n", sep = "")
cat("Training total mean: ", round(total_distribution$mean[1], 4),
    " | Test total mean: ", round(total_distribution$mean[2], 4), "\n", sep = "")
cat("Binary score reconstruction: ",
    sum(binary_reconciliation$exact_match_to_binary_score), "/15 PASS\n", sep = "")
cat("Final QC: ", sum(qc$status == "PASS"), "/", nrow(qc), " PASS\n", sep = "")
cat("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 09/11: Step_04_1_2_RWDID_C
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 10/11: Step_04_2_1_BETA_WEIGHTED_TRAINING
# Source: Step_04_2_1_beta_weighted_diet_score_training_V1_3_RCS_DATADIST_FIX(20261118-041346).R
# Provenance: Library source dated 2026-11-18; current Step_04_2_1 V1.3 RCS/DATADIST fix
# SHA256 of included block: 665f8763025cf8e2ba6dfc2789ef577a161599c75b030045deee7905cb6abc3a
########################################################################################

# =============================================================================
# Step_04_2_1_beta_weighted_diet_score_training_FINAL.R
#
# FINAL beta-weight estimation and weighted-score association.
#
# Frozen inputs:
#   - Exact Step_03_1_0 70% training cohort (seed 20261103).
#   - Final intervals, scoring rules and participant-level binary scores
#     exported by Step_04_1_1_4.
#   - Original Model 3 covariates and mutually adjusted 15-component Cox model.
#
# Separation of responsibilities:
#   - Step_04_1_1 performs RCS, interval derivation, binary scoring and Table Sx.
#   - This Step_04_2_1 script imports those frozen outputs without modification.
#   - No RCS model is fitted and no food-group interval is re-derived here.
#   - All 15 imported binary scores are independently reconciled with raw intake.
#   - Beta signs are estimated without sign constraints.
#   - Only the training cohort is used; the held-out 30% test cohort is excluded.
#
# Run from the UK Biobank project root.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(broom)
  library(forestplot)
  library(cols4all)
  library(grid)
})

options(stringsAsFactors = FALSE, scipen = 999)

# -----------------------------------------------------------------------------
# 0. Frozen inputs, output paths, and expected split metadata
# -----------------------------------------------------------------------------

SPLIT_DIR <- "./results/Step_03/ukb_70_30_split"
MAIN_SPLIT_DIR <- file.path(SPLIT_DIR, "main")

SPLIT_MANIFEST_FILE <- file.path(
  SPLIT_DIR, "ukb_70_30_split_manifest.csv"
)
TRAIN_RDS <- file.path(
  MAIN_SPLIT_DIR, "ukb_depression_training_cohort.rds"
)

# This is permanently frozen to the RCS20 workflow.
# No environment-variable switch, file search, wildcard, or fallback path is used.
V53_VARIANT <- "RCS20"

# Exact files written by Step_04_1_1:
# ./results/Step_04/Step_04_1_1_4
RCS_RANGE_FILE <- paste0(
  "./results/Step_04/Step_04_1_1_4/",
  "food_rcs_range_all_final.csv"
)
RCS_SCORE_FILE <- paste0(
  "./results/Step_04/Step_04_1_1_4/",
  "RWDID_binary_component_scores_training.csv"
)
RCS_RULE_FILE <- paste0(
  "./results/Step_04/Step_04_1_1_4/",
  "RWDID_binary_scoring_rules_15_components.csv"
)
RCS_DIAG_FILE <- paste0(
  "./results/Step_04/Step_04_1_1_4/",
  "audit/RCS_model_summary_and_range_diagnostics.csv"
)
REPAIR_AUDIT_FILE <- paste0(
  "./results/Step_04/Step_04_1_1_4/",
  "audit/two_component_interval_derivation.csv"
)

# New, separate output directory. Existing V5.1/V5.3 outputs remain untouched.
OUT_DIR <- file.path(
  "./results/Step_04",
  paste0("Step_04_2_1")
)
PLOT_DIR <- file.path(OUT_DIR, "plots")
AUDIT_DIR <- file.path(OUT_DIR, "audit")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPECTED_SPLIT_SEED <- 20261103L
EXPECTED_TRAIN_PROP <- 0.70
EXPECTED_TRAIN_N <- 76578L
EXPECTED_TEST_N <- 32819L
EXPECTED_TOTAL_N <- EXPECTED_TRAIN_N + EXPECTED_TEST_N

# Original Model 3 adjustment set retained from the uploaded script and
# Step_04_1_1 final RCS workflow.
MODEL3_COVARIATES <- c(
  "age_food",
  "sex",
  "ethnicity",
  "education_attainment",
  "twnsend",
  "career",
  "smoking_status",
  "drink_status",
  "bmi_group",
  "physical_activity_proup",
  "sedentary_group",
  "household_income",
  "energy_intake",
  "central_obesity",
  "triglycerides",
  "hdl_cholesterol",
  "t2diabetes_status_final",
  "hypertension_status_final",
  "dyslipidemia_status_final",
  "cvd_status_final",
  "all_cancer_status_final"
)

# The exact 15-component order from the final Step_04_1_1 RCS script.
# The legacy Green leafy vegetables and High-fat dairy components are removed;
# Fruit juices and Water are included.
FOOD_MAP <- tibble(
  group = c(
    "Citrus fruits",
    "Red/orange vegetables",
    "Other vegetables",
    "Potatoes",
    "Refined grains",
    "Eggs",
    "Low fat dairy",
    "Red meats",
    "Processed meats",
    "Fish and seafood",
    "Coffee",
    "Fruit juices",
    "Sweetened beverages",
    "Water",
    "Wine"
  ),
  food = c(
    "citrus_fruits",
    "red_orange_veg",
    "other_veg",
    "potatoes",
    "refined_grains",
    "eggs",
    "low_fat_dairy",
    "red_meats",
    "processed_meats",
    "fish_seafood",
    "coffee",
    "fruit_juices",
    "sweetened_bev",
    "water",
    "wine"
  ),
  display_name = c(
    "Citrus fruits",
    "Red and orange vegetables",
    "Other vegetables",
    "Potatoes",
    "Refined grains",
    "Eggs",
    "Low-fat dairy",
    "Red meats",
    "Processed meats",
    "Fish and seafood",
    "Coffee",
    "Fruit juices",
    "Sweetened beverages",
    "Water",
    "Wine"
  ),
  score_col = paste0(
    "score_",
    c(
      "citrus_fruits", "red_orange_veg", "other_veg", "potatoes",
      "refined_grains", "eggs", "low_fat_dairy", "red_meats",
      "processed_meats", "fish_seafood", "coffee", "fruit_juices",
      "sweetened_bev", "water", "wine"
    )
  ),
  # Keep the original *_cat naming convention for the beta model terms.
  model_var = paste0(
    c(
      "citrus_fruits", "red_orange_veg", "other_veg", "potatoes",
      "refined_grains", "eggs", "low_fat_dairy", "red_meats",
      "processed_meats", "fish_seafood", "coffee", "fruit_juices",
      "sweetened_bev", "water", "wine"
    ),
    "_cat"
  )
)

ALL_FOOD_VARS <- FOOD_MAP$model_var
ALL_SCORE_COLS <- FOOD_MAP$score_col

# V5.3 uses all 15 audited components. The *_PRIMARY* aliases are retained only
# for compatibility with the existing downstream Step_04_2_2 file structure;
# they do not represent an eligibility filter in this V5.1 workflow.
FOOD_VARS <- ALL_FOOD_VARS
SCORE_COLS <- ALL_SCORE_COLS
ACTIVE_FOOD_MAP <- FOOD_MAP
ACTIVE_SCORE_COLS <- ALL_SCORE_COLS
N_COMPONENTS <- nrow(FOOD_MAP)
N_PRIMARY_COMPONENTS <- N_COMPONENTS
COVARS <- MODEL3_COVARIATES

# Only these two components are direction-repaired in this sensitivity run.
REPAIR_FOODS <- c("potatoes", "red_meats")
EXPECTED_REPAIR_INTERVALS <- tibble(
  food = c("potatoes", "red_meats"),
  expected_lower = c(1.00, 1.00),
  expected_upper = c(1.50, 1.75)
)

if (nrow(FOOD_MAP) != 15L ||
    anyDuplicated(FOOD_MAP$food) > 0L ||
    anyDuplicated(FOOD_MAP$model_var) > 0L ||
    anyDuplicated(FOOD_MAP$score_col) > 0L) {
  stop("STOP: Internal final 15-component map is invalid.")
}
if (any(c("green_leafy_veg", "high_fat_dairy") %in% FOOD_MAP$food)) {
  stop("STOP: Excluded legacy components entered the final component map.")
}
if (!all(c("fruit_juices", "water") %in% FOOD_MAP$food)) {
  stop("STOP: Fruit juices and Water must be included in the final 15.")
}

# -----------------------------------------------------------------------------
# 1. Utility functions
# -----------------------------------------------------------------------------

stop_now <- function(...) {
  stop(paste0(...), call. = FALSE)
}

fmt_p <- function(p) {
  if (length(p) != 1L || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  formatC(p, format = "f", digits = 3)
}

fmt_num <- function(x, digits = 3) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

fmt_hr_ci <- function(hr, lo, hi) {
  if (any(!is.finite(c(hr, lo, hi)))) return(NA_character_)
  paste0(
    formatC(hr, format = "f", digits = 2),
    " (",
    formatC(lo, format = "f", digits = 2),
    ", ",
    formatC(hi, format = "f", digits = 2),
    ")"
  )
}

write_csv_utf8 <- function(x, file) {
  utils::write.csv(
    x, file = file, row.names = FALSE, na = "",
    fileEncoding = "UTF-8"
  )
}

# Exact independent implementation of the imported V5.1 scoring logic.
# In particular, V5.1 high-direction intervals have an open upper tail:
# score=1 for intake >= value1, with no finite value2 cap.
score_one_component <- function(x, direction, value1, value2, food_name) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  score <- rep(NA_integer_, length(x))
  observed <- !is.na(x)

  if (any(!is.finite(x[observed]))) {
    stop_now(
      "Non-finite non-missing exposure values found for ", food_name, "."
    )
  }
  if (any(x[observed] < -1e-8)) {
    stop_now("Negative intake values found for ", food_name, ".")
  }
  if (!any(observed)) return(score)

  scale_value <- max(
    1,
    abs(value1),
    ifelse(is.finite(value2), abs(value2), 0)
  )
  tol <- max(1e-8, scale_value * 1e-8)

  if (direction == "middle") {
    if (!is.finite(value1) || !is.finite(value2) || value2 < value1) {
      stop_now("Invalid middle interval for ", food_name, ".")
    }
    score[observed] <- as.integer(
      x[observed] >= value1 - tol &
        x[observed] <= value2 + tol
    )
  } else if (direction == "low") {
    if (!is.finite(value1)) {
      stop_now("Invalid low interval for ", food_name, ".")
    }
    score[observed] <- as.integer(x[observed] <= value1 + tol)
  } else if (direction == "high") {
    if (!is.finite(value1)) {
      stop_now("Invalid high interval for ", food_name, ".")
    }
    # V5 explicitly opens the upper tail for high-intake patterns.
    score[observed] <- as.integer(x[observed] >= value1 - tol)
  } else if (direction == "point") {
    if (!is.finite(value1)) {
      stop_now("Invalid point interval for ", food_name, ".")
    }
    score[observed] <- as.integer(
      abs(x[observed] - value1) <= tol
    )
  } else {
    stop_now(
      "Unsupported V5 interval direction for ", food_name,
      ": ", direction
    )
  }

  score
}

make_survival_formula <- function(exposure, covars = COVARS) {
  stats::as.formula(
    paste0(
      "survival::Surv(depression_time, depression_status) ~ ",
      paste(c(exposure, covars), collapse = " + ")
    )
  )
}

fit_cox_model <- function(formula, data, model_name) {
  fit <- tryCatch(
    survival::coxph(
      formula,
      data = data,
      ties = "efron",
      x = TRUE,
      y = TRUE,
      model = TRUE,
      na.action = na.fail
    ),
    error = function(e) {
      stop_now(
        "Cox model failed (", model_name, "): ", conditionMessage(e)
      )
    }
  )

  if (is.null(fit$coefficients) ||
      any(!is.finite(fit$coefficients)) ||
      anyNA(fit$coefficients)) {
    stop_now(
      "Cox model contains non-finite or non-estimable coefficients: ",
      model_name, "."
    )
  }
  if (!is.finite(as.numeric(tail(fit$loglik, 1)))) {
    stop_now("Cox model log-likelihood is non-finite: ", model_name, ".")
  }

  fit
}

extract_tidy_term <- function(fit, term) {
  # FIX_TERM_V1_2: `td` has a column also named `term`.
  # In dplyr::filter(), `.data$term == term` would compare that column with
  # itself and select all coefficients (27 rows in the reported Cox model).
  # `.env$term` explicitly refers to this function's requested coefficient.
  if (!is.character(term) || length(term) != 1L ||
      is.na(term) || !nzchar(trimws(term))) {
    stop_now("extract_tidy_term requires one non-empty coefficient name.")
  }

  td <- broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE,
    conf.level = 0.95
  )
  if (!"term" %in% names(td)) {
    stop_now("The Cox tidy result has no 'term' column.")
  }

  # IMPORTANT: exact term match, not prefix / substring matching.
  out <- td %>% dplyr::filter(.data$term == .env$term)
  if (nrow(out) != 1L) {
    stop_now(
      "Expected exactly one coefficient for term '", term,
      "', but found ", nrow(out), ". Available coefficients: ",
      paste(as.character(td$term), collapse = ", ")
    )
  }
  out
}

make_result_row <- function(variable, term, fit, data,
                            count_data = data, reference = FALSE) {
  n_total <- nrow(count_data)
  n_events <- sum(count_data$depression_status == 1, na.rm = TRUE)
  events_total <- paste0(
    formatC(n_events, format = "d", big.mark = ","),
    "/",
    formatC(n_total, format = "d", big.mark = ",")
  )

  if (reference) {
    return(tibble(
      variable = variable,
      Events_Total = events_total,
      hr_ci = "Ref.",
      p_value = "Ref.",
      hazard_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value_numeric = NA_real_,
      term = NA_character_
    ))
  }

  z <- extract_tidy_term(fit, term)
  tibble(
    variable = variable,
    Events_Total = events_total,
    hr_ci = fmt_hr_ci(z$estimate[1], z$conf.low[1], z$conf.high[1]),
    p_value = fmt_p(z$p.value[1]),
    hazard_ratio = as.numeric(z$estimate[1]),
    ci_low = as.numeric(z$conf.low[1]),
    ci_high = as.numeric(z$conf.high[1]),
    p_value_numeric = as.numeric(z$p.value[1]),
    term = term
  )
}

safe_cox_zph <- function(fit, model_name) {
  z <- tryCatch(
    survival::cox.zph(fit, transform = "km"),
    error = function(e) NULL
  )
  if (is.null(z)) {
    return(tibble(
      model = model_name,
      term = "NOT_ESTIMABLE",
      chisq = NA_real_,
      df = NA_real_,
      p_value = NA_real_
    ))
  }

  tab <- as.data.frame(z$table)
  tibble(
    model = model_name,
    term = rownames(tab),
    chisq = as.numeric(tab[, 1]),
    df = as.numeric(tab[, 2]),
    p_value = as.numeric(tab[, 3])
  )
}

# -----------------------------------------------------------------------------
# 2. Validate the exact frozen Step_03_1_0 70/30 split
# -----------------------------------------------------------------------------

required_inputs <- c(
  SPLIT_MANIFEST_FILE,
  TRAIN_RDS,
  RCS_RANGE_FILE,
  RCS_SCORE_FILE
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop_now(
    "Required frozen input file(s) are missing:\n",
    paste(missing_inputs, collapse = "\n"),
    "\nRun Step_03_1_0 and the final Step_04_1_1 RCS script first. ",
    "This script will not rebuild or resplit the cohort."
  )
}

split_manifest <- read.csv(
  SPLIT_MANIFEST_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(
    participant_id = "character",
    split = "character",
    id_column = "character",
    split_seed = "integer",
    training_fraction = "numeric"
  )
)

required_manifest_cols <- c(
  "participant_id", "split", "id_column",
  "split_seed", "training_fraction"
)
if (!all(required_manifest_cols %in% names(split_manifest))) {
  stop_now("Unexpected Step_03_1_0 split-manifest structure.")
}
if (anyNA(split_manifest$participant_id) ||
    any(trimws(split_manifest$participant_id) == "") ||
    anyDuplicated(split_manifest$participant_id) > 0L) {
  stop_now("Invalid or duplicated participant IDs in the split manifest.")
}
if (!all(split_manifest$split %in% c("training", "test"))) {
  stop_now("Invalid split labels in the frozen manifest.")
}
if (length(unique(split_manifest$id_column)) != 1L ||
    anyNA(split_manifest$id_column)) {
  stop_now("The manifest must contain exactly one non-missing ID column.")
}
if (length(unique(split_manifest$split_seed)) != 1L ||
    anyNA(split_manifest$split_seed) ||
    unique(split_manifest$split_seed) != EXPECTED_SPLIT_SEED) {
  stop_now("Split seed does not match the frozen Step_03_1_0 seed.")
}
if (length(unique(split_manifest$training_fraction)) != 1L ||
    anyNA(split_manifest$training_fraction) ||
    abs(unique(split_manifest$training_fraction) - EXPECTED_TRAIN_PROP) > 1e-12) {
  stop_now("Training fraction does not match the frozen 70% split.")
}

ID_COL <- unique(split_manifest$id_column)

training_ids <- sort(
  trimws(as.character(
    split_manifest$participant_id[split_manifest$split == "training"]
  )),
  method = "radix"
)
test_ids <- sort(
  trimws(as.character(
    split_manifest$participant_id[split_manifest$split == "test"]
  )),
  method = "radix"
)

if (length(training_ids) != EXPECTED_TRAIN_N ||
    length(test_ids) != EXPECTED_TEST_N ||
    nrow(split_manifest) != EXPECTED_TOTAL_N) {
  stop_now(
    "Frozen split counts do not match the expected cohort: training=",
    length(training_ids), " (expected ", EXPECTED_TRAIN_N,
    "); test=", length(test_ids), " (expected ", EXPECTED_TEST_N,
    "); total=", nrow(split_manifest), " (expected ", EXPECTED_TOTAL_N, ")."
  )
}

training_data <- readRDS(TRAIN_RDS)
if (!is.data.frame(training_data)) {
  stop_now("The frozen training RDS is not a data.frame.")
}
if (!ID_COL %in% names(training_data)) {
  stop_now("Manifest ID column is absent from the training RDS: ", ID_COL)
}

training_data[[ID_COL]] <- trimws(as.character(training_data[[ID_COL]]))
training_ids_observed <- sort(training_data[[ID_COL]], method = "radix")

if (anyNA(training_ids_observed) ||
    any(training_ids_observed == "") ||
    anyDuplicated(training_ids_observed) > 0L) {
  stop_now("Invalid or duplicated participant IDs in the training RDS.")
}
if (!identical(training_ids_observed, training_ids)) {
  stop_now(
    "Training RDS IDs do not exactly match the Step_03_1_0 manifest."
  )
}
if (length(intersect(training_ids_observed, test_ids)) != 0L) {
  stop_now("Test participants are present in the training RDS.")
}
if (nrow(training_data) != EXPECTED_TRAIN_N) {
  stop_now("Training RDS row count does not match the frozen manifest.")
}
if ("analysis_split" %in% names(training_data) &&
    any(as.character(training_data$analysis_split) != "training")) {
  stop_now("Non-training rows were found in the frozen training RDS.")
}

required_analysis_cols <- c(
  ID_COL, "depression_time", "depression_status",
  FOOD_MAP$food, COVARS
)
missing_analysis_cols <- setdiff(required_analysis_cols, names(training_data))
if (length(missing_analysis_cols) > 0L) {
  stop_now(
    "Training RDS is missing required analysis variables: ",
    paste(missing_analysis_cols, collapse = ", ")
  )
}

if (anyNA(training_data$depression_time) ||
    any(!is.finite(training_data$depression_time)) ||
    any(training_data$depression_time <= 0)) {
  stop_now("Invalid depression_time values are present in the training RDS.")
}
if (anyNA(training_data$depression_status) ||
    !all(training_data$depression_status %in% c(0, 1))) {
  stop_now("depression_status must be complete and coded 0/1.")
}

# -----------------------------------------------------------------------------
# 3. Import and strictly validate the FINAL RCS outputs.
#
# All RCS inputs below are fixed literal paths to files written by Step_04_1_1.
# file.exists() validates those exact paths only; it does not search for files.
#
# This is deliberately a beta-weight-only step:
#   - no interval derivation or modification is performed here;
#   - no RCS models are fitted here;
#   - all 15 intervals, rules and participant-level binary scores are read from
#     the completed Step_04_1_1_4 output.
# -----------------------------------------------------------------------------

# The four files below are the authoritative inputs required for beta scoring:
# final intervals, participant-level component scores, scoring rules, and RCS
# diagnostics. The two-component interval-derivation audit is supplementary
# provenance only; its absence must not block an otherwise valid beta run.
required_v531_inputs <- c(
  RCS_RANGE_FILE, RCS_SCORE_FILE, RCS_RULE_FILE, RCS_DIAG_FILE
)
missing_v531_inputs <- required_v531_inputs[!file.exists(required_v531_inputs)]
if (length(missing_v531_inputs) > 0L) {
  stop_now(
    "Required frozen Step_04_1_1 output(s) are missing:\\n",
    paste(missing_v531_inputs, collapse = "\\n"),
    "\\nRun Step_04_1_1_RCS_analysis_between_food_and_depression_FINAL.R first."
  )
}

rcs_ranges <- read.csv(
  RCS_RANGE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
required_range_cols <- c(
  "food", "value1", "value2", "range_direction", "range_source",
  "score_lower", "score_upper", "score_rule"
)
missing_range_cols <- setdiff(required_range_cols, names(rcs_ranges))
if (length(missing_range_cols) > 0L) {
  stop_now("Final range CSV is missing: ",
           paste(missing_range_cols, collapse = ", "))
}
if (nrow(rcs_ranges) != 15L ||
    anyDuplicated(as.character(rcs_ranges$food)) > 0L ||
    !setequal(as.character(rcs_ranges$food), FOOD_MAP$food)) {
  stop_now("Final Step_04_1_1 range file does not contain exactly the frozen 15 foods.")
}
rcs_ranges <- rcs_ranges[
  match(FOOD_MAP$food, as.character(rcs_ranges$food)), , drop = FALSE
]
if (!identical(as.character(rcs_ranges$food), FOOD_MAP$food)) {
  stop_now("Could not reconcile final range component order.")
}
rcs_ranges$value1 <- as.numeric(rcs_ranges$value1)
rcs_ranges$value2 <- as.numeric(rcs_ranges$value2)
rcs_ranges$score_lower <- as.numeric(rcs_ranges$score_lower)
rcs_ranges$score_upper <- as.numeric(rcs_ranges$score_upper)
if (anyNA(rcs_ranges$range_direction) ||
    any(!rcs_ranges$range_direction %in% c("middle", "low", "high", "point")) ||
    any(!is.finite(rcs_ranges$value1)) ||
    anyNA(rcs_ranges$range_source) ||
    any(trimws(rcs_ranges$range_source) == "") ||
    anyNA(rcs_ranges$score_rule) ||
    any(rcs_ranges$score_rule == "not_defined")) {
  stop_now("Final interval definitions contain missing or invalid fields.")
}
if (any(rcs_ranges$range_direction == "high" &
        !is.infinite(rcs_ranges$score_upper))) {
  stop_now("Final high-direction score intervals must retain an open upper tail.")
}
middle_rows <- rcs_ranges$range_direction == "middle"
if (any(!is.finite(rcs_ranges$value2[middle_rows])) ||
    any(rcs_ranges$value2[middle_rows] < rcs_ranges$value1[middle_rows])) {
  stop_now("One or more final middle intervals are invalid.")
}

# Import the model-derived RCS diagnostic table from the exact Step_04_1_1 path.
# This is an audit import only: Step_04_2_1 does not refit RCS or change ranges.
rcs_diagnostics <- read.csv(
  RCS_DIAG_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
required_rcs_diag_cols <- c(
  "food", "reference_HR", "p10", "p50", "p90",
  "rcs_knot_1", "rcs_knot_2", "rcs_knot_3",
  "plateau_20_low", "plateau_20_high", "plateau_20_threshold_hr",
  "p_overall", "p_non_linear", "model_n", "model_events"
)
missing_rcs_diag_cols <- setdiff(required_rcs_diag_cols, names(rcs_diagnostics))
if (length(missing_rcs_diag_cols) > 0L) {
  stop_now(
    "Step_04_1_1 RCS diagnostic CSV is missing required columns: ",
    paste(missing_rcs_diag_cols, collapse = ", ")
  )
}
if (nrow(rcs_diagnostics) != 15L ||
    anyDuplicated(as.character(rcs_diagnostics$food)) > 0L ||
    !setequal(as.character(rcs_diagnostics$food), FOOD_MAP$food)) {
  stop_now("Step_04_1_1 RCS diagnostic CSV does not contain exactly the frozen 15 foods.")
}
rcs_diagnostics <- rcs_diagnostics[
  match(FOOD_MAP$food, as.character(rcs_diagnostics$food)), , drop = FALSE
]
if (!identical(as.character(rcs_diagnostics$food), FOOD_MAP$food)) {
  stop_now("Could not reconcile Step_04_1_1 RCS diagnostic component order.")
}
diag_numeric_cols <- setdiff(required_rcs_diag_cols, "food")
for (nm in diag_numeric_cols) {
  rcs_diagnostics[[nm]] <- suppressWarnings(as.numeric(rcs_diagnostics[[nm]]))
}
if (any(!is.finite(rcs_diagnostics$reference_HR)) ||
    any(abs(rcs_diagnostics$reference_HR - 1) > 1e-6) ||
    any(!is.finite(as.matrix(rcs_diagnostics[, c(
      "p10", "p50", "p90", "rcs_knot_1", "rcs_knot_2", "rcs_knot_3"
    ), drop = FALSE]))) ||
    any(!is.finite(rcs_diagnostics$model_n)) ||
    any(!is.finite(rcs_diagnostics$model_events)) ||
    any(rcs_diagnostics$model_n <= 0) ||
    any(rcs_diagnostics$model_events < 0)) {
  stop_now("Step_04_1_1 RCS diagnostic values failed numeric/reference QC.")
}
if (any(rcs_diagnostics$rcs_knot_1 >= rcs_diagnostics$rcs_knot_2) ||
    any(rcs_diagnostics$rcs_knot_2 >= rcs_diagnostics$rcs_knot_3)) {
  stop_now("Step_04_1_1 RCS diagnostic knots are not strictly increasing.")
}

rcs_diagnostic_import_qc <- tibble(
  check = c(
    "exactly_15_RCS_diagnostic_rows",
    "RCS_diagnostic_order_matches_frozen_15",
    "all_RCS_reference_HR_equal_1",
    "P10_P50_P90_and_three_knots_are_finite",
    "all_three_RCS_knots_are_strictly_increasing",
    "model_N_and_events_are_valid"
  ),
  status = c(
    nrow(rcs_diagnostics) == 15L,
    identical(as.character(rcs_diagnostics$food), FOOD_MAP$food),
    all(abs(rcs_diagnostics$reference_HR - 1) <= 1e-6),
    all(is.finite(as.matrix(rcs_diagnostics[, c(
      "p10", "p50", "p90", "rcs_knot_1", "rcs_knot_2", "rcs_knot_3"
    ), drop = FALSE]))),
    all(rcs_diagnostics$rcs_knot_1 < rcs_diagnostics$rcs_knot_2) &&
      all(rcs_diagnostics$rcs_knot_2 < rcs_diagnostics$rcs_knot_3),
    all(is.finite(rcs_diagnostics$model_n)) &&
      all(is.finite(rcs_diagnostics$model_events)) &&
      all(rcs_diagnostics$model_n > 0) &&
      all(rcs_diagnostics$model_events >= 0)
  )
)
if (anyNA(rcs_diagnostic_import_qc$status) ||
    any(!rcs_diagnostic_import_qc$status)) {
  print(rcs_diagnostic_import_qc)
  stop_now("Imported Step_04_1_1 RCS diagnostic QC failed.")
}

v5_rules <- read.csv(
  RCS_RULE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
required_rule_cols <- c(
  "food", "range_direction", "value1", "value2", "score_lower",
  "score_upper", "score_rule", "score_1", "score_0", "missing_rule"
)
missing_rule_cols <- setdiff(required_rule_cols, names(v5_rules))
if (length(missing_rule_cols) > 0L) {
  stop_now("Final scoring-rule CSV is missing: ",
           paste(missing_rule_cols, collapse = ", "))
}
if (nrow(v5_rules) != 15L ||
    anyDuplicated(as.character(v5_rules$food)) > 0L ||
    !setequal(as.character(v5_rules$food), FOOD_MAP$food)) {
  stop_now("Final scoring-rule CSV does not contain exactly 15 foods.")
}
v5_rules <- v5_rules[
  match(FOOD_MAP$food, as.character(v5_rules$food)), , drop = FALSE
]

same_numeric <- function(a, b, tolerance = 1e-10) {
  a <- suppressWarnings(as.numeric(a))
  b <- suppressWarnings(as.numeric(b))
  if (length(a) != length(b)) return(FALSE)
  same_na <- is.na(a) & is.na(b)
  same_finite <- is.finite(a) & is.finite(b) & abs(a - b) <= tolerance
  same_infinite <- is.infinite(a) & is.infinite(b) & sign(a) == sign(b)
  all(same_na | same_finite | same_infinite)
}
source_rule_interval_match <- identical(
  as.character(v5_rules$range_direction),
  as.character(rcs_ranges$range_direction)
) &&
  identical(as.character(v5_rules$score_rule), as.character(rcs_ranges$score_rule)) &&
  same_numeric(v5_rules$value1, rcs_ranges$value1) &&
  same_numeric(v5_rules$value2, rcs_ranges$value2) &&
  same_numeric(v5_rules$score_lower, rcs_ranges$score_lower) &&
  same_numeric(v5_rules$score_upper, rcs_ranges$score_upper)
if (!source_rule_interval_match) {
  stop_now("Final Step_04_1_1 range CSV and scoring-rule CSV do not reconcile.")
}
source_rcs_ranges <- rcs_ranges
source_v5_rules <- v5_rules

# The two-component interval-derivation audit is supplementary. The final
# range CSV, scoring-rule CSV, and participant-level score file above remain
# the authoritative frozen inputs and are reconciled independently below.
# If the optional audit exists, import and verify it; if absent, skip it without
# stopping the beta model or the newly added RCS/Table 1 modules.
repair_audit_available <- file.exists(REPAIR_AUDIT_FILE)
repair_audit_match <- TRUE
repair_audit_detail <- if (repair_audit_available) {
  "Supplementary two-component interval audit was imported and reconciled."
} else {
  "Supplementary two-component interval audit was not supplied; optional audit skipped. Final range/rule/score files remain authoritative."
}

repair_intervals <- tibble(
  food = character(),
  final_interval_low = numeric(),
  final_interval_high = numeric(),
  candidate_interval_low = numeric(),
  candidate_interval_high = numeric()
)
expected_match <- tibble(
  food = character(),
  candidate_interval_low = numeric(),
  candidate_interval_high = numeric(),
  expected_lower = numeric(),
  expected_upper = numeric()
)

if (repair_audit_available) {
  repair_intervals <- read.csv(
    REPAIR_AUDIT_FILE, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA")
  )
  required_repair_cols <- c("food", "final_interval_low", "final_interval_high")
  if (!all(required_repair_cols %in% names(repair_intervals)) ||
      nrow(repair_intervals) != 2L ||
      !setequal(as.character(repair_intervals$food), REPAIR_FOODS)) {
    stop_now("The supplied optional Step_04_1_1 two-component interval audit is invalid.")
  }
  repair_intervals <- repair_intervals %>%
    mutate(
      candidate_interval_low = as.numeric(final_interval_low),
      candidate_interval_high = as.numeric(final_interval_high)
    )
  expected_match <- repair_intervals %>%
    left_join(EXPECTED_REPAIR_INTERVALS, by = "food")
  if (anyNA(expected_match$expected_lower) ||
      any(abs(expected_match$candidate_interval_low - expected_match$expected_lower) > 1e-6) ||
      any(abs(expected_match$candidate_interval_high - expected_match$expected_upper) > 1e-6)) {
    print(expected_match)
    stop_now("The supplied optional interval audit differs from the frozen expected intersections.")
  }
  for (i in seq_len(nrow(expected_match))) {
    rr <- rcs_ranges[match(expected_match$food[i], rcs_ranges$food), , drop = FALSE]
    if (rr$range_direction[1] != "middle" ||
        abs(rr$value1[1] - expected_match$candidate_interval_low[i]) > 1e-8 ||
        abs(rr$value2[1] - expected_match$candidate_interval_high[i]) > 1e-8) {
      stop_now(
        "Final range CSV does not contain the interval recorded in the optional audit for ",
        expected_match$food[i], "."
      )
    }
  }
}

# Import the final Step_04_1_1 participant-level scores.
rcs_score_data <- read.csv(
  RCS_SCORE_FILE, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA")
)
if (!ID_COL %in% names(rcs_score_data) ||
    !all(SCORE_COLS %in% names(rcs_score_data))) {
  stop_now("Final Step_04_1_1 score file is missing participant ID or component scores.")
}
required_v531_score_summary_cols <- c(
  "n_components_scored", "RWDID_total_score_complete15", "analysis_split"
)
missing_score_summary_cols <- setdiff(
  required_v531_score_summary_cols, names(rcs_score_data)
)
if (length(missing_score_summary_cols) > 0L) {
  stop_now("Final Step_04_1_1 score file is missing: ",
           paste(missing_score_summary_cols, collapse = ", "))
}
rcs_score_data[[ID_COL]] <- trimws(as.character(rcs_score_data[[ID_COL]]))
rcs_score_ids <- sort(rcs_score_data[[ID_COL]], method = "radix")
if (nrow(rcs_score_data) != EXPECTED_TRAIN_N ||
    anyNA(rcs_score_data[[ID_COL]]) ||
    any(rcs_score_data[[ID_COL]] == "") ||
    anyDuplicated(rcs_score_data[[ID_COL]]) > 0L ||
    !identical(rcs_score_ids, training_ids)) {
  stop_now("Final Step_04_1_1 score file does not exactly match frozen training IDs.")
}
if (any(as.character(rcs_score_data$analysis_split) != "training")) {
  stop_now("Final Step_04_1_1 score file contains non-training rows.")
}
score_match <- match(training_data[[ID_COL]], rcs_score_data[[ID_COL]])
if (anyNA(score_match)) stop_now("Failed to align final Step_04_1_1 scores to training participants.")
rcs_score_aligned <- rcs_score_data[score_match, , drop = FALSE]
rcs_score_source_aligned <- rcs_score_aligned

for (score_col in SCORE_COLS) {
  z <- suppressWarnings(as.numeric(rcs_score_aligned[[score_col]]))
  if (any(!is.na(z) & !(z %in% c(0, 1)))) {
    stop_now("Non-binary values found in final score column: ", score_col)
  }
  rcs_score_aligned[[score_col]] <- as.integer(z)
}
source_score_matrix <- as.matrix(rcs_score_source_aligned[, SCORE_COLS, drop = FALSE])
storage.mode(source_score_matrix) <- "numeric"
source_n_scored <- as.integer(rowSums(!is.na(source_score_matrix)))
if (!identical(as.integer(rcs_score_source_aligned$n_components_scored), source_n_scored)) {
  stop_now("Final Step_04_1_1 n_components_scored does not match score columns.")
}
source_complete15 <- source_n_scored == 15L
source_total_recomputed <- rep(NA_integer_, nrow(rcs_score_source_aligned))
if (any(source_complete15)) {
  source_total_recomputed[source_complete15] <- as.integer(rowSums(
    source_score_matrix[source_complete15, , drop = FALSE], na.rm = FALSE
  ))
}
source_total_imported <- suppressWarnings(
  as.numeric(rcs_score_source_aligned$RWDID_total_score_complete15)
)
source_same_total <- (
  is.na(source_total_imported) & is.na(source_total_recomputed)
) | (
  !is.na(source_total_imported) & !is.na(source_total_recomputed) &
    source_total_imported == source_total_recomputed
)
source_score_summary_match <- all(source_same_total)
if (!source_score_summary_match) {
  stop_now("Final Step_04_1_1 complete15 total does not match its component scores.")
}

# Independently recompute all 15 scores from raw training intake and frozen rules.
score_reconciliation <- vector("list", nrow(FOOD_MAP))
for (i in seq_len(nrow(FOOD_MAP))) {
  food_name <- FOOD_MAP$food[i]
  score_col <- FOOD_MAP$score_col[i]
  rr <- rcs_ranges[i, , drop = FALSE]
  recomputed_score <- score_one_component(
    x = training_data[[food_name]],
    direction = as.character(rr$range_direction[1]),
    value1 = as.numeric(rr$value1[1]),
    value2 = as.numeric(rr$value2[1]),
    food_name = food_name
  )
  imported_score <- rcs_score_aligned[[score_col]]
  same_score <- (
    is.na(imported_score) & is.na(recomputed_score)
  ) | (
    !is.na(imported_score) & !is.na(recomputed_score) &
      imported_score == recomputed_score
  )
  if (!all(same_score)) {
    stop_now(
      "Final Step_04_1_1 binary scores do not reconcile with raw training intake for ",
      food_name, "; mismatches=", sum(!same_score), "."
    )
  }
  score_reconciliation[[i]] <- tibble(
    group = FOOD_MAP$group[i],
    food = food_name,
    display_name = FOOD_MAP$display_name[i],
    score_column = score_col,
    final_range_direction = as.character(rr$range_direction[1]),
    final_value1 = as.numeric(rr$value1[1]),
    final_value2 = as.numeric(rr$value2[1]),
    final_range_source = as.character(rr$range_source[1]),
    final_score_rule = as.character(rr$score_rule[1]),
    direction_repair_applied = food_name %in% REPAIR_FOODS,
    n_training = length(imported_score),
    n_missing_imported = sum(is.na(imported_score)),
    n_missing_recomputed = sum(is.na(recomputed_score)),
    n_score_1_imported = sum(imported_score == 1, na.rm = TRUE),
    n_score_1_recomputed = sum(recomputed_score == 1, na.rm = TRUE),
    n_mismatch = sum(!same_score),
    exact_match_to_RCS_frozen_score = all(same_score)
  )
}
score_reconciliation <- bind_rows(score_reconciliation)
unchanged_reconciliation <- score_reconciliation %>%
  filter(!direction_repair_applied)
if (nrow(score_reconciliation) != 15L ||
    !all(score_reconciliation$exact_match_to_RCS_frozen_score) ||
    nrow(unchanged_reconciliation) != 13L ||
    !all(unchanged_reconciliation$exact_match_to_RCS_frozen_score)) {
  stop_now("Final Step_04_1_1 scores failed the independent 15-component reconciliation.")
}

# Refresh score summary columns from the imported final score matrix.
v5_score_matrix <- as.matrix(rcs_score_aligned[, SCORE_COLS, drop = FALSE])
storage.mode(v5_score_matrix) <- "numeric"
v5_n_scored_recomputed <- as.integer(rowSums(!is.na(v5_score_matrix)))
v5_complete15_recomputed <- v5_n_scored_recomputed == 15L
v5_total_recomputed <- rep(NA_integer_, nrow(rcs_score_aligned))
if (any(v5_complete15_recomputed)) {
  v5_total_recomputed[v5_complete15_recomputed] <- as.integer(rowSums(
    v5_score_matrix[v5_complete15_recomputed, , drop = FALSE], na.rm = FALSE
  ))
}
rcs_score_aligned$n_components_scored <- v5_n_scored_recomputed
rcs_score_aligned$RWDID_total_score_complete15 <- v5_total_recomputed
rcs_score_aligned$analysis_split <- "training"
new_component_count_match <- identical(
  as.integer(rcs_score_aligned$n_components_scored), v5_n_scored_recomputed
)
new_total_check <- v5_total_recomputed
new_total_imported <- as.numeric(rcs_score_data$RWDID_total_score_complete15[score_match])
same_total <- (
  is.na(new_total_imported) & is.na(new_total_check)
) | (
  !is.na(new_total_imported) & !is.na(new_total_check) &
    new_total_imported == new_total_check
)
if (!new_component_count_match || !all(same_total)) {
  stop_now("Imported final score summaries failed reconciliation.")
}

# Re-export the imported frozen interval/rule/score files into this beta-run
# directory as provenance copies; no interval or score is modified here.
write_csv_utf8(rcs_ranges, file.path(OUT_DIR, "food_rcs_range_all.csv"))
write_csv_utf8(v5_rules, file.path(OUT_DIR, "RWDID_binary_scoring_rules_15_components.csv"))
write_csv_utf8(rcs_score_aligned, file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv"))
write_csv_utf8(rcs_diagnostics, file.path(OUT_DIR, "RCS_model_summary_and_range_diagnostics_imported.csv"))
write_csv_utf8(rcs_diagnostic_import_qc, file.path(AUDIT_DIR, "RCS_diagnostic_import_QC.csv"))
write_csv_utf8(repair_intervals, file.path(AUDIT_DIR, "two_component_interval_derivation.csv"))
write_csv_utf8(score_reconciliation, file.path(AUDIT_DIR, "participant_score_reconciliation.csv"))

write_csv_utf8(
  tibble(
    check = c(
      "final_RCS_range_and_rule_files_match",
      "final_RCS_score_summary_matches_source_score_matrix",
      "optional_repair_audit_matches_frozen_expected_intersections_if_supplied",
      "all_15_final_RCS_scores_reproduce_from_raw_training_intake",
      "n_components_scored_matches_recomputed_scores",
      "complete15_total_matches_recomputed_scores"
    ),
    status = c(
      source_rule_interval_match,
      source_score_summary_match,
      repair_audit_match,
      nrow(score_reconciliation) == 15L &&
        all(score_reconciliation$exact_match_to_RCS_frozen_score),
      new_component_count_match,
      all(same_total)
    ),
    detail = c(
      "Step_04_1_1 final range and rule exports reconcile",
      "Imported final total score summary reconciles with component score columns",
            repair_audit_detail,
      paste0("exact final score matches=", sum(score_reconciliation$exact_match_to_RCS_frozen_score), "/15"),
      "n_components_scored independently recomputed",
      "complete15 total independently recomputed"
    )
  ),
  file.path(AUDIT_DIR, "input_reconciliation_QC.csv")
)

# -----------------------------------------------------------------------------
# 4. Assemble the training analysis dataset using all 15 binary scores
# -----------------------------------------------------------------------------

analysis_data <- training_data
for (i in seq_len(nrow(FOOD_MAP))) {
  analysis_data[[FOOD_MAP$model_var[i]]] <-
    as.integer(rcs_score_aligned[[FOOD_MAP$score_col[i]]])
}
analysis_data$analysis_split <- "training"

score_matrix <- as.matrix(
  analysis_data[, FOOD_VARS, drop = FALSE]
)
storage.mode(score_matrix) <- "numeric"

analysis_data$n_components_scored <- as.integer(
  rowSums(!is.na(score_matrix))
)
complete_primary <- analysis_data$n_components_scored == N_COMPONENTS

# Backward-compatible field names are retained for the existing validation
# workflow. All 15 V5 components are active; there is no eligibility subset.
analysis_data$n_V5_score_components <- N_COMPONENTS
analysis_data$n_primary_score_components <- N_COMPONENTS
analysis_data$RWDID_total_score_complete15 <- NA_integer_
if (any(complete_primary)) {
  analysis_data$RWDID_total_score_complete15[complete_primary] <-
    as.integer(rowSums(
      score_matrix[complete_primary, , drop = FALSE],
      na.rm = FALSE
    ))
}
analysis_data$RWDID_total_score_complete_primary <-
  analysis_data$RWDID_total_score_complete15

if (any(
  analysis_data$RWDID_total_score_complete15[complete_primary] < 0 |
    analysis_data$RWDID_total_score_complete15[complete_primary] > 15L,
  na.rm = TRUE
)) {
  stop_now("The complete 15-component unweighted total is outside 0-15.")
}
if (sum(complete_primary) == 0L) {
  stop_now("No participant has complete binary scores across all 15 V5 components.")
}

# -----------------------------------------------------------------------------
# 5. Estimate mutually adjusted Model 3 Cox beta coefficients for all 15 V5.3.1
#    components in the frozen TRAINING cohort only
# -----------------------------------------------------------------------------

message("\n[1/7] Estimating mutually adjusted Model 3 beta coefficients for all 15 components ...")

beta_model_vars <- c(
  "depression_time", "depression_status", FOOD_VARS, COVARS
)
beta_complete <- complete.cases(
  analysis_data[, beta_model_vars, drop = FALSE]
)
beta_model_data <- analysis_data[beta_complete, , drop = FALSE]

if (nrow(beta_model_data) == 0L ||
    sum(beta_model_data$depression_status == 1, na.rm = TRUE) == 0L) {
  stop_now("No complete-case training participants/events for the 15-component beta model.")
}

component_variation <- tibble(
  term = FOOD_VARS,
  n_zero = vapply(
    FOOD_VARS,
    function(v) sum(beta_model_data[[v]] == 0, na.rm = TRUE),
    numeric(1)
  ),
  n_one = vapply(
    FOOD_VARS,
    function(v) sum(beta_model_data[[v]] == 1, na.rm = TRUE),
    numeric(1)
  )
) %>%
  left_join(
    FOOD_MAP %>% select(model_var, group, food),
    by = c("term" = "model_var")
  ) %>%
  select(group, food, term, n_zero, n_one)

if (any(component_variation$n_zero == 0L) ||
    any(component_variation$n_one == 0L)) {
  print(component_variation)
  stop_now(
    "At least one V5 binary component has no variation in the beta-estimation ",
    "sample. Inspect the component variation audit; do not force an ",
    "unidentifiable coefficient."
  )
}

beta_formula <- make_survival_formula(FOOD_VARS, COVARS)
fit_beta <- fit_cox_model(
  beta_formula,
  beta_model_data,
  "mutually_adjusted_Model3_Cox_beta_model_all_15_V5_components"
)

beta_tidy_all <- broom::tidy(
  fit_beta,
  exponentiate = FALSE,
  conf.int = TRUE,
  conf.level = 0.95
)

beta_results <- beta_tidy_all %>%
  filter(.data$term %in% FOOD_VARS) %>%
  transmute(
    term = as.character(.data$term),
    estimate = as.numeric(.data$estimate),
    std.error = as.numeric(.data$std.error),
    statistic = as.numeric(.data$statistic),
    p.value = as.numeric(.data$p.value),
    conf.low = as.numeric(.data$conf.low),
    conf.high = as.numeric(.data$conf.high),
    hazard_ratio = exp(.data$estimate),
    hr_ci_low = exp(.data$conf.low),
    hr_ci_high = exp(.data$conf.high)
  ) %>%
  mutate(
    p_fdr_BH = p.adjust(p.value, method = "BH"),
    p_bonf = p.adjust(p.value, method = "bonferroni"),
    beta_order = match(term, FOOD_VARS)
  ) %>%
  left_join(
    FOOD_MAP %>% select(model_var, group, food, display_name),
    by = c("term" = "model_var")
  ) %>%
  arrange(beta_order) %>%
  select(
    beta_order, group, food, display_name, term,
    estimate, std.error, statistic, p.value,
    p_fdr_BH, p_bonf,
    hazard_ratio, hr_ci_low, hr_ci_high,
    conf.low, conf.high
  )

if (nrow(beta_results) != N_COMPONENTS ||
    anyDuplicated(beta_results$term) > 0L ||
    !setequal(beta_results$term, FOOD_VARS) ||
    any(!is.finite(beta_results$estimate)) ||
    anyNA(beta_results$estimate)) {
  stop_now(
    "The mutually adjusted Cox model did not return one finite beta for ",
    "each of the 15 components."
  )
}

beta_vector <- beta_results$estimate[
  match(FOOD_VARS, beta_results$term)
]
names(beta_vector) <- FOOD_VARS

beta_sum <- sum(beta_vector)
beta_abs_sum <- sum(abs(beta_vector))
beta_sum_tolerance <- max(1e-10, 1e-8 * beta_abs_sum)

beta_sign_pattern <- if (all(beta_vector < 0)) {
  "all_negative_inverse_associations"
} else if (all(beta_vector > 0)) {
  "all_positive_direct_associations"
} else if (any(beta_vector < 0) && any(beta_vector > 0)) {
  "mixed_positive_and_negative_associations"
} else {
  "includes_zero_beta"
}

diet_score_direction_note <- if (all(beta_vector < 0)) {
  "Higher values reflect greater weighted adherence to score-1 intervals with inverse associations."
} else {
  "Beta signs are not uniformly negative; do not automatically describe higher scores as protective."
}

if (!is.finite(beta_sum) || abs(beta_sum) <= beta_sum_tolerance) {
  stop_now(
    "The sum of the 15 beta coefficients is numerically zero or unstable. ",
    "The original diet_score = weighted_sum / sum(beta) formula cannot be ",
    "applied safely. Inspect beta_weight_audit.csv."
  )
}

beta_weight_audit <- beta_results %>%
  transmute(
    beta_order,
    group,
    food,
    term,
    beta = estimate,
    beta_direction = case_when(
      estimate > 0 ~ "positive",
      estimate < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    hazard_ratio,
    p_value = p.value,
    p_fdr_BH,
    p_bonf
  )

beta_weight_audit <- bind_rows(
  beta_weight_audit,
  tibble(
    beta_order = N_COMPONENTS + 1L,
    group = "SUM OF ALL 15 COMPONENT BETA COEFFICIENTS",
    food = NA_character_,
    term = "sum_beta",
    beta = beta_sum,
    beta_direction = ifelse(beta_sum > 0, "positive_sum", "negative_sum"),
    hazard_ratio = NA_real_,
    p_value = NA_real_,
    p_fdr_BH = NA_real_,
    p_bonf = NA_real_
  )
)

write_csv_utf8(
  beta_results,
  file.path(OUT_DIR, "beta_coefficients_training.csv")
)

beta_sign_audit <- beta_results %>%
  transmute(
    group, food, display_name, term,
    beta = estimate,
    hazard_ratio,
    beta_sign = case_when(
      estimate < 0 ~ "negative_inverse_association",
      estimate > 0 ~ "positive_direct_association",
      TRUE ~ "zero"
    ),
    p_value = p.value,
    p_fdr_BH,
    p_bonf
  )
write_csv_utf8(
  beta_sign_audit,
  file.path(OUT_DIR, "beta_sign_audit.csv")
)

component_manifest_used <- FOOD_MAP %>%
  left_join(
    rcs_ranges %>%
      select(
        food, value1, value2, range_direction, range_source,
        score_lower, score_upper, score_rule
      ),
    by = "food"
  ) %>%
  left_join(
    v5_rules %>%
      select(food, interval_source, score_1, score_0, missing_rule),
    by = "food"
  ) %>%
  arrange(match(food, FOOD_MAP$food))

write_csv_utf8(
  component_manifest_used,
  file.path(OUT_DIR, "component_manifest_used.csv")
)
write_csv_utf8(
  beta_weight_audit,
  file.path(OUT_DIR, "beta_weight_audit.csv")
)
write_csv_utf8(
  component_variation,
  file.path(AUDIT_DIR, "beta_component_variation.csv")
)
save(
  beta_results,
  file = file.path(OUT_DIR, "beta_results_depression_food_group.RData")
)

beta_model_sample_summary <- tibble(
  analysis = "Mutually adjusted Model 2 Cox beta estimation for all 15 components",
  cohort = "Frozen Step_03_1_0 70% training cohort",
  n_training_rds = nrow(analysis_data),
  n_complete_case_beta_model = nrow(beta_model_data),
  events_beta_model = sum(beta_model_data$depression_status == 1, na.rm = TRUE),
  censored_beta_model = sum(beta_model_data$depression_status == 0, na.rm = TRUE),
  number_of_components = N_COMPONENTS,
  number_of_beta_coefficients = nrow(beta_results),
  sum_beta = beta_sum,
  sum_absolute_beta = beta_abs_sum,
  denominator_tolerance = beta_sum_tolerance,
  denominator_sign = ifelse(beta_sum > 0, "positive", "negative"),
  beta_sign_pattern = beta_sign_pattern,
  diet_score_direction_note = diet_score_direction_note
)
write_csv_utf8(
  beta_model_sample_summary,
  file.path(OUT_DIR, "beta_model_sample_summary.csv")
)

# 6. Calculate the original beta-weighted diet_score
# -----------------------------------------------------------------------------

message("\n[2/7] Constructing beta-weighted diet_score ...")

# This is the original user-specified formula:
# weighted_sum = sum(component_binary_score * beta_i)
# diet_score   = weighted_sum / sum(beta_i)
#
# Incomplete 15-component rows remain NA. No missing component is imputed;
# a complete score across all 15 V5 components is required.
weighted_sum <- rep(NA_real_, nrow(analysis_data))
diet_score <- rep(NA_real_, nrow(analysis_data))

if (any(complete_primary)) {
  complete_matrix <- score_matrix[complete_primary, , drop = FALSE]
  weighted_sum[complete_primary] <- as.numeric(
    complete_matrix %*% beta_vector
  )
  diet_score[complete_primary] <- weighted_sum[complete_primary] / beta_sum
}

analysis_data$weighted_sum <- weighted_sum
analysis_data$diet_score <- diet_score

if (any(!is.finite(analysis_data$diet_score[complete_primary]))) {
  stop_now("Non-finite beta-weighted diet_score values were generated.")
}
if (sum(complete_primary) == 0L) {
  stop_now("No participant has complete data for all 15 V5 binary components.")
}

# Training-derived distribution and tertile thresholds. The tertile cutpoints
# are frozen here and will be applied unchanged to the held-out test cohort.
score_values <- analysis_data$diet_score[
  is.finite(analysis_data$diet_score)
]
score_min <- min(score_values)
score_max <- max(score_values)
score_mean <- mean(score_values)
score_sd <- stats::sd(score_values)

if (!is.finite(score_min) || !is.finite(score_max) ||
    !is.finite(score_mean) || !is.finite(score_sd) ||
    score_max <= score_min || score_sd <= 0) {
  stop_now("The training beta-weighted diet_score has no usable variation.")
}

tertile_quantiles <- as.numeric(stats::quantile(
  score_values,
  probs = c(1 / 3, 2 / 3),
  na.rm = TRUE,
  names = FALSE,
  type = 7
))
tertile_q1 <- tertile_quantiles[1]
tertile_q2 <- tertile_quantiles[2]

if (!all(is.finite(c(tertile_q1, tertile_q2))) ||
    tertile_q1 >= tertile_q2) {
  stop_now(
    "Training beta-weighted diet_score does not support three distinct ",
    "tertiles without splitting tied values. Inspect its distribution."
  )
}

analysis_data$diet_score_tertile <- NA_character_
analysis_data$diet_score_tertile[complete_primary] <- as.character(cut(
  analysis_data$diet_score[complete_primary],
  breaks = c(-Inf, tertile_q1, tertile_q2, Inf),
  include.lowest = TRUE,
  right = TRUE,
  labels = c("T1", "T2", "T3")
))
analysis_data$diet_score_tertile <- factor(
  analysis_data$diet_score_tertile,
  levels = c("T1", "T2", "T3")
)
analysis_data$diet_score_tertile_num <- as.numeric(
  analysis_data$diet_score_tertile
)

# Preserve the original definition: 20% of the observed training score range.
range_20 <- 0.20 * (score_max - score_min)
if (!is.finite(range_20) || range_20 <= 0) {
  stop_now("The training 20%-of-range increment is invalid.")
}
analysis_data$diet_score_20 <- NA_real_
analysis_data$diet_score_20[complete_primary] <-
  (analysis_data$diet_score[complete_primary] - score_mean) / range_20

# Retain the original quartile variable for the supplementary baseline table;
# quartiles are descriptive only and are not part of the requested Cox models.
quartile_quantiles <- as.numeric(stats::quantile(
  score_values,
  probs = c(0.25, 0.50, 0.75),
  na.rm = TRUE,
  names = FALSE,
  type = 7
))
quartile_available <- all(is.finite(quartile_quantiles)) &&
  length(unique(quartile_quantiles)) == 3L

if (quartile_available) {
  analysis_data$diet_score_quartile <- NA_character_
  analysis_data$diet_score_quartile[complete_primary] <- as.character(cut(
    analysis_data$diet_score[complete_primary],
    breaks = c(-Inf, quartile_quantiles, Inf),
    include.lowest = TRUE,
    right = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4")
  ))
  analysis_data$diet_score_quartile <- factor(
    analysis_data$diet_score_quartile,
    levels = c("Q1", "Q2", "Q3", "Q4")
  )
} else {
  analysis_data$diet_score_quartile <- NA
  warning(
    "Quartile cutpoints contain ties; quartile baseline table will be skipped. ",
    "The requested tertile association analysis is unaffected."
  )
}

# Freeze all parameters needed for Step_04_2_2. No parameter is estimated in
# the held-out test cohort.
scaling_parameters <- tibble(
  parameter = c(
    "sum_beta",
    "sum_absolute_beta",
    "diet_score_min_training",
    "diet_score_max_training",
    "diet_score_mean_training",
    "diet_score_sd_training",
    "tertile_cutpoint_1_training",
    "tertile_cutpoint_2_training",
    "per20_increment_training",
    "quantile_type",
    "n_training_rds",
    "n_V5_score_components",
    "n_primary_score_components",
    "n_training_complete_primary",
    "n_training_association_model"
  ),
  value = c(
    beta_sum,
    beta_abs_sum,
    score_min,
    score_max,
    score_mean,
    score_sd,
    tertile_q1,
    tertile_q2,
    range_20,
    7,
    nrow(analysis_data),
    N_COMPONENTS,
    N_COMPONENTS,
    sum(complete_primary),
    NA_real_
  ),
  definition = c(
    "Algebraic sum of mutually adjusted Cox beta coefficients for 15 V5 components",
    "Sum of absolute values of beta coefficients for 15 V5 components",
    "Minimum beta-weighted score among participants complete on all 15 V5 components",
    "Maximum beta-weighted score among participants complete on all 15 V5 components",
    "Mean beta-weighted score among participants complete on all 15 V5 components",
    "SD beta-weighted score among participants complete on all 15 V5 components",
    "Type-7 empirical 33.333rd percentile; upper boundary of T1",
    "Type-7 empirical 66.667th percentile; upper boundary of T2",
    "0.20 multiplied by the observed training score range",
    "R quantile type used for tertile thresholds",
    "Number of participants in the frozen Step_03_1_0 training RDS",
    "Number of components in the V5 RWDID; all 15 components are included and no additional V6.1 eligibility gate is applied",
    "Backward-compatible alias; equals the 15 V5 components, not an eligibility subset",
    "Number with all 15 V5 binary component scores observed; missing values are not imputed",
    "Filled after the common complete-case association cohort is defined"
  )
)

# Participant-level training score output. Component scores are copied from
# Step_04_1_1 V5 and have passed independent V5 interval/rule reconciliation.
participant_score_output <- analysis_data %>%
  select(
    all_of(ID_COL),
    all_of(ALL_FOOD_VARS),
    n_components_scored,
    n_V5_score_components,
    n_primary_score_components,
    RWDID_total_score_complete_primary,
    RWDID_total_score_complete15,
    weighted_sum,
    diet_score,
    diet_score_tertile,
    diet_score_tertile_num,
    diet_score_20,
    analysis_split
  )

write_csv_utf8(
  participant_score_output,
  file.path(OUT_DIR, "beta_weighted_diet_score_training.csv")
)

write_csv_utf8(
  component_manifest_used,
  file.path(OUT_DIR, "RCS_interval_manifest_used.csv")
)
write_csv_utf8(
  v5_rules,
  file.path(OUT_DIR, "binary_scoring_rules_used.csv")
)

# -----------------------------------------------------------------------------
# 7. Define one common complete-case association cohort and fit Cox models
# -----------------------------------------------------------------------------

message("\n[3/7] Fitting continuous, tertile, trend, and per-20% Cox models ...")

association_vars <- c(
  "depression_time",
  "depression_status",
  "diet_score",
  "diet_score_tertile",
  "diet_score_tertile_num",
  "diet_score_20",
  COVARS
)

association_complete <- complete.cases(
  analysis_data[, association_vars, drop = FALSE]
)
association_data <- analysis_data[association_complete, , drop = FALSE]

if (nrow(association_data) == 0L ||
    sum(association_data$depression_status == 1, na.rm = TRUE) == 0L) {
  stop_now("No eligible complete-case training participants/events for association models.")
}
if (!identical(
  sort(association_data[[ID_COL]], method = "radix"),
  sort(beta_model_data[[ID_COL]], method = "radix")
)) {
  stop_now(
    "The beta-estimation cohort and diet-score association cohort differ. ",
    "The score construction/complete-case handling requires inspection."
  )
}
if (!all(c("T1", "T2", "T3") %in%
         as.character(unique(association_data$diet_score_tertile)))) {
  stop_now("All three training tertiles must be represented in the association cohort.")
}

tertile_counts <- association_data %>%
  group_by(diet_score_tertile) %>%
  summarise(
    total = n(),
    events = sum(depression_status == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_fmt = formatC(total, format = "d", big.mark = ","),
    events_fmt = formatC(events, format = "d", big.mark = ","),
    Events_Total = paste0(events_fmt, "/", total_fmt)
  )

if (nrow(tertile_counts) != 3L || any(tertile_counts$total == 0L)) {
  stop_now("Tertile distribution does not contain exactly three non-empty groups.")
}

fit_continuous <- fit_cox_model(
  make_survival_formula("diet_score", COVARS),
  association_data,
  "continuous_diet_score"
)

fit_tertile <- fit_cox_model(
  make_survival_formula("diet_score_tertile", COVARS),
  association_data,
  "diet_score_tertiles"
)

fit_trend <- fit_cox_model(
  make_survival_formula("diet_score_tertile_num", COVARS),
  association_data,
  "diet_score_tertile_trend"
)

fit_per20 <- fit_cox_model(
  make_survival_formula("diet_score_20", COVARS),
  association_data,
  "diet_score_per_20_percent_training_range"
)

# Continuous score: per one-unit increase.
continuous_row <- make_result_row(
  variable = "RWDID (per score increase)",
  term = "diet_score",
  fit = fit_continuous,
  data = association_data
)

# Tertile Cox model: T1 is the reference group, exactly as in the original code.
t1_row <- make_result_row(
  variable = "T1",
  term = NA_character_,
  fit = fit_tertile,
  data = association_data %>% filter(diet_score_tertile == "T1"),
  count_data = association_data %>% filter(diet_score_tertile == "T1"),
  reference = TRUE
)
t2_row <- make_result_row(
  variable = "T2",
  term = "diet_score_tertileT2",
  fit = fit_tertile,
  data = association_data %>% filter(diet_score_tertile == "T2"),
  count_data = association_data %>% filter(diet_score_tertile == "T2")
)
t3_row <- make_result_row(
  variable = "T3",
  term = "diet_score_tertileT3",
  fit = fit_tertile,
  data = association_data %>% filter(diet_score_tertile == "T3"),
  count_data = association_data %>% filter(diet_score_tertile == "T3")
)

# P for trend: T1/T2/T3 are coded 1/2/3, preserving the original approach.
trend_tidy <- extract_tidy_term(fit_trend, "diet_score_tertile_num")
trend_row <- tibble(
  variable = "P for trend",
  Events_Total = "",
  hr_ci = "",
  p_value = fmt_p(trend_tidy$p.value[1]),
  hazard_ratio = NA_real_,
  ci_low = NA_real_,
  ci_high = NA_real_,
  p_value_numeric = as.numeric(trend_tidy$p.value[1]),
  term = "diet_score_tertile_num"
)

# Per 20% increment: 20% of the observed training score range.
per20_row <- make_result_row(
  variable = "Per 20% increment",
  term = "diet_score_20",
  fit = fit_per20,
  data = association_data
)

diet_score_all <- bind_rows(
  continuous_row,
  t1_row,
  t2_row,
  t3_row,
  trend_row,
  per20_row
) %>%
  mutate(
    model_n = nrow(association_data),
    model_events = sum(association_data$depression_status == 1, na.rm = TRUE),
    analysis_cohort = "Frozen 70% training"
  )

# Update the frozen scaling file with the final association sample size.
scaling_parameters$value[
  scaling_parameters$parameter == "n_training_association_model"
] <- nrow(association_data)

write_csv_utf8(
  scaling_parameters,
  file.path(OUT_DIR, "diet_score_scaling_parameters.csv")
)
write_csv_utf8(
  diet_score_all,
  file.path(OUT_DIR, "diet_score_association_results.csv")
)

# -----------------------------------------------------------------------------
# 8. Forest plot: retain the original forestplot layout and aesthetics
# -----------------------------------------------------------------------------

message("\n[4/7] Creating the original-style forest plot ...")

forest_final_result <- diet_score_all

mycolors <- cols4all::c4a("friendly13", 7)

# The original five-colour cycle and circle confidence-interval markers are
# retained. The P-for-trend row and T1 reference row have no CI marker.
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker) {
    i <<- i + 1
    if (i %% 5 == 1) {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5, clr.line = mycolors[1], clr.marker = mycolors[1]
      )
    } else if (i %% 5 == 2) {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5, clr.line = mycolors[2], clr.marker = mycolors[2]
      )
    } else if (i %% 5 == 3) {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5, clr.line = mycolors[3], clr.marker = mycolors[3]
      )
    } else if (i %% 5 == 4) {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5, clr.line = mycolors[4], clr.marker = mycolors[4]
      )
    } else {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5, clr.line = mycolors[5], clr.marker = mycolors[5]
      )
    }
  }
})

tabletext_HR <- cbind(
  c("Groups", as.character(forest_final_result$variable)),
  c("Events/Total", as.character(forest_final_result$Events_Total)),
  c("HR (95% CI)", as.character(forest_final_result$hr_ci)),
  c("P value", as.character(forest_final_result$p_value))
)

model1_forestplot <- forestplot::forestplot(
  labeltext = tabletext_HR,
  graph.pos = 4,
  fn.ci_norm = fn,
  boxsize = 0.2,
  mean = c(NA_real_, forest_final_result$hazard_ratio),
  lower = c(NA_real_, forest_final_result$ci_low),
  upper = c(NA_real_, forest_final_result$ci_high),
  xlog = TRUE,
  col = forestplot::fpColors(zero = "black"),
  zero = 1,
  graphwidth = grid::unit(0.10, "npc"),
  lty.zero = 3,
  lwd.zero = 1.5,
  xlab = "HR (95% CI)",
  lwd.xaxis = 1.5,
  xticks = c(0.12, 0.25, 0.5, 1, 2),
  txt_gp = forestplot::fpTxtGp(
    label = grid::gpar(cex = 1.1),
    ticks = grid::gpar(cex = 0.8),
    xlab = grid::gpar(cex = 0.9)
  )
)

grDevices::pdf(
  file = file.path(OUT_DIR, "diet_score_forestplot.pdf"),
  width = 15,
  height = 3
)
print(model1_forestplot)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# 9. Beta coefficient heatmap: diverging palette by estimated beta sign
# -----------------------------------------------------------------------------

message("
[5/7] Creating the beta coefficient heatmap ...")

beta_results_plot <- beta_results %>%
  mutate(
    variable = factor(
      display_name,
      levels = rev(FOOD_MAP$display_name)
    ),
    dummy = "estimate",
    estimate_label = formatC(estimate, format = "f", digits = 3)
  )

# Beta signs are estimated, not constrained. Negative beta means the
# component's score=1 category is associated with lower hazard conditional
# on other active components and Model 3 covariates; positive beta means the
# opposite. A diverging palette prevents a one-sided visual interpretation.
if (!all(is.finite(beta_results_plot$estimate))) {
  stop_now("Non-finite beta estimates cannot be plotted.")
}
beta_abs_limit <- max(abs(beta_results_plot$estimate))
if (!is.finite(beta_abs_limit) || beta_abs_limit <= 0) {
  stop_now("The active beta coefficients have no non-zero range.")
}
beta_limits <- c(-beta_abs_limit, beta_abs_limit)

p_heatmap <- ggplot(
  beta_results_plot,
  aes(x = dummy, y = variable, fill = estimate)
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = estimate_label),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "PiYG")[6:11]),
    name = "Estimate"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Mutually adjusted beta coefficients",
    subtitle = "Green: inverse association; magenta: positive association"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", colour = "black"),
    axis.title.x = element_text(face = "bold", colour = "black"),
    axis.title.y = element_text(face = "bold", colour = "black"),
    panel.grid = element_blank()
  )

grDevices::pdf(
  file.path(OUT_DIR, "beta_results.pdf"),
  width = 4.2,
  height = 4.2
)
print(p_heatmap)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# 10. Score distributions and cohort summaries
# -----------------------------------------------------------------------------

message("\n[6/7] Exporting score distributions, model diagnostics, and tables ...")

diet_score_distribution <- tibble(
  cohort = "Frozen 70% training",
  n_training_rds = nrow(analysis_data),
  n_components_in_V5_score = N_COMPONENTS,
  n_complete_15_component_scores = sum(complete_primary),
  n_incomplete_15_component_scores = sum(!complete_primary),
  n_association_model = nrow(association_data),
  events_association_model = sum(
    association_data$depression_status == 1, na.rm = TRUE
  ),
  mean = mean(score_values),
  sd = stats::sd(score_values),
  median = stats::median(score_values),
  p25 = as.numeric(stats::quantile(score_values, 0.25, names = FALSE)),
  p75 = as.numeric(stats::quantile(score_values, 0.75, names = FALSE)),
  minimum = min(score_values),
  maximum = max(score_values),
  sum_beta = beta_sum,
  range_20 = range_20
)

write_csv_utf8(
  diet_score_distribution,
  file.path(OUT_DIR, "diet_score_distribution.csv")
)

diet_score_tertile_distribution <- analysis_data %>%
  filter(!is.na(diet_score_tertile)) %>%
  group_by(diet_score_tertile) %>%
  summarise(
    total = n(),
    events = sum(depression_status == 1, na.rm = TRUE),
    mean_diet_score = mean(diet_score),
    median_diet_score = median(diet_score),
    min_diet_score = min(diet_score),
    max_diet_score = max(diet_score),
    .groups = "drop"
  ) %>%
  mutate(
    Events_Total = paste0(
      formatC(events, format = "d", big.mark = ","),
      "/",
      formatC(total, format = "d", big.mark = ",")
    )
  )

write_csv_utf8(
  diet_score_tertile_distribution,
  file.path(OUT_DIR, "diet_score_tertile_distribution.csv")
)

component_score_distribution <- bind_rows(lapply(seq_len(nrow(FOOD_MAP)), function(i) {
  v <- FOOD_MAP$model_var[i]
  x <- analysis_data[[v]]
  n_obs <- sum(!is.na(x))
  n0 <- sum(x == 0, na.rm = TRUE)
  n1 <- sum(x == 1, na.rm = TRUE)

  tibble(
    group = FOOD_MAP$group[i],
    food = FOOD_MAP$food[i],
    display_name = FOOD_MAP$display_name[i],
    range_direction = rcs_ranges$range_direction[
      match(FOOD_MAP$food[i], rcs_ranges$food)
    ],
    range_source = rcs_ranges$range_source[
      match(FOOD_MAP$food[i], rcs_ranges$food)
    ],
    score_rule = rcs_ranges$score_rule[
      match(FOOD_MAP$food[i], rcs_ranges$food)
    ],
    n_training = length(x),
    n_missing = sum(is.na(x)),
    n_scored = n_obs,
    n_score_0 = n0,
    n_score_1 = n1,
    pct_score_0_among_scored = ifelse(n_obs > 0, 100 * n0 / n_obs, NA_real_),
    pct_score_1_among_scored = ifelse(n_obs > 0, 100 * n1 / n_obs, NA_real_)
  )
}))

write_csv_utf8(
  component_score_distribution,
  file.path(OUT_DIR, "binary_component_score_distribution.csv")
)

# Cox proportional-hazards diagnostics for the beta model and four association
# models. These diagnostics do not change the model specifications.
ph_diagnostics <- bind_rows(
  safe_cox_zph(fit_beta, "mutually_adjusted_beta_model"),
  safe_cox_zph(fit_continuous, "continuous_diet_score"),
  safe_cox_zph(fit_tertile, "diet_score_tertiles"),
  safe_cox_zph(fit_trend, "diet_score_tertile_trend"),
  safe_cox_zph(fit_per20, "diet_score_per_20_percent_range")
)

write_csv_utf8(
  ph_diagnostics,
  file.path(OUT_DIR, "cox_ph_diagnostics.csv")
)

# Baseline characteristics tables from the original script are retained for
# the training cohort. They are descriptive only and do not alter the Cox
# models. If optional table packages or source variables are unavailable, the
# V5 score/association outputs remain valid and the omission is recorded.
table1_status <- "NOT_RUN"
table1_reason <- ""

table1_vars <- c(
  "age_food", "follow_up_years", "sex", "ethnicity",
  "education_attainment", "career", "twnsend", "smoking_status",
  "drink_status", "bmi_group", "household_income",
  "physical_activity_proup", "sedentary_group", "energy_intake",
  "central_obesity", "triglycerides", "hdl_cholesterol",
  "hypertension_status_final", "t2diabetes_status_final",
  "dyslipidemia_status_final", "cvd_status_final",
  "all_cancer_status_final"
)

if (requireNamespace("gtsummary", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE) &&
    all(table1_vars %in% names(association_data))) {

  table_data <- association_data %>%
    mutate(
      depression_status_label = factor(
        ifelse(depression_status == 1, "Depression", "No depression"),
        levels = c("No depression", "Depression")
      )
    )

  if (all(as.character(table_data$sex) %in% c("Male", "Female"))) {
    table_data$sex <- factor(
      as.character(table_data$sex),
      levels = c("Male", "Female")
    )
  }
  if (all(as.character(table_data$bmi_group) %in% c(
    "Underweight (<18.5)",
    "Normal weight (18.5 to <25)",
    "Overweight (25 to <30)",
    "Obese (>=30)"
  ))) {
    table_data$bmi_group <- factor(
      as.character(table_data$bmi_group),
      levels = c(
        "Underweight (<18.5)",
        "Normal weight (18.5 to <25)",
        "Overweight (25 to <30)",
        "Obese (>=30)"
      )
    )
  }
  if (all(as.character(table_data$household_income) %in%
          c("Low", "Middle", "High"))) {
    table_data$household_income <- factor(
      as.character(table_data$household_income),
      levels = c("Low", "Middle", "High")
    )
  }

  table1_labels <- list(
    age_food ~ "Age, year, median (IQR)",
    follow_up_years ~ "Follow-up, years, median (IQR)",
    sex ~ "Sex, no. (%)",
    ethnicity ~ "Ethnicity, no. (%)",
    education_attainment ~ "Education attainment",
    career ~ "Employment status, no. (%)",
    twnsend ~ "Townsend deprivation index",
    energy_intake ~ "Energy intake, KJ/day",
    smoking_status ~ "Smoking status, no. (%)",
    drink_status ~ "Drinking status, no. (%)",
    bmi_group ~ "BMI category (kg/m2), no. (%)",
    household_income ~ "Family income, no. (%)",
    physical_activity_proup ~ "Physical activity, no. (%)",
    sedentary_group ~ "Sedentary behavior, no. (%)",
    central_obesity ~ "Central obesity, no. (%)",
    triglycerides ~ "High triglycerides, no. (%)",
    hdl_cholesterol ~ "Low HDL, no. (%)",
    hypertension_status_final ~ "Hypertension, no. (%)",
    t2diabetes_status_final ~ "Diabetes, no. (%)",
    dyslipidemia_status_final ~ "Dyslipidemia, no. (%)",
    cvd_status_final ~ "CVD, no. (%)",
    all_cancer_status_final ~ "Cancer, no. (%)"
  )

  make_table1 <- function(data, by_var) {
    # Retain the original Table 1 continuous2 specification.
    type_spec <- list(
      c(age_food, follow_up_years, education_attainment, twnsend, energy_intake) ~
        "continuous2"
    )

    tbl <- gtsummary::tbl_summary(
      data = data,
      by = by_var,
      include = all_of(table1_vars),
      label = table1_labels,
      type = type_spec,
      digits = all_continuous() ~ 2,
      statistic = list(
        all_continuous() ~ c("{mean} ({sd})", "{median} ({p25}, {p75})")
      ),
      missing = "no"
    ) %>%
      gtsummary::add_overall() %>%
      gtsummary::add_n(
        statistic = "{N_nonmiss}",
        col_label = "**N**",
        footnote = TRUE
      )
    tbl
  }

  tryCatch({
    tbl_depression <- make_table1(table_data, "depression_status_label")
    ft_depression <- gtsummary::as_flex_table(tbl_depression)
    flextable::save_as_docx(
      "Table 1" = ft_depression,
      path = file.path(OUT_DIR, "Table_1_baseline_by_depression.docx")
    )

    tbl_tertile <- make_table1(table_data, "diet_score_tertile")
    ft_tertile <- gtsummary::as_flex_table(tbl_tertile)
    flextable::save_as_docx(
      "Table 1" = ft_tertile,
      path = file.path(OUT_DIR, "Table_1_baseline_by_tertile.docx")
    )

    if (quartile_available) {
      tbl_quartile <- make_table1(table_data, "diet_score_quartile")
      ft_quartile <- gtsummary::as_flex_table(tbl_quartile)
      flextable::save_as_docx(
        "Table 1" = ft_quartile,
        path = file.path(OUT_DIR, "Table_1_baseline_by_quartile.docx")
      )
    }

    table1_status <- "COMPLETED"
  }, error = function(e) {
    table1_status <<- "FAILED_OPTIONAL_TABLES"
    table1_reason <<- conditionMessage(e)
    warning(
      "Baseline Table 1 export failed; core Cox and forestplot outputs ",
      "have already been generated. Reason: ", table1_reason
    )
  })

} else {
  table1_status <- "SKIPPED_OPTIONAL_TABLES"
  table1_reason <- paste0(
    "Missing gtsummary/flextable package or baseline variables: ",
    paste(setdiff(table1_vars, names(association_data)), collapse = ", ")
  )
}


# -----------------------------------------------------------------------------
# ADD-ON A. RWDID restricted cubic spline (RCS) association curve
#
# This is an additional descriptive/model-based analysis. It does not alter
# beta estimation, the four prespecified Cox models, score construction,
# tertile cutpoints, or any existing result.
# The 3 knots are estimated in the training cohort and exported for Step_04_2_2.
# -----------------------------------------------------------------------------

message("\n[ADD-ON 1/2] Fitting the original-style 3-knot RWDID RCS model ...")

if (!requireNamespace("rms", quietly = TRUE)) {
  stop_now("Package 'rms' is required for the additional RWDID RCS analysis.")
}

fit_rwdid_rcs_addon <- function(data, knots = NULL, reference_value,
                                cohort_label, file_prefix) {
  if (!is.data.frame(data) ||
      !all(c("depression_time", "depression_status", "diet_score", COVARS) %in%
           names(data))) {
    stop_now("RCS input data are missing the outcome, RWDID, or Model 3 covariates.")
  }
  if (nrow(data) < 100L ||
      sum(data$depression_status == 1, na.rm = TRUE) < 20L) {
    stop_now("Insufficient observations/events for the additional RCS model.")
  }
  if (!is.finite(reference_value)) {
    stop_now("The frozen RWDID RCS reference value is not finite.")
  }

  # Build datadist only from variables actually used in this RCS model.
  # This avoids warnings from unrelated constant audit/metadata columns.
  rcs_datadist_vars <- unique(c("diet_score", COVARS))
  ddist_data <- data[, rcs_datadist_vars, drop = FALSE]
  ddist <- rms::datadist(ddist_data)

  if (!"diet_score" %in% colnames(ddist$limits)) {
    stop_now("rms::datadist did not retain the RWDID variable.")
  }
  ddist$limits["Adjust to", "diet_score"] <- reference_value

  # rms::Design() resolves options(datadist = "<object name>") from an
  # environment on the search path. A function-local object is not sufficient.
  # Register a temporary, uniquely named object in .GlobalEnv and restore the
  # prior option/object on exit, including when model fitting raises an error.
  datadist_object_name <- ".rwdid_rcs_datadist_addon"
  had_previous_datadist_object <- exists(
    datadist_object_name, envir = .GlobalEnv, inherits = FALSE
  )
  if (had_previous_datadist_object) {
    previous_datadist_object <- get(
      datadist_object_name, envir = .GlobalEnv, inherits = FALSE
    )
  }

  assign(datadist_object_name, ddist, envir = .GlobalEnv)
  previous_datadist_option <- options(datadist = datadist_object_name)

  on.exit({
    options(previous_datadist_option)
    if (had_previous_datadist_object) {
      assign(
        datadist_object_name, previous_datadist_object,
        envir = .GlobalEnv
      )
    } else if (exists(datadist_object_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = datadist_object_name, envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (is.null(knots)) {
    rcs_term <- "rms::rcs(diet_score, 3)"
  } else {
    knots <- as.numeric(knots)
    if (length(knots) != 3L ||
        any(!is.finite(knots)) ||
        is.unsorted(knots, strictly = TRUE)) {
      stop_now("The frozen RCS knots must be three finite, strictly increasing values.")
    }
    knots_text <- paste(
      format(knots, digits = 17, scientific = FALSE, trim = TRUE),
      collapse = ", "
    )
    rcs_term <- paste0("rms::rcs(diet_score, c(", knots_text, "))")
  }

  rcs_formula <- stats::as.formula(
    paste0(
      "survival::Surv(depression_time, depression_status == 1) ~ ",
      rcs_term,
      " + ",
      paste(COVARS, collapse = " + ")
    )
  )

  fit <- rms::cph(
    rcs_formula,
    data = data,
    x = TRUE,
    y = TRUE
  )

  if (is.null(knots)) {
    if (is.null(fit$Design$parms$diet_score)) {
      stop_now("Could not recover the three RCS knots from the training rms model.")
    }
    knots <- as.numeric(fit$Design$parms$diet_score)
  }
  if (length(knots) != 3L ||
      any(!is.finite(knots)) ||
      is.unsorted(knots, strictly = TRUE)) {
    stop_now("The fitted RCS model did not yield three valid RWDID knots.")
  }

  anova_tab <- as.data.frame(stats::anova(fit))
  anova_rows <- trimws(rownames(anova_tab))
  p_col <- names(anova_tab)[tolower(names(anova_tab)) == "p"]
  if (length(p_col) != 1L) {
    p_col <- grep("p", names(anova_tab), ignore.case = TRUE, value = TRUE)[1]
  }
  if (is.na(p_col) || !nzchar(p_col)) {
    stop_now("Could not identify the P-value column in the rms ANOVA table.")
  }

  overall_idx <- which(anova_rows == "diet_score")
  if (!length(overall_idx)) {
    overall_idx <- which(
      grepl("diet_score", anova_rows, fixed = TRUE) &
        !grepl("Nonlinear", anova_rows, ignore.case = TRUE)
    )
  }
  nonlinear_idx <- which(grepl("Nonlinear", anova_rows, ignore.case = TRUE))

  if (!length(overall_idx)) {
    stop_now("Could not identify the overall RWDID term in the rms ANOVA table.")
  }

  p_overall <- as.numeric(anova_tab[overall_idx[1], p_col])
  p_nonlinear <- if (length(nonlinear_idx)) {
    as.numeric(anova_tab[nonlinear_idx[1], p_col])
  } else {
    NA_real_
  }

  pred <- as.data.frame(
    rms::Predict(
      fit,
      diet_score,
      fun = exp,
      type = "predictions",
      ref.zero = TRUE
    )
  )
  required_prediction_cols <- c("diet_score", "yhat", "lower", "upper")
  if (!all(required_prediction_cols %in% names(pred)) ||
      nrow(pred) < 10L ||
      any(!is.finite(pred$yhat)) ||
      any(!is.finite(pred$lower)) ||
      any(!is.finite(pred$upper))) {
    stop_now("The RWDID RCS prediction table is incomplete or non-finite.")
  }

  pred$cohort <- cohort_label
  write_csv_utf8(
    pred,
    paste0(file_prefix, "_curve_data.csv")
  )

  x_range <- range(pred$diet_score, na.rm = TRUE)
  y_upper <- max(pred$upper, na.rm = TRUE)
  y_lower <- min(pred$lower, na.rm = TRUE)
  x_annotate <- x_range[1] + 0.025 * diff(x_range)
  y_annotate <- y_upper * 1.13
  p_label <- paste0(
    "P for overall ", fmt_p(p_overall),
    "\nP for non-linear ", fmt_p(p_nonlinear)
  )

  friendly_colours <- cols4all::c4a("friendly13")
  rcs_plot <- ggplot2::ggplot(
    pred,
    ggplot2::aes(x = diet_score, y = yhat)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.10,
      fill = friendly_colours[1]
    ) +
    ggplot2::geom_line(
      linewidth = 1,
      alpha = 0.70,
      colour = friendly_colours[1]
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = 2,
      linewidth = 0.5,
      colour = friendly_colours[2],
      alpha = 0.50
    ) +
    ggplot2::annotate(
      "text",
      x = x_annotate,
      y = y_annotate,
      label = p_label,
      hjust = 0,
      vjust = 1,
      size = 4.5
    ) +
    ggplot2::scale_x_continuous(name = "RWDID") +
    ggplot2::scale_y_continuous(
      name = "Hazard ratio (95% CI)",
      expand = ggplot2::expansion(mult = c(0.03, 0.18))
    ) +
    ggplot2::coord_cartesian(
      ylim = c(max(0, y_lower * 0.90), y_annotate * 1.03)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(
        colour = "black", size = 18, face = "bold"
      ),
      axis.title = ggplot2::element_text(
        colour = "black", size = 18, face = "bold"
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 18
      ),
      axis.title.x = ggplot2::element_text(
        colour = "black", size = 18, face = "bold"
      ),
      axis.title.y = ggplot2::element_text(
        colour = "black", size = 18, face = "bold"
      )
    ) +
    ggplot2::ggtitle(cohort_label)

  ggplot2::ggsave(
    filename = paste0(file_prefix, ".pdf"),
    plot = rcs_plot,
    width = 10,
    height = 6,
    device = "pdf"
  )
  ggplot2::ggsave(
    filename = paste0(file_prefix, "_600dpi.png"),
    plot = rcs_plot,
    width = 10,
    height = 6,
    dpi = 600,
    bg = "white"
  )

  model_summary <- tibble(
    cohort = cohort_label,
    n = nrow(data),
    events = sum(data$depression_status == 1, na.rm = TRUE),
    reference_value_training_median = reference_value,
    knot_1 = knots[1],
    knot_2 = knots[2],
    knot_3 = knots[3],
    p_overall = p_overall,
    p_nonlinear = p_nonlinear
  )
  write_csv_utf8(
    model_summary,
    paste0(file_prefix, "_model_summary.csv")
  )

  list(
    fit = fit,
    prediction = pred,
    plot = rcs_plot,
    knots = knots,
    reference_value = reference_value,
    summary = model_summary
  )
}

RWDID_RCS_REFERENCE_TRAINING <- stats::median(
  association_data$diet_score,
  na.rm = TRUE
)
rcs_training_addon <- fit_rwdid_rcs_addon(
  data = association_data,
  knots = NULL,
  reference_value = RWDID_RCS_REFERENCE_TRAINING,
  cohort_label = "UK Biobank - Training cohort",
  file_prefix = file.path(PLOT_DIR, "RWDID_RCS_training")
)

RWDID_RCS_KNOTS_TRAINING <- rcs_training_addon$knots
RCS_KNOTS_FILE <- file.path(OUT_DIR, "RWDID_RCS_knots_training.csv")
write_csv_utf8(
  tibble(
    knot_order = seq_along(RWDID_RCS_KNOTS_TRAINING),
    knot_value = RWDID_RCS_KNOTS_TRAINING,
    reference_value_training = RWDID_RCS_REFERENCE_TRAINING
  ),
  RCS_KNOTS_FILE
)
saveRDS(
  rcs_training_addon$fit,
  file.path(OUT_DIR, "RWDID_RCS_training_model.rds")
)

# -----------------------------------------------------------------------------
# ADD-ON B. Publication-ready Table 1 by the three RWDID tertiles
#
# The summary is restricted to the same complete-case association cohort used
# by the prespecified Cox models. Follow-up time is stored in days and is
# converted to years using 365.25 days/year for display only.
# -----------------------------------------------------------------------------

message("\n[ADD-ON 2/2] Creating the publication-ready training Table 1 ...")

if (!requireNamespace("gtsummary", quietly = TRUE) ||
    !requireNamespace("flextable", quietly = TRUE)) {
  stop_now(
    "Packages 'gtsummary' and 'flextable' are required for the additional ",
    "publication-ready Table 1."
  )
}

make_rwdid_publication_table1 <- function(data, output_prefix, cohort_note) {
  table1_vars <- c(
    "age_food", "follow_up_years", "sex", "ethnicity",
    "education_attainment", "career", "twnsend", "smoking_status",
    "drink_status", "bmi_group", "household_income",
    "physical_activity_proup", "sedentary_group", "energy_intake",
    "central_obesity", "triglycerides", "hdl_cholesterol",
    "hypertension_status_final", "t2diabetes_status_final",
    "dyslipidemia_status_final", "cvd_status_final",
    "all_cancer_status_final"
  )
  continuous_vars <- c(
    "age_food", "follow_up_years", "education_attainment",
    "twnsend", "energy_intake"
  )
  binary_vars <- c(
    "central_obesity", "triglycerides", "hdl_cholesterol",
    "hypertension_status_final", "t2diabetes_status_final",
    "dyslipidemia_status_final", "cvd_status_final",
    "all_cancer_status_final"
  )

  required <- c(
    "depression_time", "diet_score", "diet_score_tertile",
    table1_vars[table1_vars != "follow_up_years"]
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop_now(
      "Table 1 source data are missing variables: ",
      paste(missing, collapse = ", ")
    )
  }

  table_data <- as.data.frame(data)
  table_data$follow_up_years <- as.numeric(table_data$depression_time) / 365.25
  if (any(!is.finite(table_data$follow_up_years)) ||
      any(table_data$follow_up_years <= 0)) {
    stop_now("Derived follow-up years contain invalid values.")
  }

  table_data$diet_score_tertile <- factor(
    as.character(table_data$diet_score_tertile),
    levels = c("T1", "T2", "T3")
  )
  if (anyNA(table_data$diet_score_tertile)) {
    stop_now("Table 1 contains missing or invalid RWDID tertile assignments.")
  }

  order_factor <- function(x, preferred_levels) {
    z <- trimws(as.character(x))
    observed <- unique(z[!is.na(z) & nzchar(z)])
    factor(
      z,
      levels = c(
        preferred_levels[preferred_levels %in% observed],
        setdiff(observed, preferred_levels)
      )
    )
  }

  table_data$sex <- order_factor(
    table_data$sex, c("Male", "Female")
  )
  table_data$ethnicity <- order_factor(
    table_data$ethnicity, c("Non-white", "White")
  )
  table_data$career <- order_factor(
    table_data$career, c("Employed", "Retired/other")
  )
  table_data$smoking_status <- order_factor(
    table_data$smoking_status, c("Current", "Never", "Previous")
  )
  table_data$drink_status <- order_factor(
    table_data$drink_status, c("Current", "Never", "Previous")
  )

  bmi_values <- trimws(as.character(table_data$bmi_group))
  bmi_values[bmi_values == "Obese (>=30)"] <- "Obese (≥30)"
  table_data$bmi_group <- order_factor(
    bmi_values,
    c(
      "Underweight (<18.5)",
      "Normal weight (18.5 to <25)",
      "Overweight (25 to <30)",
      "Obese (≥30)"
    )
  )
  table_data$household_income <- order_factor(
    table_data$household_income, c("Low", "Middle", "High")
  )
  table_data$physical_activity_proup <- order_factor(
    table_data$physical_activity_proup,
    c("High physical activity", "Low physical activity")
  )
  table_data$sedentary_group <- order_factor(
    table_data$sedentary_group, c("Not sedentary", "Sedentary")
  )

  # Convert conventional 0/1 and Yes/No indicators to explicit No/Yes factors.
  # Unexpected source labels are retained as categorical levels rather than
  # silently recoded or discarded.
  normalize_binary <- function(x) {
    z <- trimws(as.character(x))
    z_lower <- tolower(z)
    map <- c(
      "0" = "No", "1" = "Yes",
      "no" = "No", "yes" = "Yes",
      "false" = "No", "true" = "Yes"
    )
    observed <- unique(z_lower[!is.na(z_lower) & nzchar(z_lower)])
    if (all(observed %in% names(map))) {
      factor(unname(map[z_lower]), levels = c("No", "Yes"))
    } else {
      factor(z, levels = unique(z[!is.na(z)]))
    }
  }
  for (v in binary_vars) {
    table_data[[v]] <- normalize_binary(table_data[[v]])
  }

  for (v in continuous_vars) {
    original <- table_data[[v]]
    converted <- suppressWarnings(as.numeric(as.character(original)))
    if (any(!is.na(original) & is.na(converted))) {
      stop_now("A Table 1 continuous variable is not numeric: ", v)
    }
    table_data[[v]] <- converted
  }

  dichotomous_vars <- binary_vars[vapply(
    binary_vars,
    function(v) {
      observed <- unique(as.character(table_data[[v]][!is.na(table_data[[v]])]))
      length(observed) > 0L && all(observed %in% c("No", "Yes"))
    },
    logical(1)
  )]

  type_spec <- c(
    list(
      c(
        age_food, follow_up_years, education_attainment,
        twnsend, energy_intake
      ) ~ "continuous"
    ),
    lapply(
      dichotomous_vars,
      function(v) stats::as.formula(paste0(v, " ~ 'dichotomous'"))
    )
  )
  value_spec <- lapply(
    dichotomous_vars,
    function(v) stats::as.formula(paste0(v, " ~ 'Yes'"))
  )

  table1_labels <- list(
    age_food ~ "Age, year, median (IQR)",
    follow_up_years ~ "Follow-up, years, median (IQR)",
    sex ~ "Sex, no. (%)",
    ethnicity ~ "Ethnicity, no. (%)",
    education_attainment ~ "Education attainment, mean (SD)",
    career ~ "Employment status, no. (%)",
    twnsend ~ "Townsend deprivation index, mean (SD)",
    smoking_status ~ "Smoking status, no. (%)",
    drink_status ~ "Drinking status, no. (%)",
    bmi_group ~ "BMI category (kg/m2), no. (%)",
    household_income ~ "Family income, no. (%)",
    physical_activity_proup ~ "Physical activity, no. (%)",
    sedentary_group ~ "Sedentary behavior, no. (%)",
    energy_intake ~ "Energy intake, KJ/day, mean (SD)",
    central_obesity ~ "Central obesity, no. (%)",
    triglycerides ~ "High triglycerides, no. (%)",
    hdl_cholesterol ~ "Low HDL, no. (%)",
    hypertension_status_final ~ "Hypertension, no. (%)",
    t2diabetes_status_final ~ "Diabetes, no. (%)",
    dyslipidemia_status_final ~ "Dyslipidemia, no. (%)",
    cvd_status_final ~ "CVD, no. (%)",
    all_cancer_status_final ~ "Cancer, no. (%)"
  )

  tbl <- gtsummary::tbl_summary(
    data = table_data,
    by = all_of("diet_score_tertile"),
    include = all_of(table1_vars),
    label = table1_labels,
    type = type_spec,
    value = value_spec,
    statistic = list(
      age_food ~ "{median} ({p25}, {p75})",
      follow_up_years ~ "{median} ({p25}, {p75})",
      education_attainment ~ "{mean} ({sd})",
      twnsend ~ "{mean} ({sd})",
      energy_intake ~ "{mean} ({sd})",
      gtsummary::all_categorical() ~ "{n} ({p}%)",
      gtsummary::all_dichotomous() ~ "{n} ({p}%)"
    ),
    digits = list(
      age_food ~ 0,
      follow_up_years ~ 2,
      education_attainment ~ 2,
      twnsend ~ 2,
      energy_intake ~ 0,
      gtsummary::all_categorical() ~ c(0, 1),
      gtsummary::all_dichotomous() ~ c(0, 1)
    ),
    missing = "no"
  ) %>%
    gtsummary::add_overall(last = FALSE)

  group_n <- table_data %>%
    dplyr::count(diet_score_tertile, .drop = FALSE)
  if (nrow(group_n) != 3L || any(group_n$n <= 0L)) {
    stop_now("The publication Table 1 requires non-empty T1, T2, and T3 groups.")
  }

  tb <- as.data.frame(tbl$table_body)
  if (!all(c("variable", "row_type", "label") %in% names(tb))) {
    stop_now("Unexpected gtsummary table_body structure.")
  }
  stat_cols <- grep("^stat_[0-9]+$", names(tb), value = TRUE)
  if (!identical(stat_cols, c("stat_0", "stat_1", "stat_2", "stat_3"))) {
    stop_now(
      "Expected Total plus three RWDID tertile columns in Table 1; found: ",
      paste(stat_cols, collapse = ", ")
    )
  }

  tb$variable <- as.character(tb$variable)
  tb$row_type <- as.character(tb$row_type)
  tb$label <- as.character(tb$label)
  for (v in stat_cols) tb[[v]] <- as.character(tb[[v]])

  insert_section_row <- function(x, before_variable, section_label) {
    idx <- which(x$variable == before_variable)[1]
    if (is.na(idx)) {
      stop_now("Could not insert Table 1 section header before ", before_variable)
    }
    new_row <- x[NA_integer_, , drop = FALSE]
    new_row$variable <- paste0(".section_", make.names(section_label))
    new_row$row_type <- "label"
    new_row$label <- section_label
    for (v in stat_cols) new_row[[v]] <- ""
    dplyr::bind_rows(
      x[seq_len(idx - 1L), , drop = FALSE],
      new_row,
      x[seq.int(idx, nrow(x)), , drop = FALSE]
    )
  }

  tb <- insert_section_row(
    tb, "central_obesity", "Metabolic syndrome"
  )
  tb <- insert_section_row(
    tb, "dyslipidemia_status_final", "Other diseases"
  )

  table_export <- tb %>%
    dplyr::select(
      Characteristic = label,
      dplyr::all_of(stat_cols)
    )

  fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
  total_label <- paste0("Total (n = ", fmt_n(nrow(table_data)), ")")
  t1_label <- paste0("T1 (n = ", fmt_n(group_n$n[group_n$diet_score_tertile == "T1"]), ")")
  t2_label <- paste0("T2 (n = ", fmt_n(group_n$n[group_n$diet_score_tertile == "T2"]), ")")
  t3_label <- paste0("T3 (n = ", fmt_n(group_n$n[group_n$diet_score_tertile == "T3"]), ")")

  ft <- flextable::flextable(table_export)
  ft <- flextable::set_header_labels(
    ft,
    Characteristic = "Characteristic",
    stat_0 = total_label,
    stat_1 = t1_label,
    stat_2 = t2_label,
    stat_3 = t3_label
  )
  ft <- flextable::add_header_row(
    ft,
    values = c("Characteristic", total_label, "RWDID"),
    colwidths = c(1, 1, 3),
    top = TRUE
  )
  ft <- flextable::merge_v(
    ft,
    j = c("Characteristic", "stat_0"),
    part = "header"
  )
  section_rows <- which(
    table_export$Characteristic %in% c("Metabolic syndrome", "Other diseases")
  )
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
  ft <- flextable::fontsize(ft, size = 8.5, part = "all")
  ft <- flextable::bold(ft, part = "header", bold = TRUE)
  if (length(section_rows)) {
    ft <- flextable::bold(
      ft, i = section_rows, j = "Characteristic",
      bold = TRUE, part = "body"
    )
    ft <- flextable::bg(
      ft, i = section_rows, bg = "#F2F2F2", part = "body"
    )
  }
  ft <- flextable::align(ft, j = "Characteristic", align = "left", part = "all")
  ft <- flextable::align(ft, j = stat_cols, align = "center", part = "all")
  ft <- flextable::valign(ft, valign = "center", part = "all")
  ft <- flextable::padding(ft, padding = 2, part = "all")
  ft <- flextable::width(ft, j = "Characteristic", width = 2.25)
  ft <- flextable::width(ft, j = stat_cols, width = 1.05)
  ft <- flextable::set_table_properties(ft, layout = "fixed", width = 1)
  ft <- flextable::add_footer_lines(
    ft,
    values = paste0(
      "Continuous variables are presented as median (IQR) or mean (SD), ",
      "as labelled; categorical variables are presented as n (%). ",
      "Percentages are column percentages. ", cohort_note
    )
  )

  csv_out <- table_export
  names(csv_out) <- c(
    "Characteristic", total_label, t1_label, t2_label, t3_label
  )
  utils::write.csv(
    csv_out,
    file = paste0(output_prefix, ".csv"),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
  flextable::save_as_docx(
    "Table 1" = ft,
    path = paste0(output_prefix, ".docx")
  )
  saveRDS(
    tbl,
    file = paste0(output_prefix, "_gtsummary.rds")
  )

  list(
    table = tbl,
    flextable = ft,
    n = nrow(table_data),
    tertile_n = group_n,
    output_docx = paste0(output_prefix, ".docx"),
    output_csv = paste0(output_prefix, ".csv")
  )
}

table1_training_addon <- make_rwdid_publication_table1(
  data = association_data,
  output_prefix = file.path(
    OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_training"
  ),
  cohort_note = "Training cohort; RWDID tertile cutpoints were estimated in the training cohort."
)

addon_qc_training <- tibble(
  check = c(
    "training_RWDID_RCS_model_fitted",
    "training_RWDID_RCS_knots_exported",
    "training_RWDID_RCS_curve_pdf_exported",
    "training_RWDID_RCS_curve_600dpi_png_exported",
    "training_RWDID_RCS_prediction_data_exported",
    "training_RWDID_Table1_word_exported",
    "training_RWDID_Table1_csv_exported",
    "training_RWDID_Table1_contains_Total_T1_T2_T3"
  ),
  status = c(
    inherits(rcs_training_addon$fit, "cph"),
    file.exists(RCS_KNOTS_FILE),
    file.exists(file.path(PLOT_DIR, "RWDID_RCS_training.pdf")),
    file.exists(file.path(PLOT_DIR, "RWDID_RCS_training_600dpi.png")),
    file.exists(file.path(PLOT_DIR, "RWDID_RCS_training_curve_data.csv")),
    file.exists(file.path(OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_training.docx")),
    file.exists(file.path(OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_training.csv")),
    nrow(table1_training_addon$tertile_n) == 3L
  ),
  detail = c(
    paste0("N=", nrow(association_data),
           "; events=", sum(association_data$depression_status == 1)),
    RCS_KNOTS_FILE,
    file.path(PLOT_DIR, "RWDID_RCS_training.pdf"),
    file.path(PLOT_DIR, "RWDID_RCS_training_600dpi.png"),
    file.path(PLOT_DIR, "RWDID_RCS_training_curve_data.csv"),
    table1_training_addon$output_docx,
    table1_training_addon$output_csv,
    paste0(
      "Total N=", table1_training_addon$n,
      "; T1/T2/T3 N=",
      paste(table1_training_addon$tertile_n$n, collapse = "/")
    )
  )
)
write_csv_utf8(
  addon_qc_training,
  file.path(OUT_DIR, "RCS_Table1_addon_QC.csv")
)
if (any(!addon_qc_training$status)) {
  print(addon_qc_training %>% dplyr::filter(!status))
  stop_now("The additional training RCS/Table 1 QC did not pass.")
}


# Save compact analysis data and model objects for exact reproducibility.
analysis_data_export <- association_data %>%
  select(
    all_of(ID_COL),
    depression_time,
    depression_status,
    all_of(ALL_FOOD_VARS),
    all_of(COVARS),
    n_components_scored,
    n_V5_score_components,
    n_primary_score_components,
    RWDID_total_score_complete_primary,
    RWDID_total_score_complete15,
    weighted_sum,
    diet_score,
    diet_score_tertile,
    diet_score_tertile_num,
    diet_score_20
  )

saveRDS(
  analysis_data_export,
  file.path(OUT_DIR, "diet_score_analysis_data.rds")
)

saveRDS(
  list(
    beta_model = fit_beta,
    continuous_model = fit_continuous,
    tertile_model = fit_tertile,
    trend_model = fit_trend,
    per20_model = fit_per20,
    beta_results = beta_results,
    beta_vector = beta_vector,
    beta_sum = beta_sum,
    tertile_cutpoints = c(T1_upper = tertile_q1, T2_upper = tertile_q2),
    per20_increment = range_20,
    score_mean = score_mean,
    score_min = score_min,
    score_max = score_max,
    id_column = ID_COL,
    food_map = FOOD_MAP,
    active_food_map = FOOD_MAP,
    rcs_range_manifest = rcs_ranges %>%
      select(
        food, value1, value2, range_direction, range_source,
        score_lower, score_upper, score_rule
      ),
    v5_binary_scoring_rules = v5_rules,
    covariates = COVARS
  ),
  file.path(OUT_DIR, "diet_score_models.rds")
)

# -----------------------------------------------------------------------------
# 11. Final QC gate and run summary
# -----------------------------------------------------------------------------

message("\n[7/7] Running final V5 QC gate ...")

final_qc <- tibble(
  check = c(
    "frozen_manifest_seed_matches_20261103",
    "frozen_manifest_training_fraction_matches_70_percent",
    "frozen_manifest_training_n_matches_76578",
    "frozen_manifest_test_n_matches_32819",
    "training_RDS_exactly_matches_manifest_training_IDs",
    "test_participants_excluded_from_training_RDS",
    "RCS_interval_file_contains_exact_15_components",
    "binary_rule_file_contains_exact_15_components",
    "rule_file_matches_V5_interval_file",
    "score_file_exactly_matches_training_IDs",
    "original_V5_3_score_summary_matches_source_score_matrix",
    "optional_repair_audit_matches_expected_intersections_if_supplied",
    "all_13_unchanged_components_match_original_V5_3_scores",
    "n_components_scored_matches_recomputed_scores",
    "complete15_total_matches_recomputed_scores",
    "all_component_scores_are_0_1_or_NA",
    "all_15_components_in_mutually_adjusted_beta_model",
    "all_15_beta_coefficients_are_finite",
    "all_15_components_have_score_variation",
    "beta_signs_are_reported_without_forcing_a_direction",
    "original_sum_beta_denominator_is_numerically_stable",
    "beta_weighted_score_only_calculated_for_complete_15_component_rows",
    "beta_model_and_association_model_use_same_training_participants",
    "training_tertiles_have_three_nonempty_groups",
    "per20_training_increment_is_positive_and_finite",
    "continuous_tertile_trend_and_per20_models_fitted",
    "forestplot_file_exists",
    "beta_heatmap_file_exists",
    "beta_coefficients_file_exists",
    "participant_beta_weighted_score_file_exists",
    "scaling_parameters_file_exists",
    "component_manifest_file_exists",
    "analysis_data_RDS_exists",
    "model_RDS_exists",
    "interval_export_exists",
    "rule_export_exists",
    "participant_score_export_exists"
  ),
  status = c(
    unique(split_manifest$split_seed) == EXPECTED_SPLIT_SEED,
    abs(unique(split_manifest$training_fraction) - EXPECTED_TRAIN_PROP) <= 1e-12,
    length(training_ids) == EXPECTED_TRAIN_N,
    length(test_ids) == EXPECTED_TEST_N,
    identical(training_ids_observed, training_ids),
    length(intersect(training_ids_observed, test_ids)) == 0L,
    nrow(rcs_ranges) == 15L &&
      identical(as.character(rcs_ranges$food), FOOD_MAP$food),
    nrow(v5_rules) == 15L &&
      identical(as.character(v5_rules$food), FOOD_MAP$food),
    source_rule_interval_match,
    identical(rcs_score_ids, training_ids),
    source_score_summary_match,
    repair_audit_match,
    nrow(unchanged_reconciliation) == 13L &&
      all(unchanged_reconciliation$exact_match_to_RCS_frozen_score),
    new_component_count_match,
    all(same_total),
    all(vapply(
      SCORE_COLS,
      function(v) all(is.na(rcs_score_aligned[[v]]) |
                        rcs_score_aligned[[v]] %in% c(0, 1)),
      logical(1)
    )),
    nrow(beta_results) == 15L &&
      setequal(beta_results$term, FOOD_VARS),
    nrow(beta_results) == 15L &&
      all(is.finite(beta_results$estimate)),
    all(component_variation$n_zero > 0L &
          component_variation$n_one > 0L),
    file.exists(file.path(OUT_DIR, "beta_sign_audit.csv")) &&
      nrow(beta_sign_audit) == 15L,
    is.finite(beta_sum) && abs(beta_sum) > beta_sum_tolerance,
    all(is.na(analysis_data$diet_score[!complete_primary])) &&
      all(is.finite(analysis_data$diet_score[complete_primary])),
    identical(
      sort(association_data[[ID_COL]], method = "radix"),
      sort(beta_model_data[[ID_COL]], method = "radix")
    ),
    nrow(tertile_counts) == 3L && all(tertile_counts$total > 0L),
    is.finite(range_20) && range_20 > 0,
    all(c(
      "diet_score", "diet_score_tertile", "diet_score_tertile_num",
      "diet_score_20"
    ) %in% names(association_data)),
    file.exists(file.path(OUT_DIR, "diet_score_forestplot.pdf")),
    file.exists(file.path(OUT_DIR, "beta_results.pdf")),
    file.exists(file.path(OUT_DIR, "beta_coefficients_training.csv")),
    file.exists(file.path(OUT_DIR, "beta_weighted_diet_score_training.csv")),
    file.exists(file.path(OUT_DIR, "diet_score_scaling_parameters.csv")),
    file.exists(file.path(OUT_DIR, "component_manifest_used.csv")),
    file.exists(file.path(OUT_DIR, "diet_score_analysis_data.rds")),
    file.exists(file.path(OUT_DIR, "diet_score_models.rds")),
    file.exists(file.path(OUT_DIR, "food_rcs_range_all.csv")),
    file.exists(file.path(OUT_DIR, "RWDID_binary_scoring_rules_15_components.csv")),
    file.exists(file.path(OUT_DIR, "RWDID_binary_component_scores_training.csv"))
  ),
  detail = c(
    paste0("seed=", unique(split_manifest$split_seed)),
    paste0("training_fraction=", unique(split_manifest$training_fraction)),
    paste0("training_n=", length(training_ids)),
    paste0("test_n=", length(test_ids)),
    paste0("training_RDS_n=", nrow(training_data)),
    paste0("training/test ID overlap=", length(intersect(training_ids_observed, test_ids))),
    paste0("Final Step_04_1_1 interval rows=", nrow(rcs_ranges)),
    paste0("Final Step_04_1_1 scoring-rule rows=", nrow(v5_rules)),
    "range direction, values, score bounds, and score-rule text match",
    paste0("Step_04_1_1 frozen score rows=", nrow(rcs_score_data)),
    "Step_04_1_1 score summary independently reconciled",
        repair_audit_detail,
    paste0(
      "unchanged component exact matches=",
      sum(unchanged_reconciliation$exact_match_to_RCS_frozen_score),
      "/13; direction-repaired components=",
      paste(REPAIR_FOODS, collapse = ", ")
    ),
    "n_components_scored independently reconciled",
    "complete15 total independently reconciled",
    "All recomputed component score values are 0/1/NA",
    paste0("beta model components=", length(FOOD_VARS), "/15"),
    paste0("estimated finite beta coefficients=", nrow(beta_results), "/15"),
    paste0(
      "components with both score values present=",
      sum(component_variation$n_zero > 0 & component_variation$n_one > 0),
      "/15"
    ),
    paste0(
      "negative=", sum(beta_vector < 0),
      "; zero=", sum(beta_vector == 0),
      "; positive=", sum(beta_vector > 0),
      "; no sign constraint applied"
    ),
    paste0("sum_beta=", format(beta_sum, digits = 12)),
    paste0(
      "complete_15=", sum(complete_primary),
      "; incomplete_15=", sum(!complete_primary)
    ),
    paste0(
      "beta model n=", nrow(beta_model_data),
      "; association n=", nrow(association_data)
    ),
    paste0("T1/T2/T3 n=", paste(tertile_counts$total, collapse = "/")),
    paste0("range_20=", format(range_20, digits = 12)),
    "All four requested association model specifications are present",
    "diet_score_forestplot.pdf",
    "beta_results.pdf",
    "beta_coefficients_training.csv",
    "beta_weighted_diet_score_training.csv",
    "diet_score_scaling_parameters.csv",
    "component_manifest_used.csv",
    "diet_score_analysis_data.rds",
    "diet_score_models.rds",
    "food_rcs_range_all.csv",
    "RWDID_binary_scoring_rules_15_components.csv",
    "RWDID_binary_component_scores_training.csv"
  )
)

write_csv_utf8(final_qc, file.path(OUT_DIR, "analysis_QC.csv"))
if (any(!final_qc$status)) {
  print(final_qc %>% filter(!status))
  stop_now("STOP: One or more final beta-weighted score QC checks failed.")
}

# README for the corresponding held-out validation analysis.
readme_lines <- c(
  "STEP_04_2_1 TRAINING-SET BETA-WEIGHTED RWDID ANALYSIS - V5.3.1",
  "",
  "Cohort and frozen inputs:",
  paste0("- Split manifest: ", SPLIT_MANIFEST_FILE),
  paste0("- Frozen training RDS: ", TRAIN_RDS),
  paste0("- Frozen V5.3 RCS20 intervals: ", RCS_RANGE_FILE),
  paste0("- Frozen RCS diagnostic file: ", RCS_DIAG_FILE),
  paste0("- Frozen V5.3 binary scoring rules: ", RCS_RULE_FILE),
  paste0("- Original V5.3 participant-level training scores: ", RCS_SCORE_FILE),
  paste0("- Split seed: ", EXPECTED_SPLIT_SEED),
  paste0("- Training participants: ", nrow(training_data)),
  paste0("- Held-out test participants: ", length(test_ids)),
  "",
  "Scoring and beta estimation:",
  "- All 15 audited food components are included in the binary RWDID.",
  "- No additional V6.1 shape-specific eligibility gate is applied; the 15",
  "  component set is retained, with only the two documented intervals revised.",
  "- Green leafy vegetables and high-fat dairy are excluded from the frozen 15.",
  "- Fruit juices and water are included.",
  "- The 13 non-target components retain their V5.3 RCS20 intervals and scores.",
  "- Final intervals and binary component scores are imported from Step_04_1_1.",
  "- The two-component interval-derivation audit is supplementary and optional;",
  "  the final range/rule/participant-score files are the authoritative inputs.",
  "- Potatoes and red meats were already repaired in Step_04_1_1 by the",
  "  frozen RCS20 20% plateau x Model 2 lowest-risk-category intersection.",
  "- This beta step does not alter interval boundaries or component scores.",
  "- All 15 imported binary scores are independently recomputed from raw intake",
  "  and checked against the frozen Step_04_1_1 participant-level score export.",
  "- Middle intervals are inclusive at both boundaries; low intervals use",
  "  intake <= upper boundary; high intervals use intake >= lower boundary",
  "  with an open upper tail; point intervals use the V5 numeric tolerance.",
  "- Values outside the defined score-1 interval receive 0; missing intake",
  "  remains NA and is not imputed.",
  "- Beta coefficients are estimated for all 15 binary components in one",
  "  mutually adjusted Model 3 Cox model in the frozen 70% training cohort.",
  "- Beta signs are estimated and reported without sign constraints.",
  "- FDR-BH and Bonferroni-adjusted P values are reported across the 15 betas.",
  "- The original formula is retained:",
  "  weighted_sum = sum(component_binary_score * beta_i)",
  "  diet_score = weighted_sum / sum(beta_i).",
  "- A participant receives a beta-weighted score only when all 15 component",
  "  scores are observed; no missing component is imputed.",
  "",
  "Association analyses:",
  "- Continuous RWDID (per score increase).",
  "- Tertiles T1/T2/T3; T1 is the reference.",
  "- P for trend uses tertile codes 1/2/3.",
  "- Per 20% increment is 20% of the observed training score range.",
  "- All association models use the original 21 Model 3 covariates.",
  "- The forestplot layout, markers, axis ticks, and dimensions retain the",
  "  uploaded V1.4 script's established style.",
  "- The beta heatmap uses a diverging palette: green indicates negative beta",
  "  and magenta indicates positive beta.",
  paste0("- Beta sign pattern: ", beta_sign_pattern, "."),
  paste0("- Score interpretation: ", diet_score_direction_note),
  "",
  "Frozen parameters for Step_04_2_2:",
  "- beta_coefficients_training.csv",
  "- diet_score_scaling_parameters.csv",
  "- V5_RCS_interval_manifest_used.csv",
  "- V5_binary_scoring_rules_used.csv",
  "- These parameters must be applied unchanged to the held-out test cohort.",
  "- Do not re-estimate beta coefficients, tertile cutpoints, or the per-20%",
  "  increment in Step_04_2_2.",
  "",
  paste0("Baseline Table 1 status: ", table1_status),
  if (nzchar(table1_reason)) paste0("Baseline Table 1 note: ", table1_reason) else "",
  "",
  "Core outputs:",
  "- food_rcs_range_all.csv",
  "- RWDID_binary_scoring_rules_15_components.csv",
  "- RWDID_binary_component_scores_training.csv",
  "- beta_coefficients_training.csv",
  "- beta_sign_audit.csv",
  "- V5_component_manifest_used.csv (intervals)",
  "- beta_weight_audit.csv",
  "- beta_weighted_diet_score_training.csv",
  "- diet_score_association_results.csv",
  "- diet_score_forestplot.pdf",
  "- beta_results.pdf",
  "- diet_score_scaling_parameters.csv",
  "- diet_score_distribution.csv",
  "- diet_score_tertile_distribution.csv",
  "- binary_component_score_distribution.csv",
  "- cox_ph_diagnostics.csv",
  "- diet_score_analysis_data.rds",
  "- diet_score_models.rds",
  "- analysis_QC.csv"
)
writeLines(
  readme_lines[nzchar(readme_lines)],
  con = file.path(OUT_DIR, "README.txt"),
  useBytes = TRUE
)

cat("\n============================================================\n")
cat("STEP_04_2_1 TRAINING BETA-WEIGHTED RWDID DIRECTION-REPAIRED COMPLETED\n")
cat("Cohort: exact frozen Step_03_1_0 70% training RDS\n")
cat("Training N: ", nrow(training_data), "\n", sep = "")
cat("Beta-model N: ", nrow(beta_model_data), "\n", sep = "")
cat(
  "Beta-model events: ",
  sum(beta_model_data$depression_status == 1, na.rm = TRUE),
  "\n",
  sep = ""
)
cat("components in score: 15 (two repaired; 13 unchanged)\n")
cat("Beta coefficients: ", length(beta_vector),
    ", mutually adjusted Model 3 Cox, training only\n", sep = "")
cat("Sum beta denominator: ", format(beta_sum, digits = 12), "\n", sep = "")
cat("Beta sign pattern: ", beta_sign_pattern, "\n", sep = "")
cat("Score interpretation: ", diet_score_direction_note, "\n", sep = "")
cat("Complete 15-component scores: ", sum(complete_primary), "\n", sep = "")
cat("Association-model N: ", nrow(association_data), "\n", sep = "")
cat("Association-model events: ",
    sum(association_data$depression_status == 1, na.rm = TRUE), "\n", sep = "")
cat("Analyses: continuous / tertiles / P for trend / per 20% increment\n")
cat("Validation parameters frozen: beta / tertile cutpoints / per20 scale\n")
cat("Forest plot: ", file.path(OUT_DIR, "diet_score_forestplot.pdf"), "\n", sep = "")
cat("Beta heatmap: ", file.path(OUT_DIR, "beta_results.pdf"), "\n", sep = "")
cat("QC: ", file.path(OUT_DIR, "analysis_QC.csv"), "\n", sep = "")
cat("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 10/11: Step_04_2_1_BETA_WEIGHTED_TRAINING
########################################################################################


########################################################################################
# BEGIN ANALYSIS BLOCK 11/11: Step_04_2_2_HELDOUT_VALIDATION
# Source: Step_04_2_2_beta_weighted_diet_score_test_validation_V1_4_RCS_DATADIST_FIX(6).R
# Provenance: Library source dated 2026-11-16; current Step_04_2_2 V1.4 RCS/DATADIST fix
# SHA256 of included block: 09e36e0da1323add9cacf9c2221440eaa0ba67952f517c356f0da5980bd7f922
########################################################################################

# =============================================================================
# Step_04_2_2_RWDID_HELD_OUT_TEST_SET_VALIDATION_V1_1.R
#
# HELD-OUT 30% TEST-SET VALIDATION OF THE FROZEN RWDID
#
# Design is intentionally locked to the current Step_04_2_1 training outputs:
#   beta_coefficients_training.csv
#   diet_score_scaling_parameters.csv
#   RWDID_binary_scoring_rules_15_components.csv
#
# This script:
#   1) reads the exact frozen Step_03 30% test cohort and split manifest;
#   2) imports the 15 beta coefficients and scaling parameters exported by the
#      current Step_04_2_1 training script, plus its RCS rules and the source
#      ranges exported by Step_04_1_1_4;
#   3) applies those parameters unchanged to the held-out test cohort;
#   4) fits the same four Model 3 Cox association specifications as Step_04_2_1;
#   5) exports validation results, diagnostics, figures, provenance, and QC.
#
# It DOES NOT:
#   - create or modify a cohort split;
#   - estimate or update any beta coefficient;
#   - re-select food components or re-estimate RCS intervals;
#   - calculate new tertile cutpoints in the test set;
#   - re-estimate the per-20%-of-range increment or training centering value;
#   - impute missing component scores;
#   - clip test scores to the training score range.
#
# Run from the UK Biobank project root.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(broom)
  library(forestplot)
  library(cols4all)
  library(grid)
})

# -----------------------------------------------------------------------------
# 0. Frozen paths and cohort metadata
# -----------------------------------------------------------------------------

PROJECT_ROOT <- "."

SPLIT_DIR <- file.path(
  PROJECT_ROOT, "results", "Step_03", "ukb_70_30_split"
)
MAIN_SPLIT_DIR <- file.path(SPLIT_DIR, "main")
SPLIT_MANIFEST_FILE <- file.path(
  SPLIT_DIR, "ukb_70_30_split_manifest.csv"
)
TEST_RDS <- file.path(
  MAIN_SPLIT_DIR, "ukb_depression_test_cohort.rds"
)

TRAINING_OUTPUT_DIR <- file.path(
  PROJECT_ROOT, "results", "Step_04", "Step_04_2_1"
)

# Exact outputs written by the current Step_04_2_1 training script.
BETA_FILE <- file.path(
  TRAINING_OUTPUT_DIR, "beta_coefficients_training.csv"
)
SCALING_FILE <- file.path(
  TRAINING_OUTPUT_DIR, "diet_score_scaling_parameters.csv"
)
RCS_RULE_FILE <- file.path(
  TRAINING_OUTPUT_DIR, "RWDID_binary_scoring_rules_15_components.csv"
)
RCS_RANGE_FILE <- file.path(
  TRAINING_OUTPUT_DIR, "food_rcs_range_all.csv"
)
TRAIN_ASSOC_FILE <- file.path(
  TRAINING_OUTPUT_DIR, "diet_score_association_results.csv"
)
TRAIN_ANALYSIS_RDS <- file.path(
  TRAINING_OUTPUT_DIR, "diet_score_analysis_data.rds"
)
RCS_KNOTS_FILE <- file.path(
  TRAINING_OUTPUT_DIR, "RWDID_RCS_knots_training.csv"
)

# Independent source-of-truth range exported by the final RCS/binary-scoring
# script. Step_04_2_1 imports and re-exports this file without changing ranges.
RCS_SOURCE_DIR <- file.path(
  PROJECT_ROOT, "results", "Step_04", "Step_04_1_1_4"
)
RCS_SOURCE_RANGE_FILE <- file.path(
  RCS_SOURCE_DIR, "food_rcs_range_all_final.csv"
)

OUT_DIR <- file.path(
  PROJECT_ROOT, "results", "Step_04", "Step_04_2_2"
)
PLOT_DIR <- file.path(OUT_DIR, "plots")
AUDIT_DIR <- file.path(OUT_DIR, "audit")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPECTED_SPLIT_SEED <- 20261103L
EXPECTED_TRAIN_PROP <- 0.70
EXPECTED_TRAIN_N <- 76578L
EXPECTED_TEST_N <- 32819L
EXPECTED_TOTAL_N <- EXPECTED_TRAIN_N + EXPECTED_TEST_N

# Same 21 Model 3 covariates as the current Step_04_2_1 training script.
MODEL3_COVARIATES <- c(
  "age_food",
  "sex",
  "ethnicity",
  "education_attainment",
  "twnsend",
  "career",
  "smoking_status",
  "drink_status",
  "bmi_group",
  "physical_activity_proup",
  "sedentary_group",
  "household_income",
  "energy_intake",
  "central_obesity",
  "triglycerides",
  "hdl_cholesterol",
  "t2diabetes_status_final",
  "hypertension_status_final",
  "dyslipidemia_status_final",
  "cvd_status_final",
  "all_cancer_status_final"
)

# Exact 15-component order from Step_04_2_1 current frozen training.
FOOD_MAP <- tibble(
  beta_order = seq_len(15L),
  group = c(
    "Citrus fruits",
    "Red/orange vegetables",
    "Other vegetables",
    "Potatoes",
    "Refined grains",
    "Eggs",
    "Low fat dairy",
    "Red meats",
    "Processed meats",
    "Fish and seafood",
    "Coffee",
    "Fruit juices",
    "Sweetened beverages",
    "Water",
    "Wine"
  ),
  food = c(
    "citrus_fruits",
    "red_orange_veg",
    "other_veg",
    "potatoes",
    "refined_grains",
    "eggs",
    "low_fat_dairy",
    "red_meats",
    "processed_meats",
    "fish_seafood",
    "coffee",
    "fruit_juices",
    "sweetened_bev",
    "water",
    "wine"
  ),
  display_name = c(
    "Citrus fruits",
    "Red and orange vegetables",
    "Other vegetables",
    "Potatoes",
    "Refined grains",
    "Eggs",
    "Low-fat dairy",
    "Red meats",
    "Processed meats",
    "Fish and seafood",
    "Coffee",
    "Fruit juices",
    "Sweetened beverages",
    "Water",
    "Wine"
  ),
  score_col = paste0(
    "score_",
    c(
      "citrus_fruits", "red_orange_veg", "other_veg", "potatoes",
      "refined_grains", "eggs", "low_fat_dairy", "red_meats",
      "processed_meats", "fish_seafood", "coffee", "fruit_juices",
      "sweetened_bev", "water", "wine"
    )
  ),
  model_var = paste0(
    c(
      "citrus_fruits", "red_orange_veg", "other_veg", "potatoes",
      "refined_grains", "eggs", "low_fat_dairy", "red_meats",
      "processed_meats", "fish_seafood", "coffee", "fruit_juices",
      "sweetened_bev", "water", "wine"
    ),
    "_cat"
  )
)

FOOD_VARS <- FOOD_MAP$model_var
COVARS <- MODEL3_COVARIATES

# -----------------------------------------------------------------------------
# 1. Utility functions
# -----------------------------------------------------------------------------

stop_now <- function(...) {
  stop(paste0(...), call. = FALSE)
}

write_csv_utf8 <- function(x, file) {
  utils::write.csv(
    x,
    file = file,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

fmt_p <- function(p) {
  if (length(p) != 1L || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  formatC(p, format = "f", digits = 3)
}

fmt_hr_ci <- function(hr, lo, hi) {
  if (any(!is.finite(c(hr, lo, hi)))) return(NA_character_)
  paste0(
    formatC(hr, format = "f", digits = 2),
    " (",
    formatC(lo, format = "f", digits = 2),
    ", ",
    formatC(hi, format = "f", digits = 2),
    ")"
  )
}

get_frozen_parameter <- function(parameter_table, parameter_name) {
  hit <- parameter_table[
    as.character(parameter_table$parameter) == parameter_name,
    ,
    drop = FALSE
  ]
  if (nrow(hit) != 1L) {
    stop_now(
      "Expected exactly one frozen scaling parameter '",
      parameter_name, "', found ", nrow(hit), "."
    )
  }
  value <- suppressWarnings(as.numeric(hit$value[1]))
  if (length(value) != 1L || !is.finite(value)) {
    stop_now("Frozen parameter is non-finite: ", parameter_name, ".")
  }
  value
}

score_one_component <- function(x, direction, value1, value2, food_name) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  score <- rep(NA_integer_, length(x))
  observed <- !is.na(x)

  if (any(!is.finite(x[observed]))) {
    stop_now(
      "Non-finite non-missing intake values found for ", food_name, "."
    )
  }
  if (any(x[observed] < -1e-8)) {
    stop_now(
      "Negative intake values found for ", food_name,
      ". Check the test-cohort source coding before scoring."
    )
  }
  if (!any(observed)) return(score)

  scale_value <- max(
    1,
    abs(value1),
    ifelse(is.finite(value2), abs(value2), 0)
  )
  tol <- max(1e-8, scale_value * 1e-8)

  if (direction == "middle") {
    if (!is.finite(value1) || !is.finite(value2) || value2 < value1) {
      stop_now("Invalid frozen middle interval for ", food_name, ".")
    }
    score[observed] <- as.integer(
      x[observed] >= value1 - tol &
        x[observed] <= value2 + tol
    )
  } else if (direction == "low") {
    if (!is.finite(value1)) {
      stop_now("Invalid frozen low threshold for ", food_name, ".")
    }
    score[observed] <- as.integer(x[observed] <= value1 + tol)
  } else if (direction == "high") {
    if (!is.finite(value1)) {
      stop_now("Invalid frozen high threshold for ", food_name, ".")
    }
    score[observed] <- as.integer(x[observed] >= value1 - tol)
  } else if (direction == "point") {
    if (!is.finite(value1)) {
      stop_now("Invalid frozen point threshold for ", food_name, ".")
    }
    score[observed] <- as.integer(
      abs(x[observed] - value1) <= tol
    )
  } else {
    stop_now(
      "Unsupported frozen RCS interval direction for ",
      food_name, ": ", direction
    )
  }

  score
}

make_survival_formula <- function(exposure, covars = COVARS) {
  stats::as.formula(
    paste0(
      "survival::Surv(depression_time, depression_status) ~ ",
      paste(c(exposure, covars), collapse = " + ")
    )
  )
}

MODEL_WARNINGS <- list()

fit_cox_model <- function(formula, data, model_name) {
  warning_messages <- character(0)

  fit <- withCallingHandlers(
    tryCatch(
      survival::coxph(
        formula,
        data = data,
        ties = "efron",
        x = TRUE,
        y = TRUE,
        model = TRUE,
        na.action = na.fail
      ),
      error = function(e) {
        stop_now(
          "Test-set Cox model failed (", model_name, "): ",
          conditionMessage(e)
        )
      }
    ),
    warning = function(w) {
      warning_messages <<- c(
        warning_messages,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )

  MODEL_WARNINGS[[model_name]] <<- warning_messages

  if (is.null(fit$coefficients) ||
      any(!is.finite(fit$coefficients)) ||
      anyNA(fit$coefficients)) {
    stop_now(
      "Test-set Cox model contains non-finite or non-estimable coefficients: ",
      model_name,
      ". Inspect covariate variation and event distribution."
    )
  }
  if (!is.finite(as.numeric(tail(fit$loglik, 1)))) {
    stop_now("Test-set Cox log-likelihood is non-finite: ", model_name, ".")
  }

  fit
}

extract_tidy_term <- function(fit, term) {
  if (!is.character(term) || length(term) != 1L ||
      is.na(term) || !nzchar(trimws(term))) {
    stop_now("extract_tidy_term requires one non-empty coefficient name.")
  }

  td <- broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE,
    conf.level = 0.95
  )

  if (!"term" %in% names(td)) {
    stop_now("The Cox tidy result has no 'term' column.")
  }

  # Exact match; .env$term avoids the V1_1/V1_2 dplyr scoping bug.
  out <- td %>% dplyr::filter(.data$term == .env$term)

  if (nrow(out) != 1L) {
    stop_now(
      "Expected exactly one coefficient for test-set term '", term,
      "', found ", nrow(out), ". Available terms: ",
      paste(as.character(td$term), collapse = ", ")
    )
  }
  out
}

make_result_row <- function(variable, term, fit, count_data,
                            reference = FALSE) {
  n_total <- nrow(count_data)
  n_events <- sum(count_data$depression_status == 1, na.rm = TRUE)
  events_total <- paste0(
    formatC(n_events, format = "d", big.mark = ","),
    "/",
    formatC(n_total, format = "d", big.mark = ",")
  )

  if (reference) {
    return(tibble(
      variable = variable,
      Events_Total = events_total,
      hr_ci = "Ref.",
      p_value = "Ref.",
      hazard_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value_numeric = NA_real_,
      term = NA_character_
    ))
  }

  z <- extract_tidy_term(fit, term)
  tibble(
    variable = variable,
    Events_Total = events_total,
    hr_ci = fmt_hr_ci(z$estimate[1], z$conf.low[1], z$conf.high[1]),
    p_value = fmt_p(z$p.value[1]),
    hazard_ratio = as.numeric(z$estimate[1]),
    ci_low = as.numeric(z$conf.low[1]),
    ci_high = as.numeric(z$conf.high[1]),
    p_value_numeric = as.numeric(z$p.value[1]),
    term = term
  )
}

safe_cox_zph <- function(fit, model_name) {
  z <- tryCatch(
    survival::cox.zph(fit, transform = "km"),
    error = function(e) e
  )

  if (inherits(z, "error")) {
    return(tibble(
      model = model_name,
      term = "NOT_ESTIMABLE",
      chisq = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      diagnostic_message = conditionMessage(z)
    ))
  }

  tab <- as.data.frame(z$table)
  tibble(
    model = model_name,
    term = rownames(tab),
    chisq = as.numeric(tab[, 1]),
    df = as.numeric(tab[, 2]),
    p_value = as.numeric(tab[, 3]),
    diagnostic_message = ""
  )
}

# -----------------------------------------------------------------------------
# 2. Validate the authoritative frozen 70/30 split and test RDS
# -----------------------------------------------------------------------------

message("\n[1/8] Validating the frozen Step_03_1_0 test cohort ...")

required_inputs <- c(
  SPLIT_MANIFEST_FILE,
  TEST_RDS,
  BETA_FILE,
  SCALING_FILE,
  RCS_RULE_FILE,
  RCS_RANGE_FILE,
  RCS_SOURCE_RANGE_FILE,
  RCS_KNOTS_FILE,
  TRAIN_ANALYSIS_RDS
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0L) {
  stop_now(
    "Required frozen input file(s) are missing:\n",
    paste(missing_inputs, collapse = "\n"),
    "\nRun Step_03_1_0 and current Step_04_2_1 training first."
  )
}

split_manifest <- read.csv(
  SPLIT_MANIFEST_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(
    participant_id = "character",
    split = "character",
    id_column = "character",
    split_seed = "integer",
    training_fraction = "numeric"
  )
)

required_manifest_cols <- c(
  "participant_id", "split", "id_column",
  "split_seed", "training_fraction"
)
if (!all(required_manifest_cols %in% names(split_manifest))) {
  stop_now("Unexpected Step_03_1_0 split-manifest structure.")
}

split_manifest$participant_id <- trimws(
  as.character(split_manifest$participant_id)
)
split_manifest$split <- as.character(split_manifest$split)

if (anyNA(split_manifest$participant_id) ||
    any(split_manifest$participant_id == "") ||
    anyDuplicated(split_manifest$participant_id) > 0L) {
  stop_now("Invalid or duplicated participant IDs in the frozen manifest.")
}
if (!all(split_manifest$split %in% c("training", "test"))) {
  stop_now("Invalid split labels in the frozen manifest.")
}
if (length(unique(split_manifest$id_column)) != 1L ||
    anyNA(split_manifest$id_column) ||
    !nzchar(unique(split_manifest$id_column)[1])) {
  stop_now("The manifest must contain exactly one non-missing ID column.")
}
if (length(unique(split_manifest$split_seed)) != 1L ||
    anyNA(split_manifest$split_seed) ||
    unique(split_manifest$split_seed) != EXPECTED_SPLIT_SEED) {
  stop_now("Split seed does not match the frozen Step_03_1_0 seed.")
}
if (length(unique(split_manifest$training_fraction)) != 1L ||
    anyNA(split_manifest$training_fraction) ||
    abs(unique(split_manifest$training_fraction) - EXPECTED_TRAIN_PROP) > 1e-12) {
  stop_now("Training fraction does not match the frozen 70% split.")
}

ID_COL <- unique(split_manifest$id_column)

training_ids <- sort(
  split_manifest$participant_id[split_manifest$split == "training"],
  method = "radix"
)
test_ids <- sort(
  split_manifest$participant_id[split_manifest$split == "test"],
  method = "radix"
)

if (length(training_ids) != EXPECTED_TRAIN_N ||
    length(test_ids) != EXPECTED_TEST_N ||
    nrow(split_manifest) != EXPECTED_TOTAL_N) {
  stop_now(
    "Frozen split counts differ from the final Step_04_2_1 run: training=",
    length(training_ids), "; test=", length(test_ids),
    "; total=", nrow(split_manifest), "."
  )
}

test_data <- readRDS(TEST_RDS)
if (!is.data.frame(test_data)) {
  stop_now("The frozen test RDS is not a data.frame.")
}
if (!ID_COL %in% names(test_data)) {
  stop_now("Manifest ID column is absent from the test RDS: ", ID_COL)
}

test_data[[ID_COL]] <- trimws(as.character(test_data[[ID_COL]]))
test_ids_observed <- sort(test_data[[ID_COL]], method = "radix")

if (anyNA(test_ids_observed) ||
    any(test_ids_observed == "") ||
    anyDuplicated(test_ids_observed) > 0L) {
  stop_now("Invalid or duplicated participant IDs in the test RDS.")
}
if (!identical(test_ids_observed, test_ids)) {
  stop_now(
    "Test RDS IDs do not exactly match the authoritative Step_03_1_0 ",
    "test IDs in the split manifest."
  )
}
if (length(intersect(test_ids_observed, training_ids)) != 0L) {
  stop_now("Training participants were found in the held-out test RDS.")
}
if (nrow(test_data) != EXPECTED_TEST_N) {
  stop_now("Test RDS row count does not match the frozen manifest.")
}
if ("analysis_split" %in% names(test_data) &&
    any(as.character(test_data$analysis_split) != "test")) {
  stop_now("Non-test rows were found in the frozen test RDS.")
}

required_analysis_cols <- c(
  ID_COL,
  "depression_time",
  "depression_status",
  FOOD_MAP$food,
  COVARS
)
missing_analysis_cols <- setdiff(required_analysis_cols, names(test_data))
if (length(missing_analysis_cols) > 0L) {
  stop_now(
    "Test RDS is missing required variables: ",
    paste(missing_analysis_cols, collapse = ", ")
  )
}

test_data$depression_time <- suppressWarnings(
  as.numeric(as.character(test_data$depression_time))
)
test_data$depression_status <- suppressWarnings(
  as.numeric(as.character(test_data$depression_status))
)

if (anyNA(test_data$depression_time) ||
    any(!is.finite(test_data$depression_time)) ||
    any(test_data$depression_time <= 0)) {
  stop_now("Invalid depression_time values are present in the test RDS.")
}
if (anyNA(test_data$depression_status) ||
    !all(test_data$depression_status %in% c(0, 1))) {
  stop_now("depression_status in the test RDS must be complete and coded 0/1.")
}

# -----------------------------------------------------------------------------
# 3. Import and strictly validate the frozen beta, scaling, and RCS parameters
# -----------------------------------------------------------------------------

message("\n[2/8] Importing the frozen Step_04_2_1 parameters ...")

beta_results <- read.csv(
  BETA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_beta_cols <- c(
  "beta_order", "group", "food", "display_name", "term", "estimate"
)
if (!all(required_beta_cols %in% names(beta_results))) {
  stop_now(
    "Frozen beta file is missing columns: ",
    paste(setdiff(required_beta_cols, names(beta_results)), collapse = ", ")
  )
}

beta_results <- beta_results %>%
  arrange(.data$beta_order)

if (nrow(beta_results) != 15L ||
    anyDuplicated(beta_results$term) > 0L ||
    !identical(as.integer(beta_results$beta_order), seq_len(15L)) ||
    !identical(as.character(beta_results$food), FOOD_MAP$food) ||
    !identical(as.character(beta_results$term), FOOD_MAP$model_var) ||
    any(!is.finite(as.numeric(beta_results$estimate))) ||
    anyNA(beta_results$estimate)) {
  stop_now(
    "Frozen beta file does not exactly match the final 15-component order ",
    "or contains non-finite/duplicated estimates."
  )
}

beta_vector <- as.numeric(beta_results$estimate)
names(beta_vector) <- as.character(beta_results$term)

# The current Step_04_2_1 training run is the source of truth for beta values.
# Do not hard-code coefficients from an older current frozen training archive: the current model
# estimates all 15 coefficients from the frozen training cohort, without sign
# constraints. Their exact values are reconciled with sum_beta below.
beta_sum_from_file <- sum(beta_vector)
beta_abs_sum_from_file <- sum(abs(beta_vector))

if (!is.finite(beta_sum_from_file) ||
    !is.finite(beta_abs_sum_from_file) ||
    beta_abs_sum_from_file <= 0) {
  stop_now("Imported beta coefficients have an invalid or degenerate sum.")
}

beta_sign_pattern <- if (all(beta_vector < 0)) {
  "all_negative"
} else if (all(beta_vector > 0)) {
  "all_positive"
} else if (any(beta_vector < 0) && any(beta_vector > 0)) {
  "mixed_positive_negative"
} else {
  "contains_zero"
}

scaling_parameters <- read.csv(
  SCALING_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(c("parameter", "value", "definition") %in%
         names(scaling_parameters))) {
  stop_now("Frozen scaling parameter file has an unexpected structure.")
}
if (anyDuplicated(as.character(scaling_parameters$parameter)) > 0L) {
  stop_now("Duplicate parameter names found in the frozen scaling file.")
}

frozen_sum_beta <- get_frozen_parameter(
  scaling_parameters, "sum_beta"
)
frozen_sum_abs_beta <- get_frozen_parameter(
  scaling_parameters, "sum_absolute_beta"
)
training_score_min <- get_frozen_parameter(
  scaling_parameters, "diet_score_min_training"
)
training_score_max <- get_frozen_parameter(
  scaling_parameters, "diet_score_max_training"
)
training_score_mean <- get_frozen_parameter(
  scaling_parameters, "diet_score_mean_training"
)
training_score_sd <- get_frozen_parameter(
  scaling_parameters, "diet_score_sd_training"
)
training_t1_upper <- get_frozen_parameter(
  scaling_parameters, "tertile_cutpoint_1_training"
)
training_t2_upper <- get_frozen_parameter(
  scaling_parameters, "tertile_cutpoint_2_training"
)
training_per20_increment <- get_frozen_parameter(
  scaling_parameters, "per20_increment_training"
)
training_quantile_type <- get_frozen_parameter(
  scaling_parameters, "quantile_type"
)
frozen_training_n <- get_frozen_parameter(
  scaling_parameters, "n_training_rds"
)
frozen_training_v5_components_n <- get_frozen_parameter(
  scaling_parameters, "n_V5_score_components"
)
frozen_training_primary_components_n <- get_frozen_parameter(
  scaling_parameters, "n_primary_score_components"
)
frozen_training_complete15_n <- get_frozen_parameter(
  scaling_parameters, "n_training_complete_primary"
)
frozen_training_association_n <- get_frozen_parameter(
  scaling_parameters, "n_training_association_model"
)

parameter_tolerance <- 1e-10
if (!is.finite(frozen_sum_beta) ||
    abs(frozen_sum_beta) <= parameter_tolerance ||
    !is.finite(frozen_sum_abs_beta) ||
    frozen_sum_abs_beta <= 0) {
  stop_now("Frozen beta denominator or absolute-beta sum is invalid.")
}
if (abs(beta_sum_from_file - frozen_sum_beta) >
    parameter_tolerance * max(1, abs(frozen_sum_beta))) {
  stop_now("Imported beta coefficients do not reconcile to the frozen sum_beta.")
}
if (abs(beta_abs_sum_from_file - frozen_sum_abs_beta) >
    parameter_tolerance * max(1, abs(frozen_sum_abs_beta))) {
  stop_now("Imported beta coefficients do not reconcile to sum_absolute_beta.")
}

if (frozen_training_v5_components_n != 15 ||
    frozen_training_primary_components_n != 15) {
  stop_now("Current Step_04_2_1 scaling metadata does not declare all 15 components.")
}

if (frozen_training_n != EXPECTED_TRAIN_N ||
    frozen_training_complete15_n < 1 ||
    frozen_training_complete15_n > frozen_training_n ||
    frozen_training_association_n < 1 ||
    frozen_training_association_n > frozen_training_n) {
  stop_now("Frozen training sample-size parameters are invalid or inconsistent.")
}

if (training_quantile_type != 7) {
  stop_now("The frozen tertile quantile type is not the current frozen training type-7 definition.")
}
if (!is.finite(training_score_min) ||
    !is.finite(training_score_max) ||
    training_score_max <= training_score_min ||
    !is.finite(training_score_mean) ||
    !is.finite(training_score_sd) ||
    training_score_sd <= 0 ||
    training_t1_upper <= training_score_min ||
    training_t2_upper <= training_t1_upper ||
    training_t2_upper >= training_score_max ||
    training_per20_increment <= 0) {
  stop_now("One or more frozen training scaling parameters are invalid.")
}
expected_per20 <- 0.20 * (training_score_max - training_score_min)
if (abs(training_per20_increment - expected_per20) >
    1e-10 * max(1, abs(expected_per20))) {
  stop_now(
    "Frozen per20 increment is inconsistent with 20% of the training score range."
  )
}

rcs_rules <- read.csv(
  RCS_RULE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
required_rule_cols <- c(
  "group", "display_name", "food",
  "value1", "value2", "range_direction", "score_rule"
)
if (!all(required_rule_cols %in% names(rcs_rules))) {
  stop_now(
    "Frozen RCS scoring-rule file is missing columns: ",
    paste(setdiff(required_rule_cols, names(rcs_rules)), collapse = ", ")
  )
}
if (nrow(rcs_rules) != 15L ||
    anyDuplicated(as.character(rcs_rules$food)) > 0L ||
    !setequal(as.character(rcs_rules$food), FOOD_MAP$food)) {
  stop_now("Frozen RCS scoring-rule file does not contain exactly the final 15 foods.")
}

rcs_rules <- rcs_rules[
  match(FOOD_MAP$food, as.character(rcs_rules$food)),
  ,
  drop = FALSE
]
if (!identical(as.character(rcs_rules$food), FOOD_MAP$food) ||
    !identical(as.character(rcs_rules$group), FOOD_MAP$group)) {
  stop_now("Frozen RCS rules do not match the current frozen training component order and labels.")
}
if (anyNA(rcs_rules$range_direction) ||
    any(!rcs_rules$range_direction %in% c("middle", "low", "high", "point")) ||
    anyNA(rcs_rules$value1) ||
    anyNA(rcs_rules$score_rule) ||
    any(trimws(as.character(rcs_rules$score_rule)) == "")) {
  stop_now("At least one frozen RCS scoring rule is missing or unsupported.")
}

# Reconcile the current binary scoring rules with both the Step_04_2_1
# provenance copy and the final Step_04_1_1_4 RCS export. This replaces stale
# hard-coded current frozen training cutpoints and preserves the current V5.3.1 interval rules.
training_range <- read.csv(
  RCS_RANGE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
source_range <- read.csv(
  RCS_SOURCE_RANGE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

required_range_cols <- c(
  "food", "value1", "value2", "range_direction", "score_rule"
)
for (obj_name in c("training_range", "source_range")) {
  obj <- get(obj_name)
  if (!all(required_range_cols %in% names(obj))) {
    stop_now(
      "Range file ", obj_name, " is missing required columns: ",
      paste(setdiff(required_range_cols, names(obj)), collapse = ", ")
    )
  }
  if (nrow(obj) != 15L ||
      anyDuplicated(as.character(obj$food)) > 0L ||
      !setequal(as.character(obj$food), FOOD_MAP$food)) {
    stop_now("Range file ", obj_name, " does not contain exactly the final 15 foods.")
  }
}

training_range <- training_range[
  match(FOOD_MAP$food, as.character(training_range$food)), , drop = FALSE
]
source_range <- source_range[
  match(FOOD_MAP$food, as.character(source_range$food)), , drop = FALSE
]

same_numeric_with_na <- function(x, y, tolerance = 1e-10) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) return(FALSE)
  same_na <- is.na(x) == is.na(y)
  both <- !is.na(x) & !is.na(y)
  all(same_na) && all(abs(x[both] - y[both]) <= tolerance)
}

if (!identical(as.character(rcs_rules$food), FOOD_MAP$food) ||
    !identical(as.character(rcs_rules$group), FOOD_MAP$group)) {
  stop_now("Frozen RCS rules do not match the current 15-component order and labels.")
}

rule_vs_training_range <- identical(
  as.character(rcs_rules$range_direction),
  as.character(training_range$range_direction)
) &&
  same_numeric_with_na(rcs_rules$value1, training_range$value1) &&
  same_numeric_with_na(rcs_rules$value2, training_range$value2) &&
  identical(as.character(rcs_rules$score_rule),
            as.character(training_range$score_rule))

training_vs_source_range <- identical(
  as.character(training_range$range_direction),
  as.character(source_range$range_direction)
) &&
  same_numeric_with_na(training_range$value1, source_range$value1) &&
  same_numeric_with_na(training_range$value2, source_range$value2) &&
  identical(as.character(training_range$score_rule),
            as.character(source_range$score_rule))

if (!rule_vs_training_range || !training_vs_source_range) {
  stop_now(
    "Current RCS scoring rules, Step_04_2_1 range copy, and ",
    "Step_04_1_1_4 source range do not reconcile."
  )
}

frozen_beta_import_audit <- beta_results %>%
  transmute(
    beta_order = as.integer(.data$beta_order),
    group = as.character(.data$group),
    food = as.character(.data$food),
    term = as.character(.data$term),
    beta = as.numeric(.data$estimate),
    hazard_ratio = exp(.data$estimate),
    beta_direction = dplyr::case_when(
      .data$estimate < 0 ~ "negative",
      .data$estimate > 0 ~ "positive",
      TRUE ~ "zero"
    ),
    applied_without_reestimation = TRUE
  )

write_csv_utf8(
  frozen_beta_import_audit,
  file.path(AUDIT_DIR, "frozen_beta_import_audit.csv")
)

frozen_input_files <- c(
  split_manifest = SPLIT_MANIFEST_FILE,
  test_cohort_rds = TEST_RDS,
  beta_coefficients = BETA_FILE,
  scaling_parameters = SCALING_FILE,
  rcs_scoring_rules = RCS_RULE_FILE,
  training_range_copy = RCS_RANGE_FILE,
  final_rcs_source_range = RCS_SOURCE_RANGE_FILE
)
frozen_input_provenance <- tibble(
  input_name = names(frozen_input_files),
  path = unname(frozen_input_files),
  md5 = unname(tools::md5sum(unname(frozen_input_files))),
  size_bytes = as.numeric(file.info(unname(frozen_input_files))$size)
)
write_csv_utf8(
  frozen_input_provenance,
  file.path(AUDIT_DIR, "frozen_input_provenance.csv")
)

# -----------------------------------------------------------------------------
# 4. Apply the frozen RCS binary scoring rules and frozen beta vector to TEST
# -----------------------------------------------------------------------------

message("\n[3/8] Applying frozen RCS rules and beta coefficients to the test cohort ...")

for (i in seq_len(nrow(FOOD_MAP))) {
  food_var <- FOOD_MAP$food[i]
  model_var <- FOOD_MAP$model_var[i]
  rule_row <- rcs_rules[i, , drop = FALSE]

  test_data[[model_var]] <- score_one_component(
    x = test_data[[food_var]],
    direction = as.character(rule_row$range_direction[1]),
    value1 = as.numeric(rule_row$value1[1]),
    value2 = as.numeric(rule_row$value2[1]),
    food_name = food_var
  )
}

score_matrix <- as.matrix(
  test_data[, FOOD_VARS, drop = FALSE]
)
storage.mode(score_matrix) <- "numeric"

test_data$n_components_scored <- as.integer(
  rowSums(!is.na(score_matrix))
)
complete_15 <- test_data$n_components_scored == 15L

test_data$RWDID_total_score_complete15 <- NA_integer_
if (any(complete_15)) {
  test_data$RWDID_total_score_complete15[complete_15] <-
    as.integer(rowSums(
      score_matrix[complete_15, , drop = FALSE],
      na.rm = FALSE
    ))
}

if (sum(complete_15) == 0L) {
  stop_now("No test participant has all 15 RCS-derived component scores observed.")
}
if (any(
  test_data$RWDID_total_score_complete15[complete_15] < 0 |
    test_data$RWDID_total_score_complete15[complete_15] > 15,
  na.rm = TRUE
)) {
  stop_now("The complete 15-component unweighted total is outside 0-15.")
}

weighted_sum <- rep(NA_real_, nrow(test_data))
diet_score <- rep(NA_real_, nrow(test_data))

if (any(complete_15)) {
  weighted_sum[complete_15] <- as.numeric(
    score_matrix[complete_15, , drop = FALSE] %*% beta_vector
  )
  diet_score[complete_15] <- weighted_sum[complete_15] / frozen_sum_beta
}

test_data$weighted_sum <- weighted_sum
test_data$diet_score <- diet_score

if (any(!is.finite(test_data$diet_score[complete_15]))) {
  stop_now("Non-finite test RWDID values were generated.")
}
# With uniformly signed beta coefficients, the normalized score lies in
# [0,1]. If the mutually adjusted model yields mixed signs, the mathematically
# valid attainable bounds are derived from the frozen beta vector itself.
beta_sum_denominator <- sum(beta_vector)
possible_score_endpoints <- c(
  sum(pmin(beta_vector, 0)) / beta_sum_denominator,
  sum(pmax(beta_vector, 0)) / beta_sum_denominator
)
theoretical_score_lower <- min(possible_score_endpoints)
theoretical_score_upper <- max(possible_score_endpoints)

if (any(
  test_data$diet_score[complete_15] <
    theoretical_score_lower - 1e-10 |
    test_data$diet_score[complete_15] >
    theoretical_score_upper + 1e-10
)) {
  stop_now(
    "Test RWDID values fall outside the theoretical bounds implied by ",
    "the frozen beta vector and binary component scores."
  )
}

# Apply the TRAINING cutpoints unchanged. Do not recalculate test quantiles.
test_data$diet_score_tertile <- NA_character_
test_data$diet_score_tertile[complete_15] <- as.character(cut(
  test_data$diet_score[complete_15],
  breaks = c(-Inf, training_t1_upper, training_t2_upper, Inf),
  include.lowest = TRUE,
  right = TRUE,
  labels = c("T1", "T2", "T3")
))
test_data$diet_score_tertile <- factor(
  test_data$diet_score_tertile,
  levels = c("T1", "T2", "T3")
)
test_data$diet_score_tertile_num <- as.numeric(
  test_data$diet_score_tertile
)

# Preserve the current frozen training definition: center by training mean and divide by 20%
# of the observed TRAINING score range. No test-derived rescaling or clipping.
test_data$diet_score_20 <- NA_real_
test_data$diet_score_20[complete_15] <-
  (test_data$diet_score[complete_15] - training_score_mean) /
  training_per20_increment

if (any(!is.finite(test_data$diet_score_20[complete_15]))) {
  stop_now("Non-finite per-20%-training-range test score values were generated.")
}

test_data$analysis_split <- "test"

# Component-level scoring audit.
component_score_audit <- vector("list", nrow(FOOD_MAP))
for (i in seq_len(nrow(FOOD_MAP))) {
  v <- FOOD_MAP$model_var[i]
  z <- test_data[[v]]
  component_score_audit[[i]] <- tibble(
    beta_order = FOOD_MAP$beta_order[i],
    group = FOOD_MAP$group[i],
    food = FOOD_MAP$food[i],
    model_var = v,
    score_rule = as.character(rcs_rules$score_rule[i]),
    n_test = length(z),
    n_score_0 = sum(z == 0L, na.rm = TRUE),
    n_score_1 = sum(z == 1L, na.rm = TRUE),
    n_missing = sum(is.na(z)),
    pct_score_1 = ifelse(
      sum(!is.na(z)) > 0L,
      100 * mean(z == 1L, na.rm = TRUE),
      NA_real_
    )
  )
}
component_score_audit <- bind_rows(component_score_audit)
write_csv_utf8(
  component_score_audit,
  file.path(OUT_DIR, "test_component_score_distribution.csv")
)

# Participant-level scores are retained for reproducibility and downstream QC.
participant_score_output <- test_data %>%
  select(
    all_of(ID_COL),
    all_of(FOOD_MAP$food),
    all_of(FOOD_VARS),
    n_components_scored,
    RWDID_total_score_complete15,
    weighted_sum,
    diet_score,
    diet_score_tertile,
    diet_score_tertile_num,
    diet_score_20,
    analysis_split
  )

write_csv_utf8(
  participant_score_output,
  file.path(OUT_DIR, "RWDID_participant_scores_test.csv")
)

# -----------------------------------------------------------------------------
# 5. Define the common test association cohort and fit the same four Cox models
# -----------------------------------------------------------------------------

message("\n[4/8] Defining the test association cohort and fitting Model 3 Cox models ...")

association_vars <- c(
  "depression_time",
  "depression_status",
  "diet_score",
  "diet_score_tertile",
  "diet_score_tertile_num",
  "diet_score_20",
  COVARS
)

association_complete <- complete.cases(
  test_data[, association_vars, drop = FALSE]
)
association_data <- test_data[association_complete, , drop = FALSE]

if (nrow(association_data) == 0L ||
    sum(association_data$depression_status == 1, na.rm = TRUE) == 0L) {
  stop_now("No eligible complete-case test participants/events for association models.")
}

if (!all(c("T1", "T2", "T3") %in%
         as.character(unique(association_data$diet_score_tertile)))) {
  stop_now(
    "The held-out test association cohort does not contain all three ",
    "training-defined RWDID tertiles. Do not recalculate test cutpoints; ",
    "inspect the test score distribution."
  )
}

tertile_counts <- association_data %>%
  group_by(diet_score_tertile) %>%
  summarise(
    total = n(),
    events = sum(depression_status == 1, na.rm = TRUE),
    mean_diet_score = mean(diet_score, na.rm = TRUE),
    median_diet_score = median(diet_score, na.rm = TRUE),
    min_diet_score = min(diet_score, na.rm = TRUE),
    max_diet_score = max(diet_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_fmt = formatC(total, format = "d", big.mark = ","),
    events_fmt = formatC(events, format = "d", big.mark = ","),
    Events_Total = paste0(events_fmt, "/", total_fmt)
  )

if (nrow(tertile_counts) != 3L ||
    any(tertile_counts$total == 0L) ||
    any(tertile_counts$events == 0L)) {
  print(tertile_counts)
  stop_now(
    "At least one training-defined test tertile is empty or has no events. ",
    "The test-set tertile Cox contrast is not reliably estimable."
  )
}

# Remove unused factor levels after the test complete-case restriction.
# This does not alter any participant's covariate value or the RWDID score.
association_data <- droplevels(association_data)

covariate_variation <- purrr::map_dfr(COVARS, function(v) {
  x <- association_data[[v]]
  n_unique <- dplyr::n_distinct(x[!is.na(x)])
  tibble(
    covariate = v,
    class = paste(class(x), collapse = "/"),
    n_unique_test = n_unique,
    has_variation = n_unique >= 2L
  )
})

write_csv_utf8(
  covariate_variation,
  file.path(AUDIT_DIR, "test_covariate_variation_audit.csv")
)

if (any(!covariate_variation$has_variation)) {
  stop_now(
    "At least one Model 3 covariate has no variation in the held-out test ",
    "association cohort: ",
    paste(
      covariate_variation$covariate[!covariate_variation$has_variation],
      collapse = ", "
    ),
    ". Do not silently remove covariates; inspect the test cohort."
  )
}

fit_continuous <- fit_cox_model(
  make_survival_formula("diet_score", COVARS),
  association_data,
  "test_continuous_RWDID"
)

fit_tertile <- fit_cox_model(
  make_survival_formula("diet_score_tertile", COVARS),
  association_data,
  "test_RWDID_tertiles"
)

fit_trend <- fit_cox_model(
  make_survival_formula("diet_score_tertile_num", COVARS),
  association_data,
  "test_RWDID_tertile_trend"
)

fit_per20 <- fit_cox_model(
  make_survival_formula("diet_score_20", COVARS),
  association_data,
  "test_RWDID_per_20_percent_training_range"
)

continuous_row <- make_result_row(
  variable = "RWDID (per score increase)",
  term = "diet_score",
  fit = fit_continuous,
  count_data = association_data
)

t1_data <- association_data %>%
  filter(diet_score_tertile == "T1")
t2_data <- association_data %>%
  filter(diet_score_tertile == "T2")
t3_data <- association_data %>%
  filter(diet_score_tertile == "T3")

t1_row <- make_result_row(
  variable = "T1",
  term = NA_character_,
  fit = fit_tertile,
  count_data = t1_data,
  reference = TRUE
)
t2_row <- make_result_row(
  variable = "T2",
  term = "diet_score_tertileT2",
  fit = fit_tertile,
  count_data = t2_data
)
t3_row <- make_result_row(
  variable = "T3",
  term = "diet_score_tertileT3",
  fit = fit_tertile,
  count_data = t3_data
)

trend_tidy <- extract_tidy_term(
  fit_trend, "diet_score_tertile_num"
)
trend_row <- tibble(
  variable = "P for trend",
  Events_Total = "",
  hr_ci = "",
  p_value = fmt_p(trend_tidy$p.value[1]),
  hazard_ratio = NA_real_,
  ci_low = NA_real_,
  ci_high = NA_real_,
  p_value_numeric = as.numeric(trend_tidy$p.value[1]),
  term = "diet_score_tertile_num"
)

per20_row <- make_result_row(
  variable = "Per 20% increment",
  term = "diet_score_20",
  fit = fit_per20,
  count_data = association_data
)

test_association_results <- bind_rows(
  continuous_row,
  t1_row,
  t2_row,
  t3_row,
  trend_row,
  per20_row
) %>%
  mutate(
    model_n = nrow(association_data),
    model_events = sum(
      association_data$depression_status == 1,
      na.rm = TRUE
    ),
    analysis_cohort = "Frozen 30% held-out test"
  )

write_csv_utf8(
  test_association_results,
  file.path(OUT_DIR, "test_association_results.csv")
)

# Optional descriptive comparison with apparent training-set estimates.
# This does not change the test analysis and is not a formal interaction test.
training_test_comparison_status <- "NOT_CREATED_OPTIONAL_TRAINING_RESULT_FILE_MISSING"
if (file.exists(TRAIN_ASSOC_FILE)) {
  training_association_results <- read.csv(
    TRAIN_ASSOC_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_training_result_cols <- c(
    "variable", "Events_Total", "hr_ci", "p_value",
    "hazard_ratio", "ci_low", "ci_high", "p_value_numeric",
    "model_n", "model_events"
  )
  if (!all(required_training_result_cols %in%
           names(training_association_results))) {
    stop_now(
      "Optional Step_04_2_1 association-results file has an unexpected structure."
    )
  }

  expected_result_rows <- c(
    "RWDID (per score increase)",
    "T1", "T2", "T3",
    "P for trend",
    "Per 20% increment"
  )
  if (!identical(
        as.character(training_association_results$variable),
        expected_result_rows
      )) {
    stop_now(
      "Step_04_2_1 training association rows do not match the current frozen training order."
    )
  }

  training_comparison <- training_association_results %>%
    transmute(
      variable = as.character(.data$variable),
      training_Events_Total = as.character(.data$Events_Total),
      training_hr_ci = as.character(.data$hr_ci),
      training_p_value = as.character(.data$p_value),
      training_HR = as.numeric(.data$hazard_ratio),
      training_CI_low = as.numeric(.data$ci_low),
      training_CI_high = as.numeric(.data$ci_high),
      training_p_value_numeric = as.numeric(.data$p_value_numeric),
      training_N = as.integer(.data$model_n),
      training_events = as.integer(.data$model_events)
    ) %>%
    inner_join(
      test_association_results %>%
        transmute(
          variable = as.character(.data$variable),
          test_Events_Total = as.character(.data$Events_Total),
          test_hr_ci = as.character(.data$hr_ci),
          test_p_value = as.character(.data$p_value),
          test_HR = as.numeric(.data$hazard_ratio),
          test_CI_low = as.numeric(.data$ci_low),
          test_CI_high = as.numeric(.data$ci_high),
          test_p_value_numeric = as.numeric(.data$p_value_numeric),
          test_N = as.integer(.data$model_n),
          test_events = as.integer(.data$model_events)
        ),
      by = "variable"
    ) %>%
    mutate(
      test_minus_training_logHR = ifelse(
        is.finite(test_HR) & is.finite(training_HR),
        log(test_HR) - log(training_HR),
        NA_real_
      ),
      comparison_type = "Descriptive only; no formal interaction test"
    )

  if (nrow(training_comparison) != length(expected_result_rows) ||
      !identical(training_comparison$variable, expected_result_rows)) {
    stop_now("Could not align all six training and held-out test result rows.")
  }

  write_csv_utf8(
    training_comparison,
    file.path(OUT_DIR, "training_vs_test_association_comparison.csv")
  )
  training_test_comparison_status <- "CREATED"
} else {
  warning(
    "Optional training association results not found; ",
    "test-set analysis will proceed without the descriptive comparison table."
  )
}

write_csv_utf8(
  tertile_counts %>%
    select(
      diet_score_tertile, total, events,
      mean_diet_score, median_diet_score,
      min_diet_score, max_diet_score, Events_Total
    ),
  file.path(OUT_DIR, "test_tertile_distribution.csv")
)

# -----------------------------------------------------------------------------
# 6. Test-set forest plot and RWDID distribution plot
# -----------------------------------------------------------------------------

message("\n[5/8] Creating test-set forest plot and RWDID distribution plot ...")

forest_final_result <- test_association_results
mycolors <- cols4all::c4a("friendly13", 7)

fn <- local({
  i <- 0
  function(..., clr.line, clr.marker) {
    i <<- i + 1
    if (i %% 5 == 1) {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5,
        clr.line = mycolors[1], clr.marker = mycolors[1]
      )
    } else if (i %% 5 == 2) {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5,
        clr.line = mycolors[2], clr.marker = mycolors[2]
      )
    } else if (i %% 5 == 3) {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5,
        clr.line = mycolors[3], clr.marker = mycolors[3]
      )
    } else if (i %% 5 == 4) {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5,
        clr.line = mycolors[4], clr.marker = mycolors[4]
      )
    } else {
      forestplot::fpDrawCircleCI(
        ..., lwd = 1.5,
        clr.line = mycolors[5], clr.marker = mycolors[5]
      )
    }
  }
})

tabletext_HR <- cbind(
  c("Groups", as.character(forest_final_result$variable)),
  c("Events/Total", as.character(forest_final_result$Events_Total)),
  c("HR (95% CI)", as.character(forest_final_result$hr_ci)),
  c("P value", as.character(forest_final_result$p_value))
)

test_forestplot <- forestplot::forestplot(
  labeltext = tabletext_HR,
  graph.pos = 4,
  fn.ci_norm = fn,
  boxsize = 0.2,
  mean = c(NA_real_, forest_final_result$hazard_ratio),
  lower = c(NA_real_, forest_final_result$ci_low),
  upper = c(NA_real_, forest_final_result$ci_high),
  xlog = TRUE,
  col = forestplot::fpColors(zero = "black"),
  zero = 1,
  graphwidth = grid::unit(0.10, "npc"),
  lty.zero = 3,
  lwd.zero = 1.5,
  xlab = "HR (95% CI)",
  lwd.xaxis = 1.5,
  xticks = c(0.12, 0.25, 0.5, 1, 2),
  txt_gp = forestplot::fpTxtGp(
    label = grid::gpar(cex = 1.1),
    ticks = grid::gpar(cex = 0.8),
    xlab = grid::gpar(cex = 0.9)
  )
)

grDevices::pdf(
  file = file.path(OUT_DIR, "RWDID_test_forestplot.pdf"),
  width = 15,
  height = 3
)
print(test_forestplot)
grDevices::dev.off()

test_score_values <- test_data$diet_score[
  is.finite(test_data$diet_score)
]

if (length(test_score_values) != nrow(test_data) ||
    any(!is.finite(test_score_values)) ||
    any(test_score_values < -1e-10 | test_score_values > 1 + 1e-10)) {
  stop_now(
    "The held-out test RWDID distribution contains an unexpected number ",
    "of values or scores outside the theoretical 0-1 range."
  )
}

# Use the exact green-only PiYG palette construction adopted by the
# current Step_04_2_1 training beta-coefficient heatmap.
# A medium green is used for the histogram; the darkest green marks the
# two frozen training-derived tertile cutpoints.
piyg_full <- RColorBrewer::brewer.pal(11, "PiYG")
piyg_green_half <- rev(piyg_full[6:11])
rwdid_histogram_fill <- piyg_green_half[3]
rwdid_cutpoint_color <- piyg_green_half[1]

test_n_for_plot <- length(test_score_values)
cutpoint_1_label <- formatC(
  training_t1_upper, digits = 3, format = "f"
)
cutpoint_2_label <- formatC(
  training_t2_upper, digits = 3, format = "f"
)

p_test_distribution <- ggplot(
  tibble(diet_score = test_score_values),
  aes(x = diet_score)
) +
  geom_histogram(
    bins = 40,
    fill = rwdid_histogram_fill,
    color = "white",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = c(training_t1_upper, training_t2_upper),
    color = rwdid_cutpoint_color,
    linetype = "longdash",
    linewidth = 1.0
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = ggplot2::expansion(mult = c(0, 0.08))
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    x = "RWDID",
    y = "Participants",
    title = "RWDID distribution in the held-out test set",
    subtitle = paste0(
      "Dashed lines indicate frozen training tertile cutpoints: ",
      cutpoint_1_label, " and ", cutpoint_2_label
    ),
    caption = paste0(
      "N = ", formatC(test_n_for_plot, format = "d", big.mark = ","),
      " participants. Cutpoints were estimated in the 70% training set ",
      "and applied unchanged."
    )
  ) +
  theme_classic(
    base_size = 14,
    base_family = "sans"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16.5,
      colour = "#111111",
      hjust = 0,
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = element_text(
      face = "bold",
      size = 11.5,
      colour = "#333333",
      hjust = 0,
      margin = ggplot2::margin(b = 12)
    ),
    plot.caption = element_text(
      size = 9.5,
      colour = "#444444",
      hjust = 0,
      margin = ggplot2::margin(t = 10)
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 13.5,
      colour = "#111111",
      margin = ggplot2::margin(t = 9)
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 13.5,
      colour = "#111111",
      margin = ggplot2::margin(r = 9)
    ),
    axis.text = element_text(
      face = "bold",
      size = 11,
      colour = "#222222"
    ),
    axis.line = element_line(
      colour = "#222222",
      linewidth = 0.75
    ),
    axis.ticks = element_line(
      colour = "#222222",
      linewidth = 0.65
    ),
    axis.ticks.length = grid::unit(2.2, "mm"),
    panel.grid = element_blank(),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.margin = ggplot2::margin(
      t = 12, r = 16, b = 12, l = 12
    )
  )

# Publication-ready exports:
# - PDF and SVG: vector formats for manuscript production/editing.
# - TIFF: 600 dpi with LZW compression for journal submission.
# - PNG: 300 dpi preview for quick inspection.
RWDID_FIGURE_WIDTH_IN <- 7.5
RWDID_FIGURE_HEIGHT_IN <- 5.6

grDevices::pdf(
  file = file.path(OUT_DIR, "RWDID_test_distribution.pdf"),
  width = RWDID_FIGURE_WIDTH_IN,
  height = RWDID_FIGURE_HEIGHT_IN,
  family = "Helvetica",
  pointsize = 12,
  useDingbats = FALSE,
  bg = "white"
)
print(p_test_distribution)
grDevices::dev.off()

grDevices::svg(
  filename = file.path(OUT_DIR, "RWDID_test_distribution.svg"),
  width = RWDID_FIGURE_WIDTH_IN,
  height = RWDID_FIGURE_HEIGHT_IN,
  pointsize = 12,
  bg = "white"
)
print(p_test_distribution)
grDevices::dev.off()

grDevices::tiff(
  filename = file.path(OUT_DIR, "RWDID_test_distribution_600dpi.tiff"),
  width = RWDID_FIGURE_WIDTH_IN,
  height = RWDID_FIGURE_HEIGHT_IN,
  units = "in",
  res = 600,
  compression = "lzw",
  bg = "white"
)
print(p_test_distribution)
grDevices::dev.off()

grDevices::png(
  filename = file.path(OUT_DIR, "RWDID_test_distribution_preview.png"),
  width = RWDID_FIGURE_WIDTH_IN,
  height = RWDID_FIGURE_HEIGHT_IN,
  units = "in",
  res = 300,
  bg = "white"
)
print(p_test_distribution)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# 7. Test-set summaries and model diagnostics
# -----------------------------------------------------------------------------

message("\n[6/8] Exporting test-set summaries and Cox diagnostics ...")

n_test_events_all <- sum(test_data$depression_status == 1, na.rm = TRUE)
n_test_events_association <- sum(
  association_data$depression_status == 1,
  na.rm = TRUE
)

test_score_distribution <- tibble(
  cohort = "Frozen 30% held-out test",
  n_test_rds = nrow(test_data),
  n_test_events_all = n_test_events_all,
  n_complete_15_components = sum(complete_15),
  n_incomplete_15_components = sum(!complete_15),
  n_association_model = nrow(association_data),
  events_association_model = n_test_events_association,
  mean = mean(test_score_values),
  sd = stats::sd(test_score_values),
  median = stats::median(test_score_values),
  p25 = as.numeric(stats::quantile(
    test_score_values, 0.25, names = FALSE, type = 7
  )),
  p75 = as.numeric(stats::quantile(
    test_score_values, 0.75, names = FALSE, type = 7
  )),
  minimum = min(test_score_values),
  maximum = max(test_score_values),
  training_mean_reference = training_score_mean,
  training_sd_reference = training_score_sd,
  training_min_reference = training_score_min,
  training_max_reference = training_score_max,
  training_t1_upper_frozen = training_t1_upper,
  training_t2_upper_frozen = training_t2_upper,
  training_per20_increment_frozen = training_per20_increment,
  n_below_training_min = sum(test_score_values < training_score_min),
  n_above_training_max = sum(test_score_values > training_score_max),
  pct_outside_training_range = 100 * mean(
    test_score_values < training_score_min |
      test_score_values > training_score_max
  ),
  sum_beta_frozen = frozen_sum_beta
)

write_csv_utf8(
  test_score_distribution,
  file.path(OUT_DIR, "test_score_distribution.csv")
)

cox_ph_diagnostics <- bind_rows(
  safe_cox_zph(fit_continuous, "test_continuous_RWDID"),
  safe_cox_zph(fit_tertile, "test_RWDID_tertiles"),
  safe_cox_zph(fit_trend, "test_RWDID_tertile_trend"),
  safe_cox_zph(fit_per20, "test_RWDID_per_20_percent_training_range")
)
write_csv_utf8(
  cox_ph_diagnostics,
  file.path(OUT_DIR, "test_cox_ph_diagnostics.csv")
)

model_warning_summary <- tibble(
  model = names(MODEL_WARNINGS),
  n_warnings = vapply(MODEL_WARNINGS, length, integer(1)),
  warning_messages = vapply(
    MODEL_WARNINGS,
    function(x) paste(unique(x), collapse = " | "),
    character(1)
  )
)
write_csv_utf8(
  model_warning_summary,
  file.path(AUDIT_DIR, "test_cox_model_warnings.csv")
)

# Save all-test score data and the exact common association sample.
test_analysis_data_export <- association_data %>%
  select(
    all_of(ID_COL),
    depression_time,
    depression_status,
    all_of(FOOD_MAP$food),
    all_of(FOOD_VARS),
    all_of(COVARS),
    n_components_scored,
    RWDID_total_score_complete15,
    weighted_sum,
    diet_score,
    diet_score_tertile,
    diet_score_tertile_num,
    diet_score_20
  )

saveRDS(
  test_data,
  file.path(OUT_DIR, "test_score_all_participants.rds")
)
saveRDS(
  test_analysis_data_export,
  file.path(OUT_DIR, "test_analysis_data.rds")
)
saveRDS(
  list(
    continuous_model = fit_continuous,
    tertile_model = fit_tertile,
    trend_model = fit_trend,
    per20_model = fit_per20,
    beta_vector_frozen = beta_vector,
    beta_sum_frozen = frozen_sum_beta,
    beta_sign_pattern_frozen = beta_sign_pattern,
    theoretical_score_lower_frozen = theoretical_score_lower,
    theoretical_score_upper_frozen = theoretical_score_upper,
    training_mean_frozen = training_score_mean,
    training_sd_frozen = training_score_sd,
    training_min_frozen = training_score_min,
    training_max_frozen = training_score_max,
    tertile_cutpoints_frozen = c(
      T1_upper = training_t1_upper,
      T2_upper = training_t2_upper
    ),
    per20_increment_frozen = training_per20_increment,
    rcs_rules_frozen = rcs_rules,
    food_map = FOOD_MAP,
    covariates = COVARS,
    id_column = ID_COL
  ),
  file.path(OUT_DIR, "test_cox_models_and_frozen_parameters.rds")
)


# -----------------------------------------------------------------------------
# ADD-ON A. Held-out RWDID RCS curve using TRAINING-derived knots/reference
#
# The three spline knots and reference value are imported from Step_04_2_1.
# The test cohort estimates only the validation association coefficients;
# neither knots nor the RWDID reference value are re-estimated in the test set.
# -----------------------------------------------------------------------------

message("\n[ADD-ON 1/2] Fitting the held-out test RWDID RCS model ...")

if (!requireNamespace("rms", quietly = TRUE)) {
  stop_now("Package 'rms' is required for the additional test RWDID RCS analysis.")
}

frozen_rcs_knots_table <- read.csv(
  RCS_KNOTS_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(c("knot_order", "knot_value", "reference_value_training") %in%
         names(frozen_rcs_knots_table)) ||
    nrow(frozen_rcs_knots_table) != 3L) {
  stop_now("The frozen training RCS knot file has an unexpected structure.")
}
frozen_rcs_knots_table <- frozen_rcs_knots_table %>%
  dplyr::arrange(knot_order)
frozen_rcs_knots <- as.numeric(frozen_rcs_knots_table$knot_value)
frozen_rcs_reference <- unique(
  as.numeric(frozen_rcs_knots_table$reference_value_training)
)
if (length(frozen_rcs_reference) != 1L ||
    !is.finite(frozen_rcs_reference) ||
    length(frozen_rcs_knots) != 3L ||
    any(!is.finite(frozen_rcs_knots)) ||
    is.unsorted(frozen_rcs_knots, strictly = TRUE)) {
  stop_now("The frozen training RCS knots/reference are invalid.")
}

fit_rwdid_rcs_addon <- function(data, knots, reference_value,
                                cohort_label, file_prefix) {
  if (!is.data.frame(data) ||
      !all(c("depression_time", "depression_status", "diet_score", COVARS) %in%
           names(data))) {
    stop_now("Test RCS input data are missing the outcome, RWDID, or Model 3 covariates.")
  }
  if (!is.finite(reference_value) ||
      length(knots) != 3L ||
      any(!is.finite(knots)) ||
      is.unsorted(knots, strictly = TRUE)) {
    stop_now("Invalid frozen training RCS knots/reference supplied to the test model.")
  }

  # Build datadist only from variables actually used in this RCS model.
  # This avoids warnings from unrelated constant audit/metadata columns.
  rcs_datadist_vars <- unique(c("diet_score", COVARS))
  ddist_data <- data[, rcs_datadist_vars, drop = FALSE]
  ddist <- rms::datadist(ddist_data)

  if (!"diet_score" %in% colnames(ddist$limits)) {
    stop_now("rms::datadist did not retain the test RWDID variable.")
  }
  # Use the training median as the common HR=1 reference in both cohorts.
  ddist$limits["Adjust to", "diet_score"] <- reference_value

  # rms::Design() resolves options(datadist = "<object name>") from an
  # environment on the search path. A function-local object is not sufficient.
  # Register a temporary object in .GlobalEnv and restore prior state on exit.
  datadist_object_name <- ".rwdid_rcs_datadist_addon"
  had_previous_datadist_object <- exists(
    datadist_object_name, envir = .GlobalEnv, inherits = FALSE
  )
  if (had_previous_datadist_object) {
    previous_datadist_object <- get(
      datadist_object_name, envir = .GlobalEnv, inherits = FALSE
    )
  }

  assign(datadist_object_name, ddist, envir = .GlobalEnv)
  previous_datadist_option <- options(datadist = datadist_object_name)

  on.exit({
    options(previous_datadist_option)
    if (had_previous_datadist_object) {
      assign(
        datadist_object_name, previous_datadist_object,
        envir = .GlobalEnv
      )
    } else if (exists(datadist_object_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = datadist_object_name, envir = .GlobalEnv)
    }
  }, add = TRUE)

  knots_text <- paste(
    format(knots, digits = 17, scientific = FALSE, trim = TRUE),
    collapse = ", "
  )
  rcs_formula <- stats::as.formula(
    paste0(
      "survival::Surv(depression_time, depression_status == 1) ~ ",
      "rms::rcs(diet_score, c(", knots_text, ")) + ",
      paste(COVARS, collapse = " + ")
    )
  )

  fit <- rms::cph(
    rcs_formula,
    data = data,
    x = TRUE,
    y = TRUE
  )

  anova_tab <- as.data.frame(stats::anova(fit))
  anova_rows <- trimws(rownames(anova_tab))
  p_col <- names(anova_tab)[tolower(names(anova_tab)) == "p"]
  if (length(p_col) != 1L) {
    p_col <- grep("p", names(anova_tab), ignore.case = TRUE, value = TRUE)[1]
  }
  if (is.na(p_col) || !nzchar(p_col)) {
    stop_now("Could not identify the P-value column in the test rms ANOVA table.")
  }

  overall_idx <- which(anova_rows == "diet_score")
  if (!length(overall_idx)) {
    overall_idx <- which(
      grepl("diet_score", anova_rows, fixed = TRUE) &
        !grepl("Nonlinear", anova_rows, ignore.case = TRUE)
    )
  }
  nonlinear_idx <- which(grepl("Nonlinear", anova_rows, ignore.case = TRUE))
  if (!length(overall_idx)) {
    stop_now("Could not identify the overall RWDID term in the test rms ANOVA table.")
  }

  p_overall <- as.numeric(anova_tab[overall_idx[1], p_col])
  p_nonlinear <- if (length(nonlinear_idx)) {
    as.numeric(anova_tab[nonlinear_idx[1], p_col])
  } else {
    NA_real_
  }

  pred <- as.data.frame(
    rms::Predict(
      fit,
      diet_score,
      fun = exp,
      type = "predictions",
      ref.zero = TRUE
    )
  )
  required_prediction_cols <- c("diet_score", "yhat", "lower", "upper")
  if (!all(required_prediction_cols %in% names(pred)) ||
      nrow(pred) < 10L ||
      any(!is.finite(pred$yhat)) ||
      any(!is.finite(pred$lower)) ||
      any(!is.finite(pred$upper))) {
    stop_now("The test RWDID RCS prediction table is incomplete or non-finite.")
  }

  pred$cohort <- cohort_label
  write_csv_utf8(
    pred,
    paste0(file_prefix, "_curve_data.csv")
  )

  x_range <- range(pred$diet_score, na.rm = TRUE)
  y_upper <- max(pred$upper, na.rm = TRUE)
  y_lower <- min(pred$lower, na.rm = TRUE)
  x_annotate <- x_range[1] + 0.025 * diff(x_range)
  y_annotate <- y_upper * 1.13
  p_label <- paste0(
    "P for overall ", fmt_p(p_overall),
    "\nP for non-linear ", fmt_p(p_nonlinear)
  )

  friendly_colours <- cols4all::c4a("friendly13")
  rcs_plot <- ggplot2::ggplot(
    pred,
    ggplot2::aes(x = diet_score, y = yhat)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.10,
      fill = friendly_colours[1]
    ) +
    ggplot2::geom_line(
      linewidth = 1,
      alpha = 0.70,
      colour = friendly_colours[1]
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = 2,
      linewidth = 0.5,
      colour = friendly_colours[2],
      alpha = 0.50
    ) +
    ggplot2::annotate(
      "text",
      x = x_annotate,
      y = y_annotate,
      label = p_label,
      hjust = 0,
      vjust = 1,
      size = 4.5
    ) +
    ggplot2::scale_x_continuous(name = "RWDID") +
    ggplot2::scale_y_continuous(
      name = "Hazard ratio (95% CI)",
      expand = ggplot2::expansion(mult = c(0.03, 0.18))
    ) +
    ggplot2::coord_cartesian(
      ylim = c(max(0, y_lower * 0.90), y_annotate * 1.03)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(
        colour = "black", size = 18, face = "bold"
      ),
      axis.title = ggplot2::element_text(
        colour = "black", size = 18, face = "bold"
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 18
      ),
      axis.title.x = ggplot2::element_text(
        colour = "black", size = 18, face = "bold"
      ),
      axis.title.y = ggplot2::element_text(
        colour = "black", size = 18, face = "bold"
      )
    ) +
    ggplot2::ggtitle(cohort_label)

  ggplot2::ggsave(
    filename = paste0(file_prefix, ".pdf"),
    plot = rcs_plot,
    width = 10,
    height = 6,
    device = "pdf"
  )
  ggplot2::ggsave(
    filename = paste0(file_prefix, "_600dpi.png"),
    plot = rcs_plot,
    width = 10,
    height = 6,
    dpi = 600,
    bg = "white"
  )

  model_summary <- tibble(
    cohort = cohort_label,
    n = nrow(data),
    events = sum(data$depression_status == 1, na.rm = TRUE),
    reference_value_training_median = reference_value,
    knot_1_training = knots[1],
    knot_2_training = knots[2],
    knot_3_training = knots[3],
    p_overall = p_overall,
    p_nonlinear = p_nonlinear
  )
  write_csv_utf8(
    model_summary,
    paste0(file_prefix, "_model_summary.csv")
  )

  list(
    fit = fit,
    prediction = pred,
    plot = rcs_plot,
    knots = knots,
    reference_value = reference_value,
    summary = model_summary
  )
}

rcs_test_addon <- fit_rwdid_rcs_addon(
  data = association_data,
  knots = frozen_rcs_knots,
  reference_value = frozen_rcs_reference,
  cohort_label = "UK Biobank - Held-out test cohort",
  file_prefix = file.path(PLOT_DIR, "RWDID_RCS_test")
)
saveRDS(
  rcs_test_addon$fit,
  file.path(OUT_DIR, "RWDID_RCS_test_model.rds")
)

# -----------------------------------------------------------------------------
# ADD-ON B. Publication-ready Table 1 by RWDID tertiles
#
# Two cohort-specific tables are produced (test only and pooled total cohort).
# The pooled table combines the saved training association sample with this
# held-out test association sample; it retains the training-defined T1/T2/T3
# assignments and does not re-estimate tertile cutpoints.
# -----------------------------------------------------------------------------

message("\n[ADD-ON 2/2] Creating test and full-cohort publication Table 1 ...")

if (!requireNamespace("gtsummary", quietly = TRUE) ||
    !requireNamespace("flextable", quietly = TRUE)) {
  stop_now(
    "Packages 'gtsummary' and 'flextable' are required for the additional ",
    "publication-ready Table 1."
  )
}

make_rwdid_publication_table1 <- function(data, output_prefix, cohort_note) {
  table1_vars <- c(
    "age_food", "follow_up_years", "sex", "ethnicity",
    "education_attainment", "career", "twnsend", "smoking_status",
    "drink_status", "bmi_group", "household_income",
    "physical_activity_proup", "sedentary_group", "energy_intake",
    "central_obesity", "triglycerides", "hdl_cholesterol",
    "hypertension_status_final", "t2diabetes_status_final",
    "dyslipidemia_status_final", "cvd_status_final",
    "all_cancer_status_final"
  )
  continuous_vars <- c(
    "age_food", "follow_up_years", "education_attainment",
    "twnsend", "energy_intake"
  )
  binary_vars <- c(
    "central_obesity", "triglycerides", "hdl_cholesterol",
    "hypertension_status_final", "t2diabetes_status_final",
    "dyslipidemia_status_final", "cvd_status_final",
    "all_cancer_status_final"
  )

  required <- c(
    "depression_time", "diet_score", "diet_score_tertile",
    table1_vars[table1_vars != "follow_up_years"]
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop_now(
      "Table 1 source data are missing variables: ",
      paste(missing, collapse = ", ")
    )
  }

  table_data <- as.data.frame(data)
  table_data$follow_up_years <- as.numeric(table_data$depression_time) / 365.25
  if (any(!is.finite(table_data$follow_up_years)) ||
      any(table_data$follow_up_years <= 0)) {
    stop_now("Derived follow-up years contain invalid values.")
  }

  table_data$diet_score_tertile <- factor(
    as.character(table_data$diet_score_tertile),
    levels = c("T1", "T2", "T3")
  )
  if (anyNA(table_data$diet_score_tertile)) {
    stop_now("Table 1 contains missing or invalid training-defined RWDID tertiles.")
  }

  order_factor <- function(x, preferred_levels) {
    z <- trimws(as.character(x))
    observed <- unique(z[!is.na(z) & nzchar(z)])
    factor(
      z,
      levels = c(
        preferred_levels[preferred_levels %in% observed],
        setdiff(observed, preferred_levels)
      )
    )
  }

  table_data$sex <- order_factor(
    table_data$sex, c("Male", "Female")
  )
  table_data$ethnicity <- order_factor(
    table_data$ethnicity, c("Non-white", "White")
  )
  table_data$career <- order_factor(
    table_data$career, c("Employed", "Retired/other")
  )
  table_data$smoking_status <- order_factor(
    table_data$smoking_status, c("Current", "Never", "Previous")
  )
  table_data$drink_status <- order_factor(
    table_data$drink_status, c("Current", "Never", "Previous")
  )

  bmi_values <- trimws(as.character(table_data$bmi_group))
  bmi_values[bmi_values == "Obese (>=30)"] <- "Obese (≥30)"
  table_data$bmi_group <- order_factor(
    bmi_values,
    c(
      "Underweight (<18.5)",
      "Normal weight (18.5 to <25)",
      "Overweight (25 to <30)",
      "Obese (≥30)"
    )
  )
  table_data$household_income <- order_factor(
    table_data$household_income, c("Low", "Middle", "High")
  )
  table_data$physical_activity_proup <- order_factor(
    table_data$physical_activity_proup,
    c("High physical activity", "Low physical activity")
  )
  table_data$sedentary_group <- order_factor(
    table_data$sedentary_group, c("Not sedentary", "Sedentary")
  )

  normalize_binary <- function(x) {
    z <- trimws(as.character(x))
    z_lower <- tolower(z)
    map <- c(
      "0" = "No", "1" = "Yes",
      "no" = "No", "yes" = "Yes",
      "false" = "No", "true" = "Yes"
    )
    observed <- unique(z_lower[!is.na(z_lower) & nzchar(z_lower)])
    if (all(observed %in% names(map))) {
      factor(unname(map[z_lower]), levels = c("No", "Yes"))
    } else {
      factor(z, levels = unique(z[!is.na(z)]))
    }
  }
  for (v in binary_vars) {
    table_data[[v]] <- normalize_binary(table_data[[v]])
  }

  for (v in continuous_vars) {
    original <- table_data[[v]]
    converted <- suppressWarnings(as.numeric(as.character(original)))
    if (any(!is.na(original) & is.na(converted))) {
      stop_now("A Table 1 continuous variable is not numeric: ", v)
    }
    table_data[[v]] <- converted
  }

  dichotomous_vars <- binary_vars[vapply(
    binary_vars,
    function(v) {
      observed <- unique(as.character(table_data[[v]][!is.na(table_data[[v]])]))
      length(observed) > 0L && all(observed %in% c("No", "Yes"))
    },
    logical(1)
  )]

  type_spec <- c(
    list(
      c(
        age_food, follow_up_years, education_attainment,
        twnsend, energy_intake
      ) ~ "continuous"
    ),
    lapply(
      dichotomous_vars,
      function(v) stats::as.formula(paste0(v, " ~ 'dichotomous'"))
    )
  )
  value_spec <- lapply(
    dichotomous_vars,
    function(v) stats::as.formula(paste0(v, " ~ 'Yes'"))
  )

  table1_labels <- list(
    age_food ~ "Age, year, median (IQR)",
    follow_up_years ~ "Follow-up, years, median (IQR)",
    sex ~ "Sex, no. (%)",
    ethnicity ~ "Ethnicity, no. (%)",
    education_attainment ~ "Education attainment, mean (SD)",
    career ~ "Employment status, no. (%)",
    twnsend ~ "Townsend deprivation index, mean (SD)",
    smoking_status ~ "Smoking status, no. (%)",
    drink_status ~ "Drinking status, no. (%)",
    bmi_group ~ "BMI category (kg/m2), no. (%)",
    household_income ~ "Family income, no. (%)",
    physical_activity_proup ~ "Physical activity, no. (%)",
    sedentary_group ~ "Sedentary behavior, no. (%)",
    energy_intake ~ "Energy intake, KJ/day, mean (SD)",
    central_obesity ~ "Central obesity, no. (%)",
    triglycerides ~ "High triglycerides, no. (%)",
    hdl_cholesterol ~ "Low HDL, no. (%)",
    hypertension_status_final ~ "Hypertension, no. (%)",
    t2diabetes_status_final ~ "Diabetes, no. (%)",
    dyslipidemia_status_final ~ "Dyslipidemia, no. (%)",
    cvd_status_final ~ "CVD, no. (%)",
    all_cancer_status_final ~ "Cancer, no. (%)"
  )

  tbl <- gtsummary::tbl_summary(
    data = table_data,
    by = all_of("diet_score_tertile"),
    include = all_of(table1_vars),
    label = table1_labels,
    type = type_spec,
    value = value_spec,
    statistic = list(
      age_food ~ "{median} ({p25}, {p75})",
      follow_up_years ~ "{median} ({p25}, {p75})",
      education_attainment ~ "{mean} ({sd})",
      twnsend ~ "{mean} ({sd})",
      energy_intake ~ "{mean} ({sd})",
      gtsummary::all_categorical() ~ "{n} ({p}%)",
      gtsummary::all_dichotomous() ~ "{n} ({p}%)"
    ),
    digits = list(
      age_food ~ 0,
      follow_up_years ~ 2,
      education_attainment ~ 2,
      twnsend ~ 2,
      energy_intake ~ 0,
      gtsummary::all_categorical() ~ c(0, 1),
      gtsummary::all_dichotomous() ~ c(0, 1)
    ),
    missing = "no"
  ) %>%
    gtsummary::add_overall(last = FALSE)

  group_n <- table_data %>%
    dplyr::count(diet_score_tertile, .drop = FALSE)
  if (nrow(group_n) != 3L || any(group_n$n <= 0L)) {
    stop_now("The publication Table 1 requires non-empty T1, T2, and T3 groups.")
  }

  tb <- as.data.frame(tbl$table_body)
  if (!all(c("variable", "row_type", "label") %in% names(tb))) {
    stop_now("Unexpected gtsummary table_body structure.")
  }
  stat_cols <- grep("^stat_[0-9]+$", names(tb), value = TRUE)
  if (!identical(stat_cols, c("stat_0", "stat_1", "stat_2", "stat_3"))) {
    stop_now(
      "Expected Total plus three RWDID tertile columns in Table 1; found: ",
      paste(stat_cols, collapse = ", ")
    )
  }

  tb$variable <- as.character(tb$variable)
  tb$row_type <- as.character(tb$row_type)
  tb$label <- as.character(tb$label)
  for (v in stat_cols) tb[[v]] <- as.character(tb[[v]])

  insert_section_row <- function(x, before_variable, section_label) {
    idx <- which(x$variable == before_variable)[1]
    if (is.na(idx)) {
      stop_now("Could not insert Table 1 section header before ", before_variable)
    }
    new_row <- x[NA_integer_, , drop = FALSE]
    new_row$variable <- paste0(".section_", make.names(section_label))
    new_row$row_type <- "label"
    new_row$label <- section_label
    for (v in stat_cols) new_row[[v]] <- ""
    dplyr::bind_rows(
      x[seq_len(idx - 1L), , drop = FALSE],
      new_row,
      x[seq.int(idx, nrow(x)), , drop = FALSE]
    )
  }

  tb <- insert_section_row(
    tb, "central_obesity", "Metabolic syndrome"
  )
  tb <- insert_section_row(
    tb, "dyslipidemia_status_final", "Other diseases"
  )

  table_export <- tb %>%
    dplyr::select(
      Characteristic = label,
      dplyr::all_of(stat_cols)
    )

  fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
  total_label <- paste0("Total (n = ", fmt_n(nrow(table_data)), ")")
  t1_label <- paste0("T1 (n = ", fmt_n(group_n$n[group_n$diet_score_tertile == "T1"]), ")")
  t2_label <- paste0("T2 (n = ", fmt_n(group_n$n[group_n$diet_score_tertile == "T2"]), ")")
  t3_label <- paste0("T3 (n = ", fmt_n(group_n$n[group_n$diet_score_tertile == "T3"]), ")")

  ft <- flextable::flextable(table_export)
  ft <- flextable::set_header_labels(
    ft,
    Characteristic = "Characteristic",
    stat_0 = total_label,
    stat_1 = t1_label,
    stat_2 = t2_label,
    stat_3 = t3_label
  )
  ft <- flextable::add_header_row(
    ft,
    values = c("Characteristic", total_label, "RWDID"),
    colwidths = c(1, 1, 3),
    top = TRUE
  )
  ft <- flextable::merge_v(
    ft,
    j = c("Characteristic", "stat_0"),
    part = "header"
  )
  section_rows <- which(
    table_export$Characteristic %in% c("Metabolic syndrome", "Other diseases")
  )
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
  ft <- flextable::fontsize(ft, size = 8.5, part = "all")
  ft <- flextable::bold(ft, part = "header", bold = TRUE)
  if (length(section_rows)) {
    ft <- flextable::bold(
      ft, i = section_rows, j = "Characteristic",
      bold = TRUE, part = "body"
    )
    ft <- flextable::bg(
      ft, i = section_rows, bg = "#F2F2F2", part = "body"
    )
  }
  ft <- flextable::align(ft, j = "Characteristic", align = "left", part = "all")
  ft <- flextable::align(ft, j = stat_cols, align = "center", part = "all")
  ft <- flextable::valign(ft, valign = "center", part = "all")
  ft <- flextable::padding(ft, padding = 2, part = "all")
  ft <- flextable::width(ft, j = "Characteristic", width = 2.25)
  ft <- flextable::width(ft, j = stat_cols, width = 1.05)
  ft <- flextable::set_table_properties(ft, layout = "fixed", width = 1)
  ft <- flextable::add_footer_lines(
    ft,
    values = paste0(
      "Continuous variables are presented as median (IQR) or mean (SD), ",
      "as labelled; categorical variables are presented as n (%). ",
      "Percentages are column percentages. ", cohort_note
    )
  )

  csv_out <- table_export
  names(csv_out) <- c(
    "Characteristic", total_label, t1_label, t2_label, t3_label
  )
  utils::write.csv(
    csv_out,
    file = paste0(output_prefix, ".csv"),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
  flextable::save_as_docx(
    "Table 1" = ft,
    path = paste0(output_prefix, ".docx")
  )
  saveRDS(
    tbl,
    file = paste0(output_prefix, "_gtsummary.rds")
  )

  list(
    table = tbl,
    flextable = ft,
    n = nrow(table_data),
    tertile_n = group_n,
    output_docx = paste0(output_prefix, ".docx"),
    output_csv = paste0(output_prefix, ".csv")
  )
}

table1_test_addon <- make_rwdid_publication_table1(
  data = association_data,
  output_prefix = file.path(
    OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_test"
  ),
  cohort_note = "Held-out test cohort; T1/T2/T3 cutpoints were imported unchanged from the training cohort."
)

# Create the requested full-cohort table (Total n = training + test) by
# combining the already saved training association cohort and current test
# association cohort. Both cohorts already carry training-defined tertiles.
training_table_source <- readRDS(TRAIN_ANALYSIS_RDS)
if (!is.data.frame(training_table_source)) {
  stop_now("The saved Step_04_2_1 training analysis RDS is not a data.frame.")
}
pooled_required <- c(
  ID_COL, "depression_time", "depression_status",
  "diet_score", "diet_score_tertile", COVARS
)
missing_pooled <- setdiff(pooled_required, names(training_table_source))
if (length(missing_pooled)) {
  stop_now(
    "The saved training analysis RDS is missing fields for the pooled Table 1: ",
    paste(missing_pooled, collapse = ", ")
  )
}
pooled_table_data <- dplyr::bind_rows(
  training_table_source,
  association_data
)
if (anyDuplicated(as.character(pooled_table_data[[ID_COL]])) > 0L) {
  stop_now("Duplicate participant IDs were found in the pooled Table 1 source.")
}
if (nrow(pooled_table_data) !=
    nrow(training_table_source) + nrow(association_data)) {
  stop_now("Pooled Table 1 row count does not equal training plus test rows.")
}
if (nrow(pooled_table_data) > EXPECTED_TOTAL_N) {
  stop_now("Pooled Table 1 contains more rows than the frozen total cohort.")
}

table1_pooled_addon <- make_rwdid_publication_table1(
  data = pooled_table_data,
  output_prefix = file.path(
    OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_full_cohort"
  ),
  cohort_note = "Full 70/30 cohort; tertile cutpoints were derived in the training cohort and applied unchanged to the held-out test cohort."
)

addon_qc_test <- tibble(
  check = c(
    "test_RWDID_RCS_model_fitted",
    "test_RWDID_RCS_curve_pdf_exported",
    "test_RWDID_RCS_curve_600dpi_png_exported",
    "test_RWDID_RCS_prediction_data_exported",
    "test_RWDID_RCS_knots_equal_training_knots",
    "test_Table1_word_exported",
    "test_Table1_csv_exported",
    "pooled_full_cohort_Table1_word_exported",
    "pooled_full_cohort_Table1_csv_exported",
    "pooled_Table1_uses_training_defined_tertiles"
  ),
  status = c(
    inherits(rcs_test_addon$fit, "cph"),
    file.exists(file.path(PLOT_DIR, "RWDID_RCS_test.pdf")),
    file.exists(file.path(PLOT_DIR, "RWDID_RCS_test_600dpi.png")),
    file.exists(file.path(PLOT_DIR, "RWDID_RCS_test_curve_data.csv")),
    isTRUE(all.equal(rcs_test_addon$knots, frozen_rcs_knots)) &&
      isTRUE(all.equal(rcs_test_addon$reference_value, frozen_rcs_reference)),
    file.exists(file.path(OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_test.docx")),
    file.exists(file.path(OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_test.csv")),
    file.exists(file.path(OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_full_cohort.docx")),
    file.exists(file.path(OUT_DIR, "Table_1_baseline_by_RWDID_tertiles_full_cohort.csv")),
    nrow(table1_pooled_addon$tertile_n) == 3L &&
      nrow(pooled_table_data) ==
        nrow(training_table_source) + nrow(association_data)
  ),
  detail = c(
    paste0("N=", nrow(association_data),
           "; events=", sum(association_data$depression_status == 1)),
    file.path(PLOT_DIR, "RWDID_RCS_test.pdf"),
    file.path(PLOT_DIR, "RWDID_RCS_test_600dpi.png"),
    file.path(PLOT_DIR, "RWDID_RCS_test_curve_data.csv"),
    paste0("knots=", paste(frozen_rcs_knots, collapse = " / "),
           "; training reference=", frozen_rcs_reference),
    table1_test_addon$output_docx,
    table1_test_addon$output_csv,
    table1_pooled_addon$output_docx,
    table1_pooled_addon$output_csv,
    paste0(
      "pooled N=", nrow(pooled_table_data),
      "; expected frozen total N=", EXPECTED_TOTAL_N,
      "; T1/T2/T3 N=",
      paste(table1_pooled_addon$tertile_n$n, collapse = "/")
    )
  )
)
write_csv_utf8(
  addon_qc_test,
  file.path(OUT_DIR, "RCS_Table1_addon_QC.csv")
)
if (any(!addon_qc_test$status)) {
  print(addon_qc_test %>% dplyr::filter(!status))
  stop_now("The additional test RCS/Table 1 QC did not pass.")
}


# -----------------------------------------------------------------------------
# 8. Final QC gate, README, and run summary
# -----------------------------------------------------------------------------

message("\n[7/8] Running final test-set QC gate ...")

expected_output_files <- c(
  "audit/frozen_beta_import_audit.csv",
  "audit/frozen_input_provenance.csv",
  "audit/test_covariate_variation_audit.csv",
  "audit/test_cox_model_warnings.csv",
  "test_component_score_distribution.csv",
  "RWDID_participant_scores_test.csv",
  "test_association_results.csv",
  "test_tertile_distribution.csv",
  "RWDID_test_forestplot.pdf",
  "RWDID_test_distribution.pdf",
  "RWDID_test_distribution.svg",
  "RWDID_test_distribution_600dpi.tiff",
  "RWDID_test_distribution_preview.png",
  "test_score_distribution.csv",
  "test_cox_ph_diagnostics.csv",
  "test_score_all_participants.rds",
  "test_analysis_data.rds",
  "test_cox_models_and_frozen_parameters.rds"
)

qc_checks <- c(
  "manifest_seed_matches_20261103",
  "manifest_training_fraction_matches_70_percent",
  "manifest_training_n_matches_76578",
  "manifest_test_n_matches_32819",
  "test_RDS_exactly_matches_manifest_test_IDs",
  "test_participants_excluded_from_training_IDs",
  "test_RDS_contains_only_test_rows",
  "test_outcome_time_and_status_are_valid",
  "frozen_beta_file_contains_exactly_15_ordered_components",
  "all_15_frozen_beta_values_are_finite_and_sum_is_stable",
  "sum_of_frozen_beta_matches_frozen_scaling_parameter",
  "frozen_training_scaling_parameters_are_valid",
  "frozen_per20_increment_matches_20_percent_training_range",
  "frozen_RCS_rules_contain_exactly_15_components",
  "RCS_rules_match_Step_04_2_1_range_copy",
  "Step_04_2_1_range_copy_matches_final_Step_04_1_1_4_export",
  "all_15_test_component_scores_are_0_1_or_NA",
  "test_RWDID_only_calculated_for_complete_15_components",
  "test_RWDID_is_within_beta_implied_theoretical_bounds",
  "test_tertiles_use_frozen_training_cutpoints",
  "all_three_training_defined_test_tertiles_are_present",
  "all_three_test_tertiles_have_events",
  "all_21_Model3_covariates_have_test_sample_variation",
  "continuous_tertile_trend_and_per20_models_fitted",
  "all_four_models_use_the_same_test_association_cohort",
  "test_forestplot_exists",
  "test_RWDID_distribution_plot_exists",
  "test_RWDID_distribution_svg_exists",
  "test_RWDID_distribution_600dpi_tiff_exists",
  "test_RWDID_distribution_preview_png_exists",
  "test_association_results_file_exists",
  "test_PH_diagnostics_file_exists",
  "test_analysis_RDS_exists",
  "test_model_RDS_exists",
  "all_core_validation_outputs_exist"
)

qc_status <- c(
  unique(split_manifest$split_seed) == EXPECTED_SPLIT_SEED,
  abs(unique(split_manifest$training_fraction) - EXPECTED_TRAIN_PROP) <= 1e-12,
  length(training_ids) == EXPECTED_TRAIN_N,
  length(test_ids) == EXPECTED_TEST_N,
  identical(test_ids_observed, test_ids),
  length(intersect(test_ids_observed, training_ids)) == 0L,
  !("analysis_split" %in% names(test_data)) ||
    all(as.character(test_data$analysis_split) == "test"),
  all(is.finite(test_data$depression_time)) &&
    all(test_data$depression_time > 0) &&
    all(test_data$depression_status %in% c(0, 1)),
  nrow(beta_results) == 15L &&
    identical(as.character(beta_results$term), FOOD_VARS),
  all(is.finite(beta_vector)) &&
    abs(sum(beta_vector)) > parameter_tolerance,
  abs(beta_sum_from_file - frozen_sum_beta) <=
    parameter_tolerance * max(1, abs(frozen_sum_beta)),
  all(is.finite(c(
    training_score_min, training_score_max, training_score_mean,
    training_score_sd, training_t1_upper, training_t2_upper,
    training_per20_increment
  ))) &&
    training_score_max > training_score_min &&
    training_score_sd > 0 &&
    training_t1_upper < training_t2_upper,
  abs(training_per20_increment - expected_per20) <=
    1e-10 * max(1, abs(expected_per20)),
  nrow(rcs_rules) == 15L &&
    identical(as.character(rcs_rules$food), FOOD_MAP$food),
  isTRUE(rule_vs_training_range),
  isTRUE(training_vs_source_range),
  all(vapply(
    FOOD_VARS,
    function(v) all(is.na(test_data[[v]]) |
                      test_data[[v]] %in% c(0L, 1L)),
    logical(1)
  )),
  all(is.na(test_data$diet_score[!complete_15])) &&
    all(is.finite(test_data$diet_score[complete_15])),
  all(
    test_data$diet_score[complete_15] >= theoretical_score_lower - 1e-10 &
      test_data$diet_score[complete_15] <= theoretical_score_upper + 1e-10
  ),
  all(is.na(test_data$diet_score_tertile[!complete_15])) &&
    all(is.na(test_data$diet_score_tertile_num[!complete_15])) &&
    all(is.na(test_data$diet_score_20[!complete_15])) &&
    identical(
      levels(test_data$diet_score_tertile),
      c("T1", "T2", "T3")
    ),
  all(c("T1", "T2", "T3") %in%
        as.character(unique(association_data$diet_score_tertile))),
  all(tertile_counts$events > 0L),
  all(covariate_variation$has_variation),
  length(list(
    fit_continuous, fit_tertile, fit_trend, fit_per20
  )) == 4L &&
    nrow(test_association_results) == 6L,
  all(vapply(
    list(fit_continuous, fit_tertile, fit_trend, fit_per20),
    function(x) {
      identical(as.integer(x$n), as.integer(nrow(association_data))) &&
        identical(as.integer(x$nevent), as.integer(n_test_events_association))
    },
    logical(1)
  )),
  file.exists(file.path(OUT_DIR, "RWDID_test_forestplot.pdf")),
  file.exists(file.path(OUT_DIR, "RWDID_test_distribution.pdf")),
  file.exists(file.path(OUT_DIR, "RWDID_test_distribution.svg")),
  file.exists(file.path(OUT_DIR, "RWDID_test_distribution_600dpi.tiff")),
  file.exists(file.path(OUT_DIR, "RWDID_test_distribution_preview.png")),
  file.exists(file.path(OUT_DIR, "test_association_results.csv")),
  file.exists(file.path(OUT_DIR, "test_cox_ph_diagnostics.csv")),
  file.exists(file.path(OUT_DIR, "test_analysis_data.rds")),
  file.exists(file.path(OUT_DIR, "test_cox_models_and_frozen_parameters.rds")),
  all(file.exists(file.path(OUT_DIR, expected_output_files)))
)

qc_detail <- c(
  paste0("seed=", unique(split_manifest$split_seed)),
  paste0("training_fraction=", unique(split_manifest$training_fraction)),
  paste0("training_n=", length(training_ids)),
  paste0("test_n=", length(test_ids)),
  paste0(
    "test_RDS_n=", nrow(test_data),
    "; exact manifest ID match=", identical(test_ids_observed, test_ids)
  ),
  paste0("training/test ID overlap=", length(intersect(test_ids_observed, training_ids))),
  "Test RDS analysis_split is test when the column is present",
  paste0(
    "test_events_all=", n_test_events_all,
    "; valid time/status rows=", nrow(test_data)
  ),
  paste0("frozen beta coefficients=", nrow(beta_results)),
  paste0(
    "beta range=", min(beta_vector), "/", max(beta_vector),
    "; sign pattern=", beta_sign_pattern
  ),
  paste0("sum_beta=", format(frozen_sum_beta, digits = 12)),
  paste0(
    "frozen training mean=", format(training_score_mean, digits = 10),
    "; SD=", format(training_score_sd, digits = 10),
    "; min/max=", format(training_score_min, digits = 8), "/",
    format(training_score_max, digits = 8)
  ),
  paste0("frozen per20 increment=", format(training_per20_increment, digits = 12)),
  paste0("frozen RCS rules=", nrow(rcs_rules)),
  "RCS binary rules exactly match Step_04_2_1 range export",
  "Step_04_2_1 range export exactly matches Step_04_1_1_4 final RCS export",
  paste0(
    "test component scores: complete15=", sum(complete_15),
    "; incomplete15=", sum(!complete_15)
  ),
  "Only complete 15-component test scores receive RWDID",
  paste0(
    "test RWDID range=",
    format(min(test_score_values), digits = 8), " to ",
    format(max(test_score_values), digits = 8),
    "; theoretical bounds=",
    format(theoretical_score_lower, digits = 8), " to ",
    format(theoretical_score_upper, digits = 8)
  ),
  "Training tertile cutpoints applied unchanged; test tertiles were not re-estimated",
  paste0(
    "test tertile n=",
    paste(tertile_counts$total, collapse = "/")
  ),
  paste0(
    "test tertile events=",
    paste(tertile_counts$events, collapse = "/")
  ),
  paste0(
    "non-constant Model 3 covariates=",
    sum(covariate_variation$has_variation), "/",
    length(COVARS)
  ),
  "Four prespecified current frozen training Cox association specifications were fitted",
  paste0(
    "common test association N=", nrow(association_data),
    "; all four models have the same N and event count"
  ),
  "RWDID_test_forestplot.pdf",
  "RWDID_test_distribution.pdf",
  "RWDID_test_distribution.svg",
  "RWDID_test_distribution_600dpi.tiff",
  "RWDID_test_distribution_preview.png",
  "test_association_results.csv",
  "test_cox_ph_diagnostics.csv",
  "test_analysis_data.rds",
  "test_cox_models_and_frozen_parameters.rds",
  paste0(
    "existing core outputs=",
    sum(file.exists(file.path(OUT_DIR, expected_output_files))),
    "/", length(expected_output_files)
  )
)

final_qc <- tibble(
  check = qc_checks,
  status = as.logical(qc_status),
  detail = qc_detail
)

write_csv_utf8(
  final_qc,
  file.path(OUT_DIR, "test_analysis_QC.csv")
)

if (any(!final_qc$status)) {
  print(final_qc %>% filter(!status))
  stop_now(
    "STOP: One or more final Step_04_2_2 validation QC checks failed. ",
    "Inspect test_analysis_QC.csv; do not interpret the run as finalized."
  )
}

readme_lines <- c(
  "STEP_04_2_2 HELD-OUT TEST-SET VALIDATION OF FROZEN RWDID",
  "",
  "Design:",
  "- Uses the exact frozen 30% Step_03_1_0 test cohort.",
  "- Split seed = 20261103; training fraction = 0.70.",
  "- No resampling, cohort rebuilding, or participant reassignment.",
  "- No beta coefficients are re-estimated in the test cohort.",
  "- No RCS cutpoints or binary scoring rules are re-estimated.",
  "- No tertile thresholds or per-20% scaling parameters are recalculated.",
  "",
  "Frozen current Step_04_2_1 training inputs:",
  "- The imported 15 beta estimates are used unchanged and reconciled to the frozen sum_beta and sum_absolute_beta parameters.",
  "- Training-derived scaling parameters are imported unchanged and checked for internal consistency.",
  "- RCS rules are reconciled against the Step_04_2_1 range copy and the Step_04_1_1_4 final RCS range export.",
  paste0("- Beta coefficients: ", BETA_FILE),
  paste0("- Training scaling parameters: ", SCALING_FILE),
  paste0("- RCS binary scoring rules: ", RCS_RULE_FILE),
  paste0("- Step_04_2_1 range copy: ", RCS_RANGE_FILE),
  paste0("- Step_04_1_1_4 final RCS range: ", RCS_SOURCE_RANGE_FILE),
  paste0("- Split manifest: ", SPLIT_MANIFEST_FILE),
  paste0("- Held-out test RDS: ", TEST_RDS),
  "",
  "Scoring:",
  "- The exact final 15 components and order from the current Step_04_1_1_4 / Step_04_2_1 workflow are used.",
  "- The original formula is retained:",
  "  weighted_sum = sum(binary_component_score * frozen_beta_i)",
  "  RWDID = weighted_sum / frozen_sum_beta.",
  "- Missing component scores remain NA; no component score is imputed.",
  "- RWDID is calculated only when all 15 component scores are observed.",
  "- Test scores are not clipped to the training score range.",
  "- T1/T2/T3 use the two frozen training cutpoints without test-set recalculation.",
  "- Per-20% scaling uses the frozen training mean and training range increment.",
  "",
  "Association analyses:",
  "- Continuous RWDID (per one score-unit increase).",
  "- T2 and T3 versus T1 using training-defined tertiles.",
  "- P for trend using training-defined tertile codes 1/2/3.",
  "- Per 20% increment of the training score range.",
  "- All four models use the same 21 Model 3 covariates and common test sample.",
  "- Cox ties = Efron; 95% confidence intervals.",
  "- Schoenfeld residual proportional-hazards diagnostics are exported.",
  "",
  "Figure styling (V1_1):",
  "- The RWDID distribution histogram uses the same green-only PiYG palette construction as the Step_04_2_1 beta heatmap.",
  "- Main title, subtitle, axis titles, and axis tick labels use enlarged/bold typography on a clean white background.",
  "- Frozen training tertile cutpoints remain unchanged and are shown as dark-green dashed lines.",
  "- Figure exports: vector PDF, vector SVG, 600-dpi LZW-compressed TIFF, and 300-dpi PNG preview.",
  "- This update changes figure presentation only; no cohort, score, beta, cutpoint, scaling, or Cox model calculation is changed.",
  "",
  "Important interpretation:",
  "- These are held-out test-set association estimates, not beta-development estimates.",
  "- Test tertiles may have unequal sample sizes because training cutpoints are fixed.",
  "- The test set must not be used to revise the score or its frozen parameters.",
  "",
  "Core outputs:",
  "- test_association_results.csv",
  "- training_vs_test_association_comparison.csv (if the Step_04_2_1 association CSV is present)",
  "- RWDID_test_forestplot.pdf",
  "- RWDID_test_distribution.pdf (vector)",
  "- RWDID_test_distribution.svg (vector)",
  "- RWDID_test_distribution_600dpi.tiff (600 dpi, LZW)",
  "- RWDID_test_distribution_preview.png (300 dpi)",
  "- test_score_distribution.csv",
  "- test_tertile_distribution.csv",
  "- test_component_score_distribution.csv",
  "- test_cox_ph_diagnostics.csv",
  "- test_analysis_QC.csv",
  "- RWDID_participant_scores_test.csv",
  "- test_analysis_data.rds",
  "- test_cox_models_and_frozen_parameters.rds",
  "- audit/frozen_beta_import_audit.csv",
  "- audit/frozen_input_provenance.csv",
  "- audit/test_covariate_variation_audit.csv",
  "- audit/test_cox_model_warnings.csv"
)

writeLines(
  readme_lines,
  con = file.path(OUT_DIR, "README.txt"),
  useBytes = TRUE
)

cat("\n============================================================\n")
cat("STEP_04_2_2 HELD-OUT TEST-SET VALIDATION COMPLETED\n")
cat("Frozen test cohort N: ", nrow(test_data), "\n", sep = "")
cat("Test events (all): ", n_test_events_all, "\n", sep = "")
cat("Complete 15-component test scores: ", sum(complete_15), "\n", sep = "")
cat("Test association N: ", nrow(association_data), "\n", sep = "")
cat("Test association events: ", n_test_events_association, "\n", sep = "")
cat("Training-vs-test comparison: ", training_test_comparison_status, "\n", sep = "")
cat("Frozen sum(beta): ", format(frozen_sum_beta, digits = 12), "\n", sep = "")
cat(
  "Frozen training tertile cutpoints: ",
  format(training_t1_upper, digits = 8), " / ",
  format(training_t2_upper, digits = 8), "\n",
  sep = ""
)
cat("Test tertile N: ", paste(tertile_counts$total, collapse = " / "), "\n", sep = "")
cat("Test tertile events: ", paste(tertile_counts$events, collapse = " / "), "\n", sep = "")
cat("RWDID distribution figure: PDF / SVG / 600-dpi TIFF / 300-dpi PNG\n")
cat("QC: ", file.path(OUT_DIR, "test_analysis_QC.csv"), "\n", sep = "")
cat("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")

########################################################################################
# END ANALYSIS BLOCK 11/11: Step_04_2_2_HELDOUT_VALIDATION
########################################################################################
