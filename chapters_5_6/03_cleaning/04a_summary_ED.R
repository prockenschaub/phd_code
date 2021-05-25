###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     04a_summary_ED.R
# Date:     25/02/2020
# Task:     Create summaries for the ED dataset
#
###########################################################################

source("00_init.R")

library(tableone)
library(kableExtra)
library(ggpubr)
library(ggplot2)
library(corrplot)
library(mice)

.dir_imp <- file.path("02_imported", vers)
.dir_root <- "03_cleaning"
.dir_der <- file.path(.dir_root, "01_derived")
.dir_data <- file.path(.dir_root, "02_analysis")


# Load the data -----------------------------------------------------------

# Datasets
full_data <- readRDS(file.path(.dir_data, "all_pat_ed.rds"))
urine <- readRDS(file.path(.dir_der, "urine.rds"))
sens <- readRDS(file.path(.dir_imp, "sens.rds"))


# Combine variables for easier display
full_data[, ethnicity := fct_collapse(ethnicity, other = c("mixed", "other"))]

age_map <- c(rep("18-64", 5), rep("65+", 4))
names(age_map) <- levels(full_data$age)
full_data[, age_bin := fct_relabel(age, function(x) age_map[as.character(x)])]

full_data[, cci_tri := cut(cci, c(-Inf, 0, 2, Inf), c("0", "1-2", "$geq$3"))]

full_data[, hosp_12m := fct_yesno(hosp_n_12m > 0)]
#full_data[, ddd_yes := fct_yesno(ddd_12m > 0)]

# Limit only to patients whose samples were cultured
full_data %<>% .[results == "yes"]

# Limit to main analysis observations
data <- full_data[arrival_date >= ymd("2011-11-01") & hosp_uti_30d == "no"]



# Table 1a ----------------------------------------------------------------

create_table1 <- function(total, strat, shell, nonnormal = NULL, 
                          catDigits = 1, contDigits = 2, pDigits = 3){
  # Make a latex version of tables 1a and 1b
  
  shell %<>% .[!is.na(label)]
  
  # Extract the table from the `tableone` object
  tbl_1 <- 
    list(total, strat) %>% 
    map(print, dropEqual = TRUE, printToggle = FALSE, nonnormal = nonnormal,
        catDigits = catDigits, contDigits = contDigits, pDigits = pDigits) %>% 
    reduce(cbind) %>%
    as.data.table(keep.rownames = "var") %>% 
    .[, !("test")]
  
  tbl_1 %<>% .[, map(.SD, str_trim)]
  
  # Remove the headers
  tbl_1[ , p := if_else(p == "" & lag(p) != "", lag(p), p, p)] 
  tbl_1 %<>% .[Overall != ""]
  
  # Merge to the table shell
  tbl_1 %<>% .[shell, on = "var", nomatch = 0]
  tbl_1[, label := str_replace(label, "\\%", "\\\\%")]
  tbl_1[, label := str_replace(label, "\\&", "\\\\&")]
  
  # Bold those columns without headers
  tbl_1[is.na(header), label := cell_spec(label, "latex", bold = TRUE, escape = FALSE)]
  
  # Render the table in latex
  tbl_1_rend <- 
    tbl_1[, c("label", "Overall", "yes", "no", "p")] %>% 
    kable(
      "latex", 
      booktabs = TRUE,
      linesep = "",
      col.names = c("", "Overall", "Yes", "No", "p-value"),
      escape = FALSE,
      caption = "Table 1" 
    )
  
  # Add row headers
  rhead_1 <- tbl_1 %>% 
    .[, .(header, row = 1:.N)] %>% 
    .[!is.na(header), .(start = min(row), end = max(row)), by = header]
  
  for(i in 1:nrow(rhead_1)){
    tbl_1_rend %<>% 
      pack_rows(rhead_1[i]$header, rhead_1[i]$start, rhead_1[i]$end)
  }
  
  # Add column headers
  chead_1 <- c(2, 2, 1)
  names(chead_1) <- c(" ", linebreak("Bacterial growth", "c", T), " ")
  
  tbl_1_rend %<>% add_header_above(chead_1, escape = FALSE)
  
  # Print the table
  tbl_1_rend
}

to_row_percent <- function(table){
  # Convert percentages in tableone (which only provides
  # column percentages) to row percentages
  #
  # Parameters
  # ----------
  # table : TableOne object
  #
  # Returns
  # -------
  # TableOne with row percent
  
  for(v in names(table$CatTable$no)){
    no <- table$CatTable$no[[v]]$freq
    yes <- table$CatTable$yes[[v]]$freq
    
    tot <- yes + no
    
    table$CatTable$no[[v]]$percent <- no / tot * 100
    table$CatTable$yes[[v]]$percent <- yes / tot * 100
  }
  
  table
}



# Table 1a: patient characteristics and medical history
covar <- c("age_bin", "sex", "ethnicity", 
           "cci_tri", "cancer", "immuno", "renal", "uro", "renal_surg", 
           "hosp_12m", "hosp_uti_12m", "urine_12m", "pos_12m", "abx_12m")


tbl_1a_total <- CreateTableOne(covar, data = data)
tbl_1a_strat <- CreateTableOne(covar, strata = "growth", data = data, )

shell_1a <- read_excel(file.path(.dir_root, "04a_table_shells.xlsx"), 
                       sheet = "table_1a")
setDT(shell_1a)

create_table1(tbl_1a_total, to_row_percent(tbl_1a_strat), shell_1a)


# Table 1b: ED diagnosis
covar <- c("susp")

tbl_1b_total <- CreateTableOne(covar, data = data)
tbl_1b_strat <- CreateTableOne(covar, strata = "growth", data = data)

shell_1b <- read_excel(file.path(.dir_root, "04a_table_shells.xlsx"), 
                       sheet = "table_1b")
setDT(shell_1b)

create_table1(tbl_1b_total, to_row_percent(tbl_1b_strat), shell_1b)

for(l in shell_1b$var[-1]){
  cat(l, ":", round(chisq.test(data$susp == l, data$growth)$p.value, 3), "\n")
}



# Table 1c: clinical presentation

covar <- c("ua_bacteria", "ua_epithelial", "ua_rbc_total", "ua_wbc",
           "hr", "rr", "temp", "bp", "o2",
           "crp", "wcc", "plats", "creat", "bili", "alp")

# Scale the variables for easier display
scl_data <- copy(data)
scl_data[, ua_bacteria := ua_bacteria / 1000]


tbl_1c_total <- CreateTableOne(covar, data = scl_data)
tbl_1c_strat <- CreateTableOne(covar, strata = "growth", data = scl_data)

           
shell_1c <- read_excel(file.path(.dir_root, "04a_table_shells.xlsx"), 
                       sheet = "table_1c")
setDT(shell_1c)        
           
create_table1(tbl_1c_total, tbl_1c_strat, shell_1c, covar, contDigits = 1)



# Plot distribution across years ------------------------------------------

prop_by_quarter <- full_data %>% 
  .[!year %in% c("2019")] %>% 
  .[, quarter := year(arrival_date) + 0.125 + 0.25 * (quarter(arrival_date) - 1) ] %>%
  .[, .(.N, prop = mean(growth == "yes")), by = quarter] %>% 
  .[, sd := sqrt(prop * (1-prop) / N)] %>% 
  .[]

total_cultured <- full_data %>% 
  .[, .(year = year(arrival_date))] %>% 
  .[, .N, by = year]

total_urines <- urine %>% 
  .[, .(year = year(collect_date))] %>% 
  .[, .N, by = year]

g1 <- ggplot(full_data, aes(as.numeric(as.character(year)) + 0.5)) + 
  geom_bar(width = 0.8, alpha = 0.8) + 
  geom_segment(aes(x = year + 0.1, xend = year + 0.9, y = N, yend = N),
               data = total_urines, 
               colour = "#414487FF", size = 1) + 
  geom_vline(xintercept = 2011.83) +
  geom_vline(xintercept = 2015.75) +
  scale_x_continuous(breaks = 2010 + 0:9, labels = c(2010 + 0:8, ""), expand = expansion(0, 0)) + 
  scale_y_continuous(limits = c(0, 4000), expand = expansion(0, 0)) + 
  annotate("segment", x = 2011.83, xend = 2012.03, y = 3850, yend = 3850) + 
  annotate("text", x = 2012.1, y = 3850, label = "ED diagnoses recorded", hjust = 0) + 
  annotate("text", x = 2012.1, y = 3650, label = "(November 2011)", size = 2, hjust = 0) + 
  annotate("segment", x = 2015.75, xend = 2015.95, y = 3850, yend = 3850) + 
  annotate("text", x = 2016, y = 3850, label = "Culture thresholds change", hjust = 0) + 
  annotate("text", x = 2016, y = 3650, label = "(October 2015)", size = 2, hjust = 0) + 
  labs(x = "",
       y = "Number of urine samples\n") +
  coord_cartesian(xlim = c(2010, 2019)) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(hjust = -1), 
    panel.grid.minor = element_blank()
  )


g2 <- ggplot(prop_by_quarter, aes(x = quarter, y = prop, 
                         ymin = prop - 1.96 * sd,
                         ymax = prop + 1.96 * sd,
                         group = quarter > 2015.75)) + 
  geom_vline(xintercept = 2011.83) +
  geom_vline(xintercept = 2015.75) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, colour = "#414487FF") +
  scale_x_continuous(breaks = 2010 + 0:9, 
                     labels = c(2010 + 0:8, ""), 
                     expand = expansion(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.5),
                     expand = expansion(0, 0), 
                     labels = scales::percent) + 
  annotate("segment", x = 2011.83, xend = 2012.03, y = 0.085, yend = 0.085) + 
  annotate("text", x = 2012.1, y = 0.085, label = "ED diagnoses recorded", hjust = 0) + 
  annotate("text", x = 2012.1, y = 0.06, label = "(November 2011)", size = 2, hjust = 0) + 
  annotate("segment", x = 2015.75, xend = 2015.95, y = 0.085, yend = 0.085) + 
  annotate("text", x = 2016, y = 0.085, label = "Culture thresholds change", hjust = 0) + 
  annotate("text", x = 2016, y = 0.06, label = "(October 2015)", size = 2, hjust = 0) + 
  coord_cartesian(xlim = c(2010, 2019)) + 
  labs(x = "\n Year",
       y = "Proportion of culture-positive samples\n") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = -1), 
        panel.grid.minor = element_blank())



ggarrange(g1, g2, labels = c("A", "B"), ncol = 1, nrow = 2, align = "v")
ggsave(file.path(.dir_root, "04_summary", "samples_by_year.png"), 
       width = 7, height = 9, dpi = 600)


# Plot missing pattern ----------------------------------------------------

# Get missingness pattern
mispat <- data[, .SD, .SDcols = covar] %>% 
  md.pattern(plot = FALSE) %>% 
  as.data.table(keep.rownames = TRUE)

mispat %<>% 
  .[order(-as.numeric(rn))] %>% 
  .[, .SD, .SDcols = c("rn", covar)]

# Label and bring counts in right order
mispat[, order_y := 1:.N]
mispat[, cumprop := .(cumsum(as.integer(rn)) / sum(as.integer(rn), 
                                                   na.rm = TRUE))]
mispat[, label := str_c(scales::number(as.integer(rn), big.mark = ","), 
                        " (", round(cumprop * 100), "%)")]

# Bring into long format for plotting
var_order <- seq_along(covar)
names(var_order) <- covar

mispat %<>% 
  select(-rn, -cumprop) %>% 
  gather("var", "obs", -order_y, -label)
  
mispat %<>% mutate(order_x = var_order[var])
  
# Distinguish variable types  
mispat %<>% mutate(
  type = case_when(
    str_starts(var, "ua_") ~ "Flow cytometry", 
    var %in% c("hr", "rr", "temp", "bp", "o2", "avpu") ~ "Vital signs", 
    TRUE ~ "Blood biomarkers"
  )
)

# Also get the missing percent per var across all patients
misvar <- data[, map(.SD, ~ mean(!is.na(.))), .SDcols = covar]
misvar %<>%
  unlist() %>% 
  (function(x) round(x * 100)) %>% 
  str_c("%")
  
# Define plot  
g <- 
  ggplot(NULL, aes(
    x = order_x, 
    y = as.character(order_y), 
    fill = type,
    alpha = obs * 0.3 + 0.7)
  ) + 
  geom_tile(color = "black") +
  scale_x_continuous(
    position = "top",
    breaks = var_order,
    labels = c("BAC", "EPC", "RBC", "WBC", 
              "HR", "RR", "TMP", "SBP", "O2",
               "CRP", "WBC", "PLT", "CRE", "BIL", "ALP"),
    sec.axis = sec_axis(
      ~ .,
      breaks = var_order, 
      labels = misvar
    )
  ) + 
  scale_y_discrete(
    limits = as.character(10:1), 
    labels = mispat$label[10:1]
  ) + 
  scale_fill_viridis_d(
    "Predictor type:",
    begin = 0.2,
    limits = c("Flow cytometry", "Vital signs", "Blood biomarkers")
  ) + 
  scale_alpha(range = c(0, 1)) + 
  guides(
    fill = guide_legend(title = ""),
    alpha = FALSE
  ) + 
  labs(
    x = "Predictors\n",
    y = "Number (cum-%) of patients with a measurement\n"
  ) + 
  coord_cartesian(expand = FALSE) + 
  theme_bw() + 
  theme(
    axis.title.x.top = element_text(),
    axis.text.x.bottom = element_text(margin = margin(t = 10)),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

g %+% (mispat %>% filter(order_y <= 10))

ggsave(file.path(.dir_root, "04_summary", "missing_pattern.png"), 
       width = 7, height = 4, dpi = 600)



# Correlations ------------------------------------------------------------

.variables <- c(
  AGE = "age", GND = "sex",
  TMP = "temp", RR = "rr", HR = "hr", O2 = "o2", BP = "bp", SEWS = "sews",
  UAB = "ua_bacteria", UACT = "ua_casts_total", UACP = "ua_casts_prcnt",
  UACO = "ua_conductivity", UAE = "ua_epithelial", UAO = "ua_other", 
  UART = "ua_rbc_total", UARP = "ua_rbc_prcnt", UAS = "ua_sml_rnd_cells",
  UAW = "ua_wbc", UACY = "ua_crystals",
  WBC = "wcc", PLT = "plats", CRP = "crp", CREA = "creat", ALP = "alp", 
  BILI = "bili",
  CCI = "cci", CANC = "cancer", REN = "renal", URO = "uro", 
  SURG = "renal_surg",
  UTI = "hosp_uti_12m", UTIN = "hosp_uti_n_24m", 
  EUTI = "ed_uti_12m", EUTIN = "ed_uti_n_24m", 
  URIN = "urine_12m", POS = "pos_12m", 
  HOSP = "hosp_n_12m", HD12 = "hosp_days_12m", H7 = "hosp_7d", 
  ED = "ed_n_12m", ABX = "abx" 
)

corr_data <- data[, map(.SD, as.numeric), .SDcols = unlist(.variables)]
setnames(corr_data, .variables, names(.variables))
corr <- cor(corr_data, 
            method = "spearman",
            use = "pairwise.complete.obs")


for(i in 2:ncol(corr)){
  for(j in 1:(i-1)){
    if(abs(corr[i, j]) > 0.8){
      cat(rownames(corr)[i], "-", colnames(corr)[j], ":", corr[i, j], "\n")
    }
  }
}


corrplot(corr, 
         method = "color",
         type = "full",
         order = "hclust",
         hclust.method = "complete",
         tl.col = "black",
         tl.cex = 0.7)


# Additional calculations -------------------------------------------------

# Isolated organisms
orgs <- sens %>% 
  .[, .(pat_id, analysis_id, org)] %>% 
  unique() %>% 
  .[(data[, .(pat_id, sex, older = fct_yesno(as.numeric(age) >= 6), idx_urine)]), on = .(pat_id, analysis_id = idx_urine)] %>% 
  .[!is.na(org) & org != ""]

tail(orgs$org %>% table() %>% sort())
tail(orgs$org %>% table() %>% prop.table() %>%  sort())

tail(orgs[sex == "female"]$org %>% table() %>% prop.table() %>%  sort())
tail(orgs[sex == "male"]$org %>% table() %>% prop.table() %>%  sort())

tail(orgs[older == "no"]$org %>% table() %>% prop.table() %>%  sort())
tail(orgs[older == "yes"]$org %>% table() %>% prop.table() %>%  sort())


