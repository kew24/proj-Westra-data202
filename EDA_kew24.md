Exploratory Data Analysis
================
Kaitlyn Westra, partner: Danish Ansari
November 2020

Background Information
----------------------

-   *Topic or Overall Question:* Genetic Mutations arising from Cancer Therapies

-   *Data sources I'm considering:* <https://github.com/papaemmelab/bolton_NG_CH>

-   *Why I'm interested:*
    -   *Nature Genetics* is one of the most well-known and well-respected journals in the field of genetics. The work that these articles are based on is of high quality, and presents some of the latest and greatest ideas in the field. Due to this, many experts frequently read these recently published articles, and regard getting published in the journal as a great achievement.
    -   Over the past 6 months, I've worked as a research intern at Grand Rapids' Van Andel Institute, in a [lab](https://braslab.vai.org) that studies the genetics of neurodegenerative diseases, like Parkinson's Disease, Alzhiemer's Disease, and dementia with Lewy bodies. Throughout this internship, I have been able to learn about the intersection of genetics, disease, and data science, using R and other genomics tools to analyze genetic mutation data and make discoveries that have an impact on people's lives.
    -   Because of this introduction I've had to the field of genetics annd bioinformatics, I have been curious what *other* labs and groups do with similar data. This lead me to browse through the [cBioPortal for Cancer Genomics](https://www.cbioportal.org), which I had heard about from a series of online [Dataviz + Cancer Microlabs](https://apply.hub.ki/cancerplusviz/). The dataset I had found was used in a recent Nature Genetics article, "[Cancer therapy shapes the fitness landscape of clonal hematopoiesis](https://www.nature.com/articles/s41588-020-00710-0)" so I decided this would bee perfect for my final project.

Setup
-----

Read in Data from repo (code copied from author's .Rmd):

``` r
#Input data
M_long = suppressWarnings(data.table::fread('M_long.txt', sep = '\t', header = T)) %>%
  as.data.frame()

unzip(zipfile = 'M_wide_all.txt.zip')
M_wide_all = suppressWarnings(data.table::fread('M_wide_all.txt', sep = '\t', header = T)) %>%
  as.data.frame()

M2_all = suppressWarnings(data.table::fread('M2_all.txt', sep = '\t', header = T))  %>%
  as.data.frame()

M_tmn_wide_st = suppressWarnings(data.table::fread('M_tmn_wide_st.txt', sep = '\t', header = T)) %>% 
  as.data.frame()

M_tmn_h = suppressWarnings(data.table::fread('M_tmn_h.txt', sep = '\t', header = T))  %>%
  as.data.frame()

M_tmn_long = suppressWarnings(data.table::fread('M_tmn_long.txt', sep = '\t', header = T))  %>%
  as.data.frame()

P_serial = suppressWarnings(data.table::fread('P_serial.txt', sep = '\t', header = T))  %>%
  as.data.frame()

T_seer = read_tsv('breast_prop_chemo.txt', col_types = cols()) 

tmn_inc = read_tsv('CumIncByYear-20190716_KB_Wolff.txt', col_types = cols()) %>%
    dplyr::select(Year, yearly_inc) %>% as.data.frame %>% unname %>% as.matrix

mort_inc = read_tsv('overall_mort_50_75.txt', col_types = cols()) %>%
    dplyr::select(time, year_inc) %>% as.data.frame %>% unname %>% as.matrix


#Dictionaries:
class_dict <- read.csv("chemoclass_jan2019_revised.csv", stringsAsFactors = F)

drugsets_dict <- read.csv("top_sets_dec2018.csv", header = T, stringsAsFactors = F)

class_dict <- class_dict %>% filter(drugclass_c=="cytotoxic_therapy")
```

Exploratory Data Analysis
-------------------------

-   Specific questions I'm thinking of asking:
-   which genes (Gene) predict specific types of tumors (generaltumortype)?
-   which genes (Gene) have more mutations in smokers (smoke)?
-   which genes (Gene) have more mutations in patients that have undergone cancer therapy (ind\_target\_therapy)(ind\_immune\_therapy)(XRT)(ind\_cytotoxic\_therapy)?
-   which variant class (silent, missense, splice, frameshift, etc) or variant type (deletion, insertion, SNV) leads
-   Characteristics of the data that will allow me to answer those questions: \_\_\_\_\_ (e.g., names of specific columns in specific datasets that have that data) -columns include: "generaltumortype", "Gene", "VariantClass", "variant\_type", "myeloid\_gene", "smoke"/"smoke\_bin", "therapy\_binary", "ind\_anychemo", "ind\_radiotherapy", "ind\_cytotoxic\_therapy".
