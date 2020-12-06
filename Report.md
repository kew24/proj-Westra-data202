Genetic Mutations from Cancer Therapies
================
Kaitlyn Westra
Fall 2020

-   [About](#about)
-   [Introduction](#introduction)
-   [Question](#question)
    -   [Background](#background)
    -   [Project](#project)
-   [Dataset](#dataset)
-   [(?) Data Wrangling](#data-wrangling)
-   [Exploratory Data Analysis](#exploratory-data-analysis)
-   [Modeling](#modeling)
-   [Findings:](#findings)
-   [Limitations:](#limitations)
-   [Future Directions:](#future-directions)
-   [Conclusions](#conclusions)

### About

This is an outline of the things that need to be included for this final project. Below, I've copied & pasted the project outline from Professor Arnold's project description. Interspersed is bullet points of specific rubric objectives (copied & pasted from the "3 points" column in Moodle). This hopefully makes it easier to see what progress has been made & what we have left to do. At the very end, when these sections are filled in, the filler words can be deleted.

Introduction
------------

*Nature Genetics* is one of the most well-known and well-respected journals in the field of genetics. The work that these articles are based on reflect years of rigorous research, and presents some of the latest and greatest ideas in the field. Due to this, many experts regard getting published in the journal as a great achievement and frequently read these published articles.

Over the past 6 months, I've worked as a research intern at Grand Rapids' Van Andel Institute, in a [lab](https://braslab.vai.org) that studies the genetics of neurodegenerative diseases, like Parkinson's Disease, Alzhiemer's Disease, and dementia with Lewy bodies. Throughout this internship, I have been able to learn about the intersection of genetics, disease, and data science, using R and other genomics tools to analyze genetic mutation data and make discoveries that have an impact on people's lives.

Because of this introduction I've had to the field of genetics annd bioinformatics, I have been curious what *other* labs and groups do with similar data. This lead me to browse through the [cBioPortal for Cancer Genomics](https://www.cbioportal.org), which I had heard about from a series of online [Dataviz + Cancer Microlabs](https://apply.hub.ki/cancerplusviz/) that occured earlier this year. The dataset I had found was used in a recent *Nature Genetics* article, "[Cancer therapy shapes the fitness landscape of clonal hematopoiesis](https://www.nature.com/articles/s41588-020-00710-0)" so I decided this would be perfect for my final project.

Question
--------

#### Background

The background of this article follows first author Kelly Bolton, a physician scientist (MD-PhD) and first year fellow in medical oncology at Memorial Sloan Kettering, the world’s oldest and largest private cancer center. Bolton was new to studying clonal hematopoiesis (CH), but jumped on the opportunity to research it, especially to understand if specific therapy types resulted in a higher frequency of CH. With her knowledge of epidemiology, she realized that she could use data from both electronic health records and sequential samples from patients, to understand how therapy could be promoting pre-existing CH or inducing new mutations.

CH essentially happens when a hematopoietic stem cell, which can develop into different types of blood cells, starts making cells with the same genetic mutation in individuals without a blood disease. As such, clonal hematopoiesis [includes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7065480/) the entire spectrum of premalignant conditions related to somatic mutations in genes associated with myeloid disorders. CH remains benign in most people, but sometimes progresses to malignancy.

Previously, studies have shown that certain types of chemotherapy lead to a higher risk of developing myeloid neoplasms (a class of clonal hematopoietic stem cell disorders), and also that the presence of clonal hematopoiesis increases the risk of developing these therapy-related myleoid neoplasm diseases (tMNs). Thus, the relationship between specific therapies and CH is important to investigate. Understanding this relationship can aid in the understanding of the mechanisms by which therapy-related myeloid diseases occur, and from there, can help develop interventions.

#### Project

With this project, I aim to replicate this lab's findings that mutations in genes are enriched based on specific exposures. Specifically, I will set out to verify that "**mutations in *ASXL1* are enriched in current or former smokers, whereas cancer therapy with radiation, platinum and topoisomerase II inhibitors preferentially selects for mutations in DNA damage response genes (*TP53*, *PPM1D*, *CHEK2*).**"

To do so, I will use exactly the same data the authors used, with analyses using a combination of their published R code and some of my own original code. In this process, I will be thoroughly narrating my steps to give insight into the process and decisions that were made. This demonstrates the knowledge I've gained this semester in DATA-202, and that I am equipped to practically apply this in new situations.

> *Much of the technical background information, including specific phrases, from this section was found in Bolton's [Nature Research blog post](https://cancercommunity.nature.com/posts/cancer-therapy-shapes-the-fitness-landscape-of-clonal-hematopoiesis).*

Dataset
-------

An analysis of the appropriateness of your dataset for addressing these questions.

A brief (2-4 sentences) verbal description of the dataset: what is the dataset about?

-   **High-level dataset description:** Scholarly-level background about the data is provided. *(more than "Data Clearly Described")*

A description of its provenance: where did the data come from originally? Where did you download it from? And (as much as you can tell or speculate) how did it end up available there?

-   **Data provenance**: Thoughtful consideration of provenance, e.g., selecting between different datasets based on provenance considerations, or drawing implications for results based on provenance
    -   more than: *Answers basic provenance questions: where did the data come from originally? Where did you download it from? And (as much as you can tell or speculate) how did it end up available there? What does all of this say about whether the data is reliable for addressing the real-world questions of interest? (2 points)*

The number of records in the dataset, and what each one represents

A list of the features in the dataset and their types

-   **Data strcuture**: Data structure clearly identified: number of records, what each record is, what each feature is and what type it is, with examples (2 points)
    -   \*\*WITH:\* ... with thoughtful and critical reflection, e.g., consideration of the appropriateness of different data types (full 3 points)
-   **Approach** *(not entirely sure where this fits in)*: Approach clearly stated and thoughtfully motivated by the intersection of the real-world question and the available data
    -   more than: *Approach stated clearly in terms that map clearly to the data used*
-   **Data Appropriateness**: Particularly thoughtful consideration of data appropriateness
    -   more than: *Reflected on how the provenance and structure of the data make it appropriate (or not appropriate) for the task*

(?) Data Wrangling
------------------

\[?\] *Not included explictly on outline, but I think this section would be helpful...*

-   Thoughtful reflection on choices made in data wrangling and/or their implications on subsequent analysis
    -   more than: *Correctly performed needed data wrangling, e.g., merging multiple data sources, aggregating, re-coding data, identifying and dealing with missing data*

Exploratory Data Analysis
-------------------------

Show plots illustrating the distribution of at least 4 variables in your dataset. Comment on anything interesting you observe.

-   **Univariate EDA**: Includes insightful / thoughtful commentary on implications of the distributions plotted (perhaps how that informs further analysis)
    -   more than: *Includes plots illustrating the distribution of more than one variable in the dataset, with some commentary.*

Show plots illustrating bivariate relationships for at least 2 pairs of variables. Explain what you observe (e.g., positive/negative correlation, no correlation, etc.).

-   **Bivariate EDA**: Includes insightful / thoughtful commentary on the implication of the relationships plotted
    -   more than: *Includes plots illustrating bivariate relationships for at least 2 pairs of variables and some description of observations (e.g., strength of relationship, dependence on other factors).*

Modeling
--------

Clearly state what is the target variable you are trying to predict, which variables (features) you are using to predict it, and why you chose those features.

-   **Modeling setup/ task description**: Modeling task is both clearly described and well motivated (e.g., alternative decisions are identified and considered)
    -   more than: *Clear description of the model setup. e.g., for a predictive model, what are the features and targets, why those were chosen, what validation method was employed and why, etc.*

Fit a basic predictive model using one of the techniques we discussed in class (regression: Nearest Neighbors or Linear Regression, classification: Nearest Neighbors or Logistic Regression; other choices such as Decision Trees are also fine)

Report the results of your basic predictive model on held-out data or via cross-validation.

-   **Modeling baselie results**: A baseline was thoughtfully chosen and evaluated.
    -   more than: *Results of a correctly applied baseline model are reported correctly*

Make one or more changes to the predictive model to improve the accuracy. Discuss what changes you made, why you made them, and what the results were.

-   **Modeling refinements**: Particular thought or insight was given to the choice of changes to make and/or analysis of their implications.
    -   more than: *Reports what changes were made, why they were made, and what the results were.*

**Alternative**: instead of a supervised prediction task, you can define an unsupervised learning task and use clustering. In this case, clearly state what you want to understand through the clustering, and report your observations.

\[FIXTHIS: pick up on rubric here... **Discussion of findings**\]

Findings:
---------

Summarize the analyses you performed and what the results told you. What do your findings say about the real-world and prediction (or clustering) questions you posed?

Limitations:
------------

What are some limitations of your analyses and potential biases of the data you used?

Future Directions:
------------------

What new questions came up following your exploration of this data? Describe at least one question that could not be answered using your data alone, and specify what additional data you would collect to address it.

Conclusions
-----------

> this is my own added section. see if I want to include this or not. if not, put this paragraph somewhere else

In addressing this question, I gained not only an understanding of how data science is used in my chosen career field, but also an appreciation for the hard work done by women in STEM -- as, to my surprise, first author Kelly Bolton MD PhD, as well as lab PI, Elli Papaemmanuil, PhD, are both women.
