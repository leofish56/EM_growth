###Order of the file####
To ensure optimal code performance, we recommend first cleaning your data (clean_data.R), then executing the function (final_function.R), and finally plotting your results (plots.R).
###Dataset format:###
The function has been created for a specific type of dataset. This dataset are available at https://shiny.marine.ie/igfs/.
###Context:###
This project regroup the R code used for my master thesis on the implications of omitting unclassified individuals from sex-specific growth mode. 
The codes as been created following methodology of the article Minto et al. (2018) (DOI:10.1139/cjfas-2016-0450)
###Abstarct of the thesis:###
To fit sex-specific growth model, the sex of the fish must apparently be known for each observation. However, in some species, the sex of juveniles cannot be visually identified 
generating a bias in the growth fit at early ages. In this study, a model using a latent variable model is tested on several species. The model estimates sex-specific von Bertalanffy 
growth model using an EM algorithm to include the individuals having unclassified sex observations. The species analysed are exploited species having sexually dimorphic growth and 
partially unclassified observations. Significant differences in parameter estimates were found between the traditional approach of omitting unclassified observations and the EM algorithm.
The parameters K and t_0 of the von Bertalanffy growth model appear to be the most affected. These parameters are the one affecting the most the fit of the growth at young ages. These 
results show the potential of the model and the implication of omitting unclassified sex observation for fisheries management. However, further research must be done for the improvement
and generalisation of the model; investigating the unclassified individuals, studying the difference between immature and mature growth, implementing the EM algorithm on other growth 
models, and further research in statistical tests.
