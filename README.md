# human-footprint-index-VBD
> Code for: *Human footprint is associated with shifts in assemblages of major vector-borne diseases.* 2022. Skinner EB, Glidden CK, MacDonald AJ, Mordecai EA. *Nature Sustainability*: accepted.
> Code by Caroline K. Glidden & Eloise B. Skinner

This github contains the following R files:

1. spatiotemporal_cv.R : code to perform spatiotemporal cross-validation for each model included in the manuscript
2. bootstrapping_pdps.R: code to build partial dependence plots and calculate inflection points, with variablility in pdps and confidence intervals calculated via bootstrapping
3. variable_importance.R: code to build variable importance plots and tables for each model in the manuscript

The data (full_data.csv) includes the following columns: 

* Incidence for the vector-borne disease analyzed in the manuscript
  * Dengue, Chikungunya, Zika virus, cutaneous leishmaniasis, visceral leishmaniasis incidence was downloaded from the Brazilian national disease surveillance system (SINAN: http://portalsinan.saude.gov.br/)
  * Malaria incidence was downloaded from the Brazil Epidemiological Surveillance Information System for Malaria (SIVEP-MALARIA: https://public.tableau.com/app/profile/mal.ria.brasil/viz/Dadosparacidado_201925_03_2020/Incio)
  * Human footprint index was downloaded from Williams et al. (2020) (https://doi.org/10.1016/j.oneear.2020.08.009) and Keys et al. 2021 (https://doi.org/10.1088/1748-9326/abe00a)
  * Climate data was downloaded from the Climate Research Unit: https://crudata.uea.ac.uk/cru/data/hrg/
  * Land class data was downloaded from MAPBIOMAS: https://mapbiomas.org/
 
 Brazil municipality shapefiles for data cleaning and building figures were downloaded from: https://www.ibge.gov.br/en/geosciences/territorial-organization/territorial-meshes/18890-municipal-mesh.html?=&t=acesso-ao-produto 
  

