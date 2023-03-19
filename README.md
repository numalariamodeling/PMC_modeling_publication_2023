## Perennial Malaria Chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis

### Objective and summary
Perennial malaria chemoprevention (PMC), previously known as intermittent preventive treatment in infants, aims to reduce malaria burden during the first few years of life. The intervention has been recently relabeled by WHO to allow for more flexibility in deployment and several pilot implementation studies are ongoing to gather evidence needed for policy adoption in countries. 
The recent endorsement of the malaria vaccine RTS,S further adds to the complexity of malaria prevention options and combinations in young children that need to be addressed.

In this work, we used a long-standing individual based malaria model ([EMOD](https://docs.idmod.org/projects/emod-malaria/en/latest/index.html)) to simulate the impact of varying deployment schemes of PMC with or without malaria vaccination on malaria burden in young children across different transmission settings. 
We incorporated two coverage scenarios, one at a target coverage of 80%, and n operational coverage corresponding to Southern Nigeria.  

We find that the most impactful timing and number of PMC doses highly depends on the disease metric and transmission intensity. 
The results further highlight the complementary effects of PMC and the malaria vaccine as their combination is more effective than either intervention alone. 
Similarly, PMC with more doses is more impactful and at operational coverage additional doses might be considered to offset lower coverage.

Our analysis provides a simplified approach for exploring varying PMC schedules across different transmission settings. 

----------------------------------
_Manuscript under review, analysis scripts can be found in the branch `prepublication`_

----------------------------------

### Repository structure:  
1. `data_files` stores input and reference data used in analysis*
2. `simulation_inputs` stores scenario csv files *
3. `simulations` script to run EMOD simulations (requires HPC)
4. `simulation_output` and `postprocessed_output`: processed EMOD ouput files in csv format*
5. `analysis`: generates result figures, written into `figures`

*) Data and simulation files stored on Zenodo (DOI: 10.5281/zenodo.7750298)

### Requirements:

**EMOD**
https://docs.idmod.org/projects/emod-malaria/en/latest/index.html

**Python**   
dtk==0.2  
malaria.egg==info  
numpy==1.19.5  
pandas==1.1.5  
requests==2.25.1  
scipy==1.5.4  
simtools==0.2  
simulations==0.1.1.dev36  


**R - RVersion: 4.0.2**
dplyr  
data.table  
wrapr  
dplyr  
tidyr  
lubridate  
zoo  
ggplot2  
ggthemes  
cowplot  
scales  
RColorBrewer  
see  
pracma  



