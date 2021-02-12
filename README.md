# SHM_HCV_incidence

The following documents in this reprository will allow you to conduct a similar analysis as the one found in the following article: 

HCV micro-elimination in individuals with HIV in the Netherlands 4 years after universal access to direct-acting antivirals: a retrospective cohort study.
Smit C, Boyd A, Rijnders BJA, van de Laar TJW, Leyten EM, Bierman WF, Brinkman K, Claassen MAA, den Hollander J, Boerekamps A, Newsum AM, Schinkel J, Prins M, Arends JE, Op de Coul ELM, van der Valk M, Reiss P; ATHENA observational cohort.
Lancet HIV. 2021 Feb;8(2):e96-e105. doi: 10.1016/S2352-3018(20)30301-5. 


Given the highly confidential nature of our data, we cannot publically release the dataset used during the analysis of this article. Instead, I have generated a simulated dataset with all the essential variables used in analysis. This dataset was generated using Stata (v15.1, College Station, TX) and the code can be found in the do-file entitled: 

smit_hcv_incidentie_data_sample.do

This do-file provides step-by-step instructions of how variables were constructed. Inputs from the ATHENA dataset were partially used to generate variables, meaning that the assumptions I used here should NOT be used to get back to the original results in the manuscript. The content of these variables are only to be used for illustrative purposes.

I cannot stress enough that this is a simulated dataset. Variables were simulated to the minimum assumptions needed to loosely replicate the statistical models and figures found in the manuscript. Any use of this dataset to answer research questions other than those addressed in the manuscript is strongly discouraged (unless the author is 100% certain of what they are doing). 


The analysis for HCV primary infection is given in the following do-file: 

smit_hcv_incidentie_analyse_sample.do

The do-file provides step-by-step instructions of the procedures and models used to generate the more complex figures in the manuscript. Anyone should be able to run the code from start to finish and to replicate a set of figures similar to those in the manuscript. 


NOTE: I have only provided simulated data for primary incidence calculations. The exact same models were used to calculate incidence of HCV re-infection. (The only difference is that the HCV re-infection dataset will contain indiviudals who re-enter follow-up after SVR or spontaneous clearance, see the manuscript for specific details on the left-censoring used for individuals with re-re-infection.) 


NOTE: As with anything in life, mistakes can happen. The goal is to get the correct informaton out there. I am hoping that providing this code as open source will help others provide input on ways in which this analysis could be improved or if any mistakes occurred. I will ensure that any errors be checked in the original code and if present, be reported to the editors of the Lancet HIV in due course. You can contact me at the following email address: a.c.boyd@amsterdamumc.nl
