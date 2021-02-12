/***********************************************
Re: Program to produce simulated data to mimic 
analysis published in Smit C, Boyd A, et al. the 
Lancet HIV, 2021
(only concerns primary HCV infection)

Author: Anders Boyd
Date created: 29 September 2020 
************************************************/


					/* NOTE: If you havent read the READ ME file accompanying the repo 
							(or more importantly, if you havent read the article), then dont bother 
							looking at this do-file. */
							
					/* NOTE: All binomial variables in the dataset are coded as: 
							1 = positive/yes, 0 = negative/no */
							

							
clear all
set seed 218423


							/* insert location of file */ 
cd ""


	/* setting number of individuals included */ 
set obs 23590


	/* patient id number */ 
gen s = "S"
gen n = _n 
egen id = concat(s n)
drop s n
	

	/* generating key group membership */ 
gen group = rbinomial(1, 0.63)							/* 1 = MSM */ 
replace group = 2*rbinomial(1, 0.95) if group == 0   	/* 2 = PWID */ 
replace group = 3 if group == 0							/* 3 = Hetero/Other */ 


	/* setting up follow-up */ 
				/* generating start year */ 
gen start_2000 = rbinomial(1, 0.35)	/* roughly 35% of the cohort was already included in 2000 */ 
gen year_start = 2000 if start_2000 == 1
replace year_start = 2000 + round(rbeta(1, 1.25)*19, 1) if start_2000 == 0 /* inclusions based on Beta(1, 1.25) distribued rate after 2000 */
tab year_start

				/* generating last year */ 
gen fu_complete = rbinomial(1, 0.8) if year_start == 2000	/* roughly 80% of the cohort included in 2000 is still in follow-up */ 
foreach y of numlist 2001/2019 {							/* loss to follow-up (LTFU) decreased at a constant rate depending on year of inclusion */   	
replace fu_complete = rbinomial(1, 0.8 + (1-((2020-`y')/20))*0.2) if year_start == `y'	
}
gen year_end = 2019 if fu_complete == 1						/* assigning year based on whether individual is LTFU */ 
foreach y of numlist 2000/2017 {
replace year_end = runiformint(`y'+1, 2018) if fu_complete == 0 & year_start == `y'
}
replace year_end = 2018 if fu_complete == 0 & year_start == 2018

drop start_2000 fu_complete


	/* expanding data for follow-up */ 
gen n = year_end - year_start
foreach n of numlist 1/19 {
expand `n'+1 if n == `n'
}
drop n 

sort id, stable
by id: gen year = year_start + (_n - 1)


	/* generating variable for person-years of observation */ 
gen pyo = 1 if year > year_start & year < year_end  
replace pyo = runiform(1, 365)/365 if year == year_start
replace pyo = runiform(1, 365)/365 if year == year_end & year_end != 2019
replace pyo = rigaussian(0.83, 8) if year == year_end & year_end == 2019
replace pyo = 1 if pyo > 1 & year_end == 2019


	
	/* generating tested during the year */ 
gen hcv_test = rbinomial(1, 0.20 + runiform(-0.20, 0.20)) if group == 1 & year <= 2006		/* testing in MSM */ 
replace hcv_test = rbinomial(1, 0.28 + runiform(-0.20, 0.20)) if group == 1 & year == 2007
replace hcv_test = rbinomial(1, 0.37 + runiform(-0.20, 0.20)) if group == 1 & year > 2007 & year <= 2016
replace hcv_test = rbinomial(1, 0.40 + runiform(-0.20, 0.20)) if group == 1 & year > 2016 
replace hcv_test = rbinomial(1, 0.20 + runiform(-0.20, 0.20)) if group != 1 				/* testing for all other groups */ 
								/* NOTE: this simulation assumes random distribution of HCV testing moments within individuals -  
										considering certain individuals will be more frequently engaging in risk activities 
										associated with HCV infection, they will likely be more often tested, which then violates
										this assumption. */ 
								/* NOTE: the "runiform(-0.20, 0.20)" bit was added as a noise factor to bring some variability into the mix. */

								
	/* generating covariates for the model: age, HIV-RNA and CD4 markers */ 
								/* NOTE: these variables were created before the test results because this endpoint 
										is dependent on age. */
					/* age */ 
gen birth_y = .
foreach y of numlist 2000/2019 /* birth year needs to be generated per year coming into the cohort */ {
by id: replace birth_y = round(rnormal(`y'-38 + runiform(-5, 5), 10), 1) if year_start == `y' & _n == 1
by id: replace birth_y = `y'-18 if birth_y > `y'-18 & year_start == `y' & _n == 1   /* individuals had to be >18 years old to be included */					
}	
by id: replace birth_y = birth_y[_n-1] if birth_y == . 
gen age = year - birth_y


					/* HIV-RNA viral load */
gen hivrna_udt = .
foreach y of numlist 2000/2019 {
by id: replace hivrna_udt = rbinomial(1, 0.50 + (.015*(`y'-2000))) if year == `y' & year < 2015
by id: replace hivrna_udt = rbinomial(1, 0.80 + runiform(-0.2, 0.2)) if year == `y' & year >= 2015				
}
gen hivrna_log = rnormal(5.00, 0.08) if hivrna_udt == 0 & year < 2015
replace hivrna_log = rnormal(2.30, 0.03) if hivrna_udt == 0 & year >= 2015
replace hivrna_log = log10(50) if hivrna_udt == 1
								/* NOTE: Simulations were done based on the distributions at each year 
										and were not based on individual trajectors. As a result, some of 
										the patients will have very unrealistic blips or patterns in HIV RNA 
										replication. Normally it should not matter since HCV incidence 
										did not seem to change at different levels of HIV RNA viral loads. */ 
										
									
					/* CD4+ count */ 
gen cd4 = .
foreach y of numlist 2000/2019 {
by id: replace cd4 = round(rnormal(326 + (15*(`y'-2000)), 250), 1) if year == `y'				
}
by id: replace cd4 = 0 if cd4 < 0 
								/* NOTE: Again, simulations were done based on the distributions at each year 
										and were not based on individual trajectors. As a result, some of 
										the patients will have very unrealistic increases / descreases in 
										CD4 cell count. Normally it should not matter since HCV incidence 
										did not seem to change at different levels of CD4 cell counts. */ 


	
	/* generating test result */ 
								/* NOTE: This is quite painful to simulate with multivariate distributions 
										(since HCV incidence is in function of transmission group, year, 
										whether someone got tested and age), so I simulated the HCV test 
										result outcomes for certain easily identifiable strata. */
								/* NOTE: The positive test result is loosely simulated from the incidence 
										observed in the stratum (after the if statement) divided by the proprotion 
										of individuals tested (as only tested individuals can have a result). */ 
gen event = rbinomial(1, 0.0005/0.20) if year == 2000 & hcv_test == 1	& group == 1 & age <= 30
replace event = rbinomial(1, 0.002/0.20) if year == 2001 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.002/0.20) if year == 2002 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.002/0.20) if year == 2003 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.002/0.20) if year == 2004 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.002/0.20) if year == 2005 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.005/0.20) if year == 2006 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.010/0.28) if year == 2007 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.014/0.37) if year == 2008 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.006/0.37) if year == 2009 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.016/0.37) if year == 2010 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.008/0.37) if year == 2011 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.008/0.37) if year == 2012 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.015/0.37) if year == 2013 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.015/0.37) if year == 2014 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.012/0.37) if year == 2015 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.003/0.37) if year == 2016 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.003/0.40) if year == 2017 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.003/0.40) if year == 2018 & hcv_test == 1 & group == 1 & age <= 30
replace event = rbinomial(1, 0.003/0.40) if year == 2019 & hcv_test == 1 & group == 1 & age <= 30

replace event = rbinomial(1, 0.002/0.20) if year == 2000 & hcv_test == 1 & group == 1 & age > 30 & age <= 55 
replace event = rbinomial(1, 0.002/0.20) if year == 2001 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.002/0.20) if year == 2002 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.002/0.20) if year == 2003 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.002/0.20) if year == 2004 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.002/0.20) if year == 2005 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.005/0.20) if year == 2006 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.017/0.28) if year == 2007 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.014/0.37) if year == 2008 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.011/0.37) if year == 2009 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.011/0.37) if year == 2010 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.015/0.37) if year == 2011 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.014/0.37) if year == 2012 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.011/0.37) if year == 2013 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.012/0.37) if year == 2014 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.012/0.37) if year == 2015 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.003/0.37) if year == 2016 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.003/0.40) if year == 2017 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.003/0.40) if year == 2018 & hcv_test == 1 & group == 1 & age > 30 & age <= 55
replace event = rbinomial(1, 0.003/0.40) if year == 2019 & hcv_test == 1 & group == 1 & age > 30 & age <= 55

replace event = rbinomial(1, 0.0001/0.20) if year == 2000 & hcv_test == 1 & group == 1 & age > 55 
replace event = rbinomial(1, 0.0001/0.20) if year == 2001 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0001/0.20) if year == 2002 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0002/0.20) if year == 2003 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0001/0.20) if year == 2004 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0001/0.20) if year == 2005 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0003/0.20) if year == 2006 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.004/0.28) if year == 2007 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.005/0.37) if year == 2008 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.003/0.37) if year == 2009 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.003/0.37) if year == 2010 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.005/0.37) if year == 2011 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.006/0.37) if year == 2012 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.004/0.37) if year == 2013 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.003/0.37) if year == 2014 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.006/0.37) if year == 2015 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0001/0.37) if year == 2016 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0001/0.40) if year == 2017 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0001/0.40) if year == 2018 & hcv_test == 1 & group == 1 & age > 55
replace event = rbinomial(1, 0.0001/0.40) if year == 2019 & hcv_test == 1 & group == 1 & age > 55

replace event = rbinomial(1, 0/0.20) if year == 2000 & hcv_test == 1	& group == 2
replace event = rbinomial(1, 0/0.20) if year == 2001 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2002 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0.001/0.20) if year == 2003 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2004 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0.00005/0.20) if year == 2005 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0.0003/0.20) if year == 2006 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0.0001/0.20) if year == 2007 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0.00005/0.20) if year == 2008 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2009 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0.001/0.20) if year == 2010 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2011 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2012 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2013 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0.001/0.20) if year == 2014 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2015 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2016 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2017 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2018 & hcv_test == 1 & group == 2
replace event = rbinomial(1, 0/0.20) if year == 2019 & hcv_test == 1 & group == 2	

replace event = rbinomial(1, 0.0024/0.20) if year == 2000 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0015/0.20) if year == 2001 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0008/0.20) if year == 2002 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0015/0.20) if year == 2003 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0006/0.20) if year == 2004 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0006/0.20) if year == 2005 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0008/0.20) if year == 2006 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0012/0.20) if year == 2007 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0010/0.20) if year == 2008 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0010/0.20) if year == 2009 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0014/0.20) if year == 2010 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0012/0.20) if year == 2011 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0008/0.20) if year == 2012 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0009/0.20) if year == 2013 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0016/0.20) if year == 2014 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0011/0.20) if year == 2015 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0010/0.20) if year == 2016 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0010/0.20) if year == 2017 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0007/0.20) if year == 2018 & hcv_test == 1 & group == 3
replace event = rbinomial(1, 0.0008/0.20) if year == 2019 & hcv_test == 1 & group == 3				

by id: replace event = 0 if event[1] == .			/* last observation carried forward */ 
by id: replace event = event[_n-1] if event == . 

by id: gen t_y = pyo 								/* individuals with an event are right-censored */
by id: replace t_y = t_y + t_y[_n-1] if _n != 1
stset t_y, f(event == 1) id(id)						 
drop if _d == .
drop _* t_y
													
replace pyo = runiform(1, 365)/365 if event == 1	/* need to redress follow-up at last visit / HCV positive result */





	/* generating treatment data */ 
					/* ever treated */ 
gen tx = rbinomial(1, 0.8 + runiform(-0.05, 0.10)) if year < 2017 & event == 1
replace tx = rbinomial(1, 0.6 + runiform(-0.25, 0)) if year >= 2017 & event == 1


					/* time from HCV diagnosis to treatment */
gen hcvtx_delai = rigaussian(45, 10) if year == 2005 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(44, 10) if year == 2006 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(39, 10) if year == 2007 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(39, 8) if year == 2008 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(23, 8) if year == 2009 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(18, 8) if year == 2010 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(12, 6) if year == 2011 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(11, 6) if year == 2012 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(14, 6) if year == 2013 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(12, 6) if year == 2014 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(7, 2) if year == 2015 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(6, 2) if year == 2016 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(5, 2) if year == 2017 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = rigaussian(4, 2) if year == 2018 & group == 1 & event == 1 & tx == 1
replace hcvtx_delai = 1 if hcvtx_delai < 1

	
					/* year of being treated with HCV */ 
gen year_hcvtx = round(year + hcvtx_delai/12, 1) 	

													/* need to correct for those with treatment years > 2018 */ 
replace hcvtx_delai = (2018 - year)*12 if year_hcvtx > 2018 & year_hcvtx != .
replace year_hcvtx = round(year + hcvtx_delai/12, 1) 	
	

					/* treatment uptake */	
					
								/* within 12 months */ 
gen hcvtx_uptake12 = 1 if hcvtx_delai < 13 & hcvtx_delai != . 
replace hcvtx_uptake = 0 if (hcvtx_delai >= 13 & hcvtx_delai != .) | tx == 0

								/* within 6 months */ 
gen hcvtx_uptake6 = 1 if hcvtx_delai < 7 & hcvtx_delai != . 
replace hcvtx_uptake6 = 0 if (hcvtx_delai >= 7 & hcvtx_delai != .) | tx == 0



					/* sustained virological response */
gen svr = rbinomial(1, 0.5) if year_hcvtx == 2005 
replace svr = rbinomial(1, 0.5) if year_hcvtx == 2006 
replace svr = rbinomial(1, 0.55) if year_hcvtx == 2007 
replace svr = rbinomial(1, 0.74) if year_hcvtx == 2008 
replace svr = rbinomial(1, 0.70) if year_hcvtx == 2009 
replace svr = rbinomial(1, 0.70) if year_hcvtx == 2010 
replace svr = rbinomial(1, 0.75) if year_hcvtx == 2011 
replace svr = rbinomial(1, 0.80) if year_hcvtx == 2012 
replace svr = rbinomial(1, 0.79) if year_hcvtx == 2013 
replace svr = rbinomial(1, 0.98) if year_hcvtx == 2014 
replace svr = rbinomial(1, 0.98) if year_hcvtx == 2015 
replace svr = 1 if year_hcvtx == 2016 
replace svr = 1 if year_hcvtx == 2017 
replace svr = 1 if year_hcvtx == 2018 




	/* dropping unnecessary variables */ 
drop birth_y year_start year_end tx
	
	
	
	/* saving dataset */ 
save "data_primary_inc.dta", replace	

				
/****************** end of do-file ******************/
