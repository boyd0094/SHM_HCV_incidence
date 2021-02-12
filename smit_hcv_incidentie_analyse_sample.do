/***********************************************
Re: Sample program to mimic analysis published 
in Smit C, Boyd A, et al. the Lancet HIV, 2021
(only concerns primary HCV infection)

Database: Simulated 
Author: Anders Boyd
Date created: 29 September 2020 
************************************************/


					/* NOTE: If you havent read the READ ME file accompanying the repo 
							(or more importantly, if you havent read the article), then dont bother 
							looking at this do-file. */



/* specifications */ 

set more off
capture log close

							/* insert location of file */ 
cd ""





/* data */ 

					/* NOTE: The data used in this do-file are simulated from the inputs found in the article. 
							IN NO WAY do they represent ANY individual in the ATHENA cohort. Hence, use of these 
							data are not meant for any secondary analyses or individual-level meta-analysis. 
							Your use of these data are entirely your responsiblity. */ 

use "data_primary_inc.dta", clear 		






/* Figure 1A. Testing per calendar year */ 
			
		/* HCV testing during the year */
tabstat hcv_test, by(year) save				

gen g_year = . 								/* extracting information to be used in graph */
foreach n of numlist 1/20 {
replace g_year = 1999+`n' if _n == `n'
}

gen g_mean1 = .
tabstat hcv_test if group == 1, by(year) save /* getting basic information - MSM only */ 
foreach n of numlist 1/20 {
mat def S`n' = r(Stat`n')
replace g_mean1 = S`n'[1,1] if _n == `n'
}											/* modeling endpoint with year as cubic splines */ 
mkspline year_sp = year if year >= 2000 & group == 1, nk(5) cubic				
glm hcv_test year_sp* if year >= 2000 & group == 1, f(bi) l(logit)
											/* obtaining prediction */
predict yhat_1 if year >= 2000 & group == 1, xb
replace yhat_1 = invlogit(yhat_1)			/* need to bring it back to probabilities */ 

gen g_mean2 = .								/* same procedure with PWID */
tabstat hcv_test if group == 2, by(year) save
foreach n of numlist 1/20 {
mat def S`n' = r(Stat`n')
replace g_mean2 = S`n'[1,1] if _n == `n'
}
drop year_sp*								
mkspline year_sp = year if year >= 2000 & group == 2, nk(5) cubic	
glm hcv_test year_sp* if year >= 2000 & group == 2, f(bi) l(logit)
predict yhat_2 if year >= 2000 & group == 2, xb	 
replace yhat_2 = invlogit(yhat_2)			

gen g_mean3 = .								/* same procedure with hetero/other */
tabstat hcv_test if group == 3, by(year) save
foreach n of numlist 1/20 {
mat def S`n' = r(Stat`n')
replace g_mean3 = S`n'[1,1] if _n == `n'
}
drop year_sp*
mkspline year_sp = year if year >= 2000 & group == 3, nk(5) cubic	
glm hcv_test year_sp* if year >= 2000 & group == 3, f(bi) l(logit)
predict yhat_3 if year >= 2000 & group == 3, xb
replace yhat_3 = invlogit(yhat_3)

											/* plotting */ 
twoway (mspline yhat_1 year if year >= 2000, lcolor(black)) ///
		(mspline yhat_2 year if year >= 2000, lcolor(black) lpattern(longdash)) ///
		(mspline yhat_3 year if year >= 2000, lcolor(black) lpattern(longdash_dot)) ///
		(scatter g_mean1 g_year if g_year >= 2000, mcolor(black) msize(medlarge) msymbol(circle_hollow)) ///
		(scatter g_mean2 g_year if g_year >= 2000, mcolor(black) msize(medlarge) msymbol(triangle_hollow)) ///
		(scatter g_mean3 g_year if g_year >= 2000, mcolor(black) msize(medlarge) msymbol(square_hollow)), ///
				ytitle("Proportion with" "HCV testing during the year") ///
				ytitle(, size(large)) ///
				ylabel(0(0.2)1, labsize(large) angle(horizontal) format(%2.1f) nogrid) ///
				xtitle(Year) xtitle(, size(large)) ///
				xlabel(2001(2)2019, labsize(large) angle(forty_five)) ///
				legend(order(1 "MSM" 2 "PWID" 3 "Hetero/Other" 4 "" 5 "" 6 "") rows(2) size(large)) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))		
			
						/* NOTE: I used the mspline command to get smoothed functions of the 
							predicted probability from the model. There is a default number of bands
							that you can use with this type of graph. Please keep in mind that that 
							this number can drastically change the shape of the figure, especially 
							at lower numbers! */
			
			
		

		
/* Figure 2B. Total number of infections per risk group */		
		
		
		
		/* description of incident cases */
stset pyo, f(event == 1)

by id: gen last = 1 if _n == _N
tabstat _t if last == 1, stat(p25 p50 p75 n)
stsum
strate, per(1000)
strate year, per(1000)
tab year if _d == 1

	
	
	
		/* transmission risk factors at incident infection */
tab year group if _d == 1, row

drop g_*
recode group (1=1) (2/3=0), gen(group1) 
tabstat group1 if _d == 1, by(year) save
gen g_year = . 
gen g_mean1 = .
foreach n of numlist 1/20 {
replace g_year = 1999+`n' if _n == `n'
mat def S`n' = r(Stat`n')
replace g_mean1 = S`n'[1,1] if _n == `n'
}
recode group (1=0) (2=1) (3=0), gen(group2)
tabstat group2 if _d == 1, by(year) save
gen g_mean2bis = .
foreach n of numlist 1/20 {
mat def S`n' = r(Stat`n')
replace g_mean2bis = S`n'[1,1] if _n == `n'
}
recode group (1/2=0) (3=1), gen(group3)
tabstat group3 if _d == 1, by(year) save
gen g_mean3bis = .
foreach n of numlist 1/20 {
mat def S`n' = r(Stat`n')
replace g_mean3bis = S`n'[1,1] if _n == `n'
}

gen g_mean2 = g_mean1 + g_mean2bis /* to get percentages to total 1 for the bar chart */
gen g_mean3 = g_mean1 + g_mean2bis + g_mean3bis 
		
twoway (bar g_mean3 g_year, fcolor(black) lcolor(black) barwidth(0.8)) ///
		(bar g_mean2 g_year, fcolor(gs4) lcolor(black) barwidth(0.8)) ///
		(bar g_mean1 g_year, fcolor(gs12) lcolor(black) barwidth(0.8)) ///
					if g_year >= 2005, ///
				ytitle("Proportion belonging" "to risk group") ///
				ytitle(, size(large)) ///
				ylabel(0(0.2)1, labsize(large) angle(horizontal) format(%3.2f) nogrid) ///
				xtitle(Year) ///
				xtitle(, size(large)) ///
				xlabel(2005(1)2019, labsize(large) angle(forty_five)) ///
				legend(order(3 "MSM" 2 "PWID" 1 "Hetero/Other") rows(1) size(large)) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))

				
				
		
		
		
		
				
/* NB : the large majority of HCV re-infections come from MSM; this will be the focus 
		group from this point onwards */

		
keep if group == 1		 





/* Figure 3A. Incidence rate per year */ 


	/* incidence rates over time */

		/* direct calculation */

								/* declaring data to be survival time */							
stset pyo, f(event == 1)
tab event if pyo == 0 							
stsum															
tab year if _d == 1

								/* direct calculation (using this setup) */
strate, per(1000)
strate year, per(1000)






		/* fixed-effect Poisson regression */

					/* need to get indicator variables for each year */
tab year, gen(year)

								/* direct calculation (using this setup) */
table year, contents(sum pyo sum event) 

								
								/* models */
												/* Model 1: to obtain estimates on calendar year (unadjusted) */
glm event year1-year20, f(poisson) l(log) exposure(pyo) eform noconst irls

mat def b_v = r(table)
mat list b_v
drop g_year
gen g_year = . 
gen ir = .
gen ir_lb = .
gen ir_ub = .
foreach n of numlist 1/20 {
replace g_year = 1999+`n' if _n == `n'
replace ir = b_v[1,`n'] if _n == `n'
replace ir_lb = b_v[5,`n'] if _n == `n'
replace ir_ub = b_v[6,`n'] if _n == `n'
}
export delim g_year ir ir_lb ir_ub using "fig_hcv_msm_ir", delim(";") replace


save "temp.dta", replace


				
								/* plotting */	
import delim using "fig_hcv_msm_ir.csv", ///
		delimiter(";") clear

foreach v of varlist ir* {		
replace `v' = `v'*1000	
}

twoway (scatter ir g_year, mcolor(black) msymbol(circle_hollow)) ///
		(rcap ir_lb ir_ub g_year, lcolor(gs6) msize(medlarge)), ///
				ytitle(Incidence (per 1000 py)) ///
				ytitle(, size(large)) ///
				ylabel(0(5)20, labsize(large) angle(horizontal) format(%2.1f) nogrid) ///
				xtitle(Year) ///
				xtitle(, size(large)) ///
				xlabel(2001(2)2019, /*angle(forty_five)*/ labsize(large)) ///
				legend(off) ///
				xsize(3.5) ///
				ysize(2) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))								
	
	
	
	
/* supplementary figure S2 */	
	
	
												/* Model 2: to obtain estimates on age */
use "temp.dta", clear
erase "temp.dta"
												
mkspline age_sp = age, nk(3) cubic
glm event year1-year20 age_sp*, f(poisson) l(log) exposure(pyo) eform noconst irls

predict yhat, xb
predict error, stdp
gen lb = yhat - invnormal(0.975)*error
gen ub = yhat + invnormal(0.975)*error				
gen plb = 1000*invlogit(lb)
gen pub = 1000*invlogit(ub)				
replace yhat = 1000*invlogit(yhat)

twoway (mspline yhat age, bands(10) lcolor(black)) ///
		(mspline plb age, bands(10) lcolor(black) lpattern(longdash)) ///
		(mspline pub age, bands(10) lcolor(black) lpattern(longdash)), ///
				ytitle("Incidence (per 1000 person-years)") ///
				ytitle(, size(large)) ///
				ylabel(0(2)8, labsize(large) angle(horizontal) format(%2.1f) nogrid) ///
				xtitle(Age (years)) xtitle(, size(large) margin(medsmall)) ///
				xlabel(20(10)90, labsize(large)) ///
				legend(off) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))

			

												/* Model 3: to obtain estimates on hiv rna viral load and cd4 count */
glm event year1-year20 hivrna_log cd4, f(poisson) l(log) exposure(pyo) eform noconst irls

drop yhat error lb ub plb pub
predict yhat, xb
predict error, stdp
gen lb = yhat - invnormal(0.975)*error
gen ub = yhat + invnormal(0.975)*error				
gen plb = 1000*invlogit(lb)
gen pub = 1000*invlogit(ub)				
replace yhat = 1000*invlogit(yhat)

twoway (mspline yhat cd4, bands(8) lcolor(black)) ///
		(mspline plb cd4, bands(8) lcolor(black) lpattern(longdash)) ///
		(mspline pub cd4, bands(8) lcolor(black) lpattern(longdash)), ///
				ytitle("Incidence (per 1000 person-years)") ///
				ytitle(, size(large)) ///
				ylabel(0(2)10, labsize(large) angle(horizontal) format(%2.1f) nogrid) ///
				xtitle(CD4+ cell count (/mm3)) xtitle(, size(large) margin(medsmall)) ///
				xlabel(200(250)1500, labsize(large)) ///
				legend(off) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))

twoway (mspline yhat hivrna_log, bands(8) lcolor(black)) ///
		(mspline plb hivrna_log, bands(8) lcolor(black) lpattern(longdash)) ///
		(mspline pub hivrna_log, bands(8) lcolor(black) lpattern(longdash)), ///
				ytitle("Incidence (per 1000 person-years)") ///
				ytitle(, size(large)) ///
				ylabel(0(3)15, labsize(large) angle(horizontal) format(%2.1f) nogrid) ///
				xtitle(HIV RNA (log10 copies/mL)) xtitle(, size(large) margin(medsmall)) ///
				xlabel(2(1)5, labsize(large) format(%2.1f)) ///
				legend(off) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))

								
								
		/* mean age at incident infection */
drop g_year g_mean*
tabstat age if _d == 1, stat(p25 p50 p75 n mean sd) by(year) save
						
								/* inputting data for figure */	
gen g_year = . 
gen g_mean = .
foreach n of numlist 1/19 {
replace g_year = 2000+`n' if _n == `n'
mat def S`n' = r(Stat`n')
replace g_mean = S`n'[5,1] if _n == `n'
}

								/* figure with fitted versus observed average age */
drop year1 year2 year3 year4 yhat error lb ub
mkspline year = year if year >= 2005, nk(5) cubic
glm age year1 year2 year3 year4 if _d == 1 & year >= 2005, f(gau) l(id)
predict yhat if _d == 1, xb
predict error if _d == 1, stdp
gen lb = yhat - invnormal(0.975)*error
gen ub = yhat + invnormal(0.975)*error

twoway (mspline yhat year if _d == 1 & year >= 2005, lcolor(black)) ///
		(mspline lb year if _d == 1 & year >= 2005, lcolor(black) lpattern(longdash)) ///
		(mspline ub year if _d == 1 & year >= 2005, lcolor(black) lpattern(longdash)) ///
		(scatter g_mean g_year if g_year >= 2005, mcolor(black) msize(medlarge) msymbol(circle_hollow)), ///
				ytitle(Age (years)) ///
				ytitle(, size(large)) ///
				ylabel(20(5)50, labsize(large) angle(horizontal) nogrid) ///
				xtitle(Year) xtitle(, size(large)) ///
				xlabel(2005(1)2019, labsize(large) angle(forty_five)) ///
				legend(off) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))

				
				
		/* mean cd4 at incident infection */
drop g_year g_mean		
tabstat cd4 if _d == 1, stat(p25 p50 p75 n mean sd) by(year) save 

gen g_year = . 
gen g_mean = .
foreach n of numlist 1/19 {
replace g_year = 2000+`n' if _n == `n'
mat def S`n' = r(Stat`n')
replace g_mean = S`n'[5,1] if _n == `n'
}
		
glm cd4 year1 year2 year3 year4 if _d == 1 & year >= 2005, f(gau) l(id)
drop yhat error lb ub
predict yhat if _d == 1, xb
predict error if _d == 1, stdp
gen lb = yhat - invnormal(0.975)*error
gen ub = yhat + invnormal(0.975)*error

twoway (mspline yhat year if _d == 1 & year >= 2005, lcolor(black)) ///
		(mspline lb year if _d == 1 & year >= 2005, lcolor(black) lpattern(longdash)) ///
		(mspline ub year if _d == 1 & year >= 2005, lcolor(black) lpattern(longdash)) ///
		(scatter g_mean g_year if g_year >= 2005, mcolor(black) msize(medlarge) msymbol(circle_hollow)), ///
				ytitle(CD4+ T cell (/mm3)) ///
				ytitle(, size(large)) ///
				ylabel(300(100)800, labsize(large) angle(horizontal) nogrid) ///
				xtitle(Year) xtitle(, size(large)) ///
				xlabel(2005(1)2019, labsize(large) angle(forty_five)) ///
				legend(off) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))


		/* proportion with detectable HIV-RNA at incident infection */				
drop g_year g_mean		
tabstat hivrna_udt if _d == 1, by(year) save

gen g_year = . 
gen g_mean = .
foreach n of numlist 1/19 {
replace g_year = 2000+`n' if _n == `n'
mat def S`n' = r(Stat`n')
replace g_mean = S`n'[1,1] if _n == `n'
}
		
glm hivrna_udt year1 year2 year3 year4 if _d == 1 & year >= 2005, f(bi) l(logit)
drop yhat error lb ub plb pub
predict yhat if _d == 1, xb
predict error if _d == 1, stdp
gen lb = yhat - invnormal(0.975)*error
gen ub = yhat + invnormal(0.975)*error				
gen plb = invlogit(lb)
gen pub = invlogit(ub)				
replace yhat = invlogit(yhat)

twoway (mspline yhat year if _d == 1 & year >= 2005, lcolor(black)) ///
		(mspline plb year if _d == 1 & year >= 2005, lcolor(black) lpattern(longdash)) ///
		(mspline pub year if _d == 1 & year >= 2005, lcolor(black) lpattern(longdash)) ///
		(scatter g_mean g_year if g_year >= 2005, mcolor(black) msize(medlarge) msymbol(circle_hollow)), ///
				ytitle("Proportion with" "HIV RNA <200 copies/mL") ///
				ytitle(, size(large)) ///
				ylabel(0(0.2)1, labsize(large) angle(horizontal) format(%2.1f) nogrid) ///
				xtitle(Year) xtitle(, size(large)) ///
				xlabel(2005(1)2019, labsize(large) angle(forty_five)) ///
				legend(off) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))

			
				
				
				
					
				
				
				
				
/* Figure 4. Treatment outcomes */ 				
		
		
		/* median (IQR) time from diagnosis (primary infection) to treatment */
graph box hcvtx_delai if hcvtx_delai > 0 & year >= 2005 & year < 2019, ///
		over(year, label(angle(forty_five) labsize(large))) ///
		box(1, fcolor(gs12) lcolor(black)) ///
		marker(1, mcolor(gs8) msize(medium)) ///
		ytitle(HCV diagnosis to treatment (months)) ///
		ytitle(, size(large)) ///
		yscale(log) ///
		ylabel(1 12 24 36 60 120, labsize(large) angle(horizontal) nogrid) ///
		graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
		plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))


				
 		/* treatment uptake for primary infection */
drop g_year g_mean*		

gen g_year = . 
foreach n of numlist 1/18 {
replace g_year = 2000+`n' if _n == `n'
}

tabstat hcvtx_uptake12 if _d == 1, stat(mean n) by(year) save 
gen g_mean1 = .
foreach n of numlist 1/18 {
mat def S`n' = r(Stat`n')
replace g_mean1 = S`n'[1,1] if _n == `n'
}		
mkspline yearbis = year if year >= 2005 & year < 2019, nk(5) cubic
glm hcvtx_uptake12 yearbis1 yearbis2 yearbis3 yearbis4 if _d == 1 & year >= 2005 & year < 2019, f(bi) l(logit)
predict yhat1 if _d == 1 & year >= 2005 & year < 2019, xb
predict error1 if _d == 1 & year >= 2005 & year < 2019, stdp
gen lb1 = yhat1 - invnormal(0.975)*error1
gen ub1 = yhat1 + invnormal(0.975)*error1				
gen plb1 = invlogit(lb1)
gen pub1 = invlogit(ub1)				
replace yhat1 = invlogit(yhat1)

tabstat hcvtx_uptake6 if _d == 1, stat(mean n) by(year) save 
gen g_mean2 = .
foreach n of numlist 1/18 {
mat def S`n' = r(Stat`n')
replace g_mean2 = S`n'[1,1] if _n == `n'
}		
glm hcvtx_uptake6 yearbis1 yearbis2 yearbis3 yearbis4 if _d == 1 & year >= 2005 & year < 2019, f(bi) l(logit)
predict yhat2 if _d == 1 & year >= 2005 & year < 2019, xb
predict error2 if _d == 1 & year >= 2005 & year < 2019, stdp
gen lb2 = yhat2 - invnormal(0.975)*error2
gen ub2 = yhat2 + invnormal(0.975)*error2				
gen plb2 = invlogit(lb2)
gen pub2 = invlogit(ub2)				
replace yhat2 = invlogit(yhat2)

twoway (scatter g_mean1 g_year if g_year >= 2005 & g_year < 2019, mcolor(black) msize(medlarge) msymbol(circle_hollow)) ///
		(scatter g_mean2 g_year if g_year >= 2005 & g_year < 2019, mcolor(black) msize(medlarge) msymbol(triangle_hollow)), ///
				ytitle("Proportion treatment uptake") ///
				ytitle(, size(large)) ///
				ylabel(0(0.2)1, labsize(large) angle(horizontal) format(%2.1f) nogrid) ///
				xtitle(Year of diagnosis) xtitle(, size(large)) ///
				xlabel(2005(1)2018, labsize(large) angle(forty_five)) ///
				legend(order(1 "<=12 months" 2 "<=6 months") rows(1) size(large)) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))				


		
				
 		/* svr (in those treated) for primary infection */
drop g_*		

tabstat svr if _d == 1, stat(mean n) by(year_hcvtx) save

gen g_year = . 
gen g_mean = .
foreach n of numlist 1/15 {
replace g_year = 2004+`n' if _n == `n'
mat def S`n' = r(Stat`n')
replace g_mean = S`n'[1,1] if _n == `n'
}

mkspline year_tx = year_hcvtx if year >= 2006, nk(3) cubic
glm svr year_tx1 year_tx2 if _d == 1 & year >= 2006, f(bi) l(logit)
drop yhat error lb ub plb pub
predict yhat if _d == 1 & year >= 2006, xb
predict error if _d == 1 & year >= 2006, stdp
gen lb = yhat - invnormal(0.975)*error
gen ub = yhat + invnormal(0.975)*error				
gen plb = invlogit(lb)
gen pub = invlogit(ub)				
replace yhat = invlogit(yhat)

twoway (mspline yhat year_hcvtx if _d == 1 & year_hcvtx >= 2006 & year_hcvtx < 2019, lcolor(black)) ///
		(mspline plb year_hcvtx if _d == 1 & year_hcvtx >= 2006 & year_hcvtx < 2019, lcolor(black) lpattern(longdash)) ///
		(mspline pub year_hcvtx if _d == 1 & year_hcvtx >= 2006 & year_hcvtx < 2019, lcolor(black) lpattern(longdash)) ///
		(scatter g_mean g_year if g_year >= 2006 & g_year < 2019, mcolor(black) msize(medlarge) msymbol(circle_hollow)), ///
				ytitle("Proportion with SVR") ///
				ytitle(, size(large)) ///
				ylabel(0(0.2)1, labsize(large) angle(horizontal) format(%2.1f) nogrid) ///
				xtitle(Year initiating anti-HCV treatment) xtitle(, size(large)) ///
				xlabel(2006(1)2018, labsize(large) angle(forty_five)) ///
				legend(off) ///
				graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
				plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))






				
/****************** end of do-file ******************/
