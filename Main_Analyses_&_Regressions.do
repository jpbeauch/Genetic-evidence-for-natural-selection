************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
* PRELIMINARIES

clear all
set maxvar 12000

cap log close
set more off
cd "PATH"

* Open log file
set logtype text
*log using "Main analyses -- study sample", replace
*log using "Main regressions -- HRS0 cohort only", replace
*log using "Main regressions -- HRS1 cohort only", replace
*log using "Main regressions -- HRS2 cohort only", replace
*log using "Main regressions -- HRS3 cohort only", replace
*log using "Main regressions -- study sample, with LRS (instead of rLRS)", replace
*log using "Main regressions -- study sample, age 50-70 for females and 55-70 for males", replace
*log using "Main regressions -- all cohorts (HRS0-3) together", replace

************************************************************
************************************************************
************************************************************
************************************************************
************************************************************









************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
* LOAD THE DATA AND DEFINE THE SAMPLE AND THE OPTIONS FOR THE ANALYSES
************************************************************

use dataset_master.dta, clear

************************************************************
************************************************************

* Define relative LRS as our ReproSuccess variable (for the baseline analyses)
gen ReproSuccess = relLRS

************************************************************
************************************************************

* To keep only the genotyped individuals 
* (For the analyses with the scores, only genotyped individuals are used by default;
*  but for the analyses with the phenotypic variables, we can either use the sample of all individuals (the default) 
*  or only use the sample of genotyped individuals. For the latter, unstar the following line:
*keep if hwe_eur_keep==1


************************************************************
************************************************************
* For the robustness checks:

* 1. With LRS instead of rLRS as the ReproSuccess variable
*replace ReproSuccess = absLRS

* 2. With the PLINK2 scores (already outputted below -- no need to alter anything in the code)

* 3. With individuals who were 
*	(i) 70 or less in 2008 (born >= 1938) at the time of the last genotyping wave (to see if selection on survival after 70 years old biases results)
* 	(ii) age_childrenMeas>=50 for females and >=55 for males
/***
keep if byear >=38  & byear <= 53 
drop if (age_childrenMeas < 55) & sex == 0 
drop if (age_childrenMeas < 50) & sex == 1
***/

* 4. With the study sample (the HRS1-3 cohorts) together with the HRS0 cohort (to run this, comment out the line "keep if study_sample==1" below)

************************************************************
************************************************************
* Defining the sample for the analyses

* To use the study sample (this sample is used for the main analyses and combines the HRS1-3 cohorts):
keep if study_sample==1

* To use only the HRS0 cohort:
*keep if byear >=24  & byear <= 30

* To use only the HRS1 cohort:
*keep if byear >=31  & byear <= 41

* To use only the HRS2 cohort:
*keep if byear >=42  & byear <= 47 

* To use only the HRS3 cohort:
*keep if byear >=48  & byear <= 53 
* NOTE: there are no observations for pheno_TC for the males in the HRS3 cohort => set pheno_TC to 1 when we need to run the other analyses for males in the HRS3 cohort
*replace pheno_TC = 1

************************************************************
************************************************************
************************************************************
************************************************************
************************************************************







************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
* RUNNING THE MAIN REGRESSION ANALYSES (FOR TABLES 1 AND 2)
************************************************************

cap drop PGS
gen PGS = .
cap drop Pheno
gen Pheno = .
cap drop PhenoXgenotyped 
gen PhenoXgenotyped = .

************************************************************
************************************************************
* REGRESSION ON THE PGS -- females only 

* Regress ReproSuccess on the LDpred PGS ("PGS"), on the PLINK PGS ("PGSraw"), and on the phenotypic variables:
foreach var in  BMI EA3 GLU HGT SCZ2 TC AAM  {
	qui replace PGS = PGSldp_`var'
	qui reg ReproSuccess PGS i.byear pc*     i.hacohort  if sex == 1
	est store `var'_ldp
	qui replace PGS = PGSraw_`var'
	qui reg ReproSuccess PGS 	i.byear pc*	i.hacohort  if sex == 1
	est store `var'_raw
	qui replace Pheno = pheno_`var'
	qui reg ReproSuccess Pheno i.byear			i.hacohort    if sex == 1
	est store `var'_pheno
}
* Estimate table for each trait separately:
estimates table BMI_ldp EA3_ldp GLU_ldp HGT_ldp SCZ2_ldp TC_ldp  AAM_ldp, keep(PGS) b se t p stats(N)	
estimates table BMI_raw EA3_raw GLU_raw HGT_raw SCZ2_raw TC_raw  AAM_raw, keep(PGS) b se t p stats(N)
estimates table BMI_pheno EA3_pheno GLU_pheno  HGT_pheno SCZ2_pheno TC_pheno  AAM_pheno, keep(Pheno) b se t p stats(N)

************************************************************
************************************************************
* REGRESSION ON THE PGS -- males only (=> sex == 0) 

* Regress ReproSuccess on the LDpred PGS ("PGS"), on the PLINK PGS ("PGSraw"), and on the phenotypic variables:
foreach var in    BMI EA3 GLU HGT SCZ2 TC {
	qui replace PGS = PGSldp_`var'
	qui reg ReproSuccess PGS i.byear pc*  i.hacohort    if sex == 0
	est store `var'_ldp
	qui replace PGS = PGSraw_`var'
	qui reg ReproSuccess PGS 	i.byear pc*  i.hacohort	if sex == 0
	est store `var'_raw
	qui replace Pheno = pheno_`var'
	qui reg ReproSuccess Pheno i.byear 			i.hacohort    if sex == 0
	est store `var'_pheno
}
* Estimate table for each trait separately:
estimates table BMI_ldp EA3_ldp GLU_ldp HGT_ldp SCZ2_ldp TC_ldp, keep(PGS) b se t p stats(N)	
estimates table BMI_raw EA3_raw GLU_raw HGT_raw SCZ2_raw TC_raw, keep(PGS) b se t p stats(N)
estimates table BMI_pheno EA3_pheno GLU_pheno HGT_pheno SCZ2_pheno TC_pheno, keep(Pheno) b se t p stats(N)

************************************************************
************************************************************
* REGRESSION ON THE PGS -- males & females pooled 

* NOTE: As described in the Supporting Information of my article, to ensure that my regressions are well-specified 
*	(given that females and males in the same household will typically have the same number of children),
*	I keep only one person per household for these regressions.

* Create a variable ("HH_PN") that identifies the respondent with the lower PN in each household
sort genotyped HHID PN
bysort genotyped HHID: gen HH_PN = 1 if HHID[_n]~=HHID[_n-1]
replace HH_PN = 0 if missing(HH_PN)

* Regress ReproSuccess on the LDpred PGS ("PGS"):
foreach var in  BMI EA3 GLU HGT SCZ2 TC AAM  {
	qui replace PGS = PGSldp_`var'
	qui reg ReproSuccess PGS i.byear pc*  sex    i.hacohort  if  HH_PN == 1
	est store `var'_ldp
	qui replace PGS = PGSraw_`var'
}
* Estimate table for each trait separately:
estimates table BMI_ldp EA3_ldp GLU_ldp HGT_ldp SCZ2_ldp TC_ldp , keep(PGS) b se t p stats(N)	


************************************************************
************************************************************
************************************************************
************************************************************
************************************************************













***
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
* CALCULATING (1) THE RESCALED DIRECTIONAL SELECTION DIFFERENTIAL FOR THE PGS OF EA (IN YEARS OF EA)
*		  AND (2) THE DIRECTIONAL SELECTION DIFFERENTIAL FOR EA (IN YEARS OF EA)
* 		AND THEIR SE'S/95% CI'S

********************************
 
 
* 0. Defining a program for the R2 of the PGS of EA
capture program drop PROG_R2_PGSofEA
program PROG_R2_PGSofEA, rclass
	qui reg pheno_EA3	sex i.byear pc*
	scalar R2o=e(r2)
	qui reg  pheno_EA3 PGSldp_EA3 sex i.byear pc*
	return scalar dR2=e(r2)-R2o
	* To test:
	scalar dR2=e(r2)-R2o
	*scalar dir dR2
end

* Call the program (to test it)
*PROG_R2_PGSofEA

*****
* 1. Defining a program for (1) the rescaled directional selection differential for the PGS of EA, in years of EA
*	 	The program also calls PROG_R2_PGSofEA and returns its output, so we can use the bootstrap command only once for both quantities

capture program drop PROG_dz_PGS_yrsEA
program PROG_dz_PGS_yrsEA, rclass
	* i. First, estimate dz_PGS in Haldanes (this is the coefficient from the main regressions)
	qui reg ReproSuccess PGSldp_EA3 i.byear pc*     i.hacohort
	mat b=e(b)
	scalar dz_PGS = b[1,1]
	* ii. Second, estimate sigma_PGS_yrsEA (=SD(EA)*R2(PGSofEA))
	PROG_R2_PGSofEA
	scalar dR2 = r(dR2)
	qui sum pheno_EA3
	scalar sigma_PGS_yrsEA = r(sd)*sqrt(dR2)
	* Calulate and return dz_PGS_yrsEA
	return scalar dz_PGS_yrsEA = dz_PGS * sigma_PGS_yrsEA
	return scalar dR2 = dR2
	* To test:
	scalar dz_PGS_yrsEA = dz_PGS * sigma_PGS_yrsEA
	*scalar dir
end

* Call the program (to test it)
*PROG_dz_PGS_yrsEA

*****
* 2. Defining a program for (2) the directional selection differential for EA (in years of EA)
* 		The program also calls PROG_dz_PGS_yrsEA and returns its output, so we can use the bootstrap command only once for both quantities

* Assumed heritability of EA:
scalar h2_EA = 0.40

capture program drop PROG_dz_EA_yrsEA
program PROG_dz_EA_yrsEA, rclass
	PROG_dz_PGS_yrsEA
	scalar dz_PGS_yrsEA=r(dz_PGS_yrsEA)
	scalar dR2 = r(dR2)
	*** 
	* Calulate and return dz_EA_yrsEA
	return scalar dz_EA_yrsEA = dz_PGS_yrsEA * h2_EA / dR2
	return scalar dz_PGS_yrsEA = dz_PGS_yrsEA
	* To test:
	scalar dz_EA_yrsEA = dz_PGS_yrsEA * h2_EA / dR2
	scalar dir
end

* Call the program (to test it)
*PROG_dz_EA_yrsEA


*****

* Bootstrap for females:
preserve
keep if  ~missing(pheno_EA3, ReproSuccess, pc1, byear, hacohort) & sex == 1  
tab sex
bootstrap dz_PGS_yrsEA_females=r(dz_PGS_yrsEA) dz_EA_yrsEA_females=r(dz_EA_yrsEA) ///
			if ~missing(pheno_EA3, ReproSuccess, pc1, byear, hacohort) & sex == 1, reps(1000): PROG_dz_EA_yrsEA
estat bootstrap, all
restore

* Bootstrap for males:
preserve
keep if  ~missing(pheno_EA3, ReproSuccess, pc1, byear, hacohort) & sex == 0 
tab sex
bootstrap dz_PGS_yrsEA_males=r(dz_PGS_yrsEA) dz_EA_yrsEA_males=r(dz_EA_yrsEA) ///
			if ~missing(pheno_EA3, ReproSuccess, pc1, byear, hacohort) & sex == 0, reps(1000): PROG_dz_EA_yrsEA
estat bootstrap, all			
restore
			

************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***/








************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
* 2. TESTING FOR NONLINEAR SELECTION (LANDE & ARNOLD REGRESSION) 

* Generate variables for the score interactions
scalar count1 = 0
foreach pheno1 in BMI EA3 GLU HGT SCZ2 TC AAM {
	gen PGS_`pheno1' = PGSldp_`pheno1'
	scalar count1 = count1 + 1
	scalar count2 = 0
	foreach pheno2 in BMI EA3 GLU HGT SCZ2 TC AAM  {
		scalar count2 = count2 + 1
		if count2>=count1 {
			gen PGSx_`pheno1'_`pheno2' = PGSldp_`pheno1'*PGSldp_`pheno2'
		}
	}
}

* For females (with AAM)
reg ReproSuccess   PGS_BMI PGS_EA3 PGS_GLU PGS_HGT PGS_SCZ2 PGS_TC PGS_AAM  ///
			PGSx_BMI_BMI PGSx_EA3_EA3 PGSx_GLU_GLU PGSx_HGT_HGT PGSx_SCZ2_SCZ2 PGSx_TC_TC PGSx_AAM_AAM ///
			PGSx_BMI_EA3 PGSx_BMI_GLU PGSx_BMI_HGT PGSx_BMI_SCZ2 PGSx_BMI_TC PGSx_BMI_AAM ///
			PGSx_EA3_GLU PGSx_EA3_HGT PGSx_EA3_SCZ2 PGSx_EA3_TC PGSx_EA3_AAM ///
			PGSx_GLU_HGT PGSx_GLU_SCZ2 PGSx_GLU_TC PGSx_GLU_AAM ///
			PGSx_HGT_SCZ2 PGSx_HGT_TC PGSx_HGT_AAM ///
			PGSx_SCZ2_TC PGSx_SCZ2_AAM ///
			PGSx_TC_AAM ///
			i.byear pc*	i.hacohort if sex == 1
			
* For males (without AAM)
foreach pheno in PGS_AAM PGSx_AAM_AAM PGSx_BMI_AAM PGSx_EA3_AAM PGSx_GLU_AAM PGSx_HGT_AAM PGSx_SCZ2_AAM PGSx_TC_AAM  {
	replace `pheno' = 0
}
reg ReproSuccess   PGS_BMI PGS_EA3 PGS_GLU PGS_HGT PGS_SCZ2 PGS_TC PGS_AAM  ///
			PGSx_BMI_BMI PGSx_EA3_EA3 PGSx_GLU_GLU PGSx_HGT_HGT PGSx_SCZ2_SCZ2 PGSx_TC_TC PGSx_AAM_AAM ///
			PGSx_BMI_EA3 PGSx_BMI_GLU PGSx_BMI_HGT PGSx_BMI_SCZ2 PGSx_BMI_TC PGSx_BMI_AAM ///
			PGSx_EA3_GLU PGSx_EA3_HGT PGSx_EA3_SCZ2 PGSx_EA3_TC PGSx_EA3_AAM ///
			PGSx_GLU_HGT PGSx_GLU_SCZ2 PGSx_GLU_TC PGSx_GLU_AAM ///
			PGSx_HGT_SCZ2 PGSx_HGT_TC PGSx_HGT_AAM ///
			PGSx_SCZ2_TC PGSx_SCZ2_AAM ///
			PGSx_TC_AAM ///
			i.byear pc*	i.hacohort if sex == 0


************************************************************
************************************************************
************************************************************
************************************************************
************************************************************






 
* Close logs
cap log close
