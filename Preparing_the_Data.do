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
*log using "Preparing_the_Data_Log.txt", replace

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
* PREPARING THE DATA
************************************************************

* The file "dataset_master.dta" contains the following variables:
* - the polygenic scores 
*		(NOTE: the PGSraw_XX variables are the PLINK scores; the PGSldpPXX_XX are the LDpred scores)
* - the top principal components of the genetic relatedness matrix
* - the following variables from the The RAND HRS Data File (v.N) ("rndhrs_n.dta"):
* 		hhidpn s1hhidpn h*hhid raevbrn raevbrnf r*bmi r*height r*psyche r*diabe raedyrs inw* hacohort racohbyr raracem rahispan rabyear ragender r*iwstat
*		(NOTE: rabyear and ragender have been renamed "byear" and "sex")	
* - the HRS variable cholesterol_w1_v435 
* - the variable saliva_consent (which equals 1 if HRS variables LI913 or KI913 equal 1) 
* - the variable age_childrenMeas (which indicates the age at which number of children was measured) 
*		(constructed from RAND HRS variable rabyear and from HRS variables V145, D668, E668, F1006, G1093, HB033, JB033, KB033, LB033, MB033, NB033)
* - the variables z_cognition and z_neuroticism (constructed as described in section 6.4 of the Supporting Information of Okbay et al Nature 2016)
use dataset_master_initial.dta, clear

************************************************************
************************************************************
* Dropping individuals who will not be used for any analyses

drop if byear < 24 | byear > 53
drop if (age_childrenMeas < 50 ) & sex == 0 
drop if (age_childrenMeas < 45 ) & sex == 1 
gen year_childrenMeas = byear + age_childrenMeas
keep if year_childrenMeas <= 2008
keep if ( raracem == 1 & rahispan == 0 ) | (hwe_eur_keep == 1)

gen study_sample=0
replace study_sample=1 if byear >=31  & byear <= 53 

************************************************************
************************************************************
* Constructing the rLRS variable (which equals LRS divided by the mean LRS for the previous, current, and subsequent years)
* 	(NOTE: I do this separately for the HRS0 cohort and for the HRS1-3 cohorts togeter, since selection bias may affect the HRS0 cohort
*	 Also, for the "corner" years (i.e., 1924 and 1930 for the HRS0 cohort and 1931 and 1953 for the HRS1-3 cohorts), I use 2-year windows instead of 3-year windows) 

gen raevbrn_mean3yr = .
gen raevbrn_sd3yr = .
gen raevbrn_n3yr= .

sort sex byear

forval sex = 0/1 {
* The "inside" years for the HRS0 cohort [born 1924-29]:
forval byr = 25/29 {
	gen temp1 = raevbrn if byear>=`byr'-1 & byear<=`byr'+1 & sex==`sex'
	egen temp_mean = mean(temp1)
	egen temp_sd = sd(temp1)
	egen temp_n = count(temp1)
	qui replace raevbrn_mean3yr=temp_mean if byear==`byr' & sex==`sex'
	* We generate also the sd and n over 3 years -- needed for the plots:
	qui replace raevbrn_sd3yr=temp_sd if byear==`byr' & sex==`sex'
	qui replace raevbrn_n3yr=temp_n if byear==`byr' & sex==`sex'
	drop temp*
}
* The "inside" years for the HRS1-3 cohorts [born 1931-53]
forval byr = 32/52 {
	gen temp1 = raevbrn if byear>=`byr'-1 & byear<=`byr'+1 & sex==`sex'
	egen temp_mean = mean(temp1)
	egen temp_sd = sd(temp1)
	egen temp_n = count(temp1)
	qui replace raevbrn_mean3yr=temp_mean if byear==`byr' & sex==`sex'
	* We generate also the sd and n over 3 years -- needed for the plots:
	qui replace raevbrn_sd3yr=temp_sd if byear==`byr' & sex==`sex'
	qui replace raevbrn_n3yr=temp_n if byear==`byr' & sex==`sex'	
	drop temp*
}
* The "lower corner" years for the HRS0 and the HRS1-3 cohorts:
forval byr = 24(7)31 {
	gen temp1 = raevbrn if byear>=`byr' & byear<=`byr'+1 & sex==`sex'
	egen temp_mean = mean(temp1)
	egen temp_sd = sd(temp1)
	egen temp_n = count(temp1)	
	qui replace raevbrn_mean3yr=temp_mean if byear==`byr' & sex==`sex'
	* We generate also the sd and n over 3 years -- needed for the plots:
	qui replace raevbrn_sd3yr=temp_sd if byear==`byr' & sex==`sex'
	qui replace raevbrn_n3yr=temp_n if byear==`byr' & sex==`sex'
	drop temp*
}
* The "upper corner" years for the HRS0 and the HRS1-3 cohorts:
forval byr = 30(23)53 {
	gen temp1 = raevbrn if byear>=`byr'-1 & byear<=`byr' & sex==`sex'
	egen temp_mean = mean(temp1)
	egen temp_sd = sd(temp1)
	egen temp_n = count(temp1)
	qui replace raevbrn_mean3yr=temp_mean if byear==`byr' & sex==`sex'
	* We generate also the sd over and n 3 years -- needed for the plots:
	qui replace raevbrn_sd3yr=temp_sd if byear==`byr' & sex==`sex'
	qui replace raevbrn_n3yr=temp_n if byear==`byr' & sex==`sex'
	drop temp*
}
}

* Defining the rLRS variable:
gen relLRS = raevbrn/raevbrn_mean3yr
* Define also absolute LRS:
gen absLRS = raevbrn
* Also define childlessness (needed for summary stats table)
gen childless=0 if absLRS>0 & absLRS<100
replace childless=1 if absLRS==0
tab absLRS childless

************************************************************
************************************************************
* Defining the height variable

* Generate the mean residualized height separately for males and females:
cap drop resid*
forval w = 1/11 {
	qui reg r`w'height i.byear if sex == 0
	predict resid0w`w' if sex == 0, resid
	qui reg r`w'height i.byear if sex == 1
	predict resid1w`w' if sex == 1, resid 
}
egen height0 = rowmean(resid0w*)
egen height1 = rowmean(resid1w*)
gen height_resid = height0  if sex == 0
replace height_resid = height1 if sex == 1
drop resid* height0 height1

* Generating the mean height across all individuals and all waves
egen height0_mean = rowmean(r*height) if sex == 0
sum height0_mean
scalar height0_mean_sca = r(mean)
egen height1_mean = rowmean(r*height) if sex == 1
sum height1_mean
scalar height1_mean_sca = r(mean)

* Adding the resulting height means to the height residual, to obtain our height phenotype
gen height_clean = height0_mean_sca + height_resid if sex == 0
replace height_clean = height1_mean_sca + height_resid if sex == 1

* NOTE: height_resid and height_clean are perfectly correlated by construction, so can use either in our analyses. But need height_clean for sum stats table
bysort sex: corr height_resid height_clean

************************************************************
************************************************************
* Defining the BMI variable

* Generate the mean residualized BMI separately for males and females:
cap drop resid*
forval w = 1/11 {
	qui reg r`w'bmi i.byear if sex == 0
	predict resid0w`w' if sex == 0, resid
	qui reg r`w'bmi i.byear if sex == 1
	predict resid1w`w' if sex == 1, resid 
}
egen bmi0 = rowmean(resid0w*)
egen bmi1 = rowmean(resid1w*)
gen bmi_resid = bmi0  if sex == 0
replace bmi_resid = bmi1 if sex == 1
drop resid* bmi0 bmi1

* Generating the mean bmi across all individuals and all waves
egen bmi0_mean = rowmean(r*bmi) if sex == 0
sum bmi0_mean
scalar bmi0_mean_sca = r(mean)
egen bmi1_mean = rowmean(r*bmi) if sex == 1
sum bmi1_mean
scalar bmi1_mean_sca = r(mean)

* Adding the resulting bmi means to the bmi residual, to obtain our bmi phenotype
gen bmi_clean = bmi0_mean_sca + bmi_resid if sex == 0
replace bmi_clean = bmi1_mean_sca + bmi_resid if sex == 1

* NOTE: bmi_resid and bmi_clean are perfectly correlated by construction, so can use either in our analyses. But need bmi_clean for sum stats table
bysort sex: corr bmi_resid bmi_clean

************************************************************
************************************************************

* Generate an indicator equal to the last r*diabe variable (each of which is a 0/1 variable for have you ever had diabetes or high blood sugar):
gen diabe = r11diabe
forval w = 1/10 {
	local wloc = 11-`w'
	replace diabe = r`wloc'diabe if missing(diabe)
}

************************************************************
************************************************************

sort IID 
compress
order IID FID hhidpn s1hhidpn LOCAL_ID sex byear byear2
save dataset_master.dta, replace

sum

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
* CHOOSING THE GAUSSIAN MIXTURE WEIGHT FOR THE LDPRED SCORES FOR EACH PHENOTYPE
* As mentioned in the Supporting Information of my article: 
* 	"For each phenotype, LDpred scores were constructed for each of the following Gaussian mixture weights: 
*	0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, and 1. 
* 	For BMI, EA, HGT, and TC—-for which there are phenotypic variables in the HRS—-I selected the weights that maximize 
*	the incremental R2 of each score in an OLS regression of the phenotypic variable on the score and on variables 
*	for sex, birth year, birth year squared, and the top 20 principal components of the genetic relatedness matrix. 
*	For each of GLU, SCZ, and AAM—-for which there are no phenotypic variables in the HRS—-I selected the weights that
*	maximize the correlations between the score and known correlates of the phenotype, controlling for sex, birth year, 
*	birth year squared, and the top 20 principal compo- nents. 
*	For GLU, I selected the weight that maximizes the correlation with a variable indicating if an individual ever had diabetes or high blood sugar; 
*	for SCZ, I selected the weight that maximizes the correlations with neuroticism (ref. 60) and cognitive ability (ref. 61); 
*	for AAM, I selected the weight that maximizes the correlations with HGT and BMI (ref. 62)." 
************************************************************

use dataset_master.dta, clear

**************************************************
**************************************************
***

* PGS FOR BMI (ON BMI)

scalar drop _all

qui reg bmi_resid 			sex byear byear2 pc*   if study_sample==1
scalar R2o=e(r2)
qui reg  bmi_resid PGSraw_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_raw=e(r2)-R2o
qui reg bmi_resid PGSldpP1Em4_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em4=e(r2)-R2o
qui reg bmi_resid PGSldpP3Em4_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em4=e(r2)-R2o
qui reg bmi_resid PGSldpP1Em3_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em3=e(r2)-R2o
qui reg bmi_resid PGSldpP3Em3_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em3=e(r2)-R2o
qui reg bmi_resid PGSldpP1Em2_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em2=e(r2)-R2o
qui reg bmi_resid PGSldpP3Em2_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em2=e(r2)-R2o
qui reg bmi_resid PGSldpP1Em1_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em1=e(r2)-R2o
qui reg bmi_resid PGSldpP3Em1_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em1=e(r2)-R2o
qui reg bmi_resid PGSldpP1_BMI sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1=e(r2)-R2o

scalar dir

* The LDpred score with the Gaussian mixture weight P=1Em1 (i.e., 0.1) has the highest incremental R2, with dR2_ldpP1Em1 =  .08904619

************************************************************
************************************************************

* PGS FOR EA

qui reg raedyrs 			sex byear byear2 pc*   if study_sample==1
scalar R2o=e(r2)
qui reg raedyrs PGSraw_EA3 sex byear byear2 pc*   if study_sample==1   
scalar dR2_raw=e(r2)-R2o
qui reg raedyrs PGSldpP1Em4_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em4=e(r2)-R2o
qui reg raedyrs PGSldpP3Em4_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em4=e(r2)-R2o
qui reg raedyrs PGSldpP1Em3_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em3=e(r2)-R2o
qui reg raedyrs PGSldpP3Em3_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em3=e(r2)-R2o
qui reg raedyrs PGSldpP1Em2_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em2=e(r2)-R2o
qui reg raedyrs PGSldpP3Em2_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em2=e(r2)-R2o
reg raedyrs PGSldpP1Em1_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em1=e(r2)-R2o
qui reg raedyrs PGSldpP3Em1_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em1=e(r2)-R2o
qui reg raedyrs PGSldpP1_EA3 sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1=e(r2)-R2o

scalar dir

* The LDpred score with the Gaussian mixture weight P=1Em1 (i.e., 0.1) has the highest incremental R2, with dR2_ldpP1Em1 =  .07427988

**************************************************
**************************************************

* PGS FOR GLUCOSE


scalar drop _all

qui reg diabe 			sex byear byear2 pc*   if study_sample==1 
scalar R2o=e(r2)
reg  diabe PGSraw_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_raw=e(r2)-R2o
reg diabe PGSldpP1Em4_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_ldpP1Em4=e(r2)-R2o
reg diabe PGSldpP3Em4_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_ldpP3Em4=e(r2)-R2o
reg diabe PGSldpP1Em3_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_ldpP1Em3=e(r2)-R2o
reg diabe PGSldpP3Em3_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_ldpP3Em3=e(r2)-R2o
reg diabe PGSldpP1Em2_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_ldpP1Em2=e(r2)-R2o
reg diabe PGSldpP3Em2_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_ldpP3Em2=e(r2)-R2o
reg diabe PGSldpP1Em1_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_ldpP1Em1=e(r2)-R2o
reg diabe PGSldpP3Em1_GLU sex byear byear2 pc*   if study_sample==1 
scalar dR2_ldpP3Em1=e(r2)-R2o
reg diabe PGSldpP1_GLU sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1=e(r2)-R2o

scalar dir

* The LDpred score with the Gaussian mixture weight P=3Em2 (i.e., 0.03) has the highest incremental R2, with dR2_ldpP1Em1 =  .0077873

**************************************************
**************************************************

* PGS FOR HEIGHT 

scalar drop _all

qui reg height_resid 			sex byear byear2 pc*   if study_sample==1
scalar R2o=e(r2)
qui reg  height_resid PGSraw_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_raw=e(r2)-R2o
qui reg height_resid PGSldpP1Em4_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em4=e(r2)-R2o
qui reg height_resid PGSldpP3Em4_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em4=e(r2)-R2o
qui reg height_resid PGSldpP1Em3_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em3=e(r2)-R2o
qui reg height_resid PGSldpP3Em3_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em3=e(r2)-R2o
qui reg height_resid PGSldpP1Em2_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em2=e(r2)-R2o
qui reg height_resid PGSldpP3Em2_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em2=e(r2)-R2o
qui reg height_resid PGSldpP1Em1_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em1=e(r2)-R2o
qui reg height_resid PGSldpP3Em1_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em1=e(r2)-R2o
qui reg height_resid PGSldpP1_HGT sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1=e(r2)-R2o

scalar dir

* The LDpred score with the Gaussian mixture weight P=1 has the highest incremental R2, with dR2_ldpP1Em1 =  .17370011

**************************************************
**************************************************

* PGS FOR SCZ 

foreach P in 1Em4 3Em4 1Em3 3Em3 1Em2 3Em2 1Em1 3Em1 1 {
	display "----- P = `P' -----"
	reg PGSldpP`P'_SCZ2 z_cognition z_neuroticism sex	byear byear2 pc*    if study_sample==1
}
display "----- Raw score -----"
	reg PGSraw_SCZ2  	z_cognition z_neuroticism sex	byear byear2 pc*	  if study_sample==1 

* The LDpred score the Gaussian mixture weight P=3m1 (i.e., 0.3) has the highest correlations with z_cognition and z_neuroticism (with the correct signs)

**************************************************
**************************************************

* PGS FOR TC (TOTAL CHOLESTEROL)

scalar drop _all

qui reg cholesterol_w1_v435 			sex byear byear2 pc*   if study_sample==1
scalar R2o=e(r2)
qui reg  cholesterol_w1_v435 PGSraw_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_raw=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP1Em4_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em4=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP3Em4_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em4=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP1Em3_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em3=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP3Em3_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em3=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP1Em2_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em2=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP3Em2_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em2=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP1Em1_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1Em1=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP3Em1_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP3Em1=e(r2)-R2o
qui reg cholesterol_w1_v435 PGSldpP1_TC sex byear byear2 pc*   if study_sample==1
scalar dR2_ldpP1=e(r2)-R2o

scalar dir

* The LDpred score with the Gaussian mixture weight P=3Em1 (i.e., 0.3) has the highest incremental R2, with dR2_ldpP3Em1 =  .01226485

**************************************************
**************************************************

* PGS FOR MENARCHE (for females only)

foreach P in 1Em4 3Em4 1Em3 3Em3 1Em2 3Em2 1Em1 3Em1 1 {
	display "----- P = `P' -----"
	reg PGSldpP`P'_AAM height_resid bmi_resid byear byear2 pc*  if sex == 1 & study_sample==1
}	
display "----- Raw score -----"
	reg PGSraw_AAM height_resid bmi_resid byear byear2 pc*  if sex == 1 & study_sample==1

* The LDpred score with the Gaussian mixture weight P=3m1 (i.e., 0.3) has the highest correlations with height and bmi (with the correct signs)

**************************************************
**************************************************

* Defining the PGS variables
gen PGSldp_EA3 = PGSldpP1Em1_EA3
gen PGSldp_HGT = PGSldpP1_HGT
gen PGSldp_BMI = PGSldpP1Em1_BMI
gen PGSldp_SCZ2 = PGSldpP3Em1_SCZ2
gen PGSldp_AAM = PGSldpP3Em1_AAM
gen PGSldp_GLU = PGSldpP3Em2_GLU
gen PGSldp_TC = PGSldpP3Em1_TC


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

* Defining the phenotype variables for the main analyses
* 	(Note: we set pheno_ to 1 for SCZ2, menarche, Glucose, for then Stata will simply omit these in the regressions in the for loop)
gen pheno_EA3 = raedyrs
* Express height in cm rather than meters:
gen pheno_HGT = height_clean*100
gen pheno_BMI = bmi_clean
gen pheno_SCZ2 = 1
gen pheno_AAM = 1
gen pheno_GLU = 1
gen pheno_TC = cholesterol_w1_v435

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
* SAVE THE MASTER DATASET

sort IID 
compress
save dataset_master.dta, replace
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************







 
* Close log file
cap log close
