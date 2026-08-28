/**************************************************************************
 REPLICATION RUN INSTRUCTIONS

 Set Stata's working directory to the Ricci_replication folder before
 running this do-file. The code will then use that folder as ROOT and keep
 the original relative folder structure (raw/, output_*/, calibration_*).

Required user-written Stata command:
  - boottest

 No analytical commands, parameter values, sample definitions, estimation
 choices, or output construction have been changed from the authors'
 original file. Only the machine-specific ROOT path has been made portable.
**************************************************************************/

/****************************************************************************************
 EXTEND COAUTHOR'S RICCI/PWT CALIBRATION TO OBSERVED WIOD YEARS, 1995-2009

 This do-file uses the two Table-2 datasets produced by:
   ricci_all_years_with_ROW_estimate_v2.do

 It does the coauthor's calibration exercise twice:
   (A) Ricci-comparable 12-region Table 2 (no ROW)
   (B) 12 Ricci regions + estimated ROW

 It:
   1. Aggregates PWT 11.0 to Ricci's regions for every year 1995-2009.
   2. Reproduces the attached calib_inputs.csv construction.
   3. Exports two 1995-2009 calib-input CSVs with the SAME 14 columns as the
      coauthor's file.
   4. Estimates theta in Stata using the transformed form of equation (14).
   5. Uses aggregate capital Omega = K_core/K_region as the PRIMARY power measure; per-worker power is retained as a robustness check.
   6. Repeats the lambda/b/delta grid robustness exercise using aggregate Omega.
   7. Recreates the empirical extraction-schedule figure in Stata using aggregate Omega and ONLY
      regions with net transfer < 0 in EVERY year of 1995-2009.
   8. Omits the Hickel et al. band from the figure.

 IMPORTANT SIGN CONVENTION
   Ricci Table 2 reports providers/outflows as NEGATIVE net transfers.
   The model's extraction rate is therefore:
       e = - net_transfer_pct_va / 100
   Thus a region with an outflow has e > 0.

 BASELINE MODEL PARAMETERS, matching the coauthor exercise:
   lambda = 0.5, b = 0.5, delta = 0.5

 Equation (14):
   e = m / (chi0*Omega^(-theta) - lambda*delta)
   m = 1 + lambda*(1-b)

 Linearized estimating equation:
   ln(m/e + lambda*delta) = ln(chi0) - theta*ln(Omega)

 The attached 1995/2000/2007 calib_inputs.csv gives theta ~= 0.353 with
 SE ~= 0.116 under this transformed OLS specification.
****************************************************************************************/

version 17.0
clear all
set more off
set varabbrev off

********************************************************************************
* 0. PATHS
********************************************************************************

global ROOT "`c(pwd)'"
global RAW  "$ROOT/raw"
global OUT  "$ROOT/output_all_years"
global CAL  "$ROOT/calibration_1995_2009_aggregate"

capture mkdir "$CAL"

global PWT "$RAW/pwt110.dta"
global T2REG "$OUT/table2_12regions_all_years_long.dta"
global T2ROW "$OUT/table2_13regions_ROW_est_all_years_long.dta"

foreach f in "$PWT" "$T2REG" "$T2ROW" {
    capture confirm file `f'
    if _rc {
        di as error "Required file not found: `f'"
        di as error "Place pwt110.dta in $RAW and run ricci_all_years_with_ROW_estimate_v2.do first."
        exit 601
    }
}

local FIRSTYEAR 1995
local LASTYEAR  2009
local NYEARS = `LASTYEAR' - `FIRSTYEAR' + 1

* Coauthor baseline parameters.
local lambda 0.5
local b      0.5
local delta  0.5
local m = 1 + `lambda'*(1-`b')

* The original Ricci 40-country sample.
local CTRY40 "AUS AUT BEL BGR BRA CAN CHN CYP CZE DEU DNK ESP EST FIN FRA GBR GRC HUN IDN IND IRL ITA JPN KOR LTU LUX LVA MEX MLT NLD POL PRT ROU RUS SVK SVN SWE TUR TWN USA"

********************************************************************************
* 1. BUILD PWT 11.0 COUNTRY-YEAR DATA AND RICCI REGION CROSSWALK
********************************************************************************

use "$PWT", clear
keep if inrange(year,`FIRSTYEAR',`LASTYEAR')
keep countrycode country year cn emp cgdpo labsh csh_i
rename countrycode iso3

* The 40 Ricci/WIOD countries must be complete for these calibration inputs.
gen byte in_wiod40 = strpos(" `CTRY40' ", " " + iso3 + " ") > 0
count if in_wiod40 & missing(cn,emp,cgdpo,labsh,csh_i)
assert r(N)==0

* Ricci's 12-region classification.
gen str24 region = ""
replace region = "North America"   if inlist(iso3,"CAN","USA")
replace region = "North EMU"       if inlist(iso3,"AUT","BEL","DEU","FIN","FRA","LUX","NLD")
replace region = "South EMU"       if inlist(iso3,"CYP","ESP","GRC","IRL","ITA","MLT","PRT")
replace region = "North Europe"    if inlist(iso3,"DNK","GBR","SWE")
replace region = "East Europe"     if inlist(iso3,"BGR","CZE","EST","HUN","LTU")
replace region = "East Europe"     if inlist(iso3,"LVA","POL","ROU","SVK","SVN")
replace region = "Latin America"   if inlist(iso3,"BRA","MEX")
replace region = "China"           if iso3=="CHN"
replace region = "India"           if iso3=="IND"
replace region = "North East Asia" if inlist(iso3,"JPN","KOR")
replace region = "Other Asia"      if inlist(iso3,"IDN","TUR","TWN")
replace region = "Russia"          if iso3=="RUS"
replace region = "Australia"       if iso3=="AUS"

* Variables used by the coauthor's aggregation.
gen double labor_income_proxy = labsh*cgdpo
gen double investment_proxy   = csh_i*cgdpo

* Fixed core bloc from the paper: North America, North EMU, North Europe,
* and North East Asia. Keep this definition fixed through time so the extension
* is directly comparable to the original calibration.
gen byte core_country = inlist(region,"North America","North EMU","North Europe","North East Asia")

* Save country-level PWT data for later ROW construction.
tempfile pwt_country
save `pwt_country'

********************************************************************************
* 2. CORE BLOC AGGREGATES, 1995-2009
********************************************************************************

preserve
    keep if core_country
    collapse (sum) K_core=cn L_core=emp Y_core=cgdpo ///
                   labor_core=labor_income_proxy inv_core=investment_proxy, by(year)

    gen double alpha_core = labor_core/Y_core
    gen double s_core     = inv_core/Y_core

    keep year K_core L_core Y_core alpha_core s_core
    isid year
    count
    assert r(N)==`NYEARS'
    tempfile core_agg
    save `core_agg'
restore

********************************************************************************
* 3. PWT AGGREGATES FOR ALL 12 RICCI REGIONS
********************************************************************************

preserve
    keep if in_wiod40
    assert region!=""

    collapse (sum) K_per=cn L_per=emp Y_per=cgdpo ///
                   labor_per=labor_income_proxy inv_per=investment_proxy, by(year region)

    gen double alpha_per = labor_per/Y_per
    gen double s_per     = inv_per/Y_per

    merge m:1 year using `core_agg'
    assert _merge==3
    drop _merge

    gen double Omega_agg   = K_core/K_per
    gen double Omega_pw    = (K_core/L_core)/(K_per/L_per)
    gen double Y_ratio_agg = Y_core/Y_per
    gen double Y_ratio_pw  = (Y_core/L_core)/(Y_per/L_per)
    gen double sigma       = s_core/s_per

    keep year region Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
         alpha_core alpha_per sigma K_core K_per L_core L_per

    sort region year
    isid year region
    count
    assert r(N)==12*`NYEARS'

    tempfile pwt12
    save `pwt12'
restore

********************************************************************************
* 4. PWT PROXY FOR REST OF WORLD
********************************************************************************

* WIOD's ROW does not have a direct PWT aggregate. For the ROW calibration row,
* approximate its capital/labor side with all PWT economies outside Ricci's 40
* that have the five required PWT variables in a given year.
use `pwt_country', clear
keep if !in_wiod40
keep if !missing(cn,emp,cgdpo,labsh,csh_i)

collapse (sum) K_per=cn L_per=emp Y_per=cgdpo ///
               labor_per=labor_income_proxy inv_per=investment_proxy ///
         (count) row_pwt_n_countries=cn, by(year)

gen str24 region = "Rest of World (est.)"
gen double alpha_per = labor_per/Y_per
gen double s_per     = inv_per/Y_per

merge 1:1 year using `core_agg'
assert _merge==3
drop _merge

gen double Omega_agg   = K_core/K_per
gen double Omega_pw    = (K_core/L_core)/(K_per/L_per)
gen double Y_ratio_agg = Y_core/Y_per
gen double Y_ratio_pw  = (Y_core/L_core)/(Y_per/L_per)
gen double sigma       = s_core/s_per

keep year region Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
     alpha_core alpha_per sigma K_core K_per L_core L_per row_pwt_n_countries

sort year
isid year
count
assert r(N)==`NYEARS'

tempfile pwtrow
save `pwtrow'

* Append the ROW proxy to the 12 Ricci-region PWT file.
use `pwt12', clear
append using `pwtrow'
sort region year
tempfile pwt13
save `pwt13'

********************************************************************************
* 5. CALIBRATION INPUTS: REGULAR TABLE 2 (12 REGIONS)
********************************************************************************

use "$T2REG", clear
isid year region
keep if inrange(year,`FIRSTYEAR',`LASTYEAR')

merge 1:1 year region using `pwt12'
assert _merge==3
drop _merge

* Negative Ricci net transfer = provider/outflow = positive model extraction rate.
gen double e = -net_transfer_pct_va/100

* Flags are retained in the DTA diagnostic file but omitted from the exact-format CSV.
gen byte outflow = net_transfer<0
gen byte core_component = inlist(region,"North America","North EMU","North Europe","North East Asia")
bysort region: egen byte all_year_outflow = min(outflow)
bysort region: egen int  n_region_years = count(year)
replace all_year_outflow = 0 if n_region_years!=`NYEARS'

order year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
      alpha_core alpha_per sigma K_core K_per L_core L_per ///
      net_transfer net_transfer_pct_va outflow all_year_outflow ///
      core_component
sort region year

save "$CAL/calib_inputs_1995_2009_table2_reg.dta", replace

preserve
    keep year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
         alpha_core alpha_per sigma K_core K_per L_core L_per
    order year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
          alpha_core alpha_per sigma K_core K_per L_core L_per
    sort year region
    export delimited using "$CAL/calib_inputs_1995_2009_table2_reg.csv", replace
restore

preserve
    keep year region net_transfer net_transfer_pct_va e outflow all_year_outflow ///
         core_component
    sort region year
    export delimited using "$CAL/calib_inputs_1995_2009_table2_reg_diagnostics.csv", replace
restore

********************************************************************************
* 6. CALIBRATION INPUTS: TABLE 2 + ESTIMATED ROW (13 REGIONS)
********************************************************************************

use "$T2ROW", clear
isid year region
keep if inrange(year,`FIRSTYEAR',`LASTYEAR')

merge 1:1 year region using `pwt13'
assert _merge==3
drop _merge

gen double e = -net_transfer_pct_va/100

gen byte outflow = net_transfer<0
gen byte core_component = inlist(region,"North America","North EMU","North Europe","North East Asia")
gen byte row_estimated_calibration = region=="Rest of World (est.)"
bysort region: egen byte all_year_outflow = min(outflow)
bysort region: egen int n_region_years = count(year)
replace all_year_outflow = 0 if n_region_years!=`NYEARS'

order year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
      alpha_core alpha_per sigma K_core K_per L_core L_per ///
      net_transfer net_transfer_pct_va outflow all_year_outflow ///
      core_component row_estimated_calibration
sort region year

save "$CAL/calib_inputs_1995_2009_table2_ROW.dta", replace

preserve
    keep year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
         alpha_core alpha_per sigma K_core K_per L_core L_per
    order year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
          alpha_core alpha_per sigma K_core K_per L_core L_per
    sort year region
    export delimited using "$CAL/calib_inputs_1995_2009_table2_ROW.csv", replace
restore

preserve
    keep year region net_transfer net_transfer_pct_va e outflow all_year_outflow ///
         core_component row_estimated_calibration
    sort region year
    export delimited using "$CAL/calib_inputs_1995_2009_table2_ROW_diagnostics.csv", replace
restore

********************************************************************************
* 7. ESTIMATE THETA: REGULAR TABLE 2
*    PRIMARY SPECIFICATION = AGGREGATE CAPITAL POWER
********************************************************************************

use "$CAL/calib_inputs_1995_2009_table2_reg.dta", clear
keep if all_year_outflow==1
assert e>0
assert Omega_agg>0
assert Omega_pw>0

* Transformed equation (14):
* ln(m/e + lambda*delta) = ln(chi0) - theta*ln(Omega)
egen region_id=group(region)

gen double lhs = ln(`m'/e + `lambda'*`delta')

* ---------------------------------------------------------------------------
* PRIMARY: aggregate capital, matching the model's Omega = K_core / K_periphery
* ---------------------------------------------------------------------------
gen double lnOmega_agg = ln(Omega_agg)
reg lhs lnOmega_agg
scalar theta_reg_agg    = -_b[lnOmega_agg]
scalar se_reg_agg       = _se[lnOmega_agg]
scalar chi0_reg_agg     = exp(_b[_cons])
scalar r2_reg_agg       = e(r2)
scalar n_reg_agg        = e(N)
scalar theta_reg_agg_lo = theta_reg_agg - invnormal(.975)*se_reg_agg
scalar theta_reg_agg_hi = theta_reg_agg + invnormal(.975)*se_reg_agg

di as result ""
di as result "TABLE 2 REGULAR -- PRIMARY: aggregate capital power"
di as result "theta = " %9.4f scalar(theta_reg_agg) "   SE = " %9.4f scalar(se_reg_agg) ///
    "   95% CI = [" %9.4f scalar(theta_reg_agg_lo) ", " %9.4f scalar(theta_reg_agg_hi) "]"
di as result "chi0 = " %9.4f scalar(chi0_reg_agg) "   R2 = " %9.4f scalar(r2_reg_agg) ///
    "   N = " %9.0f scalar(n_reg_agg)

* Repeated-region robustness: clustered SE. Theta point estimate is identical to OLS.
reg lhs lnOmega_agg, vce(cluster region_id)

scalar se_reg_agg_cl = _se[lnOmega_agg]

scalar df_reg_agg_cl = e(df_r)
scalar tcrit_reg_agg_cl = invttail(df_reg_agg_cl,.025)

scalar lo_reg_agg_cl = ///
    theta_reg_agg - tcrit_reg_agg_cl*se_reg_agg_cl

scalar hi_reg_agg_cl = ///
    theta_reg_agg + tcrit_reg_agg_cl*se_reg_agg_cl


********************************************************************************
* Wild-cluster bootstrap: aggregate Omega
********************************************************************************

boottest lnOmega_agg, ///
    cluster(region_id) ///
    weight(webb) ///
    reps(9999) ///
    seed(12345) ///
    level(95) ///
	jk
	
scalar p_reg_agg_wild = r(p)

matrix CI_reg_agg_wild = r(CI)

scalar lo_reg_agg_wild = -CI_reg_agg_wild[1,2]
scalar hi_reg_agg_wild = -CI_reg_agg_wild[1,1]

di as result "REGULAR SAMPLE: wild-cluster p-value = " ///
    %9.4f scalar(p_reg_agg_wild)

di as result "REGULAR SAMPLE: wild-cluster theta 95% CI = [" ///
    %9.4f scalar(lo_reg_agg_wild) ", " ///
    %9.4f scalar(hi_reg_agg_wild) "]"
* ---------------------------------------------------------------------------
* ROBUSTNESS: capital per worker, matching the coauthor's original figure
* ---------------------------------------------------------------------------
gen double lnOmega_pw = ln(Omega_pw)
reg lhs lnOmega_pw
scalar theta_reg_pw    = -_b[lnOmega_pw]
scalar se_reg_pw       = _se[lnOmega_pw]
scalar chi0_reg_pw     = exp(_b[_cons])
scalar r2_reg_pw       = e(r2)
scalar n_reg_pw        = e(N)
scalar theta_reg_pw_lo = theta_reg_pw - invnormal(.975)*se_reg_pw
scalar theta_reg_pw_hi = theta_reg_pw + invnormal(.975)*se_reg_pw

di as result "TABLE 2 REGULAR -- robustness: capital per worker"
di as result "theta = " %9.4f scalar(theta_reg_pw) "   SE = " %9.4f scalar(se_reg_pw) ///
    "   95% CI = [" %9.4f scalar(theta_reg_pw_lo) ", " %9.4f scalar(theta_reg_pw_hi) "]"

********************************************************************************
* 8. ESTIMATE THETA: ROW-INCLUSIVE TABLE 2
*    PRIMARY SPECIFICATION = AGGREGATE CAPITAL POWER
********************************************************************************

use "$CAL/calib_inputs_1995_2009_table2_ROW.dta", clear
keep if all_year_outflow==1
assert e>0
assert Omega_agg>0
assert Omega_pw>0

gen double lhs = ln(`m'/e + `lambda'*`delta')

* PRIMARY: aggregate capital.
gen double lnOmega_agg = ln(Omega_agg)
reg lhs lnOmega_agg
scalar theta_row_agg    = -_b[lnOmega_agg]
scalar se_row_agg       = _se[lnOmega_agg]
scalar chi0_row_agg     = exp(_b[_cons])
scalar r2_row_agg       = e(r2)
scalar n_row_agg        = e(N)
scalar theta_row_agg_lo = theta_row_agg - invnormal(.975)*se_row_agg
scalar theta_row_agg_hi = theta_row_agg + invnormal(.975)*se_row_agg

di as result ""
di as result "TABLE 2 + ROW -- PRIMARY: aggregate capital power"
di as result "theta = " %9.4f scalar(theta_row_agg) "   SE = " %9.4f scalar(se_row_agg) ///
    "   95% CI = [" %9.4f scalar(theta_row_agg_lo) ", " %9.4f scalar(theta_row_agg_hi) "]"
di as result "chi0 = " %9.4f scalar(chi0_row_agg) "   R2 = " %9.4f scalar(r2_row_agg) ///
    "   N = " %9.0f scalar(n_row_agg)

* ROBUSTNESS: capital per worker.
gen double lnOmega_pw = ln(Omega_pw)
reg lhs lnOmega_pw
scalar theta_row_pw    = -_b[lnOmega_pw]
scalar se_row_pw       = _se[lnOmega_pw]
scalar chi0_row_pw     = exp(_b[_cons])
scalar r2_row_pw       = e(r2)
scalar n_row_pw        = e(N)
scalar theta_row_pw_lo = theta_row_pw - invnormal(.975)*se_row_pw
scalar theta_row_pw_hi = theta_row_pw + invnormal(.975)*se_row_pw

di as result "TABLE 2 + ROW -- robustness: capital per worker"
di as result "theta = " %9.4f scalar(theta_row_pw) "   SE = " %9.4f scalar(se_row_pw) ///
    "   95% CI = [" %9.4f scalar(theta_row_pw_lo) ", " %9.4f scalar(theta_row_pw_hi) "]"

********************************************************************************
* 9. EXPORT THETA RESULTS
********************************************************************************

clear
set obs 4
gen str20 sample = ""
gen str16 power_measure = ""
gen byte primary = .
gen double theta = .
gen double se = .
gen double ci_lo = .
gen double ci_hi = .
gen double ci_lo_cluster=.
gen double ci_hi_cluster=.
gen double p_wild_cluster=.
gen double ci_lo_wild_cluster=.
gen double ci_hi_wild_cluster=.
gen double cluster_df=.
gen double chi0 = .
gen double r2 = .
gen double N = .

* Aggregate specifications first because they are the paper's preferred measure.
replace sample = "Table 2 regular" in 1
replace power_measure = "aggregate" in 1
replace primary = 1 in 1
replace theta = scalar(theta_reg_agg) in 1
replace se = scalar(se_reg_agg) in 1
replace ci_lo = scalar(theta_reg_agg_lo) in 1
replace ci_hi = scalar(theta_reg_agg_hi) in 1
replace p_wild_cluster=scalar(p_reg_agg_wild) in 1
replace ci_lo_wild_cluster=scalar(lo_reg_agg_wild) in 1
replace ci_hi_wild_cluster=scalar(hi_reg_agg_wild) in 1
replace cluster_df=scalar(df_reg_agg_cl) in 1
replace ci_lo_cluster=scalar(lo_reg_agg_cl) in 1
replace ci_hi_cluster=scalar(hi_reg_agg_cl) in 1
replace chi0 = scalar(chi0_reg_agg) in 1
replace r2 = scalar(r2_reg_agg) in 1
replace N = scalar(n_reg_agg) in 1

replace sample = "Table 2 + ROW" in 2
replace power_measure = "aggregate" in 2
replace primary = 1 in 2
replace theta = scalar(theta_row_agg) in 2
replace se = scalar(se_row_agg) in 2
replace ci_lo = scalar(theta_row_agg_lo) in 2
replace ci_hi = scalar(theta_row_agg_hi) in 2
replace chi0 = scalar(chi0_row_agg) in 2
replace r2 = scalar(r2_row_agg) in 2
replace N = scalar(n_row_agg) in 2

replace sample = "Table 2 regular" in 3
replace power_measure = "per_worker" in 3
replace primary = 0 in 3
replace theta = scalar(theta_reg_pw) in 3
replace se = scalar(se_reg_pw) in 3
replace ci_lo = scalar(theta_reg_pw_lo) in 3
replace ci_hi = scalar(theta_reg_pw_hi) in 3
replace chi0 = scalar(chi0_reg_pw) in 3
replace r2 = scalar(r2_reg_pw) in 3
replace N = scalar(n_reg_pw) in 3

replace sample = "Table 2 + ROW" in 4
replace power_measure = "per_worker" in 4
replace primary = 0 in 4
replace theta = scalar(theta_row_pw) in 4
replace se = scalar(se_row_pw) in 4
replace ci_lo = scalar(theta_row_pw_lo) in 4
replace ci_hi = scalar(theta_row_pw_hi) in 4
replace chi0 = scalar(chi0_row_pw) in 4
replace r2 = scalar(r2_row_pw) in 4
replace N = scalar(n_row_pw) in 4

format theta se ci_lo ci_hi chi0 r2 %9.4f
sort primary sample
gsort -primary sample

save "$CAL/theta_estimates_1995_2009.dta", replace
export delimited using "$CAL/theta_estimates_1995_2009.csv", replace

********************************************************************************
* 10. PARAMETER-GRID ROBUSTNESS: lambda x b x delta
*     Uses AGGREGATE Omega, the preferred power measure.
********************************************************************************

tempname gridpost
postfile `gridpost' str20 sample double lambda b delta theta se r2 N using ///
    "$CAL/theta_parameter_grid_1995_2009.dta", replace

foreach SAMPLE in reg row {
    if "`SAMPLE'"=="reg" use "$CAL/calib_inputs_1995_2009_table2_reg.dta", clear
    if "`SAMPLE'"=="row" use "$CAL/calib_inputs_1995_2009_table2_ROW.dta", clear

    keep if all_year_outflow==1
    assert e>0 & Omega_agg>0
    gen double lnOmega_agg = ln(Omega_agg)

    foreach L in 0.2 0.5 0.8 {
        foreach B in 0.3 0.5 0.9 {
            foreach D in 0.3 0.5 0.9 {
                local MM = 1 + `L'*(1-`B')
                capture drop lhs_grid
                gen double lhs_grid = ln(`MM'/e + `L'*`D')
                quietly reg lhs_grid lnOmega_agg
                local TH = -_b[lnOmega_agg]
                local SE = _se[lnOmega_agg]
                local R2 = e(r2)
                local NN = e(N)
                post `gridpost' ("`SAMPLE'") (`L') (`B') (`D') (`TH') (`SE') (`R2') (`NN')
            }
        }
    }
}
postclose `gridpost'

use "$CAL/theta_parameter_grid_1995_2009.dta", clear
sort sample lambda b delta
export delimited using "$CAL/theta_parameter_grid_1995_2009.csv", replace

bysort sample: summarize theta

********************************************************************************
* 11. FIGURE A: REGULAR TABLE 2
*     1995-2009; persistent outflow regions only; AGGREGATE Omega; no Hickel band
********************************************************************************

use "$CAL/calib_inputs_1995_2009_table2_reg.dta", clear
keep if all_year_outflow==1
assert e>0 & Omega_agg>0
sort region year

* Re-estimate the aggregate specification immediately before prediction.
gen double lnOmega_agg = ln(Omega_agg)
gen double lhs = ln(`m'/e + `lambda'*`delta')
reg lhs lnOmega_agg
scalar theta_graph = -_b[lnOmega_agg]

quietly summarize Omega_agg
scalar XMIN = r(min)*0.85
scalar XMAX = r(max)*1.15

tempfile obs_reg grid_reg
save `obs_reg'

clear
set obs 250
gen byte grid = 1
gen double lnOmega_agg = ln(scalar(XMIN)) + (_n-1)/(250-1)*(ln(scalar(XMAX))-ln(scalar(XMIN)))
gen double Omega_agg = exp(lnOmega_agg)

predict double lhs_hat, xb
predict double lhs_se, stdp

gen double lhs_lo = lhs_hat - invnormal(.975)*lhs_se
gen double lhs_hi = lhs_hat + invnormal(.975)*lhs_se

* Back-transform. Since e is decreasing in lhs, the CI endpoints reverse.
gen double e_hat   = `m'/(exp(lhs_hat)-`lambda'*`delta')
gen double e_ci_lo = `m'/(exp(lhs_hi)-`lambda'*`delta')
gen double e_ci_hi = `m'/(exp(lhs_lo)-`lambda'*`delta')

gen str24 region = ""
gen int year = .
gen double e = .
save `grid_reg'

use `obs_reg', clear
gen byte grid = 0
append using `grid_reg'
sort grid region year

local theta_txt : display %4.2f scalar(theta_graph)

#delimit ;
twoway
    (rarea e_ci_hi e_ci_lo Omega_agg if grid==1, sort color(gs10%25) lcolor(gs8) lpattern(dash))
    (line e_hat Omega_agg if grid==1, sort lcolor(black) lwidth(medthick))
    (connected e Omega_agg if grid==0 & region=="East Europe",     pstyle(p1) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Latin America",   pstyle(p2) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="China",           pstyle(p3) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="India",           pstyle(p4) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Other Asia",      pstyle(p5) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="South EMU",       pstyle(p6) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Russia",          pstyle(p7) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Australia",       pstyle(p8) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="North America",   pstyle(p9) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="North EMU",       pstyle(p10) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="North Europe",    pstyle(p11) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="North East Asia", pstyle(p12) msymbol(O) msize(vsmall) lwidth(vthin))
    (scatter e Omega_agg if grid==0 & year==`LASTYEAR', msymbol(none) mlabel(region) mlabsize(small) mlabposition(3))
    ,
    xscale(log)
    xlabel(2 3 5 10 20 30, labsize(small))
    ytitle("Extraction rate  e")
    xtitle("Relative power  {&Omega} (core / region, aggregate capital)")
    legend(order(1 "95% CI" 2 "fitted, {&theta} = `theta_txt'") rows(1) position(6) ring(1) size(small))
    graphregion(color(white))
    plotregion(color(white))
    name(fig_reg, replace)
;
#delimit cr

graph export "$CAL/fig_emp_schedule_1995_2009_aggregate_table2_reg.pdf", replace
graph export "$CAL/fig_emp_schedule_1995_2009_aggregate_table2_reg.png", width(2200) replace

********************************************************************************
* 12. FIGURE B: TABLE 2 + ROW
*     1995-2009; persistent outflow regions only; AGGREGATE Omega; no Hickel band
********************************************************************************

use "$CAL/calib_inputs_1995_2009_table2_ROW.dta", clear
keep if all_year_outflow==1
assert e>0 & Omega_agg>0
sort region year

gen double lnOmega_agg = ln(Omega_agg)
gen double lhs = ln(`m'/e + `lambda'*`delta')
reg lhs lnOmega_agg
scalar theta_graph = -_b[lnOmega_agg]

quietly summarize Omega_agg
scalar XMIN = r(min)*0.85
scalar XMAX = r(max)*1.15

tempfile obs_row grid_row
save `obs_row'

clear
set obs 250
gen byte grid = 1
gen double lnOmega_agg = ln(scalar(XMIN)) + (_n-1)/(250-1)*(ln(scalar(XMAX))-ln(scalar(XMIN)))
gen double Omega_agg = exp(lnOmega_agg)

predict double lhs_hat, xb
predict double lhs_se, stdp

gen double lhs_lo = lhs_hat - invnormal(.975)*lhs_se
gen double lhs_hi = lhs_hat + invnormal(.975)*lhs_se
gen double e_hat   = `m'/(exp(lhs_hat)-`lambda'*`delta')
gen double e_ci_lo = `m'/(exp(lhs_hi)-`lambda'*`delta')
gen double e_ci_hi = `m'/(exp(lhs_lo)-`lambda'*`delta')

gen str24 region = ""
gen int year = .
gen double e = .
save `grid_row'

use `obs_row', clear
gen byte grid = 0
append using `grid_row'
sort grid region year

local theta_txt : display %4.2f scalar(theta_graph)

#delimit ;
twoway
    (rarea e_ci_hi e_ci_lo Omega_agg if grid==1, sort color(gs10%25) lcolor(gs8) lpattern(dash))
    (line e_hat Omega_agg if grid==1, sort lcolor(black) lwidth(medthick))
    (connected e Omega_agg if grid==0 & region=="East Europe",           pstyle(p1) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Latin America",         pstyle(p2) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="China",                 pstyle(p3) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="India",                 pstyle(p4) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Other Asia",            pstyle(p5) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="South EMU",             pstyle(p6) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Russia",                pstyle(p7) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Australia",             pstyle(p8) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="North America",         pstyle(p9) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="North EMU",             pstyle(p10) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="North Europe",          pstyle(p11) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="North East Asia",       pstyle(p12) msymbol(O) msize(vsmall) lwidth(vthin))
    (connected e Omega_agg if grid==0 & region=="Rest of World (est.)",  pstyle(p13) msymbol(O) msize(vsmall) lwidth(vthin))
    (scatter e Omega_agg if grid==0 & year==`LASTYEAR', msymbol(none) mlabel(region) mlabsize(small) mlabposition(3))
    ,
    xscale(log)
    xlabel(2 3 5 10 20 30, labsize(small))
    ytitle("Extraction rate  e")
    xtitle("Relative power  {&Omega} (core / region, aggregate capital)")
    legend(order(1 "95% CI" 2 "fitted, {&theta} = `theta_txt'") rows(1) position(6) ring(1) size(small))
    graphregion(color(white))
    plotregion(color(white))
    name(fig_row, replace)
;
#delimit cr

graph export "$CAL/fig_emp_schedule_1995_2009_aggregate_table2_ROW.pdf", replace
graph export "$CAL/fig_emp_schedule_1995_2009_aggregate_table2_ROW.png", width(2200) replace

********************************************************************************
* 13. LIST REGIONS ACTUALLY USED IN EACH THETA/FIGURE SAMPLE
********************************************************************************

tempname eligpost
postfile `eligpost' str20 sample str24 region using "$CAL/eligible_regions_1995_2009.dta", replace

use "$CAL/calib_inputs_1995_2009_table2_reg.dta", clear
levelsof region if all_year_outflow==1, local(RREG)
foreach r of local RREG {
    post `eligpost' ("Table 2 regular") ("`r'")
}

use "$CAL/calib_inputs_1995_2009_table2_ROW.dta", clear
levelsof region if all_year_outflow==1, local(RROW)
foreach r of local RROW {
    post `eligpost' ("Table 2 + ROW") ("`r'")
}
postclose `eligpost'

use "$CAL/eligible_regions_1995_2009.dta", clear
sort sample region
export delimited using "$CAL/eligible_regions_1995_2009.csv", replace

********************************************************************************
* 14. FINAL MESSAGE
********************************************************************************

di as result ""
di as result "1995-2009 aggregate-power calibration complete."
di as result "Outputs are in: $CAL"
di as result ""
di as result "Main CSVs:"
di as result "  calib_inputs_1995_2009_table2_reg.csv"
di as result "  calib_inputs_1995_2009_table2_ROW.csv"
di as result "  theta_estimates_1995_2009.csv"
di as result "  theta_parameter_grid_1995_2009.csv"
di as result "  eligible_regions_1995_2009.csv"
di as result ""
di as result "Figures (aggregate capital Omega; no Hickel band):"
di as result "  fig_emp_schedule_1995_2009_aggregate_table2_reg.pdf/png"
di as result "  fig_emp_schedule_1995_2009_aggregate_table2_ROW.pdf/png"
