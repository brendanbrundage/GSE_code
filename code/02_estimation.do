/****************************************************************************************
 02_estimation.do

 Global Stratification Economics and International Extraction

 PURPOSE
 -------
 ONE referee-facing estimation file.

 This file begins from the compact Stata calibration datasets produced by
 01_clean_all_big_data.do. It DOES NOT read the original large WIOD, SEA, PPP,
 or PWT source files.

 It runs the original:
   - theta estimates
   - conventional and region-clustered inference
   - Webb-weight wild-cluster bootstrap
   - aggregate-capital main specification
   - capital-per-worker robustness
   - lambda/b/delta parameter-grid robustness
   - empirical extraction-schedule figures
   - optional old-vs-new WIOD overlap check

 ANALYTICAL CODE
 ---------------
 The estimation blocks are copied from the two original calibration do-files.
 No estimating equations, parameters, sample restrictions, bootstrap settings,
 robustness grids, or figure calculations have been changed.

 REQUIREMENT
 -----------
 The original code uses the user-written Stata command boottest.

 Run from the ROOT of the replication repository:
      do "02_estimation.do"
****************************************************************************************/

version 17.0
clear all
set more off
set varabbrev off

global ROOT "`c(pwd)'"


********************************************************************************
* PART 1 SETUP: ORIGINAL WIOD RELEASE, 1995-2009
********************************************************************************

global CAL "$ROOT/calibration_1995_2009_aggregate"

foreach f in ///
    "$CAL/calib_inputs_1995_2009_table2_reg.dta" ///
    "$CAL/calib_inputs_1995_2009_table2_ROW.dta" {
    capture confirm file `f'
    if _rc {
        di as error "Required compact replication file not found: `f'"
        di as error "Run 01_clean_all_big_data.do first, or place the supplied .dta files in the calibration folder."
        exit 601
    }
}

local FIRSTYEAR 1995
local LASTYEAR  2009
local NYEARS = `LASTYEAR' - `FIRSTYEAR' + 1

local lambda 0.5
local b      0.5
local delta  0.5
local m = 1 + `lambda'*(1-`b')

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


/****************************************************************************************
 PART 2 SETUP: WIOD NOVEMBER 2016 RELEASE, 2000-2014
****************************************************************************************/

clear all
set more off
set varabbrev off

global ROOT "`c(pwd)'"
global CAL    "$ROOT/calibration_2016release_2000_2014"
global OUTNEW "$ROOT/output_2016release_2000_2014"
global OUTOLD "$ROOT/output_all_years"

global T2REG    "$OUTNEW/table2_13regions_all43_2000_2014_long.dta"
global T2ROW    "$OUTNEW/table2_14regions_ROW_est_2000_2014_long.dta"
global T2OLDREG "$OUTOLD/table2_12regions_all_years_long.dta"

foreach f in ///
    "$CAL/calib_inputs_2000_2014_table2_reg.dta" ///
    "$CAL/calib_inputs_2000_2014_table2_ROW.dta" {
    capture confirm file `f'
    if _rc {
        di as error "Required compact replication file not found: `f'"
        di as error "Run 01_clean_all_big_data.do first, or place the supplied .dta files in the calibration folder."
        exit 601
    }
}

local FIRSTYEAR 2000
local LASTYEAR  2014
local NYEARS = `LASTYEAR'-`FIRSTYEAR'+1

local lambda 0.5
local b      0.5
local delta  0.5
local m = 1 + `lambda'*(1-`b')

********************************************************************************
* 8. THETA: REGULAR TABLE 2 -- AGGREGATE OMEGA MAIN
********************************************************************************

use "$CAL/calib_inputs_2000_2014_table2_reg.dta", clear
keep if theta_eligible
assert e>0 & Omega_agg>0 & Omega_pw>0

egen region_id=group(region)
gen double lhs=ln(`m'/e + `lambda'*`delta')
gen double lnOmega_agg=ln(Omega_agg)
gen double lnOmega_pw=ln(Omega_pw)

* MAIN: aggregate power, matching Omega = K_C/K_P in the model.
reg lhs lnOmega_agg
scalar theta_reg_agg=-_b[lnOmega_agg]
scalar se_reg_agg=_se[lnOmega_agg]
scalar chi0_reg_agg=exp(_b[_cons])
scalar r2_reg_agg=e(r2)
scalar n_reg_agg=e(N)
scalar lo_reg_agg=theta_reg_agg-invnormal(.975)*se_reg_agg
scalar hi_reg_agg=theta_reg_agg+invnormal(.975)*se_reg_agg

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

* Per-worker robustness.
reg lhs lnOmega_pw
scalar theta_reg_pw=-_b[lnOmega_pw]
scalar se_reg_pw=_se[lnOmega_pw]
scalar chi0_reg_pw=exp(_b[_cons])
scalar r2_reg_pw=e(r2)
scalar n_reg_pw=e(N)
scalar lo_reg_pw=theta_reg_pw-invnormal(.975)*se_reg_pw
scalar hi_reg_pw=theta_reg_pw+invnormal(.975)*se_reg_pw

di as result "NEW WIOD 2000-2014, TABLE 2 REGULAR"
di as result "Aggregate theta = " %9.4f scalar(theta_reg_agg) ///
    "  conventional SE = " %9.4f scalar(se_reg_agg) ///
    "  region-clustered SE = " %9.4f scalar(se_reg_agg_cl)
di as result "Per-worker theta = " %9.4f scalar(theta_reg_pw)

********************************************************************************
* 9. THETA: ROW-INCLUSIVE TABLE 2 -- AGGREGATE OMEGA MAIN
********************************************************************************

use "$CAL/calib_inputs_2000_2014_table2_ROW.dta", clear
keep if theta_eligible
assert e>0 & Omega_agg>0 & Omega_pw>0

egen region_id=group(region)
gen double lhs=ln(`m'/e + `lambda'*`delta')
gen double lnOmega_agg=ln(Omega_agg)
gen double lnOmega_pw=ln(Omega_pw)

reg lhs lnOmega_agg
scalar theta_row_agg=-_b[lnOmega_agg]
scalar se_row_agg=_se[lnOmega_agg]
scalar chi0_row_agg=exp(_b[_cons])
scalar r2_row_agg=e(r2)
scalar n_row_agg=e(N)
scalar lo_row_agg=theta_row_agg-invnormal(.975)*se_row_agg
scalar hi_row_agg=theta_row_agg+invnormal(.975)*se_row_agg

reg lhs lnOmega_agg, vce(cluster region_id)
scalar se_row_agg_cl=_se[lnOmega_agg]
scalar lo_row_agg_cl=theta_row_agg-invnormal(.975)*se_row_agg_cl
scalar hi_row_agg_cl=theta_row_agg+invnormal(.975)*se_row_agg_cl

reg lhs lnOmega_pw
scalar theta_row_pw=-_b[lnOmega_pw]
scalar se_row_pw=_se[lnOmega_pw]
scalar chi0_row_pw=exp(_b[_cons])
scalar r2_row_pw=e(r2)
scalar n_row_pw=e(N)
scalar lo_row_pw=theta_row_pw-invnormal(.975)*se_row_pw
scalar hi_row_pw=theta_row_pw+invnormal(.975)*se_row_pw

di as result "NEW WIOD 2000-2014, TABLE 2 + ROW"
di as result "Aggregate theta = " %9.4f scalar(theta_row_agg) ///
    "  conventional SE = " %9.4f scalar(se_row_agg) ///
    "  region-clustered SE = " %9.4f scalar(se_row_agg_cl)
di as result "Per-worker theta = " %9.4f scalar(theta_row_pw)

********************************************************************************
* 10. EXPORT THETA RESULTS
********************************************************************************

clear
set obs 4
gen str20 sample=""
gen str18 power_measure=""
gen double theta=.
gen double se=.
gen double se_cluster_region=.
gen double ci_lo=.
gen double ci_hi=.
gen double ci_lo_cluster=.
gen double ci_hi_cluster=.
gen double p_wild_cluster=.
gen double ci_lo_wild_cluster=.
gen double ci_hi_wild_cluster=.
gen double cluster_df=.
gen double chi0=.
gen double r2=.
gen double N=.

replace sample="Table 2 regular" in 1
replace power_measure="aggregate" in 1
replace theta=scalar(theta_reg_agg) in 1
replace se=scalar(se_reg_agg) in 1
replace se_cluster_region=scalar(se_reg_agg_cl) in 1
replace ci_lo=scalar(lo_reg_agg) in 1
replace ci_hi=scalar(hi_reg_agg) in 1
replace p_wild_cluster=scalar(p_reg_agg_wild) in 1
replace ci_lo_wild_cluster=scalar(lo_reg_agg_wild) in 1
replace ci_hi_wild_cluster=scalar(hi_reg_agg_wild) in 1
replace cluster_df=scalar(df_reg_agg_cl) in 1
replace ci_lo_cluster=scalar(lo_reg_agg_cl) in 1
replace ci_hi_cluster=scalar(hi_reg_agg_cl) in 1
replace chi0=scalar(chi0_reg_agg) in 1
replace r2=scalar(r2_reg_agg) in 1
replace N=scalar(n_reg_agg) in 1

replace sample="Table 2 regular" in 2
replace power_measure="per_worker" in 2
replace theta=scalar(theta_reg_pw) in 2
replace se=scalar(se_reg_pw) in 2
replace ci_lo=scalar(lo_reg_pw) in 2
replace ci_hi=scalar(hi_reg_pw) in 2
replace chi0=scalar(chi0_reg_pw) in 2
replace r2=scalar(r2_reg_pw) in 2
replace N=scalar(n_reg_pw) in 2

replace sample="Table 2 + ROW" in 3
replace power_measure="aggregate" in 3
replace theta=scalar(theta_row_agg) in 3
replace se=scalar(se_row_agg) in 3
replace se_cluster_region=scalar(se_row_agg_cl) in 3
replace ci_lo=scalar(lo_row_agg) in 3
replace ci_hi=scalar(hi_row_agg) in 3
replace ci_lo_cluster=scalar(lo_row_agg_cl) in 3
replace ci_hi_cluster=scalar(hi_row_agg_cl) in 3
replace chi0=scalar(chi0_row_agg) in 3
replace r2=scalar(r2_row_agg) in 3
replace N=scalar(n_row_agg) in 3

replace sample="Table 2 + ROW" in 4
replace power_measure="per_worker" in 4
replace theta=scalar(theta_row_pw) in 4
replace se=scalar(se_row_pw) in 4
replace ci_lo=scalar(lo_row_pw) in 4
replace ci_hi=scalar(hi_row_pw) in 4
replace chi0=scalar(chi0_row_pw) in 4
replace r2=scalar(r2_row_pw) in 4
replace N=scalar(n_row_pw) in 4

save "$CAL/theta_estimates_2000_2014.dta", replace
export delimited using "$CAL/theta_estimates_2000_2014.csv", replace

********************************************************************************
* 11. PARAMETER GRID: AGGREGATE OMEGA
********************************************************************************

tempname gridpost
postfile `gridpost' str20 sample double lambda b delta theta se r2 N using ///
    "$CAL/theta_parameter_grid_2000_2014.dta", replace

foreach SAMPLE in reg row {
    if "`SAMPLE'"=="reg" use "$CAL/calib_inputs_2000_2014_table2_reg.dta", clear
    if "`SAMPLE'"=="row" use "$CAL/calib_inputs_2000_2014_table2_ROW.dta", clear

    keep if theta_eligible
    gen double lnOmega_agg=ln(Omega_agg)

    foreach L in 0.2 0.5 0.8 {
        foreach B in 0.3 0.5 0.9 {
            foreach D in 0.3 0.5 0.9 {
                local MM=1+`L'*(1-`B')
                capture drop lhs_grid
                gen double lhs_grid=ln(`MM'/e + `L'*`D')
                quietly reg lhs_grid lnOmega_agg

                local TH=-_b[lnOmega_agg]
                local SE=_se[lnOmega_agg]
                local R2=e(r2)
                local NN=e(N)

                post `gridpost' ("`SAMPLE'") (`L') (`B') (`D') ///
                    (`TH') (`SE') (`R2') (`NN')
            }
        }
    }
}
postclose `gridpost'

use "$CAL/theta_parameter_grid_2000_2014.dta", clear
sort sample lambda b delta
export delimited using "$CAL/theta_parameter_grid_2000_2014.csv", replace
bysort sample: summarize theta

********************************************************************************
* 12. FIGURE A: REGULAR TABLE 2 -- AGGREGATE OMEGA, NO HICKEL BAND
********************************************************************************

use "$CAL/calib_inputs_2000_2014_table2_reg.dta", clear
keep if theta_eligible
assert e>0 & Omega_agg>0
sort region year

gen double lnOmega_agg=ln(Omega_agg)
gen double lhs=ln(`m'/e + `lambda'*`delta')
reg lhs lnOmega_agg

scalar theta_graph=-_b[lnOmega_agg]

quietly summarize Omega_agg
scalar XMIN=r(min)*0.85
scalar XMAX=r(max)*1.15

tempfile obs_reg grid_reg
save `obs_reg'

clear
set obs 250
gen byte grid=1
gen double lnOmega_agg=ln(scalar(XMIN))+(_n-1)/(250-1)*(ln(scalar(XMAX))-ln(scalar(XMIN)))
gen double Omega_agg=exp(lnOmega_agg)

predict double lhs_hat, xb
predict double lhs_se, stdp

gen double lhs_lo=lhs_hat-invnormal(.975)*lhs_se
gen double lhs_hi=lhs_hat+invnormal(.975)*lhs_se

gen double e_hat=`m'/(exp(lhs_hat)-`lambda'*`delta')
gen double e_ci_lo=`m'/(exp(lhs_hi)-`lambda'*`delta')
gen double e_ci_hi=`m'/(exp(lhs_lo)-`lambda'*`delta')

gen str24 region=""
gen int year=.
gen double e=.
save `grid_reg'

use `obs_reg', clear
gen byte grid=0
append using `grid_reg'

levelsof region if grid==0, local(regions)
local regionplots ""
local p=1
foreach r of local regions {
    local regionplots `"`regionplots' (connected e Omega_agg if grid==0 & region=="`r'", pstyle(p`p') msymbol(O) msize(vsmall) lwidth(vthin))"'
    local ++p
}

local theta_txt : display %4.2f scalar(theta_graph)

gen yearlabel = ""
replace yearlabel = "2000" if year==2000
replace yearlabel = "2014" if year==2014

#delimit ;
twoway
    (rarea e_ci_hi e_ci_lo Omega_agg if grid==1,
        sort color(gs10%25) lcolor(gs8) lpattern(dash))

    (line e_hat Omega_agg if grid==1,
        sort lcolor(black) lwidth(medthick))

    `regionplots'

    /* Label beginning and end years */
    (scatter e Omega_agg if grid==0 & inlist(year,2000,2014),
        msymbol(O)
        msize(vsmall)
        mlabel(yearlabel)
        mlabsize(vsmall)
        mlabposition(12))

    /* Region name at 2014 endpoint */
    (scatter e Omega_agg if grid==0 & year==`LASTYEAR',
        msymbol(none)
        mlabel(region)
        mlabsize(small)
        mlabposition(3))
    ,
    xscale(log)
    ytitle("Extraction rate  e")
    xtitle("Relative power {&Omega} (core / region, aggregate capital)")
    legend(order(1 "95% CI" 2 "fitted, {&theta} = `theta_txt'") rows(1) size(small))
    graphregion(color(white))
    plotregion(margin(small))
    name(fig_reg_newwiod, replace)
;
#delimit cr

graph export "$CAL/fig_emp_schedule_2000_2014_table2_reg_aggregate.pdf", replace
graph export "$CAL/fig_emp_schedule_2000_2014_table2_reg_aggregate.png", width(2200) replace

********************************************************************************
* 13. FIGURE B: TABLE 2 + ROW -- AGGREGATE OMEGA, NO HICKEL BAND
********************************************************************************

use "$CAL/calib_inputs_2000_2014_table2_ROW.dta", clear
keep if theta_eligible
assert e>0 & Omega_agg>0
sort region year

gen double lnOmega_agg=ln(Omega_agg)
gen double lhs=ln(`m'/e + `lambda'*`delta')
reg lhs lnOmega_agg

scalar theta_graph=-_b[lnOmega_agg]

quietly summarize Omega_agg
scalar XMIN=r(min)*0.85
scalar XMAX=r(max)*1.15

tempfile obs_row grid_row
save `obs_row'

clear
set obs 250
gen byte grid=1
gen double lnOmega_agg=ln(scalar(XMIN))+(_n-1)/(250-1)*(ln(scalar(XMAX))-ln(scalar(XMIN)))
gen double Omega_agg=exp(lnOmega_agg)

predict double lhs_hat, xb
predict double lhs_se, stdp

gen double lhs_lo=lhs_hat-invnormal(.975)*lhs_se
gen double lhs_hi=lhs_hat+invnormal(.975)*lhs_se

gen double e_hat=`m'/(exp(lhs_hat)-`lambda'*`delta')
gen double e_ci_lo=`m'/(exp(lhs_hi)-`lambda'*`delta')
gen double e_ci_hi=`m'/(exp(lhs_lo)-`lambda'*`delta')

gen str24 region=""
gen int year=.
gen double e=.
save `grid_row'

use `obs_row', clear
gen byte grid=0
append using `grid_row'

levelsof region if grid==0, local(regions)
local regionplots ""
local p=1
foreach r of local regions {
    local regionplots `"`regionplots' (connected e Omega_agg if grid==0 & region=="`r'", pstyle(p`p') msymbol(O) msize(vsmall) lwidth(vthin))"'
    local ++p
}

local theta_txt : display %4.2f scalar(theta_graph)

gen yearlabel = ""
replace yearlabel = "2000" if year==2000
replace yearlabel = "2014" if year==2014

#delimit ;
twoway
    (rarea e_ci_hi e_ci_lo Omega_agg if grid==1,
        sort color(gs10%25) lcolor(gs8) lpattern(dash))

    (line e_hat Omega_agg if grid==1,
        sort lcolor(black) lwidth(medthick))

    `regionplots'

    /* Label beginning and end years */
    (scatter e Omega_agg if grid==0 & inlist(year,2000,2014),
        msymbol(O)
        msize(vsmall)
        mlabel(yearlabel)
        mlabsize(vsmall)
        mlabposition(12))

    /* Region name at 2014 endpoint */
    (scatter e Omega_agg if grid==0 & year==`LASTYEAR',
        msymbol(none)
        mlabel(region)
        mlabsize(small)
        mlabposition(3))
    ,
    xscale(log)
    ytitle("Extraction rate  e")
    xtitle("Relative power {&Omega} (core / region, aggregate capital)")
    legend(order(1 "95% CI" 2 "fitted, {&theta} = `theta_txt'") rows(1) size(small))
    graphregion(color(white))
    plotregion(margin(small))
    name(fig_row_newwiod, replace)
;
#delimit cr

graph export "$CAL/fig_emp_schedule_2000_2014_table2_ROW_aggregate.pdf", replace
graph export "$CAL/fig_emp_schedule_2000_2014_table2_ROW_aggregate.png", width(2200) replace

********************************************************************************
* 14. OPTIONAL OLD-vs-NEW WIOD OVERLAP CHECK, 2000-2009
********************************************************************************

capture confirm file "$T2OLDREG"
if !_rc {
    use "$T2OLDREG", clear
    keep if inrange(year,2000,2009)
    keep year region net_transfer net_transfer_pct_va
    rename net_transfer old_net_transfer
    rename net_transfer_pct_va old_net_transfer_pct_va
    tempfile old_overlap
    save `old_overlap'

    use "$T2REG", clear
    keep if inrange(year,2000,2009)
    drop if region=="Additional Core"
    keep year region net_transfer net_transfer_pct_va
    rename net_transfer new_net_transfer
    rename net_transfer_pct_va new_net_transfer_pct_va

    merge 1:1 year region using `old_overlap'
    assert _merge==3
    drop _merge

    gen double old_e=-old_net_transfer_pct_va/100
    gen double new_e=-new_net_transfer_pct_va/100
    gen double diff_e=new_e-old_e
    gen double abs_diff_e=abs(diff_e)
    gen byte same_transfer_sign=sign(new_net_transfer)==sign(old_net_transfer)

    sort region year
    save "$CAL/old_vs_new_extraction_2000_2009.dta", replace
    export delimited using "$CAL/old_vs_new_extraction_2000_2009.csv", replace

    quietly corr old_e new_e
   matrix C = r(C)
scalar overlap_corr = C[1,2]
    quietly summarize abs_diff_e, meanonly
    scalar overlap_mae=r(mean)
    quietly summarize same_transfer_sign, meanonly
    scalar overlap_sign_agreement=r(mean)

    di as result "OLD vs NEW WIOD, 2000-2009:"
    di as result "Correlation of regional extraction rates = " %9.4f scalar(overlap_corr)
    di as result "Mean absolute extraction-rate difference = " %9.4f scalar(overlap_mae)
    di as result "Transfer-sign agreement = " %9.4f scalar(overlap_sign_agreement)
}
else {
    di as text "Old-WIOD overlap file not found; skipping optional 2000-2009 release comparison."
}

********************************************************************************
* 15. FINAL MESSAGE
********************************************************************************

di as result "Completed WIOD-2016 calibration for 2000-2014."
di as result "MAIN theta uses aggregate Omega = K_core/K_region."
di as result "Core includes original four core regions PLUS CHE, HRV, and NOR."
di as result "CHE, HRV, NOR are grouped as Additional Core and never treated as ROW."
di as result "Only non-core regions with negative net transfers in every year 2000-2014 enter theta/figures."
di as result "No Hickel band is included."
di as result "Outputs: $CAL"
