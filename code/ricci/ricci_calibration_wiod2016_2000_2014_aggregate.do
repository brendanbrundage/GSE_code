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
 RICCI / PWT CALIBRATION USING WIOD NOVEMBER 2016 RELEASE, 2000-2014

 MAIN DESIGN
   - Uses ONLY the WIOD 2016 release for 2000-2014. No hybrid splice.
   - Main power measure matches the model:
         Omega = K_core / K_region
   - Per-worker Omega is retained as a robustness check.
   - Extraction:
         e = - net_transfer_pct_va / 100
   - Estimation sample contains only NON-CORE regions with net_transfer < 0
     in EVERY year 2000-2014.
   - Figures contain only those persistent extraction-provider regions.
   - No Hickel band.

 CORE
   Original coauthor core:
       North America + North EMU + North Europe + North East Asia
   plus the three additional countries separately identified by WIOD 2016:
       Switzerland (CHE) + Croatia (HRV) + Norway (NOR)

   CHE + HRV + NOR are grouped as "Additional Core" in the Table 2/PWT crosswalk.
   They are NEVER treated as ROW and NEVER enter the peripheral theta sample.

 MODEL / TRANSFORMED ESTIMATING EQUATION
   e = m / (chi0*Omega^(-theta) - lambda*delta)
   m = 1 + lambda*(1-b)

   ln(m/e + lambda*delta) = ln(chi0) - theta*ln(Omega)

 INPUTS
   output_2016release_2000_2014/table2_13regions_all43_2000_2014_long.dta
   output_2016release_2000_2014/table2_14regions_ROW_est_2000_2014_long.dta
   raw/pwt110.dta

 OPTIONAL OVERLAP CHECK
   If output_all_years/table2_12regions_all_years_long.dta exists, the code
   also compares old-WIOD and new-WIOD extraction rates over 2000-2009 for
   the 12 common Ricci regions.

 OUTPUTS
   calibration_2016release_2000_2014/
       calib_inputs_2000_2014_table2_reg.csv
       calib_inputs_2000_2014_table2_ROW.csv
       eligible_regions_2000_2014.csv
       theta_estimates_2000_2014.csv
       theta_parameter_grid_2000_2014.csv
       fig_emp_schedule_2000_2014_table2_reg_aggregate.pdf/png
       fig_emp_schedule_2000_2014_table2_ROW_aggregate.pdf/png
       old_vs_new_extraction_2000_2009.csv   [if old output exists]
****************************************************************************************/

version 17.0
clear all
set more off
set varabbrev off

********************************************************************************
* 0. PATHS / PARAMETERS
********************************************************************************

global ROOT "`c(pwd)'"
global RAW    "$ROOT/raw"
global OUTNEW "$ROOT/output_2016release_2000_2014"
global OUTOLD "$ROOT/output_all_years"
global CAL    "$ROOT/calibration_2016release_2000_2014"

capture mkdir "$CAL"

global PWT      "$RAW/pwt110.dta"
global T2REG    "$OUTNEW/table2_13regions_all43_2000_2014_long.dta"
global T2ROW    "$OUTNEW/table2_14regions_ROW_est_2000_2014_long.dta"
global T2OLDREG "$OUTOLD/table2_12regions_all_years_long.dta"

foreach f in "$PWT" "$T2REG" "$T2ROW" {
    capture confirm file `f'
    if _rc {
        di as error "Required file not found: `f'"
        di as error "Run ricci_wiod2016_2000_2014_all43.do first."
        exit 601
    }
}

local FIRSTYEAR 2000
local LASTYEAR  2014
local NYEARS = `LASTYEAR'-`FIRSTYEAR'+1

* Coauthor baseline parameters.
local lambda 0.5
local b      0.5
local delta  0.5
local m = 1 + `lambda'*(1-`b')

* All 43 countries separately observed in WIOD 2016.
local CTRY43 "AUS AUT BEL BGR BRA CAN CHE CHN CYP CZE DEU DNK ESP EST FIN FRA GBR GRC HRV HUN IDN IND IRL ITA JPN KOR LTU LUX LVA MEX MLT NLD NOR POL PRT ROU RUS SVK SVN SWE TUR TWN USA"

********************************************************************************
* 1. PWT 11.0 COUNTRY DATA + WIOD-2016 REGION CROSSWALK
********************************************************************************

use "$PWT", clear
keep if inrange(year,`FIRSTYEAR',`LASTYEAR')
keep countrycode country year cn emp cgdpo labsh csh_i
rename countrycode iso3

gen byte in_wiod43 = strpos(" `CTRY43' ", " " + iso3 + " ")>0
count if in_wiod43 & missing(cn,emp,cgdpo,labsh,csh_i)
assert r(N)==0

gen str24 region=""
replace region="North America"   if inlist(iso3,"CAN","USA")
replace region="North EMU"       if inlist(iso3,"AUT","BEL","DEU","FIN","FRA","LUX","NLD")
replace region="South EMU"       if inlist(iso3,"CYP","ESP","GRC","IRL","ITA","MLT","PRT")
replace region="North Europe"    if inlist(iso3,"DNK","GBR","SWE")
replace region="East Europe"     if inlist(iso3,"BGR","CZE","EST","HUN","LTU")
replace region="East Europe"     if inlist(iso3,"LVA","POL","ROU","SVK","SVN")
replace region="Latin America"   if inlist(iso3,"BRA","MEX")
replace region="China"           if iso3=="CHN"
replace region="India"           if iso3=="IND"
replace region="North East Asia" if inlist(iso3,"JPN","KOR")
replace region="Other Asia"      if inlist(iso3,"IDN","TUR","TWN")
replace region="Russia"          if iso3=="RUS"
replace region="Australia"       if iso3=="AUS"

* User-specified treatment: CHE, HRV, NOR are observed CORE, not ROW.
replace region="Additional Core" if inlist(iso3,"CHE","HRV","NOR")

gen double labor_income_proxy=labsh*cgdpo
gen double investment_proxy=csh_i*cgdpo

gen byte core_country = ///
    inlist(region,"North America","North EMU","North Europe","North East Asia") | ///
    region=="Additional Core"

tempfile pwt_country
save `pwt_country'

********************************************************************************
* 2. CORE AGGREGATES: ORIGINAL CORE + CHE/HRV/NOR
********************************************************************************

preserve
    keep if core_country
    collapse (sum) K_core=cn L_core=emp Y_core=cgdpo ///
                   labor_core=labor_income_proxy inv_core=investment_proxy, by(year)
    gen double alpha_core=labor_core/Y_core
    gen double s_core=inv_core/Y_core
    keep year K_core L_core Y_core alpha_core s_core
    isid year
    count
    assert r(N)==`NYEARS'
    tempfile core_agg
    save `core_agg'
restore

********************************************************************************
* 3. PWT AGGREGATES: 13 OBSERVED TABLE-2 REGIONS
********************************************************************************

preserve
    keep if in_wiod43
    assert region!=""

    collapse (sum) K_per=cn L_per=emp Y_per=cgdpo ///
                   labor_per=labor_income_proxy inv_per=investment_proxy, by(year region)

    gen double alpha_per=labor_per/Y_per
    gen double s_per=inv_per/Y_per

    merge m:1 year using `core_agg'
    assert _merge==3
    drop _merge

    gen double Omega_agg=K_core/K_per
    gen double Omega_pw=(K_core/L_core)/(K_per/L_per)
    gen double Y_ratio_agg=Y_core/Y_per
    gen double Y_ratio_pw=(Y_core/L_core)/(Y_per/L_per)
    gen double sigma=s_core/s_per

    keep year region Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
         alpha_core alpha_per sigma K_core K_per L_core L_per
    isid year region
    count
    assert r(N)==13*`NYEARS'

    tempfile pwt13
    save `pwt13'
restore

********************************************************************************
* 4. PWT PROXY FOR CURRENT WIOD ROW
********************************************************************************

use `pwt_country', clear
keep if !in_wiod43
keep if !missing(cn,emp,cgdpo,labsh,csh_i)

collapse (sum) K_per=cn L_per=emp Y_per=cgdpo ///
               labor_per=labor_income_proxy inv_per=investment_proxy ///
         (count) row_pwt_n_countries=cn, by(year)

gen str24 region="Rest of World (est.)"
gen double alpha_per=labor_per/Y_per
gen double s_per=inv_per/Y_per

merge 1:1 year using `core_agg'
assert _merge==3
drop _merge

gen double Omega_agg=K_core/K_per
gen double Omega_pw=(K_core/L_core)/(K_per/L_per)
gen double Y_ratio_agg=Y_core/Y_per
gen double Y_ratio_pw=(Y_core/L_core)/(Y_per/L_per)
gen double sigma=s_core/s_per

keep year region Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
     alpha_core alpha_per sigma K_core K_per L_core L_per row_pwt_n_countries

isid year
count
assert r(N)==`NYEARS'

tempfile pwtrow
save `pwtrow'

use `pwt13', clear
append using `pwtrow'
isid year region
tempfile pwt14
save `pwt14'

********************************************************************************
* 5. CALIBRATION INPUTS: REGULAR TABLE 2
********************************************************************************

use "$T2REG", clear
keep if inrange(year,`FIRSTYEAR',`LASTYEAR')
merge 1:1 year region using `pwt13'
assert _merge==3
drop _merge

gen double e=-net_transfer_pct_va/100
gen byte outflow=net_transfer<0

gen byte core_component = ///
    inlist(region,"North America","North EMU","North Europe","North East Asia") | ///
    region=="Additional Core"

bysort region: egen byte all_year_outflow=min(outflow)
bysort region: egen int n_region_years=count(year)
replace all_year_outflow=0 if n_region_years!=`NYEARS'

* A core component never enters the peripheral calibration even if its transfer sign is negative.
gen byte theta_eligible=(all_year_outflow==1 & core_component==0)

order year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
      alpha_core alpha_per sigma K_core K_per L_core L_per ///
      net_transfer net_transfer_pct_va outflow all_year_outflow ///
      theta_eligible core_component
sort region year

save "$CAL/calib_inputs_2000_2014_table2_reg.dta", replace

preserve
    keep year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
         alpha_core alpha_per sigma K_core K_per L_core L_per
    sort year region
    export delimited using "$CAL/calib_inputs_2000_2014_table2_reg.csv", replace
restore

********************************************************************************
* 6. CALIBRATION INPUTS: TABLE 2 + ESTIMATED ROW
********************************************************************************

use "$T2ROW", clear
keep if inrange(year,`FIRSTYEAR',`LASTYEAR')
merge 1:1 year region using `pwt14'
assert _merge==3
drop _merge

gen double e=-net_transfer_pct_va/100
gen byte outflow=net_transfer<0

gen byte core_component = ///
    inlist(region,"North America","North EMU","North Europe","North East Asia") | ///
    region=="Additional Core"

gen byte row_estimated_calibration=region=="Rest of World (est.)"

bysort region: egen byte all_year_outflow=min(outflow)
bysort region: egen int n_region_years=count(year)
replace all_year_outflow=0 if n_region_years!=`NYEARS'
gen byte theta_eligible=(all_year_outflow==1 & core_component==0)

order year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
      alpha_core alpha_per sigma K_core K_per L_core L_per ///
      net_transfer net_transfer_pct_va outflow all_year_outflow ///
      theta_eligible core_component row_estimated_calibration
sort region year

save "$CAL/calib_inputs_2000_2014_table2_ROW.dta", replace

preserve
    keep year region e Omega_agg Omega_pw Y_ratio_agg Y_ratio_pw ///
         alpha_core alpha_per sigma K_core K_per L_core L_per
    sort year region
    export delimited using "$CAL/calib_inputs_2000_2014_table2_ROW.csv", replace
restore

********************************************************************************
* 7. ELIGIBLE REGION LIST
********************************************************************************

use "$CAL/calib_inputs_2000_2014_table2_reg.dta", clear
preserve
    keep if theta_eligible
    keep region
    duplicates drop
    sort region
    gen str20 sample="Table 2 regular"
    tempfile eligreg
    save `eligreg'
restore

use "$CAL/calib_inputs_2000_2014_table2_ROW.dta", clear
keep if theta_eligible
keep region
duplicates drop
sort region
gen str20 sample="Table 2 + ROW"
append using `eligreg'
sort sample region
save "$CAL/eligible_regions_2000_2014.dta", replace
export delimited using "$CAL/eligible_regions_2000_2014.csv", replace

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
