/****************************************************************************************
 01_clean_all_big_data.do

 Global Stratification Economics and International Extraction

 PURPOSE
 -------
 ONE local data-construction file.

 Run this file on the computer that contains the original large WIOD / SEA /
 PPP / PWT source files. It reproduces the original Ricci-style data
 construction for both WIOD vintages and then creates the compact calibration
 .dta files used by the estimation do-file.

 NO ESTIMATION IS RUN IN THIS FILE.

 IMPORTANT
 ---------
 The analytical code below is copied from the original four do-files. No
 country definitions, industry definitions, sample restrictions, imputations,
 formulas, transfer calculations, PWT aggregation rules, or calibration-data
 construction rules have been changed.

 The only portability change is:
     global ROOT "`c(pwd)'"

 Therefore, before running this file:
   1. Put this do-file in the ROOT of the Ricci_replication folder.
   2. Put the original large source files in ROOT/raw using the same filenames
      expected by the original code.
   3. In Stata, cd to ROOT.
   4. Run:
          do "01_clean_all_big_data.do"

 The compact files needed by 02_estimation.do are created in:
   calibration_1995_2009_aggregate/
      calib_inputs_1995_2009_table2_reg.dta
      calib_inputs_1995_2009_table2_ROW.dta

   calibration_2016release_2000_2014/
      calib_inputs_2000_2014_table2_reg.dta
      calib_inputs_2000_2014_table2_ROW.dta
      eligible_regions_2000_2014.dta

 Those compact .dta files, not the huge WIOD/SEA files, are the data files to
 upload with the referee-facing replication package.
****************************************************************************************/



/****************************************************************************************
 PART 1. ORIGINAL WIOD RELEASE: RAW WIOD/SEA/PPP -> RICCI TABLE 2 DATA
 Source: ricci_all_years_with_ROW_estimate_v2.do
****************************************************************************************/

/****************************************************************************************
 Ricci (2019), "Unequal Exchange in the Age of Globalization"
 ALL-YEARS WIOD extension: 1995-2011
 Includes:
   A. Ricci-comparable 40-country results for every WIOD year
   B. Extended 40 countries + estimated Rest of World (ROW)

 INPUTS
   raw/wiot_full.dta
   raw/Socio_Economic_Accounts_July14.xlsx
   raw/P_Data_Extract_From_World_Development_Indicators.xlsx
   raw/dataset_2026-08-04T19_36_35.328252551Z_DEFAULT_INTEGRATION_IMF.RES_WEO_9.0.0.csv

 IMPORTANT
   - 1995-2009 use observed SEA LAB and H_EMP except structural-zero cells.
   - 2010-2011 are not fully covered by SEA labor variables. Where H_EMP is missing,
     hours are estimated as current EMP times the country's 2009 industry hours/worker.
     Where LAB is missing, labor compensation is estimated as current VA times the
     country's 2009 industry labor-compensation share.
   - ROW has WIOT trade, value added and gross output but no SEA labor accounts.
     Therefore ROW transfer rates are ESTIMATED, by year and industry, as the
     gross-output-weighted average transfer-per-dollar rate of the observed 40 countries.
     This keeps the extended 41-region bilateral accounting zero-sum while leaving the
     original Ricci 40-country coefficients untouched.
   - Industry 35 (private households with employed persons) is excluded, as in Ricci.

 MAIN OUTPUT
   output/ricci_all_years_1995_2011.xlsx
****************************************************************************************/

version 17.0
clear all
set more off
set varabbrev off
set maxvar 32767

********************************************************************************
* 0. USER PATHS
********************************************************************************

global ROOT "`c(pwd)'"

global RAW  "$ROOT/raw"
global WORK "$ROOT/work_all_years"
global OUT  "$ROOT/output_all_years"

capture mkdir "$WORK"
capture mkdir "$OUT"

global WIOT     "$RAW/wiot_full.dta"
global SEA_XLSX "$RAW/Socio_Economic_Accounts_July14.xlsx"
global PPP_XLSX "$RAW/P_Data_Extract_From_World_Development_Indicators.xlsx"
global IMF_CSV  "$RAW/dataset_2026-08-04T19_36_35.328252551Z_DEFAULT_INTEGRATION_IMF.RES_WEO_9.0.0.csv"

foreach f in "$WIOT" "$SEA_XLSX" "$PPP_XLSX" "$IMF_CSV" {
    capture confirm file `f'
    if _rc {
        di as error "Required input file not found: `f'"
        exit 601
    }
}

********************************************************************************
* 1. SAMPLE AND YEARS
********************************************************************************

local CTRY40 "AUS AUT BEL BGR BRA CAN CHN CYP CZE DEU DNK ESP EST FIN FRA GBR GRC HUN IDN IND IRL ITA JPN KOR LTU LUX LVA MEX MLT NLD POL PRT ROU RUS SVK SVN SWE TUR TWN USA"
local CTRY41 "`CTRY40' ROW"

local FIRSTYEAR 1995
local LASTYEAR  2011

********************************************************************************
* 2. REDUCE WIOT TO 1995-2011
********************************************************************************

use year row_country col_country row_item col_item value ///
    if inrange(year,`FIRSTYEAR',`LASTYEAR') using "$WIOT", clear

* Standardize WIOD's mixed-case Rest-of-World code.
replace row_country = "ROW" if row_country=="RoW"
replace col_country = "ROW" if col_country=="RoW"

compress
save "$WORK/wiot_1995_2011.dta", replace

********************************************************************************
* 3. COUNTRY-INDUSTRY VALUE ADDED AND GROSS OUTPUT: 40 + ROW
********************************************************************************

use "$WORK/wiot_1995_2011.dta", clear
keep if inlist(row_country,"VA","GO")

gen byte keep_country = strpos(" `CTRY41' ", " " + col_country + " ") > 0
keep if keep_country
keep if inrange(col_item,1,34)

keep year row_country col_country col_item value
rename col_country country
rename col_item industry

isid year country industry row_country
reshape wide value, i(year country industry) j(row_country) string
rename valueVA va_usd
rename valueGO go_usd

label var va_usd "Value added, million current US$"
label var go_usd "Gross output, million current US$"

assert !missing(va_usd)
assert !missing(go_usd)

compress
save "$WORK/wiod_accounts_41.dta", replace

preserve
    keep if country!="ROW"
    isid year country industry
    save "$WORK/wiod_accounts_40.dta", replace
restore

********************************************************************************
* 4. BILATERAL TRADE: 40 + ROW
********************************************************************************

use "$WORK/wiot_1995_2011.dta", clear

gen byte keep_exporter = strpos(" `CTRY41' ", " " + row_country + " ") > 0
gen byte keep_importer = strpos(" `CTRY41' ", " " + col_country + " ") > 0
keep if keep_exporter & keep_importer
keep if inrange(row_item,1,34)

* Intermediate-use columns are 1-35. Final-demand columns in this WIOD Stata
* file are 37, 38, 39, 41, and 42.
gen byte valid_use = inrange(col_item,1,35) | ///
                     inlist(col_item,37,38,39,41,42)
keep if valid_use

drop if row_country == col_country

collapse (sum) trade_usd=value, ///
    by(year row_country col_country row_item)

rename row_country exporter
rename col_country importer
rename row_item industry

label var trade_usd "Bilateral trade, million current US$"
isid year exporter importer industry
compress
save "$WORK/bilateral_trade_41.dta", replace

preserve
    keep if exporter!="ROW" & importer!="ROW"
    save "$WORK/bilateral_trade_40.dta", replace
restore

********************************************************************************
* 5. WIOD INDUSTRY CROSSWALK
********************************************************************************

clear
input str8 sea_code byte industry str70 industry_name
"AtB"   1  "Agriculture, hunting, forestry and fishing"
"C"     2  "Mining and quarrying"
"15t16" 3  "Food, beverages and tobacco"
"17t18" 4  "Textiles and textile products"
"19"    5  "Leather, leather products and footwear"
"20"    6  "Wood and products of wood and cork"
"21t22" 7  "Pulp, paper, printing and publishing"
"23"    8  "Coke, refined petroleum and nuclear fuel"
"24"    9  "Chemicals and chemical products"
"25"   10  "Rubber and plastics"
"26"   11  "Other non-metallic mineral products"
"27t28" 12 "Basic metals and fabricated metal"
"29"   13  "Machinery, nec"
"30t33" 14 "Electrical and optical equipment"
"34t35" 15 "Transport equipment"
"36t37" 16 "Manufacturing nec; recycling"
"E"    17  "Electricity, gas and water supply"
"F"    18  "Construction"
"50"   19  "Motor vehicle trade and repair; retail fuel"
"51"   20  "Wholesale trade"
"52"   21  "Retail trade and repair of household goods"
"H"    22  "Hotels and restaurants"
"60"   23  "Inland transport"
"61"   24  "Water transport"
"62"   25  "Air transport"
"63"   26  "Supporting transport and travel agencies"
"64"   27  "Post and telecommunications"
"J"    28  "Financial intermediation"
"70"   29  "Real estate activities"
"71t74" 30 "Renting and other business activities"
"L"    31  "Public administration and defence"
"M"    32  "Education"
"N"    33  "Health and social work"
"O"    34  "Other community, social and personal services"
"P"    35  "Private households with employed persons"
end

isid sea_code
save "$WORK/wiod_industry_map.dta", replace
export delimited using "$OUT/wiod_industry_map.csv", replace

********************************************************************************
* 6. IMPORT SEA: VA, LAB, H_EMP, AND EMP FOR ALL YEARS
********************************************************************************

import excel using "$SEA_XLSX", sheet("DATA") allstring clear

rename (A B C D E F G H I J K L M N O P Q R S T U) ///
       (country variable description sea_code ///
        sea1995 sea1996 sea1997 sea1998 sea1999 sea2000 sea2001 sea2002 ///
        sea2003 sea2004 sea2005 sea2006 sea2007 sea2008 sea2009 sea2010 sea2011)

replace country  = upper(strtrim(country))
replace variable = upper(strtrim(variable))
replace sea_code = strtrim(sea_code)

keep if inlist(variable,"VA","LAB","H_EMP","EMP")
gen byte keep_country = strpos(" `CTRY40' ", " " + country + " ") > 0
keep if keep_country

foreach y of numlist 1995/2011 {
    destring sea`y', replace force
}

merge m:1 sea_code using "$WORK/wiod_industry_map.dta", keepusing(industry)
keep if _merge==3
drop _merge
keep if inrange(industry,1,34)

reshape long sea, i(country variable industry) j(year)
rename sea sea_value
keep if inrange(year,`FIRSTYEAR',`LASTYEAR')

isid year country industry variable
reshape wide sea_value, i(year country industry) j(variable) string

rename sea_valueVA    va_ncu
rename sea_valueLAB   lab_ncu
rename sea_valueH_EMP hours
rename sea_valueEMP   employment

label var va_ncu     "SEA value added, million national currency"
label var lab_ncu    "SEA labour compensation, million national currency"
label var hours      "SEA total hours worked by persons engaged, millions"
label var employment "SEA persons engaged, thousands"

isid year country industry
compress
save "$WORK/sea_all_years_raw.dta", replace

********************************************************************************
* 7. MERGE SEA TO WIOT; HANDLE STRUCTURAL ZEROS; IMPUTE 2010-2011 SEA GAPS
********************************************************************************

use "$WORK/wiod_accounts_40.dta", clear
merge 1:1 year country industry using "$WORK/sea_all_years_raw.dta"
assert _merge==3
drop _merge

* Structural-zero cells: set missing SEA values to zero only if both WIOT GO and VA
* are zero. This handles the persistent CHN/IDN zero cells without inventing activity.
foreach v in va_ncu lab_ncu hours employment {
    replace `v' = 0 if missing(`v') & abs(go_usd)<1e-10 & abs(va_usd)<1e-10
}

* Construct 2009 country-industry benchmarks used ONLY where 2010-2011 labor data
* are absent. EMP and VA are substantially more complete in 2010-2011.
gen double hpw_2009_tmp = hours/employment if year==2009 & hours>0 & employment>0
bysort country industry: egen double hpw_2009 = max(hpw_2009_tmp)

gen double labshare_2009_tmp = lab_ncu/va_ncu if year==2009 & lab_ncu>=0 & va_ncu>0
bysort country industry: egen double labshare_2009 = max(labshare_2009_tmp)

gen byte hours_imputed = 0
replace hours_imputed = 1 if inrange(year,2010,2011) & missing(hours) & ///
                             employment>0 & !missing(hpw_2009)
replace hours = employment*hpw_2009 if hours_imputed

gen byte lab_imputed = 0
replace lab_imputed = 1 if inrange(year,2010,2011) & missing(lab_ncu) & ///
                           va_ncu>0 & !missing(labshare_2009)
replace lab_ncu = va_ncu*labshare_2009 if lab_imputed

* Export every imputed cell so 2010-2011 can always be separated from the fully
* observed 1995-2009 results.
preserve
    keep if hours_imputed | lab_imputed
    keep year country industry hours_imputed lab_imputed ///
         hours employment hpw_2009 lab_ncu va_ncu labshare_2009
    sort year country industry
    save "$OUT/sea_2010_2011_imputations.dta", replace
    export delimited using "$OUT/sea_2010_2011_imputations.csv", replace
restore

drop hpw_2009_tmp labshare_2009_tmp

* Nothing active should remain MISSING after the transparent 2010-2011 imputation.
count if go_usd>1e-10 & missing(va_ncu)
assert r(N)==0
count if go_usd>1e-10 & missing(hours)
assert r(N)==0
count if go_usd>1e-10 & missing(lab_ncu)
assert r(N)==0

* The all-year WIOD panel contains a handful of tiny/statistical source-sector
* anomalies with positive gross output but nonpositive value added and/or zero
* recorded hours. These do not occur in Ricci's three benchmark years, except
* outside those benchmark observations. Ricci's homogeneous-labor coefficient
* is not economically defined for nonpositive VA, so DO NOT force positive
* hours or VA. Flag and export these cells, then exclude their exporter-sector
* flows from the transfer-rate calculation below.
gen byte ricci_ineligible = go_usd>1e-10 & (va_usd<=1e-10 | hours<=0)

preserve
    keep if ricci_ineligible
    keep year country industry go_usd va_usd va_ncu hours employment lab_ncu ///
         hours_imputed lab_imputed
    sort year country industry
    save "$OUT/ricci_ineligible_source_cells.dta", replace
    export delimited using "$OUT/ricci_ineligible_source_cells.csv", replace
restore

* For all otherwise usable active sectors, Ricci requires positive VA and labor.
assert hours>0 if go_usd>1e-10 & !ricci_ineligible
assert va_usd>0 if go_usd>1e-10 & !ricci_ineligible

* Derive WIOD's current market exchange rate, LCU per current US dollar.
bysort year country: egen double total_va_ncu = total(va_ncu)
bysort year country: egen double total_va_usd = total(va_usd)
gen double xr_lcu_per_usd = total_va_ncu/total_va_usd
assert xr_lcu_per_usd>0 & !missing(xr_lcu_per_usd)

* Diagnostic only. Do not stop the whole all-year exercise over tiny vintage/rounding
* inconsistencies; export them if they exceed 0.1 percent.
gen double xr_industry = va_ncu/va_usd if va_ncu>1e-10 & va_usd>1e-10
gen double xr_rel_error = abs(xr_industry/xr_lcu_per_usd-1) if !missing(xr_industry)

preserve
    keep if xr_rel_error>=0.001 & !missing(xr_rel_error)
    count
    if r(N)>0 {
        sort year country industry
        export delimited year country industry va_ncu va_usd ///
            xr_lcu_per_usd xr_industry xr_rel_error using ///
            "$OUT/exchange_rate_diagnostics.csv", replace
    }
restore

gen double labcomp_usd = lab_ncu/xr_lcu_per_usd
assert labcomp_usd>=0 & !missing(labcomp_usd)

preserve
    keep year country xr_lcu_per_usd total_va_ncu total_va_usd
    duplicates drop
    isid year country
    save "$WORK/wiod_implied_exchange_rates.dta", replace
    export delimited using "$OUT/wiod_implied_exchange_rates.csv", replace
restore

drop total_va_ncu total_va_usd xr_industry xr_rel_error
compress
save "$WORK/wiod_accounts_sea_40.dta", replace

********************************************************************************
* 8. PPP: WORLD BANK FOR 39 COUNTRIES + IMF TAIWAN, 1995-2011
********************************************************************************

import excel using "$PPP_XLSX", sheet("Data") allstring clear

rename A year_text
rename D country
rename E ppp_text

keep year_text country ppp_text
replace country = upper(strtrim(country))
destring year_text, gen(year) force
destring ppp_text, gen(ppp_lcu_per_intl) force

gen byte keep_country = strpos(" `CTRY40' ", " " + country + " ") > 0
keep if keep_country & inrange(year,`FIRSTYEAR',`LASTYEAR')
drop if country=="TWN"

keep year country ppp_lcu_per_intl
sort year country
isid year country
count
assert r(N)==663
assert ppp_lcu_per_intl>0 & !missing(ppp_lcu_per_intl)

tempfile wdi_ppp
save `wdi_ppp'

* Taiwan PPP from IMF WEO.
import delimited using "$IMF_CSV", clear varnames(1) case(lower) ///
    colrange(1:7) encoding("UTF-8")

keep if strtrim(country)=="Taiwan Province of China"
keep if strpos(strtrim(indicator), ///
    "Rate, Domestic currency per international dollar in PPP terms")==1

capture confirm numeric variable time_period
if _rc destring time_period, replace force
capture confirm numeric variable obs_value
if _rc destring obs_value, replace ignore(",") force

keep if frequency=="Annual" & inrange(time_period,`FIRSTYEAR',`LASTYEAR')
keep time_period obs_value
rename time_period year
rename obs_value ppp_lcu_per_intl
gen str3 country = "TWN"
order year country ppp_lcu_per_intl

assert ppp_lcu_per_intl>0 & !missing(ppp_lcu_per_intl)
isid year country
count
assert r(N)==17

sort year
export delimited using "$OUT/taiwan_ppp_from_imf_1995_2011.csv", replace

append using `wdi_ppp'
sort year country
isid year country
count
assert r(N)==680
assert ppp_lcu_per_intl>0 & !missing(ppp_lcu_per_intl)

label var ppp_lcu_per_intl "GDP PPP conversion factor, LCU per international $"
compress
save "$WORK/ppp_clean_40_1995_2011.dta", replace
export delimited using "$OUT/ppp_clean_40_1995_2011.csv", replace

********************************************************************************
* 9. RICCI COEFFICIENTS FOR THE 40 OBSERVED COUNTRIES
********************************************************************************

use "$WORK/wiod_accounts_sea_40.dta", clear
merge m:1 year country using "$WORK/ppp_clean_40_1995_2011.dta"
assert _merge==3
drop _merge

gen double ppp_usd_factor = xr_lcu_per_usd/ppp_lcu_per_intl
assert ppp_usd_factor>0 & !missing(ppp_usd_factor)

* Ricci equation (3): industry-normalized PPP.
gen double va_ppp_raw = va_usd*ppp_usd_factor
bysort year industry: egen double world_va_usd_j = total(va_usd)
bysort year industry: egen double world_va_ppp_raw_j = total(va_ppp_raw)

gen double industry_ppp_scale = world_va_usd_j/world_va_ppp_raw_j
gen double ppp_usd_factor_ind = ppp_usd_factor*industry_ppp_scale
gen double va_ppp = va_usd*ppp_usd_factor_ind
gen double erdi = 1/ppp_usd_factor_ind

bysort year industry: egen double world_va_ppp_j = total(va_ppp)
gen double ppp_norm_error = abs(world_va_ppp_j-world_va_usd_j) / ///
                            max(abs(world_va_usd_j),1e-12)
summarize ppp_norm_error, meanonly
assert r(max)<1e-8

bysort year industry: egen double world_hours_j = total(hours)
bysort year industry: egen double world_labcomp_j = total(labcomp_usd)
assert world_hours_j>0

gen double labor_h = (va_ppp/world_va_usd_j)*world_hours_j
assert !missing(labor_h)

* A valid Ricci source-sector coefficient requires positive gross output,
* positive value added, positive observed/imputed hours, and therefore positive
* homogeneous labor. Nonpositive-VA cells are retained in WIOD aggregates but
* are not assigned exporter transfer rates.
gen byte ricci_eligible = go_usd>1e-10 & va_usd>1e-10 & hours>0 & labor_h>0
assert labor_h>0 if go_usd>1e-10 & !ricci_ineligible

* Ricci equation (2): MEV.
bysort year: egen double world_va = total(va_usd)
bysort year: egen double world_hours = total(hours)
bysort year: egen double world_labcomp = total(labcomp_usd)
gen double mev = world_va/world_hours

gen double va_per_h_world_ind   = world_va_usd_j/world_hours_j
gen double wage_per_h_world_ind = world_labcomp_j/world_hours_j
gen double wage_per_h_world     = world_labcomp/world_hours
gen double wage_per_h_nat_ind   = labcomp_usd/labor_h if labor_h>0

* Equations (7)-(9).
gen double coef_inter = va_per_h_world_ind-mev
gen double coef_intra = (erdi-1)*va_per_h_world_ind
gen double coef_total = coef_inter+coef_intra

gen double coef_wage_inter = wage_per_h_world_ind-wage_per_h_world
gen double coef_wage_intra = wage_per_h_nat_ind-wage_per_h_world_ind if labor_h>0
gen double coef_wage = coef_wage_inter+coef_wage_intra if labor_h>0

gen double coef_profit_inter = coef_inter-coef_wage_inter
gen double coef_profit_intra = coef_intra-coef_wage_intra if labor_h>0
gen double coef_profit = coef_total-coef_wage if labor_h>0

* Direct equation-(6) identity check.
gen double coef_total_direct = va_usd/labor_h-mev if labor_h>0
gen double coef_error = abs(coef_total-coef_total_direct) / ///
                        max(abs(coef_total_direct),1e-12) if labor_h>0
summarize coef_error, meanonly
assert r(max)<1e-8

* Transfer rate per dollar of the exporter's gross output. Multiplying this by
* bilateral trade is algebraically identical to coef * (X/Q) * L_h.
gen double rate_total        = coef_total*labor_h/go_usd if ricci_eligible
gen double rate_inter        = coef_inter*labor_h/go_usd if ricci_eligible
gen double rate_intra        = coef_intra*labor_h/go_usd if ricci_eligible
gen double rate_wage         = coef_wage*labor_h/go_usd if ricci_eligible
gen double rate_wage_inter   = coef_wage_inter*labor_h/go_usd if ricci_eligible
gen double rate_wage_intra   = coef_wage_intra*labor_h/go_usd if ricci_eligible
gen double rate_profit       = coef_profit*labor_h/go_usd if ricci_eligible
gen double rate_profit_inter = coef_profit_inter*labor_h/go_usd if ricci_eligible
gen double rate_profit_intra = coef_profit_intra*labor_h/go_usd if ricci_eligible

keep year country industry va_usd go_usd hours labcomp_usd ///
     hours_imputed lab_imputed ricci_ineligible ricci_eligible ///
     xr_lcu_per_usd ppp_lcu_per_intl ppp_usd_factor ///
     ppp_usd_factor_ind erdi labor_h mev ///
     coef_total coef_inter coef_intra ///
     coef_wage coef_wage_inter coef_wage_intra ///
     coef_profit coef_profit_inter coef_profit_intra ///
     rate_total rate_inter rate_intra ///
     rate_wage rate_wage_inter rate_wage_intra ///
     rate_profit rate_profit_inter rate_profit_intra

isid year country industry
compress
save "$WORK/ricci_coefficients_40_all_years.dta", replace

********************************************************************************
* 10. RICCI-COMPARABLE 40-COUNTRY BILATERAL TRANSFERS, 1995-2011
********************************************************************************

use "$WORK/ricci_coefficients_40_all_years.dta", clear
rename country exporter
keep year exporter industry go_usd va_usd hours ricci_ineligible ricci_eligible rate_*
save "$WORK/transfer_rates_40.dta", replace

use "$WORK/bilateral_trade_40.dta", clear
merge m:1 year exporter industry using "$WORK/transfer_rates_40.dta"

* Save positive trade flows whose exporter-sector cannot receive a valid Ricci
* rate (zero GO, nonpositive VA, or zero labor). They are omitted rather than
* silently imputed. In this release the nonpositive-VA/zero-labor cases are very
* sparse and are separately listed in ricci_ineligible_source_cells.csv.
preserve
    keep if _merge==3 & missing(rate_total) & abs(trade_usd)>1e-8
    count
    if r(N)>0 {
        sort year exporter industry importer
        export delimited using "$OUT/ineligible_source_positive_exports_40.csv", replace
    }
restore

keep if _merge==3
drop _merge
drop if missing(rate_total)

gen double tr_total        = rate_total*trade_usd
gen double tr_inter        = rate_inter*trade_usd
gen double tr_intra        = rate_intra*trade_usd
gen double tr_wage         = rate_wage*trade_usd
gen double tr_wage_inter   = rate_wage_inter*trade_usd
gen double tr_wage_intra   = rate_wage_intra*trade_usd
gen double tr_profit       = rate_profit*trade_usd
gen double tr_profit_inter = rate_profit_inter*trade_usd
gen double tr_profit_intra = rate_profit_intra*trade_usd

gen double check_source = abs(tr_total-tr_inter-tr_intra)
gen double check_income = abs(tr_total-tr_wage-tr_profit)
summarize check_source, meanonly
assert r(max)<1e-7
summarize check_income, meanonly
assert r(max)<1e-7
drop check_source check_income

keep year exporter importer industry trade_usd tr_*
compress
save "$WORK/bilateral_value_transfers_40.dta", replace

********************************************************************************
* 11. COUNTRY NET TRANSFERS: 40-COUNTRY RICCI-COMPARABLE SYSTEM
********************************************************************************

use "$WORK/bilateral_value_transfers_40.dta", clear

preserve
    keep year exporter tr_*
    rename exporter country
    collapse (sum) tr_*, by(year country)
    tempfile export_credits40
    save `export_credits40'
restore

keep year importer tr_*
rename importer country
foreach v of varlist tr_* {
    replace `v' = -`v'
}
collapse (sum) tr_*, by(year country)
append using `export_credits40'
collapse (sum) tr_*, by(year country)

rename tr_total        net_transfer
rename tr_inter        net_inter
rename tr_intra        net_intra
rename tr_wage         net_wage
rename tr_wage_inter   net_wage_inter
rename tr_wage_intra   net_wage_intra
rename tr_profit       net_profit
rename tr_profit_inter net_profit_inter
rename tr_profit_intra net_profit_intra

preserve
    use "$WORK/wiod_accounts_40.dta", clear
    collapse (sum) va_usd, by(year country)
    tempfile country_va40
    save `country_va40'
restore

merge 1:1 year country using `country_va40'
assert _merge==3
drop _merge

gen double net_transfer_pct_va = 100*net_transfer/va_usd

* Zero-sum check across the 40-country system.
bysort year: egen double net_sum_check = total(net_transfer)
assert abs(net_sum_check)<1e-5
drop net_sum_check

sort country year
compress
save "$OUT/table3_40countries_all_years_long.dta", replace

********************************************************************************
* 12. EXPORT TABLE 3, 40 COUNTRIES
********************************************************************************

preserve
    keep country year net_transfer net_transfer_pct_va
    sort country year
    export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
        sheet("T3_40_long") firstrow(variables) replace

    keep country year net_transfer
    reshape wide net_transfer, i(country) j(year)
    sort country
    export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
        sheet("T3_40_wide") firstrow(variables) sheetmodify
restore

********************************************************************************
* 13. TABLE 2, 12 RICCI REGIONS, ALL YEARS
********************************************************************************

gen str20 region = ""
replace region = "North America"   if inlist(country,"CAN","USA")
replace region = "North EMU"       if inlist(country,"AUT","BEL","DEU","FIN","FRA","LUX","NLD")
replace region = "South EMU"       if inlist(country,"CYP","ESP","GRC","IRL","ITA","MLT","PRT")
replace region = "North Europe"    if inlist(country,"DNK","GBR","SWE")
replace region = "East Europe"     if inlist(country,"BGR","CZE","EST","HUN","LTU")
replace region = "East Europe"     if inlist(country,"LVA","POL","ROU","SVK","SVN")
replace region = "Latin America"   if inlist(country,"BRA","MEX")
replace region = "China"           if country=="CHN"
replace region = "India"           if country=="IND"
replace region = "North East Asia" if inlist(country,"JPN","KOR")
replace region = "Other Asia"      if inlist(country,"IDN","TUR","TWN")
replace region = "Russia"          if country=="RUS"
replace region = "Australia"       if country=="AUS"
assert region!=""

collapse (sum) net_transfer net_inter net_intra ///
               net_wage net_wage_inter net_wage_intra ///
               net_profit net_profit_inter net_profit_intra va_usd, ///
    by(year region)

gen double net_transfer_pct_va = 100*net_transfer/va_usd
sort region year
compress
save "$OUT/table2_12regions_all_years_long.dta", replace

export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
    sheet("T2_12reg_long") firstrow(variables) sheetmodify

preserve
    keep region year net_transfer net_transfer_pct_va
    reshape wide net_transfer net_transfer_pct_va, i(region) j(year)
    sort region
    export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
        sheet("T2_12reg_wide") firstrow(variables) sheetmodify
restore

********************************************************************************
* 14. TABLE 2 BOTTOM PANEL, 12-REGION RICCI-COMPARABLE SYSTEM
********************************************************************************

bysort year: egen double global_va = total(va_usd)
keep if net_transfer>0

collapse (sum) total_net_transfer=net_transfer ///
               total_inter=net_inter ///
               total_intra=net_intra ///
               total_wage=net_wage ///
               total_profit=net_profit ///
         (max) global_va, by(year)

gen double total_net_pct_global_va = 100*total_net_transfer/global_va
gen double industrial_specialization_pct = 100*total_inter/total_net_transfer
gen double absolute_rent_pct = 100*total_intra/total_net_transfer
gen double wage_differential_pct = 100*total_wage/total_net_transfer
gen double profit_differential_pct = 100*total_profit/total_net_transfer

gen byte sea_partly_imputed = inrange(year,2010,2011)

order year sea_partly_imputed total_net_transfer total_net_pct_global_va ///
      industrial_specialization_pct absolute_rent_pct ///
      wage_differential_pct profit_differential_pct

sort year
save "$OUT/table2_12regions_decomposition_all_years.dta", replace
export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
    sheet("T2_12reg_decomp") firstrow(variables) sheetmodify

********************************************************************************
* 15. ESTIMATE ROW TRANSFER RATES FROM OBSERVED 40 COUNTRIES
********************************************************************************

* ROW lacks SEA hours and labor compensation. We estimate its transfer-per-dollar
* rate separately by year and industry as the gross-output-weighted mean of the
* 40 observed countries' rates. This is a neutral WIOD-only estimate.
use "$WORK/ricci_coefficients_40_all_years.dta", clear
keep if go_usd>1e-10 & !missing(rate_total)

local RATES "rate_total rate_inter rate_intra rate_wage rate_wage_inter rate_wage_intra rate_profit rate_profit_inter rate_profit_intra"

foreach v of local RATES {
    gen double num_`v' = `v'*go_usd
}

collapse (sum) num_rate_total num_rate_inter num_rate_intra ///
               num_rate_wage num_rate_wage_inter num_rate_wage_intra ///
               num_rate_profit num_rate_profit_inter num_rate_profit_intra ///
               observed_output=go_usd, by(year industry)

foreach v of local RATES {
    gen double `v' = num_`v'/observed_output
}

gen str3 exporter = "ROW"
gen byte row_rate_estimated = 1

keep year exporter industry observed_output row_rate_estimated `RATES'
isid year exporter industry

save "$OUT/row_transfer_rates_estimated.dta", replace
export delimited using "$OUT/row_transfer_rates_estimated.csv", replace

* Append estimated ROW rates to actual observed-country rates.
tempfile rowrates
save `rowrates'

use "$WORK/transfer_rates_40.dta", clear
gen byte row_rate_estimated = 0
gen double observed_output = .
append using `rowrates'
isid year exporter industry
save "$WORK/transfer_rates_41_with_ROW_est.dta", replace

********************************************************************************
* 16. EXTENDED BILATERAL TRANSFERS: 40 COUNTRIES + ESTIMATED ROW
********************************************************************************

use "$WORK/bilateral_trade_41.dta", clear
merge m:1 year exporter industry using "$WORK/transfer_rates_41_with_ROW_est.dta"

* Save unmatched/zero-output positive flows. ROW itself should match all industries;
* unmatched observed-country flows are the same zero-output anomalies as before.
preserve
    keep if (_merge!=3 | missing(rate_total)) & abs(trade_usd)>1e-8
    count
    if r(N)>0 {
        sort year exporter industry importer
        export delimited using "$OUT/unusable_positive_trade_41.csv", replace
    }
restore

keep if _merge==3
drop _merge
drop if missing(rate_total)

gen double tr_total        = rate_total*trade_usd
gen double tr_inter        = rate_inter*trade_usd
gen double tr_intra        = rate_intra*trade_usd
gen double tr_wage         = rate_wage*trade_usd
gen double tr_wage_inter   = rate_wage_inter*trade_usd
gen double tr_wage_intra   = rate_wage_intra*trade_usd
gen double tr_profit       = rate_profit*trade_usd
gen double tr_profit_inter = rate_profit_inter*trade_usd
gen double tr_profit_intra = rate_profit_intra*trade_usd

gen double check_source = abs(tr_total-tr_inter-tr_intra)
gen double check_income = abs(tr_total-tr_wage-tr_profit)
summarize check_source, meanonly
assert r(max)<1e-7
summarize check_income, meanonly
assert r(max)<1e-7
drop check_source check_income

keep year exporter importer industry trade_usd row_rate_estimated tr_*
compress
save "$WORK/bilateral_value_transfers_41_ROW_est.dta", replace

********************************************************************************
* 17. COUNTRY/ROW NET TRANSFERS IN THE EXTENDED 41-REGION SYSTEM
********************************************************************************

use "$WORK/bilateral_value_transfers_41_ROW_est.dta", clear

preserve
    keep year exporter tr_*
    rename exporter country
    collapse (sum) tr_*, by(year country)
    tempfile export_credits41
    save `export_credits41'
restore

keep year importer tr_*
rename importer country
foreach v of varlist tr_* {
    replace `v' = -`v'
}
collapse (sum) tr_*, by(year country)
append using `export_credits41'
collapse (sum) tr_*, by(year country)

rename tr_total        net_transfer
rename tr_inter        net_inter
rename tr_intra        net_intra
rename tr_wage         net_wage
rename tr_wage_inter   net_wage_inter
rename tr_wage_intra   net_wage_intra
rename tr_profit       net_profit
rename tr_profit_inter net_profit_inter
rename tr_profit_intra net_profit_intra

preserve
    use "$WORK/wiod_accounts_41.dta", clear
    collapse (sum) va_usd, by(year country)
    tempfile country_va41
    save `country_va41'
restore

merge 1:1 year country using `country_va41'
assert _merge==3
drop _merge

gen double net_transfer_pct_va = 100*net_transfer/va_usd
gen byte row_estimated = country=="ROW"
gen byte sea_partly_imputed = inrange(year,2010,2011)

* Extended system must also be zero-sum because every bilateral transfer is credited
* to the exporter and debited from the importer.
bysort year: egen double net_sum_check = total(net_transfer)
assert abs(net_sum_check)<1e-5
drop net_sum_check

sort country year
compress
save "$OUT/table3_41regions_ROW_est_all_years_long.dta", replace

********************************************************************************
* 18. EXPORT EXTENDED TABLE 3 INCLUDING ROW
********************************************************************************

preserve
    keep country year row_estimated sea_partly_imputed ///
         net_transfer net_transfer_pct_va
    sort country year
    export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
        sheet("T3_41_ROW_long") firstrow(variables) sheetmodify

    keep country year net_transfer
    reshape wide net_transfer, i(country) j(year)
    sort country
    export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
        sheet("T3_41_ROW_wide") firstrow(variables) sheetmodify
restore

********************************************************************************
* 19. EXTENDED TABLE 2: 12 RICCI REGIONS + ROW
********************************************************************************

gen str24 region = ""
replace region = "North America"   if inlist(country,"CAN","USA")
replace region = "North EMU"       if inlist(country,"AUT","BEL","DEU","FIN","FRA","LUX","NLD")
replace region = "South EMU"       if inlist(country,"CYP","ESP","GRC","IRL","ITA","MLT","PRT")
replace region = "North Europe"    if inlist(country,"DNK","GBR","SWE")
replace region = "East Europe"     if inlist(country,"BGR","CZE","EST","HUN","LTU")
replace region = "East Europe"     if inlist(country,"LVA","POL","ROU","SVK","SVN")
replace region = "Latin America"   if inlist(country,"BRA","MEX")
replace region = "China"           if country=="CHN"
replace region = "India"           if country=="IND"
replace region = "North East Asia" if inlist(country,"JPN","KOR")
replace region = "Other Asia"      if inlist(country,"IDN","TUR","TWN")
replace region = "Russia"          if country=="RUS"
replace region = "Australia"       if country=="AUS"
replace region = "Rest of World (est.)" if country=="ROW"
assert region!=""

collapse (sum) net_transfer net_inter net_intra ///
               net_wage net_wage_inter net_wage_intra ///
               net_profit net_profit_inter net_profit_intra va_usd, ///
    by(year region)

gen double net_transfer_pct_va = 100*net_transfer/va_usd
gen byte row_estimated = region=="Rest of World (est.)"
gen byte sea_partly_imputed = inrange(year,2010,2011)

sort region year
compress
save "$OUT/table2_13regions_ROW_est_all_years_long.dta", replace

export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
    sheet("T2_13_ROW_long") firstrow(variables) sheetmodify

preserve
    keep region year net_transfer net_transfer_pct_va
    reshape wide net_transfer net_transfer_pct_va, i(region) j(year)
    sort region
    export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
        sheet("T2_13_ROW_wide") firstrow(variables) sheetmodify
restore

********************************************************************************
* 20. EXTENDED TABLE 2 BOTTOM PANEL INCLUDING ESTIMATED ROW
********************************************************************************

bysort year: egen double global_va = total(va_usd)
keep if net_transfer>0

collapse (sum) total_net_transfer=net_transfer ///
               total_inter=net_inter ///
               total_intra=net_intra ///
               total_wage=net_wage ///
               total_profit=net_profit ///
         (max) global_va, by(year)

gen double total_net_pct_global_va = 100*total_net_transfer/global_va
gen double industrial_specialization_pct = 100*total_inter/total_net_transfer
gen double absolute_rent_pct = 100*total_intra/total_net_transfer
gen double wage_differential_pct = 100*total_wage/total_net_transfer
gen double profit_differential_pct = 100*total_profit/total_net_transfer
gen byte row_included_estimate = 1
gen byte sea_partly_imputed = inrange(year,2010,2011)

order year row_included_estimate sea_partly_imputed ///
      total_net_transfer total_net_pct_global_va ///
      industrial_specialization_pct absolute_rent_pct ///
      wage_differential_pct profit_differential_pct

sort year
save "$OUT/table2_13regions_ROW_est_decomposition.dta", replace
export excel using "$OUT/ricci_all_years_1995_2011.xlsx", ///
    sheet("T2_13_ROW_decomp") firstrow(variables) sheetmodify

********************************************************************************
* 21. FINAL MESSAGE
********************************************************************************

di as result "All-years extension completed: 1995-2011."
di as result "1995-2009: observed SEA labor inputs (apart from structural-zero cells)."
di as result "2010-2011: missing SEA hours/LAB transparently imputed from 2009 benchmarks."
di as result "ROW: estimated using output-weighted industry transfer rates of the 40 observed countries."
di as result "Workbook: $OUT/ricci_all_years_1995_2011.xlsx"
di as result "40-country results remain separate from the ROW-inclusive extension."


/****************************************************************************************
 PART 2. WIOD NOVEMBER 2016 RELEASE: RAW WIOD/SEA/PPP/PWT -> RICCI TABLE 2 DATA
 Source: ricci_wiod2016_2000_2014_all43.do
****************************************************************************************/

/****************************************************************************************
 Ricci-style unequal-exchange extension using the WIOD November 2016 release, 2000-2014

 PURPOSE
   Reproduce Ricci-style Table 2 / Table 3 accounting for every year 2000-2014
   using the October 2016 WIOT files and November 2016 SEA.

 COUNTRY COVERAGE
   The 2016 WIOD release separately identifies 43 countries plus ROW.
   This version KEEPS Switzerland (CHE), Croatia (HRV), and Norway (NOR) as
   separately observed countries. They are NOT folded into ROW.

   For regional reporting, CHE + HRV + NOR are grouped as:
       "Additional Core"
   This leaves the original 12 Ricci regions unchanged for clean comparison,
   while keeping the three extra economies inside the observed core bloc.

   Therefore this file produces:
     A. 43 observed-country results (ROW excluded)
     B. 43 observed countries + WIOD ROW, with ROW's transfer coefficients estimated
     C. 13-region Table 2 = original 12 Ricci regions + Additional Core
     D. 14-region Table 2 = those 13 + Rest of World (est.)

 INDUSTRIES / LABOR
   - WIOD 2016 has 56 industries. We use industries 1-54 and exclude:
       55 = household services
       56 = extraterritorial organizations
   - The 2016 SEA reports H_EMPE (hours worked by employees), not H_EMP.
     Total hours for persons engaged are reconstructed as:
       H = H_EMPE * EMP / EMPE
     where possible.
   - China has no usable H_EMPE/EMPE series in this SEA release. For those cells
     only, PWT 11.0 average annual hours (avh) are multiplied by SEA industry EMP.
     Every fallback is flagged and exported.

 PPP
   - World Bank WDI GDP PPP conversion factor for all observed countries except Taiwan.
   - IMF WEO PPP conversion factor for Taiwan.

 REQUIRED FILES in $RAW
   WIOT2000_October16_ROW.dta ... WIOT2014_October16_ROW.dta
   (the code also recognizes downloaded suffixes such as "(1)" or "(2)")
   Socio_Economic_Accounts.xlsx
   P_Data_Extract_From_World_Development_Indicators.xlsx
   dataset_2026-08-04T19_36_35.328252551Z_DEFAULT_INTEGRATION_IMF.RES_WEO_9.0.0.csv
   pwt110.dta

 MAIN OUTPUTS
   output_2016release_2000_2014/ricci_2000_2014_wiod2016_all43.xlsx
   output_2016release_2000_2014/table2_13regions_all43_2000_2014_long.dta
   output_2016release_2000_2014/table2_14regions_ROW_est_2000_2014_long.dta
   output_2016release_2000_2014/table3_43countries_2000_2014_long.dta
   output_2016release_2000_2014/table3_44regions_ROW_est_2000_2014_long.dta
****************************************************************************************/

version 17.0
clear all
set more off
set varabbrev off
set maxvar 32767

********************************************************************************
* 0. PATHS
********************************************************************************

global ROOT "`c(pwd)'"
global RAW  "$ROOT/raw"
global WORK "$ROOT/work_2016release_2000_2014"
global OUT  "$ROOT/output_2016release_2000_2014"

capture mkdir "$WORK"
capture mkdir "$OUT"

global SEA_XLSX "$RAW/Socio_Economic_Accounts.xlsx"
global PPP_XLSX "$RAW/P_Data_Extract_From_World_Development_Indicators.xlsx"
global IMF_CSV  "$RAW/dataset_2026-08-04T19_36_35.328252551Z_DEFAULT_INTEGRATION_IMF.RES_WEO_9.0.0.csv"
global PWT      "$RAW/pwt110.dta"

* Resolve WIOT filenames once. The code accepts clean names and common
* duplicate-download suffixes "(1)" and "(2)".
tempname wiotmanifest
postfile `wiotmanifest' int year str244 wiotfile using "$WORK/wiot_file_manifest.dta", replace

foreach y of numlist 2000/2014 {
    local wiotfile "$RAW/WIOT`y'_October16_ROW.dta"
    capture confirm file "`wiotfile'"
    if _rc {
        local wiotfile "$RAW/WIOT`y'_October16_ROW(1).dta"
        capture confirm file "`wiotfile'"
    }
    if _rc {
        local wiotfile "$RAW/WIOT`y'_October16_ROW(2).dta"
        capture confirm file "`wiotfile'"
    }
    if _rc {
        di as error "Required WIOT file not found for year `y'."
        di as error "Expected WIOT`y'_October16_ROW.dta, ...ROW(1).dta, or ...ROW(2).dta in $RAW"
        exit 601
    }
    post `wiotmanifest' (`y') ("`wiotfile'")
}
postclose `wiotmanifest'

foreach f in "$SEA_XLSX" "$PPP_XLSX" "$IMF_CSV" "$PWT" {
    capture confirm file `f'
    if _rc {
        di as error "Required input file not found: `f'"
        exit 601
    }
}

local FIRSTYEAR 2000
local LASTYEAR  2014

* WIOD 2016's 43 separately identified countries (all kept observed).
local CTRY43 "AUS AUT BEL BGR BRA CAN CHE CHN CYP CZE DEU DNK ESP EST FIN FRA GBR GRC HRV HUN IDN IND IRL ITA JPN KOR LTU LUX LVA MEX MLT NLD NOR POL PRT ROU RUS SVK SVN SWE TUR TWN USA"

* WIOD 2016 separately identifies these 43 countries plus current ROW.
local CTRY44 "AUS AUT BEL BGR BRA CAN CHE CHN CYP CZE DEU DNK ESP EST FIN FRA GBR GRC HRV HUN IDN IND IRL ITA JPN KOR LTU LUX LVA MEX MLT NLD NOR POL PRT ROU RUS SVK SVN SWE TUR TWN USA ROW"

* Stubs used by the wide WIOT files: vAUS1-vAUS61, ..., vROW1-vROW61.
local STUBS ""
foreach c of local CTRY44 {
    local STUBS "`STUBS' v`c'"
}

capture program drop resolve_wiot
program define resolve_wiot, rclass
    syntax , YEAR(integer)
    local f "$RAW/WIOT`year'_October16_ROW.dta"
    capture confirm file "`f'"
    if _rc {
        local f "$RAW/WIOT`year'_October16_ROW(1).dta"
        capture confirm file "`f'"
    }
    if _rc {
        local f "$RAW/WIOT`year'_October16_ROW(2).dta"
        capture confirm file "`f'"
    }
    if _rc {
        di as error "Could not resolve WIOT file for `year'."
        exit 601
    }
    return local file "`f'"
end

********************************************************************************
* 1. INDUSTRY MAP FROM THE 2016 WIOT
********************************************************************************

resolve_wiot, year(2014)
local wiot2014 "`r(file)'"
use IndustryCode IndustryDescription Country RNr Year using "`wiot2014'", clear
keep if Country=="AUS" & inrange(RNr,1,54)
keep RNr IndustryCode IndustryDescription
rename RNr industry
rename IndustryCode sea_code
rename IndustryDescription industry_name
isid industry
isid sea_code
sort industry
save "$WORK/wiod2016_industry_map.dta", replace
export delimited using "$OUT/wiod2016_industry_map.csv", replace

********************************************************************************
* 2. VALUE ADDED AND GROSS OUTPUT: 43 OBSERVED COUNTRIES + CURRENT ROW
********************************************************************************

local first_accounts 1
tempfile accounts_stack

foreach y of numlist 2000/2014 {
    resolve_wiot, year(`y')
    local wiotfile "`r(file)'"
    use IndustryCode Country Year vAUS1-vROW61 using "`wiotfile'", clear
    keep if Country=="TOT" & inlist(IndustryCode,"VA","GO")

    * First reshape columns 1-61 within each destination-country stub.
    reshape long `STUBS', i(Year IndustryCode) j(industry)
    keep if inrange(industry,1,54)

    * Then stack destination countries.
    reshape long v, i(Year IndustryCode industry) j(country) string
    rename v value
    rename Year year

    reshape wide value, i(year country industry) j(IndustryCode) string
    rename valueVA va_usd
    rename valueGO go_usd

    isid year country industry
    if `first_accounts' {
        save `accounts_stack', replace
        local first_accounts 0
    }
    else {
        append using `accounts_stack'
        save `accounts_stack', replace
    }
}

use `accounts_stack', clear
sort year country industry
isid year country industry
assert !missing(va_usd,go_usd)
label var va_usd "Value added, million current US$"
label var go_usd "Gross output, million current US$"
save "$WORK/wiod2016_accounts_44.dta", replace

preserve
    keep if country!="ROW"
    gen byte keep43 = strpos(" `CTRY43' ", " " + country + " ")>0
    keep if keep43
    drop keep43
    isid year country industry
    save "$WORK/wiod2016_accounts_43.dta", replace
restore

********************************************************************************
* 3. BILATERAL TRADE: 43 OBSERVED COUNTRIES + CURRENT ROW
********************************************************************************

local first_trade 1
tempfile trade_stack

foreach y of numlist 2000/2014 {
    resolve_wiot, year(`y')
    local wiotfile "`r(file)'"
    use Country RNr Year vAUS1-vROW61 using "`wiotfile'", clear

    gen byte keep_source = strpos(" `CTRY44' ", " " + Country + " ")>0
    keep if keep_source & inrange(RNr,1,54)
    drop keep_source

    * Stack the 61 intermediate/final-use columns within each destination country,
    * then sum them to bilateral exporter-industry -> importer trade.
    reshape long `STUBS', i(Country RNr Year) j(usecol)
    collapse (sum) `STUBS', by(Country RNr Year)
    reshape long v, i(Country RNr Year) j(importer) string

    rename Country exporter
    rename RNr industry
    rename Year year
    rename v trade_usd

    collapse (sum) trade_usd, by(year exporter importer industry)
    drop if exporter==importer
    drop if abs(trade_usd)<1e-12

    isid year exporter importer industry
    if `first_trade' {
        save `trade_stack', replace
        local first_trade 0
    }
    else {
        append using `trade_stack'
        save `trade_stack', replace
    }
}

use `trade_stack', clear
sort year exporter importer industry
isid year exporter importer industry
label var trade_usd "Bilateral trade, million current US$"
save "$WORK/bilateral_trade_44.dta", replace

preserve
    keep if exporter!="ROW" & importer!="ROW"
    gen byte keep_exporter = strpos(" `CTRY43' ", " " + exporter + " ")>0
    gen byte keep_importer = strpos(" `CTRY43' ", " " + importer + " ")>0
    keep if keep_exporter & keep_importer
    drop keep_exporter keep_importer
    isid year exporter importer industry
    save "$WORK/bilateral_trade_43.dta", replace
restore

********************************************************************************
* 4. IMPORT 2016 SEA AND RECONSTRUCT TOTAL HOURS WORKED
********************************************************************************

import excel using "$SEA_XLSX", sheet("DATA") allstring clear

rename (A B C D E F G H I J K L M N O P Q R S) ///
       (country variable description sea_code ///
        sea2000 sea2001 sea2002 sea2003 sea2004 sea2005 sea2006 sea2007 ///
        sea2008 sea2009 sea2010 sea2011 sea2012 sea2013 sea2014)

replace country  = upper(strtrim(country))
replace variable = upper(strtrim(variable))
replace sea_code = strtrim(sea_code)

keep if inlist(variable,"VA","LAB","COMP","EMP","EMPE","H_EMPE")
gen byte keep_country = strpos(" `CTRY43' ", " " + country + " ")>0
keep if keep_country
drop keep_country

foreach y of numlist 2000/2014 {
    destring sea`y', replace force
}

merge m:1 sea_code using "$WORK/wiod2016_industry_map.dta", keepusing(industry)
keep if _merge==3
drop _merge
keep if inrange(industry,1,54)

reshape long sea, i(country variable industry) j(year)
rename sea sea_value
keep if inrange(year,`FIRSTYEAR',`LASTYEAR')

isid year country industry variable
reshape wide sea_value, i(year country industry) j(variable) string

rename sea_valueVA     va_ncu
rename sea_valueLAB    lab_raw_ncu
rename sea_valueCOMP   comp_ncu
rename sea_valueEMP    employment
rename sea_valueEMPE   employees
rename sea_valueH_EMPE h_empe

label var va_ncu     "SEA value added, million national currency"
label var lab_raw_ncu "SEA labour compensation (LAB), million national currency"
label var comp_ncu   "SEA employee compensation (COMP), million national currency"
label var employment "SEA persons engaged, thousands"
label var employees  "SEA employees, thousands"
label var h_empe     "SEA hours worked by employees, millions"

* Reconstruct hours worked by all persons engaged using average employee hours.
* H_EMPE is in millions of hours and EMPE/EMP are in thousands, so
* H_EMPE/EMPE is numerically thousands of hours per employee.
gen double hours = h_empe*(employment/employees) if h_empe>0 & employees>0 & employment>=0

* First fallback: country-year SEA average hours per employee.
gen byte hours_countryavg = 0
bysort year country: egen double cy_h_empe = total(h_empe)
bysort year country: egen double cy_employees = total(employees)
gen double cy_hours_per_employee = cy_h_empe/cy_employees if cy_employees>0 & cy_h_empe>0
replace hours_countryavg = 1 if missing(hours) & employment>0 & cy_hours_per_employee>0
replace hours = employment*cy_hours_per_employee if hours_countryavg

* The 2016 SEA has no H_EMPE/EMPE series for China. Use PWT 11.0 country-year
* average annual hours only where SEA cannot provide any hours ratio.
preserve
    use "$PWT", clear
    keep if inrange(year,`FIRSTYEAR',`LASTYEAR')
    keep countrycode year avh
    rename countrycode country
    keep if strpos(" `CTRY43' ", " " + country + " ")>0
    keep if !missing(avh)
    isid year country
    tempfile pwt_avh
    save `pwt_avh'
restore

merge m:1 year country using `pwt_avh', keep(master match) nogen
gen byte hours_pwt = 0
replace hours_pwt = 1 if missing(hours) & employment>0 & avh>0
* EMP is thousands and avh is hours/person/year; divide by 1000 -> million hours.
replace hours = employment*(avh/1000) if hours_pwt
replace hours = 0 if employment==0 & missing(hours)

* LAB is preferred. If LAB is missing but employee compensation is observed,
* scale COMP by EMP/EMPE as a transparent fallback for self-employment labor.
gen double lab_ncu = lab_raw_ncu
gen byte lab_from_comp = 0
replace lab_from_comp = 1 if missing(lab_ncu) & comp_ncu>=0 & employees>0 & employment>=0
replace lab_ncu = comp_ncu*(employment/employees) if lab_from_comp
replace lab_ncu = 0 if employment==0 & missing(lab_ncu)

keep year country industry va_ncu lab_ncu lab_raw_ncu comp_ncu ///
     employment employees h_empe hours avh hours_countryavg hours_pwt lab_from_comp
isid year country industry
compress
save "$WORK/sea2016_2000_2014_clean.dta", replace

********************************************************************************
* 5. MERGE SEA TO WIOT ACCOUNTS; DIAGNOSTICS; CURRENT EXCHANGE RATES
********************************************************************************

use "$WORK/wiod2016_accounts_43.dta", clear
merge 1:1 year country industry using "$WORK/sea2016_2000_2014_clean.dta"

* There should be a matched SEA row for every 43 observed-country-industry-year.
count if _merge!=3
if r(N)>0 {
    preserve
        keep if _merge!=3
        export delimited using "$OUT/sea_wiot_unmatched_cells.csv", replace
    restore
    di as error "Unmatched WIOT/SEA cells exist. See sea_wiot_unmatched_cells.csv"
    exit 459
}
drop _merge

* Structural zeros only.
foreach v in va_ncu lab_ncu hours employment employees h_empe {
    replace `v'=0 if missing(`v') & abs(go_usd)<1e-10 & abs(va_usd)<1e-10
}

* Active cells must have the inputs needed for world labor/value calculations.
foreach v in va_ncu hours lab_ncu {
    count if go_usd>1e-10 & missing(`v')
    if r(N)>0 {
        preserve
            keep if go_usd>1e-10 & missing(`v')
            export delimited using "$OUT/missing_active_`v'_cells.csv", replace
        restore
        di as error "Active cells with missing `v'. See missing_active_`v'_cells.csv"
        exit 459
    }
}

* Flag nonpositive-VA/zero-hours source sectors instead of forcing values.
gen byte ricci_ineligible = go_usd>1e-10 & (va_usd<=1e-10 | hours<=0)
preserve
    keep if ricci_ineligible
    sort year country industry
    export delimited using "$OUT/ricci_ineligible_source_cells_2000_2014.csv", replace
restore

* Save the reconstruction/fallback flags.
preserve
    keep if hours_countryavg | hours_pwt | lab_from_comp
    sort year country industry
    export delimited using "$OUT/sea_hours_lab_fallbacks_2000_2014.csv", replace
restore

* Current market exchange rate: LCU per current US$, inferred within the SAME
* 2016 release from aggregate SEA VA / WIOT VA.
bysort year country: egen double total_va_ncu = total(va_ncu)
bysort year country: egen double total_va_usd = total(va_usd)
gen double xr_lcu_per_usd = total_va_ncu/total_va_usd
assert xr_lcu_per_usd>0 & !missing(xr_lcu_per_usd)

gen double labcomp_usd = lab_ncu/xr_lcu_per_usd
assert labcomp_usd>=0 & !missing(labcomp_usd)

preserve
    keep year country xr_lcu_per_usd total_va_ncu total_va_usd
    duplicates drop
    isid year country
    save "$WORK/wiod2016_implied_exchange_rates.dta", replace
    export delimited using "$OUT/wiod2016_implied_exchange_rates.csv", replace
restore

drop total_va_ncu total_va_usd
compress
save "$WORK/wiod2016_accounts_sea_43.dta", replace

********************************************************************************
* 6. PPP: WORLD BANK FOR 42 COUNTRIES + IMF TAIWAN, 2000-2014
********************************************************************************

import excel using "$PPP_XLSX", sheet("Data") allstring clear
rename A year_text
rename D country
rename E ppp_text

keep year_text country ppp_text
replace country=upper(strtrim(country))
destring year_text, gen(year) force
destring ppp_text, gen(ppp_lcu_per_intl) force

gen byte keep_country = strpos(" `CTRY43' ", " " + country + " ")>0
keep if keep_country & inrange(year,`FIRSTYEAR',`LASTYEAR')
drop if country=="TWN"

keep year country ppp_lcu_per_intl
sort year country
isid year country
count
assert r(N)==630
assert ppp_lcu_per_intl>0 & !missing(ppp_lcu_per_intl)

tempfile wdi_ppp
save `wdi_ppp'

* Taiwan: IMF WEO PPP conversion factor.
import delimited using "$IMF_CSV", clear varnames(1) case(lower) ///
    colrange(1:7) encoding("UTF-8")

keep if strtrim(country)=="Taiwan Province of China"
keep if strpos(strtrim(indicator), ///
    "Rate, Domestic currency per international dollar in PPP terms")==1

capture confirm numeric variable time_period
if _rc destring time_period, replace force
capture confirm numeric variable obs_value
if _rc destring obs_value, replace ignore(",") force

keep if frequency=="Annual" & inrange(time_period,`FIRSTYEAR',`LASTYEAR')
keep time_period obs_value
rename time_period year
rename obs_value ppp_lcu_per_intl
gen str3 country="TWN"
order year country ppp_lcu_per_intl

assert ppp_lcu_per_intl>0 & !missing(ppp_lcu_per_intl)
isid year country
count
assert r(N)==15

sort year
export delimited using "$OUT/taiwan_ppp_from_imf_2000_2014.csv", replace

append using `wdi_ppp'
sort year country
isid year country
count
assert r(N)==645
assert ppp_lcu_per_intl>0 & !missing(ppp_lcu_per_intl)

save "$WORK/ppp_clean_43_2000_2014.dta", replace
export delimited using "$OUT/ppp_clean_43_2000_2014.csv", replace

********************************************************************************
* 7. RICCI COEFFICIENTS FOR THE 43 OBSERVED COUNTRIES
********************************************************************************

use "$WORK/wiod2016_accounts_sea_43.dta", clear
merge m:1 year country using "$WORK/ppp_clean_43_2000_2014.dta"
assert _merge==3
drop _merge

gen double ppp_usd_factor = xr_lcu_per_usd/ppp_lcu_per_intl
assert ppp_usd_factor>0 & !missing(ppp_usd_factor)

* Ricci equation (3): industry-normalized PPP.
gen double va_ppp_raw = va_usd*ppp_usd_factor
bysort year industry: egen double world_va_usd_j = total(va_usd)
bysort year industry: egen double world_va_ppp_raw_j = total(va_ppp_raw)

gen double industry_ppp_scale = world_va_usd_j/world_va_ppp_raw_j
gen double ppp_usd_factor_ind = ppp_usd_factor*industry_ppp_scale
gen double va_ppp = va_usd*ppp_usd_factor_ind
gen double erdi = 1/ppp_usd_factor_ind

bysort year industry: egen double world_va_ppp_j = total(va_ppp)
gen double ppp_norm_error = abs(world_va_ppp_j-world_va_usd_j) / max(abs(world_va_usd_j),1e-12)
summarize ppp_norm_error, meanonly
assert r(max)<1e-8

bysort year industry: egen double world_hours_j = total(hours)
bysort year industry: egen double world_labcomp_j = total(labcomp_usd)
assert world_hours_j>0

gen double labor_h = (va_ppp/world_va_usd_j)*world_hours_j
gen byte ricci_eligible = go_usd>1e-10 & va_usd>1e-10 & hours>0 & labor_h>0

* Ricci equation (2): MEV.
bysort year: egen double world_va = total(va_usd)
bysort year: egen double world_hours = total(hours)
bysort year: egen double world_labcomp = total(labcomp_usd)
gen double mev = world_va/world_hours

gen double va_per_h_world_ind   = world_va_usd_j/world_hours_j
gen double wage_per_h_world_ind = world_labcomp_j/world_hours_j
gen double wage_per_h_world     = world_labcomp/world_hours
gen double wage_per_h_nat_ind   = labcomp_usd/labor_h if labor_h>0

* Equations (7)-(9).
gen double coef_inter = va_per_h_world_ind-mev
gen double coef_intra = (erdi-1)*va_per_h_world_ind
gen double coef_total = coef_inter+coef_intra

gen double coef_wage_inter = wage_per_h_world_ind-wage_per_h_world
gen double coef_wage_intra = wage_per_h_nat_ind-wage_per_h_world_ind if labor_h>0
gen double coef_wage = coef_wage_inter+coef_wage_intra if labor_h>0

gen double coef_profit_inter = coef_inter-coef_wage_inter
gen double coef_profit_intra = coef_intra-coef_wage_intra if labor_h>0
gen double coef_profit = coef_total-coef_wage if labor_h>0

* Equation-(6) identity check.
gen double coef_total_direct = va_usd/labor_h-mev if labor_h>0
gen double coef_error = abs(coef_total-coef_total_direct) / max(abs(coef_total_direct),1e-12) if labor_h>0
summarize coef_error, meanonly
assert r(max)<1e-8

* Transfer per dollar of exporter gross output.
gen double rate_total        = coef_total*labor_h/go_usd if ricci_eligible
gen double rate_inter        = coef_inter*labor_h/go_usd if ricci_eligible
gen double rate_intra        = coef_intra*labor_h/go_usd if ricci_eligible
gen double rate_wage         = coef_wage*labor_h/go_usd if ricci_eligible
gen double rate_wage_inter   = coef_wage_inter*labor_h/go_usd if ricci_eligible
gen double rate_wage_intra   = coef_wage_intra*labor_h/go_usd if ricci_eligible
gen double rate_profit       = coef_profit*labor_h/go_usd if ricci_eligible
gen double rate_profit_inter = coef_profit_inter*labor_h/go_usd if ricci_eligible
gen double rate_profit_intra = coef_profit_intra*labor_h/go_usd if ricci_eligible

keep year country industry va_usd go_usd hours labcomp_usd ///
     hours_countryavg hours_pwt lab_from_comp ricci_ineligible ricci_eligible ///
     xr_lcu_per_usd ppp_lcu_per_intl ppp_usd_factor ppp_usd_factor_ind ///
     erdi labor_h mev coef_* rate_*

isid year country industry
compress
save "$WORK/ricci_coefficients_43_2000_2014.dta", replace

********************************************************************************
* 8. 43-COUNTRY BILATERAL VALUE TRANSFERS
********************************************************************************

use "$WORK/ricci_coefficients_43_2000_2014.dta", clear
rename country exporter
keep year exporter industry go_usd va_usd hours ricci_ineligible ricci_eligible rate_*
save "$WORK/transfer_rates_43_2000_2014.dta", replace

use "$WORK/bilateral_trade_43.dta", clear
merge m:1 year exporter industry using "$WORK/transfer_rates_43_2000_2014.dta"

preserve
    keep if _merge==3 & missing(rate_total) & abs(trade_usd)>1e-8
    count
    if r(N)>0 {
        sort year exporter industry importer
        export delimited using "$OUT/ineligible_source_positive_exports_43_2000_2014.csv", replace
    }
restore

keep if _merge==3
drop _merge
drop if missing(rate_total)

gen double tr_total        = rate_total*trade_usd
gen double tr_inter        = rate_inter*trade_usd
gen double tr_intra        = rate_intra*trade_usd
gen double tr_wage         = rate_wage*trade_usd
gen double tr_wage_inter   = rate_wage_inter*trade_usd
gen double tr_wage_intra   = rate_wage_intra*trade_usd
gen double tr_profit       = rate_profit*trade_usd
gen double tr_profit_inter = rate_profit_inter*trade_usd
gen double tr_profit_intra = rate_profit_intra*trade_usd

assert abs(tr_total-tr_inter-tr_intra)<1e-6
assert abs(tr_total-tr_wage-tr_profit)<1e-6

keep year exporter importer industry trade_usd tr_*
compress
save "$WORK/bilateral_value_transfers_43_2000_2014.dta", replace

********************************************************************************
* 9. TABLE 3: COUNTRY NET TRANSFERS, 43-COUNTRY SYSTEM
********************************************************************************

use "$WORK/bilateral_value_transfers_43_2000_2014.dta", clear
preserve
    keep year exporter tr_*
    rename exporter country
    collapse (sum) tr_*, by(year country)
    tempfile export_credits43
    save `export_credits43'
restore

keep year importer tr_*
rename importer country
foreach v of varlist tr_* {
    replace `v'=-`v'
}
collapse (sum) tr_*, by(year country)
append using `export_credits43'
collapse (sum) tr_*, by(year country)

rename tr_total        net_transfer
rename tr_inter        net_inter
rename tr_intra        net_intra
rename tr_wage         net_wage
rename tr_wage_inter   net_wage_inter
rename tr_wage_intra   net_wage_intra
rename tr_profit       net_profit
rename tr_profit_inter net_profit_inter
rename tr_profit_intra net_profit_intra

preserve
    use "$WORK/wiod2016_accounts_43.dta", clear
    collapse (sum) va_usd, by(year country)
    tempfile country_va43
    save `country_va43'
restore
merge 1:1 year country using `country_va43'
assert _merge==3
drop _merge

gen double net_transfer_pct_va=100*net_transfer/va_usd
gen byte wiod_release2016=1

bysort year: egen double net_sum_check=total(net_transfer)
assert abs(net_sum_check)<1e-5
drop net_sum_check

sort country year
isid year country
count
assert r(N)==43*15
save "$OUT/table3_43countries_2000_2014_long.dta", replace

********************************************************************************
* 10. TABLE 2: ORIGINAL 12 RICCI REGIONS + ADDITIONAL CORE, 2000-2014
********************************************************************************

gen str20 region=""
replace region="North America"   if inlist(country,"CAN","USA")
replace region="North EMU"       if inlist(country,"AUT","BEL","DEU","FIN","FRA","LUX","NLD")
replace region="South EMU"       if inlist(country,"CYP","ESP","GRC","IRL","ITA","MLT","PRT")
replace region="North Europe"    if inlist(country,"DNK","GBR","SWE")
replace region="East Europe"     if inlist(country,"BGR","CZE","EST","HUN","LTU")
replace region="East Europe"     if inlist(country,"LVA","POL","ROU","SVK","SVN")
replace region="Latin America"   if inlist(country,"BRA","MEX")
replace region="China"           if country=="CHN"
replace region="India"           if country=="IND"
replace region="North East Asia" if inlist(country,"JPN","KOR")
replace region="Additional Core" if inlist(country,"CHE","HRV","NOR")
replace region="Other Asia"      if inlist(country,"IDN","TUR","TWN")
replace region="Russia"          if country=="RUS"
replace region="Australia"       if country=="AUS"
assert region!=""

collapse (sum) net_transfer net_inter net_intra ///
               net_wage net_wage_inter net_wage_intra ///
               net_profit net_profit_inter net_profit_intra va_usd, by(year region)

gen double net_transfer_pct_va=100*net_transfer/va_usd
gen byte wiod_release2016=1
sort region year
isid year region
count
assert r(N)==13*15
save "$OUT/table2_13regions_all43_2000_2014_long.dta", replace

preserve
    bysort year: egen double global_va=total(va_usd)
    keep if net_transfer>0
    collapse (sum) total_net_transfer=net_transfer total_inter=net_inter ///
                   total_intra=net_intra total_wage=net_wage total_profit=net_profit ///
             (max) global_va, by(year)
    gen double total_net_pct_global_va=100*total_net_transfer/global_va
    gen double industrial_specialization_pct=100*total_inter/total_net_transfer
    gen double absolute_rent_pct=100*total_intra/total_net_transfer
    gen double wage_differential_pct=100*total_wage/total_net_transfer
    gen double profit_differential_pct=100*total_profit/total_net_transfer
    gen byte wiod_release2016=1
    save "$OUT/table2_13regions_all43_decomposition_2000_2014.dta", replace
restore

********************************************************************************
* 11. ESTIMATE CURRENT ROW TRANSFER RATES FROM THE 43 OBSERVED COUNTRIES
********************************************************************************

use "$WORK/ricci_coefficients_43_2000_2014.dta", clear
keep if go_usd>1e-10 & !missing(rate_total)
local RATES "rate_total rate_inter rate_intra rate_wage rate_wage_inter rate_wage_intra rate_profit rate_profit_inter rate_profit_intra"
foreach v of local RATES {
    gen double num_`v'=`v'*go_usd
}
collapse (sum) num_rate_total num_rate_inter num_rate_intra ///
               num_rate_wage num_rate_wage_inter num_rate_wage_intra ///
               num_rate_profit num_rate_profit_inter num_rate_profit_intra ///
               observed_output=go_usd, by(year industry)
foreach v of local RATES {
    gen double `v'=num_`v'/observed_output
}
gen str3 exporter="ROW"
gen byte row_rate_estimated=1
keep year exporter industry observed_output row_rate_estimated `RATES'
isid year exporter industry
save "$OUT/row_transfer_rates_estimated_2000_2014.dta", replace
export delimited using "$OUT/row_transfer_rates_estimated_2000_2014.csv", replace

tempfile rowrates
save `rowrates'

use "$WORK/transfer_rates_43_2000_2014.dta", clear
gen byte row_rate_estimated=0
gen double observed_output=.
append using `rowrates'
isid year exporter industry
save "$WORK/transfer_rates_44_ROW_est_2000_2014.dta", replace

********************************************************************************
* 12. 43 OBSERVED COUNTRIES + ESTIMATED CURRENT ROW BILATERAL TRANSFERS
********************************************************************************

use "$WORK/bilateral_trade_44.dta", clear
merge m:1 year exporter industry using "$WORK/transfer_rates_44_ROW_est_2000_2014.dta"

preserve
    keep if (_merge!=3 | missing(rate_total)) & abs(trade_usd)>1e-8
    count
    if r(N)>0 {
        sort year exporter industry importer
        export delimited using "$OUT/unusable_positive_trade_44_2000_2014.csv", replace
    }
restore

keep if _merge==3
drop _merge
drop if missing(rate_total)

gen double tr_total        = rate_total*trade_usd
gen double tr_inter        = rate_inter*trade_usd
gen double tr_intra        = rate_intra*trade_usd
gen double tr_wage         = rate_wage*trade_usd
gen double tr_wage_inter   = rate_wage_inter*trade_usd
gen double tr_wage_intra   = rate_wage_intra*trade_usd
gen double tr_profit       = rate_profit*trade_usd
gen double tr_profit_inter = rate_profit_inter*trade_usd
gen double tr_profit_intra = rate_profit_intra*trade_usd

assert abs(tr_total-tr_inter-tr_intra)<1e-6
assert abs(tr_total-tr_wage-tr_profit)<1e-6

keep year exporter importer industry trade_usd row_rate_estimated tr_*
save "$WORK/bilateral_value_transfers_44_ROW_est_2000_2014.dta", replace

********************************************************************************
* 13. TABLE 3 INCLUDING CURRENT WIOD ROW (TRANSFER RATE ESTIMATED)
********************************************************************************

use "$WORK/bilateral_value_transfers_44_ROW_est_2000_2014.dta", clear
preserve
    keep year exporter tr_*
    rename exporter country
    collapse (sum) tr_*, by(year country)
    tempfile export_credits44
    save `export_credits44'
restore

keep year importer tr_*
rename importer country
foreach v of varlist tr_* {
    replace `v'=-`v'
}
collapse (sum) tr_*, by(year country)
append using `export_credits44'
collapse (sum) tr_*, by(year country)

rename tr_total        net_transfer
rename tr_inter        net_inter
rename tr_intra        net_intra
rename tr_wage         net_wage
rename tr_wage_inter   net_wage_inter
rename tr_wage_intra   net_wage_intra
rename tr_profit       net_profit
rename tr_profit_inter net_profit_inter
rename tr_profit_intra net_profit_intra

preserve
    use "$WORK/wiod2016_accounts_44.dta", clear
    collapse (sum) va_usd, by(year country)
    tempfile country_va44
    save `country_va44'
restore
merge 1:1 year country using `country_va44'
assert _merge==3
drop _merge

gen double net_transfer_pct_va=100*net_transfer/va_usd
gen byte row_estimated=country=="ROW"
gen byte wiod_release2016=1

bysort year: egen double net_sum_check=total(net_transfer)
assert abs(net_sum_check)<1e-5
drop net_sum_check

sort country year
isid year country
count
assert r(N)==44*15
save "$OUT/table3_44regions_ROW_est_2000_2014_long.dta", replace

********************************************************************************
* 14. TABLE 2 INCLUDING CURRENT WIOD ROW (TRANSFER RATE ESTIMATED)
********************************************************************************

gen str24 region=""
replace region="North America"   if inlist(country,"CAN","USA")
replace region="North EMU"       if inlist(country,"AUT","BEL","DEU","FIN","FRA","LUX","NLD")
replace region="South EMU"       if inlist(country,"CYP","ESP","GRC","IRL","ITA","MLT","PRT")
replace region="North Europe"    if inlist(country,"DNK","GBR","SWE")
replace region="East Europe"     if inlist(country,"BGR","CZE","EST","HUN","LTU")
replace region="East Europe"     if inlist(country,"LVA","POL","ROU","SVK","SVN")
replace region="Latin America"   if inlist(country,"BRA","MEX")
replace region="China"           if country=="CHN"
replace region="India"           if country=="IND"
replace region="North East Asia" if inlist(country,"JPN","KOR")
replace region="Additional Core" if inlist(country,"CHE","HRV","NOR")
replace region="Other Asia"      if inlist(country,"IDN","TUR","TWN")
replace region="Russia"          if country=="RUS"
replace region="Australia"       if country=="AUS"
replace region="Rest of World (est.)" if country=="ROW"
assert region!=""

collapse (sum) net_transfer net_inter net_intra ///
               net_wage net_wage_inter net_wage_intra ///
               net_profit net_profit_inter net_profit_intra va_usd, by(year region)

gen double net_transfer_pct_va=100*net_transfer/va_usd
gen byte row_estimated=region=="Rest of World (est.)"
gen byte wiod_release2016=1
sort region year
isid year region
count
assert r(N)==14*15
save "$OUT/table2_14regions_ROW_est_2000_2014_long.dta", replace

preserve
    bysort year: egen double global_va=total(va_usd)
    keep if net_transfer>0
    collapse (sum) total_net_transfer=net_transfer total_inter=net_inter ///
                   total_intra=net_intra total_wage=net_wage total_profit=net_profit ///
             (max) global_va, by(year)
    gen double total_net_pct_global_va=100*total_net_transfer/global_va
    gen double industrial_specialization_pct=100*total_inter/total_net_transfer
    gen double absolute_rent_pct=100*total_intra/total_net_transfer
    gen double wage_differential_pct=100*total_wage/total_net_transfer
    gen double profit_differential_pct=100*total_profit/total_net_transfer
    gen byte row_included_estimate=1
    gen byte wiod_release2016=1
    save "$OUT/table2_14regions_ROW_est_decomposition_2000_2014.dta", replace
restore

********************************************************************************
* 15. EXCEL EXPORTS
********************************************************************************

local XLSX "$OUT/ricci_2000_2014_wiod2016_all43.xlsx"

use "$OUT/table3_43countries_2000_2014_long.dta", clear
export excel using "`XLSX'", sheet("T3_43_long") firstrow(variables) replace
preserve
    keep country year net_transfer net_transfer_pct_va
    reshape wide net_transfer net_transfer_pct_va, i(country) j(year)
    export excel using "`XLSX'", sheet("T3_43_wide") firstrow(variables) sheetmodify
restore

use "$OUT/table2_13regions_all43_2000_2014_long.dta", clear
export excel using "`XLSX'", sheet("T2_13reg_long") firstrow(variables) sheetmodify
preserve
    keep region year net_transfer net_transfer_pct_va
    reshape wide net_transfer net_transfer_pct_va, i(region) j(year)
    export excel using "`XLSX'", sheet("T2_13reg_wide") firstrow(variables) sheetmodify
restore
use "$OUT/table2_13regions_all43_decomposition_2000_2014.dta", clear
export excel using "`XLSX'", sheet("T2_13reg_decomp") firstrow(variables) sheetmodify

use "$OUT/table3_44regions_ROW_est_2000_2014_long.dta", clear
export excel using "`XLSX'", sheet("T3_44_ROW_long") firstrow(variables) sheetmodify
preserve
    keep country year net_transfer net_transfer_pct_va
    reshape wide net_transfer net_transfer_pct_va, i(country) j(year)
    export excel using "`XLSX'", sheet("T3_44_ROW_wide") firstrow(variables) sheetmodify
restore

use "$OUT/table2_14regions_ROW_est_2000_2014_long.dta", clear
export excel using "`XLSX'", sheet("T2_14_ROW_long") firstrow(variables) sheetmodify
preserve
    keep region year net_transfer net_transfer_pct_va
    reshape wide net_transfer net_transfer_pct_va, i(region) j(year)
    export excel using "`XLSX'", sheet("T2_14_ROW_wide") firstrow(variables) sheetmodify
restore
use "$OUT/table2_14regions_ROW_est_decomposition_2000_2014.dta", clear
export excel using "`XLSX'", sheet("T2_14_ROW_decomp") firstrow(variables) sheetmodify

********************************************************************************
* 16. FINAL MESSAGE
********************************************************************************

di as result "Completed Ricci-style WIOD 2016 exercise for 2000-2014 using all 43 observed countries."
di as result "Observed sample: all 43 WIOD countries; CHE, HRV, and NOR remain separately observed."
di as result "ROW = current WIOD ROW only; its transfer rates are estimated from the 43 observed-country industry rates."
di as result "2016 SEA hours use H_EMPE*EMP/EMPE; PWT avh is used only where SEA has no hours series (China)."
di as result "Workbook: `XLSX'"


/****************************************************************************************
 PART 3. ORIGINAL WIOD RELEASE: TABLE 2 + PWT -> COMPACT CALIBRATION DATA
 Source: ricci_calibration_1995_2009_aggregate.do
 Data-construction sections only; estimation begins in 02_estimation.do.
****************************************************************************************/

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


/****************************************************************************************
 PART 4. WIOD 2016 RELEASE: TABLE 2 + PWT -> COMPACT CALIBRATION DATA
 Source: ricci_calibration_wiod2016_2000_2014_aggregate.do
 Data-construction sections only; estimation begins in 02_estimation.do.
****************************************************************************************/

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


/****************************************************************************************
 DATA CONSTRUCTION COMPLETE

 Upload these compact Stata inputs with the public replication package:

   calibration_1995_2009_aggregate/calib_inputs_1995_2009_table2_reg.dta
   calibration_1995_2009_aggregate/calib_inputs_1995_2009_table2_ROW.dta

   calibration_2016release_2000_2014/calib_inputs_2000_2014_table2_reg.dta
   calibration_2016release_2000_2014/calib_inputs_2000_2014_table2_ROW.dta
   calibration_2016release_2000_2014/eligible_regions_2000_2014.dta

 The large WIOD/SEA source files are not required by 02_estimation.do.
****************************************************************************************/
