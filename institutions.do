clear 
cd "/Users/brendanbrundage/Library/CloudStorage/OneDrive-MorehouseCollege/Desktop/Research/BrundageDarityTavani"

import excel using "P_Institutional Quality", firstrow
drop TimeCode
rename Time year
destring year ControlofCorruptionGovernan GovernmentEffectivenessGover PoliticalStabilityGovernance RegulatoryQualityGovernance RuleofLawGovernancescore VoiceandAccountabilityGover, replace force
keep if year==2024
rename CountryName country
rename CountryCode code
rename ControlofCorruptionGovernan corruption
rename GovernmentEffectivenessGover gov
rename PoliticalStabilityGovernance stabilitu
rename stabilitu stability
rename RegulatoryQualityGovernance regulatory
rename RuleofLawGovernancescore rol
rename VoiceandAccountabilityGover voice

save WGI, replace 

clear

import excel using "soveriegn", firstrow 

replace status="1" if status=="Independent"
replace status="0" if status!="1"

drop dep 

destring status, replace force


save status, replace

clear

use WGI, clear

merge 1:1 code using status

keep if _merge==3

drop _merge

keep if status==1

save WGIstatus, replace

clear

import delimited world-regions-according-to-the-world-bank.csv, varnames(1)
drop if code==""

drop year

save region, replace

use WGIstatus, clear

merge 1:1 code using region

keep if _merge==3

rename worldregionsaccordingtowb region

drop _merge

gen region2=region

drop if code=="SMR"

replace region2="Caribbean" if code=="JAM" | code=="SUR" | code=="GUY" | code=="KNA" | code=="BRB" | code=="CUB" | code=="HTI" | code=="DOM" | code=="LCA" | code=="BHS" | code=="GRD" | code=="VCT" | code=="TTO" | code=="DMA" | code=="ATG" 

replace region2="Western Countries" if code=="USA" | code=="CAN" | code=="ISL" | code=="GBR" | code=="ESP" | code=="CHE" | code=="BEL" | code=="ITA" | code=="CYP" | code=="GRC" | code=="POL" | code=="IRL" | code=="FIN" | code=="SWE" | code=="FRA" | code=="MCO" | code=="LIE" | code=="" | code=="AUT" | code=="LUX" | code=="PRT" | code=="DNK" | code=="NLD" | code=="DEU" | code=="NOR" | code=="NZL" | code=="AUS" | code=="AND" | code=="EST" | code=="" | code=="" | code=="LTU" | code=="ROU" | code=="SVN" | code=="LVA" | code=="SVK" | code=="HUN" | code=="GRL" | code=="MLT" | code=="HRV" | code=="BGR" | code=="CZE"

replace region2="South Asia (WB)" if code=="PAK" | code=="AFG"

pca corruption gov stability regulatory rol voice
predict institution

foreach var of varlist corruption gov stability regulatory rol voice status{
	bysort region2: egen r_`var'=mean(`var')
}


bysort region2: egen r_inst=mean(institution)

bysort region2: gen keep=_n

keep if keep==1

gsort -r_inst

foreach var of varlist r_corruption r_gov r_stability r_regulatory r_rol r_voice r_inst{
	replace `var'=round(`var', .01)
}

save institutions_complete, replace

drop year country code corruption gov stability regulatory rol voice status entity region institution keep r_status

texsave region2-r_inst using institutions.tex, title(Institutional Rankings by Region) replace

clear

import excel using gdp, firstrow

drop GDPpercapitaPPPconstant20 GDPpercapitaPPPcurrentint TimeCode 

rename CountryName country

rename CountryCode code 

rename Time year

rename GiniindexSIPOVGINI gini
rename GDPpercapitaconstant2015US gdp_pc
rename GDPpercapitaconstantLCUN gdp_pc_LCU
rename GDPpercapitagrowthannual pc_growth

destring gini gdp_pc gdp_pc_LCU pc_growth, replace force

bysort country: egen has_pre1971 = max(year < 1971 & !missing(gdp_pc))
drop if has_pre1971 == 0
drop has_pre1971

sort country year

drop if year<1970

save gdp, replace

clear

use region, clear

merge 1:m code using gdp

drop if _merge==1

sort country year

rename worldregionsaccordingtowb region

drop _merge

gen region2=region

drop if code=="SMR"

replace region2="Caribbean" if code=="JAM" | code=="SUR" | code=="GUY" | code=="KNA" | code=="BRB" | code=="CUB" | code=="HTI" | code=="DOM" | code=="LCA" | code=="BHS" | code=="GRD" | code=="VCT" | code=="TTO" | code=="DMA" | code=="ATG" 

replace region2="Western Countries" if code=="USA" | code=="CAN" | code=="ISL" | code=="GBR" | code=="ESP" | code=="CHE" | code=="BEL" | code=="ITA" | code=="CYP" | code=="GRC" | code=="POL" | code=="IRL" | code=="FIN" | code=="SWE" | code=="FRA" | code=="MCO" | code=="LIE" | code=="" | code=="AUT" | code=="LUX" | code=="PRT" | code=="DNK" | code=="NLD" | code=="DEU" | code=="NOR" | code=="NZL" | code=="AUS" | code=="AND" | code=="EST" | code=="" | code=="" | code=="LTU" | code=="ROU" | code=="SVN" | code=="LVA" | code=="SVK" | code=="HUN" | code=="GRL" | code=="MLT" | code=="HRV" | code=="BGR" | code=="CZE"

replace region2="South Asia (WB)" if code=="PAK" | code=="AFG"

replace region2="Middle East, North Africa, Afghanistan and Pakistan (WB)" if code=="TUR" | code=="GEO"

drop if code=="BMU"

foreach var of varlist gini gdp_pc gdp_pc_LCU {
	bysort region2 year: egen r_`var'=mean(`var')
}


bysort region2: egen r_pc_growth=mean(pc_growth)

foreach var of varlist r_gini r_gdp_pc r_gdp_pc_LCU r_pc_growth{
	replace `var'=round(`var', .01)
}

bysort code: egen c_pc_growth=mean(pc_growth)

bysort region2 year: gen keep=_n

keep if keep==1

drop entity country code region keep pc_growth gini gdp_pc gdp_pc_LCU r_gini r_gdp_pc_LCU

sort year region2

drop if year==2025

egen region_id=group(region2)

label define REGION 1 "Caribbean (1.42%)" 2 "East Asia and Pacific (1.71%)" 3  "Latin America (1.60%)" 4 "Middle East (0.03%)" 5 "South Asia (3.35%)" 6 "Sub-Saharan Africa (1.12%)" 7 "Western Countries (1.79%)"

label values region_id REGION




xtset region_id year

xtline r_gdp_pc, overlay ///
    legend(cols(2) size(small)) ///
    ytitle("GDP per capita") ///
    xtitle("Year") ///
    title("GDP per Capita by Region")  yscale(log) ylabel(1875 3750 7500 15000 30000 60000)
	
graph export "gdp_pc_regions.png", width(2400) replace	


***CAGR	
*Caribbean=1.42, East Asia=1.71, LA=1.60, ME=0.03, South Asia=3.35, SSA=1.12, Western=1.79.



preserve
drop r_pc_growth c_pc_growth
keep if inlist(year, 1970, 2024)
reshape wide r_gdp_pc, i(region2) j(year)
gen growth_1970_2024 = 100 * ((r_gdp_pc2024 / r_gdp_pc1970)^(1/54) - 1)

restore
*twoway line r_gdp_pc year, ///
    ytitle("GDP per capita") ///
    xtitle("Year") ///
    title("GDP per Capita by Region") ///
    legend(cols(3))
