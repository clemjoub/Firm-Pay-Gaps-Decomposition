* File paths
	global path=""
	global varlist1 ="   meanend_d         meanend_s        meanendfit          sdendfit            nend_d          nend_s         muentry_d       	muentry_s        muentrydif          nentry_d          nentry_s          muexit_d          muexit_s         muexitdif           nexit_d           nexit_s         meanint_d       meanint_s         meanintfit"

* Number of bootstrap replications
	global nbbootstraps = 1 //Set to one when not doing bootstrapping
		if $nbbootstraps >10 {
			set matsize $nbbootstraps
		}

* Number of gaps to decompose 
	global nbscenarios = 1 // Run through several decompositions

* Decomposition order (default is all possible orders)
	global PrefOrder=1 // Do decomposition for the preferred counterfactual order only
	
********************************************************************************

forvalues scenario=1(1)$nbscenarios {

	* Matrix containing the decomposition results
	mat Bootstrap = J($nbbootstraps,5,.)
	matrix colnames Bootstrap = attrition salgrowth entrysal gradecomp pre87 

	local i=1

	forvalues boot = 1(1) $nbbootstraps { 
			local rownames=""

			
			* Define the gaps to be decomposed (between a group of interest and a reference group in a given year)
			* groups are defined in terms of "types" which identify gender, nationality and entry grade (see bottom of dofile)

			if `scenario'==1 {
				loc scenario_name="AggGender"
				loc refgroup="(type>0 & type<13)"
				loc group="(type>12 & type<25)"
				local simpath="${path}\Simulations\AggregateGap\"
				local outpath = "$path\Decomposition\AggregateGap"
				local gap=1 // 1 is gender, 2 is nationality
				local caption = "Aggregate Gender Gap decomposition, all grades"
				global year = 2015
				
			}
			else if `scenario'==2 {
				loc scenario_name="AggGenderGE+"
				loc refgroup=	"(type>2  & type<7)  | (type>8  & type<13)"
				loc group=		"(type>14 & type<19) | (type>20 & type<25)"
				local simpath="${path}\Simulations\AggregateGap\"
				local outpath = "$path\Decomposition\AggregateGap"
				local gap=1 // 1 is gender, 2 is nationality
				local caption = "Aggregate Gender Gap decomposition, grades GE and above"
				global year = 2015
			}
			else if `scenario'==3 {
				loc scenario_name="AggCountry"
				loc refgroup= "(type>0 & type<7) | (type>12 & type<19)"
				loc group=  "(type>6 & type<13) | (type>18 & type<25) "
				local simpath="${path}\Simulations\AggregateGap\"
				local gap=2 // 1 is gender, 2 is nationality
				local caption = "Aggregate Nationality Gap decomposition, all grades"
				global year = 2015
			}
			else if `scenario'==4 {
				loc scenario_name="AggCountryGE+"
				loc refgroup="(type>2 & type<7) | (type>14 & type<19)"
				loc group="(type>8 & type<13) | (type>20 & type<25) "
				local simpath="${path}\Simulations\AggregateGap\"
				local outpath = "$path\Decomposition\AggregateGap"
				local gap=2 // 1 is gender, 2 is nationality
				local caption = "Aggregate Nationality Gap decomposition, grades GE and above"
				global year = 2015
			}

					
			// Loop through possible orders of counterfactuals in the decomposition
			// Note: Order of counterfactuals in simulation file names is
			// nolegacy(10000), nocompgap(1000), noentrygap(100), nogrowthgap(10), noexitgap(1) 
			
			// n1, n2, n3, n4, n5 denote the order of counterfactuals  where 5 is nolegacy, 1 is no exit
			// e.g. n1=1, n2=2, n3=3, n4=4, n5=5 denotes the following order of decomposition:
			// exit, growth, entry salaries, composition, legacy ("preferred" order)
			
			//

			forvalues n1=1(1)5 { // first counterfactual 
			forvalues n2=1(1)5 { // second counterfactual
			if `n2'==`n1' {
				continue
				}
			forvalues n3=1(1)5 { // third counterfactual
			if `n3' == `n1' | `n3' == `n2' {
				continue
				}
			forvalues n4=1(1)5 { // fourth counterfactual
			if `n4' == `n1' | `n4' == `n2' | `n4'==`n3' {
				continue
				}
			forvalues n5=1(1)5 { // fifth counterfactual
			if `n5' == `n1' | `n5' == `n2' | `n5'==`n3' | `n5'==`n4' {
				continue
				}
			if $PrefOrder==1 & (`n1' != 1 | `n2' != 2 | `n3' != 3 | `n4' != 4) { //Only do preferred order
				continue
			}

			*Create filenames to be used in each counterfactual
			local c0=10^5*`gap' // 100000 or 200000: baseline scenario (no counterfactual)
			local c1=`c0'+10^(`n1'-1) //e.g. for preferred order: 100000 + 10^(1-1) = 100001 => counterfactual in which exit gap is shut down
			local c2=`c1'+10^(`n2'-1) //e.g. for preferred order: 100001+ 10^(2-1) = 100011 => counterfactual in which exit gap and growth gap are shut down
			local c3=`c2'+10^(`n3'-1) //etc.
			local c4=`c3'+10^(`n4'-1)
			local c5=`c4'+10^(`n5'-1)

			local list = "`c0'" + " `c1'"+ " `c2'"+ " `c3'"+ " `c4'"+ " `c5'"
			di "`list'"
			
			* matrix row names indicate the order of counterfactuals ("preferred"= 12345)
			local rownames="`rownames'" + " `n1'" + "`n2'" +"`n3'" +"`n4'" + "`n5'"
			di "`rownames'"

			loc pctgap_cum=1
			local c=1
			local lagcounterfactual=100000*`gap'
			foreach counterfactual in `list' {
					di "previous counterfactual: `lagcounterfactual'"
					di "current counterfactual: `counterfactual'"
					local factor= 1+log(`counterfactual'-`lagcounterfactual')/log(10) //identifies the factor that is being shifted in this counterfactual
					di "Shifting factor number: `factor'"
					 

					*load simulated salary under the counterfactual  
					if $nbbootstraps==1 {
						quietly infile   year  type $varlist1  using "`simpath'\simulationfit`counterfactual'.txt", clear
					}
					else if $nbbootstraps>1 {
						quietly infile   year  type $varlist1  using "`simpath'\simulationfit`counterfactual'_bt`boot'.txt", clear
						sort year type
					}
					
					
					if `scenario'<5 & ${year}<2016{
						sort year
						merge year using "C:\Users\wb452275\OneDrive - WBG\HR Longitudinal Database\PolicyPaperFiles_Spring2017\DATA\CPI.dta"
						keep if _merge==3
						drop _merge

						foreach var of varlist meanint_d  meanint_s  meanintfit  meanend_d meanend_s meanendfit sdendfit {
							replace `var'=`var'*CPI
						}
					}
					
					keep if year==$year
					keep year type meanend_s nend_s meanend_d nend_d
					g aggsalary=meanend_s*nend_s
					
					
					*compute average salary for the group of interest (e.g. women)
					
						
					g temp = aggsalary if `group'
						egen num = total(temp) 
						drop temp
					g temp = nend_s if `group'
						egen denom = total(temp)
						drop temp
					g meansal_group=num/denom			
						drop num denom	

					*compute average salary for the group of reference (e.g men)

					g temp = aggsalary if `refgroup'
						egen num = total(temp) 
						drop temp
					g temp = nend_s if `refgroup'
						egen denom = total(temp)
						drop temp
					g meansal_refgroup=num/denom
						drop num denom

					*compute the gap between group of interest and of reference
					
					g gap=meansal_refgroup-meansal_group
					if `counterfactual'==`lagcounterfactual' {
						local baselinegap=gap[1]
					}
					local gap`counterfactual'=gap[1]		
					local pctgap=`pctgap_cum'-gap[1]/`baselinegap'
					
					di "% gap: `pctgap'"
					di "c: `c'"
					di "i: `i'"
					di "factor: `factor'"
					*store results in a matrix
					
					if `c' >1 {
					
						mat Bootstrap[`i',`factor']=int((`pctgap_cum'-gap[1]/`baselinegap')*1000)/10
					}
					local pctgap_cum=gap[1]/`baselinegap'
					local lagcounterfactual = `counterfactual'
					local c=`c'+1
			}
			}
			}
			}
			}
			}
	matlist Bootstrap[1..`i',.]
	local i = `i'+1
	}

	matrix rownames Bootstrap = `rownames' 

	xsvmat Bootstrap, saving("`outpath'\BTSP_`scenario_name'.dta", replace) names(col)
	mat2txt , matrix(Bootstrap) saving("`outpath'\BTSP_`scenario_name'") replace


}

* Type definitions:
* type = 1 Male Part I Ungraded
* type = 2 Male Part I GA-GD
* type = 3 Male Part I GE
* type = 4 Male Part I GF
* type = 5 Male Part I GG
* type = 6 Male Part I GH+
* type = 7 Male Part II Ungraded
* type = 8 Male Part II GA-GD
* type = 9 Male Part II GE
* type = 10 Male Part II GF
* type = 11 Male Part II GG
* type = 12 Male Part II GH+
* type = 13 Female Part I Ungraded
* type = 14 Female Part I GA-GD
* type = 15 Female Part I GE
* type = 16 Female Part I GF
* type = 17 Female Part I GG
* type = 18 Female Part I GH+
* type = 19 Female Part II Ungraded
* type = 20 Female Part II GA-GD
* type = 21 Female Part II GE
* type = 22 Female Part II GF
* type = 23 Female Part II GG
* type = 24 Female Part II GH+
