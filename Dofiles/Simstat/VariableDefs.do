
* Use Gspells instead of spells (i.e treat ungraded years as separate spells)

	replace newspell=newGspell
	replace endspell=endGspell
	replace endspell=0 if year==2015

*Check cell sizes
	*codebook type
	g entryexit=1*(endspell!=1 & newspell!=1)+ 2*(endspell==1)+3*(newspell==1 & endspell!=1)
	g freq=1
	g exits=endspell==1 & active==1
	g entries=newspell==1 & active==1
	g eop = endspell!=1 & active==1
	g inter = newspell!=1 &  active==1	
	table year type entryexit, c(n upi)

* Salary histogram
	g salbin=int(salarycurr000/$binsize)+1
	g salgroup10=10*int(salarycurr000/10)
	
*Trim salary increases at p99.5
	g logincr_tr=logincr
	quietly su logincr if year>1987, d
	quietly su logincr if logincr>r(p99) & year>1987, d 
	replace logincr_tr=. if  logincr>r(p50) & year>1987
	*tabstat logincr_tr, by(typeyear) s(min p10 p90 p99 max)
	
*Set counter (years or tenure)
	g ct = year
	
*De-mean salaries
	bysort type ct: egen temp=mean(salarycurr000)
	g salarycurr000centered=salarycurr-temp

