cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition\stata_meeting_presentation"
set seed 1234

cap ado uninstall nopo
net install nopo, replace from("https://raw.githubusercontent.com/mhamjediers/nopo_decomposition/master/")

cap progrdem drop logout
run "logout.ado"
cap program drop noheadlog
run "noheadlog.ado"

cap log close
set logmsg off


*After other projekt 02_gen.do
use "real_example.dta", clear

quietly { // preparation of dataset 
bys pid (syear): gen r=runiform()
bys pid (syear): egen maxr = max(r)
keep if r==maxr
drop r maxr
tab syear

// group var
gen grp = woman + immi*2 + 1
lab var grp "Group indicator"
lab def grp 1 "Native Men" 2 "Native Women" 3 "Immigrant Men" 4 "Immigrant Women"
lab val grp grp

// dummies for descriptive table
tab ger_skills, gen(ger_skills_) // dummies for descriptive table

// dummies/quadratic for KBO
for any age_c edu isco2d lmexp: tab X, gen(X_)
gen lmexpc2 = lmexpc*lmexpc

keep if inlist(grp,1,4) 
}

global pred age_c married edu lmexp parttime isco2d	
	
********
* Stand alone
********
	
*First output
logout using "log/standalone.log", cmd: nopo decomp wage age_c married edu lmexp parttime isco2d, by(grp)

*Further options
log using "log/main_options.log", replace
nopo decomp wage ${pred}, by(grp) bref(grp == 1) swap normalize
log close
noheadlog using "log/main_options.log"


*other matching approaches
log using "log/matchings.log", replace
qui: nopo decomp wage ${pred}, by(grp)
qui: est store em
qui: nopo decomp wage ${pred}, by(grp) kmatch(ps)
qui: est store ps
qui: nopo decomp wage ${pred}, by(grp) kmatch(md)
qui: est store md 
qui: nopo decomp wage ${pred}, by(grp) kmatch(ps) kmopt(pscmd(probit) bw(0.0001))
qui: est store ps_probbw
log close
noheadlog using "log/matchings.log"

log using "log/match_table.log", replace
esttab em ps md ps_probbw, se nonumbers nonotes ///
	 mtitles("exact" "prop. score" "multi. dist." "probit ps") ///
	stats(nA mshareuwA nB mshareuwB, label("N(A)" "% matched A" "N(B)" "% matched B"))
log close
noheadlog using "log/match_table.log"

********
* Postestimation to kmatch
********
log using "log/after_kmatch.log", replace
qui: kmatch ps grp ${pred} (wage), ///
	tval(4) atc att bw(0.001) pscmd(probit) generate wgenerate replace 
nopo decomp
log close
noheadlog using "log/after_kmatch.log"


*Read out the used kmatch-command for potential adjustment or kmatch postestimation
log using "log/readout_kmatch.log", replace
qui: nopo decomp wage ${pred}, by(grp) kmatch(ps)
display "`e(kmatch_cmdline)'"
log close


********
* Postestimation
********

*Summarize table
log using "log/nopo_summarize.log", replace
qui: nopo decomp wage ${pred}, by(grp)
nopo summarize wage age married edu_1 edu_2 edu_3, label
log close
noheadlog using "log/nopo_summarize.log"


*contribution to DA/DB
logout using "log/dadb.log", cmd: nopo dadb edu
graph export "log/dadb.pdf", replace


*Over the distribution
logout using "log/gapoverdist.log", cmd: nopo gapoverdist
graph export "log/gapoverdist.pdf", replace




