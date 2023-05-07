//
// postestimation: plot gap component over dv distribution
//

/*

At the moment, the plot shows the gaps by comparing the outcome means in each quantile between 
groups. So, the mean across all these comparisons is the same as the decomposition component values
produced by nopo. But that also means that:

- At each quantile, the single component values do not add up to d
- D_A and D_B are scaled as the nopo decomp, so higher values can mean more people or larger gaps
  (though the factor is always the same due to the quantile logic: n0/(n0+n1) is the same for each
  quantile). The `rawumdiff' option circumvents the scaling and shows the absolute differences,
  but then the values do not sum up do D_A/D_B

Does that sound sensible?

*/

// plotting wrapper: 5 plots needed (one for each gap component)
cap program drop nopoplot
program define nopoplot
syntax [if] [aweight], ///
	[NQuantiles(integer 100)] ///
	[QMIN(integer 1)] ///
	[QMAX(integer 100)] ///
	[RAWUMdiff] ///
	[twtype(string)] ///
	[twopts(string asis)] ///
	[twoptsd(string asis)] ///
	[twoptsd0(string asis)] ///
	[twoptsdx(string asis)] ///
	[twoptsda(string asis)] ///
	[twoptsdb(string asis)] ///
	[SAVE(string asis)]

	quietly {

		// check if prior command was nopo
		if ("`e(cmd)'" != "nopodecomp") {
			noisily dis as error "Previous command was not nopo."
			error 301
			exit
		}

		// set input from syntax and nopo returns
		// sample
		tempvar touse
		mark `touse' `if'
		replace `touse' = 0 if !e(sample)
		// depvar
		local depvar "`e(depvar)'"
		// strata
		local strata "`e(prefix)'_strata"
		// treatment indicator (fix to 0/1)
		tempvar treat
		gen `treat' = 1 if `e(ref)'
		replace `treat' = 0 if `treat' != 1 & !mi(`e(by)')
		// ref group switch
		if ("`e(swap)'" == "1") local bref = abs(`e(ref)' - 1) // reference group for returns
			else local bref = `e(ref)'
		// support
		local support "`e(prefix)'_matched"
		// matching weights
		local mweight "`e(prefix)'_weights"
		// aweight
		if ("`e(weight)'" != "") {
			local weightexp "[aw = `e(weight)']"
		}
		else {
			tempvar w1
			gen `w1' = 1
			local weightexp "[aw = `w1']"
		}

		// set defaults
		if ("`twtype'" == "") local twtype "line"
		if (`"`twopts'"' == "") local twopts `" legend(order(1 "D" 2 "DX" 3 "D0" 4 "DA" 5 "DB") rows(1) span) yline(0) scheme(s1mono) "'
		if ("`twtype'" == "line") {
			if (`"`twoptsd'"' == "") local twoptsd "lp(solid) lw(0.5)"
			if (`"`twoptsd0'"' == "") local twoptsd0 "lp(shortdash)"
			if (`"`twoptsdx'"' == "") local twoptsdx "lp(dash)"
			if (`"`twoptsda'"' == "") local twoptsda "lp(dash_dot)"
			if (`"`twoptsdb'"' == "") local twoptsdb "lp(longdash_dot)"
		}

		// abort if quantiles > depvar groups
		qui levelsof `depvar' if `touse', local(depvarlvls)
		local nqlvls : word count `depvarlvls'
		if (`nquantiles' > `nqlvls') {
			noisily dis as error "Groups in `depvar' < quantiles requested (`nquantiles')."
			error 148
			exit
		}

		// options passthru
		local opts `"nq(`nquantiles') qmin(`qmin') qmax(`qmax') `revsign' `relative'"'

		// create plot values for each component
		// D
		tempfile d
		nopo_gapdist `depvar' if `touse' `weightexp', by(`treat') `opts' save(`d')
		// D0
		tempfile d0
		nopo_gapdist `depvar' if `touse' & `support' [aw = `mweight'] ///
			, by(`treat') `opts' save(`d0')
		// DX (requires new var containing stratum specific wages for treatment group)
		tempfile dx
		tempvar depvar_strata_mt
		bys `strata' `support': egen `depvar_strata_mt' = mean(`depvar') ///
			if `support' & `treat' == `bref'
		tempvar depvar_strata_m
		bys `strata' `support': egen `depvar_strata_m' = max(`depvar_strata_mt')
		replace `depvar_strata_m' = . if !`support'
		nopo_gapdist `depvar_strata_m' if `touse' & `support' `weightexp' ///
			, by(`treat') `opts' save(`dx')
		// DA
		tempfile da
		nopo_gapdist `depvar' if `touse' & `treat' == 1 `weightexp' ///
			, by(`support') comp(da) `rawumdiff' `opts' save(`da')
		// DB
		tempfile db
		nopo_gapdist `depvar' if `touse' & `treat' == 0 `weightexp' ///
			, by(`support') comp(db) `rawumdiff' `opts' save(`db')

		// plot
		preserve
			use "`d'", clear
			rename diff d
			foreach c in d0 dx da db {
				merge 1:1 q using "``c''", nogen
				rename diff `c'
			}
			// checking
			noisily {
				dis "Check means over component quantiles (some precision lost):"
				sum d d0 dx da db
			}
			twoway ///
				(`twtype' d q, `twoptsd') ///
				(`twtype' d0 q, `twoptsd0') ///
				(`twtype' dx q, `twoptsdx') ///
				(`twtype' da q, `twoptsda') ///
				(`twtype' db q, `twoptsdb') ///
				, `twopts'
			// save if requested
			if (`"`save'"' != "") noisily save `save', replace
		restore
	}

end

// get nopo decomposition component values over distribution of depvar
cap program drop nopo_gapdist 
program define nopo_gapdist 
syntax varname [if] [aweight], ///
	by(varlist max=1) ///
	[comp(string)] /// gap components, only relevant for D_A, D_B
	[NQuantiles(integer 100)] ///
	[QMIN(integer 1)] ///
	[QMAX(integer 100)] ///
	[RAWUMdiff] /// do not scale by N_A/N_B
	[SAVE(string asis)] * 

quietly {
	preserve

		// weight
		local weightvar = (subinstr("`exp'","=","",.))

		// subset
		if "`if'" != "" keep `if'
		
		// allow for quantile mean (or any other stat) estimation'
		// xtile aggregates quantiles if they contain constant values. Fill up to avoid empty cells.
		tempvar quantile
		gen `quantile' = .
		lab var `quantile' "Compared `varlist' quantile between groups (component-specific)"
		tempvar totweight
		forvalues i = 0/1 {
			// quantiles
		    xtile `quantile'_`i' = `varlist' if `by' == `i' [`weight'`exp'], nquantiles(`nquantiles')
			bys `by': egen `totweight' = total(`weightvar')
			bys `by' (`quantile'_`i'): replace `quantile'_`i' = ///
				ceil(sum(`weightvar') * `nquantiles' / `totweight') if `by' == `i'
		    replace `quantile' = `quantile'_`i' if `by' == `i'
			drop `totweight'
		}
		
		// collapse, use sum of weights (passed via `exp') as N
		tempvar meanq
		tempvar nq
		collapse ///
            (mean) `meanq' = `varlist' ///
            (sum) `nq' = `weightvar' ///
            if !mi(`varlist') [`weight'`exp'] ///
            , by(`by' `quantile')
		reshape wide `meanq' `nq', i(`quantile') j(`by')
		
		// gen diff
		tempvar diff
		if ("`comp'" == "da") gen `diff' = `meanq'1 - `meanq'0
			else gen `diff' = `meanq'0 - `meanq'1
		lab var `diff' "Gap"

		// scale D_A/D_B if not otherwise requested
		if (inlist("`comp'", "da", "db") & "`rawumdiff'" == "") {
			replace `diff' = `diff' * (`nq'0/(`nq'0+`nq'1))
		}

		// save temp data
		if ("`save'" != "") {
			drop if !inrange(`quantile', `qmin', `qmax')
			keep `diff' `quantile' // PERHAPS keep group values per quantile with correct labels to allow for manual processing
			rename `diff' diff
			rename `quantile' q		
			save `save'
		}
	restore
}

end

//
// postestimation: plot contribution to DA/DB by specified variable
//
/*
 Plot the contribution of each X level to DA/DB (weighted) and absolute difference in the outcome.
*/

cap program drop nopoplot2
program define nopoplot2
syntax varname [aweight], ///
	[NOSORT] /// do not sort by depvar
	[DESCending] /// sort descending (as opposed to ascending if nosort is not specified)
	[KEEPALLlevels] /// keep all levels of var (if cond. ignored)
	[FORCE] /// do not check for no. of levels in var
	[nmin(real 1)] /// minimum number of unmatched weighted obs per cat. be printed
	[twopts(string asis)] ///
	[twoptsbar(string asis)] ///
	[twoptsscatter(string asis)] ///
	[twoptsby(string asis)] ///
	[SAVE(string asis)]

	quietly {

		// check if prior command was nopo
		if ("`e(cmd)'" != "nopodecomp") {
			noisily dis as error "Previous command was not nopo."
			error 301
			exit
		}
		
		// plotvar
		tempvar plotby
		clonevar `plotby' = `varlist'
		local plotbyname `varlist' // for renaming tempvar upon save
		local plotbylbl : variable label `plotbyname'

		// check if plotvar is in matching set (otherwise it does not make much sense)
		if (ustrregexm("`e(match_set)'", "\b`plotbyname'\b") == 0) {
			noisily dis as error "Variable `varlist' not in matching set (`e(match_set)')."
			error 321
			exit
		}

		// set input from syntax and nopo returns
		// sample
		tempvar touse
		mark `touse' if e(sample)
		// depvar
		local depvar "`e(depvar)'"
		// strata
		local strata "`e(prefix)'_strata"
		// treatment indicator (fix to 0/1, relabel accordingly)
		tempvar treat
		gen `treat' = 1 if `e(ref)'
		replace `treat' = 0 if `treat' != 1 & !mi(`e(by)')
		local treatname `e(by)' // for renaming tempvar upon save
		local treatlbl : variable label `treatname'
		lab var `treat' `treatlbl'
		local vallbl : value label `treatname'
		if ("`vallbl'" != "") {
			label list `vallbl'
			levelsof `e(by)', local(bylvls)
			levelsof `e(by)' if `e(ref)', local(reflvl)
			foreach lvl in `bylvls' {
				dis "`lvl'"
				local lbl : label `vallbl' `lvl'
				dis "`lbl'"
				if (`lvl' == `reflvl') lab def _bylbl 1 "`lbl'", modify
					else lab def _bylbl 0 "`lbl'", modify
			}
			lab val `treat' _bylbl
		}
		// ref group switch
		if ("`e(swap)'" == "1") local bref = abs(`e(ref)' - 1) // reference group for returns
			else local bref = `e(ref)'
		// support
		local support "`e(prefix)'_matched"
		// matching weights
		local mweight "`e(prefix)'_weights"
		// aweight
		if ("`e(weight)'" != "") {
			local weight = `e(weight)'
			local weightexp "[aw = `e(weight)']"
		}
		else {
			tempvar weight
			gen `weight' = 1
			local weightexp "[aw = `weight']"
		}

		// save levels of plotby
		if ("`keepalllevels'" != "") levelsof `plotby', local(plotbylvls)
			else levelsof `plotby' if `touse', local(plotbylvls)

		// check if levels sensible
		local nplotbylvls : word count `plotbylvls'
		if (`"force"' == "" & `nplotbylvls' > 30) {
			noisily dis as error "`plotby' has more than 30 levels. Specify 'force' option to override."
			error 134
			exit
		}

		// relevel plotby: gen variable with no gaps
		// default: sort by depvar
		tempvar plotbyreleveled
		gen `plotbyreleveled' = .
		// gen means
		if ("`nosort'" == "") {	
			preserve
				collapse (mean) `depvar' `weightexp', by(`plotby')
				sort `depvar'
				drop `depvar'
				tempvar sorter
				if ("`descending'" != "") gen `sorter' = _n
					else gen `sorter' = _N - _n + 1
 				tempfile sorted
				save `sorted'
			restore
			merge m:1 `plotby' using `sorted', nogen
		}	
		// reorder and relabel
		local plotbyvallbl : value label `plotby'
		local s = 0
		foreach l in `plotbylvls' {
			if ("`nosort'" == "") levelsof `sorter' if `plotby' == `l', local(s)
				else local ++s
			replace `plotbyreleveled' = `s' if `plotby' == `l'
			if ("`plotbyvallbl'" != "") {
				local lbl : label `varlist' `l'
				lab def _releveledlbl `s' "`lbl'", modify
			}
		}
		lab val `plotbyreleveled' _releveledlbl
		if ("`nosort'" == "") {
			drop `sorter'
			lab var `plotbyreleveled' "`plotbylbl' releveled by mean `depvar' `descending'"
		} 
		else {
			lab var `plotbyreleveled' "`plotbylbl' releveled"
		}
		
		// build plot components
		preserve

			// calculate values by levels of plotby
			cap drop mdepvar_matched
			gen mdepvar_matched = .
			cap drop mdepvar_diff
			gen mdepvar_diff = .
			cap drop mdepvar_diff_weighted
			gen mdepvar_diff_weighted = .
			gen n_weighted = .

			foreach l in `plotbylvls' {
				// DA/DB
				forvalues t = 0/1 {
					sum `depvar' `weightexp' if `treat' == `t' & `support' == 1 & `touse'
					local mdepvar_matched = r(mean)
					sum `depvar' `weightexp' if `treat' == `t' & `touse'
					local wntotal = r(sum_w)
					replace mdepvar_diff = `depvar' - `mdepvar_matched' ///
						if `plotby' == `l' & `treat' == `t' & `support' == 0 & `touse'
					replace mdepvar_diff_weighted = mdepvar_diff * (`weight'/`wntotal') ///
						if `plotby' == `l' & `treat' == `t' & `support' == 0 & `touse'
					count if `treat' == `t' & `touse'
					replace n_weighted = `weight' * (r(N)/`wntotal') ///
						if `plotby' == `l' & `treat' == `t' & `support' == 0 & `touse'
				}
			}
			replace mdepvar_diff_weighted = -mdepvar_diff_weighted if `treat' == 1 // reverse for DB

			// check if values the same as in nopo table
			noisily dis "Component sum check:"
			noisily table `treat' if `support' == 0, stat(sum mdepvar_diff_weighted)

			// keep all levels for plot? 
			// useful if plotted comparisons do not have the same plotby levels due to missings
			if ("`keepallevels'" == "") keep if `touse' 

			// collapse
			collapse ///
				(mean) mdepvar_diff (sum) mdepvar_diff_weighted ///
				(sum) n_weighted (mean) `plotby' ///
				, by(`plotbyreleveled' `treat')
			replace mdepvar_diff_weighted = . if n_weighted == 0

			// N as text: get plot area and coordinates from data
			tostring n_weighted, gen(n_weighted_str) format(%9.0f) force
			sum mdepvar_diff
			if (abs(r(max)) > abs(r(min))) local mmax = abs(r(max)) * 1.75 // make room for obs text
				else local mmax = abs(r(min)) * 1.75
			sum mdepvar_diff_weighted
			if (abs(r(max)) > abs(r(min))) local wmmax = abs(r(max)) * 1.75 // make room for obs text
				else local wmmax = abs(r(min)) * 1.75
			if (`nplotbylvls'/10 < 1) local yrangemax = `nplotbylvls' + 1
				else local yrangemax = `nplotbylvls'/10 + `nplotbylvls'
			cap drop nx
			gen nx = `mmax' // x value for n counts (added as mlabel)
			local text `" text(`yrangemax' `mmax' "N unmatched" "(weighted)", place(sw) just(right) size(small) xaxis(2)) "'
			local ysize = `nplotbylvls'/10 + 5

			// set default plot options
			#delimit ;
			if (`"`twopts'"' == "") local twopts `"
				legend(order(
					3 "Contribution of unmatched to D (top x-axis)"
					1 "Category-specific mean of unmatched - overall mean of matched (bottom x-axis)"
					) rows(2) margin(zero)  region(style(none)))
				ylabel(1(1)`nplotbylvls', valuelabel grid angle(horizontal)) 
				yscale(range(`yrangemax' 1)) ytitle("")
				xscale(range(-`wmmax' `wmmax') axis(1))
				xscale(range(-`mmax' `mmax') axis(2)) xlab(, axis(2) grid)
				xtitle("Difference in means", axis(2) margin(0 0 0 3)) 
				subtitle(, bcolor("237 237 237") margin(1 1 1 1.5))
				scheme(s1mono) xsize(9) ysize(`ysize')
				"';
			if (`"`twoptsby'"' == "") local twoptsby `" 
				ixtitle note("") b1title("") graphregion(margin(zero)) 
				"';
			if (`"`twoptsbar'"' == "") local twoptsbar `"
				fcolor(gs10%50) lcolor(gs10) lp(solid) lw(0.2)
				"';
			if (`"`twoptsscatter'"' == "") local twoptsscatter `" 
				mcolor(black) xline(0, lcolor(black) lwidth(0.2)) xaxis(1) 
				xtitle("Contribution of unmatched to D", margin(0 0 3 3)) 
				"';
			#delimit cr

			// plot
			twoway ///
				(bar mdepvar_diff `plotbyreleveled' if n_weighted >= `nmin' ///
					, horizontal xaxis(2) `text' `twoptsbar') ///
				(scatter `plotbyreleveled' mdepvar_diff_weighted if n_weighted >= `nmin' ///
					, `twoptsscatter') ///
				(scatter `plotbyreleveled' nx ///
					, xaxis(2) mcolor(none) mlabel(n_weighted_str) mlabpos(9) mlabgap(0) ///
					msize(vtiny)) ///
				, by(`treat', `twoptsby') `twopts'
				
			// save plot data?
			if (`"`save'"' != "") {
				drop nx n_weighted_str
				// order, rename and label
				order `treat' `plotby' `plotbyreleveled'
				lab var `plotby' "`plotbylbl'"
				lab val `plotby' `plotbyvallbl'
				rename `plotby' `plotbyname'
				rename `plotbyreleveled' `plotbyname'_relevel
				lab var `treat' "`treatlbl'"
				rename `treat' `treatname'
				rename mdepvar_diff `depvar'_diff
				lab var `depvar'_diff "Difference mean unmatched - overall mean of matched"
				rename mdepvar_diff_weighted `depvar'_diff_weighted
				lab var `depvar'_diff_weighted "Contribution of unmatched to D"
				lab var n_weighted "N unmatched (weighted)"
				// save
				noisily save `save', replace
			}

		restore
	}

end


//
// postestimation: balance/desc table by group/matching status/weight
//