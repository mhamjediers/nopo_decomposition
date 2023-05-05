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
syntax varname [if] [aweight], ///
	[NOSORT] /// do not sort by depvar
	[DESCending] /// sort descending (as opposed to ascending if nosort is not specified)
	[KEEPALLlevels] /// keep all levels of var (if cond. ignored)
	[FORCE] /// do not check for no. of levels in var
	[twopts(string asis)] ///
	[twoptsbar(string asis)] ///
	[twoptsscatter(string asis)] ///
	[twoptsby(string asis)] ///
	[SAVE(string asis)]

	//quietly {

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
			local weight = `e(weight)'
			local weightexp "[aw = `e(weight)']"
		}
		else {
			tempvar weight
			gen `weight' = 1
			local weightexp "[aw = `weight']"
		}
		// plotvar
		tempvar plotby
		clonevar `plotby' = `varlist'

		// save all levels of plotby
		if ("`keepalllevels'" != "") levelsof `plotby', local(plotbylvls)
			else levelsof `plotby' if `touse', local(plotbylvls)

		// check if levels sensible
		local nplotbylvls : word count `plotbylvls'
		if (`"force"' == "" & `nplotbylvls' > 30) {
			noisily dis as error "`plotby' has more than 30 levels. Specify 'force' option to override."
			error 134
			exit
		}

		// sort by dv?
		// in any case: gen variable with no gaps
		tempvar plotbysorted
		gen `plotbysorted' = .
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
		local vallbl : value label `plotby'
		local s = 0
		foreach l in `plotbylvls' {
			if ("`nosort'" == "") levelsof `sorter' if `plotby' == `l', local(s)
				else local ++s
			replace `plotbysorted' = `s' if `plotby' == `l'
			if ("`vallbl'" != "") {
				local lbl : label `varlist' `l'
				lab def _sortedlbl `s' "`lbl'", modify
			}
		}
		lab val `plotbysorted' _sortedlbl
		if ("`nosort'" == "") drop `sorter'

		preserve

			keep if `touse'

			// calculate values by levels of plotby
			cap drop mdepvar_matched
			gen mdepvar_matched = .
			cap drop mdepvar_diff
			gen mdepvar_diff = .
			cap drop mdepvar_diff_weighted
			gen mdepvar_diff_weighted = .
			gen n_weighted = .

			foreach l in `plotbylvls' {

				// DA
				sum `depvar' `weightexp' if `treat' == 0 & `support' == 1
				local w2 = r(mean)
				sum `depvar' `weightexp' if `treat' == 0
				local n2 = r(sum_w)
				replace mdepvar_matched = `w2' if `plotby' == `l' & `treat' == 0
				replace mdepvar_diff = wage - `w2' if `plotby' == `l' & `treat' == 0
				dis "`weight'"
				replace mdepvar_diff_weighted = mdepvar_diff * (`weight'/`n2') ///
					if `plotby' == `l' & `treat' == 0
				count if `treat' == 0
				replace n_weighted = `weight'*(r(N)/`n2') if `plotby' == `l' & `treat' == 0

				// DB
				sum `depvar' `weightexp' if `treat' == 1 & `support' == 1
				local w2 = r(mean)
				sum `depvar' `weightexp' if `treat' == 1
				local n2 = r(sum_w)
				replace mdepvar_matched = `w2' if `plotby' == `l' & `treat' == 1
				replace mdepvar_diff = wage - `w2' if `plotby' == `l' & `treat' == 1
				replace mdepvar_diff_weighted = -mdepvar_diff * (`weight'/`n2') ///
					if `plotby' == `l' & `treat' == 1
				count if `treat' == 1
				replace n_weighted = `weight'*(r(N)/`n2') if `plotby' == `l' & `treat' == 1

			}

			// check if values the same as in nopo table
			noisily table `treat' if `support' == 0, stat(sum mdepvar_diff_weighted)
			
			// collapse again
			collapse ///
				(mean) mdepvar_diff (sum) mdepvar_diff_weighted (sum) n_weighted if `support' == 0 ///
				, by(`plotbysorted' `treat')
			
			// rmake nice when weighted n < 1
			replace mdepvar_diff = 0 if n_weighted < 0.5
			replace mdepvar_diff_weighted = 0 if n_weighted < 0.5
			//replace n_weighted = 0 if n_weighted < 0.5
			tostring n_weighted, replace format(%9.0f) force
			
			// N as text: get coordinates from data
			sum mdepvar_diff
			local mmax = r(max)
			dis "`mmax'"
			sum mdepvar_diff_weighted
			local wmmax = r(max)
			dis "`wmmax'"
			if (`nplotbylvls'/10 < 1) local yrangemax = `nplotbylvls' + 1
				else `nplotbylvls'/10 + `nplotbylvls'
			gen mx = `mmax' * 1.5 // x value for n counts (added as mlabel)
			local text `" text(`yrangemax' `=(`mmax'*1.5)' "{bf:N unmatched}" "(weighted)", place(sw) xaxis(2)) "'
			local ysize = `nplotbylvls'/10 + 5
			// text(40.2 30 "(weighted)", place(w) xaxis(2))

			// set default plot options
			#delimit ;
			if (`"`twopts'"' == "") local twopts `"
				legend(order(
					1 "Category-specific mean of unmatched - overall mean of matched (bottom x-axis)"
					3 "Contribution of unmatched to D (top x-axis)"
					) rows(2) margin(zero)  region(style(none)))
				ylabel(1(1)`nplotbylvls', valuelabel grid angle(horizontal)) 
				yscale(range(`yrangemax' 1)) ytitle("")
				xscale(range(-`wmmax' `=(`wmmax'*1.5)') axis(1))
				xscale(range(-`mmax' `=(`mmax'*1.5)') axis(2))
				xtitle("Gap", axis(2)) subtitle(, bcolor("237 237 237") margin(1 1 1 1.5))
				scheme(s1mono) xsize(9) ysize(`ysize')
				"';
			#delimit cr
			if (`"`twoptsby'"' == "") local twoptsby `" ixtitle note("") b1title("") graphregion(margin(zero)) "'
			if (`"`twoptsbar'"' == "") local twoptsbar "fcolor(none) lcolor(gs8) lp(solid) lw(0.5)"
			if (`"`twoptsline'"' == "") local twoptsscatter `" mcolor(black) xline(0, lcolor(black) lwidth(0.2)) xaxis(1) xtitle("Contribution of unmatched to D") "'
			///if (`"`twoptsobs'"' == "") local twoptsobs 

			// plot
			twoway ///
				(bar mdepvar_diff `plotbysorted', horizontal xaxis(2) `text' `twoptsbar') ///
				(scatter `plotbysorted' mdepvar_diff_weighted, `twoptsscatter') ///
				(scatter `plotbysorted' mx ///
					, xaxis(2) mcolor(none) mlabel(n_weighted) mlabpos(9) mlabgap(0) msize(vtiny)) ///
				, by(`treat', `twoptsby') `twopts'
				
			// save data?

		restore
	//}

end


//
// postestimation: balance/desc table by group/matching status/weight
//