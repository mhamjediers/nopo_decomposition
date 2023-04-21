//
// postestimation: plot gap component over dv distribution
//

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
		    replace `quantile' = `quantile'_`i' if `by'==`i'
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
// postestimation: balance/desc table by group/matching status/weight
//