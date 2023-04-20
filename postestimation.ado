//
// estimation returns definition
//

// nopo has to assert 2 group levels
// nopo has to export treatment var
// nopo has to export treatment value
// nopo has to export normalized outcome variable
// nopo has to export tsupport variable 0/1

/* local depvar e(depvar) // points to normalized var if that was requested!
local support e(support)
local weight e(weight)
local mweight e(mweight)
local treat e(treat)
scalar treatval = e(treatval)
local strata e(strata) */


//
// postestimation: plot gap component over dv distribution
//

// plotting wrapper: 5 plots needed (one for each gap component)
cap program drop nopoplot
program define nopoplot
syntax [if] [aweight], ///
	[treat(varlist max=1)] ///
	[Quantiles(integer 100)] ///
	[QMIN(integer 1)] ///
	[QMAX(integer 100)] ///
	[Stat(string)] ///
	[RAWUMdiff] ///
	[REVsign] /// reverse sign
	[TWType(string)] ///
	[TWOpts(string asis)] ///
	[RELative]

    // check if prior command was nopo
   	/* if ("e(cmd)" != "nopo") {
		noisily dis as error "Previous command was not nopo."
		error 301
		exit
	} */

	// set input from syntax and nopo returns
	tempvar touse
    mark `touse' `if'
	//replace `touse' = 0 if !e(sample)

	local strata "strata"
	local treatval = 1
	//local local weightexp "[aw = e(weight)]"
	local support _supp
	local mweight _match
	local treatment t
	local depvar wage

	// set defaults
	local stat "mean"
	if ("`twtype'" == "") local twtype "line"
	
	/* if ("e(weight)" != "") {
		local weightexp "[aw = `e(weight)']"
	}
	else { */
		tempvar w1
		gen `w1' = 1
		local weightexp "[aw = `w1']"
	/* } */
	tempvar treat
	gen `treat' = 1 if `treatment' == `treatval'
	replace `treat' = 0 if `treatment' != `treatval' & !mi(`treatment')
	noisily tab `treat'

	// abort if quantiles > depvar groups
	qui levelsof `depvar', local(depvarlvls)
	local ndepvarlvls : word count `depvarlvls'
	if (`quantiles' > `ndepvarlvls') {
		noisily dis as error "Less groups in `depvar' than quantiles requested (`quantiles')."
		error 148
		exit
	}

	// options passthru
	local opts `"q(`quantiles') qmin(`qmin') qmax(`qmax') stat(`stat') `revsign' `relative'"'

    // create plot values for each component
    // D
	tempfile d
	nopo_gapdist `depvar' if `touse' `weightexp', by(`treat') `opts' save(`d')
	// DX (requires new var containing stratum specific wages for treatment group)
	tempfile dx
	tempvar depvar_strata_mt
	bys `strata' `support': egen `depvar_strata_mt' = mean(`depvar') ///
		if `support' & `treat' == `treatval'
	tempvar depvar_strata_m
	bys `strata' `support': egen `depvar_strata_m' = max(`depvar_strata_mt')
	replace `depvar_strata_m' = . if !`support'
	nopo_gapdist `depvar_strata_m' if `touse' & `support' `weightexp' ///
		, by(`treat') `opts' save(`dx')
	// D0
	tempfile d0
	nopo_gapdist `depvar' if `touse' & `support' [aw = `mweight'] ///
		, by(`treat') `opts' save(`d0')
	// DA
	tempfile da
	nopo_gapdist `depvar' if `touse' & `treat' == `treatval' `weightexp' ///
		, by(`support') `rawumdiff' `opts' save(`da')
	// DB
	tempfile db
	nopo_gapdist `depvar' if `touse' & `treat' != `treatval' `weightexp' ///
		, by(`support') `rawumdiff' `opts' save(`db')

	// plot
	preserve
		use "`d'", clear
		rename diff d
		foreach c in dx d0 da db {
			merge 1:1 q using "``c''", nogen
			rename diff `c'
		}
		twoway ///
			(line d q, lp(solid)) ///
			(line dx q, lp(dash)) ///
			(line d0 q, lp(shortdash)) ///
			(line da q, lp(dash_dot)) ///
			(line db q, lp(longdash_dot)), ///
			legend(order(1 "D" 2 "DX" 3 "D0" 4 "DA" 5 "DB") rows(1) span) ///
			yline(0) scheme(s1mono)
	restore

end

// get nopo decomposition component values over distribution of depvar
cap program drop nopo_gapdist 
program define nopo_gapdist 
syntax varname [if] [aweight], ///
	by(varlist max=1) ///
	[Quantiles(integer 100)] ///
	[QMIN(integer 1)] ///
	[QMAX(integer 100)] ///
	[Stat(string)] ///
	[RAWUMdiff] ///
	[revsign] /// reverse sign
	[TWType(string)] ///
	[TWOpts(string asis)] ///
	[RELative] ///
	[SAVE(string asis)]

quietly {
	preserve

		// subset
		if "`if'" != "" keep `if'
		
		// allow for quantile mean (or any other stat) estimation'
		// xtile aggregates quantiles if they contain constant values. Fill up to avoid empty cells.
		gen q = .
		lab var q "Quantile of group-specific wage distribution"
		forvalues i = 0/1 {
			// quantiles
		    xtile q_`i' = `varlist' if `by' == `i' [`weight'`exp'], nquantiles(`quantiles')
			bys `by' (q_`i'): replace q_`i' = ceil(_n * `quantiles'/_N) if `by' == `i'
		    replace q = q_`i' if `by'==`i'
		}
		
		// collapse, use sum of weights (passed via `exp') as N
		collapse ///
            (`stat') `stat'_depvar = `varlist' ///
            (sum) n = `=(subinstr("`exp'","=","",.))' ///
            if !mi(`varlist') [`weight'`exp'] ///
            , by(`by' q)
		reshape wide `stat'_depvar n, i(q) j(`by')
		
		// gen diff
		if ("`relative'" == "") {
			gen diff = `stat'_depvar1 - `stat'_depvar0
			if ("`rawumdiff'" == "") replace diff = diff * (n0/(n0+n1))
			lab var diff "Diff. in `stat' values per quantile"
		}
		else {
			gen diff = (1 - (`stat'_depvar0 / `stat'_depvar1)) * 100
			if ("`rawumdiff'" == "") replace diff = diff * (n0/(n0+n1))
			lab var diff "Relative diff. in `stat' values per quantile"
		}
		
		// reverse gap direction
		if ("`revsign'" != "") replace diff = -diff

		// save temp data
		if ("`save'" != "") {
			drop if !inrange(q, `qmin', `qmax') 		
			save `save'
		}
	restore
}

end


//
// postestimation: balance/desc table by group/matching status/weight
//