//
// wrapper
//

cap program drop nopo
program define nopo, eclass properties(svyb)
syntax [anything] [if] [in] [fweight pweight iweight] , ///
	[ 	/// standalone onlys
		by(varlist max=1) /// matching groups
		swap /// swap groups and reference vector
		bref(string) /// set coefficient reference group like group == 0
		NORMalize /// normalize by dividing by reference group mean
		KMatch(string) /// kmatch subcmd: md ps em
		KMOpts(string asis) /// pass on kmatch options
		KMPASSthru(string) /// pass on additional ereturns from kmatch to nopo decomp
		KMKEEPgen /// keep all generated variables
		KMNOISily /// show kmatch output
		/// post
		att atc /// allow for these optiosn to keep terminology consistent
		* ///
	]

	/*
	 In essence, this wrapper does two things:

	 (1) Call a Nopo (2008) style decomposition after matching via kmatch. If a varlist is
	     specified, this wrapper calls a default version of kmatch, otherwise it checks if all 
		 requirements are met by the kmatch performed before nopo. Estimates are returned for
		 the decomposition components and a few auxiliary things.
	 (2) Call postestimation stuff:
	 	 - plot gap components over the outcome distribution
		 - plot contributions to DA/DB by variable level
		 - show summary table by group:
		   A_unmatched, A_matched, A_matched_weighted/B_matched_weighted, B_matched, B_unmatched
	*/

	// tokenize; determine decomp operation
	if ("`anything'" != "") gettoken subcmd varlist : anything
		else local subcmd "decomp"
	if (!inlist("`subcmd'", "decomp", "gapoverdist", "dadb", "summarize")) {
		dis as error "nopo subcommand must be one of:"
		dis as error "'decomp', 'gapoverdist', 'dadb', 'summarize'"
		error 198
		exit
	}

	// run kmatch or check if kmatch requirements met
	if ("`subcmd'" == "decomp" & "`varlist'" != "") {
		local nvars : word count `varlist'
		if (`nvars' == 1) {
			dis as error "'nopo decomp `varlist'' invalid. Valid options are:"
			dis as error "(1) Provide at least one matching variable in addition to depvar `varlist'."
			dis as error "(2) Use only 'nopo decomp' after running 'kmatch'"
			error 100
			exit
		}
		else if ("`by'" == "") {
			dis as error "Option 'by(groupvar)' required for decomposition"
			error 102
			exit
		}
		else {
			// no strings
			cap confirm numeric variable `by'
			if (_rc == 7) {
				dis as error "`by' must be numeric"
				error 7
				exit
			}
			// assert 2 levels
			qui levelsof `by', local(_bylvls)
			if (r(r) < 2) {
				dis as error "`by' must have 2 levels"
				error 148
				exit
			} 
			else if (r(r) > 2) {
				dis as error "`by' must have 2 levels"
				error 149
				exit
			}

			//
			// convert top-level syntax into ATT/ATC logic which is used in kmatch
			//
			/*
			- treatment indicator 0/1
				- group A = treat == 0
				- group B = treat == 1
				- ref is B
			- option swap:
				- group B = treat == 0
				- group A = treat == 1
				- ref is A
			- manually specify reference group for COEFFICIENTS bref(1)
			*/
			// treatment value = group order
			if ("`swap'" == "") local _tval : word 2 of `_bylvls'
				else local _tval : word 1 of `_bylvls'
			// att/atc logic depending on beta reference vector
			if ("`bref'" != "") {
				// get numeric value
				local _bref = stritrim("`bref'")
				local _bref = ustrregexra("`_bref'", "^\w+[\s]?[=]+\s", "")
				if (`_tval' == `_bref') local _te = "att"
					else local _te = "atc"
			}
			local att // unset (in case it was passed)
			local atc // unset (in case it was passed)
			if ("`_te'" == "atc") local atc = "atc"
				else local att = "att" // no _te defaults to att (correct for swap/non-swap)
			
			//
			// Run kmatch
			//
			
			// get input
			gettoken _depvar varlist : varlist
			if ("`kmatch'" == "") local kmatch = "em"
			if ("`weight'" != "") local _weightexp "[`weight'`exp']"
			if ("`kmnoisily'" != "") local kmnoisily = "noisily"
			
			// normalize depvar if requested: this results in a new variable!
			if ("`normalize'" != "") {
				local _depvarabbrev = abbrev("`_depvar'", 26)
				tempvar touse
				mark `touse' `if' `in'
				if ("`weight'" != "") {
					if ("`weight'" == "pweight") local _sum_weightexp = "[aw`exp']"
						else local _sum_weightexp = "[`weight'`exp']"
				}
				cap drop _`_depvarabbrev'_norm
				if ("`atc'" != "") {
					qui sum `_depvar' ///
						if `by' != `_tval' & !mi(`by') & `touse' `_sum_weightexp', meanonly
				}
				else {
					qui sum `_depvar' if `by' == `_tval' & `touse' `_sum_weightexp', meanonly
				}
				qui gen _`_depvarabbrev'_norm = `_depvar' / r(mean) if `touse'
				local _depvarlbl : variable label `_depvar'			
				// set normalized var as _depvar!
				local _depvar = "_`_depvarabbrev'_norm"
				if ("`_depvarlbl'" != "") lab var `_depvar' "`_depvarlbl' (normalized)"
				dis as text "Normalized outcome generated: `_depvar'"
			}
			
			// run
			quietly {	
				`kmnoisily' kmatch `kmatch' `by' `varlist' (`_depvar') `if' `in' `_weightexp' ///
					, tval(`_tval') `att' `atc' generate wgenerate replace `kmopts' 
					// perhaps strip opts of gen commands
			}
			
			// save matching weight variable for passthru
			local _mweight = "mweight(`e(wgenerate)')" // only one possible
			// clear varlist
			local varlist
		}
	}
	else if ("`subcmd'" == "decomp" & "`varlist'" == "") {
		// perform checks on prior kmatch
		if ("`e(cmd)'" != "kmatch") {
			dis as error "Previous command was not kmatch."
			error 301
			exit
		} 
		else if (!inlist("`e(subcmd)'", "md", "ps", "em")) {
			dis as error "nopo only works after kmatch md, kmatch ps, and kmatch em"
			error 301
			exit
		}
		else if ("`atc'" != "" & "`att'" != "") {
			dis as error "Specify either 'att' or 'atc' as options, not both"
			error 198
			exit
		}
		else if ("`e(att)'" == "" & "`att'" != "") {
			dis as error "No kmatch estimates found for ATT (use kmatch option 'att')"
			error 301
			exit
		}
		else if ("`e(atc)'" == "" & "`atc'" != "") {
			dis as error "No kmatch estimates found for ATC (use kmatch option 'atc')"
			error 301
			exit
		}
		else if ("`e(att)'" == "" & "`e(atc)'" == "") {
			dis as error "kmatch estimates required for either ATT or ATC (use kmatch options 'atc'/'att')"
			error 301
			exit
		}
		else if (e(N_over) > 1 | e(N_ovars) > 1) {
			dis as error "nopo requires kmatch to be specifcied with a single outcome (ovar) and without the 'over()' option"
			error 301
			exit
		}
		else if ("`e(generate)'" == "" | "`e(wgenerate)'" == "") {
			// rerun with matching variables if missing
			dis as text "nopo requires matching variables generated by kmatch. Re-running with options 'generate' and 'wgenerate'..."
			local _cmdline `" `e(cmdline)' "'
			if ("`e(generate)'" == "") local _cmdline `" `_cmdline' generate"' 
			if ("`e(wgenerate)'" == "") local _cmdline `" `_cmdline' wgenerate"'
			if (strpos(`"`_cmdline'"', " replace") == 0) local _cmdline `" `_cmdline' replace"'
			qui `_cmdline'
		}

		// default to ATT if present and not specified; fallback to ATC
		if ("`att'" == "" & "`atc'" == "" & "`e(att)'" != "") {
			local att = "att"
			local atc
		}
		else if ("`att'" == "" & "`atc'" == "" & "`e(atc)'" != "") {
			local att
			local atc = "atc"
		}
		
		// save matching weight variable for passthru (flexible approach to prior kmatch specs)
		local _idx = 1
		if ("`e(ate)'" != "") local ++_idx // if ate present, jump to next word in weight list
		if ("`e(att)'" != "" & "`atc'" != "") local ++_idx // defaults to ATT
		local _mweight : word `_idx' of `e(wgenerate)'
		local _mweight "mweight(`_mweight')"

	}
	// set passthru
	if ("`kmpassthru'" != "") local kmpassthru "kmpassthru(`kmpassthru')"

	// run subcommand with option passthru
	nopo_`subcmd' `varlist' ///
		, `_mweight' `atc' `att' `kmpassthru' `kmkeepgen' `options'

end


//
// Nopo (2008) style decomposition
//

cap program drop nopo_decomp
program define nopo_decomp, eclass

    syntax , ///
		mweight(namelist max=1) ///
		[att atc] ///
		[kmpassthru(string)] ///
		[kmkeepgen]

	quietly {
		
		//
		// use returns of kmatch for estimations or passthru
		//

		// log kmatch cmd
		local _kmatch_subcmd = e(subcmd)
		local _kmatch_cmdline = strrtrim(stritrim(`"`e(cmdline)'"'))

		// depvar
		local _depvar = e(depvar)
		
		// treatment (fixed to 0/1)
		local _tvar = e(tvar)
		local _tval = e(tval)
		tempvar treat
		gen `treat' = 0 if !mi(`_tvar')
		replace `treat' = 1 if `_tvar' == `_tval'
		
		// determine matching set from kmatch for return passthru; drop doublettes
		local _varset "`e(xvars)' `e(emvars)' `e(emxvars)'" // varnames = tokenizable as regex words
		local _nvarset : word count `_varset'
		while (`_nvarset' > 0) {
			gettoken _word _rest : _varset
			// save first occurence
			local _matchset "`_matchset' `_word'"
			// delete the rest
			local _varset = ustrregexra("`_rest'", "\b`_word'\b", "")
			// save new list
			local _nvarset : word count `_varset'
		}
		local _matchset = strrtrim(strltrim(stritrim("`_matchset'")))
		
		// weights
		if ("`e(wtype)'" != "") {
			local _wtype = e(wtype)
			local _wexp = e(wexp)
			local _weightexp = "[`_wtype'`_wexp']"
			if ("`_wtype'" == "pweight") local _sum_weightexp = "[aw`_wexp']"
				else local _sum_weightexp = `_weightexp'
			}
		if ("`e(vce)'" == "analytic") local vce // unset for default
			else local vce = e(vce)
		
		// generated matching vars processing
		/*
		 catch all for missing gen / wgen vars: 
		 - if manually deleted before calling nopo decomp
		 - if kmatch estimates are restored after nopo decomp without option 'kmkeepgen'
		*/
		cap desc `e(generate)'
		if (_rc == 111) {
			dis as error "kmatch variables `e(generate)' necessary for nopo decomp not found."
			dis as error "Use 'nopo decomp ..., kmkeepgen' if you want to subsequently restore the estimates of kmatch and run nopo decomp again."
			error 111
			exit
		}
		cap desc `e(wgenerate)'
		if (_rc == 111) {
			dis as error "kmatch variables `e(wgenerate)' necessary for nopo decomp not found."
			dis as error "Use 'nopo decomp ..., kmkeepgen' if you want to subsequently restore the estimates of kmatch and run nopo decomp again."
			error 111
			exit
		}
		tokenize `e(generate)'
		tempvar matched
		gen `matched' = 0 if `2' == 0 | `3' == 0
		replace `matched' = 1 if (`2' > 0 & !mi(`2')) | (`3' > 0 & !mi(`3'))
		if ("`_kmatch_subcmd'" == "ps") {
			local _strata
			local _ps = "`5'"
			if ("`kmkeepgen'" == "") drop `1' `2' `3' `4' `6' 
		}
		else {
			local _strata = "`5'"
			local _ps // unset
			if ("`kmkeepgen'" == "") drop `1' `2' `3' `4'
		}
		
		// obtaining number of strata and matched strata only for exact matching
		if ("`_kmatch_subcmd'" == "em") {
			mata: st_numscalar("_nstrata", colmax(st_data(., "`_strata'")))
			mata: st_numscalar("_nmstrata", length(uniqrows(st_data(., "`_strata'","`matched'"))))
		}
		
		// sample
		tempvar sample
		gen `sample' = e(sample)
		local _Nsample = e(N)
		mat Nsupport = e(_N)

		// abort if nobody has been matched
		count if `matched' == 1 & `sample'
		if (r(N) == 0) {
			dis as error "0 observations have been matched: Nopo decomposition not possible."
			error 2000
			exit
		}
		
		// determine A, B, and bref
		levelsof `_tvar', local(_tvarlvls)
		local _cval = strrtrim(strltrim(usubinstr("`_tvarlvls'", "`_tval'", "", .)))
		local _groupA = "`_tvar' == `_cval'"
		local _groupB = "`_tvar' == `_tval'"
		if ("`atc'" != "") local bref = "`_groupA'"
			else if ("`att'" != "") local bref = "`_groupB'"
		
		// save nn / kernel bandwidth / ridge param for display
		if ("`e(nn)'" != "") local _nn = e(nn)
		if ("`e(bwidth)'" != "") local _bwidth = e(bwidth)[1, "`att'`atc'"]
		if ("`e(ridge)'" != "") local _ridge = e(ridge)

		// abort if matched without replacement (replacement necessary for correct weights)
		if ("`e(wor)'" == "wor") {
			dis as error "nopo decomp requires matching with replacement. Option wor not allowed in:"
			dis as error " `_kmatch_cmdline'"
			error 322
			exit
		}

		
		// save everything from kmatch which has been requested for passthru
		/*
		 - We exclude everything we return and stuff that does not make sense for subcmds offered
		 - But: no idea what's possible with what, so just check everything and return if not empty
		 - There are a few returns with dynamic names, these are also not possible
		*/
		if ("`kmpassthru'" != "") {
			local _escalars ///
				k_omit N_clust N_outsup df_r nn_min nn_max pm_quantile pm_factor ///
				cv_factor maxiter btolerance
			local _emacros ///
				xvars ematch emxvars psvars pscore comsup generate wgenerate dygenerate ///
				idgenerate dxgenerate cemgenerate ifgenerate metric kernel keepall ///
				pscmd psopts pspredict bw_method cv_outcome cv_weighted cv_nopenalty ///
				cv_nolimit cv_exact ebalance ebvars csonly targets covariances nconstraint ///
				fitopts att atc vce clustvar title 
			local _ematrices ///
			 	_N S cv
			foreach _r in `kmpassthru' {
				if (ustrregexm("`_escalars'", "\b`_r'\b")) {
					scalar _`_r' = e(`_r') // save as local to be able to reference in loop later
					if (!mi(_`_r')) local _kmatch_escalars = "`_kmatch_escalars' _`_r'"
				}
				else if (ustrregexm("`_emacros'", "\b`_r'\b")) {
					local _`_r' = e(`_r')
					if ("`_`_r''" != "") local _kmatch_emacros = "`_kmatch_emacros' _`_r'"
				}
				else if (ustrregexm("`_ematrices'", "\b`_r'\b")) {
					cap mat _`_r' = e(`_r')
					if (_rc != 198) local _kmatch_ematrices = "`_kmatch_ematrices' _`_r'"
				}
			}
		}
		

		//
		// gather/estimate components
		//
		
		// placeholder matrices: D, D0, DA, DB
		mat b4 = J(1, 4, .)
		matname b4 D D0 DA DB, columns(1..4) explicit
		mat V4 = J(1, 4, .)
		matname V4 D D0 DA DB, columns(1..4) explicit
		
		// D0 (directly from kmatch)
		if ("`att'" != "") local _TE "ATT"
			else local _TE "ATC"
		mat b = e(b)[1, "`_TE'"]
		mat b4[1,2] = b[1,1]
		mat V = e(V)["`_TE'", "`_TE'"]
		mat V4[1,2] = V[1,1]
		
		// D
		mean `_depvar' if `sample' `_sum_weightexp', over(`treat')
		local _meanA = e(b)[1,1]
		local _meanB = e(b)[1,2]
		reg `_depvar' i.`treat' `_weightexp' if `sample', vce(`vce')
		mat b = e(b)
		mat b4[1,1] = b[1,2]
		mat V = e(V)
		mat V4[1,1] = V[2,2]
		
		// DA
		sum `_depvar' if `treat' == 0 & `matched' == 0 & `sample' `_sum_weightexp'
		scalar _numA = r(N)
		scalar _numwA = r(sum_w)
		sum `_depvar' if `treat' == 0 & `sample' `_sum_weightexp'
		scalar _nA = r(N)
		scalar _nwA = r(sum_w)
		reg `_depvar' i.`matched' if `treat' == 0 & `sample' `_weightexp', vce(`vce')
		scalar _mgapA = _b[1.`matched']
		nlcom _b[1.`matched'] * ( _numwA / _nwA ), post
		scalar _mshareuwA = (1 - _numA / _nA ) * 100
		scalar _msharewA = (1 - _numwA / _nwA ) * 100
		mat b = e(b)
		mat b4[1,3] = b[1,1]
		mat V = e(V)
		mat V4[1,3] = V[1,1]

		// DB
		sum `_depvar' if `treat' == 1 & `matched' == 0 & `sample' `_sum_weightexp'
		scalar _numB = r(N)
		scalar _numwB = r(sum_w)
		sum `_depvar' if `treat' == 1 & `sample' `_sum_weightexp'
		scalar _nB = r(N)
		scalar _nwB = r(sum_w)
		reg `_depvar' i.`matched' if `treat' == 1 & `sample' `_weightexp', vce(`vce')
		scalar _mgapB = _b[1.`matched']
		nlcom _b[1.`matched'] * -1 * ( _numwB / _nwB ), post
		scalar _mshareuwB = (1 - _numB / _nB ) * 100
		scalar _msharewB = (1 - _numwB / _nwB ) * 100
		mat b = e(b)
		mat b4[1,4] = b[1,1]
		mat V = e(V)
		mat V4[1,4] = V[1,1]
		
		// estimate DX from other components in the same model
		mat b5 = b4[1,1], b4[1,2], ., b4[1,3], b4[1,4]
		matname b5 D D0 DX DA DB, columns(1..5) explicit
		mat V5 = V4[1,1], V4[1,2], ., V4[1,3], V4[1,4]
		matname V5 D D0 DX DA DB, columns(1..5) explicit
		mat V4 = diag(V4)
		ereturn post b4 V4
		nlcom (_b[D] - _b[D0] - _b[DA] - _b[DB]), post
		mat b = e(b)
		mat b5[1,3] = b[1,1]
		mat V = e(V)
		mat V5[1,3] = V[1,1]
		mat V5 = diag(V5)

		// return
		ereturn post b5 V5, obs(`_Nsample') esample(`sample') depname(`_depvar')
		ereturn local cmd = "nopo"
		ereturn local subcmd = "`subcmd'"
		if ("`_wtype'" != "") {
			ereturn local wtype = "`_wtype'"
			ereturn local wexp = "`_wexp'"
		}
		ereturn local teffect = "`_TE'"
		ereturn local tvar = "`_tvar'"
		ereturn local tval = "`_tval'"
		ereturn local cval = "`_cval'"
		ereturn local groupA = "`_groupA'"
		ereturn local groupB = "`_groupB'"
		ereturn local bref = "`bref'"
		ereturn local matchset = strltrim("`_matchset'")
		// nopo vars (uses kmatch gen vars: weight & strata = copies = doublettes if kmkeepgen)
		// empty locals are not returned by Stata by default
		cap drop _nopo_matched
		rename `matched' _nopo_matched
		lab var _nopo_matched "Matching indicator (dummy)"
		ereturn local matched = "_nopo_matched"
		cap drop _nopo_mweight
		gen _nopo_mweight = `mweight'
		lab var _nopo_mweight "Matching weight"
		ereturn local mweight = "_nopo_mweight"
		if ("`_kmatch_subcmd'" == "em") {
			cap drop _nopo_strata
			gen _nopo_strata = `_strata'
			lab var _nopo_strata "Matching stratum"
			ereturn local strata = "_nopo_strata"
			ereturn scalar nstrata = _nstrata
			ereturn scalar nstrata_matched = _nmstrata
		}
		if ("`_kmatch_subcmd'" == "ps") {
			cap drop _nopo_ps
			gen _nopo_ps = `_ps'
			lab var _nopo_ps "Matching propensity score"
			ereturn local ps = "_nopo_ps"
		}
		if ("`kmkeepgen'" == "") drop `mweight' `_strata' `_ps'
	
		if ("`_nn'" != "") ereturn scalar nn = `_nn'
		if ("`_bwidth'" != "") ereturn scalar bwidth = `_bwidth'
		if ("`_ridge'" != "") ereturn scalar ridge = `_ridge'

		ereturn scalar N = `_Nsample'
		ereturn scalar nA = _nA
		ereturn scalar mshareuwA = _mshareuwA // unweighted
		ereturn scalar msharewA = _msharewA // weighted
		ereturn scalar mgapA = _mgapA // raw diff by matching status
		ereturn scalar nB = _nB
		ereturn scalar mshareuwB = _mshareuwB // unweighted
		ereturn scalar msharewB = _msharewB // weighted
		ereturn scalar mgapB = _mgapB // raw diff by matching status

		ereturn local kmatch_subcmd = "`_kmatch_subcmd'"
		ereturn local kmatch_cmdline = "`_kmatch_cmdline'"
		// return from passthru
		foreach _escalar in `_kmatch_escalars' {
			ereturn scalar kmatch`_escalar' = `_escalar'
		}
		foreach _emacro in `_kmatch_emacros' {
			ereturn local kmatch`_emacro' = "``_emacro''"
		}
		foreach _ematrix in `_kmatch_ematrices' {
			ereturn matrix kmatch`_ematrix' = `_ematrix'
		}
	}

	// display general info
	if ("`_groupA'" == "`bref'") local _refA "(ref)"
		else local _refB "(ref)"
	if ("`_nn'" != "") {
		local _param "NN requested:"
		local _paramval = `_nn'
	}
	else if ("`_bwidth'" != "") {
		// always set, irrespective if ridge matching or not
		local _param "Kernel bandwidth:"
		local _paramval = `_bwidth'
	}

	// get group labels
	local _treatvallbl : value label `_tvar'
	if ("`_treatvallbl'" != "") {
		local _groupAlbl : label `_treatvallbl' `_cval'
		local _groupBlbl : label `_treatvallbl' `_tval'
	}
	
	di as text " "
	di as text "Nopo decomposition" _col(42) "N" _col(68) "= " _col(71) %8.0g `_Nsample'
	if ("`_kmatch_subcmd'" == "em") {
		di as text "Exact matching:" _col(42) "N strata" _col(68) "= " _col(71) %8.0g _nstrata
		di as text _col(42) "N matched strata" _col(68) "= " _col(71) %8.0g _nmstrata
		di as text _col(42) "(unique combinations of matching set)"
	}
	else if ("`_kmatch_subcmd'" == "ps") {
		di as text "Propensity-score matching:" _col(42) "`_param'" _col(68) "= " _col(71) %05.3g `_paramval'
	}
	else if ("`_kmatch_subcmd'" == "md") {
		di as text "Multivariate-distance matching:" _col(42) "`_param'" _col(68) "= " _col(71) %05.3g `_paramval'
	}
	if ("`_ridge'" != "") di as text _col(42) "Ridge parameter:" _col(68) "= " _col(71) %05.3g `_ridge'
	dis ""
	di as text "{hline 29}{c TT}{hline 48}"
	di as text _col(30) "{c |}" /*
		*/ _col(45) "N / %" /*
		*/ _col(71) "Mean"
	di as text _col(30) "{c |}" _col(32) "{hline 33}" _col(67) "{hline 12}"
	di as text "Group " _col(30) "{c |}" /*
		*/ _col(33) "Matched" /*
		*/ _col(44) "Unmatched" /*
		*/ _col(60) "Total"  /*
		*/ _col(67) %12s abbrev("`_depvar'", 12)
	di as text "{hline 29}{c +}{hline 48}"
	di as text "A: " abbrev("`_groupAlbl'", 25)  _col(30) "{c |}" /*
		*/ as result _col(32) %8.0f `=_nA*_mshareuwA/100' /*
		*/ _col(45) %8.0f `=_nA*(1-_mshareuwA/100)' /*
		*/ _col(57) %8.0f _nA /*
		*/ _col(71) %08.3g `_meanA'
	di as text _col(4) abbrev("`_tvar'", 8) " == `_cval' `_refA'" _col(30) "{c |}" /*
		*/ as result _col(33) %7.1f _mshareuwA /*
		*/ _col(46) %7.1f `=100-_mshareuwA'
	di as text "B: " abbrev("`_groupBlbl'", 25) _col(30) "{c |}" /*
		*/ as result _col(32) %8.0f `=_nB*_mshareuwB/100' /*
		*/ _col(45) %8.0f `=_nB*(1-_mshareuwB/100)' /*
		*/ _col(57) %8.0f _nB /*
		*/ _col(71) %08.3g `_meanB'
	di as text _col(4) abbrev("`_tvar'", 8) " == `_tval' `_refB'" _col(30) "{c |}" /*
		*/ as result _col(33) %7.1f _mshareuwB /*
		*/ _col(46) %7.1f `=100-_mshareuwB'
	di as text "{hline 29}{c BT}{hline 48}"
	if ("`_wtype'" != "") di as text "Note: N and % are unweighted." 
	dis ""

	// display estimates
	ereturn display
	
end


//
// Plot gap component over dv distribution
//

/*

At the moment, the plot shows the gaps by comparing the outcome means in each quantile between 
groups. So, the mean across all these comparisons is the same as the decomposition component values
produced by nopo. But that also means that:

- At each quantile, the single component values do not add up to d
- D_A and D_B are scaled as in the nopo decomp, so higher values can mean more people or larger gaps
  (though the factor is always the same due to the quantile logic: n0/(n0+n1) is the same for each
  quantile). The `rawumdiff' option circumvents the scaling and shows the absolute differences,
  but then the values do not sum up to D_A/D_B

Does that sound sensible?

*/

// gap over distribution plotting wrapper: 5 plots needed (one for each gap component)
cap program drop nopo_gapoverdist
program define nopo_gapoverdist
syntax [if] [in], /// might produce strange results if if/in are used
	[NQuantiles(integer 100)] ///
	[RAWUMdiff] ///
	[twtype(string)] ///
	[twopts(string asis)] ///
	[twoptsd(string asis)] ///
	[twoptsd0(string asis)] ///
	[twoptsdx(string asis)] ///
	[twoptsda(string asis)] ///
	[twoptsdb(string asis)] ///
	[nodraw] ///
	[SAVE(string asis)]

	quietly {

		// check if prior command was nopo
		if ("`e(cmd)'" != "nopo") {
			noisily dis as error "Previous command was not nopo decomp"
			error 301
			exit
		}

		// set input from syntax and nopo returns
		// sample
		tempvar touse
		mark `touse' `if' `in'
		replace `touse' = 0 if !e(sample)
		// depvar
		local _depvar = e(depvar)
		// treatment indicator (fix to 0/1)
		tempvar treat
		gen `treat' = 1 if `e(tvar)' == `e(tval)'
		replace `treat' = 0 if `treat' != 1 & !mi(`e(tvar)')
		// b reference
		if ("`e(teffect)'" == "ATC") local _bref = 1
			else local _bref = 0
		// support
		local _support = e(matched)
		// matching weights
		local _mweight = e(mweight)
		// weights
		if ("`e(wtype)'" != "") {
			local _weightexp "[`e(wtype)' = `e(wexp)']"
		}
		else {
			tempvar w1
			gen `w1' = 1
			local _weightexp "[pw = `w1']"
		}

		// options passthru
		local opts `"nq(`nquantiles') qmin(`qmin') qmax(`qmax') `revsign' `relative'"'

		// create plot values for each component
		// D
		tempfile d
		noisily nopo_gapdist `_depvar' if `touse' `_weightexp', by(`treat') comp(d) `opts' save(`d')
		// D0
		tempfile d0
		noisily nopo_gapdist `_depvar' if `touse' & `_support' [pw = `_mweight'] ///
			, by(`treat') comp(d0) `opts' save(`d0')
		// DX (requires passing of matching weight)
		tempfile dx
		noisily nopo_gapdist `_depvar' if `touse' & `_support' `_weightexp' ///
			, by(`treat') bref(`_bref') comp(dx) mweight(`_mweight') `opts' save(`dx')
		// DA
		tempfile da
		noisily nopo_gapdist `_depvar' if `touse' & `treat' == 0 `_weightexp' ///
			, by(`_support') comp(da) `rawumdiff' `opts' save(`da')
		// DB
		tempfile db
		noisily nopo_gapdist `_depvar' if `touse' & `treat' == 1 `_weightexp' ///
			, by(`_support') comp(db) `rawumdiff' `opts' save(`db')

		// output
		preserve
			use `d', clear
			rename diff d
			foreach c in d0 dx da db {
				merge 1:1 q using "``c''", nogen
				rename diff `c'
				lab var `c' "Decomposition component `c'"
			}
			
			// summary table for sensibility checks
			foreach c in d d0 dx da db {
				tabstat `c' `c'_qcntmin `c'_nmin, save
				if ("`c'" == "d") mat _M = r(StatTotal)
					else mat _M = _M \ r(StatTotal)
			}
			mat _M = e(b)' , _M
			mat rownames _M = D D0 DX DA DB
			local colnames `" "Estimate" "Sum over q" "Minimum among compared groups: Unique q values" "Minimum among compared groups: N" "'
			mat colnames _M = `colnames'
			noisily dis "Component distribution across `nquantiles' quantiles of `_depvar' requested"
			noisily matlist _M, border(all) showcoleq(combined) ///
				rspec(||&&&&|) cspec(& %3s | %14.3g & %14.3g & %18.0g & %21.0g &)
			noisily dis "Note:"
			noisily dis "- The component sum across quantiles should correspond to the estimates with"
			noisily dis "  well populated quantiles."
			if (_M[2,3] < `nquantiles' | _M[3,3] < `nquantiles' | _M[4,3] < `nquantiles' | _M[5,3] < `nquantiles') {
				noisily dis "- There are less unique quantile values than quantiles requested which means"
				noisily dis "  that across some quantiles, the value of `_depvar' does not change for"
				noisily dis "  (one of) the groups compared to estimate the component."
			}
			if (inlist(., _M[2,2], _M[3,2], _M[4,2], _M[5,2])) {
				noisily dis "- No gap over the distribution could be estimated for the components where N"
				noisily dis "  of compared groups < no. of requested quantiles."
			}
			if (`nquantiles' == 100) noisily dis "- Use the nquantiles(#) option to set the number of quantiles."
			
			// plot
			if ("`nodraw'" == "") {
				
				// defaults
				local _i = 1 
				foreach _comp in d d0 dx da db {
					local _lbl = strupper("`_comp'")
					count if !mi(`_comp')
					if (r(N) > 0) {
						local _dadblegend `" `_dadblegend' `_i' "`_lbl'" "'
					}
					local ++_i
				}
				if ("`twtype'" == "") local twtype "line"
				if (`"`twopts'"' == "") local twopts `" legend(order(`_dadblegend') rows(1) span) yline(0) scheme(s1mono) ylab(, angle(horizontal)) xlab(, grid) ylab(, grid)"'
				if ("`twtype'" == "line") {
					if (`"`twoptsd'"' == "") local twoptsd "lp(solid) lw(0.5)"
					if (`"`twoptsd0'"' == "") local twoptsd0 "lp(shortdash)"
					if (`"`twoptsdx'"' == "") local twoptsdx "lp(dash)"
					if (`"`twoptsda'"' == "") local twoptsda "lp(dash_dot)"
					if (`"`twoptsdb'"' == "") local twoptsdb "lp(longdash_dot)"
				}

				twoway ///
					(`twtype' d q, `twoptsd') ///
					(`twtype' d0 q, `twoptsd0') ///
					(`twtype' dx q, `twoptsdx') ///
					(`twtype' da q, `twoptsda') ///
					(`twtype' db q, `twoptsdb') ///
					, `twopts'
			}
			// save if requested
			if (`"`save'"' != "") noisily save `save', replace
		restore
	}

end

// get nopo decomposition component values over distribution of depvar
cap program drop nopo_gapdist
program define nopo_gapdist
syntax varname [if] [fweight pweight iweight], ///
	by(varlist max=1) ///
	[bref(integer 1)] ///
	[mweight(varlist max=1)] ///
	[comp(string)] /// gap components, used as filter
	[NQuantiles(integer 100)] ///
	[RAWUMdiff] /// do not scale by N_A/N_B
	[SAVE(string asis)] * 

quietly {
	preserve
		
		// weight
		local weightvar = (subinstr("`exp'","=","",.))

		// subset
		if ("`if'" != "") keep `if'

		// dx: for by-logic, we just expand the data for `treat' == 0 and assign `treat' == 1
		if ("`comp'" == "dx") {
			keep if `by' == `bref'
			cap drop _expanded
			expand 2, gen(_expanded)
			replace `by' = abs(1 - `bref') if _expanded == 1
			replace `weightvar' = `mweight' if _expanded == 1 // replace weight with matching weight
		}
		
		// xtile aggregates quantiles if they contain constant values. Fill up to avoid empty cells.
		tempvar quantile
		gen `quantile' = .
		tempvar totweight
		local _qsuccess = 1 // use as trigger conditional on success in quantile estimation
		forvalues i = 0/1 {
			// quantiles
		    cap xtile `quantile'_`i' = `varlist' if `by' == `i' [`weight'`exp'], nquantiles(`nquantiles')
			// handle situation where requested quantiles are more than available obs per group
			if (_rc != 198) {
				// save number of quantiles (for info if they have been expanded)
				levelsof `quantile'_`i'
				if ("`_qcntmin'" == "") {
					local _qcntmin = r(r)
				}
				else {
					if (r(r) < `_qcntmin') local _qcntmin = r(r)
				}
				// expand and adjust weights
				bys `by': egen `totweight' = total(`weightvar')
				bys `by' (`quantile'_`i'): replace `quantile'_`i' = ///
					ceil(sum(`weightvar') * `nquantiles' / `totweight') if `by' == `i'
				replace `quantile' = `quantile'_`i' if `by' == `i'
				drop `totweight'
				
			} 
			else {
				local _qsuccess = 0
			}
			// save minimum N for check display
			count if `by' == `i'
			if ("`_nmin'" == "") {
				local _nmin = r(N)
			} 
			else {
				if (r(N) < `_nmin') local _nmin = r(N)
			}
		}
		
		// do only if estimation was successful; otherwise save empty data
		if (`_qsuccess' == 1) {
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
			if ("`comp'" == "db") {
				gen `diff' = `meanq'0 - `meanq'1
			}
			else {
				gen `diff' = `meanq'1 - `meanq'0
			}

			// scale D_A/D_B if not otherwise requested
			if (inlist("`comp'", "da", "db") & "`rawumdiff'" == "") {
				replace `diff' = `diff' * (`nq'0/(`nq'0+`nq'1))
			}

			keep `diff' `quantile'
			rename `diff' diff

			// add quantile count to inform about contant values
			gen `comp'_qcntmin = `_qcntmin'

		}
		else {
			keep if _n == 1
			gen diff = .
			keep diff `quantile'
			gen `comp'_qcntmin = .
		}
		
		// save temp data
		rename `quantile' q
		lab var q "Compared `varlist' quantile between groups (component-specific)"
		lab var diff "Gap"
		lab var `comp'_qcntmin "Min. number of different quantiles across compared groups"
		gen `comp'_nmin = `_nmin'
		lab var `comp'_nmin "Min. group N of compared groups"
		if ("`save'" != "")	save `save'

	restore
}

end


//
// Plot contribution to DA/DB by specified variable
//
/*
 Plot the contribution of each X level to DA/DB (weighted) and absolute difference in the outcome.
*/

cap program drop nopo_dadb
program define nopo_dadb
syntax varname [if] [in], ///
	[NOSORT] /// do not sort by depvar
	[DESCending] /// sort descending (as opposed to ascending if nosort is not specified)
	[KEEPALLlevels] /// keep all levels of var (if cond. ignored)
	[FORCE] /// do not check for no. of levels in var
	[nmin(real 1)] /// minimum number of unmatched weighted obs per cat. be printed
	[twopts(string asis)] ///
	[twoptsbar(string asis)] ///
	[twoptsscatter(string asis)] ///
	[twoptsby(string asis)] ///
	[nodraw] ///
	[SAVE(string asis)]

	quietly {
		
		// check if prior command was nopo
		if ("`e(cmd)'" != "nopo") {
			noisily dis as error "Previous command was not nopo decomp"
			error 301
			exit
		}

		// plotvar
		tempvar plotby
		clonevar `plotby' = `varlist'
		local _plotbyname `varlist' // for renaming tempvar upon save
		local _plotbylbl : variable label `_plotbyname'

		// check if plotvar is in matching set (otherwise it does not make much sense)
		if (ustrregexm("`e(matchset)'", "\b`_plotbyname'\b") == 0) {
			noisily dis as error "Variable `varlist' not in matching set (`e(matchset)')."
			error 321
			exit
		}
		
		// set input from syntax and nopo returns
		// sample
		tempvar touse
		mark `touse' `if' `in'
		replace `touse' = 0 if !e(sample)
		// depvar
		local _depvar = e(depvar)
		// treatment indicator (fix to 0/1)
		local _tval = e(tval)
		local _cval = e(cval)
		local _treatname = e(tvar) // for renaming tempvar upon save
		tempvar treat
		gen `treat' = 1 if `_treatname' == `_tval'
		replace `treat' = 0 if `_treatname' == `_cval'
		local _treatlbl : variable label `_treatname'
		if ("`treatlbl'" != "") lab var `treat' `_treatlbl'
		// label for plot putput; revert to original bylabel when saved
		local _treatvallbl : value label `_treatname'
		if ("`_treatvallbl'" != "") {
			label list `_treatvallbl'
			levelsof `_treatname', local(_bylvls)
			levelsof `_treatname' if `treat' == 1, local(_reflvl)
			foreach _lvl in `_bylvls' {
				local _lbl : label `_treatvallbl' `_lvl'
				if (`_lvl' == `_reflvl') lab def _bylbl 1 "`_lbl'", modify
					else lab def _bylbl 0 "`_lbl'", modify
			}
			lab val `treat' _bylbl
		}

		// support
		local _support = e(matched)
		// abort if complete support
		count if `_support' == 0
		if (r(N) == 0) {
			dis as error "0 unmatched observations: Plot not possible."
			error 2000
			exit
		}
		else {
			// or drop group which has full support
			count if `_support' == 0 & `treat' == 0
			if (r(N) == 0) replace `treat' = . if `treat' == 0
			count if `_support' == 0 & `treat' == 1
			if (r(N) == 0) replace `treat' = . if `treat' == 1
		}		
		
		// weights
		if ("`e(wtype)'" != "") {
			local _weightexp "[`e(wtype)' = `e(wexp)']"
			if ("`e(wtype)'" == "pweight") local _sum_weightexp = usubinstr("`_weightexp'", "pweight", "aweight", .)
				else local _sum_weightexp = `_weightexp'
		}
		else {
			tempvar weight
			gen `weight' = 1
			local _weightexp "[aw = `weight']"
		}

		// save levels of plotby
		if ("`keepalllevels'" != "") levelsof `plotby', local(_plotbylvls)
			else levelsof `plotby' if `touse', local(_plotbylvls)

		// check if levels sensible
		local _nplotbylvls : word count `_plotbylvls'
		if ("`force'" == "" & `_nplotbylvls' > 30) {
			noisily dis as error "`plotby' has more than 30 levels. Specify 'force' option to override"
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
				collapse (mean) `_depvar' `_weightexp', by(`plotby')
				sort `_depvar'
				drop `_depvar'
				tempvar sorter
				if ("`descending'" != "") gen `sorter' = _n
					else gen `sorter' = _N - _n + 1
 				tempfile sorted
				save `sorted'
			restore
			merge m:1 `plotby' using `sorted', nogen
		}	
		// reorder and relabel
		local _plotbyvallbl : value label `plotby'
		local _s = 0
		foreach _l in `_plotbylvls' {
			if ("`nosort'" == "") levelsof `sorter' if `plotby' == `_l', local(_s)
				else local ++_s
			replace `plotbyreleveled' = `_s' if `plotby' == `_l'
			if ("`_plotbyvallbl'" != "") {
				local _lbl : label `varlist' `_l'
				lab def _releveledlbl `_s' "`_lbl'", modify
			}
		}
		lab val `plotbyreleveled' _releveledlbl
		if ("`nosort'" == "") {
			drop `sorter'
			lab var `plotbyreleveled' "`_plotbylbl' releveled by mean `_depvar' `descending'"
		} 
		else {
			lab var `plotbyreleveled' "`_plotbylbl' releveled"
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

			foreach _l in `_plotbylvls' {
				// DA/DB
				forvalues _t = 0/1 {
					sum `_depvar' `_sum_weightexp' if `treat' == `_t' & `_support' == 1 & `touse'
					local _mdepvar_matched = r(mean)
					sum `_depvar' `_sum_weightexp' if `treat' == `_t' & `touse'
					local _wntotal = r(sum_w)
					replace mdepvar_diff = `_depvar' - `_mdepvar_matched' ///
						if `plotby' == `_l' & `treat' == `_t' & `_support' == 0 & `touse'
					replace mdepvar_diff_weighted = mdepvar_diff * (`weight'/`_wntotal') ///
						if `plotby' == `_l' & `treat' == `_t' & `_support' == 0 & `touse'
					count if `treat' == `_t' & `touse'
					replace n_weighted = `weight' * (r(N)/`_wntotal') ///
						if `plotby' == `_l' & `treat' == `_t' & `_support' == 0 & `touse'
				}
			}
			replace mdepvar_diff_weighted = -mdepvar_diff_weighted if `treat' == 0 // reverse for DA

			// check if values the same as in nopo table
			/* noisily dis "Component sum check:"
			noisily table `treat' if `_support' == 0, stat(sum mdepvar_diff_weighted) */

			// keep all levels for plot? 
			// useful if plotted comparisons do not have the same plotby levels due to missings
			if ("`keepalllevels'" == "") keep if `touse' 

			// collapse
			collapse ///
				(mean) mdepvar_diff (sum) mdepvar_diff_weighted ///
				(sum) n_weighted (mean) `plotby' ///
				, by(`plotbyreleveled' `treat')
			replace mdepvar_diff_weighted = . if n_weighted == 0

			if ("`nodraw'" == "") {

				// N as text: get plot area and coordinates from data
				sum mdepvar_diff if n_weighted >= `nmin'
				if (abs(r(max)) > abs(r(min))) local _mmax = abs(r(max)) * 1.75 // room obs text
					else local _mmax = abs(r(min)) * 1.8
				sum mdepvar_diff_weighted if n_weighted >= `nmin'
				if (abs(r(max)) > abs(r(min))) local _wmmax = abs(r(max)) * 1.75 // room obs text
					else local _wmmax = abs(r(min)) * 1.75
				if (`_nplotbylvls'/5 < 1) local _yrangemax = `_nplotbylvls' + 1
					else if (`_nplotbylvls'/5 < 2) local _yrangemax = `_nplotbylvls' + 2
					else local _yrangemax = `_nplotbylvls'/5 + `_nplotbylvls'
				local _ysize = `_nplotbylvls'/5 + 5
				if (`_ysize' < 8) local _xsize = 9
					else local _xsize = 9 + `_ysize'/3

				// N as text
				tostring n_weighted, gen(n_weighted_str) format(%9.0f) force
				cap drop nx
				gen nx = `_mmax' // x value for n counts (added as mlabel)
				local _text `" text(`_yrangemax' `_mmax' "N unmatched" "(weighted)", place(sw) just(right) size(small) xaxis(2)) "'

				// set default plot options
				#delimit ;
				if (`"`twopts'"' == "") local twopts `"
					legend(order(
						3 "Contribution of unmatched to D (top x-axis)"
						1 "Category-specific mean of unmatched - overall mean of matched (bottom x-axis)"
						) rows(2) margin(zero)  region(style(none)) size(small))
					ylabel(1(1)`_nplotbylvls', valuelabel grid angle(horizontal) labsize(small))
					yscale(range(`_yrangemax' 1)) ytitle("")
					xscale(range(-`_wmmax' `_wmmax') axis(1)) xlab(#5, axis(1) grid labsize(small))
					xscale(range(-`_mmax' `_mmax') axis(2)) xlab(#5, axis(2) grid labsize(small))
					xtitle("Difference in means", axis(2) margin(0 0 0 3)) 
					subtitle(, bcolor("237 237 237") margin(1 1 1 1.5))
					scheme(s1mono) xsize(`_xsize') ysize(`_ysize')
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
						, horizontal xaxis(2) `_text' `twoptsbar') ///
					(scatter `plotbyreleveled' mdepvar_diff_weighted if n_weighted >= `nmin' ///
						, `twoptsscatter') ///
					(scatter `plotbyreleveled' nx ///
						, xaxis(2) mcolor(none) mlabel(n_weighted_str) mlabpos(9) mlabgap(0) ///
						msize(vtiny)) ///
					, by(`treat', `twoptsby') `twopts'
			}

			// save plot data?
			if (`"`save'"' != "") {
				cap drop nx n_weighted_str
				// order, rename and label
				order `treat' `plotby' `plotbyreleveled'
				lab var `plotby' "`_plotbylbl'"
				lab val `plotby' `_plotbyvallbl'
				rename `plotby' `_plotbyname'
				rename `plotbyreleveled' `_plotbyname'_relevel
				lab var `treat' "`_treatlbl'"
				rename `treat' `_treatname'
				levelsof `_treatname', local(_treatlvls) // check if one group missing to label data nicely
				if (`r(r)' == 2) recode `_treatname' (1 = `_tval') (0 = `_cval')
					else if (`_treatlvls' == 0) recode `_treatname' (. = `_tval') (0 = `_cval')
					else if (`_treatlvls' == 1) recode `_treatname' (1 = `_tval') (. = `_cval')
				if ("`_treatvallbl'" != "") lab val `_treatname' `_treatvallbl'
				rename mdepvar_diff `_depvar'_diff
				lab var `_depvar'_diff "Difference mean unmatched - overall mean of matched"
				rename mdepvar_diff_weighted `_depvar'_diff_weighted
				lab var `_depvar'_diff_weighted "Contribution of unmatched to D"
				lab var n_weighted "N unmatched (weighted)"
				// save
				noisily save `save', replace
			}

		restore
	}

end


//
// Summary table by group/matching status/weight
//

cap program drop nopo_summarize
program define nopo_summarize, rclass
syntax [varlist (default=none)] [if] [in], ///
	[STATistics(string)] /// mean mean/sd?
	[label] ///
	[SAVE(string asis)]

	quietly {
		
		// check if prior command was nopo
		if ("`e(cmd)'" != "nopo") {
			noisily dis as error "Previous command was not nopo decomp"
			error 301
			exit
		}

		// set input from syntax and nopo returns
		// sample
		tempvar touse
		mark `touse' `if' `in'
		replace `touse' = 0 if !e(sample)
		// depvar
		local _depvar = e(depvar)
		// treatment indicator (fix to 0/1)
		tempvar treat
		gen `treat' = 1 if `e(tvar)' == `e(tval)'
		replace `treat' = 0 if `treat' != 1 & !mi(`e(tvar)')
		// assign correct labels
		local _treatname = e(tvar) // for renaming tempvar upon save
		local _vallbl : value label `_treatname'
		if ("`_vallbl'" != "") {
			label list `_vallbl'
			levelsof `_treatname', local(_bylvls)
			levelsof `_treatname' if `treat' == 1, local(_reflvl)
			foreach _lvl in `_bylvls' {
				local _lbl : label `_vallbl' `_lvl'
				if (`_lvl' == `_reflvl') lab def _bylbl 1 "`_lbl'", modify
					else lab def _bylbl 0 "`_lbl'", modify
			}
			lab val `treat' _bylbl
		}
		// reference group
		if ("`e(teffect)'" == "ATC") local _bref = 1
			else local _bref = 0
		// support
		local _support = e(matched)
		// matching weights
		local _mweight = e(mweight)
		// weights
		if ("`e(wtype)'" != "") {
			local _weightexp "[`e(wtype)' = `e(wexp)']"
			if ("`e(wtype)'" == "pweight") local _weightexp = usubinstr("`_weightexp'", "pweight", "aweight", .)
		}
		
		// vars to tab; defaults to matching set
		if ("`varlist'" == "") local varlist "`_depvar' `e(matchset)'"
		if ("`statistics'" == "") local statistics "mean sd"

		//
		// Create table as matrix
		// 
		/*
		 Some N handling is done in nopo decomp, so that there are matched obs in both groups.
		 Here we only need to capture missing obs for the unmatched among A and B.

		 Table labels come from matrix equation and col/row names.
		*/

		// A_matched
		tabstat `varlist' if `treat' == 0 & `_support' == 1 `_weightexp' & `touse' ///
			, stat(`statistics') save
		mat A_matched = r(StatTotal)
		nopo_stacktbl A_matched, `label'
		mat A_matched = r(A_matched)

		// A_unmatched
		cap tabstat `varlist' if `treat' == 0 & `_support' == 0 & `touse' `_weightexp' ///
			, stat(`statistics') save
		if (_rc == 2000) {
			local _umA = "mi"
		}
		else {
			mat A_unmatched = r(StatTotal)
			nopo_stacktbl A_unmatched, `label'
			mat A_unmatched = r(A_unmatched)
		}

		// A_matched_weighted or B_matched_weighted
		if (`_bref' == 0) {
			// A: b ref group for which treat == 0
			tabstat `varlist' if `treat' == 0 & `_support' == 1 & `touse' [aw = `_mweight'] ///
				, stat(`statistics') save
			local _matweighted = "A_matched_weighted"
		}
		else if (`_bref' == 1) {
			// B: b ref group for which treat == 1
			tabstat `varlist' if `treat' == 1 & `_support' == 1  & `touse' [aw = `_mweight'] ///
				, stat(`statistics') save
			local _matweighted = "B_matched_weighted"
		}
		mat `_matweighted' = r(StatTotal)
		nopo_stacktbl `_matweighted', `label'
		mat `_matweighted' = r(`_matweighted')

		// B_matched
		tabstat `varlist' if `treat' == 1 & `_support' == 1 `_weightexp' & `touse' ///
			, stat(`statistics') save
		mat B_matched = r(StatTotal)
		nopo_stacktbl B_matched, `label'
		mat B_matched = r(B_matched)

		// B_unmatched
		cap tabstat `varlist' if `treat' == 1 & `_support' == 0 `_weightexp' & `touse' ///
			, stat(`statistics') save
		if (_rc == 2000) {
			local _umB = "mi"
		}
		else {
			mat B_unmatched = r(StatTotal)
			nopo_stacktbl B_unmatched, `label'
			mat B_unmatched = r(B_unmatched)
		}

		// combine
		if ("`_umA'" != "mi") mat _M = A_unmatched, A_matched, `_matweighted', B_matched
			else mat _M = A_matched, `_matweighted', B_matched
		if ("`_umB'" != "mi") mat _M = _M, B_unmatched

		// label (always provide headings via equations)
		local _colnames : colnames _M, quoted
		local _colnames = usubinstr(`"`_colnames'"', "A_", "A:", .)
		local _colnames = usubinstr(`"`_colnames'"', "B_", "B:", .)
		if ("`label'" != "") {
			if ("`_vallbl'" != "") {
				local _Albl : label _bylbl 0
				local _colnames = usubinstr(`"`_colnames'"', "A:", "`_Albl':", .)
				local _Blbl : label _bylbl 1
				local _colnames = usubinstr(`"`_colnames'"', "B:", "`_Blbl':", .)
			}
			local _colnames = usubinstr(`"`_colnames'"', "_weighted", " & weighted", .)
		}
		mat colnames _M = `_colnames'

		// determine column format by no. of columns
		if (colsof(_M) == 3) {
			local _twidth = 14
			local _format = "%18.3g"	
		}
		else if (colsof(_M) == 4) {
			local _twidth = 13
			local _format = "%13.3g"
		}
		else {
			local _twidth = 12
			local _format = "%10.3g"
		}
		// list and return
		noisily matlist _M, lines(columns) showcoleq(combined) twidth(`_twidth') format(`_format')
		return mat npsum = _M

	}

end

// stack tabstat results for multiple statistics into a single column
cap program drop nopo_stacktbl
program define nopo_stacktbl, rclass
syntax namelist (max=1), ///
	[label]

	mat _IN = `namelist'
	local _nrows : rowsof(_IN)
	local _rownames : rownames _IN
	local _ncols : colsof(_IN)
	local _colnames : colnames _IN
	mat _OUT = J(`_nrows' * `_ncols', 1, .)
    
	local _cidx = 0
	foreach _colname in `_colnames' {
        local ++_cidx
        local _ridx = 0
        foreach _rowname in `_rownames' {
			local ++_ridx
			// gather later row names
			if ("`label'" == "") {
				local _stackednames = `"`_stackednames' "`_colname':`_rowname'""'
			}
			else {
				local _lbl : variable label `_colname'
				if ("`_lbl'" == "") local _lbl = "`_colname'"
				local _stackednames = `"`_stackednames' "`_lbl':`_rowname'""'
			}
			// replace values in placeholder matrix
			local _nidx = `_nrows' * (`_cidx' - 1) + `_ridx'
			mat _OUT[`_nidx', 1] = _IN[`_ridx', `_cidx']
		}
	}
	mat rownames _OUT = `_stackednames'
	mat colnames _OUT = `namelist'
    
	// return
	return mat `namelist' = _OUT
	
end