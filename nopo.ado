//
// wrapper
//

cap program drop nopo
program define nopo, eclass properties(svyb)
syntax [anything] [if] [in] [fweight pweight iweight] , ///
  [ 	/// standalone onlys
    by(varlist max=1) /// matching groups
    xref(string) /// set characteristics reference group like xref(group == 1)
    bref(string) /// set coefficient reference group like bref(group == 0)
    swap /// swap groups and reference vector
    NORMalize /// normalize by dividing by reference group mean
    KMatch(string) /// kmatch subcmd: md ps em
    KMOpts(string asis) /// pass on kmatch options
    KMPASSthru(string) /// pass on additional ereturns from kmatch to nopo decomp
    KMKEEPgen /// keep all generated variables
    KMNOISily /// show kmatch output
    dtable /// do not show estimates table (makes sense for bootstrap)
    naivese /// report naive SE from weighted reg and SUEST
    /// post
    att atc /// allow for these options to keep terminology consistent
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
  if (!inlist("`subcmd'", "decomp", "gapoverdist", "dadb", "summarize", "ex")) {
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
        - xref is B; bref is A
      - option swap:
        - group B = treat == 0
        - group A = treat == 1
        - xref is B; bref is A
      - manually specify reference groups (xref, bref)
        
      */
      // treatment value = group order
      if ("`swap'" == "") {
        local _tval : word 2 of `_bylvls'
        local _cval : word 1 of `_bylvls'
      }
      else {
        local _tval : word 1 of `_bylvls'
        local _cval : word 2 of `_bylvls'
      }
      // check reference vectors
      if ("`xref'" != "") {
        // get numeric value
        local _xref = stritrim("`xref'")
        local _xref = ustrregexra("`_xref'", "^\w+[\s]?[=]+[\s]?", "")
      }
      if ("`bref'" != "") {
        // get numeric value
        local _bref = stritrim("`bref'")
        local _bref = ustrregexra("`_bref'", "^\w+[\s]?[=]+[\s]?", "")
      }
      if ("`_xref'" != "" & "`_bref'" != "" & "`_xref'" == "`_bref'") {
        // check if xref and bref make sense
        dis as error "The x and the b reference cannot be the same."
        error 198
        exit
      }
      // set att/atc logic depending on reference
      if ("`_xref'" != "" | "`_bref'" != "") {
        if ("`_tval'" == "`_xref'" | "`_cval'" == "`_bref'") local _te = "att"
          else local _te = "atc"
        if ("`atc'" != "") dis as text "Option '`atc'' ignored due to manually specified reference."
        if ("`att'" != "") dis as text "Option '`att'' ignored due to manually specified reference."
      }
      else if ("`atc'" != "" & "`att'" != "") {
        dis as error "Specify either 'att' or 'atc' as options, not both (or use xref() or bref())"
        error 198
        exit
      }
      else {
        local _te `atc' `att' // only one of both possible; set
      }
      local att // unset (in case it was passed)
      local atc // unset (in case it was passed)
      // ATT IS DEFAULT (to correspond to our expositions)
      if ("`_te'" == "atc") local atc = "atc"
        else local att = "att"
      
      //
      // Run kmatch
      //
      
      // get input
      gettoken _depvar varlist : varlist
      if ("`kmatch'" == "") local kmatch = "em"
      if ("`weight'" != "") local _weightexp "[`weight'`exp']"
      if ("`kmnoisily'" != "") local kmnoisily = "noisily"

      // clean factor notation if exact matching
      // kmatch em treats everything as factor and so does nopo_summarize after kmatch em
      if ("`kmatch'" == "em") {
        local _varlist_nf = ""
        foreach _var in `varlist' {
          local _var = ustrregexra("`_var'", "^i.*\.", "")
          local _varlist_nf "`_varlist_nf' `_var'"
        }
        local varlist = ustrltrim("`_varlist_nf'")
      }

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
        // sum: always group A = Control = T == 0
        sum `_depvar' if `by' == `_cval' & `touse' `_sum_weightexp', meanonly
        qui gen _`_depvarabbrev'_norm = `_depvar' / r(mean) if `touse'
        local _depvarlbl : variable label `_depvar'			
        // set normalized var as _depvar!
        local _depvar = "_`_depvarabbrev'_norm"
        if ("`_depvarlbl'" != "") lab var `_depvar' "`_depvarlbl' (normalized)"
        dis as text "Normalized outcome generated: `_depvar'"
      }

      // SEs supposed to be bootstrapped, so default is to not to compute standard errors
      // if not specifically requested via naivese (probably wrong)
      if ("`naivese'" == "") local nose = "nose"
        else local nose
      
      // run
      quietly {	
        `kmnoisily' kmatch `kmatch' `by' `varlist' (`_depvar') `if' `in' `_weightexp' ///
          , tval(`_tval') `att' `atc' generate wgenerate replace `kmopts' `nose'
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
    , `_mweight' `atc' `att' `kmpassthru' `kmkeepgen' `dtable' `naivese' `options'

end


//
// Nopo (2008) style decomposition
//

cap program drop nopo_decomp
program define nopo_decomp, eclass

    syntax , ///
    mweight(namelist max=1) ///
    [ ///
      att ///
      atc ///
      kmpassthru(string) ///
      kmkeepgen ///
      dtable ///
      naivese ///
    ]
  
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
    local _te = strupper("`att'`atc'")
    
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
    }
    else {
      local _wtype = "aweight"
      local _wexp = "=1"
    }
    if ("`naivese'" != "") {
      if ("`e(vce)'" == "analytic") local vce = "robust"
        else local vce = "`e(vce)' `e(clustvar)'" // cluster only alternative
    }
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
    
    // determine A, B, and xref
    levelsof `_tvar', local(_tvarlvls)
    local _cval = strrtrim(strltrim(usubinstr("`_tvarlvls'", "`_tval'", "", .)))
    local _groupA = "`_tvar' == `_cval'"
    local _groupB = "`_tvar' == `_tval'"
    if ("`atc'" != "") {
      local xref = "`_groupA'"
      local bref = "`_groupB'"
    }
    else if ("`att'" != "") {
      local xref = "`_groupB'"
      local bref = "`_groupA'"
    }

    // save treatment effect (as D0)
    scalar _d0 = e(b)[1,1]
    
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
    
    /*
     suest is used to calculate SEs accounting for the covariance between components. 
     DX is estimated as delta: = D - D0 - DA -DB. suest requires some hacking around weighting 
     restrictions:

     - D0 always needs aweights (pweights not allowed in suest)
     - if user provides fweights, we need to tell Stata that the weight used in D0 estimates is
       fweight instead of the actual aweight
     - use vce only in suest estimation

     SEs are not correct and are not returned if not specifically requested via 'naivese';
     users have to use bootstrap

     SPEED UP: since bootstrapping is necessary for anybody not requesting naivese, you can estimate
     everything without standard errors.
    */
    
    // weighting helpers
    tempvar weight_cons
    gen `weight_cons' `_wexp'
    local _wexp_cons = "=`weight_cons'"
    if ("`_wtype'" == "pweight") local _wtype_cons = "aweight"
        else local _wtype_cons = "`_wtype'"

    // D
    mean `_depvar' [`_wtype_cons' `_wexp_cons'] if `sample', over(`treat')
    scalar _meanA = e(b)[1,1]
    scalar _meanB = e(b)[1,2]
    if ("`naivese'" == "") {
      mat b = J(1, 5, .)
      mat b[1,1] = _meanB - _meanA
    }
    else {
      reg `_depvar' i.`treat' [`_wtype_cons' `_wexp_cons'] if `sample'
      estimates store d
    }
    
    // DA
    sum `_depvar' [`_wtype_cons' `_wexp_cons'] if `treat' == 0 & `matched' == 0 & `sample'
    if ("`naivese'" == "") scalar _meanumA = r(mean)
    scalar _numA = r(N)
    scalar _numwA = r(sum_w)
    sum `_depvar' [`_wtype_cons' `_wexp_cons'] if `treat' == 0 & `sample'
    scalar _nA = r(N)
    scalar _nwA = r(sum_w)
    scalar _nmA = _nA - _numA
    scalar _nmwA = _nwA - _numwA
    scalar _mshareA = (_nmA / _nA) * 100
    scalar _msharewA = (_nmwA / _nwA) * 100
    if ("`naivese'" == "") {
      sum `_depvar' [`_wtype_cons' `_wexp_cons'] if `treat' == 0 & `matched' == 1 & `sample', meanonly
      scalar _meanmA = r(mean)
      scalar _mgapA = _meanmA - _meanumA
      mat b[1,4] = _mgapA * (_numwA / _nwA)
      if (b[1,4] == .) mat b[1,4] = 0
    }
    else {
      // check if no variation in Y among unmatched: SE estimation problem with suest
      // if suest SE is missing or infinitesimally small, estimate with manually plugged in constant
      // do not check for other groups: among matched (also D0 DX) always serious estimation problem
      reg `_depvar' i.`matched' [`_wtype_cons' `_wexp_cons'] if `treat' == 0 & `sample'
      estimates store da
      suest da
      if (e(V)["mean:_cons", "mean:_cons"] < 1e-10) {
        noisily dis "No variation in `_depvar' among unmatched in group A."
        if (e(b)[1, "mean:_cons"] < 1e-6) {
          reg `_depvar' i.`matched' [`_wtype_cons' `_wexp_cons'] ///
            if `treat' == 0 & `sample', nocons
        }
        else {
          tempvar _cons
          gen double `_cons' = e(b)[1, "mean:_cons"]
          reg `_depvar' i.`matched' `_cons' [`_wtype_cons' `_wexp_cons'] ///
            if `treat' == 0 & `sample', hascons
        }
        estimates store da
      }
      if (_numA > 0) scalar _mgapA = _b[1.`matched']
        else scalar _mgapA = .
    }

    // DB
    sum `_depvar' [`_wtype_cons' `_wexp_cons'] if `treat' == 1 & `matched' == 0 & `sample'
    if ("`naivese'" == "") scalar _meanumB = r(mean)
    scalar _numB = r(N)
    scalar _numwB = r(sum_w)
    sum `_depvar' [`_wtype_cons' `_wexp_cons'] if `treat' == 1 & `sample'
    scalar _nB = r(N)
    scalar _nwB = r(sum_w)
    scalar _nmB = _nB - _numB
    scalar _nmwB = _nwB - _numwB
    scalar _mshareB = (_nmB / _nB) * 100
    scalar _msharewB = (_nmwB / _nwB) * 100
    if ("`naivese'" == "") {
      sum `_depvar' [`_wtype_cons' `_wexp_cons'] if `treat' == 1 & `matched' == 1 & `sample', meanonly
      scalar _meanmB = r(mean)
      scalar _mgapB = _meanmB - _meanumB
      mat b[1,5] = -1 * _mgapB * (_numwB / _nwB)
      if (b[1,5] == .) mat b[1,5] = 0
    }
    else {
      // check if no variation in Y among unmatched: SE estimation problem with suest
      // if suest SE is missing or infinitesimally small, estimate with manually plugged in constant
      // do not check for other groups: among matched (also D0 DX) always serious estimation problem
      reg `_depvar' i.`matched' [`_wtype_cons' `_wexp_cons'] if `treat' == 1 & `sample'
      estimates store db
      suest db
      if (e(V)["mean:_cons", "mean:_cons"] < 1e-10) {
        noisily dis "No variation in `_depvar' among unmatched in group B."
        if (e(b)[1, "mean:_cons"] < 1e-6) {
          reg `_depvar' i.`matched' [`_wtype_cons' `_wexp_cons'] ///
            if `treat' == 1 & `sample', nocons
        }
        else {
          tempvar _cons
          gen double `_cons' = e(b)[1, "mean:_cons"]
          reg `_depvar' i.`matched' `_cons' [`_wtype_cons' `_wexp_cons'] ///
            if `treat' == 1 & `sample', hascons
        }
        estimates store db
      }
      if (_numB > 0) scalar _mgapB = _b[1.`matched']
        else scalar _mgapB = .
    }

    // change to matching weight
    replace `weight_cons' = `mweight'

    // D0 (always uses aweights and the matching weight returned by kmatch)
    if ("`naivese'" == "") {
      mat b[1,2] = _d0 // scalar fetched from kmatch ereturns
    }
    else {
      reg `_depvar' i.`treat' [aw `_wexp_cons'] if `matched' == 1 & `sample'
      ereturn local wtype = "`_wtype_cons'" // tell Stata we used the originally provided weight
      estimates store d0
    }

    // DX (always uses aweights and the matching weight returned by kmatch)
    if ("`naivese'" == "") {
      mat b[1,3] = b[1,1] - b[1,2] - b[1,4] - b[1,5]
    }
    else {
      /* if ("`att'" != "") local _xref = 1
        else local _xref = 0
      tempvar sub expanded
      expand 2 if `treat' != `_xref', gen(`expanded')
      gen `sub' = `_xref' if `treat' != `_xref'
      replace `sub' = abs(1 - `_xref') if `expanded' == 1
      replace `weight_cons' `_wexp' if `expanded' == 1
      noi reg `_depvar' i.`sub' [aw `_wexp_cons'] if `matched' == 1 & `sample'
      ereturn local wtype = "`_wtype_cons'" // tell Stata we used the originally provided weight
      estimates store dx
      drop if `expanded' == 1 */

      // change to standard weight
      replace `weight_cons' `_wexp'

      // suest & nlcom
      suest d d0 da db, vce(`vce') // analytic and cluster only!
      nlcom ///
        (D: [d_mean]1.`treat') ///
        (D0: [d0_mean]1.`treat') ///
        (DX: [d_mean]1.`treat' ///
            - [d0_mean]1.`treat' ///
            - [da_mean]1.`matched' * ( _numwA / _nwA ) ///
            - [db_mean]1.`matched' * -1 * ( _numwB / _nwB )) ///
        (DA: [da_mean]1.`matched' * ( _numwA / _nwA )) ///
        (DB: [db_mean]1.`matched' * -1 * ( _numwB / _nwB )) ///
        , post

      // suest & nlcom
      /* noi suest d d0 dx da db, vce(`vce') // analytic and cluster only!
      nlcom ///
        (D: [d_mean]1.`treat') ///
        (D0: [d0_mean]1.`treat') ///
        (DX: [dx_mean]1.`sub') ///
        (DA: [da_mean]1.`matched' * ( _numwA / _nwA )) ///
        (DB: [db_mean]1.`matched' * -1 * ( _numwB / _nwB )) ///
        , post */
    }

    // return
    if (`"`naivese'"' == "") {
      // default
      mat colnames b = D D0 DX DA DB
      ereturn post b, obs(`_Nsample') esample(`sample') depname(`_depvar')
    }
    else {
      mat b = e(b)
      mat V = e(V)
      ereturn post b V, obs(`_Nsample') esample(`sample') depname(`_depvar')
    }
    ereturn local cmd = "nopo"
    ereturn local subcmd = "`subcmd'"
    if ("`_wtype'" != "" & "`_wexp'" != "=1") {
      ereturn local wtype = "`_wtype'"
      ereturn local wexp = "`_wexp'"
    }
    ereturn local teffect = "`_te'"
    ereturn local tvar = "`_tvar'"
    ereturn local tval = "`_tval'"
    ereturn local cval = "`_cval'"
    ereturn local groupA = "`_groupA'"
    ereturn local groupB = "`_groupB'"
    ereturn local xref = "`xref'"
    ereturn local bref = "`bref'"
    ereturn local by = "`_tvar'"
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
    cap drop _nopo_strata
    cap drop _nopo_ps
    if ("`_kmatch_subcmd'" == "em") {
      gen _nopo_strata = `_strata'
      lab var _nopo_strata "Matching stratum"
      ereturn local strata = "_nopo_strata"
      ereturn scalar nstrata = _nstrata
      ereturn scalar nstrata_matched = _nmstrata
    }
    if ("`_kmatch_subcmd'" == "ps") {
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
    ereturn scalar mshareA = _mshareA // unweighted
    ereturn scalar msharewA = _msharewA // weighted
    ereturn scalar mgapA = _mgapA // raw diff by matching status
    ereturn scalar nB = _nB
    ereturn scalar mshareB = _mshareB // unweighted
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
  if ("`_groupA'" == "`xref'") {
    local _refA "(xref)"
    local _refB "(bref)"
  }
  else {
    local _refA "(bref)"
    local _refB "(xref)"
  }
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
    di as text "Propensity-score matching:" _col(42) "`_param'" _col(68) "= " _col(71) %08.3g `_paramval'
  }
  else if ("`_kmatch_subcmd'" == "md") {
    di as text "Multivariate-distance matching:" _col(42) "`_param'" _col(68) "= " _col(71) %08.3g `_paramval'
  }
  if ("`_ridge'" != "") di as text _col(42) "Ridge parameter:" _col(68) "= " _col(71) %08.3g `_ridge'
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
  di as text "A: " abbrev("`_tvar'", 7) " == `_cval' `_refA'"  _col(30) "{c |}" /*
    */ as result _col(32) %8.0f _nmA /*
    */ _col(45) %8.0f _numA /*
    */ _col(57) %8.0f _nA /*
    */ _col(71) %08.3g _meanA
  di as text _col(4) abbrev("`_groupAlbl'", 25) _col(30) "{c |}" /*
    */ as result _col(33) %7.1f _mshareA /*
    */ _col(46) %7.1f `=100-_mshareA'
  di as text "B: " abbrev("`_tvar'", 7) " == `_tval' `_refB'" _col(30) "{c |}" /*
    */ as result _col(32) %8.0f _nmB /*
    */ _col(45) %8.0f _numB /*
    */ _col(57) %8.0f _nB /*
    */ _col(71) %08.3g _meanB
  di as text _col(4) abbrev("`_groupBlbl'", 25) _col(30) "{c |}" /*
    */ as result _col(33) %7.1f _mshareB /*
    */ _col(46) %7.1f `=100-_mshareB'
  di as text "{hline 29}{c BT}{hline 48}"
  if ("`_wtype'" != "" & "`_wexp'" != "=1") di as text "Note: N and % are unweighted." 
  dis ""

  // display estimates
  if ("`dtable'" == "") {
    ereturn display
  }
  
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

{p 6 8 2}
{cmdab:rawum:diff} plots the raw difference in {depvar} between matched and unmatched units of groups {it:A} and {it:B} per quantile. By default, this difference is scaled by the number of unmatched units relative to all units in each group and quantile (see {browse "https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md":documentation on github}). Thus, the raw differences do not sum to {it:DA} and {it:DB} across quantiles.

*/

// gap over distribution plotting wrapper: 5 plots needed (one for each gap component)
cap program drop nopo_gapoverdist
program define nopo_gapoverdist, rclass
syntax [if] [in], /// if/in might produce misleading results; undocumented
  [NQuantiles(integer 100)] ///
  [RAWUMdiff] /// undocumented
  [recast(string)] ///
  [twopts(string asis)] ///
  [twoptsd(string asis)] ///
  [twoptsd0(string asis)] ///
  [twoptsdx(string asis)] ///
  [twoptsda(string asis)] ///
  [twoptsdb(string asis)] ///
  [xsize(real 0)] ///
  [ysize(real 0)] ///
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
    // x reference
    if ("`e(teffect)'" == "ATT") local _xref = 1
      else local _xref = 0
    // support
    local _support = e(matched)
    // matching weights
    local _mweight = e(mweight)
    // weights
    if ("`e(wtype)'" != "") {
      local _weightexp "[`e(wtype)'`e(wexp)']"
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
      , by(`treat') xref(`_xref') comp(dx) mweight(`_mweight') `opts' save(`dx')
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
        rspec(||&&&&|) cspec(& %3s | %14.3g & %14.3g | %19.0g & %18.0g &)
      noisily dis "Note:"
      noisily dis "- The component sum across well-populated quantiles should correspond to the"
      noisily dis "  component estimates."
      if (_M[2,3] < `nquantiles' | _M[3,3] < `nquantiles' | _M[4,3] < `nquantiles' | _M[5,3] < `nquantiles') {
        noisily dis "- There are less unique quantile values than quantiles requested which means"
        noisily dis "  that across some quantiles, the value of `_depvar' does not change for"
        noisily dis "  (one of) the groups compared to the estimate the component."
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
            local _dadblegend `"`_dadblegend' `_i' "`_lbl'""'
          }
          local ++_i
        }
        if ("`recast'" == "") local recast "line"
        if (`"`twopts'"' == "") local twopts `"legend(order(`_dadblegend') rows(1) span) yline(0) scheme(s1mono) ylab(, angle(horizontal)) xlab(, grid) ylab(, grid)"'
        if (`xsize' > 0) local twopts `"`twopts' xsize(`xsize')"'
        if (`ysize' > 0) local twopts `"`twopts' ysize(`ysize')"'
        if ("`recast'" == "line") {
          if (`"`twoptsd'"' == "") local twoptsd "lp(solid) lw(0.5)"
          if (`"`twoptsd0'"' == "") local twoptsd0 "lp(shortdash)"
          if (`"`twoptsdx'"' == "") local twoptsdx "lp(dash)"
          if (`"`twoptsda'"' == "") local twoptsda "lp(dash_dot)"
          if (`"`twoptsdb'"' == "") local twoptsdb "lp(longdash_dot)"
        }

        twoway ///
          (`recast' d q, `twoptsd') ///
          (`recast' d0 q, `twoptsd0') ///
          (`recast' dx q, `twoptsdx') ///
          (`recast' da q, `twoptsda') ///
          (`recast' db q, `twoptsdb') ///
          , `twopts'

      }
      
      // save if requested
      if (`"`save'"' != "") noisily save `save', replace

      // return plotoptions
      return local recast = `"`recast'"'
      return local twopts = `"`twopts'"'
      return local twoptsd = `"`twoptsd'"'
      return local twoptsd0 = `"`twoptsd0'"'
      return local twoptsdx = `"`twoptsdx'"'
      return local twoptsda = `"`twoptsda'"'
      return local twoptsdb = `"`twoptsdb'"'
      if (`xsize' == 0 | `ysize' == 0) {
        qui gr_setscheme
        if (`xsize' == 0) local xsize = "`.__SCHEME.graphsize.x'"
        if (`ysize' == 0) local ysize = "`.__SCHEME.graphsize.y'"
      }
      return scalar xsize = `xsize'
      return scalar ysize = `ysize'
    
    restore
  }

end

// get nopo decomposition component values over distribution of depvar
cap program drop nopo_gapdist
program define nopo_gapdist
syntax varname [if] [fweight pweight iweight], ///
  by(varlist max=1) ///
  [xref(integer 1)] ///
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
      keep if `by' != `xref' & !mi(`by')
      cap drop _expanded
      expand 2, gen(_expanded)
      replace `by' = `xref' if _expanded == 1
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

 {p 6 8 2}
 {cmdab:keepall:levels} /// keep all levels of var (if cond. ignored)
*/

cap program drop nopo_dadb
program define nopo_dadb, rclass
syntax varname [if] [in], /// if/in might produce misleading results; undocumented
  [NOSORT] /// do not sort by depvar
  [DESCending] /// sort descending (as opposed to ascending if nosort is not specified)
  [KEEPALLlevels] /// keep all levels of var (if cond. ignored); undocumented
  [FORCE] /// do not check for no. of levels in var
  [nmin(real 1)] /// minimum number of unmatched weighted obs per cat. be printed
  [twopts(string asis)] ///
  [twoptsbar(string asis)] ///
  [twoptsscatter(string asis)] ///
  [twoptsn(string asis)] ///
  [twoptsby(string asis)] ///
  [xsize(real 0)] ///
  [ysize(real 0)] ///
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
    // label for plot output; revert to original bylabel when saved
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
      local _weightexp "[`e(wtype)'`e(wexp)']"
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
        
        if (`ysize' == 0) local ysize = `_nplotbylvls'/5 + 5
        if (`xsize' == 0) {
          if (`ysize' < 8) local xsize = 9
            else local xsize = 9 + `ysize'/3
        }

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
          scheme(s1mono) xsize(`xsize') ysize(`ysize')
          "';
        if (`"`twoptsbar'"' == "") local twoptsbar `"
          horizontal xaxis(2) `_text' fcolor(gs10%50) lcolor(gs10) lp(solid) lw(0.2)
          "';
        if (`"`twoptsscatter'"' == "") local twoptsscatter `" 
          mcolor(black) xline(0, lcolor(black) lwidth(0.2)) xaxis(1) 
          xtitle("Contribution of unmatched to D", margin(0 0 3 3)) 
          "';
        if (`"`twoptsn'"' == "") local twoptsn `" 
          xaxis(2) mcolor(none) mlabel(n_weighted_str) mlabpos(9) mlabgap(0) msize(vtiny))
          "';
        if (`"`twoptsby'"' == "") local twoptsby `" 
          ixtitle note("") b1title("") graphregion(margin(zero)) 
          "';
        #delimit cr

        // plot
        twoway ///
          (bar mdepvar_diff `plotbyreleveled' if n_weighted >= `nmin', `twoptsbar') ///
          (scatter `plotbyreleveled' mdepvar_diff_weighted if n_weighted >= `nmin', `twoptsscatter') ///
          (scatter `plotbyreleveled' nx, `twoptsn' ///
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

      // return plotoptions
      return local twopts =  ustrtrim(ustrregexra(`"`twopts'"', "[\s\t]+", " "))
      return local twoptsbar =  ustrtrim(ustrregexra(`"`twoptsbar'"', "[\s\t]+", " "))
      return local twoptsscatter =  ustrtrim(ustrregexra(`"`twoptsscatter'"', "[\s\t]+", " "))
      return local twoptsn =  ustrtrim(ustrregexra(`"`twoptsn'"', "[\s\t]+", " "))
      return local twoptsby =  ustrtrim(ustrregexra(`"`twoptsby'"', "[\s\t]+", " "))
      return scalar xsize = `xsize'
      return scalar ysize = `ysize'

    restore
  }

end


//
// Summary table by group/matching status/weight
//

cap program drop nopo_summarize
program define nopo_summarize, rclass
syntax [varlist (default=none fv)] [if] [in], ///
  [STATistics(string)] /// mean mean/sd?
  [label] ///
  [labelwidth(real 1)] ///
  [fvdummies] ///
  [fvpercent] ///
  [keepempty]

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
    local _treatvallbl : value label `_treatname'
    if ("`_treatvallbl'" != "") {
      levelsof `_treatname', local(_bylvls)
      levelsof `_treatname' if `treat' == 1, local(_reflvl)
      foreach _lvl in `_bylvls' {
        local _lbl : label `_treatvallbl' `_lvl'
        if (`_lvl' == `_reflvl') lab def _bylbl 1 "`_lbl'", modify
          else lab def _bylbl 0 "`_lbl'", modify
      }
      lab val `treat' _bylbl
    }
    // x reference group
    if ("`e(teffect)'" == "ATT") local _xref = 1
      else local _xref = 0
    // support
    local _support = e(matched)
    // matching weights
    local _mweight = e(mweight)
    // weights
    if ("`e(wtype)'" != "") {
      local _weightexp "[`e(wtype)'`e(wexp)']"
      if ("`e(wtype)'" == "pweight") local _weightexp = usubinstr("`_weightexp'", "pweight", "aweight", .)
    }
    // matching type
    local _kmatch = e(kmatch_subcmd)

    //
    // Create table as matrix
    // 
    /*
     Some N handling is done in nopo decomp, so that there are matched obs in both groups.
     Here we only need to capture missing obs for the unmatched among A and B.

     Table labels come from matrix equation and col/row names.

     Approach is var x comparison, and the rows are stacked. Rather inelegant, but most flexible
     to accomodate that factors always need a mean estimation and need to be split into levels,
     which just have a single output per level: the share.

     We also save every statistic in its own table to ease manual postprocessing.

     [ToDo: let a separate program do the matrix processing, lots of repetetions below]

    */

    // statistics
    if ("`statistics'" == "") local statistics = "mean sd"
    local _nstats : word count `statistics' 

    // vars to tab; defaults to matching set
    if ("`varlist'" == "") local varlist "`_depvar' `e(matchset)'"

    // loop over variables and samples; table is built row by row
    local _i = 0
    foreach _var in `varlist' {
      local ++_i

      // which vars to tabstat
      local _tabstatvars = ""

      // rownames
      local _rownames = "" // gather for full table
      local _rownames_sep = "" // gather for separate table for each statistic
      
      // factor variable processing
      if ("`_var'" != "`_depvar'" & (ustrregexm("`_var'", "^i.*\.") == 1 | "`_kmatch'" == "em")) {
        
        // set factor indicator
        local _factor = 1

        // set statistics
        local _statistics = "mean"

        // expand factors
        local _var = ustrregexra("`_var'", "^i.*\.", "")
        levelsof `_var', local(_varlvls)
        if ("`_varlvls'" == "0 1") local _dummy = 1 // used to omit base level from table
          else local _dummy = 0
        local _j = 0
        foreach _lvl in `_varlvls' {
          local ++_j
          if (`_j' == 1 & `_dummy' == 1 & "`fvdummies'" == "") continue // skip dummy base level
          // gen one variable by factor lvl; gather in local
          gen _noposum_`_j' = `_var' == `_lvl' if !mi(`_var')
          local _tabstatvars = "`_tabstatvars' _noposum_`_j'"
          // gather rownames from value labels
          if ("`label'" != "") {
            local _lbl = ""
            local _vallbl : value label `_var'
            if ("`_vallbl'" != "") {
              local _lbl : label `_var' `_lvl'
              local _lbl = abbrev(ustrregexra("`_lbl'", "\.|:", ""), 32)
              }
            if ("`_lbl'" != "") local _rownames = `" `_rownames' "`_lbl'" "' 
              else local _rownames = `" `_rownames' "`_lvl'" "' // numval as fallback  
          }
          else {
            local _rownames = `" `_rownames' "`_lvl'" "' // numval as fallback  
          }
        }
      }
      else {
        // no factors
        local _factor = 0
        local _dummy = 0
        local _statistics = "`statistics'"
        local _tabstatvars = "`_var'"
        foreach _stat in `_statistics' {
          if (!inlist("`_stat'", "sd", "SD")) local _statlbl = strproper("`_stat'")
            else local _statlbl "SD"
          local _rownames = `" `_rownames' "`_statlbl'" "'
        }
      }

      // get var label if present
      local _varlbl = ""
      if ("`label'" != "") local _varlbl : variable label `_var'
      if ("`_varlbl'" == "") local _varlbl = "`_var'"
      local _varlblabbrev = abbrev(ustrregexra("`_varlbl'", "\.|:", ""), 32)

      // build rownames as equations -> varname:stat or varname:lvl
      // if dummy or for single stat tables, do not use equation labeling
      local _rownameseq = ""
      local _rownameseq_sep = ""
      local _j = 1 // label once per eq for separate stat tables
      foreach _rowname in `_rownames' {
        if (`_dummy' == 1 & "`fvdummies'" == "") local _rownameseq = `" `_rownameseq' "`_varlblabbrev'" "'
          else local _rownameseq = `" `_rownameseq' "`_varlblabbrev':`_rowname'" "'
        if (`_factor' == 1) {
          local _rownameseq_sep = `"`_rownameseq'"'
        }
        else if (`_j' == 1) {
          local _rownameseq_sep = `" `_rownameseq_sep' "`_varlblabbrev'" "'
        } 
        local ++_j
      }
      local _rownames = `"`_rownameseq'"'
      local _rownames_sep = `"`_rownameseq_sep'"'

      //
      // estimate mean for each sample and concatenate
      //
      
      // gather colnames
      local _colnames = ""

      // A_matched
      tabstat `_tabstatvars' if `treat' == 0 & `_support' == 1 & `touse' `_weightexp' ///
        , stat(`_statistics') save
      mat _S = r(StatTotal)
      if (`_factor' == 1) mat _S = _S'
      if (`_factor' == 1 & "`fvpercent'" != "") mat _S = _S * 100
      mat _V = _S // full table element
      // stat-specific table element 
      forvalues _s = 1/`_nstats' {
        if (`_factor' == 0) {
          mat _stat`_s' = _S[`_s', 1]
        }
        else if (`_factor' == 1 & `_s' == 1) {
          mat _stat`_s' = _S[`_s'..., 1]
        } 
        else {
          mat _stat`_s' = J(rowsof(_S), 1, .)
        }
      }
      local _colnames = `" `_colnames' "A_matched" "'

      // A_matched_weighted or B_matched_weighted
      if (`_xref' == 1) {
        // A -> A^B
        tabstat `_tabstatvars' if `treat' == 0 & `_support' == 1 & `touse' [aw = `_mweight'] ///
          , stat(`_statistics') save
        local _colnames = `" `_colnames' "A_matched_weighted" "'
      }
      else if (`_xref' == 0) {
        // B -> B^A
        tabstat `_tabstatvars' if `treat' == 1 & `_support' == 1  & `touse' [aw = `_mweight'] ///
          , stat(`_statistics') save
        local _colnames = `" `_colnames' "B_matched_weighted" "'
      }
      mat _S = r(StatTotal)
      if (`_factor' == 1) mat _S = _S'
      if (`_factor' == 1 & "`fvpercent'" != "") mat _S = _S * 100
      mat _V = _V, _S // full table element
      // stat-specific table element 
      forvalues _s = 1/`_nstats' {
        if (`_factor' == 0) {
          mat _stat`_s' = _stat`_s', _S[`_s', 1]
        }
        else if (`_factor' == 1 & `_s' == 1) {
          mat _stat`_s' = _stat`_s', _S[`_s'..., 1]
        } 
        else {
          mat _stat`_s' = _stat`_s', J(rowsof(_S), 1, .)
        }
      }

      // B_matched
      /* table `_support' if `treat' == 1 & `touse' `_weightexp', stat(mean `_tabstatvars') stat(sumw)
      table `_support' if `treat' == 1 & `touse' [pw`e(wexp)'], stat(mean `_tabstatvars') stat(sumw) */
      tabstat `_tabstatvars' if `treat' == 1 & `_support' == 1 & `touse' `_weightexp' ///
        , stat(`_statistics') save
      mat _S = r(StatTotal)
      if (`_factor' == 1) mat _S = _S'
      if (`_factor' == 1 & "`fvpercent'" != "") mat _S = _S * 100
      mat _V = _V, _S
      // stat-specific table element 
      forvalues _s = 1/`_nstats' {
        if (`_factor' == 0) {
          mat _stat`_s' = _stat`_s', _S[`_s', 1]
        }
        else if (`_factor' == 1 & `_s' == 1) {
          mat _stat`_s' = _stat`_s', _S[`_s'..., 1]
        } 
        else {
          mat _stat`_s' = _stat`_s', J(rowsof(_S), 1, .)
        }
      }
      local _colnames = `" `_colnames' "B_matched" "'

      // A_unmatched
      if (`e(mshareA)' < 100) {
        tabstat `_tabstatvars' if `treat' == 0 & `_support' == 0 & `touse' `_weightexp' ///
          , stat(`_statistics') save
        mat _S = r(StatTotal)
        if (`_factor' == 1) mat _S = _S'
        if (`_factor' == 1 & "`fvpercent'" != "") mat _S = _S * 100
        mat _V = _S, _V
        // stat-specific table element 
        forvalues _s = 1/`_nstats' {
          if (`_factor' == 0) {
            mat _stat`_s' = _S[`_s', 1], _stat`_s'
          }
          else if (`_factor' == 1 & `_s' == 1) {
            mat _stat`_s' = _S[`_s'..., 1], _stat`_s'
          } 
          else {
            mat _stat`_s' = J(rowsof(_S), 1, .), _stat`_s'
          }
        }
        local _colnames = `" "A_unmatched" `_colnames' "'
      }
      else if ("`keepempty'" != "") {
        local _colnames = `" "A_unmatched" `_colnames' "' // keep in table despite no unmatched obs
      } 

      // B_unmatched
      if (`e(mshareB)' < 100) {
        tabstat `_tabstatvars' if `treat' == 1 & `_support' == 0 & `touse' `_weightexp' ///
          , stat(`_statistics') save
        mat _S = r(StatTotal)
        if (`_factor' == 1) mat _S = _S'
        if (`_factor' == 1 & "`fvpercent'" != "") mat _S = _S * 100
        mat _V = _V, _S
        // stat-specific table element 
        forvalues _s = 1/`_nstats' {
          if (`_factor' == 0) {
            mat _stat`_s' = _stat`_s', _S[`_s', 1]
          }
          else if (`_factor' == 1 & `_s' == 1) {
            mat _stat`_s' = _stat`_s', _S[`_s'..., 1]
          } 
          else {
            mat _stat`_s' = _stat`_s', J(rowsof(_S), 1, .)
          }
        }
        local _colnames = `" `_colnames' "B_unmatched" "'
      }
      else if ("`keepempty'" != "") {
        local _colnames = `" `_colnames' "B_unmatched" "' // keep in table despite no unmatched obs
      } 

      // assign row names
      mat rownames _V = `_rownames'
      forvalues _s = 1/`_nstats' {
        mat rownames _stat`_s' = `_rownames_sep'
      }

      // concatenate
      if (`_i' == 1) {
        mat _M = _V
        forvalues _s = 1/`_nstats' {
          local _stat : word `_s' of `statistics'
          mat _`_stat' = _stat`_s'
        }
      }
      else {
        mat _M = _M \ _V
        forvalues _s = 1/`_nstats' {
          local _stat : word `_s' of `statistics'
          mat _`_stat' = _`_stat' \ _stat`_s'
        }
      }

      // drop generated variables; gen general indicator that factors have been used
      if ("`_factor'" == "1") {
        drop _noposum_*
        local _factor_present = 1
      }
    }

    // expand with empty columns if requested and no unmatched
    if (`e(mshareA)' == 100 & "`keepempty'" != "") {
      foreach _m in M `statistics' {
        local _rownames : rownames _`_m', quoted
        mat _`_m' = J(rowsof(_`_m'), 1, .), _`_m'
        mat rownames _`_m' = `_rownames'
      }
    }
    if (`e(mshareB)' == 100 & "`keepempty'" != "") {
      foreach _m in M `statistics' {
        mat _`_m' = _`_m', J(rowsof(_`_m'), 1, .)
      }
    }

    // label (always provide headings via equations)
    local _colnames = usubinstr(`"`_colnames'"', "A_", "A:", .)
    local _colnames = usubinstr(`"`_colnames'"', "B_", "B:", .)
    if ("`label'" != "") {
      if ("`_treatvallbl'" != "") {
        local _Albl : label _bylbl 0
        local _colnames = usubinstr(`"`_colnames'"', "A:", "`_Albl':", .)
        local _Blbl : label _bylbl 1
        local _colnames = usubinstr(`"`_colnames'"', "B:", "`_Blbl':", .)
      }
      local _colnames = usubinstr(`"`_colnames'"', "_weighted", " & weighted", .)
    }
    foreach _m in M `statistics' {
      mat colnames _`_m' = `_colnames'
    }
    
    // determine column format by no. of columns
    if ("`labelwidth'" == "") local labelwidth = 0
    if (colsof(_M) == 3) {
      if (`labelwidth' < 14) local _twidth = 14
        else local _twidth = `labelwidth'
      local _format = "%18.3g"	
    }
    else if (colsof(_M) == 4) {
      if (`labelwidth' < 13) local _twidth = 13
        else local _twidth = `labelwidth'
      local _format = "%13.3g"
    }
    else {
      if (`labelwidth' < 12) local _twidth = 12
        else local _twidth = `labelwidth'
      local _format = "%10.3g"
    }
    // list and return main table
    noisily matlist _M, lines(columns) showcoleq(combined) twidth(`_twidth') format(`_format')
    if ("`_factor_present'" == "1") {
      if ("`fvpercent'" != "") noisily dis "Note: Factor levels are printed as shares in %."
        else noisily dis "Note: Factor levels are printed as shares."
    }
    return mat table = _M
    // return separate tables for all statistics
    foreach _stat in `statistics' {
      return mat `_stat' = _`_stat'
    }

  }

end
