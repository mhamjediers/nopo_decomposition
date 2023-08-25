clear

set seed 1234
local sim = 1000
local n = 1000
matrix est = J(1,8,.) // Empty matrix for simulation runs
qui forvalues i = 1/`sim' {
	
	set obs `n'
	gen d = runiform() > 0.8 // dummy indicator
    gen notd = 1 - d
	gen y = d + rnormal() // outcome

	mata: N = st_nobs()
	mata: D = st_data(.,"d")
    mata: NOTD = st_data(.,"notd")
	mata: Y = st_data(.,"y")

	mata: y0 = mean(Y, NOTD) // mean matched
	mata: IF_d0 = N/sum(NOTD) * NOTD :* (Y:-y0)
	
	mata: y1 = mean(Y, D) // mean unmatched
	mata: IF_d1 = N/sum(D) * D :* (Y:-y1)
	
	mata: diff = (y0 - y1) * sum(D) / N // mean difference multiplied by share of D == 1
	mata: IF_diff = (IF_d0 - IF_d1) * sum(D) / N
	
	mata: EST = y0, y1, diff, sqrt(diagonal(variance((IF_d0, IF_d1, IF_diff)) / N))' // storing results
	mata: st_matrix("EST", EST)
    
    // robust estimations
    sum d, meanonly
    scalar m = r(mean)
    
    // reg on d
    reg y d, robust
    estimates store rob
    nlcom (b: _b[d] * m * sqrt((`n'-1) / `n')), post // account for scaling issue; see p. 10
    mat EST = EST, sqrt(e(V)[1,1])
    
    // suest diff
    reg y if d == 0
    estimates store d0
    reg y if d == 1
    estimates store d1
    suest d0 d1
    nlcom ([d0_mean]_cons - [d1_mean]_cons) * m, post
    mat EST = EST, sqrt(e(V)[1,1])

	matrix est = est \ EST
	drop *
}

matrix colnames est = y0 y1 diff y0se y1se diffse robse suestse
svmat est, names(col)
collapse (sd) y0 y1 emp = diff (mean) y0se y1se mod_if = diffse mod_rob = robse mod_suest = suestse
list emp mod_if mod_rob mod_suest
