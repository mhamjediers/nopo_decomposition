cap cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"
run "nopo.ado"
kmatch md female educ_c exper_c tenure_c (lnwage), ate att atc po ifgenerate(t1_*) replace
nopo decomp, kmatchse
stop

// example data
nopo ex 1
eststo clear

// ATT

// nopo
nopo decomp lnwage educ_c exper_c tenure_c, by(female) km(md)
eststo att1
local msharewA = `e(msharewA)'/100
local msharewB = `e(msharewB)'/100
local nA = `e(nA)'
local nB = `e(nB)'
local mgapA = `e(mgapA)'
local mgapB = `e(mgapB)'
kmatch em female educ_c exper_c tenure_c (lnwage) [pw=wt], nate atc att generate wgenerate replace po ifgenerate
gen _IF_DX = _IF_Y0_ATT - _IF_Y0_ATC
total _IF_NATE _IF_ATT _IF_DX [pw=wt]
mat V = diag(vecdiag(e(V)))
mat list V
stop

mata: N = st_nobs()
mata: IF = st_data(.,"_IF_DX", .)
mata: w = st_data(., "wt")
mata: W = sum(w)
mata: sqrt(colsum(IF:^2)/(W - W/888) / W)'



/*

// nlcom
nlcom ///
    (D: _b[Y1(NATE)] - _b[Y0(NATE)]) ///
    (D0: _b[ATT]) ///
    (DX: _b[Y0(ATT)] - _b[Y0(ATC)]) ///
    (DA: (_b[Y0(ATC)] - (_b[Y0(NATE)] - (_b[Y0(ATC)] * `msharewA'))/(1-`msharewA')) * (1-`msharewA') ) ///
    (DB: (_b[Y1(ATT)] - (_b[Y1(NATE)] - (_b[Y1(ATT)] * `msharewB'))/(1-`msharewB')) * -1 * (1-`msharewB') ) ///
    , post
eststo att2

// ATC

// nopo
nopo decomp lnwage educ_c exper_c tenure_c, by(female) atc
eststo atc1
local msharewA = `e(msharewA)'/100
local msharewB = `e(msharewB)'/100
local nA = `e(nA)'
local nB = `e(nB)'
local mgapA = `e(mgapA)'
local mgapB = `e(mgapB)'
qui kmatch em female educ_c exper_c tenure_c (lnwage), nate atc att generate replace po coeflegend

// ATC
nlcom ///
    (D: _b[Y1(NATE)] - _b[Y0(NATE)]) ///
    (D0: _b[ATC]) ///
    (DX: _b[Y1(ATT)] - _b[Y1(ATC)]) ///
    (DA: (_b[Y0(ATC)] - (_b[Y0(NATE)] - (_b[Y0(ATC)] * `msharewA'))/(1-`msharewA')) * (1-`msharewA') ) ///
    (DB: (_b[Y1(ATT)] - (_b[Y1(NATE)] - (_b[Y1(ATT)] * `msharewB'))/(1-`msharewB')) * -1 * (1-`msharewB') ) ///
    , post
eststo atc2

// tab
esttab att* atc*, se nodepvars // ATC ATT are the same for nopose

*/
