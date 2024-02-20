cap program drop nopo
run "/home/max/Seafile/Projects/nopo_decomposition/nopo.ado"

// example data
//nopo ex 1

// nopo
nopo decomp lnwage educ_c exper_c tenure_c, by(female)
eststo m1
local msharewA = `e(msharewA)'/100
local msharewB = `e(msharewB)'/100
local nA = `e(nA)'
local nB = `e(nB)'
local mgapA = `e(mgapA)'
local mgapB = `e(mgapB)'
qui kmatch em female educ_c exper_c tenure_c (lnwage), nate atc att generate replace po coeflegend

// mean lnwage, over(female _nopo_matched)
// mean lnwage, over(female)
nlcom ///
    (D: _b[Y1(NATE)] - _b[Y0(NATE)]) ///
    (D0: _b[ATT]) ///
    (DX: _b[Y0(ATT)] - _b[Y0(ATC)]) ///
    (DA: (_b[Y0(ATC)] - (_b[Y0(NATE)] - (_b[Y0(ATC)] * `msharewA'))/(1-`msharewA')) * (1-`msharewA') ) ///
    (DB: (_b[Y1(ATT)] - (_b[Y1(NATE)] - (_b[Y1(ATT)] * `msharewB'))/(1-`msharewB')) * -1 * (1-`msharewB') ) ///
    , post
eststo m2

esttab m1 m2, se
