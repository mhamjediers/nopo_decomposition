cap cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"
run "nopo.ado"

// example data
nopo ex 1
eststo clear

// ATT

// nopo
nopo decomp lnwage educ exper tenure, by(female) km(md) kmatchse
ereturn list
eststo att1
local msharewA = `e(msharewA)'/100
local msharewB = `e(msharewB)'/100
local nA = `e(nA)'
local nB = `e(nB)'
local mgapA = `e(mgapA)'
local mgapB = `e(mgapB)'

// kmatch + nlcom
qui kmatch md female educ exper tenure (lnwage), nate atc att generate wgenerate replace po ifgenerate sharedbwidth
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
nopo decomp lnwage educ_c exper_c tenure_c [pw=wt], by(female) atc km(md) kmatchse
eststo atc1
local msharewA = `e(msharewA)'/100
local msharewB = `e(msharewB)'/100
local nA = `e(nA)'
local nB = `e(nB)'
local mgapA = `e(mgapA)'
local mgapB = `e(mgapB)'

// kmatch + nlcom
qui kmatch md female educ_c exper_c tenure_c (lnwage) [pw=wt], nate atc att generate replace po sharedbwidth
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


