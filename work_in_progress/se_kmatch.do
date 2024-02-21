cap cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"
run "nopo.ado"
kmatch em female educ_c exper_c tenure_c (lnwage), atc ifgenerate(t1_*) replace
nopo decomp, kmatchse
stop

// example data
//nopo ex 1
eststo clear

// ATT

// nopo
nopo decomp lnwage educ_c exper_c tenure_c, by(female) km(md) kmatchse atc
stop
eststo att1
local msharewA = `e(msharewA)'/100
local msharewB = `e(msharewB)'/100
local nA = `e(nA)'
local nB = `e(nB)'
local mgapA = `e(mgapA)'
local mgapB = `e(mgapB)'
kmatch em female educ_c exper_c tenure_c (lnwage), nate atc att generate wgenerate replace po ifgenerate
ereturn list
gen _IF_DX = _IF_Y0_ATT - _IF_Y0_ATC
total(_IF_NATE _IF_ATT _IF_DX)
mean(_IF_NATE _IF_ATT _IF_DX)

local t = "_IF_ATT _IF_Y1_ATT _IF_Y0_ATT _IF_ATC _IF_Y1_ATC _IF_Y0_ATC _IF_NATE _IF_Y1_NATE _IF_Y0_NATE"
dis ustrregexm("`t'", "\b(\S+)att\b", 1)
dis ustrregexs(1)
dis usubinstr("`=ustrregexs(0)'", "`=strupper("att")'", "", .)

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

local t = ustrregexm("lol me string _wwefwefwfwefwe_ATC", `"\S+_`=strupper("atc")'"')
dis ustrregexs(0)
