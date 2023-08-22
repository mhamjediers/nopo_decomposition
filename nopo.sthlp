{smcl}
{* *! version 0.2   22aug2023  Maximilian Sprengholz & Maik Hamjediers}{...}
{vieweralsosee "kmatch" "kmatch"}{...}
{vieweralsosee "nopomatch" "nopomatch"}{...}
{viewerjumpto "Syntax" "nopo##syntax"}{...}
{viewerjumpto "Postestimation" "nopo##postest"}{...}
{viewerjumpto "Description" "nopo##desc"}{...}
{viewerjumpto "Options" "nopo##opts"}{...}
{viewerjumpto "Examples" "nopo##ex"}{...}
{title:Title}

{phang}{hi:nopo} {hline 2} Matching-based decomposition analysis of differences in outcomes between two groups following {c N~}opo (2008)

{marker:dep}{...}
{title:Dependencies}

{phang} Requires {help kmatch:{it:kmatch}} (Jann, 2017) to be installed, which has
{help moremata:{it:moremata}} as dependency.

{marker:syntax}{...}
{title:Syntax}

    Standalone (calls {help kmatch:{it:kmatch}} internally):

        {cmd:nopo decomp} {depvar} {varlist} {ifin} {weight} {cmd:, by(}varname{cmd:)} [{help nopo##comopt:{it:options}}]
   
    After calling {help kmatch:{it:kmatch}} manually:
   
        {cmd:nopo decomp}
    {pstd} Three different {help nopo##postest:{it:postestimation}} commands are available:
        {cmd:nopo summarize}, {cmd:nopo gapoverdist}, {cmd:nopo dadb} (see below).
 

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
    Standalone
{synopt :{cmdab:by(}{help varname:{it:varname}}{cmdab:)}}Group variable by which to estimate and decompose gaps in {depvar:} (required) {p_end}{...}
{synopt :{cmdab:xref(}{varname}{it: == #}{cmdab:)}}Set reference group in terms of {it:characteristics}{p_end}{...}
{synopt :{cmdab:bref(}{varname}{it: == #}{cmdab:)}}Set reference group in terms of {it:returns}{p_end}{...}
{synopt :{cmdab:swap}}Swap groups and reference{p_end}{...}
{synopt :{cmdab:norm:alize}}Normalize estimates by mean outcome of reference group{p_end}{...}
{synopt :{cmdab:km:atch(}em|md|ps{cmdab:)}}Choose between exact (em, default), multivariate-distance (md), and propensity score (ps) matching{p_end}{...}
{synopt :{cmdab:kmo:pts(}{it:string}{cmdab:)}}Options for internal {help kmatch:{it:kmatch}} call{p_end}{...}
{synopt :{cmdab:kmnois:ily}}Show output of internal kmatch call{p_end}{...}

    Both standalone and after manual {help kmatch:{it:kmatch}}
{synopt :{cmdab:kmpass:thru(}{it:string}{cmdab:)}}kmatch returns to be passed through and returned by {cmd:nopo decomp}{p_end}{...}
{synopt :{cmdab:kmkeep:gen}}Keep matching and weight variables generated by kmatch{p_end}{...}
{synoptline}{phang}{it:fweights}, {it:pweights}, and {it:iweights} are allowed. Use the {cmd:bootstrap:} prefix for bootstrapping (see {help nopo##desc:Description} for a note on standard errors).
{p2colreset}

{marker desc}{...}
{title:Description}

{pstd}
{hi:nopo decomp} provides a {c N~}opo-style (2008) decomposition of the gap {it:D} in the average outcome {it:Y} between two groups {it:A} and {it:B} by matching them on a set of characteristics {it:X} predictive of {it: Y}.

{p 6 6 2}{it:D = YB - YA}{p_end}
{p 8 8 2}{it:  = D0 + DX + DA + DB  , where}{p_end}

{p 4 6 2}
{it:- D0} is the part of the gap not attributable to compositional differences between {it:A} and {it:B} in {it:X} among matched units (classic {it:unexplained} component){p_end}
{p 4 6 2}
{it:- DX} is the part of the gap attributable to compositional differences between {it:A} and {it:B} in {it:X} among matched units (classic {it:explained} component){p_end}
{p 4 6 2}
{it:- DA} is the part of the gap attributable to unmatched units in {it:A}{p_end}
{p 4 6 2}
{it:- DB} is the part of the gap attributable to unmatched units in {it:B}{p_end}

{pstd}    
For a detailed explanation on methodology, please consult the {browse "https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md":online documentation} [...].

{pstd}{ul:Matching approaches:} 

{pstd}{c N~}opo's original proposition used exact matching but extends to other matching approaches, two of which are additionally implemented in {cmd:nopo decomp}: multivariate-distance and propensity-score matching.

{pstd}{ul:Standard errors:} 

{pstd}Currently, the standard errors [...]. Use bootstrapping to obtain empirical standard errors.

{pstd}{ul:Generated variables:}

{p2colset 5 20 24 4}{...}
{p2col:{bf:_nopo_matched}} Matching indicator (dummy){p_end}
{p2col:{bf:_nopo_mweight}} Matching weight{p_end}
{p2col:{bf:_nopo_strata}} Matching stratum ({it:em} only){p_end}
{p2col:{bf:_nopo_ps}} Matching propensity score ({it:ps} only){p_end}
{p2col:{bf:_{depvar}_norm}} Normalized {depvar} (if {cmd:normalize} was specified){p_end}


{marker opts}{...}
{title:Options}

{dlgtab:Standalone}

{phang}
{cmdab:by(}{help varname:{it:varname}}{cmdab:)} specifies the groups between which we estimate and 
decompose the gap in {depvar:} (required). Needs to be numeric with two levels. By default, the 
first {cmd:by} group is group {it:A}, which is the reference group in terms of {it:returns}
and the second group is group {it:B}, which is the reference group in terms of {it:characteristics} 
(see {browse "https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md":online documentation} for details). Use {cmd:xref()}/{cmd:bref()} or {cmd:swap} to change. 

{phang}
{cmdab:xref(}{varname}{it: == #}{cmd:)} allows to manually set the reference in terms of {it:characteristics}. 
See {browse "https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md":online documentation} 
for a detailed explanation what that means in the matching approach. Naturally, {cmd:xref()} and 
{cmd:bref()} cannot be the same.

{phang}
{cmdab:bref(}{varname}{it: == #}{cmd:)} allows to manually set the reference in terms of {it:returns}. 
See {browse "https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md":online documentation} 
for a detailed explanation what that means in the matching approach. Naturally, {cmd:bref()} and 
{cmd:xref()} cannot be the same.

{phang}
{cmdab:swap} groups {it:A} and {it:B}, so that the sign of {it:D} is reversed and and the respective reference for characteristics and returns is switched.

{phang}
{cmdab:norm:alize} estimates by the {depvar} mean of group {it: A}. Coefficients can then be interpreted in a relative manner, e.g. that group {it:B} earns on average 19 percent lower wages compared to group {it:A}. Generates {bf:_{depvar}_norm}.

{phang}
{cmdab:km:atch(}em|md|ps{cmdab:)} lets you choose the matching approach. The default is to use 
{it:exact matching} {cmd:kmatch(em)}, in which case all variables in {varlist} are treated as 
factors. For multivariate-distance {cmd:kmatch(md)} and propensity score {cmd:kmatch(ps)} 
matching, make sure to indicate via factor notation which variables are factors and which are
continuous (everything is passed through to the internal kmatch call as is). To further tweak the
matching, you can use {cmd:kmopts()} or, for maximum flexibility, call {cmd:nopo decomp} after a
manual {cmd:kmatch} call with all the necessary options.

{phang}
{cmdab:kmo:pts(}string{cmdab:)} are options passed on to {help kmatch:{it:kmatch}} when 
called internally. See {help kmatch##goptions:{it:general_options}}, 
{help kmatch##matchoptions:{it:matching_options}} and the matching-type specific options
{help kmatch##emoptions:{it:em_options}}, 
{help kmatch##mdoptions:{it:md_options}}, or {help kmatch##psoptions:{it:ps_options}}.

{phang}
{cmdab:kmnois:ily} Show output of internal kmatch call.

{dlgtab:Both standalone and after manual kmatch}

{phang}
{cmdab:kmpass:thru(}string{cmdab:)} lets you pass through further
{help kmatch##eret:{it:kmatch returns}} to the returns of {cmd:nopo decomp}. For example:
{cmd:kmpassthru(df_r metric)}.

{phang}
{cmdab:kmkeep:gen} keeps the matching and weight variables generated by
{help kmatch##gen:{it:generated by kmatch}}, which are dropped by default. In the standalone call,
these variables are prefixed by {bf:_KM_}). Note that some of the kmatch variables contain 
the same information as the variables returned by {cmd:nopo decomp} (prefixed by {bf:_nopo_}).


{marker postest}{...}
{title:Postestimation}
  
{dlgtab:nopo summarize}

{p 4 4 2}
{cmd:nopo summarize} [{varlist}] [{cmd:, label} {cmdab:stat:istics(}{help tabstat##statname:{it:statnames}}{cmd:)}]

{p 6 6 2}
Returns a descriptive table with selected statistics by group and matching/weighting status:
{it:A unmatched}, {it:A matched}, {it:A/B matched and weighted} (depends on {cmd:xref()} if {it:A} 
or {it:B}), {it:B matched}, and {it:B unmatched}. If no variables are specified, {depvar} and 
{varlist} from {cmd: nopo decomp} are summarized. By default, means and standard deviations are
reported. For factor variables, the shares of factor levels are reported in percent (indicate factors
by  using factor notation in {varlist} of either {cmd: nopo summarize} or the prior {cmd:nopo decomp}).

{p 6 8 2}
{cmd:label} displays labels instead of names for variables and value labels instead of values for 
factor variables.

{p 6 8 2}
{cmdab:stat:istics(}{help tabstat##statname:{it:statnames}}{cmd:)} are {it:mean sd} by default, but 
you can choose any other available {help tabstat##statname:{it:statnames}}.

{dlgtab:nopo gapoverdist}

{p 4 4 2}
{cmd:nopo gapoverdist} [{cmd:,} {help nopo##comopt:{it:options}}]

{p 6 6 2}
Plots decomposition-components over the distribution of {depvar} by comparing {depvar} {it:means} 
at each quantile for the respective groups compared (e.g. the matched and unmatched of group {it:B}
for {it:DB}). So, at each quantile, the single component values {it:do not} add up to {it:D}. But
the quantile values of each component sum to the decomposition component values produced by
{cmd:nopo decomp}.

{p 6 8 2}	   
{cmdab:nq:uantiles(}#{cmd:)} specifies the number of quantiles to be plotted. Default is 100.

{p 6 8 2}
{cmd:twtype(}{it:string}{cmd:)} sets the {help twoway:{it:twoway plottype}}. Default is {it: line}.

{p 6 8 2}
{cmd:twopts(}{it:string}{cmd:)}, {cmd:twoptsd(}{it:string}{cmd:)}, {cmd:twoptsd0(}{it:string}{cmd:)},
{cmd:twoptsdx(}{it:string}{cmd:)}, {cmd:twoptsda(}{it:string}{cmd:)}, and {cmd:twoptsdb(}{it:string}{cmd:)}
set the general and subplot-specific {help twoway_options:{it:twoway_options}}. Check the returns of
{cmd:nopo gapoverdist} to see defaults (e.g. {cmd:r(twopts)}). You can also use the {cmd:save()} option 
to build the plot from the underlying data yourself.

{p 6 8 2}
{cmd:xsize(}{it:real}{cmd:)} and {cmd:ysize(}{it:real}{cmd:)} set the x and y dimensions of the plot.

{p 6 8 2}
{cmd:nodraw} specifies that the plot is not produced at all. Useful in combination with {cmd:save()}
to build your own plot from the data.

{p 6 8 2}
{cmd:save(}{help filename:{it:filename}}{cmd:)} saves the underlying plot data.

{dlgtab:nopo dadb}

{p 4 4 2}
{cmd:nopo dadb} {varname} [{cmd:,} {help nopo##comopt:{it:options}}]

{p 6 6 2} 
Creates a plot showing how the different levels of {varname} contribute to {it:DA} and {it:DB}
because these levels are associated with either many unmatched units and/or large differences in 
{depvar} by matching status within groups {it:A} and {it:B}
(see {browse "https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md":online documentation})
for details. This {it:is not} the same as a detailed decomposition in regression-based approaches (which is generally not possible with matching).

{p 6 8 2}
{cmd:nosort} specifies that the plot is not to be sorted by the mean values of {depvar} per level of
{varname}. Default is to sort in ascending order.

{p 6 8 2}
{cmdab:desc:ending} specifies that the plot is to be sorted by the mean values of {depvar} per level of
{varname} in descending order. Default is to sort in ascending order.

{p 6 8 2}
{cmd:force} skips the check for the number of levels of {varname}. By default, the program exits if
{varname} has more than 30 levels.

{p 6 8 2}
{cmd:nmin(}{it:#}{cmd:)} specifies the minimum number of unmatched, weighted observations per level
of {varname} to be plotted. Default is 1. This option is useful if, for example, if a small number of 
unmatched observations strongly deviate from the average value of {depvar} in that group, producing
large bars.

{p 6 8 2}
{cmd:twopts(}{it:string}{cmd:)}, {cmd:twoptsbar(}{it:string}{cmd:)}, 
{cmd:twoptsscatter(}{it:string}{cmd:)}, {cmd:twoptsn(}{it:string}{cmd:)}, and
{cmd:twoptsby(}{it:string}{cmd:)} set the general, element-specific, and by-specific 
{help twoway_options:{it:twoway_options}}. Check the returns of {cmd:nopo dadb} to see defaults 
(e.g. {cmd:r(twopts)}). You can also use the {cmd:save()} option to build the plot from the 
underlying data yourself.

{p 6 8 2}
{cmd:xsize(}{it:real}{cmd:)} and {cmd:ysize(}{it:real}{cmd:)} set the x and y dimensions of the plot.

{p 6 8 2}
{cmd:nodraw} specifies that the plot is not produced at all. Useful in combination with {cmd:save()}
to build your own plot from the data.

{p 6 8 2}
{cmd:save(}{help filename:{it:filename}}{cmd:)} saves the underlying plot data.


{marker examples}{...}
{title:Examples}

{pstd}Example data ({stata "nopo_ex 1":{it:click to run})}{p_end}
{phang}. {stata webuse cattaneo2, clear}{p_end}
{phang}. {stata recode mage (min/18 = 1 "-18") (19/28 = 2 "19-28") (29/38 = 3 "29-38") (39/max = 4 "39-"), gen(mage_c)}{p_end}
{phang}. {stata lab var mage_c "Mother's age"}{p_end}
{phang}. {stata recode fage (min/18 = 1 "-18") (19/28 = 2 "19-28") (29/38 = 3 "29-38") (39/max = 4 "39-"), gen(fage_c)}{p_end}
{phang}. {stata lab var fage_c "Father's age"}{p_end}

{pstd}Decomposition - Standalone{p_end}
{phang}. {stata nopo decomp bweight mage_c fage_c prenatal1 mmarried fbaby foreign alcohol deadkids, by(mbsmoke)}{p_end}

{pstd}Decomposition - After {help kmatch:{it:kmatch}}{p_end}
{phang}. {stata kmatch em mbsmoke mage_c fage_c prenatal1 mmarried fbaby foreign alcohol deadkids (bweight), att}{p_end}
{phang}. {stata nopo decomp}{p_end}

{pstd}Postestimation{p_end}
{phang}. {stata nopo summarize, label}{p_end}
{phang}. {stata nopo summarize mrace, label}{p_end}
{phang}. {stata nopo gapoverdist}{p_end}
{phang}. {stata nopo dadb mage_c}{p_end}


{title:References}

{phang}
Jann, B. (2017). kmatch: Stata module for multivariate-distance and propensity-score matching,
including entropy balancing, inverse probability weighting, (coarsened) exact matching, and
regression adjustment. Available from {browse https://ideas.repec.org/c/boc/bocode/s458346.html}.

{phang}
{c N~}opo, H. (2008). Matching as a Tool to Decompose Wage Gaps. The Review of Economics and Statistics, 90(2), 290–299. {browse "https://doi.org/10/b6tqwq"}

{title:Acknowledgements}

{pstd}
Special thanks to Carla Rowold for stress testing and many helpful comments.


{title:Authors}

{pstd}
Maximilian Sprengholz,  {browse "mailto:maximilian.sprengholz@hu-berlin.de":maximilian.sprengholz@hu-berlin.de}{p_end}
{pstd}
Maik Hamjediers, {browse "mailto:maik.hamjediers@hu-berlin.de":maik.hamjediers@hu-berlin.de}{p_end}
{pstd}
Department of Social Sciences,{p_end}
{pstd}
Humboldt-Universität zu Berlin{p_end}

