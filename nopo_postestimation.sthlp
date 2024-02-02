{smcl}
{* *! version 1.0.0   02feb2024  Maximilian Sprengholz & Maik Hamjediers}{...}
{vieweralsosee "kmatch" "kmatch"}{...}
{vieweralsosee "nopomatch" "nopomatch"}{...}
{viewerjumpto "Syntax" "nopo##syntax"}{...}
{viewerjumpto "Postestimation" "nopo##postest"}{...}
{viewerjumpto "Description" "nopo##desc"}{...}
{viewerjumpto "Options" "nopo##opts"}{...}
{viewerjumpto "Examples" "nopo##ex"}{...}
{title:Title}

{hi:nopo postestimation} {hline 2} Postestimation commands for {cmd:nopo decomp}

{marker description}{...}
{title:Postestimation commands}

{pstd}
The following postestimation commands are of special interest after {cmd:nopo decomp}: 

{synoptset 17}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
{synopt :{helpb nopo_postestimation##summarize:nopo summarize}}descriptive table by group and matching/weighting status{p_end}
{synopt :{helpb nopo_postestimation##gapoverdist:nopo gapoverdist}}plotting decomposition-components over the distribution of {depvar}{p_end}
{synopt :{helpb nopo_postestimation##dadb:nopo dadb}}plot that displays contribute to components {it:DA} and {it:DB}{p_end}
{synoptline}
{p2colreset}{...}


{marker summarize}{...}
{title:Description for nopo summarize}

{p 2 2 2}
{cmd:nopo summarize} [{varlist}] [{cmd:,} {cmdab:stat:istics(}{help tabstat##statname:{it:statnames}}{cmd:)} {it: displayoptions}]

{pstd}
Returns a descriptive table with selected statistics by group and matching/weighting status:
{it:A unmatched}, {it:A matched}, {it:A^B} or {it:B^A} (matched and weighted; depends on the matching direction), 
{it:B matched}, and {it:B unmatched}. If no variables are specified, {depvar} and 
{varlist} from {cmd: nopo decomp} are summarized. By default, means and standard deviations are
reported. For factor variables, the shares of factor levels are reported (indicate factors
by  using factor notation in {varlist} of either {cmd: nopo summarize} or the prior {cmd:nopo decomp}).

{dlgtab:Options of nopo summarize}

{phang}
{cmdab:stat:istics(}{help tabstat##statname:{it:statnames}}{cmd:)} are {it:mean sd} by default, but 
you can choose any other available {help tabstat##statname:{it:statnames}}.

{phang}
{it: displayoptions} comprise:

{p 6 8 2}
{cmd:label} displays labels instead of names for variables and value labels instead of values for 
factor variables.

{p 6 8 2}
{cmd:fvpercent} displays factor variable shares as percent.

{p 6 8 2}
{cmd:fvdummies} displays both dummy categories.

{p 6 8 2}
{cmd:keepempty} keeps empty columns for the unmatched. This is useful to keep the number of table
columns fixed irrespective of common support, which allows appending {it:r(table)} after matching 
on various specifications.

{dlgtab:Stored results of nopo summarize}

{pstd}
{cmdab:nopo summarize} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}matrix of the descriptive statistics by group and matching/weighting status{p_end}
{synopt:{cmd:r({it:statnames})}}matrix of each statistic of {cmdab:stat:istics()} by group and matching/weighting status{p_end}


{marker gapoverdist}{...}
{title:Description of nopo gapoverdist}

{p 2 2 2}
{cmd:nopo gapoverdist} [d d0 dx da db] [{cmd:,} {it:options}]

{pstd}
This command plots decomposition components based on the distribution of {depvar} by group and matching/weighting status 
(thus, {it:A unmatched}, {it:A matched}, {it:A^B} or {it:B^A} , {it:B matched}, and {it:B unmatched}). 
The same calculation is performed as in the overall decomposition,
but instead of the overall {depvar} means of each group, we plug in the component-and-group-specific {depvar} 
mean for each requested quantile. Consider {it:D0} and quartiles as an example. Here, we would take the difference 
between (1) the mean wage in Q1 of the wage distribution of the matched group {it:A} and the (2) mean wage in Q1 of
the wage distribution of the matched and reweighted group {it:B^A}; repeated for each quartile. The mean of the 
quantile-spcific values per component correspond to the estimates returned by {cmd:nopo decomp}.

{pstd}
Please note that the decomposition components might not add up to {it:D} within each quantile. This is usually the 
case without common support, simply because the {depvar} distribution might vary greatly between the matched and 
unmatched in each group. For example, when all units in {depvar} Q1 of groups {it:A} and {it:B} are matched, 
{it:D} for Q1 can only be based on the {depvar} values of the matched; and these values would always be exceeded
by the {depvar} values among the unmatched (assuming there are some).

{pstd}
You can specify which decomposition components to plot (default are all) by explicitly calling, for example,
{cmd:nopo gapoverdist d d0 dx}.

{dlgtab:Options of nopo gapoverdist}

{phang}
{cmdab:nq:uantiles(}#{cmd:)} specifies the number of quantiles to be plotted. Default is 100.

{phang}
{cmd:recast(}{it:string}{cmd:)} sets the {help twoway:{it:twoway plottype}}. Default is {it: line}.

{phang}
{cmd:twoptsd(}{it:string}{cmd:)}, {cmd:twoptsd0(}{it:string}{cmd:)},
{cmd:twoptsdx(}{it:string}{cmd:)}, {cmd:twoptsda(}{it:string}{cmd:)}, and {cmd:twoptsdb(}{it:string}{cmd:)}
set subplot-specific {help twoway_options:{it:twoway_options}} (should be according to specified {cmd:recast()-option}). 

{phang}
{cmd:twopts(}{it:string}{cmd:)} sets the general {help twoway_options:{it:twoway_options}}. 
Check the stored results of {cmd:nopo gapoverdist} to see defaults (e.g. {cmd:r(twopts)}). 
You can also use the {cmd:save()} option to build the plot from the underlying data yourself.

{phang}
{cmd:xsize(}{it:real}{cmd:)} and {cmd:ysize(}{it:real}{cmd:)} set the x and y dimensions of the plot.

{phang}
{cmd:nodraw} specifies that the plot is not produced at all. Useful in combination with {cmd:save()}
to build your own plot from the data.

{phang}
{cmd:save(}{help filename:{it:filename}}{cmd:)} saves the underlying plot data (which will be replaced if it already exists).

{dlgtab:Stored results of nopo gapoverdist}

{pstd}
{cmdab:nopo gapoverdist} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(ysize)}}value of option {cmd:ysize()} or system default (if retrievable) {p_end}
{synopt:{cmd:r(xsize)}}value of option {cmd:xsize()} or system default (if retrievable) {p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(twoptsdb)}}invoked options for the plot denoting component {it:DB}{p_end}
{synopt:{cmd:r(twoptsda)}}invoked options for the plot denoting component {it:DA}{p_end}
{synopt:{cmd:r(twoptsdx)}}invoked options for the plot denoting component {it:DX}{p_end}
{synopt:{cmd:r(twoptsd0)}}invoked options for the plot denoting component {it:D0}{p_end}
{synopt:{cmd:r(twoptsd)}}invoked options for the plot denoting overall gap {it:D}{p_end}
{synopt:{cmd:r(recast)}}invoked option for {cmd:recast()}{p_end}
{synopt:{cmd:r(twopts)}}invoked options for {cmd:twopts()}{p_end}


{marker dadb}{...}
{title:Description for nopo dadb}

{p 2 2 2}
{cmd:nopo dadb} {varname} [{cmd:,} {it:options}]

{pstd}
Creates a plot showing how the different levels of {varname} contribute to {it:DA} and {it:DB}
because these levels are associated with either many unmatched units and/or large differences in 
{depvar} by matching status within groups {it:A} and {it:B}
(see {browse "https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md":online documentation}
for details).

{pstd}
The bar-segments of the plot denote the difference between category-specific mean of unmatched units and 
the overall mean of matched units. The scatter-markers denote the contribution of this difference to {it:D}, 
which depends on the size of the difference and the number of unmatched units
in the respective category. Additionally, the number of unmatched units for each category are displayed.

{pstd}
Note that the output {it:is not} the same as a detailed decomposition in regression-based approaches 
(which is generally not possible with matching). The "contribution" to {it:D} refers only to the 
comparison between matched and {varname}-level-specific unmatched units among group {it:A} and 
{it:B}, and is interdependent with all other characteristics of the matching set besides {varname}.

{dlgtab:Options of nopo dadb}

{phang}
{cmd:nosort} specifies that the plot is not to be sorted by the mean values of {depvar} per level of
{varname}. Default is to sort in ascending order.

{phang}
{cmdab:desc:ending} specifies that the plot is to be sorted by the mean values of {depvar} per level of
{varname} in descending order. Default is to sort in ascending order.

{phang}
{cmd:force} skips the check for the number of levels of {varname}. By default, the program exits if
{varname} has more than 30 levels.

{phang}
{cmd:nmin(}{it:#}{cmd:)} specifies the minimum number of unmatched, weighted observations per level
of {varname} to be plotted. Default is 1. This option is useful if, for example, a small number of 
unmatched observations strongly deviate from the average value of {depvar} in that group, producing
large bars.

{phang}
{cmd:twopts(}{it:string}{cmd:)}, {cmd:twoptsbar(}{it:string}{cmd:)}, 
{cmd:twoptsscatter(}{it:string}{cmd:)}, {cmd:twoptsn(}{it:string}{cmd:)}, and
{cmd:twoptsby(}{it:string}{cmd:)} set the general, element-specific, and by-specific 
{help twoway_options:{it:twoway_options}}. Check the stored results of {cmd:nopo dadb} to see defaults 
(e.g. {cmd:r(twopts)}). You can also use the {cmd:save()} option to build the plot from the 
underlying data yourself.

{phang}
{cmd:xsize(}{it:real}{cmd:)} and {cmd:ysize(}{it:real}{cmd:)} set the x and y dimensions of the plot.

{phang}
{cmd:xlabfmt(}{help format:{it:%fmt}}{cmd:)} sets the format of the x-axis labels and range. Try different
formats when the gridlines of the top and bottom x-axes do not align (or use {cmd:save()} to build your own graph).

{phang}
{cmd:nodraw} specifies that the plot is not produced at all. Useful in combination with {cmd:save()}
to build your own plot from the data.

{phang}
{cmd:save(}{help filename:{it:filename}}{cmd:)} saves the underlying plot data (which will be replaced if it already exists).

{dlgtab:Stored results of nopo dadb}

{pstd}
{cmdab:nopo dadb} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(ysize)}}value of option {cmd:ysize()} or computed default {p_end}
{synopt:{cmd:r(xsize)}}value of option {cmd:xsize()} or computed default {p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(twoptsby)}}invoked options for {cmd:twoptsby()}{p_end}
{synopt:{cmd:r(twoptsn)}}invoked options for {cmd:twoptsn()}{p_end}
{synopt:{cmd:r(twoptsscatter)}}invoked options for {cmd:twoptsscatter()}{p_end}
{synopt:{cmd:r(twoptsbar)}}invoked options for {cmd:twoptsbar()}{p_end}
{synopt:{cmd:r(twopts)}}invoked options for {cmd:twopts()}{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Example data ({stata "nopo_ex 1":{it:click to run}}){p_end}
{phang}. {stata "use http://fmwww.bc.edu/RePEc/bocode/o/oaxaca.dta , clear"}{p_end}
{phang}. {stata "for any exper tenure: gen X_c = round(X,5)"}{p_end}
{phang}. {stata gen educ_c = round(educ,1)}{p_end}
{phang}. {stata lab var educ_c "years of educational attainment (rounded)"}{p_end}
{phang}. {stata lab var exper_c "years of work experience (5-year intervalls)"}{p_end}
{phang}. {stata lab var tenure_c "years of job tenure (5-year intervalls)"}{p_end}
{phang}. {stata lab def female 0 "Men" 1 "Women"}{p_end}
{phang}. {stata lab val female female}{p_end}

{pstd}Example decomposition{p_end}
{phang}. {stata nopo decomp lnwage educ_c exper_c tenure_c, by(female)}{p_end}

{pstd}Postestimation{p_end}
{phang}. {stata nopo summarize, label}{p_end}
{phang}. {stata nopo summarize i.educ_c, label}{p_end}
{phang}. {stata nopo gapoverdist}{p_end}
{phang}. {stata nopo gapoverdist d d0}{p_end}
{phang}. {stata nopo dadb tenure_c}{p_end}


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

