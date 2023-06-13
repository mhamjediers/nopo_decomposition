{smcl}
{* *! version 0.2   11june2023  Maximilian Sprengholz & Maik Hamjediers}{...}

{vieweralsosee "" "--"}{...}
{viewerjumpto "Syntax" "npop##syntax"}{...}
{viewerjumpto "Postestimation" "nopo##postest"}{...}
{viewerjumpto "Description" "nopo##description"}{...}
{viewerjumpto "Options" "nopo##comopt"}{...}
{viewerjumpto "Examples" "npoo##examples"}{...}

help for {hi:nopo}
{hline}

{title:Title}

{pstd}{hi:nopo} {hline 2} Matching based decomposition analysis of differences in outcomes between two groups

{marker:syntax}{...}
{title:Syntax}

   {cmd:nopo decomp} {depvar} {varlist} {ifin} {weight} {cmd:, by(} varname{cmd:)} [{help nopo##comopt:{it:options}}]
   
 Decomposition as post-estimation to {help kmatch:{it:kmatch}}:
   
   {cmd:nopo decomp}

Allows for three different postestmation commands (nopo summarize, nopo dadb, nopo gapoverdist: described below)
 
 
{synoptset 37 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{cmdab:by(}varname{cmdab:)}} variable by which to estimate and decompose gaps; necessary {p_end}{...}
{synopt :{cmdab:swap}} 		{p_end}{...}
{synopt :{cmdab:bref(}string{cmdab:)}}		{p_end}{...}
{synopt :{cmdab:norm:alize}}   	{p_end}{...}

    Options for inbuilt matching via {help kmatch:{it:kmatch}} ({it: not if used as {help kmatch:{it:kmatch}}-post-estimation})
{synopt :{cmdab:km:atch(}em|md|ps{cmdab:)}} Choose matching approach; either exact (em, default), multivariate-distance (md), or propensity score (ps) matching {p_end}{...}
{synopt :{cmdab:kmo:opts(}string{cmdab:)}}	Matching-approach-specific options (either {help kmatch##emopts{it:emopts}}, {help kmatch##mdoptions{it:mdoptions}}, or {help kmatch##psoptions{it:psoptions}})  {p_end}{...}
{synopt :{cmdab:kmpass:thru(}string{cmdab:)}}   	{p_end}{...}
{synopt :{cmdab:kmkeep:gen}}    	{p_end}{...}
{synopt :{cmdab:kmnois:ily}}   	{p_end}{...}

{synoptline}
{p2colreset}

{marker description}{...}
{title:Description}

{pstd}
{hi:nopo decomp} does this and this.

{pstd}
It's important to not this and this.

{pstd}
Can also be used with bootstrapping.


{marker comopt}{...}
{title:Options}

{phang}
{cmdab:by(}varname{cmdab:)} specifies to plot the relative frequencies of the levels of the plotted variables. The default is to plot absolute frequencies.

{phang}
{cmdab:swap} specifies the number of points between bars through which the flows are drawn. The default are 10 points, with a larger number increasing the smoothness and computational time.

{phang}
{cmdab:bref(}string{cmdab:)} specifies to plot missing values in {it:varlist} as well. The default is to disregard them. Note that missing values can lead to an imbalance between start and end bars or they produce flows into empty bars. Instead of using {cmd:missing}, rather code missing values as distinct values in the underyling variables.

{phang}
{cmdab:norm:alize} overrules the default of not plotting variables with more than seven unique levels. Plotting variables with many levels likely impairs readability and might take relatively long. 

{dlgtab:Options for inbuilt matching via kmatch}

{phang}
{cmdab:km:atch(}em|md|ps{cmdab:)} specifies the width of the bar. The default is {cmd:barwidth(0.1)}; values greater than 0.5 are not allowed.

{phang}
{cmdab:kmo:opts(}string{cmdab:)}}  specifies the colors and opacity of the bars and flows. List the colors in the increasing order of the values of the plotted variables.

{phang}
{cmdab:kmpass:thru(}string{cmdab:)}  allow to colors and opacity of flows, independently of the option {cmd:colors()}. List the colors in the increasing order of the values of the plotted variables.

{phang}
{cmdab:kmkeep:gen} specifies the opacity of the colors of bars and flows. The default is {cmd:opacity(80)}.

{phang}
{cmdab:kmnois:ily}} allows most of the {help barlook_options} to change the look of the bars. Note that any specified option in {cmd:baroptions()} affect all bars; specifying for instance {cmd:baroptions(color(navy))} overrides any specification of the above option {cmd:colors()} and depicts all bars in {it:navy}.


{marker postest}{...}
{title:Postestimation-commands}
  
{dlgtab:nopo summarize}

   {cmd:nopo summarize} [{varlist}] [{cmd:,} label]
   
       Reports a descriptive table with means and standard devaitions by group-indicator and matching status.
   
       If no variables are specified, {depvar} and {varlist} from previous {cmd: nopo decomp} are summarized as the default.
	   
	   Option {it:label} display variable-labels instead of varible names.
      
{dlgtab:nopo gapoverdist}
   
   {cmd:nopo gapoverdist} [{cmd:,} nquantiles(#) {help twoway_options:{it:twoway_options}}]

       Plotting decomposition-components over the distribution of {depvar}.
	   
	   {cmdab:nq:uantiles(}#{cmd:)} specifies the number of quantiles to be plotted. 
	   
	   {it:twoway_options} allows for the inclusion of any additional {help twoway_options:graphing options} such as titles, axes, added lines, etc.

{dlgtab: nopo dadb}

   {cmd:nopo dadb} {varname}  [{cmd:,} {help twoway_options:{it:twoway_options}}]
   
       Plotting contribution of {varname} to components {it:D_A} and {it:D_B}.
	   
	   {it:twoway_options} allows for the inclusion of any additional {help twoway_options:graphing options} such as titles, axes, added lines, etc.
   
   
{title:Options for postestimation-commands:}



{marker examples}{...}
{title:Examples}

{pstd}Setup ({stata "sankeyplot_eg 1":{it:click to run}}) {p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 987}{p_end}
{phang2}{cmd:. gen edu_0 = _n <= 500}{p_end}
{phang2}{cmd:. replace edu_0 = 2 if _n > 900}{p_end}
{phang2}{cmd:. gen edu_1 = runiformint(0,2)}{p_end}
{phang2}{cmd:. lab def edu 0 "Low Edu." 1 "Medium Edu." 2 "High Edu."}{p_end}
{phang2}{cmd:. lab val edu_0 edu_1 edu}{p_end}

{pstd}Simple run {p_end}
{phang2}{cmd:. {stata "sankeyplot edu_0 edu_1"}}{p_end}




{title:Acknowledgements}

Thanks to Carla Rowold for useful comments.


{title:Authors}

{pstd}
Maximilian Sprengholz,  {browse "mailto:maximilian.sprengholz@hu-berlin.de"} 
Maik Hamjediers, Department of Social Sciences, Humboldt-Universität zu Berlin. {browse "mailto:maik.hamjediers@hu-berlin.de"}
{p_end}
