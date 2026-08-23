# Difference-in-Differences with Bad Controls: A Conceptual Guide

\\\newcommand{\E}{\mathbb{E}} \newcommand{\ATT}{\text{ATT}}
\newcommand{\independent}{\mathrel{\perp\\\\\\\perp}}\\

`badcontrols` implements difference-in-differences methods for settings
in which a time-varying covariate is affected by treatment (a “bad
control”). This vignette provides an intuitive overview of bad controls
based on our paper Caetano et al.
([2026](#ref-caetano-callaway-payne-santanna-2026)).

## The bad-control tension

Difference-in-differences is often more credible after conditioning on
covariates, but it is common in empirical applications to have variables
that we would like to (in some sense) control for that themselves may be
affected by the treatment. In the application in our paper, we study the
effects of job displacement on a person’s earnings. In this context, it
seems reasonable to want to compare people that work in similar
occupations, but job displacement can itself change a person’s
occupation. In this case, occupation is a “bad control” because (1) it
is affected by the treatment but also (2) seems to be an important
covariate.

Since difference-in-differences identification strategies have
traditionally been implemented using a two-way fixed effects (TWFE)
regression like \\ Y\_{it} = \theta_t + \eta_i + \alpha D\_{it} + \beta
X\_{it} + e\_{it} \\ where \\\theta_t\\ and \\\eta_i\\ are time and unit
fixed effects, \\D\_{it}\\ is the treatment indicator, and \\X\_{it}\\
is a covariate, for a researcher that is worried about \\X\_{it}\\ being
a bad control, the natural question is whether to include \\X\_{it}\\ in
the regression or not (both of which have problems):

- **Option 1: Include the bad control:** Including the bad control
  results in conditioning on a post-treatment variable. Comparisons are
  then made among treated and untreated units with the same observed
  post-treatment covariate, which is not the right comparison when the
  covariate is affected by the treatment.

This option is not very popular in empirical work in economics—to our
knowledge, the only times empirical researchers ever purposefully
include a bad control as a regressor are in supplementary robustness
checks. It has been widely criticized, probably most notably in *Mostly
Harmless Econometrics* ([Angrist and Pischke
2008](#ref-angrist-pischke-2008)).

- **Option 2: Drop the bad control:** In other words, run the regression
  above dropping the \\\beta X\_{it}\\ term. This approach avoids
  conditioning on a post-treatment variable, but it drops the covariate
  entirely. In our application, this would mean that we would not even
  attempt to compare individuals in the same occupations at all.

In contrast to Option 1, Option 2 is very popular in empirical work. To
us, it seems like this has its own issues, that could be just as bad as
including the bad control. To drop the bad control entirely effectively
says that the bad control was never needed as a covariate at all. Or, in
words, the leading solution the bad control problem is to pretend that
the covariate was never needed in the first place.

## Alternative approaches

We provide alternative approaches in our paper. The key insight involves
(1) introducing a notion of treated and untreated potential versions of
the bad control and (2) trying to control for the untreated potential
version of the bad control.

To provide more detail, let us introduce some more notation. We will
consider the setting with two time periods: \\t^\*\\ and \\t^\*-1\\. In
the first period, no units are treated, and in the second period, some
units (the treated group) become treated while other units (the
untreated group) remain untreated. We use \\D_i\\ to be an indicator for
being in the treated group, \\Y\_{it}\\ to denote the outcome, and
\\X\_{it}\\ to denote the bad control. Finally, let \\Z_i\\ denote a
vector of other exogenous covariates, which can include time-invariant
and/or time-varying covariates. Next, let \\Y\_{it}(1)\\ and
\\Y\_{it}(0)\\ denote the treated and untreated potential outcomes,
respectively, and let \\X\_{it}(1)\\ and \\X\_{it}(0)\\ denote the
treated and untreated potential bad controls, respectively. Since no one
is treated in the first period, we have \\Y\_{it^\*-1} =
Y\_{it^\*-1}(0)\\ and \\X\_{it^\*-1} = X\_{it^\*-1}(0)\\. In the second
period, we observe \\Y\_{it^\*} = D_i Y\_{it^\*}(1) + (1-D_i)
Y\_{it^\*}(0)\\ and \\X\_{it^\*} = D_i X\_{it^\*}(1) + (1-D_i)
X\_{it^\*}(0)\\ (i.e., we observe treated potential outcomes and treated
potential bad controls for treated units, and untreated potential
outcomes and untreated potential bad controls for untreated units). This
notation explicitly allows for \\X\_{it}\\ to be a bad control as it can
be the case that \\X\_{it}(1) \neq X\_{it}(0)\\. Following the
difference-in-differences literature, we target the average treatment
effect on the treated (ATT):

\\ \ATT = \E\[Y\_{t^\*}(1) - Y\_{t^\*}(0) \mid D=1\]. \\

and we make the following parallel trends assumption:

\\ \E\[\Delta Y\_{t^\*}(0) \mid X\_{t^\*}(0), X\_{t^\*-1}, Z, D=1\] =
\E\[\Delta Y\_{t^\*}(0) \mid X\_{t^\*}(0), X\_{t^\*-1}, Z, D=0\], \\

This assumption is exactly the same as the conditional parallel trends
assumptions in the difference-in-differences literature ([Heckman et al.
1997](#ref-heckman-ichimura-todd-1997); [Callaway and Sant’Anna
2021](#ref-callaway-santanna-2021)). It says that, in the absence of the
treatment, the trend in untreated potential outcomes would have been the
same for treated and untreated units with the same covariates. Notice
that, unlike Options 1 and 2 above, this assumption uses \\X\_{it}\\ as
a genuine control while also allowing for it to be affected by the
treatment.

The challenge is that the untreated potential version of the bad
control, \\X\_{t^\*}(0)\\, is not observed for treated units. This makes
the explicit comparison between treated and untreated units with the
same \\X\_{t^\*}(0)\\ infeasible. We propose two approaches to dealing
with this challenge, though we note that any approach is going to
require some additional assumptions, and so a real possibility in a
given application is that the researcher is just not able to credibly
recover the \\\ATT\\ in the presence of a bad control.

## New Approach 1: Condition on the pre-treatment bad control

A very simple approach to dealing with a bad control is to condition
only on its value on the pre-treatment period, \\X\_{t^\*-1}\\. This is
not so straightforward when implementing difference-in-differences via a
TWFE regression, but it is straightforward with the Callaway and
Sant’Anna ([2021](#ref-callaway-santanna-2021)) estimator. In fact, in
their `R` implementation in the `did` package, the default approach to
dealing with a time-varying covariate is to include its level in the
pre-treatment period, which is exactly what we suggest here.

What additional assumptions rationalize this approach? In the paper, we
provide two cases where it will hold. One under the following
unconfoundedness assumption for the bad control (we call this assumption
“simple covariate unconfoundedness” in our paper):

\\ X\_{t^\*}(0) \independent D \mid X\_{t^\*-1}, Z. \\

This assumption implies that we can learn about how the untreated bad
control would have evolved in the absence of treatment, by comparing
untreated units with the same pre-treatment value of the bad control and
other covariates. The post-treatment bad control does not need to be
included as a covariate in the final outcome comparison as this
assumption implies that it is already balanced between the treated and
untreated group conditional on the pre-treatment bad control and other
covariates. The other case that rationalizes this approach is when the
\\X\_{t^\*}(0)\\ is redundant in the parallel trends assumption after
conditioning on \\X\_{t^\*-1}\\ and \\Z\\, i.e., it does not separately
affect the trend in untreated potential outcomes after conditioning on
the pre-treatment bad control and other covariates.

## New Approach 2: Covariate unconfoundedness for the bad control

The specific version of unconfoundedness that we discussed above leads
to a convenient estimator, but it’s not necessarily the most plausible
one. Notice that, the parallel trends assumption involved the *change*
in untreated potential outcomes, while the left hand side involved the
*level* of the bad control (rather than its change). One would expect
that the number of conditioning variables needed to make an
unconfoundedness assumption plausible is larger for a level than a
change, which suggests a modified unconfoundedness assumption (we call
this assumption “covariate unconfoundedness” in our paper):

\\ X\_{t^\*}(0) \independent D \mid X\_{t^\*-1}, W, Z. \\

What is \\W\\? The short version is that it is any additional covariates
that are needed to make unconfoundedness hold. But a leading example is
to take \\W\\ to be the pre-treatment outcome. In our application, this
would say that, in the absence of treatment, occupation would have
evolved in the same way for displaced and non-displaced workers with the
same pre-treatment occupation, pre-treatment earnings, and other
covariates. Under this assumption,

\\ \ATT = \E\[\Delta Y\_{t^\*} \mid D=1\] - \E\Big\[ \E\big\[ \E\[\Delta
Y\_{t^\*} \mid X\_{t^\*}, X\_{t^\*-1}, Z\] \bigm\| X\_{t^\*-1}, W, Z,
D=0 \big\] \Bigm\| D=1 \Big\] \\

A rough explanation for this equation (starting from the inside out) is
that we find the trend in outcomes over time as a function of
\\X\_{t^\*}, X\_{t^\*-1}, Z\\ among untreated units. Then, we account
for the bad control being affected by the treatment by using the
distribution of bad control from the untreated group given exogenous
covariates \\Z\\ and \\W\\ but using the distribution of these from the
treated group—you can see the paper for more details.

`badcontrols` implements imputation, doubly robust, and machine learning
estimators building on this result.

### Conceptual side-discussion:

The approaches that we emphasize in the paper both involve
unconfoundedness assumptions. But it’s helpful to think about dealing
with bad controls at a slightly higher level. In particular, our
approach fundamentally proposes an identification strategy for the
outcome that allows us to back out the distribution of untreated
potential bad controls for the treated group. Then, we use this
counterfactual distribution in the parallel trends assumption to recover
the ATT.

An implication of this discussion is that one could use alternative
identifying assumptions in place of unconfoundedness. A natural
alternative would be to assume parallel trends for the bad control too.
The challenge here is that parallel trends is usually useful for
recovering the mean of untreated potential outcomes rather than its
entire distribution (see, e.g., Athey and Imbens
([2006](#ref-athey-imbens-2006)); Callaway and Li
([2019](#ref-callaway-li-2019))), and we really need the entire
distribution. That said, we show in the paper that, under an additional
linearity assumption, one can use parallel trends for the bad control.
We implement this estimator in `badcontrols` as
`didbc(..., est_method = "imputation", bad_control_identification_strategy = "did")`.

## Conclusion

Our paper provides new approaches to dealing with bad controls in
difference-in-differences applications. Our approaches respect the
nature of a bad control, namely, that it is needed as a covariate for
identification but it is also affected by the treatment.

## References

Angrist, Joshua D, and Jorn-Steffen Pischke. 2008. *Mostly Harmless
Econometrics: An Empiricist’s Companion*. Princeton University Press.

Athey, Susan, and Guido Imbens. 2006. “Identification and Inference in
Nonlinear Difference-in-Differences Models.” *Econometrica* 74 (2):
431–97.

Caetano, Carolina, Brantly Callaway, Stroud Payne, and Hugo Sant’Anna.
2026. “Difference-in-Differences with Bad Controls.” *arXiv Preprint
arXiv:2608.03881*. <https://arxiv.org/abs/2608.03881>.

Callaway, Brantly, and Tong Li. 2019. “Quantile Treatment Effects in
Difference in Differences Models with Panel Data.” *Quantitative
Economics* 10 (4): 1579–618.

Callaway, Brantly, and Pedro HC Sant’Anna. 2021.
“Difference-in-Differences with Multiple Time Periods.” *Journal of
Econometrics* 225 (2): 200–230.

Heckman, James, Hidehiko Ichimura, and Petra Todd. 1997. “Matching as an
Econometric Evaluation Estimator: Evidence from Evaluating a Job
Training Programme.” *The Review of Economic Studies* 64 (4): 605–54.
