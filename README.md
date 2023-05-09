Introduction to rtreefit
================
09/05/2023

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

<!-- badges: start -->
<!-- badges: end -->

You can install rtreefit like so:

``` r
devtools::install_github("nangalialab/rtreefit",build_vignettes=TRUE)
```

## Introduction

This package estimates time-based (“ultrametric”) trees, wherein the
y-axis of phylogenetic somatic mutation trees is converted from
mutations to time. The method jointly fits wild type rates, mutant rates
and absolute time branch lengths using a Bayesian per individual
tree-based model under the assumption that the observed branch lengths
are Poisson or Negative Binomial distributed with

Mean = Duration × Sensitivity × Mutation Rate

The method works with at most one change point per branch and supports
heterochronous sampling. See the rtreefit vignette for slightly fuller
mathematical details:

``` r
browseVignettes("rtreefit") 
```

## Branch Timings and Per Driver Clade Mutation Rate

We consider a rooted tree where each edge
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
consists of an observed mutation count
![m_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m_i "m_i")
and a true duration
![t_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_i "t_i").
We refer to a given edge and its child node interchangeably by the same
label. Now let
![D(i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;D%28i%29 "D(i)")
be the set of terminal nodes (tips) that descend from node
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
and let
![A(i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%28i%29 "A(i)")
be its corresponding set of ancestral nodes excluding the root. We
assume that each tip of the tree
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
has a known corresponding time
![T_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_k "T_k")
(e.g. the post conception age in years of the patient at sampling of the
cell) and so we therefore have the following constraint:

![T_k=\sum\_{i \in A(k)}t_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_k%3D%5Csum_%7Bi%20%5Cin%20A%28k%29%7Dt_i "T_k=\sum_{i \in A(k)}t_i")

and

![T_k\> t_i \> 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_k%3E%20t_i%20%3E%200 "T_k> t_i > 0")

We incorporate this constraint by performing the optimisation over the
interior branches of the tree with reparameterised branch durations
![x_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_i "x_i")
transformed to be in the range
![0\< x_i \<1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0%3C%20x_i%20%3C1 "0< x_i <1").
If
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j")
is an edge whose parent node is the root then:

![t_j=x_j  \text{min}({T_k:k \in D(j)})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_j%3Dx_j%20%20%5Ctext%7Bmin%7D%28%7BT_k%3Ak%20%5Cin%20D%28j%29%7D%29 "t_j=x_j  \text{min}({T_k:k \in D(j)})")

For other interior edges,
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i"),
we have

![t_i=\left(\text{min}\left\\{{T_k:k \in D(i)}\right\\}-\sum\_{j\in A(i)} t_j\right)x_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_i%3D%5Cleft%28%5Ctext%7Bmin%7D%5Cleft%5C%7B%7BT_k%3Ak%20%5Cin%20D%28i%29%7D%5Cright%5C%7D-%5Csum_%7Bj%5Cin%20A%28i%29%7D%20t_j%5Cright%29x_i "t_i=\left(\text{min}\left\{{T_k:k \in D(i)}\right\}-\sum_{j\in A(i)} t_j\right)x_i")

The duration of the terminal edges is fixed by the values of
![t_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_i "t_i")
on the interior edges and the overall duration constraint:

![t_i=\text{min}\left\\{{T_k:k \in D(i)}\right\\}-\sum\_{j\in A(i)} t_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_i%3D%5Ctext%7Bmin%7D%5Cleft%5C%7B%7BT_k%3Ak%20%5Cin%20D%28i%29%7D%5Cright%5C%7D-%5Csum_%7Bj%5Cin%20A%28i%29%7D%20t_j "t_i=\text{min}\left\{{T_k:k \in D(i)}\right\}-\sum_{j\in A(i)} t_j")

We assume that there are
![p-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p-1 "p-1")
change points in the tree corresponding to the acquisition of driver
mutations. This results in
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
mutation rates
![\lambda_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_j "\lambda_j")
applying throughout the tree where we allow at most one change point per
branch and the initial ancestral (or wild type) rate is
![\lambda_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_0 "\lambda_0")
and additional rate change points occur a fraction
![\alpha_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_j "\alpha_j")
along branch
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j")
and descendent branches have the rate
![\lambda_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_j "\lambda_j")
unless there are additional change points in descendant branches. The
effective rate on branches with a change point going from
![\lambda_l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_l "\lambda_l")
to
![\lambda_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_j "\lambda_j")
is just the weighted average
![\alpha_j \lambda_l+(1-\alpha_j)\lambda_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_j%20%5Clambda_l%2B%281-%5Calpha_j%29%5Clambda_j "\alpha_j \lambda_l+(1-\alpha_j)\lambda_j")
where we use a uniform unit interval prior for the
![\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha "\alpha")’s.

### Negative Binomial Model

We assume the underlying mutation process follows a Negative Binomial
Distribution with the above piecewise constant driver specific mutation
rates, the number of mutations accrued on branch
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
in time
![t_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_i "t_i")
measured in years:

![M_i \sim \text{NB}\left(\lambda\times t_i,\lambda\times t_i\times \phi\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M_i%20%5Csim%20%5Ctext%7BNB%7D%5Cleft%28%5Clambda%5Ctimes%20t_i%2C%5Clambda%5Ctimes%20t_i%5Ctimes%20%5Cphi%5Cright%29 "M_i \sim \text{NB}\left(\lambda\times t_i,\lambda\times t_i\times \phi\right)")

The number of observed mutations is:

![m_i \sim \text{Binomial}(M_i,s_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m_i%20%5Csim%20%5Ctext%7BBinomial%7D%28M_i%2Cs_i%29 "m_i \sim \text{Binomial}(M_i,s_i)")

Where we are using a per-branch estimated sensitivity
![s_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s_i "s_i")
that indirectly depends on the depth of sample and the number of samples
sharing a branch (see ?). This is equivalent too:

![m_i \sim \text{NB}\left(\lambda\times t_i \times s_i,\lambda\times t_i\times \phi\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m_i%20%5Csim%20%5Ctext%7BNB%7D%5Cleft%28%5Clambda%5Ctimes%20t_i%20%5Ctimes%20s_i%2C%5Clambda%5Ctimes%20t_i%5Ctimes%20%5Cphi%5Cright%29 "m_i \sim \text{NB}\left(\lambda\times t_i \times s_i,\lambda\times t_i\times \phi\right)")

with priors
![1/\phi \sim \text{HalfNormal}(0,10)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F%5Cphi%20%5Csim%20%5Ctext%7BHalfNormal%7D%280%2C10%29 "1/\phi \sim \text{HalfNormal}(0,10)"),
![\lambda \sim \mathcal{N}(\Lambda,0.25 \Lambda)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda%20%5Csim%20%5Cmathcal%7BN%7D%28%5CLambda%2C0.25%20%5CLambda%29 "\lambda \sim \mathcal{N}(\Lambda,0.25 \Lambda)")
where
![\Lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CLambda "\Lambda")
is the naive estimation of a single rate
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
as the per patient median of the ratio of the root to tip mutation count
and the tip sampling age, and finally we use the weakly informative
prior for the stick breaking fractions:

![x_i \sim \text{Beta}(\alpha=p_i/(1-\sum\_{j\in A(i)}p_j),\beta=1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_i%20%5Csim%20%5Ctext%7BBeta%7D%28%5Calpha%3Dp_i%2F%281-%5Csum_%7Bj%5Cin%20A%28i%29%7Dp_j%29%2C%5Cbeta%3D1%29 "x_i \sim \text{Beta}(\alpha=p_i/(1-\sum_{j\in A(i)}p_j),\beta=1)")

where the
![p_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_i "p_i")
is an initial approximation of the duration of the branch length
expressed as a fraction of the sampling time:

![p_i=\text{min}\_{j\in D(i)}\left\\{(m_j+1)/(\sum\_{k\in A(j)}\left(m_k+1\right))\right\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_i%3D%5Ctext%7Bmin%7D_%7Bj%5Cin%20D%28i%29%7D%5Cleft%5C%7B%28m_j%2B1%29%2F%28%5Csum_%7Bk%5Cin%20A%28j%29%7D%5Cleft%28m_k%2B1%5Cright%29%29%5Cright%5C%7D "p_i=\text{min}_{j\in D(i)}\left\{(m_j+1)/(\sum_{k\in A(j)}\left(m_k+1\right))\right\}")

Note that the overdispersion parameter is rescaled so that it is
comparable across branches with different mutation burden.

### Poisson Model

Here we assume the underlying mutation process follows a Poisson
Distribution again with the above piecewise constant driver specific
mutation rates, the number of observed mutations accrued on branch
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
in time
![t_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_i "t_i")
measured in years:

![m_i \sim \text{Poisson}(\lambda\times t_i\times S_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m_i%20%5Csim%20%5Ctext%7BPoisson%7D%28%5Clambda%5Ctimes%20t_i%5Ctimes%20S_i%29 "m_i \sim \text{Poisson}(\lambda\times t_i\times S_i)")

where

![S_i \sim \text{Beta}\left(\alpha=c,\beta=c(1-s_i)/s_i\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_i%20%5Csim%20%5Ctext%7BBeta%7D%5Cleft%28%5Calpha%3Dc%2C%5Cbeta%3Dc%281-s_i%29%2Fs_i%5Cright%29 "S_i \sim \text{Beta}\left(\alpha=c,\beta=c(1-s_i)/s_i\right)")

Where we have chosen the concentration parameter
![c=100](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%3D100 "c=100").
This reflects only modest uncertainty in our estimates in sensitivity
and also allows the model to mitigate larger than expected variability
in the branch lengths. In other respects the priors are the same as for
the Negative Binomial Model.

### Examples

## Neutral Case. One Rate

First lets simulate a neutral tree using rsimpop and fit the tree..

``` r
library("rtreefit")## Loads rsimpop as well
NYEARS=25
RATE=18
get_agedf_from_sim=function(simtree){
  st=get_elapsed_time_tree(simtree)## Gets "Real Time" ultrametric tree
  nh=nodeHeights(st)
  out=data.frame(tip.label=st$tip.label,age=nh[match(1:length(st$tip.label),st$edge[,2]),2]/365)
  out$age=ifelse(out$age<1e-6,1e-6,out$age)
  out
}
testing=run_neutral_sim(0.1,1/365,nyears=NYEARS)
#> n_sim_days: 9125
#> b_stop_if_empty: 0
#> b_stop_at_pop_size: 1
#> maxt: 0
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 18250 
#> MAX_SIZE= 300003 
#> n_sim_days: 9125
#> b_stop_if_empty: 0
#> b_stop_at_pop_size: 0
#> maxt: 152.215721405264
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 18250 
#> MAX_SIZE= 300003
st=get_subsampled_tree(testing,30)
#> Starting checking the validity of tmp...
#> Found number of tips: n = 31 
#> Found number of nodes: m = 30 
#> Done.
st=get_elapsed_time_tree(st,mutrateperdivision = 0,backgroundrate = RATE/365,odf=1)

plot_tree(st)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

    #> 
    #> Phylogenetic tree with 31 tips and 30 internal nodes.
    #> 
    #> Tip labels:
    #>   s1, s2, s3, s4, s5, s6, ...
    #> 
    #> Rooted; includes branch lengths.
    st$agedf=get_agedf_from_sim(st)
    res=fit_tree(tree=st,switch_nodes = c(),xcross = c(),niter = 10000,model = "poisson_tree",early_growth_model_on = 0.0)
    #> Warning in fit_tree(tree = st, switch_nodes = c(), xcross = c(), niter =
    #> 10000, : No sensitivity supplied: assuming 99%
    #> Median lambda estimate=18.34
    print(res$lambda)
    #> $mean
    #> [1] 18.37413
    #> 
    #> $sd
    #> [1] 0.1804631
    #> 
    #> $lb
    #> [1] 18.02837
    #> 
    #> $ub
    #> [1] 18.73339
    #> 
    #> $median
    #> [1] 18.37147
    par(mfcol=c(1,2))
    ut=get_elapsed_time_tree(st)
    ut$edge.length=ut$edge.length/365
    plot_tree(ut,cex.label = 0);title("True Ultrametric Tree")
    #> 
    #> Phylogenetic tree with 31 tips and 30 internal nodes.
    #> 
    #> Tip labels:
    #>   s1, s2, s3, s4, s5, s6, ...
    #> 
    #> Rooted; includes branch lengths.
    plot_tree(res$ultratree,cex.label = 0);title("Inferred Tree")
    #> 
    #> Phylogenetic tree with 31 tips and 30 internal nodes.
    #> 
    #> Tip labels:
    #>   s1, s2, s3, s4, s5, s6, ...
    #> 
    #> Rooted; includes branch lengths.

<img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" />

## Selection Case. Two fitted rates - both the same..

``` r
NYEARS=40
RATE=15
selsim=run_selection_sim(0.1,1/365,target_pop_size = 1e5,nyears_driver_acquisition = 5,nyears = NYEARS,fitness=0.3,minprop = 0.05)
#> n_sim_days: 1825
#> b_stop_if_empty: 0
#> b_stop_at_pop_size: 1
#> maxt: 0
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 3650 
#> MAX_SIZE= 300003 
#> n_sim_days: 1825
#> b_stop_if_empty: 0
#> b_stop_at_pop_size: 0
#> maxt: 128.203400258896
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 3650 
#> MAX_SIZE= 300003 
#> No driver found: tries= 0 
#>    population val fitness id driver1
#> 1           1   0     0.0  0       0
#> 2       99950   1     0.0  0       0
#> 21          1   1     0.3  1       1
#> n_sim_days: 14600
#> b_stop_if_empty: 1
#> b_stop_at_pop_size: 0
#> maxt: 1825.00043855212
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 29200 
#> MAX_SIZE= 300003 
#> No driver found: tries= 1 
#>    population val fitness id driver1
#> 1           1   0     0.0  0       0
#> 2       99950   1     0.0  0       0
#> 21          1   1     0.3  1       1
#> n_sim_days: 14600
#> b_stop_if_empty: 1
#> b_stop_at_pop_size: 0
#> maxt: 1825.00043855212
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 29200 
#> MAX_SIZE= 300003 
#> No driver found: tries= 2 
#>    population val fitness id driver1
#> 1           1   0     0.0  0       0
#> 2       99950   1     0.0  0       0
#> 21          1   1     0.3  1       1
#> n_sim_days: 14600
#> b_stop_if_empty: 1
#> b_stop_at_pop_size: 0
#> maxt: 1825.00043855212
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 29200 
#> MAX_SIZE= 300003 
#> No driver found: tries= 3 
#>    population val fitness id driver1
#> 1           1   0     0.0  0       0
#> 2       99950   1     0.0  0       0
#> 21          1   1     0.3  1       1
#> n_sim_days: 14600
#> b_stop_if_empty: 1
#> b_stop_at_pop_size: 0
#> maxt: 1825.00043855212
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 29200 
#> MAX_SIZE= 300003 
#> No driver found: tries= 4 
#>    population val fitness id driver1
#> 1           1   0     0.0  0       0
#> 2       99950   1     0.0  0       0
#> 21          1   1     0.3  1       1
#> n_sim_days: 14600
#> b_stop_if_empty: 1
#> b_stop_at_pop_size: 0
#> maxt: 1825.00043855212
#> driver_rate_per_cell_per_day: 0
#> max_driver_count: -1
#> MAX_EVENTS= 29200 
#> MAX_SIZE= 300003
st=get_subsampled_tree(selsim,30)
#> Starting checking the validity of tmp...
#> Found number of tips: n = 31 
#> Found number of nodes: m = 30 
#> Done.
st=get_elapsed_time_tree(st,mutrateperdivision = 0,backgroundrate = RATE/365,odf=1)
st$agedf=get_agedf_from_sim(st)
node=st$events$node[which(st$events$driverid==1)]
res=fit_tree(tree=st,switch_nodes = node,xcross = c(),niter = 10000,model = "poisson_tree",early_growth_model_on = 0.0)
#> Warning in fit_tree(tree = st, switch_nodes = node, xcross = c(), niter =
#> 10000, : No sensitivity supplied: assuming 99%
#> Median lambda estimate=15.33
print(res$lambda)
#> $mean
#> lambda[1] lambda[2] 
#>  15.26276  15.72464 
#> 
#> $sd
#> lambda[1] lambda[2] 
#> 0.1466417 0.5097143 
#> 
#> $lb
#> lambda[1] lambda[2] 
#>  14.98038  14.78105 
#> 
#> $ub
#> lambda[1] lambda[2] 
#>  15.55108  16.77022 
#> 
#> $median
#> lambda[1] lambda[2] 
#>  15.26248  15.70489
ut=get_elapsed_time_tree(st)
ut$edge.length=ut$edge.length/365
par(mfcol=c(1,2))
plot_tree(ut,cex.label = 0);title("True Ultrametric Tree")
#> 
#> Phylogenetic tree with 31 tips and 30 internal nodes.
#> 
#> Tip labels:
#>   s1, s2, s3, s4, s5, s6, ...
#> 
#> Rooted; includes branch lengths.
plot_tree(res$ultratree,cex.label = 0);title("Inferred Tree")
#> 
#> Phylogenetic tree with 31 tips and 30 internal nodes.
#> 
#> Tip labels:
#>   s1, s2, s3, s4, s5, s6, ...
#> 
#> Rooted; includes branch lengths.
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
