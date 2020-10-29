---
layout: post
title: A gentle introduction to OmicsPLS
subtitle: How to do cross-validation with OmicsPLS
gh-repo: selbouhaddani/OmicsPLS
gh-badge: [star, fork, follow]
tags: [OmicsPLS, cross-validation, Joint Principal Components, data integration, O2PLS]
comments: true
---

Implementing software for omics data integration can be challenging. The statistical theory is often complex, and the code needs to run fast. (To make things worse, there are many definitions and solutions for omics data integration!) Our data integration approach in OmicsPLS has statistical theory 'under the hood' (see [the article](https://doi.org/10.1186/s12859-015-0854-z)) and can be efficiently programmed. To actually use OmicsPLS in your analysis, there is some help available via the vignettes [with an actual data example](https://github.com/selbouhaddani/OmicsPLS/blob/master/OmicsPLS_vignette_realdata.pdf) and [without](https://github.com/selbouhaddani/OmicsPLS/blob/master/vignettes/OmicsPLS_vignette.pdf), and of course [the manual](https://rdrr.io/cran/OmicsPLS/man/). Still I believe that a basic introduction to the main functions in OmicsPLS will help you with getting everything out of this package. 

In this blog, I'll go through one of the functions that you will need at the start of the analysis: the cross-validation (CV) function `crossval_o2m_adjR2`. 

# The cross-validation function

From the [manual](https://rdrr.io/cran/OmicsPLS/man/crossval_o2m_adjR2.html), we can see the full CV function and its arguments:

    crossval_o2m_adjR2(X, Y, a, ax, ay, nr_folds, nr_cores = 1, 
    stripped = TRUE, p_thresh = 3000, q_thresh = p_thresh, 
    tol = 1e-10, max_iterations = 100)

First of all, make sure your omics datasets are stored in `X` and `Y`, with subjects across the rows (sometimes called subjects, or replicates) and variables across the columns (relations are sought between the variables, e.g. genes). 

Now the question is, how many components do you need? 

_Brief recap: What is a component?_ Each component represents a combination of variables, say genes. A joint component in `X` is a weighted combination of genes that has optimal covariance with another joint component in `Y`. Sometimes it helps to see a combination of genes as being a data-derived genetic pathway. The weights indicate how important each gene is in that component. For more information, [read the vignette](https://github.com/selbouhaddani/OmicsPLS/blob/master/OmicsPLS_vignette_realdata.pdf). 

With CV, you can try out several combinations of component numbers and see which one yields the best fit. Keeping this idea in mind, you specify some integers for `a`, `ax`, and `ay`. For example, take `a=1:3`, `ax=0:2` and `ay=0:2`. A good range depends of course on the data application. 

The next parameter is `nr_folds`, specifying in how many slices the rows will be cut. If you cut a pizza in 2 slices, you have 1 for fitting and 1 for independent tasting. Note that you need 2 rounds, first you fit the first slice and taste the second. Then you fit the second slice, and taste the first. Now replace the word "pizza" with "rows/subjects" and "tasting" with "testing" and you have a description of how CV works. Say we cut the pizza in 10 slices, then we can use 9 for fitting and 1 for tasting. And you can use each slice for tasting, making the total number of fitting and tasting to be 10. Of course, you need more rows to cut all slices, and choosing 1 slice leaves you without tasting slice. So for `nr_folds`, choose something between 2 and `nrow(X)`. Typically one chooses 5 or 10, but know that more slices takes more time.

The last argument of interest is `nr_cores`. If you have a computer with multiple cores (which you have, unless you've been hiding in a bunker the past decades), you can use let each core cut its own pizza and do the tasting. Be careful when `X` and `Y` are very large, and you are on a Windows machine! These matrices will be copied, quickly filling up your RAM. 

# Summary

The other arguments are not so exciting, and you usually don't need to adjust them. So to summarize:
1. Input your data `X` and `Y`
2. Choose a couple of integers that you want to try out
3. Choose the number of slices 
4. Choose the number of parallel workers (defaults to 1)

If you run the [example for the CV function](https://rdrr.io/cran/OmicsPLS/man/crossval_o2m_adjR2.html), you'll get something like this:

    minimum is at n = 1 
    Elapsed time: 0.51 sec
           MSE n nx ny
    1 1.969997 1  1  1
    2 2.087187 2  2  1
    3 2.093943 3  2  2
    4 2.057414 4  1  1

The advice based on CV is here: take 1 joint component, and also 1 specific component for `X` and `Y`. Now you can proceed to the `o2m` function, and use these numbers of components. 


****

Questions, comments, tips? Let me know!
