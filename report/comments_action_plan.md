# Action Plan for `report/comments.md`

This checklist turns the comments into work items we can verify as we revise.
The main rule is: fix factual/reporting mismatches first, then strengthen the
methodology text, and only then treat new design-aware validation results as the
final evidence layer.

## Source Files

- [ ] Main thesis text: `report/MSE_Thesis_Decomposition_of_CI_Trees_MM.qmd`
- [ ] Results generator: `report/Generate_DRC_Results_Objects.qmd`
- [ ] Current results object: `data/drc_report_results_objects.rda`
- [ ] Plot/table helpers, if labels need changing: `R/report_helpers.R`
- [ ] Reviewer notes: `report/comments.md`

## Phase 1: Lowest-Hanging Reporting Fixes

These should be corrected promptly because they are mismatches between reported
numbers, captions, and interpretation. Most should not require new methods.

- [ ] Rebuild or reload the current results object after the `shapr` and all-row
      SHAP change, then treat the regenerated tables as the numeric source of
      truth.
- [ ] Fix the predictive-forest SHAP text in the results section. The current
      text says residence and household wealth each contribute about 12%, but
      the table/comment says residence is closer to 21.6% and household wealth
      around 11.5%. Update the sentence to match the regenerated table exactly.
- [ ] Fix the same SHAP wording in the discussion/conclusion wherever the
      "about 12%" claim appears.
- [ ] Resolve the tree validation-gain mismatch. The results text reports about
      -75.5% for CI, while the appendix/table comment reports about -68.2%.
      Pull the selected values from the current tree-selection table and update
      all text to the same numbers.
- [ ] Clarify non-integer terminal-node counts. If values such as 6.8 or 7.3
      are cross-validation averages, rename the column/caption as mean terminal
      nodes across folds. If discussing a final fitted tree, report an integer
      terminal-node count from the final fitted object.
- [ ] Separate absolute impurity from signed inequality. The tree split
      criterion uses absolute CI/CIg/CIc/L values, but the results should also
      state the sign of inequality where direction matters.
- [ ] Update tree and surrogate-tree captions so readers know whether node
      labels show absolute impurity, signed concentration index, observed
      mortality rate, or average model score.
- [ ] Soften the surrogate forest interpretation. Replace high/low "mortality"
      wording with precise language: observed mortality rate when using the
      observed outcome, average forest score when using the surrogate outcome,
      and never calibrated risk unless calibration is actually performed.
- [ ] Update the abstract and summary language so the L result is described as
      a sensitivity/descriptive stability result, not the main defensible direct
      result.

## Phase 2: Methodological Explanation and Limitations

These changes strengthen the thesis without necessarily requiring new results.
They should make clear what was estimated, what is descriptive, and what remains
limited.

- [ ] Explain exactly how the `rineq` decomposition was done.
- [ ] Verify from code/documentation whether `rineq::contribution()` used the
      logistic-regression coefficients on the link scale, marginal effects, or
      another transformation.
- [ ] Add a short formula-level description of the regression decomposition
      estimand, including the role of survey weights, the binary outcome, and
      whether the result decomposes observed mortality, fitted mortality, or a
      model-scale quantity.
- [ ] If the current `rineq` decomposition uses logit-scale coefficients without
      marginal-effect correction, mark it as a limitation or redo it using a
      defensible probability-scale decomposition.
- [ ] Explain that DHS weights were used as case weights, while clustering and
      stratification were not fully propagated through tree validation, forest
      uncertainty, or SHAP uncertainty.
- [ ] Add a limitation paragraph separating weighted descriptive analysis from
      formal design-based inference.
- [ ] Reframe L clearly. CI, CIg, and CIc are rank-dependent and appropriate for
      ranking by DHS wealth index; L uses socioeconomic levels and is weaker
      with a centered PCA wealth score. Therefore L should be presented as a
      sensitivity analysis unless a ratio-scale socioeconomic variable is used.
- [ ] Explain surrogate forests carefully. A surrogate tree is an interpretable
      approximation to the forest's fitted score, not the forest itself and not
      a calibrated mortality model.
- [ ] Strengthen SHAP limitations. State that SHAP decomposes fitted predictions,
      depends on the model and SHAP convention, and should not be interpreted as
      a causal determinant decomposition when DHS covariates are correlated.
- [ ] Revisit the predictive forest after undersampling. If the class
      distribution was changed and predictions were not recalibrated, describe
      the SHAP decomposition as a decomposition of model scores rather than
      calibrated mortality probabilities.
- [ ] Add a model-performance caveat: modest ROC AUC reduces confidence in
      detailed SHAP feature rankings, even if the decomposition algebra is valid.
- [ ] Reframe the simulation study as a sanity check rather than full method
      validation unless repeated simulation results are added.

## Phase 3: Results That Need New Computation

These are heavier tasks. They should be done only after Phase 1 text mismatches
are fixed and Phase 2 limitations are clear.

- [ ] Regenerate the DRC results object after the `shapr` replacement and
      all-sample SHAP change.
- [ ] Confirm whether the DHS data contain PSU/cluster and strata variables
      needed for design-aware validation.
- [ ] Add design-aware validation for trees, ideally by resampling primary
      sampling units within strata rather than individual rows.
- [ ] Extend the same design-aware validation idea to forests, including the
      forest validation-gain summaries.
- [ ] Add uncertainty or stability summaries for tree structure and selected
      splits under cluster/strata resampling.
- [ ] Add design-aware or cluster-aware sensitivity checks for SHAP
      contributions if computationally feasible.
- [ ] Decide whether the final headline result is still supported after
      design-aware validation. If not, revise the thesis claim to "promising
      descriptive method with unstable out-of-sample rank-dependent validation."
- [ ] If ratio-scale socioeconomic data become available, rerun the L analysis
      with that variable. Otherwise keep L as sensitivity only.

## Phase 4: Final Consistency Pass

- [ ] Check every table reference against the rendered table number/caption.
- [ ] Check all percentages in the abstract, results, discussion, and conclusion
      against the regenerated tables.
- [ ] Check that "validation gain", "relative validation gain", and
      "validation gain/root" are not mixed.
- [ ] Check that "terminal nodes" always distinguishes final fitted-tree counts
      from cross-validation averages.
- [ ] Check that all uses of "risk", "mortality rate", "prediction", "score",
      and "calibrated probability" are precise.
- [ ] Check that CI signs are shown or discussed wherever interpretation of
      pro-poor versus pro-rich inequality matters.
- [ ] Render the thesis and confirm there are no broken cross-references.
- [ ] Compare this checklist against the finished thesis before submission.

## Definition of Done

- [ ] All numeric claims match the regenerated tables.
- [ ] L is no longer positioned as the main defensible result with DHS wealth;
      it is clearly a sensitivity/descriptive stability result.
- [ ] `rineq` decomposition is explained with its exact estimand and limitation.
- [ ] Surrogate forest language distinguishes observed mortality from forest
      scores.
- [ ] Survey-design limitations are explicit.
- [ ] Any final strong claim about stability is backed by design-aware
      validation results, especially for forests.
