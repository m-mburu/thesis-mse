# Narrative Consistency Review After Analysis Changes

Source reviewed: `report/MSE_Thesis_Decomposition_of_CI_Trees_MM.tex`

Scope requested: Results first, then Discussion, then Abstract. I also note a few method/appendix items where they directly affect the results narrative.

Main framing from `report/comments.md`:

- The rank-dependent CI, CIg, and CIc trees can be described as plausible in-sample subgroup discovery, but not as a stable applied method unless the new PSU-fold and bootstrap results support that claim.
- The L result should be framed as a useful sensitivity analysis or descriptive partition because the DHS wealth index is a PCA-derived relative asset score, not a ratio-scale socioeconomic variable.
- The survey design should be separated clearly into descriptive weighted analysis, PSU-held-out validation, and bootstrap stability checks. These checks are not the same as full formal survey-design inference.
- SHAP should be described as an explanation of a fitted predictive model, not a causal or population-level decomposition of determinant effects.

## Results Section

### [ ] R1. Indirect SHAP narrative does not match Table 1 values

Location: Results, lines 1701-1715; Table `tbl-results-method-comparison`, lines 1721-1759.

Current issue:

The text says that in the predictive-forest SHAP decomposition, "Residence and household wealth each contributed about 12% under CI, CIg and CIc." The table does not show this. Under CI, CIg, and CIc, residence is about 19.4%, while household wealth is about 8.3%. For L, residence is about 24.0%, while household wealth is about 6.2%.

Suggested replacement in the current Results style:

> From the indirect decomposition, the largest CI regression contribution was mother's education, followed by father's occupation in agriculture, mother's occupation in agriculture, household wealth, residence, and skilled birth attendance (Table~\ref{tbl-results-method-comparison}). In the predictive-forest SHAP decomposition, the contributions were less concentrated in the education and occupation variables. Residence gave the largest SHAP contribution under CI, CIg, CIc, and L. Under CI, CIg, and CIc, residence contributed about 19%, followed by father's agricultural occupation, household wealth, mother's education, mother's agricultural occupation, and skilled birth attendance. Under L, residence was again the largest contributor, followed by father's agricultural occupation, mother's education, Maniema, mother's agricultural occupation, Kinshasa, and household wealth. The fitted regression model is reported in Table~\ref{tbl-appendix-regression-model-summary}. Details of the fitted predictive forest are provided in Table~\ref{tbl-appendix-ranger-sampling-summary} and Table~\ref{tbl-appendix-ranger-cross-validation}, and its SHAP contributions are shown in Figure~\ref{fig-appendix-shap-summary}.

Also check:

- The SHAP method now uses `fastshap`, not `shapr`. The package name is mainly in Methodology, but the Results section should avoid wording that implies a different SHAP engine.
- The ranger cross-validation table should be regenerated because the generator now uses PSU-held-out ranger folds.

### [ ] R2. Tree validation prose uses old numbers and conflicts with appendix tables

Location: Results, lines 1867-1887; Appendix tree tables, lines 3183-3239.

Current issue:

The Results text says:

- rank-dependent trees removed approximately 87 to 88% of root-node inequality in training;
- validation gains were -75.5% for CI and -73.1% for CIc and CIg;
- the L tree removed 70.5% of root-node inequality in held-out folds;
- the L tree had an average of 6.8 terminal nodes.

The appendix tables currently show different values:

- training gain/root is about 83.3% for CI and 81.6% for CIc and CIg;
- validation gain/root is shown as very negative for CI, CIc, and CIg;
- L validation gain/root is about 64.5%;
- mean terminal nodes are 6.9 for CI, 6.0 for CIc, 6.0 for CIg, and 7.5 for L.

Suggested replacement after rerendering the generator:

> The rank-dependent trees performed well on the training data but poorly on the PSU-held-out validation folds. The CI, CIc, and CIg trees retained large training gains, but their held-out validation gains were negative at the selected settings (Figure~\ref{fig-appendix-tree-validation-path}, Table~\ref{tbl-appendix-tree-selected-setting-summary}). In this criterion, a negative validation gain means that the fitted partition did not reduce inequality when applied to the held-out fold. Thus, although these trees identified substantial rank-based structure in the fitted sample, the same partitions did not transfer well to held-out PSUs. Unlike the rank-dependent trees, the level-dependent tree retained a positive held-out validation gain at the selected setting. The rural-urban, occupational, wealth, and regional partition identified by the L tree was therefore more stable under the validation criterion used in this thesis.

Reason for using this form:

The Results section already uses descriptive language, not heavy methodological caveats. This replacement keeps that style while avoiding stale exact values. Once the QMD is rerendered, exact numbers can be reinserted if needed.

### [ ] R3. Tree complexity prose has old terminal-node values

Location: Results, lines 1902-1909; Appendix tree selection table, lines 3227-3239.

Current issue:

The text reports average terminal nodes as 7.5 for CI, 7.3 for CIc and CIg, and 6.8 for L. The current appendix table shows 6.9 for CI, 6.0 for CIc, 6.0 for CIg, and 7.5 for L.

Suggested replacement:

> Tree complexity decreased as the minimum relative split gain increased, since larger thresholds allowed only splits that removed a larger share of parent-node impurity (Figure~\ref{fig-appendix-tree-complexity-path}). The selected trees remained small, with mean terminal-node counts between about six and eight across criteria (Table~\ref{tbl-appendix-tree-selection}, Figure~\ref{fig-appendix-tree-complexity-path}).

This avoids another round of stale exact numbers if the PSU-fold render changes the selected settings again.

### [ ] R4. Forest validation prose conflicts with appendix values and with the current analysis workflow

Location: Results, lines 1957-1985; Appendix forest tables, lines 3341-3397.

Current issue:

The forest prose says that rank-dependent forest validation gains remained negative: -25.0% for CI, -28.1% for CIc, and -34.7% for CIg. The current appendix table does not match those values. It shows raw validation gain as negative for CI but slightly positive for CIc and CIg, while validation gain/root remains negative for all three.

There is also a workflow issue: the report generator now uses PSU folds for the single trees and the ranger predictive forest, but the concentration-index forest section appears to load saved forest tuning objects. Unless those saved forest objects are regenerated using PSU fold IDs, the narrative should not say that forest validation is PSU-aware in the same way as the tree validation.

Suggested replacement if forest tuning remains from the saved non-PSU object:

> We evaluated concentration-index forests using the same validation criterion as the single trees, but these results should be read as the saved forest validation analysis rather than as the main design-aware tree validation. For the rank-dependent criteria, averaging across trees reduced the size of the validation loss compared with the single trees, but it did not produce a clearly stable rank-dependent subgroup structure. The validation-gain ratios remained unfavourable for CI, CIc, and CIg, whereas the L forest retained a positive validation-gain ratio (Figure~\ref{fig-appendix-forest-validation-path}, Table~\ref{tbl-appendix-forest-selected-setting-summary}). Thus, the forest results were consistent with the single-tree results in showing stronger stability for the level-dependent criterion.

Suggested replacement if forest tuning is rerun with PSU fold IDs:

> We evaluated concentration-index forests using the same PSU-held-out validation idea as the single trees. For the rank-dependent criteria, averaging across trees reduced the validation loss compared with the single trees, but did not fully stabilise the rank-dependent subgroup structure. The L forest retained the clearest positive validation performance (Figure~\ref{fig-appendix-forest-validation-path}, Table~\ref{tbl-appendix-forest-selected-setting-summary}). Thus, the forest results supported the single-tree finding that the level-dependent criterion produced the more stable direct partition in this application.

### [ ] R5. Surrogate forest narrative may be stale after PSU validation changes

Location: Results, lines 2003-2034 and 2038-2066.

Current issue:

The text says the selected forest surrogate retained the main structure seen in the level-dependent tree, with residence as the first split and rural-low-wealth-Katanga as the highest-mortality path. This may still be correct, but it depends on the selected forest object. Since the main generator currently loads saved forest models, this narrative must be checked after deciding whether forest tuning is rerun with PSU folds.

Suggested replacement if the surrogate still looks the same:

> To interpret the selected forest, we fitted a surrogate tree to summarise its predictions (Figure~\ref{fig-results-surrogate-tree}). The surrogate tree retained the main structure seen in the level-dependent tree. The first split was type of residence, separating urban and rural children. Among urban children, the main divisions were father's agricultural occupation and household wealth. Among rural children, the main divisions were household wealth, maternal education, and residence in Katanga.

Suggested note if the surrogate changes:

> To interpret the selected forest, we fitted a surrogate tree to summarise its predictions (Figure~\ref{fig-results-surrogate-tree}). The surrogate tree should be described from the rerendered figure rather than from the previous rural-urban narrative, because the selected forest and surrogate structure may change when validation and selection objects are regenerated.

### [ ] R6. Results section should mention bootstrap only after bootstrap results are rendered

Location: Results direct-method section, currently no bootstrap subsection.

Current issue:

The separate `report/design_aware_validation.qmd` is now bootstrap-only. The current thesis Results section does not include the bootstrap stability results. That is acceptable if the bootstrap is treated as supplementary, but if it is part of the response to the lecturer's design-aware validation concern, the Results should include a short subsection after tree validation.

Suggested new subsection after Tree complexity plots:

> \subsubsection{Bootstrap stability of selected trees}
>
> We also refitted the selected tree settings under bootstrap resampling to assess the stability of the fitted subgroup structures after model selection. The bootstrap did not retune the trees. Instead, it used the selected settings from the main generator and refitted each criterion under resampled data. The results are summarised by the frequency of successful refits, the most common root split, the most common split variables, terminal-node counts, gain distributions, and similarity of terminal-node membership to the selected reference tree. These results are reported in [insert bootstrap table/appendix reference after rendering].

Suggested interpretive sentence after the bootstrap table:

> These bootstrap results should be read as stability checks for the selected tree structures, not as p-values or formal design-based confidence intervals for individual splits.

### [ ] R7. Absolute impurity and signed inequality direction need clearer reporting

Location: Direct tree results and tree captions/labels.

Current issue:

The tree impurity uses absolute concentration-index values, which is appropriate for tree splitting. However, the results narrative should not hide the direction of inequality. A terminal node with mortality concentrated among poorer children and a terminal node with mortality concentrated among richer children can have similar absolute impurity but different substantive meanings.

Suggested addition to Direct methods section:

> The tree-growing criterion uses the absolute value of the selected concentration index as an impurity measure. This means that a split is rewarded when it reduces the magnitude of within-node inequality. For interpretation, however, the sign of the concentration index still matters. Positive and negative terminal-node values describe different directions of inequality, even when their absolute values are similar. The terminal-node summaries should therefore be read together with the signed concentration-index values reported in the tables.

### [ ] R8. Surrogate forest language should not imply calibrated mortality rates unless that is exactly what is plotted

Location: Results, lines 2003-2034.

Current issue:

The comments note that the forest is trained to reduce inequality impurity, and the surrogate tree is fitted to forest scores. The current Results language sometimes reads as if the surrogate leaves are observed mortality rates or calibrated risk estimates. That is too strong unless the plotted values are explicitly observed rates.

Suggested replacement:

> To interpret the selected forest, we fitted a surrogate tree to summarise its forest score (Figure~\ref{fig-results-surrogate-tree}). The surrogate tree should be read as an interpretable approximation to the selected forest's fitted score, not as a calibrated mortality model. Where terminal-node observed mortality is reported, it describes the observed children falling into the surrogate leaves; where the forest score is reported, it describes the fitted forest output.

## Discussion Section

### [ ] D0. Main finding should be framed as descriptive and stability-limited

Location: Discussion opening, lines 2172-2188.

Current issue:

The comments make clear that the central claim should be careful. The results support a descriptive method that finds plausible subgroup structure, but the rank-dependent trees fail validation and the L result is limited by the DHS wealth-index scale.

Suggested framing sentence to add near the start of Discussion:

> The results should therefore be read as evidence that concentration-index trees are useful for descriptive subgroup discovery in inequality analysis, not as evidence that the fitted rank-dependent trees provide a stable inferential model for DHS data. In this application, the rank-dependent trees found plausible in-sample subgroups but validated poorly, while the L tree was more stable but should be treated as a sensitivity analysis because the DHS wealth index is not a ratio-scale socioeconomic variable.

### [ ] D1. SHAP discussion repeats the wrong "about 12%" statement

Location: Discussion, lines 2190-2206.

Current issue:

The text again says residence and household wealth each contributed about 12%. This conflicts with the table values in the Results.

Suggested replacement:

> The SHAP analysis, from the fitted random forest, gave a less concentrated picture. Residence was the largest SHAP contributor, while household wealth, parental occupation, mother's education, skilled birth attendance, and some regional indicators also appeared across criteria. SHAP-based decomposition therefore adds to the regression methods by allowing a more flexible fitted model, such as a random forest, to be summarised at the determinant level. These SHAP results should be read as an explanation of the fitted predictive model rather than as causal or population-level determinant contributions.

### [ ] D2. Rank-dependent tree validation discussion uses old exact values

Location: Discussion, lines 2208-2222.

Current issue:

The discussion says rank-based methods removed about 87% to 88% of root-node inequality in training, but the current appendix table shows about 82% to 83%. The sign conclusion is still consistent, but the numbers are stale.

Suggested replacement:

> The rank-dependent trees based on CI, CIg, and CIc identified a structure involving maternal education, father's agricultural occupation, skilled birth attendance, household wealth, and Katanga. Mortality was lower among children whose mothers had some education and whose fathers were not in agriculture, and higher among children whose mothers had no education, who lived in low-wealth households, and who resided in Katanga. However, these partitions did not validate well. The rank-based methods retained large training gains, but their validation gains were negative under the held-out validation criterion, suggesting that the strong training structure did not transfer well to held-out PSUs.

### [ ] D3. Forest validation claim is too strong unless forests are rerun with PSU folds

Location: Discussion, lines 2241-2254.

Current issue:

The paragraph says the same pattern was seen in the forests and that rank-dependent forest validation gains remained below zero. This conflicts with the current forest tables if raw validation gains are read literally, and it may also overstate design-aware validation if the forest objects were not regenerated with PSU folds.

Suggested replacement:

> The differences also showed up in cross-validation. The level-dependent tree retained positive held-out validation performance, whereas the rank-dependent trees had strong training performance but negative held-out gains. The forest results were directionally similar, although they should be interpreted according to the forest validation objects used in the final render. Averaging across trees reduced the validation loss for the rank-dependent criteria, but the clearest positive validation performance remained with the level-dependent forest. Overall, the validation results suggest that CI, CIg, and CIc identified detailed in-sample subgroup rules, whereas L produced the more reproducible subgroup structure in this application.

### [ ] D4. Interpretation of L is still too confident given the wealth-index limitation

Location: Discussion, lines 2256-2274 and 2294-2304.

Current issue:

The text says the cross-validation results show that L yielded a structure transferable to held-out samples using the DHS wealth score. This is probably acceptable as a result statement, but it should be softened because the same section acknowledges that L is clearer with a ratio-scale socioeconomic variable.

Suggested replacement:

> The cross-validation results suggest that the L concentration-index impurity measure yielded the most transferable structure in this application, even though the DHS wealth index is an asset-based relative score rather than a ratio-scale socioeconomic measure. For this reason, the L result is treated as a useful sensitivity and descriptive partition, rather than as a final claim about level-dependent inequality on a ratio-scale socioeconomic variable.

### [ ] D5. Tree instability limitation should mention the new bootstrap-only notebook

Location: Discussion, lines 2306-2314.

Current issue:

The text says instability was handled partly through cross-validation and splitting controls. Since `report/design_aware_validation.qmd` is now specifically a bootstrap stability analysis, the discussion should mention bootstrap after those results are generated.

Suggested replacement:

> Decision trees are known to be unstable: small changes in the data can produce different split variables, split points, and tree structures. This comes from the recursive and greedy nature of tree construction. Once an early split changes, all later splits are fitted to different subsets of the data. In this thesis, this was handled partly through PSU-held-out cross-validation, splitting controls such as minimum node size, minimum child-node proportion, maximum depth, minimum gain, and minimum relative gain, and selected-tree bootstrap stability checks. These steps help describe how stable the selected subgroup structure is, but they do not turn the fitted trees into formal design-based inferential models.

### [ ] D5b. Survey-design limitation should be explicit about what was and was not done

Location: Methodological limitations, around lines 2306-2314.

Current issue:

The comments ask for clear separation between descriptive weighted analysis and formal survey-design inference. The current thesis mentions survey design in places, but the limitation should say plainly that the tree methods use weights and PSU-aware validation/stability checks, while not providing full replicate-weight or design-based inference for every tree split, forest, or SHAP contribution.

Suggested limitation paragraph:

> The survey design is handled only partly in the tree-based analysis. DHS sampling weights are used as case weights, and the final validation uses PSU-held-out folds to avoid placing children from the same cluster in both training and validation data. The selected-tree bootstrap also resamples PSUs to assess cluster-level stability. These steps improve the design-awareness of the descriptive analysis, but they are not the same as full survey-design inference. The thesis does not provide replicate-weight-based standard errors or formal design-based confidence intervals for individual tree splits, forest structures, or SHAP contributions. The tree and forest results should therefore be interpreted as weighted descriptive and stability-checked subgroup analyses.

### [ ] D6. SHAP limitations should be updated after PSU ranger folds and fastshap

Location: Discussion, lines 2326-2346.

Current issue:

The conceptual limitation is good, but the quoted ROC AUC around 0.61 may change after the ranger model is rerendered with PSU-held-out folds. Also, the methodology currently says shapr in the TeX, while the active generator now uses fastshap.

Suggested replacement:

> Another limitation comes from SHAP-based decomposition. This analysis was based on a predictive random forest and should be interpreted as an explanation of that fitted prediction model, not as a decomposition of population mortality risk. SHAP values depend on the fitted model and on the SHAP convention used. This is important in DHS data because maternal education, occupation, residence, and household wealth are correlated. The allocation of predictive contribution across these variables is therefore not unique, and the SHAP results should not be interpreted as causal or determinant contributions. The predictive forest was also affected by undersampling, so the fitted values should be treated as prediction scores rather than calibrated mortality probabilities. The final predictive performance statement should be updated from the PSU-held-out ranger cross-validation table after rerendering.

### [ ] D7. Simulation discussion should be reduced to a sanity-check claim

Location: Simulation results and Discussion/Future research.

Current issue:

The comments say the simulation study is useful but too small to validate the method scientifically. The thesis should not overclaim from one null and one positive-control simulation.

Suggested Discussion wording:

> The simulation study should be read as a sanity check rather than a full validation study. The null example showed that the tree-growing rule did not force a split when no inequality-relevant structure was present, and the positive-control example showed that the code could recover a planted subgroup pattern in a simple setting. A full methods validation would require repeated simulations reporting recovery probability, false-split rate, validation gain, bias in impurity reduction, and sensitivity to outcome prevalence, sample size, and survey weights.

## Abstract

### [ ] A1. Abstract direct-tree narrative may be stale after PSU-fold rerender

Location: Abstract, lines 256-268.

Current issue:

The abstract says the rank-dependent pattern did not carry well to validation folds, while the L tree and forest both kept positive validation gains. The tree part likely remains directionally consistent, but the exact forest claim should be checked because the forest validation objects may not yet be PSU-fold regenerated.

Suggested replacement before fresh forest rerender:

> The direct trees added the subgroup part of the analysis. The rank-dependent trees based on the standard concentration index (CI), generalised concentration index (CIg), and corrected concentration index (CIc) split mainly through mother's education, father's agricultural occupation, skilled birth attendance, household wealth, and Katanga. The high-mortality branch was made up of children whose mothers had no education, who lived in low-wealth households, and who resided in Katanga. This pattern was clear in the fitted data, but it did not transfer well to held-out validation folds. The level-dependent concentration index (L) method gave a different tree, beginning with residence and then splitting by occupation, wealth, and province. Because the DHS wealth index is not ratio-scale, the L result was treated as a sensitivity analysis and descriptive partition rather than as the main definitive result.

Suggested replacement after confirming both tree and forest with PSU folds:

> The direct trees added the subgroup part of the analysis. The rank-dependent trees based on the standard concentration index (CI), generalised concentration index (CIg), and corrected concentration index (CIc) split mainly through mother's education, father's agricultural occupation, skilled birth attendance, household wealth, and Katanga. The high-mortality branch was made up of children whose mothers had no education, who lived in low-wealth households, and who resided in Katanga. This pattern was clear in the fitted data, but it did not transfer well to PSU-held-out validation folds. The level-dependent concentration index (L) method gave a different tree, beginning with residence and then splitting by occupation, wealth, and province. The L-based direct models retained the strongest held-out validation performance, but because the DHS wealth index is not ratio-scale, this result was treated as a stable sensitivity analysis rather than as the main definitive result.

### [ ] A2. Abstract SHAP wording should stay broad, not "about 12%"

Location: Abstract, lines 249-254.

Current issue:

The abstract does not contain the wrong 12% number, which is good. It says the SHAP decomposition was less concentrated and lists the main variables. This is broadly consistent with the table. If the SHAP table changes after PSU ranger folds and fastshap rerendering, update the variable list.

Suggested replacement if the current SHAP ranking remains similar:

> The SHAP decomposition was less concentrated on the regression-leading determinants. Residence was the strongest SHAP contributor, while household wealth, parental occupation, mother's education, skilled birth attendance, and regional indicators also appeared across criteria.

## Cross-Cutting Items To Fix Before Final Render

### [ ] C1. Methodology still says shapr, but the active generator uses fastshap

Location: Methodology, lines 1627-1629.

Current text says shapr was used. The active generator has moved back to fastshap.

Suggested replacement:

> We then used the fastshap R package to estimate approximate SHAP contributions from the fitted forest (\citeproc{ref-fastshap2024}{Greenwell, 2024}). The SHAP analysis was used to describe how the fitted forest used the predictors, not to estimate causal or population-level contributions of predictors to under-five mortality.

### [ ] C1b. Regression decomposition clarification has been addressed in code but should remain explicit in text

Location: Methodology, lines 1596-1616.

Current status:

The generator now fits the survey-weighted logistic regression with the DHS design object and uses the margins package to obtain average marginal effects on the response scale before passing the decomposition inputs to Rineq. This directly responds to the comment about logit-scale coefficients versus probability-scale decomposition.

Suggested wording to keep:

> We then used the margins R package to obtain average marginal effects from the fitted survey-weighted model. These average marginal effects express the average change in predicted under-five mortality on the response scale. We then passed these marginal effects, together with the model matrix, the DHS wealth ranking variable, the under-five mortality outcome, and the sampling weights, into the Rineq R package to compute the concentration-index decomposition.

### [ ] C2. Methodology says trees and forests used PSU-level validation, but forest generation must be checked

Location: Methodology, lines 1650-1661.

Current issue:

The paragraph says concentration-index trees and forests were validated using PSU-level folds. The generator has been changed for CI trees and ranger predictive forest folds. The concentration-index forest section should be checked carefully: if forest tuning still loads saved `forest_tuning.rdata` that was not regenerated with PSU folds, the thesis should not claim PSU-level validation for concentration-index forests.

Suggested safe wording:

> To validate the concentration-index trees, we used cross-validation that respected the clustered structure of the DHS survey. The data are collected within primary sampling units (PSUs), so children within the same PSU may be more similar to each other than children from different PSUs. If folds are created at the child level, children from the same PSU can appear in both the training and validation data. This can give an overly optimistic view of how well a fitted tree structure transfers, because the validation data are not fully separated from the training data (\citeproc{ref-roberts2017crossvalidation}{Roberts et al., 2017}). For this reason, the main tree validation folds were created at the PSU level. The predictive random forest was evaluated using the same PSU-held-out idea. The concentration-index forest validation should be described according to the final forest tuning object used in the render.

### [ ] C3. Appendix table captions should reflect PSU-held-out validation after rerender

Locations:

- Ranger CV caption, lines 3044-3045.
- Tree selected-setting caption, lines 3183-3184.
- Forest selected-setting caption, lines 3341-3342.

Suggested captions:

- "PSU-held-out cross-validation results for the predictive ranger model"
- "Tree summary for the selected PSU-held-out validation settings"
- If forest is rerun with PSU folds: "Forest summary for the selected PSU-held-out validation settings"
- If forest is not rerun with PSU folds: "Forest summary for the saved selected validation settings"

### [ ] C4. Table captions should distinguish final fitted tree size from average validation-fold tree size

Current issue:

The comments note that terminal nodes are sometimes non-integers. A final fitted tree has an integer number of terminal nodes; non-integer values are averages across validation folds or forest component trees.

Suggested captions/labels:

- Use "mean terminal nodes across validation folds" for tuning summaries.
- Use "terminal nodes in the final fitted tree" for final fitted-tree summaries.
- Use "mean terminal nodes per tree" for forest summaries.

## Priority Order

1. Rerender `report/Generate_DRC_Results_Objects.qmd` so tree and ranger results reflect PSU fold IDs and fastshap.
2. Decide whether concentration-index forest tuning will also be rerun with PSU fold IDs. If not, do not describe forest validation as design-aware.
3. Update Results lines 1701-1715 immediately because the SHAP sentence already conflicts with the table.
4. Update Results lines 1867-1985 after the fresh render because the validation numbers are stale and partly conflict with current appendix tables.
5. Update Discussion lines 2190-2254 to match the revised Results language.
6. Add bootstrap stability results after the bootstrap render, but describe them as stability checks rather than formal inference.
7. Update Abstract last, after the final Results and Discussion narrative is settled.
