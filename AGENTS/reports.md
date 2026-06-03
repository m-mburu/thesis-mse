# Report Writing Guide

This guide describes how new report sections should be written so that they
match the tone, phrasing, and story used in `report/methodology.qmd`.

## Overall Voice

Write in a clear thesis style. The language should be explanatory, direct, and
grounded in public-health interpretation. The report should not sound like a
software manual or a machine-learning benchmark. It should explain why the
method matters for health inequality analysis and then show what the method
does in the DRC DHS application.

Prefer sentences that explain the purpose of a method before giving technical
details. The methodology often uses this pattern:

> The purpose of the inequality tree is to split the sample into subgroups
> where under-five mortality is no longer strongly concentrated by
> socioeconomic position.

Use this kind of phrasing throughout the report: first state the purpose, then
describe the statistical object, then explain how it should be interpreted.

## Storyline

Each section should follow the same broad story:

1. Begin with the public-health problem.
   Under-five mortality is not only a national average. It varies by household
   poverty, place of residence, maternal background, and province. The report
   should keep returning to the idea that subgroup structure matters because
   interventions are organised in real settings.

2. Introduce socioeconomic ranking.
   Explain that the analysis orders children or households from poorer to
   richer using the DHS wealth index. The concentration index measures how the
   health outcome is distributed along this ranking.

3. Explain the limitation of average effects.
   Classical decomposition gives determinant-level contributions, but it may
   miss thresholds, interactions, and subgroup-specific effects. Use language
   such as:

   > An important structure can be missed if mortality risk varies across
   > subgroups, if determinants act jointly, or if the relationship differs
   > above and below certain thresholds.

4. Explain the tree method as a subgroup-finding tool.
   The tree is not presented mainly as a predictive model. It is presented as a
   way to find readable subgroup rules. A split is selected because it reduces
   within-node socioeconomic inequality in the outcome.

5. Interpret results in terms of groups.
   When discussing figures or tables, say what subgroup structure is visible,
   which variables define the branches, and whether terminal nodes have lower
   remaining concentration-index values.

## Preferred Phrasing

Use phrases like:

- socioeconomic inequality in under-five mortality
- ordered from poorer to richer households
- subgroup structure
- observed covariates define each branch
- terminal node
- within-node socioeconomic inequality
- concentration-index impurity
- root-node inequality
- validation gain
- fitted mortality risk
- inequality-structured mortality score
- readable subgroup rules
- planted subgroup structure, for simulation
- positive-control simulation and negative-control/null simulation

When explaining the tree, prefer:

> The split is selected because it gives the largest reduction in
> socioeconomic inequality in under-five mortality.

Do not write as if the tree is only minimising prediction error. The report
should repeatedly distinguish the proposed method from ordinary prediction
trees:

> The outcome level in each subgroup is still reported, but it is not the
> quantity used to choose the split.

## Tone And Sentence Style

Use calm explanatory sentences. Avoid exaggerated claims. The methodology
rarely says that a method is "best" or "superior." It says what the method is
designed to reveal and what limitation it addresses.

Good style:

> The fitted tree then produces decision rules that describe population groups.
> These rules are easier to discuss in public health work than regression
> coefficients alone, because they point to groups that policy makers,
> programme teams, and public health planners can recognise and act on.

Avoid:

- "The model discovers the most important predictors."
- "The algorithm optimises performance."
- "The tree predicts mortality."
- "The method proves that the subgroup is causal."

Prefer:

- "The fitted tree identifies a subgroup structure."
- "The split reduces within-node concentration-index impurity."
- "The terminal nodes can be compared by their outcome levels and remaining
  inequality."
- "The result should be interpreted as an inequality-oriented subgroup
  summary, not as a causal effect."

## Methods Sections

When writing methods, move from simple ideas to formal notation.

A good methods paragraph should usually have this order:

1. Explain the intuitive problem.
2. Define the statistical quantity.
3. Give the formula if needed.
4. Explain each term in words.
5. Say how the quantity is used in the tree or decomposition.

For example, when introducing an impurity measure:

- First say that a tree needs a criterion for comparing candidate splits.
- Then say that ordinary CART uses outcome homogeneity.
- Then say that this thesis replaces that impurity with a concentration-index
  impurity.
- Then define the formula.
- Then explain that values close to zero mean little remaining socioeconomic
  inequality inside the subgroup.

Do not introduce formulas without interpretation. Every formula should be
followed by a short explanation of what the terms mean and why the expression
matters for the method.

## Application Sections

When writing about the DRC DHS application, keep the application connected to
the methodology. Do not only describe outputs.

A good application paragraph should answer:

- What object is being shown? A table, tree, forest surrogate, SHAP summary, or
  simulation check.
- What role does it play in the argument?
- What subgroup structure or determinant pattern does it show?
- Does it agree with or differ from the classical decomposition?
- How should the result be interpreted for inequality analysis?

Use DRC consistently for the Democratic Republic of the Congo, unless the
surrounding document explicitly uses RDC.

Use "under-five mortality" for the outcome in prose. Use `U5` only in compact
figure labels where space matters.

When referring to the ranking variable, write:

> The DHS household wealth index serves as the socioeconomic ranking variable.

When referring to weights, write:

> DHS sampling weights are rescaled and used as non-negative case weights.

## Results Writing

Results should sound like the analysis has been completed and the report is
now walking the reader through what was found. The voice should be observed,
interpretive, and grounded in the displayed evidence. It should not sound like
the methods section, where the main task is to justify why an analysis was
chosen.

Use this order for most results paragraphs:

1. State the result plainly.
   Begin with what the analysis showed, not what the output was designed to
   establish.

2. Point to the evidence.
   Refer to the figure, table, or appendix output that supports the statement.
   Do not repeat all values in the table or describe every visual detail.

3. Interpret the pattern.
   Explain what the observed pattern means for socioeconomic inequality,
   subgroup structure, or the DRC DHS application.

4. Connect to the next result.
   End by showing why the next model, table, tree, or appendix check matters.

Results should not repeat the full algorithm. They should remind the reader
only of the interpretation needed for the figure or table. A reader should be
able to understand the substantive pattern without being taken back through the
technical procedure.

When writing from figures and tables, use both main-text and appendix outputs
to tell one coherent story. Do not repeat the table or graph verbatim. Mention
the pattern that matters for the argument, then explain why that pattern is
relevant for the method. For example, an exploratory paragraph can use the main
outcome and wealth figure to introduce the analysis sample, then use appendix
baseline and province outputs to explain why subgroup rules are needed.

For results-section prose, write in an evidence-led style. Start with what the
analysis showed, then point to the figure or table, then state the modelling or
substantive implication. Avoid opening a results paragraph by saying what the
outputs are intended to "establish"; that sounds like methodology. Connect the
observed pattern to the next analysis step. This keeps the paragraph in the
results voice while still telling the reader why the output matters.

Good results openings:

> The exploratory analysis showed ...

> The selected tree split the sample first by ...

> The terminal nodes indicate that ...

> The appendix outputs show the same pattern in greater detail ...

Avoid methods-style openings in the results section:

- "This figure is used to establish ..."
- "The purpose of this output is ..."
- "This section first demonstrates why ..."
- "We fit this model because ..."

Prefer results-style versions:

- "The figure shows ..."
- "The table indicates ..."
- "The selected tree separates ..."
- "The appendix summaries add ..."

When a number is already visible in a table or plot, do not copy the table into
the prose. Use the number only when it is needed to anchor the interpretation.
Otherwise, write the direction and meaning of the result: which group is higher
or lower, where inequality remains, which split comes first, or which criterion
is most stable.

The tone should remain cautious. Results can say that a pattern is visible,
consistent, concentrated, or structured. They should not say that the analysis
proves causality, discovers the true mechanism, or identifies an intervention
effect unless the design supports that claim.

For a concentration-index tree, write around the figure like this:

> The tree shows the subgroup structure selected by the inequality-based
> splitting rule. Each terminal node represents a subgroup with its own outcome
> level and remaining within-node socioeconomic inequality.

For a selected tree, discuss:

- the first split, because early splits determine which observations are
  available for later splits;
- the variables used in the branches;
- whether the terminal nodes are interpretable as public-health groups;
- the remaining concentration-index value in the terminal nodes.

For tuning results, use the language of validation gain:

> Cross-validation is used to select the stopping-rule settings that remove the
> largest share of root-node impurity. A high validation gain indicates that the
> partition reduces inequality in held-out folds, while a negative or near-zero
> gain suggests that the split structure has not transferred well.

For simulation results, use a control-check framing:

> The null simulation is a negative-control setting. Wealth, the outcome, and
> subgroup variables are generated independently, so a meaningful tree should
> not appear. The positive-control simulation plants a subgroup mechanism, so
> the fitted tree can be checked against a known structure.

## SHAP And Forest Writing

When discussing the forest, do not present it as a black-box prediction model.
Use the methodology language:

> Since the trees are grown to reduce concentration-index impurity, not to
> minimise prediction error, we do not read the forest output as a fully
> calibrated probability of under-five death. It is a mortality-scale score
> obtained by averaging terminal-node mortality values across many
> inequality-based trees.

When discussing the surrogate tree:

> The surrogate tree does not replace the forest. Its role is to give a smaller
> set of rules that approximates the main subgroup structure learned by the
> forest.

When discussing SHAP:

> SHAP rewrites the fitted tree-based prediction as a baseline value plus
> feature contributions. This gives the additive form needed for decomposition,
> while allowing the role of a determinant to differ across children and
> subgroups.

## Captions

Captions should be informative but not too long. They should identify the
object and the criterion when relevant.

Good captions:

- "Selected concentration-index tree for under-five mortality in the DRC."
- "Root-node inequality in under-five mortality by concentration-index
  criterion."
- "Tree validation gain and complexity by minimum relative split gain."
- "Positive-control simulation tree fitted with the CIg criterion."

Avoid vague captions such as:

- "Model result"
- "Tree plot"
- "Simulation output"

## Tables

Tables should be introduced by their purpose. Do not simply say that a table
"shows results." Say what comparison the table enables.

For method-comparison tables:

> The table compares determinant-level results from the classical decomposition
> and the tree-based inequality method. The purpose is to assess whether
> variables that explain average differences in mortality are also the variables
> that structure socioeconomic inequality in mortality.

For tuning tables:

> The table reports the selected stopping-rule settings within each
> concentration-index criterion. These settings are selected using validation
> gain, so they indicate which tree complexity transferred best to held-out
> folds.

## Things To Avoid

Avoid language that implies causal identification unless the section has a
causal design. The DHS application is an inequality and subgroup analysis, not
an intervention study.

Avoid saying:

- "This variable causes mortality."
- "The tree proves that this group should be targeted."
- "The forest predicts individual mortality risk."
- "The SHAP value is the causal contribution."

Use instead:

- "This variable helps define the subgroup structure."
- "The terminal node points to a group with higher fitted mortality and
  remaining socioeconomic inequality."
- "The result is useful for describing inequality-relevant groups."
- "The SHAP contribution describes how the fitted model distributes inequality
  across determinants."

## Minimal Section Template

Use this structure when drafting new report text:

1. Context sentence:
   Explain why the output matters for inequality analysis.

2. Method sentence:
   State how the object was fitted or computed.

3. Interpretation sentence:
   Explain what the reader should look for in the table or figure.

4. Application sentence:
   Connect the result to DRC DHS subgroups, determinants, or the simulation
   control setting.

Example:

> Figure X shows the concentration-index tree fitted to the DRC DHS analysis
> sample. The tree uses the DHS wealth index as the socioeconomic ranking
> variable and selects splits that reduce within-node inequality in
> under-five mortality. Each terminal node can therefore be read as a subgroup
> with its own mortality level and remaining concentration-index value. The
> figure is used to assess whether the proposed method produces subgroup rules
> that are interpretable in the DRC public-health setting.
