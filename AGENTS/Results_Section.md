
For referencing use quarto methods


1. **RDC/DRC empirical results**
2. **Simulation results**
3. **Possible Kenya replication**

Use **DRC** if the thesis is in English. Use **RDC** only if your document uses French naming.

## Proposed Results Structure

### 4. Results
### Exploratory data analysis overall summary statistics

- a graph of the distribution of the outcome variable (e.g. under-five mortality)
- histogram of the ranking variable (e.g. wealth index)
- table of the candidate predictors (e.g. mother’s education, rural residence, province, birth order/interval)
- concentration index table by the for CI type

### 4.1 Empirical Application: DRC DHS Data

This section should show what happens when your proposed inequality-tree method is applied to the real DRC data.

#### 4.1.1 Comparison of Decomposition Methods

**Table 1: Comparison of methods**

This table should compare the main methods side by side. The lecturer’s note says:

> “Comparison methods: coefficients, contribution etc…”

So the table should include something like:

use rineq package see here scripts showing this  C:\Users\moses.mburu.FIND\Pictures\personal\mse-thesis\prev_analysis

| Variable               | Linear coefficient | Mean of variable | Concentration index of variable | Wagstaff contribution | Erreygers contribution | Shap from forestt
| ---------------------- | -----------------: | ---------------: | ------------------------------: | --------------------: | ---------------------: | 
| Mother’s education     |                    |                  |                                 |                       |                        |                 
| Rural residence        |                    |                  |                                 |                       |                        |               
| Province               |                    |                  |                                 |                       |                        |               
| Birth order / interval |                    |                  |                                 |                       |                        |               

The aim is to show whether the same determinants appear important under classical decomposition and under your tree-based approach.

The writing around this table can say:

> Table 1 compares the determinant-level results from the classical decomposition approach and the proposed tree-based inequality method. The purpose is not only to compare coefficient size, but also to assess whether variables that explain average differences in the outcome are the same variables that structure socioeconomic inequality in the outcome.

#### 4.1.2 Variable Importance Across Methods

**Table 2: Variable importance across different methods**

This table should rank variables under different approaches.

| Rank | Linear decomposition | CI-tree | CI-forest | SHAP / surrogate tree | Comment                       |
| ---: | -------------------- | ------- | --------- | --------------------- | ----------------------------- |
|    1 |                      |         |           |                       | Consistently important        |
|    2 |                      |         |           |                       | Important only in tree method |
|    3 |                      |         |           |                       |                               |

This table answers the question:

> Do the same variables matter across methods, or does the tree method reveal a different structure?

Useful wording:

> Table 2 compares variable importance across the different modelling approaches. Variables that appear consistently across methods may be interpreted as robust determinants of inequality. Variables that appear mainly in the tree-based method suggest possible interaction or subgroup-specific effects that are not well captured by the linear decomposition.

#### 4.1.3 Tree Obtained Using `ctree`

**Figure 1: Tree from `ctree`**

This should show the tree obtained using the `ctree` implementation.

Suggested caption:

> **Figure 1: Concentration-index tree for under-five mortality in the DRC using the ctree implementation.**
> The tree shows the subgroup structure selected by the inequality-based splitting rule. Each terminal node represents a subgroup with its own outcome level and within-node socioeconomic inequality.

In the text, explain:

> The first split identifies the variable that produces the largest reduction in within-node socioeconomic inequality. Subsequent splits are interpreted conditionally, because each split is made within the subgroup formed by the previous split.

Be careful with the wording. `ctree` usually refers to **conditional inference trees**, so do not write “non-conditional inference” unless that is exactly what your custom implementation does. You can write:

> ctree-based implementation

or:

> conditional-inference-tree framework

#### 4.1.4 Tree Obtained Using `rpart`

**Figure 2: Tree from `rpart`**

The lecturer’s note says:

> “Tree rpart — I will send the code to be investigated before use.”

So this figure is expected, but only after the `rpart` code is checked.

Suggested caption:

> **Figure 2: Concentration-index tree for under-five mortality in the DRC using the rpart implementation.**
> The figure provides a comparison with the ctree-based implementation and allows assessment of whether the selected subgroup structure is stable across tree-growing frameworks.

In the results text:

> Figure 2 presents the corresponding tree obtained using the rpart-based implementation. This comparison is used to assess whether the main subgroup patterns depend strongly on the tree-growing framework.

---

### 4.2 Simulation Results

This section is separate from the DRC application. Its purpose is to show whether your algorithm behaves correctly under known data-generating conditions.

#### 4.2.1 Algorithm Explanation

**Figure 3: Algorithm explained**

This is your flowchart or algorithm diagram.

Suggested caption:

> **Figure 3: Algorithm for concentration-index-based recursive partitioning.**
> At each node, socioeconomic ranks are computed, candidate splits are evaluated, and the split with the largest admissible reduction in within-node inequality is selected. The process continues until stopping rules such as minimum node size, maximum depth, and minimum gain are reached.

This figure should explain:

* start at root node;
* compute rank;
* evaluate candidate splits;
* calculate CI gain;
* apply constraints such as `minbucket`, `minprob`, `min_gain`, `min_relative_gain`;
* split if admissible;
* stop if no useful split remains.

#### 4.2.2 Trees Recovered in Simulated Data

**Figure 4: Simulation tree comparison using `ctree` and `rpart`**

The lecturer wants:

> “Tree obtained with ctree and tree obtained with rpart”

So you should show either:

* one combined figure with the two trees side by side; or
* two subfigures: Figure 4A and Figure 4B.

Suggested caption:

> **Figure 4: Trees recovered from simulated data using ctree and rpart implementations.**
> The simulation evaluates whether the proposed splitting rule can recover the expected subgroup structure when the true inequality-generating mechanism is known.

The text should say:

> The simulation results are used as a diagnostic check on the proposed method. Because the data-generating structure is known, the recovered trees can be compared against the expected subgroup pattern. Agreement between the simulated structure and the fitted tree supports the internal logic of the algorithm.

---

### 4.3 Possible Replication: Kenya DHS Data

This appears to be optional.

The lecturer’s note says:

> “Possible → Kenya to repeat what you did for RDC”

So do this only if time and data are ready.

Suggested subsection:

### 4.3 Sensitivity Application: Kenya DHS Data

Purpose:

> To assess whether the method produces interpretable subgroup structures in another country setting.

You do not need a full new thesis analysis. A smaller replication is enough:

| Output                     | Purpose                                                                       |
| -------------------------- | ----------------------------------------------------------------------------- |
| Kenya Table 1              | Compare decomposition methods                                                 |
| Kenya Table 2              | Compare variable importance                                                   |
| Kenya tree figure          | Check whether subgroup structure differs from DRC                             |
| Short comparison paragraph | Explain whether Kenya and DRC show similar or different inequality structures |

Suggested wording:

> As an additional sensitivity application, the method was applied to the Kenya DHS data using the same outcome, ranking variable, candidate predictors, and tuning rules. The aim was not to provide a full country comparison, but to assess whether the proposed method remains usable when transferred to a different DHS setting.

---

## Clean List of Deliverables

Your lecturers are asking for this:

### DRC/RDC results

1. **Table 1:** Comparison of methods
   Include coefficients, concentration indices, contributions, and possibly Wagstaff/Erreygers/tree-based results.

2. **Table 2:** Variable importance under different methods
   Compare linear decomposition, tree, forest, SHAP/surrogate tree if available.

3. **Figure 1:** Tree from `ctree`

4. **Figure 2:** Tree from `rpart`
   Only after checking the code.

### Simulation results

5. **Figure 3:** Algorithm explained
   Use the flowchart you were preparing.

6. **Figure 4:** Tree from simulated data
   Show both `ctree` and `rpart` versions.

### Optional

7. **Kenya replication**
   Repeat the DRC outputs using Kenya data if time allows.

---

## Recommended Final Results Chapter Order

Use this order:

```markdown
# 4. Results

## 4.1 Empirical Application to the DRC DHS Data

### 4.1.1 Comparison of Classical and Tree-Based Decomposition Results
Table 1

### 4.1.2 Variable Importance Across Methods
Table 2

### 4.1.3 Concentration-Index Tree Using the ctree Implementation
Figure 1

### 4.1.4 Concentration-Index Tree Using the rpart Implementation
Figure 2

## 4.2 Simulation Study


### 4.2.2 Recovery of the Simulated Subgroup Structure
Figure 4

## 4.3 Optional Country Replication: Kenya DHS Data


The key message is: **they do not only want model output; they want comparison across methods.** Your results should show whether the proposed inequality-tree approach agrees with, extends, or contradicts the classical Wagstaff/Erreygers-style decomposition.
