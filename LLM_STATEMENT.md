# Statement on the use of a large language model

A large language model (Claude, Anthropic; accessed through the claude.ai chat interface in 2026) was used during the development of this analysis code. This statement describes what it was used for, what it was not used for, and how its output was checked.

## What the model was used for

- **Restructuring code.** Adapting an earlier analysis pipeline from a related study to the present research question, removing code that belonged to other projects, and rewriting helper functions.
- **Debugging.** Interpreting R error messages reported by the authors and proposing fixes.
- **Checking variable coding against the data.** The model compared the recoding rules in the script with the value distributions in the raw data file. This identified several coding errors (see `CHANGELOG.md`, version 2.x), including an incorrect coding of sex, an unmapped religion category, a constant cohabitation variable and a PDQ-39 summary index that evaluated to zero.
- **Methodological discussion.** Discussing the consequences of the small sample size relative to the number of predictors, the handling of structurally missing ("not applicable") answers, and the choice of predictors.
- **Code review** of version 3.4 and drafting of `README.md`, `CHANGELOG.md` and this statement.
- **Drafting text.** First drafts of parts of the methods and results description, which the authors revised.

## Limitations of the model's contribution

- The model had no access to an R installation and did not run any of the analyses. All results were produced by the authors on their own computers.
- The model introduced errors of its own during development, which were found and corrected, for example a section ordering error that removed 11 participants from a descriptive table, a predictor that duplicated another one, and at least one incorrect explanation of a model output. Its suggestions were therefore not accepted without checking.
- [Authors: state here which parts of version 3 (e.g. Firth regression, noise calibration, optimism correction) were developed with or without the model.]

## Verification and responsibility

All code was read, run and checked by the authors, and recodings were verified against the raw data. The choice of research question, outcomes, primary predictors and analysis strategy, and the interpretation of the results, are the authors' own. The authors take full responsibility for the code, the analyses and the reported results.

## Data protection

[Authors: confirm and describe here. During development the raw data file was uploaded to the language model service. State whether the data were pseudonymised, whether this was covered by the study's consent and data protection approval, and the data processing terms that applied.]
