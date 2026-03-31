# SampleSizeCalculation
This project compares sample size requirements across three clinical trial designs (superiority, non-inferiority, and equivalence) and three endpoint types (continuous, binary, and time-to-event), using glycaemic control in type 2 diabetes as a case study. All calculations are implemented as R functions based on standard formulas — the two-sample t-test approach for continuous endpoints, the risk difference method for binary endpoints, and the Schoenfeld formula for time-to-event endpoints. The document illustrates how the choice of endpoint and trial design jointly determine the statistical assumptions, test structure, and required sample size.

This repository contains following files:
1) sample_size_function.R - file containing only function code with minial code defining each function
2) T2D.qmd - main project folder containing fuctions, caclulated case study and resoning.
3) T2D.html - html rendered main project file
4) Endpoint Choice and Sample Size Relationship.pdf <- pdf file with the final version of case study.
