# Draft response to Reviewer – Line 917

Thank you for this comment. We agree that the online calculator and the underlying research analysis code serve different purposes and should be clearly distinguished.

We have reorganised the public materials into two linked components. First, the RWDID website provides a user-facing calculator for the frozen 15-component RWDID and RWDID-C scoring algorithms. It uses fixed binary scoring intervals, training-derived beta coefficients, continuous component-scoring parameters and training P1–P99 support limits; it does not refit or re-derive the models.

Second, the public repository provides the core R analysis scripts and machine-readable parameter files used to construct food-group exposures, audit questionnaire-field mapping, derive RCS-based scoring intervals, construct continuous proportional component scores and estimate training-set beta coefficients.

Individual-level UK Biobank data cannot be redistributed because access is governed by the UK Biobank data-access agreement. The scripts and aggregate scoring parameters are provided so that researchers with authorised access to the corresponding data can reproduce the score derivation and calculation.
