- Python 3
- follow PEP 8
- preconditions and postconditions for all functions (Instead of PyContracts, to reduce dependencies. I don't see a clear advantage of PyContracts aside from the possibly easier syntax and better error messages, but that becomes a liability if PyContracts is not maintained. It does not appear that PyContracts has been maintained in 3 years.)
- pytest
- logging
- Use Make rather than having people guess what script to run.
- track whether uncertainty was estimated for a particular data point
- Use sets of fake data for testing purposes.

- Issue: estimating uncertainty for Kusui's turbulence intensity. Might be better to exclude Kusui and use newer data as it comes in.

- wrappers for Sqlite database with units, uncertainty, and bounds
- Bayesian linear regression
- pipe turbulence intensity correlation
- vertical jet experiment
- jet breakup data compilation
