# bioinformatics

This code aims to discover the statistical association between mutational signatures and driver mutations.
To do this I will need to:

1. Filter the COSMIC databases we are using to remove irrelevant or bad data
2. Write programs to determine which mutations are most likely to be driver mutations
3. Format the data in a way that means I can use sigprofiler to generate mutational signatures
4. Use a logistic regression (or some other method) to calculate the correlation between mutational signatures and
   driver mutations
