# Landau_cutoff_for_mean

Plot lambda_max for lambda from -3 to 8 in the standard Landau function and fit.
Using TMath::Landau and TMath::VavilovAccurate classes of the ROOT.

## Run
You need the ROOT package.
I am using ROOT 6.26/06.
'''
root -l -b -q comparison.C
'''

## Functions

### void Landau_Vavilov_comparison()

Compare Vavilov and Landau functions in the TMath class to confirm that TMath::Landau(x, 0, 1) is the standard Landau function (location parameter = 0, scale parameter = &pi;/2).

![output](Vavilov_Landau_comparison.pdf)

### void landau_mean_distribution()

Plot mean of the standard Landau function as a function of &lambda;.

![output](lambda_mean_distribution.pdf)

### void lambda_max_formula_comparison()

Plot &lambda;<sub>max</sub> as a function of &lambda;, where &lambda;<sub>max</sub> is cutoff value that gives the &lambda; as mean of the standard Landau function.

![output](lambda_max_formula_comparison.pdf)
