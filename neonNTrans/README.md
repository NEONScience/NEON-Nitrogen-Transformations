neonNTrans
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- ****** Description ****** -->

This package is for calculating soil extractable inorganic nitrogen (N)
concentrations and net N transformation rates (mineralization and
nitrification) from NEON soil KCl extracts.

<!-- ****** Usage ****** -->

## Usage

The single function contained in this package, def.calc.ntrans, has the
following purpose: (1) join variables across tables (2) calculate
blank-corrected inorganic N concentrations in KCl extracts, excluding
any data with flags specified by the user (3) convert from N
concentrations in extracts (mg N/L) to soil extractable N concentrations
(ug N/g dry soil), and (4) calculate net N mineralization and
nitrification rates using initital and incubated core pairs. See the
function help file for additional details. The general flow of using
this package is:

1.  Download data for NEON.DP1.10086, Soil physical and chemical
    properties, periodic.
2.  *Recommended Workflow* Use functions from the neonUtilities package
    (on CRAN) to download files and/or read files directly into R. See
    neonUtilities package documentation for more details.

<!-- -->

     library(neonUtilities)  
     soilData <- loadByProduct(site = c("LENO", "HARV", "WREF"), dpID = "DP1.10086.001", package = "basic", check.size = F)

3.  Alternatively, download files manually from the NEON data portal, or
    via some other mechanism.
4.  Load the def.calc.ntrans package:  

<!-- -->

     library(devtools) 
     install_github("NEONScience/NEON-Nitrogen-Transformations/neonNTrans", dependencies=TRUE)  
     library(neonNTrans)  

5.  Run the function with data loaded to R using the loadByProduct()
    example above, with option to specify flagged data to be excluded.

<!-- -->

    out <- def.calc.ntrans(kclInt = soilData$ntr_internalLab, kclIntBlank = soilData$ntr_internalLabBlanks, kclExt = soilData$ntr_externalLab, soilMoist = soilData$sls_soilMoisture, dropConditions = c("extract stored at incorrect temperature", "soil stored at incorrect temperature", "mass uncertain", "volume uncertain"))

6.  Alternatively, run the function with files for each of the four
    required input tables loaded manually, with option to turn off the
    default ‘noxCorrection’ that accounts for nitrite contamination in
    KCl batches used prior to 2022 and instead drop samples with high
    nitrite in blanks.

<!-- -->

    df1 <- "path/to/data/ntr_internalLab"
    df2 <- "path/to/data/ntr_internalLabBlanks"
    df3 <- "path/to/data/ntr_ntr_externalLab"
    df4 <- "path/to/data/sls_soilMoisture"
    out <- def.calc.ntrans(kclInt = df1, kclIntBlank = df2, kclExt = df3, soilMoist = df4, noxCorrection = F, dropNitrateFlags = "blanks exceed sample value") 

Returns a list of dataframes, at least 2 are always included and up to 5
may be returned depending on specified function parameters. data_summary
will be most useful to end users as it succinctly provides inorganic N
concentrations in ug N per g dry soil plus net N transformation rates in
ug N per g per day for incubated samples if data from both initial and
final cores are available. all_data is provided so that end users can
see all of the calculations that went in to the final estimates. The
other data frames are summaries of which records were excluded due to
conditions, flags, or missing soil moisture values. Note that the
function will not run if inorganic N concentrations from KCl extractions
are not yet available from external lab analyses.

<!-- ****** Calculation Summary ****** -->

## Calculation Summary

The first step in converting inorganic N concentrations in KCl extracts
to soil extractable inorganic N is to blank-correct the concentration
values (mg/L) for each sample. This is necessary because blanks can
contain notable levels of ammonium and nitrate, and this must be
accounted for. The def.calc.ntrans function thus estimates mean N
concentration in the set of blanks associated with each batch of samples
(typically, n = 3), then subtracts the mean blank value from measured
ammonium and nitrate+nitrite N concentrations.

Sometimes this can result in a negative number, and the function sets
all negative concentrations to 0. For subtle negatives near the method
detection limit (0.01 - 0.02 mg/L), this substitution is robust. Yet
prior to 2022, data from some sites showed strongly negative
blank-corrected nitrate+nitrite N concentrations. This was caused by
inadvertently using batches of KCl with high nitrite contamination,
since nitrite is chemodenitrified, or abiotically converted to nitrogen
gas and lost, in the presence of acidic soil but not in blanks (Homyak
et al. 2015). To account for this, the function default it to correct
any pre-2022 data using a quantile regression, essentially ‘adding back’
the lost nitrite and thus removing apparent negatives, as described by
Weintraub-Leff at al. (2023). If users do not wish to apply this
correction, they can turn it off and allow those strongly negative
concentrations to be set to 0 instead or be excluded from further
calculations. See function help for more discussion of the parameters
that enable these choices.

The next step is to convert blank-corrected concentrations from
milligrams N per liter to micrograms N per gram by dividing by the mass
of dry soil extracted and multiplying by extraction volume. To calculate
the mass of dry soil extracted, the def.calc.ntrans function multiplies
the mass of field-moist soil extracted by the dry mass fraction provided
for each sample in the soil moisture table.

Lastly, if data are available for both an initial and final core in the
pair, net N mineralization is calculated as the difference in final
minus initial inorganic N (ammonium plus nitrate+nitrite N), divided by
the incubation length. Net nitrification is calculated as the difference
between final and initial nitrate+nitrite N, divided by the incubation
length. Both rates can be either positive or negative numbers.

For a more detailed breakdown of these calculations and a full list of
all variables generated by the def.calc.ntrans function, see the NEON
Data Product User Guide for Soil inorganic nitrogen pools and
transformations (NEON.DP1.10086), available at the NEON Data Portal in
the Resources \> Document Library section or the data product landing
page: <https://data.neonscience.org/data-products/DP1.10086.001>.

Note that the def.calc.ntrans function provides a summary of how many
records are missing soil concentration estimates due to lack of dry mass
fraction data in the soil moisture table. Moreover, when users decide to
drop data based on specific sample condition values or data quality
flags (see function help file for details), the function also provides a
summary of the number of records affected.

<!-- ****** References ****** -->

## References

Homyak, P.M., Vasquez, K.T., Sickman, J.O., Parker, D.R. and Schimel,
J.P. (2015). Improving Nitrite Analysis in Soils: Drawbacks of the
Conventional 2 M KCl Extraction. Soil Science Society of America
Journal, 79: 1237‐1242. <doi:10.2136/sssaj2015.02.0061n>

Weintraub-Leff, S.R., Hall, S.J., Craig, M.E., Sihi, D., Wang, Z. and
Hart, S.C. (2023). Standardized data to improve understanding and
modeling of soil nitrogen at continental scale. Earth’s Future, 11,
e2022EF003224. doi.org/10.1029/2022EF003224

<!-- ****** Acknowledgements ****** -->

## Credits & Acknowledgements

<!-- HTML tags to produce image, resize, add hyperlink. -->

<!-- ONLY WORKS WITH HTML or GITHUB documents -->

<a href="http://www.neonscience.org/">
<img src="logo.png" width="300px" /> </a>

<!-- Acknowledgements text -->

The National Ecological Observatory Network is a project solely funded
by the National Science Foundation and managed under cooperative
agreement by Battelle. Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the author(s)
and do not necessarily reflect the views of the National Science
Foundation.

<!-- ****** License ****** -->

## License

GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

<!-- ****** Disclaimer ****** -->

## Disclaimer

*Information and documents contained within this pachage are available
as-is. Codes or documents, or their use, may not be supported or
maintained under any program or service and may not be compatible with
data currently available from the NEON Data Portal.*
