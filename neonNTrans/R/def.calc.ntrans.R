################################################################################
#' @title NEON Soil Inorganic N Concentrations and Net N Transformation Rates

#' @author Samantha Weintraub-Leff \email{sweintraub@battelleecology.org}

#' @description Calculate soil extractable inorganic nitrogen concentrations and
#' net N transformation rates for NEON L1 data. Recommend using the
#' neonUtilities package to download data from DP1.10086.001 prior to running this function. 
#' Function requires several data product tables as inputs in order to complete calculations.

#' @details The function calculates blank-corrected concentration data using the mean of all blanks
#' extracted with a batch of samples, then normalizes these per gram dry soil (and per day for net rates). 
#' 
#' Data from certain sites collected prior to 2022 exhibit frequent blank-corrected nitrate + nitrite 
#' concentrations significantly negative beyond background noise (< -0.02). This was caused by inadvertently 
#' using batches of potassium chloride (KCl) for soil extractions with high nitrite contamination. Nitrite is 
#' chemodenitrified, or abiotically converted to nitrogen gas and lost, in the presence of acidic soil 
#' but not in blanks (Homyak et al. 2015), which causes this artifact. It is possible to correct the affected 
#' data using a quantile regression, essentially 'adding back' the lost nitrite and removing apparent 
#' negatives (Weintraub-Leff at al. 2023). This correction is applied by the function to all data collected 
#' in 2021 and earlier, unless users select noxCorrection = F. Starting with the 2022 sampling season and 
#' onward, all batches of KCl are purity tested before use, greatly reducing (although not eliminating) 
#' the occurrence of blank-corrected sample concentration values < -0.02. 
#' 
#' After applying the noxCorrection (for pre-2022 data), any remaining negative 
#' blank-corrected concentration values are set to 0 by the function before proceeding with 
#' further data transformations. Concentration values that blank-correct to < -0.02 are flagged by NEON
#' in the ntr_externalLab table, thus if users would rather exclude these data points instead of using
#' the noxCorrection and/or assigning them to zero, they should specify function parameters 
#' dropAmmoniumFlags or dropNitrateFlags = "blanks exceed sample value". 

#' @param kclInt A data frame containing soil masses and kcl volumes used in KCl
#'   extractions. Data product table name is ntr_internalLab
#' @param kclIntBlank A data frame containing information needed to link KCl
#'   extraction samples to procedural blanks. Data product table name is
#'   ntr_internalLabBlanks
#' @param kclExt A data frame containing inorganic N concentrations measured in sample
#'   KCl extractions and blanks. Data product table name is ntr_externalLab
#' @param soilMoist A data frame containing soil moisture values. Data product
#'   table name is sls_soilMoisture
#' @param noxCorrection A parameter that specifies whether the quantile regression
#'   approach outlined in Weintraub-Leff et al (2023) should be used to correct for 
#'   contaminant nitrite in 2021 and earlier data, defaults to T
#' @param dropConditions An optional list of sampleCondition or dataQF values for which to exclude
#'   ammonium and nitrate + nitrite N concentration measurements (set them to NA). 
#'   See categorical codes file and data product user guide for more on these fields.
#' @param dropAmmoniumFlags An optional list of ammoniumNQF values for which to exclude
#'   ammonium N concentration measurements (set them to NA). See categorical codes file 
#'   and data product user guide for more on these fields.
#' @param dropNitrateFlags An optional list of nitrateNitriteNQF values for which to exclude
#'   nitrate+nitrite N concentration measurements (set them to NA). See categorical codes file 
#'   and data product user guide for more on these fields. If 'blanks exceed sample value'
#'   is selected along with the default noxCorrection = T, values are only dropped if
#'   they still blank-correct negative after the correction is applied (relevant 
#'   to pre-2022 data).

#' @return A list that contains at least 2 data frames, and up to 5 depending on specified parameters. 
#' all_data is a combined dataframe of external and internal lab measurements,
#' plus blank-corrected concentrations normalized by soil mass as well as net rates in T-final 
#' incubated samples. data_summary is a cleaned-up version of this information showing only
#' calculated variables and excluding lab metadata. Users will likely wish to use this table for downstream
#' analyses and can join on sampleID to link with other NEON soil measurements. dropped_condition is 
#' a dataframe that shows which samples had values set to NA due to dropCondition inputs (if applicable). 
#' dropped_flags is a dataframe that shows which samples had values set to NA due to dropAmmoniumFlags 
#' and/or dropNitrateFlags inputs (if applicable). dropped_moisture is a dataframe that shows which samples 
#' are missing soil moisture (if), thus precluding downstream calculations.

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
#' Homyak, P.M., Vasquez, K.T., Sickman, J.O., Parker, D.R. and Schimel, J.P. (2015). Improving Nitrite Analysis in Soils: Drawbacks of the Conventional 2 M KCl Extraction. Soil Science Society of America Journal, 79: 1237‐1242. doi:10.2136/sssaj2015.02.0061n
#' Weintraub-Leff, S.R., Hall, S.J., Craig, M.E., Sihi, D., Wang, Z. and Hart, S.C. (2023). Standardized data to improve understanding and modeling of soil nitrogen at continental scale. Earth's Future, 11, e2022EF003224. doi.org/10.1029/2022EF003224

#' @keywords soil, nitrogen, mineralization, nitrification

#' @examples
#' \dontrun{
#' If data are loaded to R using loadByProduct (neonUtilties)
#' 
#' soilData <- loadByProduct(
#' site = c("GUAN", "HARV", "KONZ"),
#' dpID = "DP1.10086.001", 
#' package = "basic", 
#' check.size = F)
#' 
#' out <- def.calc.ntrans(
#' kclInt = soilData$ntr_internalLab, 
#' kclIntBlank = soilData$ntr_internalLabBlanks, 
#' kclExt = soilData$ntr_externalLab, 
#' soilMoist = soilData$sls_soilMoisture, 
#' dropConditions = c("extract stored at incorrect temperature", "soil stored at incorrect temperature", 
#' "mass uncertain", "volume uncertain") 
#' )
#' 
#' # If data are downloaded to a computer, then need to be loaded manually
#' 
#' df1 <- "path/to/data/ntr_internalLab"
#' df2 <- "path/to/data/ntr_internalLabBlanks"
#' df3 <- "path/to/data/ntr_ntr_externalLab"
#' df4 <- "path/to/data/sls_soilMoisture"
#'
#' out <- def.calc.ntrans(
#' kclInt = df1, 
#' kclIntBlank = df2, 
#' kclExt =df3, 
#' soilMoist = df4,
#' noxCorrection = F, 
#' dropNitrateFlags = "blanks exceed sample value" 
#' )
#' }

#' @seealso Currently none

#' @export

# changelog and author contributions / copyrights
#   Samantha Weintraub (2017-11-22)
#     original creation
#   Samantha Weintraub (2018-04-20)
#     minor updates to allow for multiple dropConditions
#   Samantha Weintraub (2019-03-29)
#     bug fixes
#   Samantha Weintraub (2020-11-05)
#     new parameters to allow for more detailed filtering based on external lab QF fields,
#     deal with KCl extraction analytical replicates by taking mean,
#     change format of output from single dataframe to list with several tables.
#   Samantha Weintraub-Leff (2023-12-07)
#     new parameter to implement the quantile regression approach to correct for 
#     contaminant nitrite in 2021 and earlier data
################################################################################

# Function
def.calc.ntrans <- function(kclInt,
                            kclIntBlank,
                            kclExt,
                            soilMoist,
                            noxCorrection = T,
                            dropConditions,
                            dropAmmoniumFlags,
                            dropNitrateFlags
) {
  
  
  # Check for missing datasets ----
  null.check = sapply(list(kclInt, kclIntBlank, kclExt, soilMoist), is.null)
  nullDSs = c('kclInt', 'kclIntBlank', 'kclExt', 'soilMoist')[null.check]
  if (length(nullDSs) > 0) {
    stop(paste0(paste(nullDSs, collapse = ', '), ' dataset(s) missing.'))
  }
  
  # Take means of analytical reps, KCl extraction data ----
  kclExt <- kclExt %>%
    group_by(kclSampleID) %>%
    summarise_all(list( ~ if (is.numeric(.)) {
      mean(., na.rm = TRUE)
    } else {
      first(.)
    }))
  
  # Join the internal and external lab data ----
  suppressWarnings(suppressMessages(combinedDF <-
                                      left_join(
                                        kclExt,
                                        kclInt,
                                        by = c(
                                          "sampleID",
                                          "kclSampleID",
                                          "sampleCode",
                                          "kclSampleCode",
                                          "domainID",
                                          "siteID",
                                          "plotID",
                                          "namedLocation",
                                          "collectDate"
                                        ), 
                                        suffix = c(".externalLab", ".internalLab")
                                      )))
  
  # Set data values to NA based on sample condition or dataQF values (optional) ----
  if(!missing(dropConditions)) {
    # conditions and NEON quality flags in external lab data
    combinedDF$kclAmmoniumNConc[combinedDF$sampleCondition.externalLab %in% dropConditions] <- NA
    combinedDF$kclNitrateNitriteNConc[combinedDF$sampleCondition.externalLab %in% dropConditions] <- NA
    combinedDF$kclAmmoniumNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.externalLab)] <- NA
    combinedDF$kclNitrateNitriteNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.externalLab)] <- NA
    # conditions and NEON quality flags in internal lab data
    combinedDF$kclAmmoniumNConc[combinedDF$sampleCondition.internalLab %in% dropConditions] <- NA
    combinedDF$kclNitrateNitriteNConc[combinedDF$sampleCondition.internalLab %in% dropConditions] <- NA
    combinedDF$kclAmmoniumNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.internalLab)] <- NA
    combinedDF$kclNitrateNitriteNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.internalLab)] <- NA
    # Export list samples excluded due to conditions
    combinedDF_condition_dropped <- combinedDF %>% 
      filter(sampleCondition.externalLab %in% dropConditions | 
             grepl(paste(dropConditions, collapse = "|"), dataQF.externalLab) |
             sampleCondition.internalLab %in% dropConditions |
             grepl(paste(dropConditions, collapse = "|"), dataQF.internalLab)) %>%
      select(sampleID, kclSampleID, sampleCondition.externalLab, dataQF.externalLab, sampleCondition.internalLab, dataQF.internalLab)
    
    # how many values got set to NA with external lab condition filtering
    if (any(combinedDF$sampleCondition.externalLab %in% dropConditions)) {
      num1 <-
        length(combinedDF$sampleID[combinedDF$sampleCondition.externalLab %in% dropConditions])
      warning1 <-
        paste(
          'Note:',
          num1,
          'records had concentration values set to NA due to external lab sample conditions',
          sep = " "
        )
      print(warning1)
    } 
    
    # how many values got set to NA with internal lab condition filtering
    if (any(combinedDF$sampleCondition.internalLab %in% dropConditions)) {
      num1a <-
        length(combinedDF$sampleID[combinedDF$sampleCondition.internalLab %in% dropConditions])
      warning1a <-
        paste(
          'Note:',
          num1a,
          'records had concentration values set to NA due to NEON lab sample conditions',
          sep = " "
        )
      print(warning1a)
    } 
    
    # how many values got set to NA with external lab dataQF filtering
    if (any(grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.externalLab))) {
      num2 <-
        length(combinedDF$sampleID[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.externalLab)])
      warning2 <-
        paste(
          'Note:',
          num2,
          'records had concentration values set to NA due to dataQF values reported in external lab data',
          sep = " "
        )
      print(warning2)
    } 
    
    # how many values got set to NA with internal lab dataQF filtering
    if (any(grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.internalLab))) {
      num2a <-
        length(combinedDF$sampleID[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.internalLab)])
      warning2a <-
        paste(
          'Note:',
          num2a,
          'records had concentration values set to NA due to dataQF values reported in NEON lab data',
          sep = " "
        )
      print(warning2a)
    } else {
      combinedDF
    }
  }
  
  # Set concentration data to NA based on ammonium quality flags (optional) ----
  if(!missing(dropAmmoniumFlags)) {
    combinedDF$kclAmmoniumNConc[combinedDF$ammoniumNQF %in% dropAmmoniumFlags] <- NA
    
    # compile list of how many ammonium values got set to NA with filtering
    if (any(combinedDF$ammoniumNQF %in% dropAmmoniumFlags)) {
      num3 <-
        length(combinedDF$kclAmmoniumNConc[combinedDF$ammoniumNQF %in% dropAmmoniumFlags])
      warning3 <-
        paste(
          'Note:',
          num3,
          'records had ammonium concentrations set to NA due to ammoniumQF values',
          sep = " "
        )
      print(warning3)
    } 
  }
  
  # Add blank info & values, calculate blank-corrected conc including noXCorrection if relevant ----
  combinedDF_2 <- suppressWarnings(
    suppressMessages(
      combinedDF %>%
        mutate(kclReferenceID = toupper(kclReferenceID)) %>%
        left_join(
          select(
            kclIntBlank,
            kclReferenceID,
            kclBlank1ID,
            kclBlank2ID,
            kclBlank3ID), by = "kclReferenceID"
        ) %>%
        mutate(
          blank1NH4 = as.numeric(kclAmmoniumNConc[match(kclBlank1ID, kclSampleID)]),
          blank2NH4 = as.numeric(kclAmmoniumNConc[match(kclBlank2ID, kclSampleID)]),
          blank3NH4 = as.numeric(kclAmmoniumNConc[match(kclBlank3ID, kclSampleID)]),
          blank1NO3 = as.numeric(kclNitrateNitriteNConc[match(kclBlank1ID, kclSampleID)]),
          blank2NO3 = as.numeric(kclNitrateNitriteNConc[match(kclBlank2ID, kclSampleID)]),
          blank3NO3 = as.numeric(kclNitrateNitriteNConc[match(kclBlank3ID, kclSampleID)]),
          blankNH4mean = rowMeans(data.frame(blank1NH4, blank2NH4, blank3NH4), na.rm = TRUE),
          blankNO3mean = rowMeans(data.frame(blank1NO3, blank2NO3, blank3NO3), na.rm = TRUE), 
          kclAmmoniumNBlankCor = kclAmmoniumNConc - blankNH4mean))) # include NH4 up here since it's not linked to noxCorrection
      
      if(noxCorrection == T) {
        combinedDF_2 <- combinedDF_2 %>%
          mutate(kclNitrateNitriteNBlankCor = ifelse(collectDate < "2022-01-01", 
                                                     (kclNitrateNitriteNConc - blankNO3mean) + 0.77*blankNO3mean, # -0.77 is the regression slope of the 0.01 quantile for 2021 and earlier data
                                                     kclNitrateNitriteNConc - blankNO3mean), 
                 noxCorrection = ifelse(collectDate < "2022-01-01", "implemented", "not applicable"))
      }
      
      if(noxCorrection == F) {
        combinedDF_2 <- combinedDF_2 %>%
          mutate(kclNitrateNitriteNBlankCor = kclNitrateNitriteNConc - blankNO3mean,
                 noxCorrection = ifelse(collectDate < "2022-01-01", "not implemented", "not applicable"))
      }
 
    # Set concentration data to NA based on nitrate quality flags (optional) ----
    # more complex since it does retain values that blank-correct >= -0.02 after noxCorrection
    if(!missing(dropNitrateFlags)) {
        combinedDF_2 <- combinedDF_2 %>%
        mutate(kclNitrateNitriteNConc = ifelse(nitrateNitriteNQF == "blanks exceed sample value" & kclNitrateNitriteNBlankCor >= -0.02, 
                                               kclNitrateNitriteNConc, ifelse(nitrateNitriteNQF %in% dropNitrateFlags, 
                                                                              NA, kclNitrateNitriteNConc)), 
               nitrateDroppedQF = ifelse(nitrateNitriteNQF == "blanks exceed sample value" & kclNitrateNitriteNBlankCor >= -0.02, 
                                         "no", ifelse(nitrateNitriteNQF %in% dropNitrateFlags, "yes", "no")), 
               nitrateDroppedQF = ifelse(is.na(nitrateDroppedQF), "no", nitrateDroppedQF),
               kclNitrateNitriteNBlankCor = ifelse(is.na(kclNitrateNitriteNConc), NA, kclNitrateNitriteNBlankCor))
        }           
                                               

    # compile list of how many nitrate + nitrite values got set to NA with filtering
    if (sum(combinedDF_2$nitrateDroppedQF == "yes") > 0) {
      num4 <- sum(combinedDF_2$nitrateDroppedQF == "yes")
      warning4 <-
        paste(
          'Note:',
          num4,
          'records had nitrate + nitrite concentrations set to NA due to the nitrateNitriteNQF value',
          sep = " "
        )
      print(warning4)
    } 
    } 
  
  # Export list of samples excluded due to ammonium and nitrate QF values
  if(!missing(dropAmmoniumFlags) | !missing(dropNitrateFlags)){
    combinedDF_flag_dropped <- combinedDF_2 %>% 
      filter(ammoniumNQF %in% dropAmmoniumFlags |
               nitrateDroppedQF == "yes") %>%
      select(sampleID, kclSampleID, ammoniumNQF, nitrateNitriteNQF)
  }
  
  # Set remaining negatives to 0, add mass and soil moisture, then normalize per gram dry soil ----
    combinedDF_3 <- combinedDF_2 %>%
    mutate(kclAmmoniumNBlankCor = ifelse(kclAmmoniumNBlankCor < 0, 0, kclAmmoniumNBlankCor),
           kclNitrateNitriteNBlankCor = ifelse(kclNitrateNitriteNBlankCor < 0, 0, kclNitrateNitriteNBlankCor), 
           noxCorrection = ifelse(is.na(kclNitrateNitriteNBlankCor), NA, noxCorrection)) %>%
  left_join(
    select(soilMoist, sampleID, soilMoisture, dryMassFraction),
    by = "sampleID"
  ) %>%
    mutate(
      soilDryMass = soilFreshMass * dryMassFraction,
      soilAmmoniumNugPerGram = kclAmmoniumNBlankCor * (kclVolume / 1000) / soilDryMass * 1000,
      soilNitrateNitriteNugPerGram = kclNitrateNitriteNBlankCor * (kclVolume / 1000) / soilDryMass * 1000,
      # the way this var is calculated, if either NH4 or NO3 is NA, the value is NA
      soilInorganicNugPerGram = soilAmmoniumNugPerGram + soilNitrateNitriteNugPerGram
    )
  
  # Produce a table for missing moisture values ----
  samples <- combinedDF_3[!grepl("BREF", combinedDF_3$kclSampleID),]
  if (any(is.na(samples$soilMoisture))) {
    num5 <-
      sum(is.na(samples$soilMoisture))
    warning5 <-
      paste(
        'Note:',
        num5,
        'records were missing soil moisture values',
        sep = " "
      )
    print(warning5)
    
    # export DF that lists these
    combinedDF_moisture_dropped <- samples %>% 
      filter(is.na(samples$soilMoisture)) %>%
      select(sampleID, kclSampleID, soilFreshMass, soilMoisture, dryMassFraction, soilDryMass)
  }
  
  
  # Create wide (cast) version of the df in order to calculate net rates with incubationPairID and nTransBoutType ----
  combinedDFforCast <- combinedDF_3 %>%
    filter(!sampleID == "")
  
  cast1 <-
    data.table::dcast(
      data.table::setDT(combinedDFforCast),
      incubationPairID ~ nTransBoutType,
      value.var = c(
        "incubationLength",
        "soilInorganicNugPerGram",
        "soilNitrateNitriteNugPerGram"
      ),
      fun = mean, # need this in case of non-unique values, will take mean
      na.rm = T
    )
  
  # Calculate net rates ----
  cast1 <- cast1 %>%
    mutate(
      netInorganicNugPerGram = soilInorganicNugPerGram_tFinal - soilInorganicNugPerGram_tInitial,
      netNitrateNitriteNugPerGram =  soilNitrateNitriteNugPerGram_tFinal - soilNitrateNitriteNugPerGram_tInitial,
      netNminugPerGramPerDay = netInorganicNugPerGram / incubationLength_tFinal,
      netNitugPerGramPerDay = netNitrateNitriteNugPerGram / incubationLength_tFinal
    )
  
  # Attach net rates onto combined df ----
  combinedDF_3 <- suppressWarnings(suppressMessages(
    combinedDF_3 %>%
      left_join(
        select(
          cast1,
          incubationPairID,
          netNminugPerGramPerDay,
          netNitugPerGramPerDay
        ),
        by = "incubationPairID"
      ) %>%
      mutate(
        netNminugPerGramPerDay = ifelse(nTransBoutType == "tInitial", NA, netNminugPerGramPerDay),
        netNitugPerGramPerDay = ifelse(nTransBoutType == "tInitial", NA, netNitugPerGramPerDay)
      )
  )) %>%
    select(-c(ammoniumNRepNum, nitrateNitriteNRepNum, nitrateDroppedQF))
  
  # Collapse concentrations and rates onto one line for the incubation pair, samples only ----
  combinedDF_collapse <- combinedDF_3 %>%
    mutate(
      soilAmmoniumNugPerGram = ifelse(nTransBoutType == "tFinal", NA, soilAmmoniumNugPerGram), # get rid of extractable N data for t-final
      soilNitrateNitriteNugPerGram = ifelse(nTransBoutType == "tFinal", NA, soilNitrateNitriteNugPerGram),
      soilInorganicNugPerGram = ifelse(nTransBoutType == "tFinal", NA, soilInorganicNugPerGram)
    ) %>%
    filter(!sampleID == "") %>% # get rid of blanks
    select(
      incubationPairID,
      soilAmmoniumNugPerGram,
      soilNitrateNitriteNugPerGram,
      soilInorganicNugPerGram,
      netNminugPerGramPerDay,
      netNitugPerGramPerDay
    ) %>%
    group_by(incubationPairID) %>%
    summarise_all(list(~ if (is.numeric(.)) {
      mean(., na.rm = TRUE)
    } else {
      first(.)
    }))
  
  # Make nice summary table, samples only and key variables ----
  combinedDF_clean <- combinedDF_3 %>%
    filter(nTransBoutType == "tInitial") %>%
    select(sampleID, collectDate, incubationPairID) %>%
    left_join(combinedDF_collapse, by = "incubationPairID")
  
  # Final data cleaning and preparation for export list ----
  # set NaN to NA
  combinedDF_3[is.na(combinedDF_3)] <- NA
  combinedDF_clean[is.na(combinedDF_clean)] <- NA
  
  # Round all numeric variables to 3 digits
  combinedDF_3 <- mutate_if(combinedDF_3, is.numeric, round, digits = 3)
  combinedDF_clean <- mutate_if(combinedDF_clean, is.numeric, round, digits = 3)
  
  output.list <- list(
    all_data = combinedDF_3,
    data_summary = combinedDF_clean
  )
  
  if (exists("combinedDF_condition_dropped")) {
    output.list[['dropped_condition']] = combinedDF_condition_dropped
  }
  if (exists("combinedDF_flag_dropped")) {
    output.list[['dropped_flags']] = combinedDF_flag_dropped
  }
  if (exists("combinedDF_moisture_dropped")) {
    output.list[['dropped_moisture']] = combinedDF_moisture_dropped
  }
  
  return(output.list)
}
