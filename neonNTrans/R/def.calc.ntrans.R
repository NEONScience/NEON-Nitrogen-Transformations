################################################################################
#' @title NEON Soil Inorganic N Concentrations and Net N Transformation Rates

#' @author Samantha Weintraub \email{sweintraub@battelleecology.org}

#' @description Calculate soil extractable inorganic nitrogen concentrations and
#' net N transformation rates for NEON L1 data. Recommend using the
#' neonUtilities package to download data from DP1.10086.001 prior to running this function. 
#' Function requires several data product tables as inputs in order to complete calculations.
#' 
#' @details The function calculates blank-corrected concentration data using the mean of all blanks
#' extracted with a batch of samples, then normalizes these per gram dry soil (and per day for rates).
#' Any negative blank-corrected KCl extraction data are set to 0 before proceeding with 
#' further transformations. This is robust for slightly negative values near method detection limits (0.01-0.02). However, 
#' highly negative values (< -0.02) may indicate an issue. Extensive troubleshooting by NEON suggests this was prevalent 
#' in samples collected prior to 2020 due to high nitrite in the KCl powder used for extractions. This nitrite was 
#' chemodenitrified (removed as a gas) under the acidic conditions of samples but not blanks, e.g., the very 
#' negative values are most often a contamination artifact. Instead of setting these highly negative values to 0, we recommend 
#' excluding them altogether, since there may have been N in the sample but it was not detectable given method artifacts.
#' Values with blank-corrected values < -0.02 are flagged by NEON, thus users can easily exclude them using the following arguments
#' in the function: dropAmmoniumFlags = "blanks exceed sample value", dropNitrateFlags = "blanks exceed sample value". 
#' If these arguments are not used, all negative values will be converted to 0 and included in the analyses.

#' @param kclInt A data frame containing soil masses and kcl volumes used in KCl
#'   extractions. Data product table name is ntr_internalLab
#' @param kclIntBlank A data frame containing information needed to link KCl
#'   extraction samples to procedural blanks. Data product table name is
#'   ntr_internalLabBlanks
#' @param kclExt A data frame containing inorganic N concentrations measured in
#'   KCl extractions and blanks. Data product table name is ntr_externalLab
#' @param soilMoist A data frame containing soil moisture values. Data product
#'   table name is sls_soilMoisture
#' @param dropConditions An optional list of sampleCondition or dataQF values for which to exclude
#'   ammonium and nitrate N concentration measurements (set them to NA). See categorical codes file 
#'   and data product user guide for more on these fields.
#' @param dropAmmoniumFlags An optional list of ammoniumNQF values for which to exclude
#'   ammonium N concentration measurements (set them to NA). See categorical codes file 
#'   and data product user guide for more on these fields.
#' @param dropNitrateFlags An optional list of nitrateNitriteNQF values for which to exclude
#'   nitrate+nitrite N concentration measurements (set them to NA). See categorical codes file 
#'   and data product user guide for more on these fields.

#' @return A list that contains at least 2 data frames, and up to 5 depending on specified parameters. 
#' all_data is a combined dataframe of external and internal lab measurements,
#' plus blank-corrected concentrations normalized by soil mass as well as net rates in T-final 
#' incubated samples. data_summary is a cleaned-up version of this information showing only
#' calculated variables and excluding lab metadata. Users will likely wish to use this table for downstream
#' analyses and can join on sampleID to link with other NEON soil measurements. dropped_condition is 
#' a dataframe that shows which samples had values set to NA due to dropCondition inputs (if applicable). dropped_flags
#' is a dataframe that shows which samples had values set to NA due to dropAmmoniumFlags and/or
#' dropNitrateFlags inputs (if applicable). dropped_moisture is a dataframe that shows which samples are
#' missing soil moisture, thus precluding downstream calculations.


#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

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
#' dropAmmoniumFlags = "blanks exceed sample value", 
#' dropNitrateFlags = "blanks exceed sample value" 
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
#' dropConditions = c("deprecatedMethod", "other"), 
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
################################################################################

# Function
def.calc.ntrans <- function(kclInt,
                            kclIntBlank,
                            kclExt,
                            soilMoist,
                            dropConditions,
                            dropAmmoniumFlags,
                            dropNitrateFlags
) {
  
  
  # check for missing datasets
  null.check = sapply(list(kclInt, kclIntBlank, kclExt, soilMoist), is.null)
  nullDSs = c('kclInt', 'kclIntBlank', 'kclExt', 'soilMoist')[null.check]
  if (length(nullDSs) > 0) {
    stop(paste0(paste(nullDSs, collapse = ', '), ' dataset(s) missing.'))
  }
  
  # take means of analytical reps, KCl extraction data
  kclExt <- kclExt %>%
    group_by(kclSampleID) %>%
    summarise_all(list( ~ if (is.numeric(.)) {
      mean(., na.rm = TRUE)
    } else {
      first(.)
    }))
  
  # join the internal and external lab data
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
  
  # set data values to NA based on sample condition or dataQF values (optional)
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
          'records had concentration values set to NA due to dataQF values reported by the external lab',
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
          'records had concentration values set to NA due to dataQF values reported during NEON lab processing',
          sep = " "
        )
      print(warning2a)
    } else {
      combinedDF
    }
  }
  
  # set concentration data to NA based on ammonium or nitrate quality flags (optional)
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
  
  if(!missing(dropNitrateFlags)) {
      combinedDF$kclNitrateNitriteNConc[combinedDF$nitrateNitriteNQF %in% dropNitrateFlags] <- NA
    
    # compile list of how many nitrate + nitrite values got set to NA with filtering
    if (any(combinedDF$nitrateNitriteNQF %in% dropNitrateFlags)) {
      num4 <-
        length(combinedDF$kclNitrateNitriteNConc[combinedDF$nitrateNitriteNQF %in% dropNitrateFlags])
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
  
  # Export list of samples excluded due to flags
  if(!missing(dropAmmoniumFlags) | !missing(dropNitrateFlags)){
  combinedDF_flag_dropped <- combinedDF %>% 
    filter(ammoniumNQF %in% dropAmmoniumFlags |
             nitrateNitriteNQF %in% dropNitrateFlags) %>%
  select(sampleID, kclSampleID, ammoniumNQF, nitrateNitriteNQF)
  }
  
  # add blank info & values, calculate blank-corrected conc, add mass and soil moisture, then normal per gram
  combinedDF_extras <- suppressWarnings(
    suppressMessages(
      combinedDF %>%
        mutate(
          kclReferenceID = toupper(kclReferenceID),
          incubationPairID = ifelse(is.na(sampleID), NA, substr(sampleID, 1, nchar(sampleID) - 9))
        ) %>%
        left_join(
          select(
            kclIntBlank,
            kclReferenceID,
            kclBlank1ID,
            kclBlank2ID,
            kclBlank3ID
          ),
          by = "kclReferenceID"
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
          kclAmmoniumNBlankCor = ifelse(
            as.numeric(kclAmmoniumNConc) - blankNH4mean < 0,
            0,
            as.numeric(kclAmmoniumNConc) - blankNH4mean
          ),
          kclNitrateNitriteNBlankCor = ifelse(
            as.numeric(kclNitrateNitriteNConc) - blankNO3mean < 0,
            0,
            as.numeric(kclNitrateNitriteNConc) - blankNO3mean
          )
        ) %>%
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
    )
  )
  
  # count how many samples are missing moisture values
  samples <- combinedDF_extras[!grepl("BREF", combinedDF_extras$kclSampleID),]
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
  
  # create wide (cast) version of the df in order to calculate net rates with incubationPairID and nTransBoutType
  combinedDFforCast <- combinedDF_extras %>%
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
      #fun = mean,
      na.rm = T
    )
  
  # calculate net rates
  cast1 <- cast1 %>%
    mutate(
      netInorganicNugPerGram = soilInorganicNugPerGram_tFinal - soilInorganicNugPerGram_tInitial,
      netNitrateNitriteNugPerGram =  soilNitrateNitriteNugPerGram_tFinal - soilNitrateNitriteNugPerGram_tInitial,
      netNminugPerGramPerDay = netInorganicNugPerGram / incubationLength_tFinal,
      netNitugPerGramPerDay = netNitrateNitriteNugPerGram / incubationLength_tFinal
    )
  
  # attach net rates onto combined df
  combinedDF_extras <- suppressWarnings(suppressMessages(
    combinedDF_extras %>%
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
    select(-c(ammoniumNRepNum, nitrateNitriteNRepNum))
  
  # collapse concentrations and rates onto one line for the incubation pair, samples only
  combinedDF_collapse <- combinedDF_extras %>%
    mutate(
      soilAmmoniumNugPerGram = ifelse(nTransBoutType == "tFinal", NA, soilAmmoniumNugPerGram),
      soilNitrateNitriteNugPerGram = ifelse(nTransBoutType == "tFinal", NA, soilNitrateNitriteNugPerGram),
      soilInorganicNugPerGram = ifelse(nTransBoutType == "tFinal", NA, soilInorganicNugPerGram)
    ) %>%
    filter(!sampleID == "") %>%
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
  
  # make nice summary table, samples only and key variables
  combinedDF_clean <- combinedDF_extras %>%
    filter(nTransBoutType == "tInitial") %>%
    select(sampleID, collectDate, incubationPairID) %>%
    left_join(combinedDF_collapse, by = "incubationPairID")
  
  # set NaN to NA
  combinedDF_extras[is.na(combinedDF_extras)] <- NA
  combinedDF_clean[is.na(combinedDF_clean)] <- NA
  
  # Round all numeric variables to 3 digits
  combinedDF_extras <- mutate_if(combinedDF_extras, is.numeric, round, digits = 3)
  combinedDF_clean <- mutate_if(combinedDF_clean, is.numeric, round, digits = 3)
  
  output.list <- list(
    all_data = combinedDF_extras,
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
