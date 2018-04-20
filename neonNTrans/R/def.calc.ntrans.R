################################################################################
#' @title NEON Soil Inorganic N Concentrations and Net N Transformation Rates

#' @author Samantha Weintraub \email{sweintraub@battelleecology.org}

#' @description Calculate soil extractable inorganic nitrogen concentrations and
#' net N transformation rates for NEON L1 data. Can use the
#' neonUtilities package (zipsByProduct and stackByTable) to download and stack monthly files prior to running this function.

#' @param kclInt A data frame containing soil masses and kcl volumes used in kcl
#'   extractions. Data product table name is ntr_internalLab
#' @param kclIntBlank A data frame containing information needed to link kcl
#'   extraction samples to procedural blanks. Data product table name is
#'   ntr_internalLabBlanks
#' @param kclExt A data frame containing inorganic N concentrations measured in
#'   kcl extractions and blanks. Data product table name is ntr_externalLab
#' @param soilMoist A data frame containing soil moisture values. Data product
#'   table name is sls_soilMoisture
#' @param dropConditions An optional list of sampleCondition or dataQF values for which to exclude
#'   ammonium and nitrate concentration measurements
#' @param dropFlagged An option to exclude ammonium and nitrate concentration measurements for
#'   samples with external lab quality flags, defaults to F
#' @param keepAll An option to keep all variables and blank info used in calculations. 
#'   Defaults to F, meaning only sample information (not blanks), critical input variables, 
#'   and calculated outputs will be included in the output data frame. 
#' @return A data frame of soil inorganic N concentrations in micrograms per gram
#'   in t-initial and t-final soil samples, as well as net N transformation rates
#'   for t-final cores. Rows for blank samples are not included when keepAll = F.

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords soil, nitrogen, mineralization, nitrification

#' @examples
#' \dontrun{
#' df1 <- "path/to/data/ntr_internalLab"
#' df2 <- "path/to/data/ntr_internalLabBlanks"
#' df3 <- "path/to/data/ntr_ntr_externalLab"
#' df4 <- "path/to/data/sls_soilMoisture"
#'
#' out <- def.calc.ntrans(kclInt = df1, kclIntBlank = df2, kclExt =df3, soilMoist = df4,
#' dropConditions = c("deprecatedMethod", "other"), dropFlagged = T)
#' }

#' @seealso Currently none

#' @export

# changelog and author contributions / copyrights
#   Samantha Weintraub (2017-11-22)
#     original creation
#   Samantha Weintraub (2018-04-20)
#     minor updates to allow for multiple dropConditions
################################################################################

# Function
def.calc.ntrans <- function(kclInt,
                    kclIntBlank,
                    kclExt,
                    soilMoist,
                    dropConditions,
                    dropFlagged = FALSE,
                    keepAll = FALSE
                    ){

  # join the internal and external lab data
  combinedDF <- dplyr::left_join(kclExt, kclInt, by = "kclSampleID")

  # set N data to NA based on sample condition or dataQF values (optional)
    if(!missing(dropConditions)) {
      # conditions and NEON quality flags in external lab data
      combinedDF$kclAmmoniumNConc[combinedDF$sampleCondition.x %in% dropConditions] <- NA
      combinedDF$kclNitrateNitriteNConc[combinedDF$sampleCondition.x %in% dropConditions] <- NA
      combinedDF$kclAmmoniumNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.x)] <- NA
      combinedDF$kclNitrateNitriteNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.x)] <- NA
      # conditions and NEON quality flags in internal lab data
      combinedDF$kclAmmoniumNConc[combinedDF$sampleCondition.y %in% dropConditions] <- NA
      combinedDF$kclNitrateNitriteNConc[combinedDF$sampleCondition.y %in% dropConditions] <- NA
      combinedDF$kclAmmoniumNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.y)] <- NA
      combinedDF$kclNitrateNitriteNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.y)] <- NA

      # how many values got set to NA with extneral lab condition filtering
      if (any(combinedDF$sampleCondition.x %in% dropConditions)) {
        num1 <-
          length(combinedDF$sampleID.x[combinedDF$sampleCondition.x %in% dropConditions])
        warning1 <-
          paste(
            'warning:',
            num1,
            'records had concentration values set to NA due to anomolous external lab sample conditions',
            sep = " "
          )
        print(warning1)
      } 
      
      # how many values got set to NA with internal lab condition filtering
      if (any(combinedDF$sampleCondition.y %in% dropConditions)) {
        num1a <-
          length(combinedDF$sampleID.x[combinedDF$sampleCondition.y %in% dropConditions])
        warning1a <-
          paste(
            'warning:',
            num1a,
            'records had concentration values set to NA due to anomolous NEON lab sample conditions',
            sep = " "
          )
        print(warning1a)
      } 
      
      # how many values got set to NA with external lab dataQF filtering
      if (any(grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.x))) {
        num2 <-
          length(combinedDF$sampleID.x[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.x)])
        warning2 <-
          paste(
            'warning:',
            num2,
            'records (including blanks) had concentration values set to NA due to data quality issues',
            sep = " "
          )
        print(warning2)
      } 
      
      # how many values got set to NA with internal lab dataQF filtering
      if (any(grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.y))) {
        num2a <-
          length(combinedDF$sampleID.x[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.y)])
        warning2a <-
          paste(
            'warning:',
            num2a,
            'records had concentration values set to NA due to data quality issues',
            sep = " "
          )
        print(warning2a)
      } 
      } else {
      combinedDF
      }
    
  # set concentration data to NA based on ammonium or nitrate quality flags (optional)
    if(dropFlagged) {
      combinedDF$kclAmmoniumNConc[combinedDF$ammoniumNQF %in% c("1", "2")] <- NA
      combinedDF$kclNitrateNitriteNConc[combinedDF$nitrateNitriteNQF %in% c("1", "2")] <- NA

      # compile list of how many ammonium values got set to NA with filtering
      if (any(combinedDF$ammoniumNQF %in% c("1", "2"))) {
        num3 <-
          length(combinedDF$kclAmmoniumNConc[combinedDF$ammoniumNQF %in% c("1", "2")])
        warning3 <-
          paste(
            'warning:',
            num3,
            'records had ammonium concentrations set to NA due to the QF value',
            sep = " "
          )
        print(warning3)
      } 

      # compile list of how many nitrate + nitrite values got set to NA with filtering
      if (any(combinedDF$nitrateNitriteNQF %in% c("1", "2"))) {
        num4 <-
          length(combinedDF$kclNitrateNitriteNConc[combinedDF$nitrateNitriteNQF %in% c("1", "2")])
        warning4 <-
          paste(
            'warning:',
            num4,
            'records had nitrate + nitrite concentrations set to NA due to the QF value',
            sep = " "
          )
        print(warning4)
      } 
      } else {
      combinedDF
      }

  # add blank info & values
  combinedDF$blank1ID <-
    kclIntBlank$kclBlank1ID[match(toupper(combinedDF$kclReferenceID), toupper(kclIntBlank$kclReferenceID))]
  combinedDF$blank1NH4 <-
    kclExt$kclAmmoniumNConc[match(toupper(combinedDF$blank1ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank1NO3 <-
    kclExt$kclNitrateNitriteNConc[match(toupper(combinedDF$blank1ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank2ID <-
    kclIntBlank$kclBlank2ID[match(toupper(combinedDF$kclReferenceID), toupper(kclIntBlank$kclReferenceID))]
  combinedDF$blank2NH4 <-
    kclExt$kclAmmoniumNConc[match(toupper(combinedDF$blank2ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank2NO3 <-
    kclExt$kclNitrateNitriteNConc[match(toupper(combinedDF$blank2ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank3ID <-
    kclIntBlank$kclBlank3ID[match(toupper(combinedDF$kclReferenceID), toupper(kclIntBlank$kclReferenceID))]
  combinedDF$blank3NH4 <-
    kclExt$kclAmmoniumNConc[match(toupper(combinedDF$blank3ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank3NO3 <-
    kclExt$kclNitrateNitriteNConc[match(toupper(combinedDF$blank3ID), toupper(kclExt$kclSampleID))]

  # calculate mean blanks
  combinedDF$blankNH4mean <-
    rowMeans(combinedDF[, c("blank1NH4", "blank2NH4", "blank3NH4")], na.rm = TRUE)
  combinedDF$blankNO3mean <-
    rowMeans(combinedDF[, c("blank1NO3", "blank2NO3", "blank3NO3")], na.rm = TRUE)

  # derive blank-corrected values
  combinedDF$kclAmmoniumNBlankCor <-
    combinedDF$kclAmmoniumNConc - combinedDF$blankNH4mean
  combinedDF$kclAmmoniumNBlankCor <-
    dplyr::if_else(combinedDF$kclAmmoniumNBlankCor < 0,
            0,
            combinedDF$kclAmmoniumNBlankCor) # set to zero if negative
  combinedDF$kclNitrateNitriteNBlankCor <-
    combinedDF$kclNitrateNitriteNConc - combinedDF$blankNO3mean
  combinedDF$kclNitrateNitriteNBlankCor <-
    dplyr::if_else(
           combinedDF$kclNitrateNitriteNBlankCor < 0,
           0,
           combinedDF$kclNitrateNitriteNBlankCor) # set to zero if negative

  # add soil moisture and dry mass fraction
  combinedDF$soilMoisture <-
    soilMoist$soilMoisture[match(combinedDF$sampleID.x, soilMoist$sampleID)]
  combinedDF$dryMassFraction <-
    soilMoist$dryMassFraction[match(combinedDF$sampleID.x, soilMoist$sampleID)]

  # count how many samples are missing moisture values
  samples <- combinedDF[!grepl("BREF", combinedDF$kclSampleID),]
  if (any(is.na(samples$soilMoisture))) {
    num5 <-
      sum(is.na(samples$soilMoisture))
    warning5 <-
      paste(
        'warning:',
        num5,
        'records were missing soil moisture values',
        sep = " "
      )
    print(warning5)
  } else {
    next()
  }
  combinedDF

  # convert concentrations to micrograms per gram soil
  combinedDF$soilDryMass <-
    combinedDF$soilFreshMass * combinedDF$dryMassFraction
  combinedDF$soilAmmoniumNugPerGram <-
    combinedDF$kclAmmoniumNBlankCor * (combinedDF$kclVolume / 1000) / combinedDF$soilDryMass * 1000
  combinedDF$soilNitrateNitriteNugPerGram <-
    combinedDF$kclNitrateNitriteNBlankCor * (combinedDF$kclVolume / 1000) / combinedDF$soilDryMass * 1000
  combinedDF$soilInorganicNugPerGram <-
    combinedDF$soilAmmoniumNugPerGram + combinedDF$soilNitrateNitriteNugPerGram
  combinedDF

  # create wide (cast) version of the df in order to calculate net rates with incubationPairID and nTransBoutType
  cast1 <-
    data.table::dcast(
      data.table::setDT(combinedDF),
      incubationPairID ~ nTransBoutType,
      value.var = c(
        "incubationLength",
        "soilInorganicNugPerGram",
        "soilNitrateNitriteNugPerGram"
      ),
      fun = mean,
      na.rm = T
    )
  cast1 <-
    dplyr::mutate(
      cast1,
      netInorganicNugPerGram = soilInorganicNugPerGram_tFinal - soilInorganicNugPerGram_tInitial,
      netNitrateNitriteNugPerGram =  soilNitrateNitriteNugPerGram_tFinal - soilNitrateNitriteNugPerGram_tInitial,
      netNminugPerGramPerDay = netInorganicNugPerGram / incubationLength_tFinal,
      netNitugPerGramPerDay = netNitrateNitriteNugPerGram / incubationLength_tFinal
    )

  # attach net rates onto combined df
  combinedDF$netNminugPerGramPerDay <-
    ifelse(combinedDF$nTransBoutType == "tInitial",
           NA,
           cast1$netNminugPerGramPerDay[match(combinedDF$incubationPairID, cast1$incubationPairID)])
  combinedDF$netNitugPerGramPerDay <-
    ifelse(combinedDF$nTransBoutType == "tInitial",
           NA,
           cast1$netNitugPerGramPerDay[match(combinedDF$incubationPairID, cast1$incubationPairID)])
  
  # set NaN to NA
  combinedDF[is.na(combinedDF)] <- NA
  
  # determine whether to keep all variable or just a subset
  if (keepAll) {
    combinedDF
  } else {
    combinedDF <- data.table::setnames(combinedDF, 
                                        old = c("plotID.x", "collectDate.x", "sampleID.x"), 
                                        new = c("plotID", "collectDate", "sampleID"))
    combinedDF <-
      subset(
        combinedDF,
        select = c(
          "plotID",
          "collectDate",
          "nTransBoutType",
          "sampleID",
          "incubationPairID",
          "incubationLength",
          "soilFreshMass",
          "dryMassFraction",
          "soilDryMass",
          "kclVolume",
          "kclAmmoniumNBlankCor",
          "kclNitrateNitriteNBlankCor",
          "soilAmmoniumNugPerGram",
          "soilNitrateNitriteNugPerGram",
          "netNminugPerGramPerDay",
          "netNitugPerGramPerDay"
        )
      )
    combinedDF <- combinedDF[!is.na(combinedDF$nTransBoutType), ]# Drop records that have no plotID (blanks)
    combinedDF <- combinedDF[order(combinedDF$incubationPairID),]
  }
  
  # Round all numeric variables to 3 digits
  combinedDF %>% mutate_if(is.numeric, round, digits=3)
}
