################################################################################
#' @title NEON Soil Inorganic N Concentration and Net N Transformation Rates

#' @author
#' Samantha Weintraub \email{sweintraub@battelleecology.org}

#' @description
#' Calculate soil extractable inorganic nitrogen concentrations and net
#' transformation rates for NEON L1 data

#' @param kclInt A data frame containing soil masses and kcl volumes used in kcl
#'   extractions. Data product table name is ntr_internalLab
#' @param kclIntBlank A data frame containing information needed to link kcl
#'   extraction samples to procedural blanks. Data product table name is
#'   ntr_internalLabBlanks
#' @param kclExt A data frame containing inorganic N concentrations measured in
#'   kcl extractions and blanks. Data product table name is ntr_externalLab
#' @param soilMoist A data frame containing soil moisture data. Data product
#'   table name is sls_soilMoisture
#' @param toFilter An optional list of the condition values to filter on
#' @param dropFlagged Parameter to filter out concentration measurements for
#'   samples with quality flags, defaults to F
#' @param keepAll Parameter that specifies whether to keep all variables used in
#'   the calculations or only relevant one, defaults to F
#' @return A dataframe of soil inorganic N concentrations in micrograms per gram
#'   in t-initial and t-final soil cores, as well as net N transformation rates
#'   for t-final cores

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords soil, nitrogen mineralization, nitrification

#' @examples
#' \dontrun{
#' df1 <- "path/to/data/ntr_internalLab"
#' df2 <- "path/to/data/ntr_internalLabBlanks"
#' df3 <- "path/to/data/ntr_ntr_externalLab"
#' df4 <- "path/to/data/sls_soilMoisture"
#'
#' out <- soilKCl(kclInt = df1, kclIntBlank = df2, kclExt =df3, soilMoist = df4,
#' filter = T, toFilter = c("thawed & warm"), dropFlagged = T)
#' }

#' @seealso Currently none

#' @export

# changelog and author contributions / copyrights
#   Samantha Weintraub (2017-11-22)
#     original creation
################################################################################

# Function
def.calc.ntrans <- function(kclInt,
                    kclIntBlank,
                    kclExt,
                    soilMoist,
                    toFilter,
                    dropFlagged = FALSE,
                    keepAll = FALSE
                    ){

  # join the internal and external lab data
  combinedDF <- dplyr::left_join(kclExt, kclInt, by = "kclSampleID")

  # set N data to NA based on sample or received condition values (optional)
    if(!missing(toFilter)) {
      combinedDF$kclAmmoniumNConc[combinedDF$sampleCondition %in% toFilter] <- NA
      combinedDF$kclNitrateNitriteNConc[combinedDF$sampleCondition %in% toFilter] <- NA
      combinedDF$kclAmmoniumNConc[combinedDF$receivedCondition %in% toFilter] <- NA
      combinedDF$kclNitrateNitriteNConc[combinedDF$receivedCondition %in% toFilter] <- NA

      # compile list of how many values got set to NA with filtering
      if (any(combinedDF$sampleCondition %in% toFilter)) {
        num1 <-
          length(combinedDF$sampleID[combinedDF$sampleCondition %in% toFilter])
        warning1 <-
          paste(
            'warning:',
            num1,
            'records had concentration values set to NA due to anomolous conditions',
            sep = " "
          )
        print(warning1)
      } else {
        print()
      }
    } else {
      combinedDF
    }

  # set N data to NA based on ammonium or nitrate quality flags (optional)
    if(dropFlagged) {
      combinedDF$kclAmmoniumNConc[combinedDF$ammoniumNQF %in% c("1", "2")] <- NA
      combinedDF$kclNitrateNitriteNConc[combinedDF$nitrateNitriteNQF %in% c("1", "2")] <- NA
      } else {
      combinedDF
      }


  # add blank info & values
  combinedDF$blank1ID <-
    kclIntBlank$kclBlank1ID[match(combinedDF$kclReferenceID, kclIntBlank$kclReferenceID)]
  combinedDF$blank1NH4 <-
    kclExt$kclAmmoniumNConc[match(toupper(combinedDF$blank1ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank1NO3 <-
    kclExt$kclNitrateNitriteNConc[match(toupper(combinedDF$blank1ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank2ID <-
    kclIntBlank$kclBlank2ID[match(combinedDF$kclReferenceID, kclIntBlank$kclReferenceID)]
  combinedDF$blank2NH4 <-
    kclExt$kclAmmoniumNConc[match(toupper(combinedDF$blank2ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank2NO3 <-
    kclExt$kclNitrateNitriteNConc[match(toupper(combinedDF$blank2ID), toupper(kclExt$kclSampleID))]
  combinedDF$blank3ID <-
    kclIntBlank$kclBlank3ID[match(combinedDF$kclReferenceID, kclIntBlank$kclReferenceID)]
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
      combinedDF$kclNitrateNitriteNBlankCor
    ) # set to zero if negative

  # add soil moisture and dry mass fraction
  combinedDF$soilMoisture <-
    soilMoist$soilMoisture[match(combinedDF$sampleID, soilMoist$sampleID)]
  combinedDF$dryMassFraction <-
    soilMoist$dryMassFraction[match(combinedDF$sampleID, soilMoist$sampleID)]
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

  # determine whether to keep all variable or just a subset
  if (keepAll) {
    combinedDF
  } else {
    combinedDF <-
      subset(
        combinedDF,
        select = c(
          "plotID",
          "setDate",
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
    combinedDF <- combinedDF[!is.na(combinedDF$plotID), ]
    combinedDF[is.na(combinedDF)] <- NA
    combinedDF %>% dplyr::mutate_at(
      vars(
        dryMassFraction,
        soilDryMass,
        kclAmmoniumNBlankCor,
        kclNitrateNitriteNBlankCor,
        soilAmmoniumNugPerGram,
        soilNitrateNitriteNugPerGram,
        netNminugPerGramPerDay,
        netNitugPerGramPerDay
      ),
      funs(round(., 3))
    )
  }
}
