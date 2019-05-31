#' Results of Ukrainian External Independent Testing 2016
#'
#' A dataset containing results of tests in Mathematics and
#' Ukrainian language and literature at External Independent Testing
#' 2016 in Ukraine. The data on persons which passed both Math and Ukrainian
#' tests are included only.
#'
#' @format A data frame with 94680 and 3 variables. Each row corresponds
#' to one person's test results. The variables are:
#' \describe{
#'  \item{obl}{code of region of Ukraine (oblast) where the examenee
#'    finished the high school (the codes are the same as in ukrelect2014
#'    data file)}
#'   \item{ukr}{result of test on Ukrainian language and literature}
#'  \item{math}{result of test on mathematics}
#' }
#' @source Official site of Ukrainian Center for
#' Educational Quality Assessment
#' \url{https://zno.testportal.com.ua/stat/2016}
"EIT2016"
#________________________________________________________
#' Results of Ukrainian Parliament elections 2014
#'
#' A dataset containing data on proportions of different
#' electoral choices in different regions (oblast) in Ukraine
#' at Parliament (Verhovna Rada) elections in 2014.
#'
#' @format A data frame with 27 rows corresponding to different
#' regions of Ukraine and 5 variables:
#' \describe{
#'   \item{name}{name of the region}
#'   \item{code}{code of the region}
#'   \item{ProEU}{percent of registered voters that voted for
#'   parties which formed pro-EU coalition  (BPP, Batkivschyna,
#'    Narodny Front, Radicals and Samopomich)}
#'    \item{ContraEU}{percent of registered voters that voted for
#'    the Opposition block}
#'    \item{Neutral}{percent of registered voters that voted for
#'    small parties, against all or didn't take part in the voting}
#' }
#'
#' (Data on four regions are \code{NA} due to special political situation
#' at these regions)
#' @source Official site of Central Election Commission (Ukraine)
#' \url{http://www.cvk.gov.ua/vnd_2014/}
"ukrelect2014"
#___________________________________________________
#' Expression levels of genes from I. Hedenfalk at al. (2001)
#' study.
#'
#' A dataset containing lewel of expression of more than 3226 genes
#' in breast cancer tissues of 22 patients. The cancer can be caused
#' by BRCA1 or BRCA2 mutation or can be of Sporadic type.
#'
#' @format A data frame with 3226 rows corresponding to genes and 22 variables
#' corresponding to tissue specimens. Names of variables correspond to
#' the type of tissue:
#' \describe{
#'  \item{BRCA1.1}{tissue with the BRCA1 mutation, 1 is the number of
#'  specimen}
#'  \item{BRCA1.2}{tissue with the BRCA1 mutation, 2 is the number of
#'  specimen}
#'  \item{BRCA1.3}{tissue with the BRCA1 mutation, 3 is the number of
#'  specimen}
#'  \item{BRCA1.4}{tissue with the BRCA1 mutation, 4 is the number of
#'  specimen}
#'  \item{BRCA1.5}{tissue with the BRCA1 mutation, 5 is the number of
#'  specimen}
#'  \item{BRCA1.6}{tissue with the BRCA1 mutation, 6 is the number of
#'  specimen}
#'  \item{BRCA1.7}{tissue with the BRCA1 mutation, 7 is the number of
#'  specimen}
#'  \item{BRCA2.1}{tissue with the BRCA2 mutation, 1 is the number of
#'  specimen}
#'  \item{BRCA2.2}{tissues with the BRCA2 mutation, 2 is the number of
#'  specimen}
#'  \item{BRCA2.3}{tissue with the BRCA2 mutation, 3 is the number of
#'  specimen}
#'  \item{BRCA2.4}{tissues with the BRCA2 mutation, 4 is the number of
#'  specimen}
#'  \item{BRCA2.5}{tissue with the BRCA2 mutation, 5 is the number of
#'  specimen}
#'  \item{BRCA2.6}{tissues with the BRCA2 mutation, 6 is the number of
#'  specimen}
#'  \item{BRCA2.7}{tissue with the BRCA2 mutation, 7 is the number of
#'  specimen}
#'  \item{BRCA2.8}{tissue with the BRCA2 mutation, 8 is the number of
#'  specimen}
#'  \item{Sporadic.1}{tissue with sporadic tumors, 1 is the number of
#'  specimen}
#'  \item{Sporadic.2}{tissue with sporadic tumors, 2 is the number of
#'  specimen}
#'  \item{Sporadic.3}{tissue with sporadic tumors, 3 is the number of
#'  specimen}
#'  \item{Sporadic.4}{tissue with sporadic tumors, 4 is the number of
#'  specimen}
#'  \item{Sporadic.5}{tissue with sporadic tumors, 5 is the number of
#'  specimen}
#'  \item{Sporadic.6}{tissue with sporadic tumors, 6 is the number of
#'  specimen}
#'  \item{Sporadic_Meth.BRCA1}{tissue with sporadic tumors in which a
#'  methilized BRCA1 mutation is present}
#' }
#'
#' @source  Hedenfalk I. at al. Gene-expression profiles in hereditary breast cancer.
#' N Engl J Med. 2001 Feb 22;344(8):539-48.
"hedenf"
