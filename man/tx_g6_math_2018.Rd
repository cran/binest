\name{tx_g6_math_2018}
\alias{tx_g6_math_2018}
\docType{data}

\title{Texas STAAR Grade-6 Mathematics, 2017-18: District-Level Bin Counts}

\description{
  District-level counts of students in each of four proficiency
  categories on the Texas State of Texas Assessments of Academic
  Readiness (STAAR) Grade-6 mathematics test, 2017-18 administration.
  For each district the dataset also reports the average scale score
  across all tested students, which can be used as ground truth for
  evaluating estimators that recover district means from binned
  counts.
}

\usage{data(tx_g6_math_2018)}

\format{
  A data frame with 1151 rows and 8 columns:
  \describe{
    \item{district_id}{Sequential integer identifier (1 to 1151).}
    \item{district_name}{District name (\code{character}).}
    \item{n_tested}{Total students tested in the district.}
    \item{unsatisfactory}{Students with scale score below 1536
                         (proficiency category "Did Not Meet Grade Level").}
    \item{approaches}{Students with scale score in [1536, 1653)
                       ("Approaches Grade Level").}
    \item{meets}{Students with scale score in [1653, 1772)
                  ("Meets Grade Level").}
    \item{masters}{Students with scale score >= 1772
                    ("Masters Grade Level").}
    \item{reported_mean}{District average scale score, computed by the
                          Texas Education Agency from individual
                          student scores.}
  }
}

\details{
  The three published cut scores defining the bin boundaries are
  1536, 1653, and 1772. The administrative floor of the STAAR scale
  is 1062 and the ceiling is 2143. Of the 1151 districts, 1014 have
  nonzero counts in all four bins, 120 have nonzero counts in three
  bins, and 17 have nonzero counts in two bins.
}

\source{
  Texas Education Agency, Academic Performance Reports (TAPR), 2017-18.
  Compiled by D.\ J.\ Hunter and P.\ T.\ von Hippel.
}

\examples{
data(tx_g6_math_2018)
str(tx_g6_math_2018)

## Recover district means using bin means with known cutpoints.
ngk <- with(tx_g6_math_2018,
            cbind(unsatisfactory, approaches, meets, masters))
fit <- bin_means(ngk, cutpoints = c(1536, 1653, 1772))

## Correlation with reported truth on the test-score scale.
## (Districts with fewer than three populated bins are NA-coded by
## bin_means; use complete.obs for the comparison.)
cor(fit$est_raw$group_mean_mle, tx_g6_math_2018$reported_mean,
    use = "complete.obs")
}

\keyword{datasets}
