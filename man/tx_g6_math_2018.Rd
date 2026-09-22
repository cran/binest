\name{tx_g6_math_2018}
\alias{tx_g6_math_2018}
\docType{data}

\title{tx_g6_math_2018: Texas STAAR Grade 6 Mathematics Scores: District-Level Bin Counts}

\description{
  \code{tx_g6_math_2018} is a bundled dataset used in the vignette. It
  can be used to test the properties of HETOP estimators in empirical
  data that may violate some model assumptions.

  \code{tx_g6_math_2018} provides district-level counts of students in
  each of four score bins on the Texas State of Texas Assessments of
  Academic Readiness (STAAR) Grade 6 mathematics test, 2017-18
  administration. For each district the dataset also reports the
  average score reported by the Texas Education Agency, which can be
  used as ground truth for evaluating estimators that recover district
  means from binned counts.
}

\usage{data(tx_g6_math_2018)}

\format{
  A data frame with 1151 rows and 8 columns:
  \describe{
    \item{district_id}{Sequential integer identifier (1 to 1151).}
    \item{district_name}{District name (\code{character}).}
    \item{n_tested}{Total students tested in the district.}
    \item{unsatisfactory}{Students scoring below 1536
                         (proficiency category "Did Not Meet Grade Level").}
    \item{approaches}{Students scoring in [1536, 1653)
                       ("Approaches Grade Level").}
    \item{meets}{Students scoring in [1653, 1772)
                  ("Meets Grade Level").}
    \item{masters}{Students scoring >= 1772
                    ("Masters Grade Level").}
    \item{reported_mean}{District average score, computed by the
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

## Recover district means using fast_hetop with known cutpoints.
ngk <- with(tx_g6_math_2018,
            cbind(unsatisfactory, approaches, meets, masters))
## These are administrative counts covering every tested student, so
## each district's students are its whole population: scope = "population".
fit <- fast_hetop(ngk, cutpoints_known = TRUE,
                  cutpoints = c(1536, 1653, 1772),
                  scope = "population")

## Correlation with reported truth on the test-score scale.
## (The 17 districts with fewer than three populated bins cannot
## identify both a mean and an SD; by default fast_hetop() estimates
## their means using a borrowed SD.  Pass
## estimate_unidentified_districts = FALSE to get NA instead, in which
## case use = "complete.obs" is needed here.)
cor(fit$est_raw$mean, tx_g6_math_2018$reported_mean)
}

\keyword{datasets}
