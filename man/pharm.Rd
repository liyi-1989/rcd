\name{pharm}
\alias{pharm}

\title{
Pharmaceutical company stock price .
}

\description{
The stock prices for three pharmaceutical companies and a market index. They are ARIAD Pharmaceuticals (ARIA), Vertex Pharmaceuticals (VRTX), and Shire Plc (SHPG), and the Nasdaq Composite (IXIC). 
}

\usage{data(pharm)}

\details{
	This data set contains 4054 samples with 4 predictors. The date range from Jan. 3, 2000 to Feb. 12, 2006.
}

\format{
  The format is a dataframe.
  1. Date - The date range from Jan. 3, 2000 to Feb. 12, 2006.
  2. ARIA - The stock price for ARIAD Pharmaceuticals.
  3. VRTX - The stock price for Vertex Pharmaceuticals.
  4. SHPG - The stock price for Shire Plc.
  5. IXIC - The index for Nasdaq Composite.
}

\examples{
data(pharm)
pairs(pharm[,2:5])
}

\source{http://finance.yahoo.com/}