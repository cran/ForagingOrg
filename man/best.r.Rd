\name{best.r}
\alias{best.r}

\title{Best-r index for foraging organization

}
\description{
Calculates best-r indicator from a set of x and y coordinates. Best-r is defined as the best Pearson correlation (in absolute value) between the order in which targets are collected and their respective x and y coordinates (Ólafsdóttir et al., 2021; Woods et al., 2013).
}
\usage{
best.r(x,y)
}

\arguments{
  \item{x}{
  A vector containing the x coordinates from the target locations.
}
  \item{y}{
  A vector containing the y coordinates from the target locations.
  }
}

\value{

Returns the best-r as a numeric value.

}
\references{
Ólafsdóttir, I. M., Gestdóttir, S., & Kristjánsson, Á. (2021) The development of foraging organization. \emph{Attention, Perception, & Psychophysics, 83}, 2891-2904.

Woods, A. J., Göksun, T., Chatterjee, A., & Zelonis, S. (2013) The development of organized visual search. \emph{Acta Psychologica, 143}, 191-199.
}
\author{
Marcos Bella-Fernández
}

\examples{

#Simulating target coordinates

n<-20 #Number of targets collected in the simulation.

X<-runif(n,0,150) #Simulating x coordinates.
Y<-runif(n,0,100) #Simulating y coordinates.

best.r(X,Y)
}
