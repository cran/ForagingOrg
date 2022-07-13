\name{meanITD}
\alias{meanITD}

\title{
Mean Intertarget Distance (ITD) from x and y coordinates for target positions.
}
\description{
Every pair of consecutively collected targets are at a certain Euclidean distance, which can be calculated from their respective x and y coordinates. The mean ITD is simply the mean of all these distances (Ólafsdóttir et al., 2021; Woods et al., 2013).
}
\usage{
meanITD(x, y)
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
A numeric value with the Mean ITD.
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

meanITD(X,Y)
}
