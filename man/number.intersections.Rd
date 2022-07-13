\name{number.intersections}
\alias{number.intersections}

\title{
Counter for the number of intersections of non-consecutive trajectories between targets among a scanpath in a foraging task.
}
\description{
This function counts the number and rate of intersections between non-consecutive targets in a scanpath.
}
\usage{
number.intersections(x, y)
}

\arguments{
  \item{x}{
A vector with the x coordinates of the targets collected in a foraging task.
}
  \item{y}{
A vector with the y coordinates of the targets collected in a foraging task.
}
}
\details{
The number of intersections between non-consecutives scanpaths in a foraging task is taken as a foraging organization measure (Ólafsdóttir et al., 2021; Woods et al., 2013). Here we developed an algorithm based on function analysis and Bolzano's Theorem (e. g. Apostol, 1974) to calculate the number of intersections. }
\value{
A numeric value with the number of intersections in a scanpath.
}
\references{
Apostol, T. M. (1974) \emph{Mathematical Analysis (2nd Edition)}. Addison Wesley Longman.

Ólafsdóttir, I. M., Gestdóttir, S., & Kristjánsson, Á. (2021) The development of foraging organization. \emph{Attention, Perception, & Psychophysics, 83}, 2891-2904.

Woods, A. J., Göksun, T., Chatterjee, A., & Zelonis, S. (2013) The development of organized visual search. \emph{Acta Psychologica, 143}, 191-199.
}
\author{
Marcos Bella-Fernández
}
