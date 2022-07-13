\name{PAO}
\alias{PAO}

\title{Percentage Above Optimal (PAO) path as a measure of organization in a foraging task.

}
\description{
The function calculates the PAO (Ólafsdóttir et al., 2021) based on the actual scanpath trough targets and the optimal path, the shortest Hamiltonian path between the targets locations.
}
\usage{
PAO(x, y)
}

\arguments{
  \item{x}{
A vector containing the x coordinates of the locations of the collected targets.
}
  \item{y}{
A vector containing the y coordinates of the locations of the collected targets.
}
}
\details{
The function estimates the Shortest Hamiltonian path among the locations of the targets collected. This shortest Hamiltonian path is the shortest (optimal) path which passes through the locations of all targets collected. The function also calculates the actual pathlength taking into account the order in which the targets were collected. The difference between the actual and the optimal paths is expressed in terms of a percentage: 0\%\ above optimal correspond to the optimal path, 100\%\ bove optimal path is the double of the optimal, and so forth.

To calculate the shortest Hamiltonian path, this functions uses two other functions: TSP from TSP (Hashler & Hornik, 2021) package, and order_tsp from PairViz package (Harley & Oldford, 2022).
}
\value{
A numeric value with the Percentag Above Optimal path.
}
\references{
Hahsler, M., and Hornik, K. (2021). \emph{TSP: Traveling Salesperson Problem (TSP)}. R
package version 1.1-11. \url{https://CRAN.R-project.org/package=TSP}.

Hurley, C. B., and Oldford, R. W. (2022). \emph{PairViz: Visualization using Graph Traversal.} R package version 1.3.5. \url{https://CRAN.R-project.org/package=PairViz}

Ólafsdóttir, I. M., Gestdóttir, S., & Kristjánsson, Á. (2021) The development of foraging organization. \emph{Attention, Perception, & Psychophysics, 83}, 2891-2904.


}
\author{
Marcos Bella-Fernández
}
