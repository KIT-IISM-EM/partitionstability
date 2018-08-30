#' @details 
#' Algorithms 1 (\link{testStability}) and 4 (\link{geq}) are nearly one-to-one 
#' implementations of the pseudocode. One exception is that we removed the need of a
#' hash map in \link{geq} by simply using a vector of the same length as the partition
#' to build our map of cluster ids. However, this implies that cluster ids must be the 
#' same as the node ids.
#' 
#' Algorithms 2 and 3 are merged into a single function (\link{computeOrbitPartition}) 
#' as R does not easily support call be reference.
#' 
#' For details of the algorithms and the background we refer to 
#' Ball, Fabian and Geyer-Schulz, Andreas (2018), 
#' ``Symmetry-based Graph Clustering Partition Stability'', 
#' Archives of Datascience Series A <vol tba>(<nr tba>), DOI: <doi tba>
#'
#' @importFrom Rdpack reprompt
"_PACKAGE"

#' testStability
#' 
#' Implementation of Algorithm 1, which checks if a partition is stabl udner the given
#' symmetry implied by the generators of a permutation group \eqn{S}.
#' 
#' @param P A partition
#' @param S Set of generators for a group \eqn{\langle S\rangle} that acts on \eqn{P}
#'
#' @examples 
#' p1 <- c(2L, 1L, 3L, 4L, 6L, 5L, 7L, 8L, 9L, 10L)
#' p2 <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 10L, 9L)
#' p3 <- c(1L, 9L, 3L, 4L, 5L, 6L, 7L, 8L, 2L, 10L)
#' S <- list(p1, p2, p3)
#' P <- rep(1L, 10)
#' isTRUE(testStability(P, S))
#' Q <- c(2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 3L, 3L)
#' !isTRUE(testStability(Q, S))
#' 
#' @export
testStability <- function(P, S) {
  if (length(S) == 0) return(TRUE)
  
  O <- computeOrbitPartition(S, length(P))
  
  if (isTRUE(geq(P, O))) return(TRUE)
  
  for (pi in S) {
    P_pi <- vapply(P, function(e) {P[pi[e]]}, 0L)
    if (!isTRUE(geq(P, P_pi))) return(FALSE)
  }
  
  return(TRUE)
}

#' geq
#' 
#' Implementation of Algorithm 4, which tests \eqn{P \geq Q}
#' 
#' 
#' @param P A partition of length n
#' @param Q Another partition of length n#'
#'
#' @examples
#' P = c(1, 1, 1, 2, 2, 2)
#' P_prime = c(2, 2, 2, 1, 1, 1)
#' Q = c(4, 4, 1, 3, 3, 2)
#' R = c(4, 4, 4, 4, 2, 2)
#' isTRUE(geq(P, Q))
#' isTRUE(geq(P, P_prime))
#' isTRUE(geq(P, P_prime) == geq(P_prime, P))
#' !isTRUE(geq(P, R))
#' !isTRUE(geq(R, P))
#'
#' @export
geq <- function(P, Q) {
  maps <- rep(-1, length(P))
  
  for (i in 1:length(P)) {
    if (maps[Q[i]] != -1) {
      oid <- maps[Q[i]]
    } else {
      oid <- Q[i]
      maps[Q[i]] <- oid
    }
    
    if (P[i] != oid) return(FALSE)
  }
  
  return(TRUE)
}
