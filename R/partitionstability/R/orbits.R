library("sets")

.set_get <- function(s) {
  for (e in s) {
    return(e)
  }
}

#' computeOrbitPartition
#' 
#' Implementation of Algorithms 2 and 3 (combined as one function), which computes the
#' orbit partition of a permutation group from a given set of generators.
#'
#' @param S Set of generators
#' @param n Length of the set the group \eqn{\langle S\rangle} acts on
#' 
#' @examples 
#' library("sets")
#' p1 <- c(2L, 1L, 3L, 4L, 6L, 5L, 7L, 8L, 9L, 10L)
#' p2 <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 10L, 9L)
#' p3 <- c(1L, 9L, 3L, 4L, 5L, 6L, 7L, 8L, 2L, 10L)
#' S <- list(p1, p2, p3)
#' computeOrbitPartition(S, 10)
#' # [1] 1 1 3 4 5 5 7 8 1 1
#' 
#' @import sets
#' @export
computeOrbitPartition <- function(S, n) {
  O <- rep(-1L, n)
  colored <- 0L
  U = set()
  
  for (i in 1L:n) {
    if (O[i] != -1L) { next }
    
    N <- set(i)
    col <- i
    
    while (set_cardinality(N) > 0) {
      j <- .set_get(N)
      N <- N - j
      
      if (O[j] == -1L) {
        O[j] <- col
        colored <- colored + 1L
      }
      
      for (pi in S) {
        if (pair(j, pi) %e% U) { next }
        
        U <- U | set(pair(j, pi))
        k <- pi[j]
        
        while (k != j) {
          if (O[k] == -1L) {
            O[k] <- col
            N <- N | set(k)
            colored <- colored + 1L
            U <- U | set(pair(k, pi))
          }
          
          k <- pi[k]
        }
        
        if (colored == n) { return(O) }
      }
    }
  }
  
  return(O)
}
