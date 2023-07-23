#' Text based representation of a vlmc
#'
#' @param ct a fitted vlmc.
#' @param prob this parameter controls the display of node level information in
#'   the tree. The default `prob=TRUE` represents the conditional distribution
#'   of the states given the (partial) context associated to the node. Setting
#'   `prob=FALSE` replaces the conditional distribution by the frequency of the
#'   states that follow the context as in [draw.ctx_tree()]. Setting `prob=NULL`
#'   removes all additional information.
#' @examples
#' dts <- sample(c("A", "B", "C"), 500, replace = TRUE)
#' model <- vlmc(dts, alpha = 0.05)
#' drawV(model)
#' drawV(model, prob = FALSE)
#' drawV(model, prob = NULL)
#' @import igraph
#' @export
drawV <- function(ct, prob = TRUE) { # type : NULL or prob
  if (is.null(prob)) {
    drawNull(ct)
  } else {
    if (depth(ct) == 0) {
      adjM <- matrix(c(0, 0, 1, 0), 2)
      if (prob) {
        colnames(adjM) <- c("", paste("p", paste(round(ct$f_by / sum(ct$f_by), 3), collapse = " "), sep = ""))
      } else {
        colnames(adjM) <- c("", paste("p", paste(ct$f_by, collapse = " "), sep = ""))
      }
      g <- graph_from_adjacency_matrix(adjM)
    } else {
      Thing <- list()
      Thing$adjM <- matrix(rep(0, (2 * ct$nb_ctx)^2), ncol = 2 * ct$nb_ctx)
      Thing$nb_node <- 1
      Thing$vals <- ct$vals
      adjM <- Thing$adjM
      colnames(adjM) <- rep(NA, ct$nb_ctx * 2)
      from <- 1
      to <- 1
      n <- length(Thing$vals)
      m <- 0
      for (i in 1:n) { # number of children
        if (length(ct$children[[i]]) > 0) {
          m <- m + 1
        }
      }
      adjM[from, (Thing$nb_node + 1)] <- 1
      if (prob) {
        colnames(adjM)[(Thing$nb_node + 1)] <- paste("p", paste(round(ct$f_by / sum(ct$f_by), 3), collapse = " "), sep = "")
      } else {
        colnames(adjM)[(Thing$nb_node + 1)] <- paste("p", paste(ct$f_by, collapse = " "), sep = "")
      }
      if (m == n) { # number of children
        adjM <- rbind(rbind(cbind(cbind(adjM, 0), 0), 0), 0)
      }
      i <- 1
      count <- 0
      Thing$adjM <- adjM
      while (count < m) {
        if (length(ct$children[[i]]) > 0) {
          count <- count + 1
          Thing <- draw_rec(ct$children[[i]], from, Thing$nb_node, Thing, paste(Thing$vals[i]), prob)
          adjM <- Thing$adjM
        }
        i <- i + 1
      }
      g <- graph_from_adjacency_matrix(adjM, add.colnames = NULL)
    }
    V(g)$color <- c("orange", ifelse(substring(V(g)$name[2:(length(V(g)))], 1, 1) == "p", "lightblue", "orange"))
    V(g)$name <- substring(V(g)$name, 2)
    # plot(g, layout=layout_with_lgl(g, root = 1))
    plot(g, layout = layout_as_tree(g, root = 1))
    # plot(g, layout= layout_as_star(g))
    # plot(g, layout= layout_nicely(g))
  }
}

draw_rec <- function(obj, from, to, Thing, nb, prob) { # final si == 1
  n <- length(Thing$vals)
  adjM <- Thing$adjM
  Thing$nb_node <- Thing$nb_node + 1

  adjM[from, (Thing$nb_node + 1)] <- 1
  adjM[(Thing$nb_node + 1), (Thing$nb_node + 2)] <- 1
  if (prob) {
    colnames(adjM)[(Thing$nb_node + 2):(Thing$nb_node + 1)] <- c(paste("p", paste(round(obj$f_by / sum(obj$f_by), 3), collapse = " "), sep = ""), paste("v", nb, sep = ""))
  } else {
    colnames(adjM)[(Thing$nb_node + 2):(Thing$nb_node + 1)] <- c(paste("p", paste(obj$f_by, collapse = " "), sep = ""), paste("v", nb, sep = ""))
  }

  Thing$nb_node <- Thing$nb_node + 1
  from <- Thing$nb_node
  m <- 0
  for (i in 1:n) {
    if (length(obj$children[[i]]) > 0) {
      m <- m + 1
    }
  }
  if (m == n) { # expand matrix for number of context
    adjM <- rbind(rbind(cbind(cbind(adjM, 0), 0), 0), 0)
  }
  i <- 1
  count <- 0
  Thing$adjM <- adjM
  while (count < m) {
    if (length(obj$children[[i]]) > 0) {
      count <- count + 1
      Thing <- draw_rec(obj$children[[i]], from, Thing$nb_node, Thing, paste(Thing$vals[i]), prob)
      adjM <- Thing$adjM
    }
    i <- i + 1
  }
  Thing$adjM <- adjM
  return(Thing)
}

drawNull <- function(ct) {
  if (depth(ct) == 0) {
    adjM <- matrix(0)
    colnames(adjM) <- ""
    g <- graph_from_adjacency_matrix(adjM)
    plot(g)
  } else {
    Thing <- list()
    Thing$adjM <- matrix(rep(0, (ct$nb_ctx)^2), ncol = ct$nb_ctx)
    Thing$nb_node <- 1
    Thing$vals <- ct$vals
    adjM <- Thing$adjM
    colnames(adjM) <- rep(NA, ct$nb_ctx)
    from <- 1
    to <- 1
    n <- length(Thing$vals)
    m <- 0
    for (i in 1:n) { # number of children
      if (length(ct$children[[i]]) > 0) {
        m <- m + 1
      }
    }
    if (m == n) { # number of children
      adjM <- rbind(cbind(adjM, 0), 0)
    }
    i <- 1
    count <- 0
    Thing$adjM <- adjM
    while (count < m) {
      if (length(ct$children[[i]]) > 0) {
        count <- count + 1
        Thing <- draw_recNull(ct$children[[i]], from, Thing$nb_node, Thing, paste(Thing$vals[i]))
        adjM <- Thing$adjM
      }
      i <- i + 1
    }
    g <- graph_from_adjacency_matrix(adjM, add.colnames = NULL)
    # plot(g, layout=layout_with_lgl(g, root = 1))
    plot(g, layout = layout_as_tree(g, root = 1))
    # plot(g, layout= layout_as_star(g))
    # plot(g, layout= layout_nicely(g))
  }
}

draw_recNull <- function(obj, from, to, Thing, nb) { # final si == 1
  n <- length(Thing$vals)
  adjM <- Thing$adjM
  Thing$nb_node <- Thing$nb_node + 1

  adjM[from, (Thing$nb_node)] <- 1
  colnames(adjM)[(Thing$nb_node)] <- nb

  from <- Thing$nb_node
  m <- 0
  for (i in 1:n) {
    if (length(obj$children[[i]]) > 0) {
      m <- m + 1
    }
  }
  if (m == n) { # expand matrix for number of context
    adjM <- rbind(cbind(adjM, 0), 0)
  }
  i <- 1
  count <- 0
  Thing$adjM <- adjM
  while (count < m) {
    if (length(obj$children[[i]]) > 0) {
      count <- count + 1
      Thing <- draw_recNull(obj$children[[i]], from, Thing$nb_node, Thing, paste(Thing$vals[i]))
      adjM <- Thing$adjM
    }
    i <- i + 1
  }
  Thing$adjM <- adjM
  return(Thing)
}
