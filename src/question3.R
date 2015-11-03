learningDualPerceptron <- function (positive, negative) {
  positivePoints <- matrix(positive, ncol=2, byrow=TRUE)
  negativePoints <- matrix(negative, ncol=2, byrow=TRUE)
  allPoints      <- rbind(cbind(positivePoints, 1), cbind(negativePoints, -1))
  lengths        <- nrow(allPoints)
  colnames(allPoints) <- c('x1', 'x2', 'y')
  # Plots.
  plot(allPoints[], type="n")
  points(positivePoints, col="blue", pch=19)
  points(negativePoints, col="red", pch=19)
  # Learning perceptron.
  alpha <- rep(0, times=lengths)
  b     <- 0
  R     <- max(sqrt(abs(allPoints[,1])^2 + abs(allPoints[,2])^2))^2
  for (time in 1:15) {
    for (i in 1:lengths) {
      xi    <- allPoints[i, c('x1', 'x2')]
      yi    <- as.numeric(allPoints[i, 'y'])
      # Get the Sigma.
      sigma <- 0
      for (j in 1:lengths) {
        aj <- alpha[j]
        xj <- allPoints[j, c('x1', 'x2')]
        yj <- as.numeric(allPoints[j, 'y'])
        sigma <- as.numeric(as.character(sigma + aj * yj * (xi %*% xj)^2))
      }
      # If total <= 0.
      if ((yi * (sigma + b)) <= 0) {
        alpha[i] <- alpha[i] + 1
        b        <- as.numeric(as.character(b + yi * R^2))
      }
    }
  }
  # Remove the points that are alpha = 0.
  for (i in 1:lengths) {
    if (alpha[i] == 0) {
      allPoints[i, 'y'] <- 0
    }
  }
  points(allPoints[allPoints[,'y']==0,], col="white", pch=19)
  # Return results.
  return(list(allPoints=allPoints, alpha=alpha, b=b))
}

randomTestsetPlot <- function (input) {
  allPoints <- input$allPoints
  alpha     <- input$alpha
  b         <- input$b
  # Total counts.
  randomCounts <- 10000
  # Container size.
  # box <- c(-1.5, 1.5) %o% c(-1.5, 1.5)
  box <- c(-1.5, 1.5)
  randomPoints <- matrix(runif(2 * randomCounts, min(box), max(box)), ncol=2, byrow=TRUE)
  colnames(randomPoints) <- c('x1', 'x2')
  randomPoints <- t(apply(randomPoints, 1, function (point) {
    sigma <- 0
    for (i in 1:nrow(allPoints)) {
      ai    <- alpha[i]
      yi    <- as.numeric(as.character(allPoints[i,'y']))
      xi    <- as.matrix(allPoints[i, c('x1', 'x2')])
      sigma <- sigma + ai * yi * (t(xi) %*% point)^2
    }
    point['y'] <- ifelse(sign(sigma + b) > 0, 1, -1)
    return(point)
  }))
  plot(randomPoints[], type="n")
  points(randomPoints[randomPoints[,'y']==1,], col="#6699FF", pch=20)
  points(randomPoints[randomPoints[,'y']==-1,], col="#EEEEEE", pch=20)
}

learningPrimalPerceptronFrom2dTo4d <- function (positive, negative) {
  posPoints <- matrix(positive, ncol=2, byrow=TRUE)
  negPoints <- matrix(negative, ncol=2, byrow=TRUE)
  posPoints <- matrix(positive, ncol=2, byrow=TRUE)
  negPoints <- matrix(negative, ncol=2, byrow=TRUE)
  p         <- rbind(cbind(posPoints, 1), cbind(negPoints, -1))
  phi <- matrix(c(p[,1] * p[,2], p[,1]^2, p[,1] * p[,2] * -1, p[,2]^2, p[,3]), ncol=5)
  colnames(phi) <- c('z1', 'z2', 'z3', 'z4', 'y')
  lengths <- nrow(phi)
  # Plots.
  plot(phi[], type="n")
  # Learning perceptron.
  error <- TRUE
  b <- 0
  k <- 0
  w <- matrix(rep(0, times=4), ncol=4)
  R <- max(sqrt(phi[,1]^2 + phi[,2]^2 + phi[,3]^2 + phi[,4]^2))
  #
  while (error) {
    error <- FALSE
    for (i in 1:lengths) {
      xi <- as.matrix(phi[i, 1:4])
      yi <- as.numeric(phi[i, 'y'])
      # If total <= 0.
      if ((yi * (w %*% xi + b)) <= 0) {
        error <- TRUE
        w <- as.numeric(as.character(w + 1 * yi * t(xi)))
        b <- as.numeric(as.character(b + 1 * yi * R^2))
      }
    }
    k <- k + ifelse(error, 1, 0)
  }
  # Return results.
  return(list(allPoints=phi, k=k, b=b, w=w))
}

randomTestsetPlot2 <- function (input) {
  allPoints <- input$allPoints
  w         <- matrix(input$w, ncol=4)
  b         <- input$b
  # Total counts.
  randomCounts <- 10000
  # Container size.
  box <- c(-1.5, 1.5)
  randomPoints <- matrix(runif(2 * randomCounts, min(box), max(box)), ncol=2, byrow=TRUE)
  randomPoints <- t(apply(randomPoints, 1, function (p) {
    p_d4 <- matrix(c(p[1] * p[2], p[1]^2, p[1] * p[2] * -1, p[2]^2), ncol=4)
    p['y'] <- ifelse(sign((w %*% t(p_d4)) + b) > 0, 1, -1)
    return(p)
  }))
  # print(randomPoints)
  colnames(randomPoints) <- c('x1', 'x2', 'y')
  plot(randomPoints[, 1:2], type="n")
  # print(t(w))
  points(randomPoints[randomPoints[,'y']==1,], col="#6699FF", pch=20)
  points(randomPoints[randomPoints[,'y']==-1,], col="#EEEEEE", pch=20)
  points(matrix(positiveA, ncol=2, byrow=TRUE), col="#003399", pch=15)
  points(matrix(negativeA, ncol=2, byrow=TRUE), col="#FF0000", pch=15)
}

# Problem A
positiveA <- c(c(0, 0), c(0.5, 0), c(0, 0.5), c(-0.5, 0), c(0, -0.5))
negativeA <- c(c(0.5, 0.5), c(0.5, -0.5), c(-0.5, 0.5), c(-0.5, -0.5), c(1, 0), c(0, 1), c(-1, 0), c(0, -1))
PLA_A     <- learningDualPerceptron(positiveA, negativeA)

# Problem B
# set.seed(1)
randomTestsetPlot(PLA_A)

# Problem C
positiveB <- c(c(0.5, 0), c(0, 0.5), c(-0.5, 0), c(0, -0.5))
negativeB <- c(c(0.5, 0.5), c(0.5, -0.5), c(-0.5, 0.5), c(-0.5, -0.5))
PLA_B     <- learningDualPerceptron(positiveB, negativeB)
# set.seed(1)
randomTestsetPlot(PLA_B)

# Problem D
phi <- learningPrimalPerceptronFrom2dTo4d(positiveA, negativeA)
# set.seed(1)
randomTestsetPlot2(phi)

