
pick_start_k <- function(y, m=5000, a1=1,b1=1,a2=1,b2=1, tries=100){
  n <- length(y)
  C <- cumsum(y); S <- sum(y)

  ks <- sample(1:(n-1), size=tries,replace='TRUE')

  score <- numeric(tries)

  for(t in seq_along(ks)){
    k <- ks[t]
    S1 <- C[k];        N1 <- k*m
    S2 <- S - C[k];    N2 <- (n-k)*m

    p1 <- (a1 + S1) / (a1 + b1 + N1)
    p2 <- (a2 + S2) / (a2 + b2 + N2)

    score[t] <- S1*log(p1) + (N1-S1)*log(1-p1) + S2*log(p2) + (N2-S2)*log(1-p2)
  }

  k0 <- ks[which.max(score)]
  list(k0=k0, ks=ks, score=score)
}

gibbs_isochores <- function(data,a1=1,a2=1,b1=1,b2=1, m = 5000, nburn = 10000, ndraw = 1000) {

  #θα χρησιμοποιησουμε prior την beta(1,1) 
  #οριζουμε τα αθροισματα μας 
  n <- length(data)
  C <- cumsum(data)
  S <- C[n]  # συνολικες επιτυχιες
#αρχικοποιηση παραμετρων
p1 <- rbeta(1,1,1)                 
p2 <- rbeta(1,1,1)

  k <-sample(1:(n-1),size=1)
  # πινακας που κραταει τις αρχικες συνθηκες
  cat("INIT: p1=", p1, " p2=", p2, " k=", k, "\n")

  draws <- matrix( ncol = 3,nrow = ndraw)
  colnames(draws) <- c("p1", "p2", "k")
#αρχιζει το gibbs επισημως 
  it <- -nburn
  while (it < ndraw) {
    it <- it + 1

    logw <- numeric(n - 1)
    for (j in 1:(n - 1)) {
      Sj1 <- C[j]
      Nj1 <- j * m
      Sj2 <- S - C[j]
      Nj2 <- (n - j) * m

      logw[j] <- Sj1*log(p1)+(Nj1-Sj1)*log(1-p1)+Sj2*log(p2)+(Nj2-Sj2)*log(1-p2)
    }


    
    # 'ξελογκαρουμε' τις πιθανοτητες μας και κανονικοποιουμε
    # αφαιρουμε με max για να μην εχουμε θεματα NaN 
    logw <- logw - max(logw)
    w <- exp(logw)
sw <- sum(w)
  w <- w / sw


        # επιλεγουμε το k με τις καινουργιες κατανομες μας για να κατευθυνουμε το επομενο iteration
     k <- sample(1:(n - 1), size = 1, prob = w)
              S1 <- C[k]
    S2 <- S - C[k]

    N1 <- k * m
    N2 <- (n - k) * m
        # παιρνουμε την posterior πριν την αλλαγη
    p1 <- rbeta(1, a1 + S1, b1 + (N1 - S1))

    # παιρνουμε την posterior μετα την αλλαγη
    p2 <- rbeta(1, a2 + S2, b2 + (N2 - S2))
    eps <- 1e-12
p1 <- min(max(p1, eps), 1 - eps)#κατω φραγμα για να ειμαστε ασφαλης 
p2 <- min(max(p2, eps), 1 - eps)
    # aαφου τελειωσει το burn in period αρχιζουμε να κραταμε τα στοιχεια μας 
    if (it > 0) {
      draws[it, ] <- c(p1, p2, k)
    }
  }

  return(draws)
}

chr1.3=c(1776, 1929, 1885, 2058, 1907, 2143, 2206, 2112, 2131, 1948, 
2291, 2112, 2122, 2985, 2690, 2634, 2339, 2924, 2805, 2755, 2901, 
2527, 2159, 2337, 2237, 2224, 2251, 2384, 2340, 2234, 2642, 2422, 
2381, 2527, 2428, 2167, 2275, 2589, 2364, 2407, 2736, 2452, 2600, 
2554, 2640, 2396, 2488, 2660, 2512, 2746, 2559, 2399, 2770, 2870, 
2300, 2391, 1936, 2308, 2128, 1940, 1925, 2108, 2343, 2223, 2174, 
2277, 2480, 2749, 2794, 2846, 2792, 2668, 2712, 2793, 2209, 2381, 
2218, 2368, 2503, 2619, 2378, 2262, 2145, 2008, 2186, 1994, 2075, 
1933, 2214, 2408, 2049, 2110, 2391, 2218, 2337, 2808, 2244, 2213, 
2569, 2455)
draws <- gibbs_isochores(chr1.3, m = 5000, nburn = 1000, ndraw = 1000)
draws[nrow(draws), ]
