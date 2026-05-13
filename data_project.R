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
DNA_project <-function(data){
n=100
n_1=n-1
s<-rep(0,n_1)
  
  
for(i in 1:n_1) s[i]=sum(data[1:i])
sdata=sum(data)
logpr=rep(0,n_1) #διανυσμα για την posterior το k
    
for(i in 1:n_1){
    logpr[i]=lbeta(s[i]+1,i*5000-s[i]+1)+lbeta(sdata-s[i]+1,(n-i)*5000-(sdata-s[i])+1)
  #εδω εχουμε την posterior του k για καθε k
}
  
logprmax=max(logpr)
pr=exp(logpr-logprmax)
pr=pr/sum(pr) #κανονικοποιουμε

barplot(pr,names=1:n_1,main="Posterior Probability of Changepoint",xlab="place",ylab="Probability")

logev1 <- lbeta(sdata + 1, 5000*n - sdata + 1)
logev2 <- -log(n-1) + logprmax + log(sum(exp(logpr - logprmax)))

Postprob1 <- plogis(logev1-logev2)  # βρισκουμε το posterior του μοντελου 1 , χρησιμοποιω plogis γιατι μου βγαζει NaN αλλιως 
Postprob2 <- 1 - Postprob1            # και ομοιως για το 2
  draws<-rep(0,103)
  draws[1]=Postprob1
  draws[2]=Postprob2 #φτιαχνουμε τον πινακα μας 
  draws[3]=which.max(pr)
  for(i in 1:99){
    draws[i+3]=pr[i] # βοηθητικο για να δω την πιθανοτητα αλλαγης σε αλλα σημεια , σιγουρα υπαρχει καλυτερος τροπος απο αυτος
  }
  
  return(draws)}

draws<-DNA_project(chr1.2)
