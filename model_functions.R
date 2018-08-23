# Tranformation function; 
# para is . inside G(.)
# r is the parameter that determines function, as illustrated below.
G_fun <- function(para,r) {
  
  if (r <= 2) { # Box-Cox transformations; if r = 1, then PH model
  output = (1/r)*((1+para)^r-1)
  }
  
  else { # Logarithmic transfomations; if r = 3, then PO model
    output = (1/(r-2))*log(1+para*(r-2)) 
  }
  return(output)
}


# Bernstein polynomial function: j from 0 to m, where m is degree; 
# [l,u] is specified range of time;
# t for observed times
# output is a Bernstein polynomial as function of t
bern <- function(j,m,l,u,t){
  b = (factorial(m)/(factorial(j)*factorial(m-j)))*(((t-l)/(u-l))^j)*((1-(t-l)/(u-l))^(m-j))
  return(b)
}


