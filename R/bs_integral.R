

Ibs.1.fun = function(bs.knots, i, x)
{
  ((x <= bs.knots[i+1]) * pmax(x - bs.knots[i], 0)^2/
     (bs.knots[i+1]-bs.knots[i])/2
   + (x > bs.knots[i+1]) * ((bs.knots[i+2]-bs.knots[i])/2
                            - pmax(bs.knots[i+2]-x, 0)^2/
                              (bs.knots[i+2]-bs.knots[i+1])/2))
}

bs.1.fun = function(bs.knots, i, x)
{
  ((x <= bs.knots[i+1]) * pmax(x - bs.knots[i], 0)/
     (bs.knots[i+1]-bs.knots[i])
   +(x > bs.knots[i+1]) * pmax(bs.knots[i+2]-x, 0)/
     (bs.knots[i+2]-bs.knots[i+1]))
}


bs.1 = function(x, bs.knots, df)
{
  if(missing(bs.knots))
  {
    bs.knots = quantile(x, probs = seq(0,1, length.out = df+2))
  }else{
    df = length(bs.knots)-2
  }
  
  bs.x = matrix(0, length(x), df)
  Ibs.x = matrix(0, length(x), df)
  for (i in 1:df)
  {
    bs.x[,i] = bs.1.fun(bs.knots, i, x)
    
    Ibs.x[,i] = Ibs.1.fun(bs.knots, i, x)
    
  }
  
  return(list(bs = bs.x, 
              Ibs = Ibs.x, 
              knots = bs.knots))
 
}

a.liang = function(x, bs.knots, gamma.liang)
{
  a.x = rep(0, length(x))
  nknots = length(bs.knots)
  
  x.bin = which((x> 0) & (x<bs.knots[1]))
  if(length(x.bin) > 0)
  {
    a.x[x.bin] = exp(gamma.liang[1])* x[x.bin]
  }
  cum.I = exp(gamma.liang[1])*bs.knots[1]
  
  for (i in 1:(nknots-1))
  {
    a = gamma.liang[1]
    b = 0 
    if(i > 1)
    {
      a = a + bs.knots[i+1]*gamma.liang[i]/(bs.knots[i+1]-bs.knots[i])
      b = b - gamma.liang[i]/(bs.knots[i+1]-bs.knots[i])
    }
    if(i < nknots - 1)
    {
      a = a - bs.knots[i]*gamma.liang[i+1]/(bs.knots[i+1]-bs.knots[i])
      b = b + gamma.liang[i+1]/(bs.knots[i+1]-bs.knots[i])
    }
    x.bin = which((x>= bs.knots[i]) & (x<bs.knots[i+1]))
    if(length(x.bin) > 0)
    {
      a.x[x.bin] = (exp(a+b*x[x.bin]) - exp(a + b*bs.knots[i]))/b + cum.I
    }
    cum.I = cum.I + (exp(a+b*bs.knots[i+1]) - exp(a+ b*bs.knots[i]))/b
  }
  x.bin = which((x> bs.knots[nknots]))
  if(length(x.bin) > 0)
  {
    a.x[x.bin] = exp(gamma.liang[1])* (x[x.bin]-bs.knots[nknots]) + cum.I
  }
  
  return(a.x)
}


Ibs.1.pred = function(bs.knots,  x, lam)
{
  out = rep(0, length(x))
  for(i in 1:length(lam))
  {
    out = out + lam[i]* Ibs.1.fun(bs.knots,i,x)
  }
  return(out)
}


# 
# x = seq(0,1, 0.01)
# 
# plot(x, bs.x[,1], type = 'l')
# for (i in 2:df){
#   lines(x, bs.x[,i])
# }
# 
# plot(x, Ibs.x[,1], type = 'l')
# for (i in 2:df){
#   lines(x, Ibs.x[,i])
# }
# 
# bs.2 = function(x, df)
# {
#     
# }
