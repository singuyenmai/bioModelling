#########
# Thu Si Nguyen Mai - s4278836
#########

calculate_sum_of_rates <- function(pars,state,time)
{
  N0 <- state$N0
  N1 <- state$N1
  N <- state$N 
  ####
  # get parameters
  # `pars` is a vector of parameters in the following order
  l0 <- pars[1] # lambda_0: max per capita speciation rate
  a0 <- pars[2] # a_0: effect of trait 0 on per capita speciation rate
  a1 <- pars[3] # a_1: effect of trait 1 on per capita speciation rate
  mu <- pars[4] # extinction rate: constant & independent of trait
  q01 <- pars[5] # rate of switching from trait 0 to trait 1
  q10 <- pars[6] # rate of switching from trait 1 to trait 0
  
  # compute overal per capita speciation rate
  lambda <- max(0, l0*(1 - a0*N0 - a1*N1))

  sum_of_rates <- (lambda + mu)*N + q01*N0 + q10*N1
  ####
  return(sum_of_rates)
}

sample_event <- function(pars,state,time)
{
  N0 <- state$N0
  N1 <- state$N1
  N <- state$N
  ####
  # list of possible events
  events <- c('speciation in type 0', 'speciation in type 1',
              'extinction in type 0', 'extinction in type 1',
              'change in state from 0 to 1', 
              'change in state from 1 to 0')
  # get parameters
  # `pars` is a vector of parameters in the following order
  l0 <- pars[1] # lambda_0: max per capita speciation rate
  a0 <- pars[2] # a_0: effect of trait 0 on per capita speciation rate
  a1 <- pars[3] # a_1: effect of trait 1 on per capita speciation rate
  mu <- pars[4] # extinction rate: constant & independent of trait
  q01 <- pars[5] # rate of switching from trait 0 to trait 1
  q10 <- pars[6] # rate of switching from trait 1 to trait 0

  # compute overal per capita speciation rate
  lambda <- max(0, l0*(1 - a0*N0 - a1*N1))

  # event probabilities
  rates <- c(lambda*N0, lambda*N1, mu*N0, mu*N1, q01*N0, q10*N1)
  event_probs <- rates/sum(rates)

  # sample the event
  event <- sample(events, 1, prob = event_probs)
  ####
  return(event)
}

update_state <- function(state,event)
{
  N0 <- state$N0
  N1 <- state$N1
  N <- state$N
  
  #### Updating state
  if (event == 'speciation in type 0'){
    N0 <- N0 + 1
  } else if (event == 'speciation in type 1'){
    N1 <- N1 + 1
  } else if (event == 'extinction in type 0'){
    N0 <- N0 - 1
  } else if (event == 'extinction in type 1'){
    N1 <- N1 -1
  } else if (event == 'change in state from 0 to 1'){
    N0 <- N0 - 1
    N1 <- N1 + 1
  } else if (event == 'change in state from 1 to 0'){
    N0 <- N0 + 1
    N1 <- N1 - 1
  }

  N <- N0 + N1
  ####

  state <- list(N0 = N0,N1 = N1,N = N)
  return(state)
}

update_newL_and_L_and_linlist <- function(L,newL,ranL,linlist0,linlist1,time,event)
{
  if(event == 'speciation in type 0' | event == 'speciation in type 1')
  {
    newL <- newL + 1
    L <- rbind(L,c(time,ranL,sign(ranL) * newL,-1))
    if(event == 'speciation in type 0')
    {
      linlist0 <- c(linlist0,sign(ranL) * newL)
    } else
    if(event == 'speciation in type 1')
    {
      linlist1 <- c(linlist1,sign(ranL) * newL)
    }
  } else
  if(event == 'extinction in type 0' | event == 'extinction in type 1')
  {
    L[abs(ranL),4] <- time
    if(event == 'extinction in type 0')
    {
      w <- which(linlist0 == ranL)
      linlist0 <- linlist0[-w]
    } else
      if(event == 'extinction in type 1')
      {
        w <- which(linlist1 == ranL)
        linlist1 <- linlist1[-w]
      }
  } else
  if(event == 'change in state from 0 to 1')
  {
    w <- which(linlist0 == ranL)
    linlist0 <- linlist0[-w]
    linlist1 <- c(linlist1,ranL)
  } else
  if(event == 'change in state from 1 to 0')
  {
    w <- which(linlist1 == ranL)
    linlist1 <- linlist1[-w]
    linlist0 <- c(linlist0,ranL)
  }
  linlist0 <- sort(linlist0)
  linlist1 <- sort(linlist1)
  return(list(newL = newL,L = L,linlist0 = linlist0,linlist1 = linlist1))
}  

calculate_lineage_of_event <- function(linlist0,linlist1,event)
{
  if(event == 'speciation in type 0' | event == 'extinction in type 0' | event == 'change in state from 0 to 1')
  {
    ranL <- DDD::sample2(x = linlist0,size = 1)
  } else
  if(event == 'speciation in type 1' | event == 'extinction in type 1' | event == 'change in state from 1 to 0')
  {
    ranL <- DDD::sample2(x = linlist1,size = 1)
  }
  return(ranL)
}
  
model_sim <- function(pars,age)
{
  # pars contains all the parameters
  # age is the time you want to run the model for
  done <- 0
  while(done == 0)
  {
    t <- rep(0,1)
    L <- matrix(0,2,4)
    i <- 1
    t[1] <- 0
    state <- list(N0 = 1,N1 = 1,N = 2)
    L[1,1:4] <- c(0,0,-1,-1)
    L[2,1:4] <- c(0,-1,2,-1)
    linlist0 <- -1
    linlist1 <- 2
    linlist <- c(linlist0,linlist1)
    newL <- 2
    denom <- calculate_sum_of_rates(pars = pars,state = state,time = t[i])
    t[i + 1] <- t[i] + stats::rexp(1,denom)
    while(t[i + 1] <= age)
    {
      i <- i + 1
      event <- sample_event(pars = pars,state = state,time = t[i])
      ranL <- calculate_lineage_of_event(linlist0,linlist1,event)
      state <- update_state(state = state,event = event)
      newL_and_L_and_linlist <- update_newL_and_L_and_linlist(L = L,newL = newL,ranL = ranL,linlist0 = linlist0,linlist1 = linlist1,time = t[i],event = event)
      newL <- newL_and_L_and_linlist$newL
      L <- newL_and_L_and_linlist$L
      linlist0 <- newL_and_L_and_linlist$linlist0
      linlist1 <- newL_and_L_and_linlist$linlist1
      linlist <- sort(c(linlist0,linlist1))
      if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
      {
        t[i + 1] <- Inf
      } else {
        denom <- calculate_sum_of_rates(pars,state,t)
        t[i + 1] <- t[i] + stats::rexp(1,denom)
      } 
    }
    if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
    {
      done <- 0
    } else {
      done <- 1
    }
  }
  
  L[,1] <- age - c(L[,1])
  notmin1 <- which(L[,4] != -1)
  L[notmin1,4] <- age - c(L[notmin1,4])
  L[which(L[,4] == age + 1),4] <- -1
  tes <- DDD::L2phylo(L,dropextinct = T)
  tas <- DDD::L2phylo(L,dropextinct = F)
  brts <- DDD::L2brts(L,dropextinct = T)
  out <- list(tes = tes,tas = tas,L = L,brts = brts)
  return(out)
}
