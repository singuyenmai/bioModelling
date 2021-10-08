library(DDD)

dd_lamuN = function(ddmodel,pars,N)
{
    la = pars[1]
    mu = pars[2]
    K = pars[3]
    n0 = (ddmodel == 2 | ddmodel == 4)
    if(length(pars) == 4)
    {
        r = pars[4]
    }
    if(ddmodel == 1)
    {
        # linear dependence in speciation rate
        laN = max(0,la - (la - mu) * N/K)
        muN = mu
    }
    if(ddmodel == 1.3)
    {
        # linear dependence in speciation rate
        laN = max(0,la * (1 - N/K))
        muN = mu
    }
    if(ddmodel == 2 | ddmodel == 2.1 | ddmodel == 2.2)
    {
        # exponential dependence in speciation rate
        al = (log(la/mu)/log(K+n0))^(ddmodel != 2.2)
        laN = la * (N + n0)^(-al)
        muN = mu
    }
    if(ddmodel == 2.3)
    {
        # exponential dependence in speciation rate
        al = K
        laN = la * (N + n0)^(-al)
        muN = mu
    }
    if(ddmodel == 3)
    {
        # linear dependence in extinction rate
        laN = la
        muN = mu + (la - mu) * N/K
    }
    if(ddmodel == 4 | ddmodel == 4.1 | ddmodel == 4.2)
    {
        # exponential dependence in extinction rate
        al = (log(la/mu)/log(K+n0))^(ddmodel != 4.2)
        laN = la
        muN = mu * (N + n0)^al
    }
    if(ddmodel == 5)
    {
        # linear dependence in speciation rate and extinction rate
        laN = max(0,la - 1/(r+1)*(la-mu) * N/K)
        muN = mu + r/(r+1)*(la-mu)/K * N
    }
    return(c(laN,muN))
}

#' Function to simulate the diversity-dependent diversification process
#' 
#' Simulating the diversity-dependent diversification process
#' 
#' 
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to
#' lambda (speciation rate) \cr \code{pars[2]} corresponds to mu (extinction
#' rate) \cr \code{pars[3]} corresponds to K (clade-level carrying capacity)
#' @param age Sets the crown age for the simulation
#' @param ddmodel Sets the model of diversity-dependence: \cr \code{ddmodel ==
#' 1} : linear dependence in speciation rate with parameter K (= diversity
#' where speciation = extinction)\cr \code{ddmodel == 1.3} : linear dependence
#' in speciation rate with parameter K' (= diversity where speciation = 0)\cr
#' \code{ddmodel == 2} : exponential dependence in speciation rate with
#' parameter K (= diversity where speciation = extinction)\cr \code{ddmodel ==
#' 2.1} : variant of exponential dependence in speciation rate with offset at
#' infinity\cr \code{ddmodel == 2.2} : 1/n dependence in speciation rate\cr
#' \code{ddmodel == 2.3} : exponential dependence in speciation rate with
#' parameter x (= exponent)\cr \code{ddmodel == 3} : linear dependence in
#' extinction rate \cr \code{ddmodel == 4} : exponential dependence in
#' extinction rate \cr \code{ddmodel == 4.1} : variant of exponential
#' dependence in extinction rate with offset at infinity \cr \code{ddmodel ==
#' 4.2} : 1/n dependence in extinction rate with offset at infinity \cr
#' \code{ddmodel == 5} : linear dependence in speciation and extinction rate
#' @return \item{ out }{ A list with the following four elements: The first
#' element is the tree of extant species in phylo format \cr The second element
#' is the tree of all species, including extinct species, in phylo format \cr
#' The third element is a matrix of all species where \cr - the first column is
#' the time at which a species is born \cr - the second column is the label of
#' the parent of the species; positive and negative values only indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values only indicate whether the species belongs to the left or
#' right crown lineage \cr - the fourth column is the time of extinction of the
#' species. If this equals -1, then the species is still extant.\cr The fourth
#' element is the set of branching times of the tree of extant species.\cr }
#' @author Rampal S. Etienne, Si-Nguyen Mai
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#'  dd_sim_traits(c(3.0, 1.0), c(0.4, 0.4), c(0.1,20), 10) 
#' @export dd_sim_traits

update_time_probs <- function(la, trait_rates, pars, N, ddmodel, t){
    # find the time to next event
    laN = sapply(c(1,2), function(i) dd_lamuN(ddmodel, 
                                            c(la[i], pars), N[i])[1])
    muN = sapply(c(1,2), function(i) dd_lamuN(ddmodel, 
                                            c(la[i], pars), N[i])[2])
    
    denom = (laN[1] + muN[1] + trait_rates[1]) * N[1] + 
            (laN[2] + muN[2] + trait_rates[2]) * N[2]
    
    t_event = ifelse(denom == 0, Inf, t + stats::rexp(1,denom))

    rates = c(laN[1]*N[1], laN[2]*N[2], muN[1]*N[1], muN[2]*N[2], 
                trait_rates[1]*N[1], trait_rates[2]*N[2])
    event_probs = rates/sum(rates)
    
    return(list(t = t_event, p = event_probs))
}

dd_sim_traits = function(la, trait_rates, pars, age, ddmodel = 1)
{
# Simulation of diversity-dependent process
#  . start from crown age
#  . no additional species at crown node
#  . no missing species in present
# pars = [la mu K]
#  . la = speciation rate per species
#  . mu = extinction rate per species
#  . K = diversification carrying capacity
#  . r = ratio of diversity-dependence in extinction rate over speciation rate
# age = crown age
# ddmodel = mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate with parameter K
#  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 2.1: variant with offset at infinity
#  . ddmodel == 2.2: 1/n dependence in speciation rate
#  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x
#  . ddmodel == 3 : linear dependence in extinction rate
#  . ddmodel == 4 : exponential dependence in extinction rate
#  . ddmodel == 4.1: variant with offset at infinity
#  . ddmodel == 4.2: 1/n dependence in speciation rate
#  . ddmodel == 5 : linear dependence in speciation rate and in extinction rate
traits <- c(0, 1)
events_list <- c("spe_0", "spe_1", "ex_0", "ex_1", "trait_01", "trait_10")

done = 0
while(done == 0)
{
    # number of species N at time t
    # i = index running through t and N
    t = rep(0,1) # time to next event
    L = matrix(0,2,5)
    i = 1

    #N = 2
    # L = data structure for lineages,
    # . L[,1] = branching times
    # . L[,2] = index of parent species
    # . L[,3] = index of daughter species
    # . L[,4] = time of extinction
    # . L[,5] = trait of parent-daughter species
    # j = index running through L
    L[1,1:5] = c(0, 0, -1, -1, traits[1])
    L[2,1:5] = c(0,-1, 2, -1, traits[2])
    lin_list = list(c(-1), c(2))
    N = c(length(lin_list[[1]]), length(lin_list[[2]]))
    newest_L = 2
    
    # update time to next event & event probabilities
    updated = update_time_probs(la, trait_rates, pars, N, ddmodel, t[i])
    t[i + 1] = updated$t
    event_probs = updated$p

    while(t[i + 1] <= age)
    {
        i = i + 1
        event = sample(events_list, 1, prob = event_probs)
        event_id = which(events_list == event)
        k = ifelse(event_id%%2==1, 1, 2) # trait where event happens
        ranL = DDD::sample2(lin_list[[k]], 1)
        
        if (event_id <= 2){ # speciation

            N[k] = N[k] + 1
            newest_L = newest_L + 1
            L = rbind(L, 
                    c(t[i], ranL, sign(ranL)*newest_L, -1, traits[k]))
            lin_list[[k]] = c(lin_list[[k]], sign(ranL)*newest_L)

        } else if (event_id <= 4){ # extinction

            N[k] = N[k] - 1
            L[abs(ranL), 4] = t[i]
            w = which(lin_list[[k]] == ranL)
            lin_list[[k]] = sort( lin_list[[k]][-w] )

        } else { # traits switching

            h = ifelse(k==1, 2, 1) # new trait
            N[k] = N[k] - 1
            N[h] = N[h] + 1
            L[abs(ranL), 5] = traits[h]

            w = which(lin_list[[k]] == ranL)
            lin_list[[k]] = sort( lin_list[[k]][-w] )
            lin_list[[h]] = sort( c(lin_list[[h]], ranL) )

        }

        if(sum(lin_list[[1]] < 0, lin_list[[2]] < 0) == 0 | 
            sum(lin_list[[1]] > 0, lin_list[[2]] > 0) == 0)
        {
            t[i + 1] = Inf
        } else {
            updated = update_time_probs(la, trait_rates, pars, N, ddmodel, t[i])
            t[i + 1] = updated$t
            event_probs = updated$p
        } 
    }

    if(sum(lin_list[[1]] < 0, lin_list[[2]] < 0) == 0 | 
        sum(lin_list[[1]] > 0, lin_list[[2]] > 0) == 0)
    {
       done = 0
    } else {
       done = 1
    }
}

L[,1] = age - c(L[,1])
notmin1 = which(L[,4] != -1)
L[notmin1,4] = age - c(L[notmin1,4])
L[which(L[,4] == age + 1),4] = -1
tes = L2phylo(L,dropextinct = T)
tas = L2phylo(L,dropextinct = F)
brts = L2brts(L,dropextinct = T)
out = list(tes = tes,tas = tas,L = L,brts = brts)
return(out)

}
