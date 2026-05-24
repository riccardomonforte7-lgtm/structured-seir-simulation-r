# This assignment was completed collaboratively by 
# Mansi Neb (s2861549),
# Riccardo Monforte (s2882823), and
# Fakhruddin Hakim Hussain (s2845938).
# Each member contributed approximately one-third (33%) to the total effort,
# covering the setup, model implementation, simulation analysis, and documentation.

# ------------------------------------------------------------------------------#

# This code implements a refined SEIR (Susceptible, Exposed, Infectious, Recovered)
# epidemiological model that accounts for social structure. It contains functions
# to generate household groups and a contact network, which are then used
# in a simulation to study how these social structures (or lack thereof)
# influence epidemic dynamics.

# ------------------------------------------------------------------------------#

# We start by setting up the population's household structure. We generate
# a vector 'h' where each element is a household ID. People with the same
# household ID live in the same household.

n <- 1000       # Number of people in the population.
hmax <- 5       # The maximum number of people in a single household.

# Create household IDs by repeating 1 to n multiple times, sampled up to hmax
# and then truncating to the population size 'n'.
h <- rep(1:n, (sample(1:hmax, n, replace=TRUE)))[1:n]


# ------------------------------------------------------------------------------#

# Now, we create a function called 'get.net()'.

get.net <- function(beta, h, nc = 15) {

# Purpose: Generate a list-based contact network for the population.
# This network represents regular, non-household contacts based on individual
# 'sociability'.

# Inputs:
#   1) beta: numeric vector of sociability for each person 
#           (higher = more likely to contact others)
#   2) h: vector indicating household membership 
#         (people in the same household share the same value)
#   3) nc: average number of contacts per person (default = 15)

# Output:
# A list where the i-th element is a vector of the contact indices for person i.
  
  n <- length(beta)       # Number of individuals in the population.
  beta_bar <- mean(beta)  # Average sociability value across all individuals.
  
  # Calculate the probability of a contact existing for every possible pair of people.
  # The 'outer' function is a fast, vectorized way to apply the probability formula
  # to all pairs of sociability scores.
  p_mat <- outer(beta, beta, function(bi, bj) {
    nc * bi * bj / (beta_bar^2 * (n - 1))
  })
  
  # Ensure all calculated probabilities are valid (between 0 and 1).
  p_mat <- pmin(1, pmax(0, p_mat))
  
  # Create a logical matrix that identifies pairs of people who are in the same
  # household or are the same person.
  same_house_or_individual <- outer(h, h, "==")
  
  # Set the contact probability to 0 for pairs identified in the previous step.
  p_mat[same_house_or_individual] <- 0
  
  # Initialize an empty matrix to represent the network's connections,
  # where a '1' will denote a link.
  A <- matrix(0, n, n)
  
  # We only consider the upper half of the matrix because the network is symmetric
  # (if a person i is connected to person j, j is also connected to i).
  upper_tri_indices <- upper.tri(A)
  
  # Randomly generate contacts in the upper triangle based on the probability matrix.
  # A contact is formed if a random uniform value is less than the calculated probability.
  A[upper_tri_indices] <- runif(sum(upper_tri_indices)) < p_mat[upper_tri_indices]
  
  # Make the matrix symmetric. This ensures that if a link exists
  # from person i to j, it also exists from j to i, completing the network.
  A <- A | t(A)
  
  # Convert the final network matrix into a list format.
  net <- lapply(1:n, function(i) which(A[i, ] == 1))
  
  return(net)
}

# ------------------------------------------------------------------------------#

# Next, we create a function called 'nseir()'.

nseir <- function(beta, h, alink,
                  alpha = c(.1, .01, .01),  
                  delta = .2, gamma = .4, nc = 15,
                  nt = 100, pinf = .005) {

# Purpose: Run a stochastic SEIR simulation with social structure (households + network).

# Inputs:
#   beta  = individual sociability
#   h     = household IDs
#   alink = list of regular non-household contacts for each person
#   alpha = (alpha_h, alpha_c, alpha_r): daily infection probs 
#           for household, network, random contacts (default = 0.1, 0.01, 0.01)
#   gamma = daily prob E (Exposed) --> I (Infectious) (default = 0.4)
#   delta = daily prob I (Infectious) --> R (Recovered) (default = 0.2)
#   nc    = average contact per person (default = 15)
#   nt    = number of days to simulate (default = 100)
#   pinf  = proportion of infective at day 1 (default = 0.005)

# State encoding: S=0; E=1; I=2; R=3

# Output:
# A list with daily counts of S, E, I, R and the time index t.
  
  n <- length(beta)   # Number of individuals in the population.
  
  # Setting up alpha for each possible social scenario
  alpha_h <- alpha[1]
  alpha_c <- alpha[2]
  alpha_r <- alpha[3]
  
  # Initialize susceptibles
  x <-  rep(0, n)
  
  # Create initial infectives (I=2)
  nI0 <- round(pinf * n)
  x[sample.int(n, nI0)] <- 2 
  
  # Setting up storage for tracking pop state each day
  S <- E <- I <- R <-  rep(0, nt)
  
  # Initialize day 1 storage variables
  S[1] <- sum(x == 0)
  E[1] <- sum(x == 1)
  I[1] <- sum(x == 2)
  R[1] <- sum(x == 3)
  
  # Loop over days
  for (t in 2:nt) {
    
    # Set up next time step
    x_next <- x
    
    # Get position of people in each state (previous day)
    S_idx <- which(x == 0)
    E_idx <- which(x == 1)
    I_idx <- which(x == 2) 
    
    # E-->I transition (latency period ends) with probability gamma
    t_I <- runif(length(E_idx)) < gamma
    x_next[E_idx[t_I]] <- 2
    
    # I-->R transition (recovery) with probability delta
    t_R <- runif(length(I_idx)) < delta
    x_next[I_idx[t_R]] <- 3
    
    
    ##---------------- Scenario 1: infection from household--------------------##
    
    # Get number of infectives per household, then per susceptible
    hmax_val <- max(h)
    I_xfam <- tabulate(h[I_idx], nbins = hmax_val) 
    kxS_h <- I_xfam[h[S_idx]] 
    
    # Probability of being infected by a household member
    p_h <- 1 - (1-alpha_h)^kxS_h
    
    
    ##---------------- Scenario 2: infection from network ---------------------##
    
    # Flag infectives; count I among each susceptible’s contacts
    flag <- rep(FALSE, n)
    flag[I_idx] <- TRUE
    kxS_c <- sapply(alink[S_idx], function(i) sum(flag[i]))
    
    # Probability of being infected by a network member
    p_c <- 1 - (1-alpha_c)^kxS_c
    
    
    ##---------------- Scenario 3: random infection ---------------------------##
    
    # Isolate constant factor for pairwise random contacts
    beta_bar <- mean(beta)
    c0 <-  (alpha_r*nc)/(beta_bar^2*(n-1))
    
    # Compute scaled contact factor for susceptibles (includes sociability and base constant)
    beta_j_scaled <- c0*beta[S_idx] 
    
    # Per-infective sociability
    bi <- beta[I_idx] 
    
    # Pairwise transmission probs (rows = susceptibles, cols = infectives)
    M <- beta_j_scaled %o% bi
    
    # Compute log probability of "no-infection"
    # log1p(-M) is numerically stable for log(1-M)
    lsum <- rowSums(log1p(-M))
    
    # Convert cumulative "no-infection" log-probability to infection probability
    p_r <- 1-exp(lsum) 
    
    # Total probability of being infected assuming independence of the three scenarios
    p_tot <- 1- (1-p_h)*(1-p_c)*(1-p_r)
    
    # S --> E transition with total probability p_tot
    new_E <- runif(length(S_idx)) < p_tot
    x_next[S_idx[new_E]] <- 1
    
    # Commit day and update counts
    x <- x_next
    S[t]<-sum(x==0)
    E[t]<-sum(x==1)
    I[t]<-sum(x==2)
    R[t]<-sum(x==3)
  }
  
  # Return the results
  list(S=S, E=E, I=I, R=R, t=0:(nt-1) , beta = beta)
}

# ------------------------------------------------------------------------------#

# Now, we create function called 'plot_sim()' to generate 2 plots.

# Setting plot window for multiple plots.
par(mfcol=c(2,4),mar=c(4,4,3,2) + 0.1)

plot_sim <- function(seir_matrix , case){
  
  # Purpose: Generate two plots: a histogram of the beta (sociability) distribution 
  # and a time-series plot of the S, E, and I populations.
  
  # Inputs:
  #   seir_matrix: A list containing simulation results (S, E, I, R) and the 
  #                beta vector used, as returned by nseir().
  #   case: An integer (1-4) identifying the simulation case for plotting a title.
  
  # Output:
  # Two plots displayed in the current par() window.
  
  ## Defining title for each case
  title <- switch(case,
                  "1" = "Case 1: Variable Beta,\nAll Structures (Default)",
                  "2" = "Case 2: Variable Beta,\nNot including Household and regular network",
                  "3" = "Case 3: Constant Beta,\nAll Structures",
                  "4" = "Case 4: Constant Beta,\nNot including Household and regular network")
  
  ## Generating the location for the title of each case.
  at_params = (1:4 - 0.5) /4
  
  ## Plot 1: Histogram of beta distribution
  hist(seir_matrix$beta,xlab="Individual Sociability (β)", ylab="Frequency (Count)", 
       main = "" ,cex.lab = 0.8) ## beta distribution
  
  ## Adding grid lines to the histogram 
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  
  ## Plot 2: Time-series of SEIR population counts
  
  ## Plotting the S (Susceptible) vector
  plot(seir_matrix$S,ylim=c(0,max(seir_matrix$S)),xlab="Time (Days)",ylab="Population Count",
       main = "Evolution of the model" , 
       cex.main = 0.8 , cex.lab = 0.8) ## S black
  
  ## Plotting the E (Exposed) and I (Infectious) vectors
  points(seir_matrix$E,col=4);points(seir_matrix$I,col=2) ## E (blue) and I (red)
  
  ## Adding legend to the plot
  legend("topright", 
         legend = c("S", "E", "I"),
         col = c(1, 4, 2), # Black, Blue, Red
         lty = 1, # Line type
         lwd = 2, # Line width
         cex = 0.6) # Font size
  
  ## Adding grid lines to the plot
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  
  ## Adding the title for the particular case
  mtext(title, side = 3, line = -2, outer = TRUE, at = at_params[case], cex = 0.5 , font = 2)
  
}

# ------------------------------------------------------------------------------#

# Now, we create function called 'run_sim()'.

run_sim <- function(beta , h , 
                    alpha = c(0.1 , 0.01,0.01) , delta = 0.2 , gamma = 0.4 ,
                    nc = 15 , nt = 100 , pinf = 0.005){
  
  # Purpose: A function to execute the full simulation pipeline: 
  # create the network and then run the SEIR model.
  
  # Inputs:
  #   beta, h, alpha, delta, gamma, nc, nt, pinf: All parameters required by nseir().
  
  # Output:
  # The list of simulation results returned by nseir().
  
  # Build the social network using the defined function.
  # 'alink' now holds the list of regular contacts for each person.
  alink = get.net(beta , h , nc)
  
  # Running the SEIR model with the given parameters
  seir_matrix = nseir(beta , h , alink , alpha , delta , gamma , nc , nt , pinf)
  
  return (seir_matrix)
}

#--------------------------------------------------------------------#
# EXPERIMENT SETUP: Now we run the four main simulation scenarios as follows:
# (The 'case' variable is used to label the results for plotting and analysis)

## Case 1 : Default Parameters
## Generate beta from a uniform distribution and 
##keep the other parameters with default values.
case = 1
# Sociability vector (beta)
beta = runif(n , 0 , 1)
seir_matrix_case_1 = run_sim(beta = beta , h = h)
plot_sim(seir_matrix_case_1 , case)

#--------------------------------------------------------------------#

## Case 2 : Not including Household and regular network
## We set αh = αc = 0 and αr = 0.04 so that we don't include
## the Household and the regular network and using the beta generated in case 1.
case = 2
## setting αh = αc = 0 and αr = 0.04
alpha_new = c(0 , 0 , 0.04)
seir_matrix_case_2 = run_sim(beta = beta , h = h , alpha = alpha_new)
plot_sim(seir_matrix_case_2 , case)

#--------------------------------------------------------------------#

## Case 3 : Constant Beta
## In this case, we use a constant beta for every person.
## This beta is generated by taking the average of the betas
## generated in case 1 and assigning it for every person in case 3
## along with taking the default values for all the other parameters.        
case = 3
## Taking the average of the beta vector and repeating it 'n' times                    
beta_new = rep(mean(beta) , n)
seir_matrix_case_3 = run_sim(beta = beta_new , h)
plot_sim(seir_matrix_case_3 , case)

#--------------------------------------------------------------------#

## Case 4 : Constant Beta and not including Household and regular network
## This is the combination of case 2 and case 3 where we take a constant beta
## for each person(beta_new) and not include the household and regular network(alpha_new).
case = 4
seir_matrix_case_4 = run_sim(beta = beta_new , h = h , alpha = alpha_new)
plot_sim(seir_matrix_case_4 , case)

#--------------------Outcome comments---------------------------------#

## Case 1 : Default Parameters

# Households and social networks create smaller, localized clusters of infection.
# People tend to infect others within their families or regular contacts, so
# outbreaks rise and fade in small groups.         

## Case 2 : Not including Household and regular network

# Without social structure, the virus mixes more freely across the whole population,
# leading to a more synchronized and intense spread than in the default model. This
# produces a higher and sharper infection peak that also declines faster. In our runs,
# the number of susceptible people remaining at the end is often higher than when families
# and networks are present. In the structured model, households
# and networks slow down transmission by creating repeated and localized contacts,
# which reduce the effective spread of the infection. In the random-mixing scenario,
# infections reach a larger set of individuals, allowing the variability
# in beta to express its mitigating effect: populations with more heterogeneous beta
# values tend to experience a less extensive epidemic overall, leaving a larger share
# of susceptibles at the end. When social structure is present, repeated contacts within
# families and networks seem to limit this positive effect overall, leading to similar epidemic outcomes regardless
# of beta variability.

## Case 3 : Constant Beta

# When beta has no variability, everyone has the same sociability level. With
# households and networks, transmission remains localized: clusters are smaller
# and the epidemic unfolds more slowly. The infection peak is lower but lasts
# longer. The total number of infections is similar to the model with variable beta,
# showing that families and networks can significantly reduce the negative effect of low beta
# variability, even if they don’t provide much benefit when beta variance is high.

## Case 4 : Constant Beta and not including Household and regular network

# With constant beta and no social structure, the epidemic spreads rapidly through
# the entire population, as contacts are more diverse and less repetitive.
# The infection peak is higher and the final number of susceptibles is smaller
# compared to the structured model, meaning that the total number of infected
# individuals is greater. This confirms that, without clustering, homogeneous
# sociability leads to a faster and wider epidemic.
