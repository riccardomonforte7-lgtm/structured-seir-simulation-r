# SEIR Social Structure Simulation in R

This project implements a structured SEIR epidemic simulation model in R.

The aim was to explore how social structure affects epidemic dynamics, moving beyond a simple random mixing assumption. In particular, the model includes households, regular contact networks, and random population-level mixing.

## Overview

The project is based on a stochastic SEIR model, where each individual can be in one of four states:

- **S**: Susceptible
- **E**: Exposed
- **I**: Infectious
- **R**: Recovered

The model adds social structure in two ways:

1. **Households**  
   Individuals are grouped into households of different sizes. People in the same household have a higher probability of infecting each other.

2. **Regular contact network**  
   Each individual has a sociability parameter, which affects the probability of forming regular non-household contacts with other individuals.

In addition to household and network transmission, the model also includes a random mixing component, representing more general population-level contacts.

## Methodology

The code is written in base R and is organised around a few main functions:

- `get.net()`  
  Builds a regular contact network based on individual sociability values, excluding contacts between people in the same household.

- `nseir()`  
  Runs the stochastic SEIR simulation over a fixed number of days. The function tracks the number of susceptible, exposed, infectious, and recovered individuals over time.

- `plot_sim()`  
  Produces plots showing both the distribution of individual sociability values and the simulated epidemic dynamics.

- `run_sim()`  
  Combines the network generation and SEIR simulation steps into a single workflow.

The simulation considers three possible routes of infection:

- infection through household members
- infection through regular network contacts
- infection through random mixing

## Simulation scenarios

The project compares four scenarios:

1. **Variable sociability with full social structure**  
   The default model, including households, regular contact networks, and random mixing.

2. **Variable sociability with random mixing only**  
   Household and regular network transmission are removed, while random mixing is increased.

3. **Constant sociability with full social structure**  
   Every individual is assigned the same sociability value, while household and network structure remain in the model.

4. **Constant sociability with random mixing only**  
   Combines constant sociability with the removal of household and regular network transmission.

## Main insights

The simulations suggest that social structure can strongly affect epidemic dynamics.

Households and regular networks tend to create more localised clusters of infection. This can slow down transmission and produce epidemic patterns that differ from a pure random-mixing model.

When household and network structure are removed, infection spreads more freely across the population, often producing a sharper and more synchronised epidemic peak.

The comparison between variable and constant sociability also shows that individual heterogeneity can affect the final epidemic size, especially when contacts are less constrained by repeated household or network interactions.

## Limitations

This project is mainly a simulation exercise rather than an empirical modelling project.

Some important limitations are:

- the model does not use real epidemic data
- parameters are chosen for simulation purposes rather than estimated from data
- household sizes and contact networks are generated using simplified assumptions
- the model is stochastic, so results may vary between runs
- the simulation focuses on illustrating mechanisms rather than making real-world predictions

## Files

- `seir.R`: commented R script containing the full model implementation, simulation scenarios, plots, and outcome comments.
- `practical-2.pdf`: project brief describing the modelling task.

## How to run

Open the R script and run:

```r
source("seir.R")
