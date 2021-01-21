# Schistosomiasis model with an R interface running functions from Julia package Schistoxpkg


#' Setup Schistoxpkg
#'
#' This function initializes Julia and the Schistoxpkg.jl package
#' The first time will be long since it includes precompilation
#'
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @examples
#'
#' \donttest{
#' ## schistox_setup() is time-consuming and requires Julia
#'
#' }
#'
#' @export
schistox_setup <- function (...){
  julia <- JuliaCall::julia_setup(...)
  # JuliaCall::julia_install_package_if_needed("Distributions")
  # JuliaCall::julia_install_package_if_needed("Random")
  # JuliaCall::julia_install_package_if_needed("JLD")
  # JuliaCall::julia_install_package_if_needed("Schistoxpkg")
  JuliaCall::julia_library("Distributions")
  JuliaCall::julia_library("Random")
  JuliaCall::julia_library("JLD")
  JuliaCall::julia_library("Schistoxpkg")
}


#' Return the default parameters defined in the julia package
#'
#' @return parameters Julia structure
#' @examples \donttest{schistox_setup()
#' pars = default_pars()}
#' @export
default_pars <- function(){
  pars =  JuliaCall::julia_eval("Parameters()")
  return(  pars)
}



#' Save humans, environmental cercariae and mircacidia to a file along with the parameters used
#'
#' @param filename  file to save to
#' @param humans  humans variable
#' @param miracidia  environmental miracidia
#' @param cercariae  environmental miracidia
#' @param pars  parameter set
#' @examples \donttest{
#'  ## save_population_to_file(filename, humans, miracidia, cercariae, pars)
#' }
#' @export
save_population_to_file <- function(filename, humans, miracidia, cercariae, pars){
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("humans", humans)
  JuliaCall::julia_assign("miracidia", miracidia)
  JuliaCall::julia_assign("pars", pars)
  JuliaCall::julia_assign("cercariae", cercariae)

  JuliaCall::julia_eval("save_population_to_file(filename, humans, miracidia, cercariae, pars)")
}


#' output information for burden and heavy-intensity burden in school age children from julia recorded data
#'
#' @param record variable which stores the data
#'
#' @export
get_sac_data_from_record <- function( record){

  JuliaCall::julia_assign("record", record)

  x = JuliaCall::julia_eval("
    function get_sac_data_from_record(record)
    sac_burden = (p->p.sac_burden[1]).(record)
    sac_heavy_burden = (p->p.sac_burden[3]).(record)
    times = (p->p.time).(record);
    return sac_burden, sac_heavy_burden, times
  end
")

  return(x(record))
}



#' output information for burden and heavy-intensity burden in adults from julia recorded data
#'
#' @param record variable which stores the data
#'
#' @export
get_adult_data_from_record <- function( record){

  JuliaCall::julia_assign("record", record)

  x = JuliaCall::julia_eval("
    function get_adult_data_from_record(record)
    adult_burden = (p->p.adult_burden[1]).(record)
    adult_heavy_burden = (p->p.adult_burden[3]).(record)
    times = (p->p.time).(record);
    return adult_burden, adult_heavy_burden, times
  end
")

  return(x(record))
}




#' output information for burden and heavy-intensity for whole population from julia recorded data
#'
#' @param record variable which stores the data
#'
#' @export
get_pop_data_from_record <- function( record){

  JuliaCall::julia_assign("record", record)

  x = JuliaCall::julia_eval("
    function get_adult_data_from_record(record)
    pop_burden = (p->p.population_burden[1]).(record)
    pop_heavy_burden = (p->p.population_burden[3]).(record)
    times = (p->p.time).(record);
    return pop_burden, pop_heavy_burden, times
  end
")

  return(x(record))
}





#' Calculate the element wise mean of whatever Julia array is input
#'
#' @param a - input object to take element wise mean of
#'
#' @export
get_dot_mean<- function(a){
  JuliaCall::julia_assign("a", a)
  JuliaCall::julia_eval("mean.(a)")
}

#' Set the parameters to the desired values.
#' @param N size of human population (Int)
#' @param time_step length of time in days to step forward each time (Float)
#' @param N_communities number of communities sharing the same environmental source (Int)
#' @param community_probs probability of an integer being in a given community. This governs
#' the population sizes in each community. This must be an array the length of
#' N_communities (Array(Float))
#' @param community_contact_rate contact rate with the environment for each of the communities.
#' This must be an array the length of N_communities (Array(Float))
#' @param density_dependent_fecundity decrease in egg production per worm due to high density
#' of worms (Float)
#' @param average_worm_lifespan average life expectancy of a worm (Float)
#' @param max_age maximum age of human in the population (Float)
#' @param initial_worms number of worms in each person to begin with. The actual number in
#' each person if chosen from a Poisson distirubiton with this mean (Int)
#' @param initial_miracidia initial number of miracidia larvae in the environment (Int)
#' @param initial_miracidia_days miracidia will age into cercariae larvae after a specified
#' number of days. This parameter will specify how many days of initial_miracidia we
#' will have already had in the environment (Int)
#' @param init_env_cercariae initial number of cercariae larvae in the environment (Int)
#' @param worm_stages how many age stages are there for the worms. Having 1 stage will give
#' Gamma distributed death ages, while more than 1 will result in Erlang distribution (Int)
#' @param contact_rate global contact rate for the uptake of larvae from the environment (Float)
#' @param max_fecundity expected number of eggs from a single worm pair. The actual number
#' will be chosen from a distribution with this mean (Float)
#' @param age_contact_rates contact rate for chosen age groups(Array(Float))
#' @param ages_for_contacts age groups for specifying contact rates (Array(Int))
#' @param contact_rate_by_age_array array holding contact rate for each age (Array(Float))
#' @param mda_adherence proportion of people who adhere to the mda (Float)
#' @param mda_access proportion of people who have access to the mda (Float)
#' @param female_factor factor for altering the contact rate for females, if we choose to
#' have gender specific behaviour which affects contact rate (Float)
#' @param male_factor factor for altering the contact rate for males, if we choose to have
#' gender specific behaviour which affects contact rate (Float)
#' @param miracidia_maturity number of days after which miracidias will mature to cercariae (Int)
#' @param birth_rate rate of birth of humans (Float)
#' @param human_cercariae_prop proportion of cercariae which are able to infect humans (Float)
#' @param predis_aggregation aggregation for predisposition of individuals to uptake larvae.
#' This is chosen from a Gamma distribution with mean 1 for each individual and set for
#' life. If this is high, then the aggregation is low, meaning that most individuals
#' have roughly the same predisposition. If it is low, then the larvae become
#' concentrated in a few individuals. (Float)
#' @param cercariae_survival what proportion of cercariae survive from one time point
#' to the next (Float)
#' @param miracidia_survival what proportion of miracidia survive from one time point
#' to the next (Float)
#' @param death_prob_by_age for specified age range, what is the probability of dying
#' each year (Array(Float))
#' @param ages_for_death age ranges for death probabilities (Array(Float))
#' @param r aggregation parameter for negative binomially distributed egg production (Float)
#' @param vaccine_effectiveness efficacy of a vaccine if one is used(Float)
#' @param drug_effectiveness efficacy of a drug given during MDA(Float)
#' @param spec_ages number of individuals by age group which we specify if we want a
#' particular age distribution for the simulation (Array(Float))
#' @param ages_per_index how many different ages we include in the spec_ages parameter (Int)
#' @param record_frequency how often we should record the prevalence in the population
#' during simulation (Float)
#' @param use_kato_katz if 0, then don't use Kato-Katz (KK) for egg counts, if 1, use KK (Int)
#' @param kato_katz_par parameter for Gamma distribution if KK is used (Float)
#' @param heavy_burden_threshold number of eggs at which an individual is said to have a
#' heavy infection (Int)
#' @param  rate_acquired_immunity rate at which immunity will be acquired for individuals.
#' This will be multiplied by the cumulative number of worms people have had throughout
#' their life to decide the level of immunity acquired (Float)
#' @param M0 if a particular for of egg production is used, this parameter is required and is
#' a proxy for mean worm burden (Float)
#' @param human_larvae_maturity_time length of time in days after which a cercariae uptaken by
#' a human will mature into a worm (Int)
#' @param egg_sample_size the proportion of eggs which are sampled from each individual every time we
#' check their burden. This is between 0 and 1, with 1 meaning that all the eggs in the person are
#' sampled, hence giving the true burden. Typical value for a urine sample may be ~ 1/100
#' @param input_ages input ages for constructing contact array
#' @param input_contact_rates input rates for constructing contact array
#' @param scenario can be one of "low adult", "moderate adult", "high adult"
#' @export
set_pars <- function(N, time_step, N_communities, community_probs,
                     community_contact_rate, density_dependent_fecundity,
                     average_worm_lifespan, max_age, initial_worms,
                     initial_miracidia, initial_miracidia_days, init_env_cercariae,
                     worm_stages, contact_rate, max_fecundity, age_contact_rates,
                     ages_for_contacts, contact_rate_by_age_array, mda_adherence,
                     mda_access, female_factor, male_factor, miracidia_maturity,
                     birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival,
                     miracidia_survival, death_prob_by_age, ages_for_death, r,
                     vaccine_effectiveness, drug_effectiveness, spec_ages, ages_per_index,
                     record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold,
                     rate_acquired_immunity, M0, human_larvae_maturity_time, egg_sample_size, input_ages, input_contact_rates,
                     scenario){

  JuliaCall::julia_assign("N", N)
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("N_communities", N_communities)
  JuliaCall::julia_assign("community_probs", community_probs)
  JuliaCall::julia_assign("community_contact_rate", community_contact_rate)
  JuliaCall::julia_assign("density_dependent_fecundity", density_dependent_fecundity)
  JuliaCall::julia_assign("average_worm_lifespan", average_worm_lifespan)
  JuliaCall::julia_assign("max_age", max_age)
  JuliaCall::julia_assign("initial_worms", initial_worms)
  JuliaCall::julia_assign("initial_miracidia", initial_miracidia)
  JuliaCall::julia_assign("initial_miracidia_days", initial_miracidia_days)
  JuliaCall::julia_assign("init_env_cercariae", init_env_cercariae)
  JuliaCall::julia_assign("average_worm_lifespan", average_worm_lifespan)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("contact_rate", contact_rate)
  JuliaCall::julia_assign("max_fecundity", max_fecundity)
  JuliaCall::julia_assign("age_contact_rates", age_contact_rates)
  JuliaCall::julia_assign("ages_for_contacts", ages_for_contacts)
  JuliaCall::julia_assign("contact_rate_by_age_array", contact_rate_by_age_array)
  JuliaCall::julia_assign("mda_adherence", mda_adherence)
  JuliaCall::julia_assign("mda_access", mda_access)
  JuliaCall::julia_assign("female_factor", female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("miracidia_maturity", miracidia_maturity)
  JuliaCall::julia_assign("birth_rate", birth_rate)
  JuliaCall::julia_assign("human_cercariae_prop", human_cercariae_prop)
  JuliaCall::julia_assign("predis_aggregation", predis_aggregation)
  JuliaCall::julia_assign("cercariae_survival", cercariae_survival)
  JuliaCall::julia_assign("miracidia_survival", miracidia_survival)
  JuliaCall::julia_assign("death_prob_by_age", death_prob_by_age)
  JuliaCall::julia_assign("ages_for_death", ages_for_death)
  JuliaCall::julia_assign("r", r)
  JuliaCall::julia_assign("vaccine_effectiveness", vaccine_effectiveness)
  JuliaCall::julia_assign("drug_effectiveness", drug_effectiveness)
  JuliaCall::julia_assign("spec_ages", spec_ages)
  JuliaCall::julia_assign("ages_per_index", ages_per_index)
  JuliaCall::julia_assign("record_frequency", record_frequency)


  JuliaCall::julia_assign("use_kato_katz",use_kato_katz)
  JuliaCall::julia_assign("kato_katz_par",kato_katz_par)
  JuliaCall::julia_assign("heavy_burden_threshold", heavy_burden_threshold)

  JuliaCall::julia_assign("rate_acquired_immunity", rate_acquired_immunity)
  JuliaCall::julia_assign("M0", M0)

  JuliaCall::julia_assign("human_larvae_maturity_time", human_larvae_maturity_time)
  JuliaCall::julia_assign("egg_sample_size", egg_sample_size)
  JuliaCall::julia_assign("scenario", scenario)

 pars = JuliaCall::julia_eval("Parameters(N, time_step, N_communities, community_probs, community_contact_rate,
        density_dependent_fecundity, average_worm_lifespan,
        max_age, initial_worms, initial_miracidia, initial_miracidia_days, init_env_cercariae,
        worm_stages, contact_rate, max_fecundity, age_contact_rates,
        ages_for_contacts, contact_rate_by_age_array, mda_adherence, mda_access,  female_factor, male_factor, miracidia_maturity,
        birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival, miracidia_survival,
        death_prob_by_age, ages_for_death, r, vaccine_effectiveness, drug_effectiveness,
        spec_ages, ages_per_index, record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold,
        rate_acquired_immunity, M0, human_larvae_maturity_time, egg_sample_size)")



  pars = make_age_contact_rate_array(pars, scenario, input_ages, input_contact_rates)

  return(pars)


}




#' Make the contact by age array given the specified inputs
#'
#' @param pars  parameters we are using
#' @param scenario  can be one of "low adult", "moderate adult" and "high adult"
#' @param input_ages  if scenario not specified, we can input ages at which we wish contact rates to change. Ages must be integers, and input with as.integer()
#' @param input_contact_rates  the contact rates corresponding to the ages in the input_ages variable
#'
#' @return will return the parameters object
#' @examples
#' \donttest{
#' # define input ages
#' input_ages = array(data = c(as.integer(4), as.integer(9), as.integer(15), as.integer(120)))
#'
#'# define contact rates for age groups
#' input_contact_rates = array(data = c(0.64, 0.91, 1, 0.018))
#'
#' # normalize so that the sum of the contact rates is 1
#' input_contact_rates = input_contact_rates/sum(input_contact_rates)
#'
#' # make the contact by age array within the parameters object
#' #pars = default_pars()
#' #pars = make_age_contact_rate_array(pars, scenario, input_ages, input_contact_rates)
#'}
#' @export
make_age_contact_rate_array <- function(pars, scenario, input_ages, input_contact_rates){


  JuliaCall::julia_assign("pars", pars)
  JuliaCall::julia_assign("scenario", scenario)
  JuliaCall::julia_assign("input_ages", input_ages)
  JuliaCall::julia_assign("input_contact_rates", input_contact_rates)
  #JuliaCall::julia_eval("result = make_age_contact_rate_array")
  update_parameters_individually(pars, name = "age_contact_rates", value = input_contact_rates)
  result <- JuliaCall::julia_eval("make_age_contact_rate_array(pars, scenario, input_ages, input_contact_rates)")

  return(result)

}


#' Will age the human population in order to get a realistic age structure based on input death rates by age.
#'
#' @param num_steps  how many steps forward to take
#' @param humans  the humans
#' @param pars  parameters used
#' @return will return the human population container
#' @export
generate_ages_and_deaths <- function(num_steps, humans, pars){
  JuliaCall::julia_assign("num_steps", num_steps)
  JuliaCall::julia_assign("humans", humans)
  JuliaCall::julia_assign("pars", pars)

  humans <- JuliaCall::julia_eval("generate_ages_and_deaths(num_steps, humans, pars)")

  return(humans)
}


#' Create contact by age rates given a scenario
#'
#' @param scenario  one of "low adult", "moderate adult" and "high adult"
#'
#' @export
create_contact_settings <- function(scenario){
  JuliaCall::julia_assign("scenario", scenario)
  result <- JuliaCall::julia_eval("create_contact_settings(scenario)")
}


#' Create the population to begin from
#'
#' @param pars parameters input
#'
#' @return list(humans, miracidia, cercariae)
#' @export
create_population_specified_ages <- function(pars){

  JuliaCall::julia_assign("pars", pars)
  x = JuliaCall::julia_eval("create_population_specified_ages(pars)")

  env = list("humans" = x[[1]], "miracidia" = x[[2]], "cercariae" = x[[3]])
  return(env)
}


#' Update the parameter set used for the simulations This will take the existing parameters and update contact_rate, max_fecundity, predis_aggregation and density_dependent_fecundity
#'
#' @param pars  existing parameters
#' @param contact_rate  updated contact rate
#' @param max_fecundity  updated max fecundity
#' @param predis_aggregation  updated aggregation of predisposition
#' @param density_dependent_fecundity  updated density dependent fecundity
#'
#' @export
update_parameters <- function(pars, contact_rate, max_fecundity, predis_aggregation, density_dependent_fecundity){

  JuliaCall::julia_assign("pars", pars)
  JuliaCall::julia_assign("contact_rate", contact_rate)
  JuliaCall::julia_assign("max_fecundity", max_fecundity)
  JuliaCall::julia_assign("predis_aggregation", predis_aggregation)
  JuliaCall::julia_assign("density_dependent_fecundity", density_dependent_fecundity)
  x = JuliaCall::julia_eval("pars.contact_rate = contact_rate")
  x = JuliaCall::julia_eval("pars.max_fecundity = max_fecundity")
  x = JuliaCall::julia_eval("pars.density_dependent_fecundity = density_dependent_fecundity")
  x = JuliaCall::julia_eval("pars.predis_aggregation = predis_aggregation")


}


#' update a specific parameter
#'
#' @param pars the parameters object
#' @param name the name of the parameter. This must match the name of a parameter defined in the set_pars function and be input as a string.
#' @param value the value to update the parameter to
#'
#' @export
update_parameters_individually <- function(pars, name, value){

  JuliaCall::julia_assign(name, value)
  julia_code = paste("pars.", name, " = ", name, sep = "")
  x = JuliaCall::julia_eval(julia_code)
}



#' Update a number of specified parameters at the same time
#'
#' @param pars the parameters object
#' @param ... list of parameters and values for chosen parameters. Must be paired, with the name of the parameter in "" first, followed by the desired valeu.
#' For example, to update parameter N to 500 and max_fecunidty to 4  the function entry must look
#' like: update_specified_parameters(pars, "N", 500, "max_fecundity", 4)
#'
#'
#' @export
#'
update_specified_parameters <- function(pars, ...){

  for(i in seq(1,(nargs()-1),2)){
    name = list(...)[[i]]
    value = list(...)[[(i+1)]]
    JuliaCall::julia_assign(name, value)
    julia_code = paste("pars.", name, " = ", name, sep = "")
    x = JuliaCall::julia_eval(julia_code)
  }
}



#' Run the simulation of the population where we include births and deaths, but each death is matched by one birth
#'
#' @param num_time_steps  how many time steps to run forward for
#' @param humans human population
#' @param miracidia  environmental miracidia
#' @param cercariae  environmental cercariae
#' @param pars parameters
#' @param mda_info  mda's to enact
#' @param vaccine_info  vaccinations to enact
#'
#' @return A list variable is returned containing the humans array, miracidia, cercariae and the record of statistics
#' over the course of the run to equilibrium. To access the humans array, if we output from this function to a variable
#' env, then by running env$humans, we will return the humans array.
#' @export
update_env_constant_population<- function(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){
  JuliaCall::julia_assign("num_time_steps",num_time_steps)
  JuliaCall::julia_assign("humans", humans)
  JuliaCall::julia_assign("miracidia", miracidia)
  JuliaCall::julia_assign("cercariae", cercariae)
  JuliaCall::julia_assign("pars", pars)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info", vaccine_info)


  x =
    JuliaCall::julia_eval("update_env_constant_population(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)")

  env = list("humans" = x[[1]], "miracidia" = x[[2]], "cercariae" = x[[3]], "record" = x[[4]])
  return(env)
}





#' Run the simulation of the population where we include births and deaths but each death is matched by one birth. Here the egg_production function is monotonically increasing
#'
#' @param num_time_steps  how many time steps to run forward for
#' @param humans  human population
#' @param miracidia  environmental miracidia
#' @param cercariae  environmental cercariae
#' @param pars  parameters
#' @param mda_info  mda's to enact
#' @param vaccine_info  vaccinations to enact
#'
#' @return A list variable is returned containing the humans array, miracidia, cercariae and the record of statistics
#' over the course of the run to equilibrium. To access the humans array, if we output from this function to a variable
#' env, then by running env$humans, we will return the humans array.
#'
#' @export
update_env_constant_population_increasing<- function(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){
  JuliaCall::julia_assign("num_time_steps",num_time_steps)
  JuliaCall::julia_assign("humans", humans)
  JuliaCall::julia_assign("miracidia", miracidia)
  JuliaCall::julia_assign("cercariae", cercariae)
  JuliaCall::julia_assign("pars", pars)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info", vaccine_info)


  x = JuliaCall::julia_eval("update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)")
  env = list("humans" = x[[1]], "miracidia" = x[[2]], "cercariae" = x[[3]], "record" = x[[4]])
  return(env)
}



#' Run the simulation of the population where we include births and deaths, but each death is matched by one birth.
#' Here cercariae are picked up as larvae within humans, rather than immediately becoming worms.
#'
#' @param num_time_steps  how many time steps to run forward for
#' @param humans  human population
#' @param miracidia  environmental miracidia
#' @param cercariae  environmental cercariae
#' @param pars  parameters
#' @param mda_info  mda's to enact
#' @param vaccine_info  vaccinations to enact
#'
#' @return A list variable is returned containing the humans array, miracidia, cercariae and the record of statistics
#' over the course of the run to equilibrium. To access the humans array, if we output from this function to a variable
#' env, then by running env$humans, we will return the humans array.
#' @export
#'
update_env_constant_population_human_larvae <- function(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){
  JuliaCall::julia_assign("num_time_steps",num_time_steps)
  JuliaCall::julia_assign("humans", humans)
  JuliaCall::julia_assign("miracidia", miracidia)
  JuliaCall::julia_assign("cercariae", cercariae)
  JuliaCall::julia_assign("pars", pars)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info", vaccine_info)


  x = JuliaCall::julia_eval("update_env_constant_population_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)")
  env = list("humans" = x[[1]], "miracidia" = x[[2]], "cercariae" = x[[3]], "record" = x[[4]])
  return(env)
}



#' Update the contact rates for individuals when they age
#'
#' @param humans - human population
#' @param pars - parameters
#' @export
update_contact_rate <- function(humans, pars){

  JuliaCall::julia_assign("humans", humans)
  JuliaCall::julia_assign("pars", pars)

  x = JuliaCall::julia_eval("update_contact_rate(humans, pars)")

}


#' create the MDA
#'
#' @param pre_SAC_prop  proportion of pre-SAC individuals to recieve MDA
#' @param SAC_prop  proportion of SAC individuals to recieve MDA
#' @param adult_prop  proportion of adults individuals to recieve MDA
#' @param first_mda_time  specifies when this will first occur in years
#' @param last_mda_time  specifies when this will last occur in years
#' @param regularity  how regularly the mda will repeat in years
#' @param pre_SAC_gender  gender of pre-SAC individuals to recieve MDA
#' @param SAC_gender  gender of SAC individuals to recieve MDA
#' @param adult_gender  gender of adults individuals to recieve MDA
#' @param mda_effectiveness  efficacy of drug used.
#' @export
create_mda <- function(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
                       last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness){

  JuliaCall::julia_assign("pre_SAC_prop", pre_SAC_prop)
  JuliaCall::julia_assign("SAC_prop", SAC_prop)
  JuliaCall::julia_assign("adult_prop", adult_prop)
  JuliaCall::julia_assign("first_mda_time", first_mda_time)
  JuliaCall::julia_assign("last_mda_time", last_mda_time)
  JuliaCall::julia_assign("regularity", regularity)
  JuliaCall::julia_assign("pre_SAC_gender", pre_SAC_gender)
  JuliaCall::julia_assign("SAC_gender", SAC_gender)
  JuliaCall::julia_assign("adult_gender", adult_gender)
  JuliaCall::julia_assign("mda_effectiveness",  mda_effectiveness)


  mda_info = JuliaCall::julia_eval("create_mda(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
           last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness)")

  return(mda_info)
}



#' run simulations where the population doesn't change
#'
#' @param filename name of file to store population in
#' @param num_time_steps how many time steps to step forward
#' @param mda_info container for all MDAs to take place
#' @param vaccine_info container for all vaccination programs to take place
#' @param num_repeats how many times do we repeat the simulation
#' @export
run_repeated_sims_no_population_change<- function(filename, num_time_steps, mda_info, vaccine_info, num_repeats){
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info", vaccine_info)
  JuliaCall::julia_assign("num_repeats", num_repeats)

  x =
    JuliaCall::julia_eval("run_repeated_sims_no_population_change(filename, num_time_steps, mda_info, vaccine_info, num_repeats)")


  env = list("times" = x[[1]], "prev" = x[[2]], "sac_prev" = x[[3]], "high_burden" = x[[4]],
             "high_burden_sac" = x[[5]], "adult_prev" = x[[6]], "high_adult_burden" = x[[7]])

  return(env)
}


#' run simulations where the population doesn't change. Larvae are uptaken to larvae within humans.
#'
#' @param filename name of file to store population in
#' @param num_time_steps how many time steps to step forward
#' @param mda_info container for all MDAs to take place
#' @param vaccine_info container for all vaccination programs to take place
#' @param num_repeats how many times do we repeat the simulation
#' @export
run_repeated_sims_no_population_change_human_larvae <- function(filename, num_time_steps, mda_info, vaccine_info, num_repeats){
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info", vaccine_info)
  JuliaCall::julia_assign("num_repeats", num_repeats)

  x =
    JuliaCall::julia_eval("run_repeated_sims_no_population_change_human_larvae(filename, num_time_steps, mda_info, vaccine_info, num_repeats)")


  env = list("times" = x[[1]], "prev" = x[[2]], "sac_prev" = x[[3]], "high_burden" = x[[4]],
             "high_burden_sac" = x[[5]], "adult_prev" = x[[6]], "high_adult_burden" = x[[7]])

  return(env)
}





#' run simulations where the population doesn't change and the egg production function is monotonically increasing
#'
#' @param filename name of file to store population in
#' @param num_time_steps how many time steps to step forward
#' @param mda_info container for all MDAs to take place
#' @param vaccine_info container for all vaccination programs to take place
#' @param num_repeats how many times do we repeat the simulation
#'
#' @export
run_repeated_sims_no_population_change_increasing<- function(filename, num_time_steps, mda_info, vaccine_info, num_repeats){
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info", vaccine_info)
  JuliaCall::julia_assign("num_repeats", num_repeats)

  x =
    JuliaCall::julia_eval("run_repeated_sims_no_population_change_increasing(filename, num_time_steps, mda_info, vaccine_info, num_repeats)")

  env = list("times" = x[[1]], "prev" = x[[2]], "sac_prev" = x[[3]], "high_burden" = x[[4]],
             "high_burden_sac" = x[[5]], "adult_prev" = x[[6]], "high_adult_burden" = x[[7]])

  return(env)
}




#' add to MDA schedule
#'
#' @param mda_info current MDA information
#' @param mda_start_time start time for new MDAs
#' @param last_mda_time end time for new MDAs
#' @param regularity how regularly the MDA takes place in years (1 is once a year, 0.5 is twice a year)
#' @param drug_efficacy how effective the drug is
#' @param pre_SAC_prop what proportion of pre school age children receive MDA
#' @param SAC_prop what proportion of school age children receive MDA
#' @param adult_prop what proportion of adults receive MDA
#'
#' @export
add_to_mda <- function(mda_info, mda_start_time, last_mda_time,  regularity,
                       drug_efficacy, pre_SAC_prop, SAC_prop, adult_prop){

  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("mda_start_time", mda_start_time)
  JuliaCall::julia_assign("last_mda_time", last_mda_time)
  JuliaCall::julia_assign("drug_efficacy", drug_efficacy)
  JuliaCall::julia_assign("pre_SAC_prop", pre_SAC_prop)
  JuliaCall::julia_assign("SAC_prop", SAC_prop)
  JuliaCall::julia_assign("adult_prop", adult_prop)
  JuliaCall::julia_assign("regularity", regularity)

  outputs = JuliaCall::julia_eval("append!(mda_info, create_mda(pre_SAC_prop, SAC_prop, adult_prop, mda_start_time,
                               last_mda_time, regularity, [0,1], [0,1], [0,1], drug_efficacy))")


}






#' function to calculate the number of worm pairs
#'
#'
#' @param female_worms number of female worms for each individual
#' @param male_worms number of male worms for each individual
#' @export
calculate_worm_pairs <- function(female_worms, male_worms){
  f_worms = array(0, length(female_worms))
  m_worms = array(0, length(male_worms))
  worm_pairs = array(0, length(male_worms))
  for(i in 1 :length(f_worms)){
    f_worms[i] = sum(female_worms[i])
    m_worms[i] = sum(male_worms[i])
    worm_pairs[i] = min(f_worms[i], m_worms[i])
  }
  return(worm_pairs)
}



#' Get all entries of a chosen variable from the Julia humans array
#'
#' @param humans Julia array containing information about human population
#' @param name name of variable we want to return.
#' Must match exactly a name of a variable e.g. "age", "eggs", "female_worms", "male_worms"
#'
#' @return Will return an r array
#' @export
get_selected_data <- function(humans, name){
  julia_code = paste("function get_selected_data(humans)
                        x = (p->p.", name, ").(humans)
                        return x
                        end", sep = "")
  x = JuliaCall::julia_eval(julia_code)
  y = x(humans)
  g = array(0,dim = c(length(y), length(y[1])))
  for(i in 1:length(y)){
    g[i,] = y[i]
  }
  return(g)
}



#' Get the value of a chosen parameter
#'
#' @param pars the parameter object currently in use
#' @param name name of variable we want to return.
#' Must match exactly a name of a variable in parameters e.g. "contact_rate", "N" etc. See ?set_pars for list of names
#'
#' @return
#' @export
parameter_value <- function(pars, name){
  julia_code = paste("function parameter_value(pars)
                        x = pars.", name, "
                        return x
                        end", sep = "")
  x = JuliaCall::julia_eval(julia_code)
  y = x(pars)
  return(y)
}






get_ages_eggs_worms <- JuliaCall::julia_eval("
function get_ages_eggs_worms(humans)
    ages = (p->p.age).(humans)
    eggs = (p->p.eggs).(humans)
    female_worms = (p->p.female_worms).(humans)
    male_worms = (p->p.male_worms).(humans)
    return ages, eggs, female_worms, male_worms
  end
")






#' Administer drugs to specified individuals.
#'
#' @param humans human population
#' @param indices which individuals to administer drug to
#' @param drug_effectiveness how effective the drug is
#' @export
administer_drug <- function(humans, indices, drug_effectiveness){



  JuliaCall::julia_assign("humans", humans)
  JuliaCall::julia_assign("indices", indices)
  JuliaCall::julia_assign("drug_effectiveness", drug_effectiveness)

  humans = JuliaCall::julia_eval("administer_drug(humans, indices, drug_effectiveness)")
  return(humans)
}

