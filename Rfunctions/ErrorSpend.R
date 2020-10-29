# rho-family spending functions (Kim-DeMets) for alpha and beta 
ErrorSpend <- function(I,  #Information at analysis to be evaluated
              rho,  #parameter for the spending function
              beta_or_alpha,  #total alpha or beta to be spent
              Imax){  #Maximum information needed at end of trial
  beta_or_alpha*min(1,(I/Imax)^rho)
}
