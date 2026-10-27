#' data.simulation
#' 
#' This function simulates data. 
#' 
#' @param N.sim Number of simulations.
#' @param N.max   Maximum number of patients to enroll.
#' @param mu.trt Log-location parameters (and therefore log-medians) for
#'   recurrence and death in the treatment arm.
#' @param Sigma.trt Variance-covariance matrix of log event times in the treatment arm.
#' @param mu.ctrl Log-location parameters (and therefore log-medians) for
#'   recurrence and death in the control arm.
#' @param Sigma.ctrl Variance-covariance matrix of log event times in the control arm.
#' @param cens_upper Upper limit of the Uniform(0, cens_upper) censoring distribution.
#' @returns An array with N.max rows, 7 columns and N.sim slices. The columns are:\tabular{ll}{
#'    \code{recurrence_t} \tab Time to recurrence.  \cr
#'    \tab \cr
#'    \code{death_t} \tab Time to death. \cr
#'    \tab \cr
#'    \code{censor_t} \tab Time to censor. \cr
#'    \tab \cr
#'    \code{group} \tab Treatment assignment.  \cr
#'    \tab \cr
#'    \code{id} \tab Patient id, numeric starting from 1. \cr 
#'    \tab \cr
#'    \code{delta1} \tab Placeholder for the censoring indicator for recurrence, leave blank. \cr
#'    \tab \cr
#'    \code{delta2} \tab Placeholder for the censoring indicator for death, leave blank. \cr
#' }
#' @examples
#' \dontrun{
#' data.simulation(N.sim = 1000, N.max = 20,
#'                mu.trt = c(0.2,0.3), Sigma.trt = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T),
#'                mu.ctrl = c(0.2,0.3),Sigma.ctrl = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T),
#'                cens_upper = 5)
#' }

#' @export
data.simulation<- function(N.sim,N.max,mu.trt,Sigma.trt,
                           mu.ctrl,Sigma.ctrl,cens_upper
){
  stopifnot(length(mu.trt) == 2L, length(mu.ctrl) == 2L)
  stopifnot(all(dim(Sigma.trt) == c(2L, 2L)), all(dim(Sigma.ctrl) == c(2L, 2L)))
  stopifnot(N.sim > 0L, N.max > 0L, cens_upper > 0)

  myData<-array(dim=c(N.max,7,N.sim));#三维数组[N.max 行, 7 列, N.sim 层]
  for (j in 1:N.sim){
    
    # Simulate randomization groups
    rand.list = sample(c(0,1), N.max, replace = T)
    
    # Simulate correlated latent times to recurrence and death. If death occurs
    # first, the later latent recurrence is not observable; analysis scripts
    # handle this through the derived first-composite-event time.
    trt.events <- 	as.data.frame(MASS::mvrnorm(sum(rand.list), mu= mu.trt,Sigma= Sigma.trt))    
    trt.events$id =  which(rand.list==1); 
    trt.events$grp = 1; 

    ctrl.events <- as.data.frame(MASS::mvrnorm(sum(rand.list==0),
     mu= mu.ctrl,
     Sigma= Sigma.ctrl))     

    ctrl.events$id = which(rand.list==0); 
    ctrl.events$grp = 0

    # simulate time to censor
    time.censor = stats::runif(N.max, 0, cens_upper); 
    eventsall = rbind(trt.events, ctrl.events); 
    eventsall =  eventsall[order(eventsall$id),]
    myData[,1,j]<- exp(eventsall[,1])                 # time to recurrence 
    myData[,2,j]<- exp(eventsall[,2])                 # time to death
    myData[,3,j]<- time.censor                        # censoring time
    myData[,4,j]<- eventsall[,4]                      # treatment group assignments
    myData[,5,j]<-  eventsall[,3]                     # id

  }
  
  # Create dimension names
  dim_names <- list(
    NULL,
    Variables = c("recurrence_t", "death_t", "censor_t", "group", "id", "delta1", "delta2"),
    Simulations = paste("Simulation", 1:N.sim)
  )
  
  # Assign the dimension names to the array
  dimnames(myData) <- dim_names
  print("data成功生成")
  return(myData)
}
