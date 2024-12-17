
getORfromLogOR <- function( logOR, se )
{
  OR <- exp(logOR)
  OR.se <- OR * exp(se)
  return(list(OR,OR.se))
}

convert_logOR_to_OR <- function(logOR, SE_logOR) {
  OR <- exp(logOR)
  SE_OR <- OR * SE_logOR
  return(data.frame(OR = OR, SE_OR = SE_OR))
}
