# Purpose: Take the gammainc source file from https://rdrr.io/rforge/pracma/src/R/gammainc.R
#           and rewrite it using the autodiff paradigm
#
#
#
#-----------------------------------------------------------------------------------


gammainc <- function(x) {
  # I found this function online from pracma package
  # https://rdrr.io/rforge/pracma/src/R/gammainc.R
  if (!is.numeric(x[2]) || !is.numeric(x[1]))
    stop("All arguments must be real numbers.")
  if (length(x[2]) > 1 || length(x[1]) > 1)
    stop("Arguments must be of length 1; function is not vectorized.")
  if (x[2] < 0)
    stop("Argument 'x[2]' must be real and nonnegative.")
  if (x[1] == 0 && x[2] == 0)
    # return(c(lowinc = 0.0, uppinc = Inf, reginc = 0.0))
    return( uppinc = Inf)
  if (x[1] == 0)
    # return(c(lowinc = 0.0, uppinc = gamma(x[2]), reginc = 0.0))
    return(uppinc = Re(gim))

  if (x[1] > 0)  xam <- -x[1] + x[2]*log(x[1])
  else        xam <- -x[1] + x[2]*log(x[1] + 0i)
  if (abs(xam) > 700.0 || abs(x[2]) > 170.0) {
    warning("Arguments 'x[1]' and/or 'x[2]' are too large.")
    return(NA)
  }

  # Computation of the incomplete gamma function
  gin <- gim <- gip <- 0

  if (x[1] == 0.0) {
    ga <- gamma(x[2])
    gim <- ga
    gip <- 0.0
  } else if (x[1] <= 1.0 + x[2]) {
    s <- 1/x[2]
    r <- s
    for  (k in 1:60) {
      r <- r * x[1]/(x[2]+k);
      s <- s+r;
      if (abs(r/s) < 1e-15) break
    }
    gin <- exp(xam) * s
    ga <- gamma(x[2])
    gip <- gin/ga
    gim <- ga - gin
  } else if (x[1] > 1.0 + x[2]) {
    t0 <- 0
    for  (k in 60:1) {
      t0 <- (k-x[2])/(1 + k/(x[1]+t0))
    }
    gim <- exp(xam)/(x[1]+t0)
    ga <- gamma(x[2])
    gin <- ga - gim
    gip <- 1 - gim/ga
  }
  # return(c(lowinc = Re(gin), uppinc = Re(gim), reginc = Re(gip)))
  return(uppinc = Re(gim))
}

# gammainc file with autodiff paradigm and one if-statment skipped
NewGammainc_d <- function(x, a, x_d, a_d) {
  #rewrite the gammainc function using autodiff paradigm and skip on if statement.
  vneg1 = x
  vneg1_d = x_d

  v0 = a
  v0_d = a_d


  if (!is.numeric(v0) || !is.numeric(vneg1))
    stop("All arguments must be real numbers.")
  if (length(v0) > 1 || length(vneg1) > 1)
    stop("Arguments must be of length 1; function is not vectorized.")
  if (v0 < 0)
    stop("Argument 'v0' must be real and nonnegative.")
  if (vneg1 == 0 && v0 == 0){
    uppinc = Inf
    uppinc_d = Inf
    return(c(uppinc, uppinc_d))
  }
  if (vneg1 == 0){
    uppinc   = gamma(v0)
    uppinc_d = gamma(v0)*digamma(v0)*v0_d
    return(c(uppinc, uppinc_d))
  }


  if (vneg1 > 0){
    v1    = v0*log(vneg1)
    v1_d  = v0_d*log(vneg1) + v0*(1/vneg1)*vneg1_d
    xam   = -vneg1 + v1
    xam_d = -vneg1_d + v1_d
  }else{
    v2    = v0*log(vneg1 + 0i)
    v2_d  = v0_d*log(vneg1+0i) +v0*(1/(vneg1+0i) )*vneg1_d
    xam   = -vneg1 + v2
    xam_d = -vneg1_d + v2_d
  }

  if (abs(xam) > 700.0 || abs(v0) > 170.0) {
    warning("Arguments 'vneg1' and/or 'v0' are too large.")
    return(NA)
  }

  # Computation of the incomplete gamma function
  gin <- gim <- gip <- 0
  gin_d <- gim_d <- gip_d <- 0

  if (vneg1 == 0.0) {
    ga <- gamma(v0)
    ga_d <- gamma(v0)*digamma(v0)*v0_d

    gim <- ga
    gim_d = ga_d

    gip <- 0.0
    gip_d = 0

  }
  # else if
  # (vneg1 <= 1.0 + v0) {
  #
  #
  #
  #   ###################
  #   s <- 1/v0
  #   s_d = -s^2*v0_d
  #
  #   r <- s
  #   r_d = s_d
  #
  #   for  (k in 1:60) {
  #     v3   = v0 + k
  #     v3_d = v0_d
  #
  #     v4   = vneg1/v3
  #     v4_d = (vneg1_d*v3 - vneg1*v3_d)/(v3^2)
  #
  #     r    = r * v4;
  #     r_d  = r_d*v4 + r*v4_d;
  #
  #     s   = s + r;
  #     s_d = s_d + r_d;
  #
  #     if (abs(r/s) < 1e-15) break
  #   }
  #   v5   = exp(xam)
  #   v5_d = v5*xam_d
  #
  #   gin   = v5 * s
  #   gin_d = v5_d*s + v5*s_d
  #
  #   ga <- gamma(v0)
  #   ga_d = gamma(v0)*digamma(v0)*v0_d
  #
  #   gim   = ga - gin
  #   gim_d = ga_d - gin_d
  #   ###################
  # }
  # else if (vneg1 > 1.0 + v0) {
    t0 <- 0
    t0_d = 0
    for  (k in 60:1) {
      v6   = k-v0
      v6_d = -v0_d
      v7   = vneg1 + t0
      v7_d = vneg1_d + t0_d
      v8   = (1 + k/v7)
      v8_d = -k/v7^2*v7_d

      t0 <- v6/v8
      t0_d = (v6_d*v8 - v8_d*v6)/v8^2
    }
    v9 = exp(xam)
    v9_d = v9*xam_d

    v10 = (vneg1+t0)
    v10_d = vneg1_d + t0_d

    gim <- v9/v10
    gim_d = (v9_d*v10 - v10_d*v9)/v10^2

  #}

  uppinc = Re(gim)
  uppinc_d = Re(gim_d)
  return(c(uppinc, uppinc_d))
}

# gammainc file with autodiff paradigm
NewGammainc_d_2 <- function(x, a, x_d, a_d) {
  #rewrite the gammainc function using autodiff paradigm and dont skip the if statement.
  vneg1 = x
  vneg1_d = x_d

  v0 = a
  v0_d = a_d


  if (!is.numeric(v0) || !is.numeric(vneg1))
    stop("All arguments must be real numbers.")
  if (length(v0) > 1 || length(vneg1) > 1)
    stop("Arguments must be of length 1; function is not vectorized.")
  if (v0 < 0)
    stop("Argument 'v0' must be real and nonnegative.")
  if (vneg1 == 0 && v0 == 0){
    uppinc = Inf
    uppinc_d = Inf
    return(c(uppinc, uppinc_d))
  }
  if (vneg1 == 0){
    uppinc   = gamma(v0)
    uppinc_d = gamma(v0)*digamma(v0)*v0_d
    return(c(uppinc, uppinc_d))
  }


  if (vneg1 > 0){
    v1    = v0*log(vneg1)
    v1_d  = v0_d*log(vneg1) + v0*(1/vneg1)*vneg1_d
    xam   = -vneg1 + v1
    xam_d = -vneg1_d + v1_d
  }else{
    v2    = v0*log(vneg1 + 0i)
    v2_d  = v0_d*log(vneg1+0i) +v0*(1/(vneg1+0i) )*vneg1_d
    xam   = -vneg1 + v2
    xam_d = -vneg1_d + v2_d
  }

  if (abs(xam) > 700.0 || abs(v0) > 170.0) {
    warning("Arguments 'vneg1' and/or 'v0' are too large.")
    return(NA)
  }

  # Computation of the incomplete gamma function
  gin <- gim <- gip <- 0
  gin_d <- gim_d <- gip_d <- 0

  if (vneg1 == 0.0) {
    ga <- gamma(v0)
    ga_d <- gamma(v0)*digamma(v0)*v0_d

    gim <- ga
    gim_d = ga_d

    gip <- 0.0
    gip_d = 0

  }else if(vneg1 <= 1.0 + v0) {

    ###################
    s <- 1/v0
    s_d = -s^2*v0_d

    r <- s
    r_d = s_d

    for  (k in 1:60) {
      v3   = v0 + k
      v3_d = v0_d

      v4   = vneg1/v3
      v4_d = (vneg1_d*v3 - vneg1*v3_d)/(v3^2)

      r    = r * v4;
      r_d  = r_d*v4 + r*v4_d;

      s   = s + r;
      s_d = s_d + r_d;

      if (abs(r/s) < 1e-15) break
    }
    v5   = exp(xam)
    v5_d = v5*xam_d

    gin   = v5 * s
    gin_d = v5_d*s + v5*s_d

    ga <- gamma(v0)
    ga_d = gamma(v0)*digamma(v0)*v0_d

    gim   = ga - gin
    gim_d = ga_d - gin_d
    ###################
  }
  else if (vneg1 > 1.0 + v0) {
  t0 <- 0
  t0_d = 0
  for  (k in 60:1) {
    v6   = k-v0
    v6_d = -v0_d
    v7   = vneg1 + t0
    v7_d = vneg1_d + t0_d
    v8   = (1 + k/v7)
    v8_d = -k/v7^2*v7_d

    t0 <- v6/v8
    t0_d = (v6_d*v8 - v8_d*v6)/v8^2
  }
  v9 = exp(xam)
  v9_d = v9*xam_d

  v10 = (vneg1+t0)
  v10_d = vneg1_d + t0_d

  gim <- v9/v10
  gim_d = (v9_d*v10 - v10_d*v9)/v10^2

  }

  uppinc = Re(gim)
  uppinc_d = Re(gim_d)
  return(c(uppinc, uppinc_d))
}


# write the poisson cdf using the autodiff paradigm
myppois = function(x,a,x_d,a_d){

  All  = NewGammainc_d(a, x+1, a_d, x_d)
  G    = All[1]
  G_d  = All[2]

  v1   = gamma(x+1)
  v1_d = gamma(x+1)*digamma(x+1)*x_d

  z = G/v1
  z_d = (G_d*v1 - G*v1_d)/v1^2

  return(c(z,z_d))

}

