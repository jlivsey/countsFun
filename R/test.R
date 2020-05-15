      # compute unnormalized weights
      wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)

      # break if I got NA weight
      if (any(is.na(wgh[t,]))| sum(wgh[t,])<10^(-8) ){
        nloglik = 10^8
        break
      }

      # normalized weights
      wghn = wgh[t,]/sum(wgh[t,])

      # Resampling: sample indices from multinomial distribution-see Step 4 of SISR in paper
      ESS = 1/sum(wghn^2)


      if(ESS<epsilon*N){
        old_state1 <- get_rand_state()
        ind = rmultinom(1,N,wghn)
        # sample particles
        znew = rep(znew,ind)
        zhat = tht*(znew-zhat)/rt0
        set_rand_state(old_state1)
      }