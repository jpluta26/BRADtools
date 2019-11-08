library(ggplot2)

# pluta 6/7/16


x <- seq(1:5)
y <- c(2,1,2,1,3)

cor.test(x,y)

cor.test(c(x,x,x,x), c(y,y,y,y))




# get variance of a binomial variable
getBinomVar <- function(n, p)
{
  return( n * p * (1-p) )
}

# function to plot comparison of two distributions
plotDstr <- function(x1, x2)
{
  dat1 <- data.frame(x=density(x1)$x, y=density(x1)$y)
  dat2 <- data.frame(x=density(x2)$x, y=density(x2)$y)
  p <- ggplot() + geom_polygon(data=dat1, aes(x=x, y=y), alpha=0.3, color="black", fill="blue") +
    geom_polygon(dat=dat2, aes(x=x, y=y), alpha=0.5, color="black", fill="red")  
  return(p)
  
}

# p-value examples
x1 <- rnorm(100, mean=0, sd=1)
x2 <- rnorm(100, mean=1, sd=1)
r <- qnorm(0.025, 0, 1)

dat1 <- data.frame(x=density(x1)$x, y=density(x1)$y)
dat2 <- data.frame(x=density(x2)$x, y=density(x2)$y)

# plot to add rejection region
p <- ggplot() + geom_polygon(data=dat1, aes(x=x, y=y), alpha=0.3, color="black", fill="blue") +
  geom_vline(xintercept=qnorm(0.025, 0, 1), linetype="dashed") + 
  geom_vline(xintercept=-qnorm(0.025,0 ,1), linetype="dashed")
    
print(p)


t <- t.test(x1, x2, paired=FALSE, alternative = "two.sided")

p <- plotDstr(x1,x2)
print(p)



x1 <- rnorm(1000, mean=0, sd=1)
x2 <- rnorm(1000, mean=1, sd=1)

p <- plotDstr(x1,x2)
print(p)

t <- t.test(x1, x2, paired=FALSE, alternative = "two.sided")



# ----------------------------------------------------------------------------------- #
# first example: frequentist test of a fair coin
# we want to test if a coin is fair.
n = 100
p = 0.5


# suppose the coin is indeed fair, and we flipped it 100 times. this is the distribution we would expect to see:
dat <- data.frame(val   = c(dbinom(seq(1:100), n, p) ),
                  x = seq(1:100) )

# base plot of 2 distributions
p1 <- ggplot() + geom_polygon(data=dat, aes(x=x, y=val), alpha=0.3, color="black", fill="blue") +
      geom_vline(xintercept = 50, linetype="dashed", color="black") + xlab("Number of Heads") + ylab("Probability")
        
print(p1)




# lognormal distribution with filled in rejection region
dat <- data.frame(val = c(dlnorm(seq(from=0, to=25, by=.1), meanlog=1, sdlog=.5) ), x=seq(from=0, to=25, by=.1))


p1 <- ggplot() + geom_polygon(data=dat, aes(x=x, y=val), fill="blue", alpha=0.3, color="black") + xlim(0,12) +
  geom_vline(xintercept=qlnorm(0.95, meanlog=1, sdlog=.5), linetype="dashed", color="black") +
  geom_ribbon(data=dat[dat$x>qlnorm(0.95, meanlog=1, sdlog=.5),], aes(x, ymin=0, ymax=val), fill="red", alpha=0.5)
print(p1)

library(BMS)


# eg, we would expect to see around 50 heads, but anything within 35-65 would be plausible
hi <- n*p + sqrt(getBinomVar(n,p)) * 3  
lo <- n*p - sqrt(getBinomVar(n,p)) * 3 



# now lets say we gather some data. we flip the coin 100 times and record the number of heads, and repeat this
# experiment 100 times
set.seed(1234)
out.vals <- rbinom(100, 100, .65)
out <- density(out.vals)

dat.sim <- data.frame(x=out$x, y=out$y)
p2 <- ggplot() + geom_polygon(data=dat, aes(x=x, y=val), alpha=0.3, color="black", fill="blue") +
  geom_polygon(dat=dat.sim, aes(x=x, y=y), alpha=0.5, color="black", fill="red")
print(p2)

w <- chisq.test(c(65.14,34.86), p=c(.5,.5))




# bayesian
# p(D|H0) is the probability of obtaining a particular test statistic, where a p-value is the 
# probability of obtaining a test statistic greater than observed
# depends on: probability hypothesis is true P(H0)
# probability alternative is true P(H1)
# and the probability that the data would have been observed if the alternative is true, P(D|H1)
# these are all unknown-> use priors
# criticism of bayes, often very little information to inform the priors

# a bayesian hypothesis test compares hte evidence of H0 to Ha; frequentists can only reject or fail to reject H0
# in other words, a prior isn't specified for Ha.
# at what point is a BF considered evidence that H0 is true? its arbitrary
# prior distribution influneces BF


# bayesian part:
# both people believe the coin is unfair. P1 believes P(H)=0.6, P2 belives P(H)=0.7. these are called "priors".
# then we get actual data and see how it cmopares to our priors.


# plot priors
n = 100
x <- seq(1,100,by=1)
dat <- data.frame(grp = factor(   rep(  c("P1", "P2"), each=n)   ),
                  val   = c(dbinom(x, n, 0.6), dbinom(x, n, 0.7) ),
                  x = rep(seq(1:100),2) )

p1 <- ggplot() + geom_polygon(data=dat, aes(x=x, y=val, fill=grp), alpha=0.3, color="black") +
  scale_fill_manual(name="Title", values=c("red", "blue")) 
print(p1)

# suppose we actually observe 62 heads out of 100 flips. whats the evidence here?
P1 <- dbinom(62, 100, .6) # the probability of observing 62 heads under P1 assumptions
P2 <- dbinom(62, 100, .7) # the probability of observing 62 heads under P2 assumptions

# bayes factor, the ratio of the evidence. implies that the observed data favors P1's assumptions by a factor of 4.
BF <- dbinom(62, 100, .6)/dbinom(62,100,.7)

# but we cant do the same for a frequentist test. we only have a prior on H0, eg, p(heads) = 0.5




# exampe stuff here
# example of plotting true distributions
# 100 trials of flipping a coin 100 times
n = 100
x <- seq(1,100,by=1)
dat <- data.frame(grp = factor(   rep(  c("A", "B"), each=n)   ),
                  val   = c(dbinom(x, n, 0.5), dbinom(x, n, 0.62) ),
                  x = rep(seq(1:100),2) )

# base plot of 2 distributions
p1 <- ggplot() + geom_polygon(data=dat, aes(x=x, y=val, fill=grp), alpha=0.3, color="black") +
      scale_fill_manual(name="Title", values=c("red", "blue")) 
print(p1)
  

# plot with simulation data overlaid on theoretical distributions
# rbinom(# of trials, # of flips per trial, prob of success per flip)
# output is number of successes per trial
set.seed(1234)
out <- rbinom(100, 100, .5)

dat.sim <- data.frame(x=out$x, y=out$y)
p2 <- ggplot() + geom_polygon(data=dat, aes(x=x, y=val, fill=grp), alpha=0.3, color="black") +
  scale_fill_manual(name="Title", values=c("red", "blue")) + geom_polygon(dat=dat.sim, aes(x=x, y=y), alpha=0.5)
print(p2)



set.seed(1234)
out <- density(rbinom(10000, 100, .5))

dat.sim <- data.frame(x=out$x, y=out$y)
p2 <- ggplot() + geom_polygon(data=dat, aes(x=x, y=val, fill=grp), alpha=0.3, color="black") +
  scale_fill_manual(name="Title", values=c("red", "blue")) + geom_polygon(dat=dat.sim, aes(x=x, y=y), alpha=0.5)
print(p2)



# H0: coin is fair -> P(H) = 0.5
# HA: coin is not fair P(H) != 0.5
w <- chisq.test(c(62,38), p=c(.5,.5))

# cant reject null; not enough evidence to suggest coin isn't fair

