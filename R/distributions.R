# univariate truncated normal using rejection sampling
rtruncnorm <- function(n, mu = 0, sigma = 1, lb = -Inf, ub = +Inf)
{
    vec <- rnorm(n, mu, sigma)
    bad_indx <- which((vec < lb) | (vec > ub))
    if (length(bad_indx) > 0) {
        vec <- vec[-bad_indx]
        res_ratio <- n / max(length(vec),1)
        nsim <- ceiling(res_ratio * (n - length(vec) + 20))
        while (length(vec) < n) {
            tmp <- rnorm(nsim, mu, sigma)
            good_indx <- which(!((tmp < lb) | (tmp > ub)))
            vec <- c(vec,tmp[good_indx])
        }
    }
    return(vec[1:n])
}