#' @title Function \code{volcano}
#' @description Computes volcano plots by plotting probabilities of propositions vs the effect sizes
#'
#' @export
#' @return a volcano plot of probabilities of propositions vs the effect sizes
#' @param chain a \code{Chain} object
volcano = function(chain){
  e = effect_sizes(chain)
  p = probs(chain)
  me = melt(e)
  mp = melt(p)
  d = data.frame(Proposition = me$Var2, Effect = me$value, Probability = mp$value)
  ggplot(d) +
    stat_binhex(aes_string(x = "Effect", y = "Probability", fill="log(..count..)"), bins = 100) + 
      xlab("\nEstimated Effect Size") + 
      ylab("Probability\n") + 
      scale_fill_continuous(low = "darkGrey", high = "black") +
      facet_wrap(~Proposition, scales = "free_x")
}