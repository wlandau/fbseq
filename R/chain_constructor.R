#' @include samplers.R special_beta_priors.R generate_starts.R
NULL

#' @title Function \code{Chain}
#' @description Construct a \code{Chain} object, the 
#' one argument to \code{single_chain(...)}.
#' @details If \code{list} is set, all other arguments will be dropped.
#' @export
#' @return a \code{Chain} object ready to be passed to \code{single_chain()}
#'
#' @param scenario a \code{Scenario} object with count data, the design matrix, etc.
#' @param configs A \code{Configs} object of MCMC control parameters.
#' @param starts A \code{Starts} object of model parameter starting bounds.
#' @param slots a list of slots to be plugged into a \code{Chain} object
Chain = function(scenario, configs = Configs(), starts = Starts(), slots = NULL){
  chain = new("Chain")

  if(!is.null(slots) & class(slots) == "list" & !is.data.frame(slots)){
    for(n in intersect(names(slots), slotNames(chain)))
      slot(chain, n) = as(slots[[n]], class(slot(chain, n)))
    return(chain)
  }

  starts = generate_starts(scenario@counts, scenario@design, starts)
  chain = plug_in_chain(chain, scenario, configs = configs, starts = starts)
  fill_easy_gaps(chain, scenario)
}

#' @title Function \code{plug_in_chain}
#' @description Plug the slots of \code{Configs} and \code{Starts}
#' objects into the analogous slots of a \code{Chain} object. Used
#' internally
#' @export
#' @return a \code{Chain} object 
#'
#' @param chain a \code{Chain} object
#' @param scenario an \code{Scenario} object
#' @param configs A \code{Configs} object of MCMC control parameters.
#' @param starts A \code{Starts} object of model parameter starting bounds.
plug_in_chain = function(chain, scenario, configs, starts){
  chain@iterations = as.integer(configs@iterations)
  stopifnot(configs@thin > 0)
  if(any(configs@priors %in% special_beta_priors())) stopifnot("xi" %in% configs@parameter_sets_update)

  subtract = c("parameter_sets_return", "parameter_sets_update", "priors")
  for(n in setdiff(intersect(slotNames(chain), slotNames(configs)), subtract))
    slot(chain, n) = as(slot(configs, n), class(slot(chain, n)))

  for(n in c("parameter_sets_return", "parameter_sets_update")){
    slot(chain, n) = as.integer(parameters() %in% slot(configs, n))
    names(slot(chain, n)) = parameters()
  }

  pvec = 0:length(special_beta_priors())
  names(pvec) = c("normal", special_beta_priors())
  chain@priors = pvec[configs@priors]
  if(length(chain@priors) == 1) chain@priors = rep(chain@priors, ncol(scenario@design))
  stopifnot(length(chain@priors) == ncol(scenario@design))

  chain = assign_samplers(chain, configs)

  for(n in slotNames(chain)){
    if(grepl("Start", n))
      slot(chain, n) = as(slot(starts, gsub("Start", "", n)), class(slot(chain, n)))
    else if(n %in% slotNames(starts) && !(paste0(n, "Start") %in% slotNames(chain)))
      slot(chain, n) = as(slot(starts, n), class(slot(chain, n)))
  }

  chain@bounds = scenario@bounds
  chain@contrasts = unlist(scenario@contrasts)
  chain@counts = as.integer(scenario@counts)
  chain@design = as.numeric(scenario@design)
  chain@propositions = as.integer(unlist(lapply(scenario@propositions, function(x){1:length(scenario@contrasts) %in% x})))
  chain@supplement = scenario@supplement

  chain
}

#' @title Function \code{fill_easy_gaps}
#' @description Set required slots in a \code{Chain} object.
#' Used internally.
#' @export
#' @return a \code{Chain} object 
#'
#' @param chain a \code{Chain} object
#' @param scenario an \code{Scenario} object
#' Must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' sets of gene-specific variables.
fill_easy_gaps = function(chain, scenario){
  designUnique = apply(scenario@design, 2, function(x){
    out = sort(unique(x[x != 0]))
    c(out, rep(0, dim(scenario@design)[1] - length(out)))
  })

  if(!length(chain@effects_update_beta)) chain@effects_update_beta = 1:ncol(scenario@design)
  if(!length(chain@genes_return)) chain@genes_return = sample.int(nrow(scenario@counts), 3, replace = T)
  if(!length(chain@genes_return_epsilon)) chain@genes_return_epsilon = sample.int(nrow(scenario@counts), 3, replace = T)
  if(!length(chain@libraries_return)) chain@libraries_return = sample.int(ncol(scenario@counts), 3, replace = T)
  if(!length(chain@libraries_return_epsilon)) chain@libraries_return_epsilon = sample.int(ncol(scenario@counts), 3, replace = T)

  chain@bound_names = as.character(names(scenario@bounds))
  chain@contrast_names = as.character(names(scenario@contrasts))
  chain@gene_names = as.character(rownames(scenario@counts))
  chain@library_names = as.character(colnames(scenario@counts))
  chain@proposition_names = as.character(names(scenario@propositions))

  chain@C = length(scenario@contrasts)
  chain@counts = as.integer(scenario@counts)
  chain@countSums_g = as.integer(apply(scenario@counts, 1, sum))
  chain@countSums_n = as.integer(apply(scenario@counts, 2, sum))
  chain@designUnique = as.numeric(designUnique)
  chain@designUniqueN = as.integer(apply(scenario@design, 2, function(x){length(unique(x[x != 0]))}))
  chain@genes_return = sort(chain@genes_return)
  chain@genes_return_epsilon = sort(chain@genes_return_epsilon)
  chain@G = G = nrow(scenario@counts)
  chain@Greturn = Greturn = length(chain@genes_return)
  chain@GreturnEpsilon = GreturnEpsilon = length(chain@genes_return_epsilon)
  chain@libraries_return = sort(chain@libraries_return)
  chain@libraries_return_epsilon = sort(chain@libraries_return_epsilon)
  chain@L = L = ncol(scenario@design)
  chain@Lupdate_beta = length(chain@effects_update_beta)
  chain@N = N = nrow(scenario@design)
  chain@Nreturn = Nreturn = length(chain@libraries_return)
  chain@NreturnEpsilon = NreturnEpsilon = length(chain@libraries_return_epsilon)
  chain@P = length(scenario@propositions)
  chain@probs = rep(0, chain@P * chain@G)

  lengths = c(
    beta = L*Greturn,
    epsilon = NreturnEpsilon*GreturnEpsilon,
    gamma = Greturn,
    nu = 1,
    sigmaSquared = L,
    tau = 1,
    theta = L,
    xi = L*Greturn)

  for(s in names(lengths))
    if(as.logical(chain@parameter_sets_return[s]))
      slot(chain, s) = rep(0, chain@iterations*lengths[s])

  lengths = c(
    beta = L*G,
    epsilon = N*G,
    gamma = G,
    nu = 1,
    sigmaSquared = L,
    tau = 1,
    theta = L,
    xi = L*G)

  for(suffix in c("PostMean", "PostMeanSquare", "Tune", "TuneAux"))
    for(s in names(lengths)){
      n = paste0(s, suffix)
      slot(chain, n) = rep(0, lengths[s])
    }

  chain
}
