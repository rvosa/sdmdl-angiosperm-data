forward.model.selection <- function(tree,data,dependent,model.template = NULL,predictors = NULL) {
    
    # create the template in the first pass, initialize once
    if ( is.null(model.template) ) {
        model.template <- paste(dependent, '~ 1 + %s')
    }
    
    # candidate predictors, initialize once
    if ( is.null(predictors) ) {
        predictors <- names(data)
        predictors <- predictors[!predictors %in% c(dependent)]
    }
    
    # named vector of AICs in current iteration
    aics <- vector( mode = "numeric", length = length(predictors))
    names(aics) <- predictors
    
    # iterate over predictors
    for ( i in 1:length(predictors) ) {
        myModel <- sprintf(model.template, predictors[i])
        result <- phyloglm(myModel, data, tree, method = "logistic_MPLE", btol = 100)
        aics[i] <- result$aic
    }
    
    # pick the best result
    aics <- sort(aics)
    best.predictor <- names(aics)[1]
    predictors <- predictors[!predictors %in% c(best.predictor)]
    best.aic <- aics[1]
    best.model <- sprintf(model.template, best.predictor)
    message(sprintf('ADD: model="%s", aic="%f"',best.model,best.aic))
    
    # now try to prune
    parsed.model <- as.formula(best.model)
    predictors.grown <- labels(terms(parsed.model))
    aics.pruned <- vector( mode = "numeric", length = length(predictors.grown))
    for ( j in 1:length(predictors.grown) ) {
        predictors.pruned <- predictors.grown[-j]
        myModel.pruned <- sprintf(
            '%s ~ 1 + %s', 
            dependent, 
            paste(predictors.pruned, collapse = '+')
        )
        message(predictors)
        result.pruned <- phyloglm(myModel.pruned, data, tree, method = "logistic_MPLE", btol = 100)
        aics.pruned[i] <- result.pruned$aic
    }
    aics.pruned <- sort(aics.pruned)
    best.pruned.aic <- aics.pruned[1]
    if ( best.pruned.aic < best.aic ) {
        pruned.predictor <- names(aics.pruned)[1]
        predictors <- predictors.grown[!predictors.grown %in% c(pruned.predictor)]
        best.model <- sprintf(
            '%s ~ 1 + %s', 
            dependent, 
            paste(predictors, collapse = '+')
        )
        message(sprintf('PRUNED: model="%s", aic="%f"',best.model,best.pruned.aic))
    }
    
    # expand the model template, shrink predictors, recurse further
    model.template <- paste(best.model, '+ %s')
    forward.model.selection(tree,data,dependent,model.template,predictors)
}