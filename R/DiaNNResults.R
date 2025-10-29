# 3. Create results class
DiaNNResults <- R6Class(
  public = list(
    precursors = NULL,
    proteins = NULL,
    genes = NULL,
    metadata = NULL,
    
    initialize = function(processor) { ... },
    summary = function() { ... },
    plot = function() { ... },
    export = function() { ... }
  )
)