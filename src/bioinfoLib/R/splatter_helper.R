library(SingleCellExperiment)
library(SummarizedExperiment)
library(splatter)

splatSimulateLoop <- function(params = newSplatParams(), nSteps, sparsify = TRUE, verbose = TRUE, ...) {
  checkmate::assertClass(params, "SplatParams")

  if (verbose) {
    message("Getting parameters...")
  }
  params <- setParams(params, ...)
  params <- splatter:::expandParams(params)
  validObject(params)

  # Set random seed
  seed <- getParam(params, "seed")
  # Get the parameters we are going to use
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  nBatches <- getParam(params, "nBatches")
  batch.cells <- getParam(params, "batchCells")
  nGroups <- getParam(params, "nGroups")
  group.prob <- getParam(params, "group.prob")

  # Set up name vectors
  if (verbose) {
    message("Creating simulation object...")
  }
  cell.names <- paste0("Cell", seq_len(nCells))
  gene.names <- paste0("Gene", seq_len(nGenes))
  batch.names <- paste0("Batch", seq_len(nBatches))

  # Create SingleCellExperiment to store simulation
  cells <- data.frame(Cell = cell.names)
  rownames(cells) <- cell.names
  features <- data.frame(Gene = gene.names)
  rownames(features) <- gene.names
  sim <- SingleCellExperiment(
    rowData = features, colData = cells,
    metadata = list(Params = params)
  )

  # Make batches vector which is the index of param$batchCells repeated
  # params$batchCells[index] times
  batches <- lapply(seq_len(nBatches), function(i, b) {
    rep(i, b[i])
  }, b = batch.cells)
  batches <- unlist(batches)
  colData(sim)$Batch <- batch.names[batches]

  withr::with_seed(seed, {
    colData(sim)$Group <- factor(rep("loop", nCells))
    if (verbose) {
      message("Simulating library sizes...")
    }
    sim <- splatter:::splatSimLibSizes(sim, params)
    if (verbose) {
      message("Simulating gene means...")
    }
    sim <- splatter:::splatSimGeneMeans(sim, params)
    if (nBatches > 1) {
      if (verbose) {
        message("Simulating batch effects...")
      }
      sim <- splatter:::splatSimBatchEffects(sim, params)
    }
    sim <- splatter:::splatSimBatchCellMeans(sim, params)

    if (verbose) {
      message("Simulating path steps...")
    }
    sim <- splatSimLoopCellMean(sim, params, nSteps)

    if (verbose) {
      message("Simulating BCV...")
    }
    sim <- splatter:::splatSimBCVMeans(sim, params)
    if (verbose) {
      message("Simulating counts...")
    }
    sim <- splatter:::splatSimTrueCounts(sim, params)
    if (verbose) {
      message("Simulating dropout (if needed)...")
    }
    sim <- splatter:::splatSimDropout(sim, params)
  })

  if (sparsify) {
    if (verbose) {
      message("Sparsifying assays...")
    }
    assays(sim) <- splatter:::sparsifyMatrices(
      assays(sim),
      auto = TRUE,
      verbose = verbose
    )
  }

  if (verbose) {
    message("Done!")
  }
  return(sim)
}

splatSimLoopCellMean <- function(sim, params, nSteps) {
  nGenes <- getParam(params, "nGenes")
  nCells <- length(colData(sim)$Cell)
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  # path.nSteps <- getParam(params, "path.nSteps")
  path.nSteps <- nSteps
  path.nonlinearProb <- getParam(params, "path.nonlinearProb")
  path.sigmaFac <- getParam(params, "path.sigmaFac")
  path.from <- getParam(params, "path.from")
  exp.lib.sizes <- colData(sim)$ExpLibSize
  nGroups <- getParam(params, "nGroups")
  path.skew <- getParam(params, "path.skew")
  groups <- colData(sim)$Group
  batch.means.cell <- assays(sim)$BatchCellMeans

  # Generate non-linear path factors
  for (idx in seq_along(path.from)) {
    # Select genes to follow a non-linear path
    is.nonlinear <- as.logical(rbinom(nGenes, 1, path.nonlinearProb))
    sigma.facs <- rep(0, nGenes)
    sigma.facs[is.nonlinear] <- path.sigmaFac
    rowData(sim)[[paste0("SigmaFacPath", idx)]] <- sigma.facs
  }

  # Generate paths. Each path is a matrix with path.nSteps columns and
  # nGenes rows where the expression from each genes changes along the path.
  path.steps <- lapply(seq_along(path.from), function(idx) {
    from <- 0
    # Find the factors at the starting position
    if (from == 0) {
      facs.start <- rep(1, nGenes)
    } else {
      facs.start <- rowData(sim)[[paste0("DEFacPath", from)]]
    }
    # Find the factors at the end position
    # Return to start
    facs.end <- facs.start

    # Get the non-linear factors
    sigma.facs <- rowData(sim)[[paste0("SigmaFacPath", idx)]]
    print(path.nSteps)
    # Build Brownian bridges from start to end
    steps <- splatter:::buildBridges(
      facs.start,
      facs.end,
      n = path.nSteps[idx],
      sigma.fac = sigma.facs
    )

    return(t(steps))
  })

  # Randomly assign a position in the appropriate path to each cell
  path.probs <- lapply(seq_len(nGroups), function(idx) {
    probs <- seq(
      path.skew[idx], 1 - path.skew[idx],
      length = path.nSteps[idx]
    )
    probs <- probs / sum(probs)
    return(probs)
  })

  steps <- vapply(factor(groups), function(path) {
    step <- sample(seq_len(path.nSteps[path]), 1, prob = path.probs[[path]])
  }, c(Step = 0))

  # Collect the underlying expression levels for each cell
  cell.facs.gene <- lapply(seq_len(nCells), function(idx) {
    path <- factor(groups)[idx]
    step <- steps[idx]
    cell.means <- path.steps[[path]][, step]
  })
  cell.facs.gene <- do.call(cbind, cell.facs.gene)

  # Adjust expression based on library size
  cell.means.gene <- batch.means.cell * cell.facs.gene
  cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
  base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)

  colnames(base.means.cell) <- cell.names
  rownames(base.means.cell) <- gene.names

  colData(sim)$Step <- steps
  assays(sim)$BaseCellMeans <- base.means.cell

  return(sim)
}
