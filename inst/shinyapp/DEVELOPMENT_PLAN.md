# ZinaSuite Shiny Application Development Plan

## Project Overview

This document outlines the comprehensive development plan for ZinaSuite Shiny application, ensuring 100% feature parity with UCSCXenaShiny while implementing modern architecture patterns.

**Target**: Complete feature coverage of UCSCXenaShiny's 69 modules in ZinaSuite's modular architecture.

---

## Architecture Comparison

### UCSCXenaShiny Structure
- **69 module files** across 8 categories
- Traditional `callModule` pattern
- Mixed UI/business logic
- Limited testing coverage

### ZinaSuite Modern Architecture
- **15 modules** currently implemented
- Modern `moduleServer` pattern
- Separated business logic in `R/shiny-logic.R`
- 319 automated tests
- Async computing with `mirai`
- TTL caching system

---

## Development Phases

### Phase 1: Core Functionality (P0 Priority) 🎯

**Goal**: Establish foundation with essential analysis features

#### 1.1 Extract General Analysis Functions
- [ ] `analyze_scatter_correlation()` - Scatter correlation analysis
- [ ] `analyze_matrix_correlation()` - Matrix correlation analysis  
- [ ] `analyze_group_comparison()` - Group comparison analysis
- [ ] `analyze_survival()` - Survival analysis
- [ ] `analyze_dimension_reduction()` - Dimension reduction analysis

**Files to create/modify**:
- `R/analysis-general.R` (new)
- `tests/testthat/test-analysis-general.R` (new)

#### 1.2 Enhance TCGA Quick Analysis Module
- [ ] Tumor vs Normal comparison (mod_quick_tcga_TN)
- [ ] Anatomy visualization (mod_quick_tcga_anatomy)
- [ ] Molecular correlation (mod_quick_tcga_cor)
- [ ] TIL analysis (mod_quick_tcga_til)
- [ ] Immune features (mod_quick_tcga_immune)
- [ ] Index analysis (TMB/Stemness/MSI)
- [ ] Pathway analysis (mod_quick_tcga_pw)
- [ ] Mutation analysis (mod_quick_tcga_mut)
- [ ] KM survival (mod_quick_tcga_km)
- [ ] Cox regression (mod_quick_tcga_cox)
- [ ] Dimension reduction (mod_quick_tcga_dim)

#### 1.3 Data Query Module Enhancement
- [ ] Gene expression query
- [ ] miRNA query
- [ ] Methylation query
- [ ] CNV query
- [ ] Protein query
- [ ] Mutation query

#### 1.4 Pan-Cancer Analysis Module
- [ ] Cross-cancer expression distribution
- [ ] Cross-cancer comparison
- [ ] Cross-cancer correlation
- [ ] Cross-cancer survival

**Deliverables**:
- All P0 functions exported with documentation
- Unit tests for all functions
- Shiny modules integrated
- shinytest2 tests for user flows

---

### Phase 2: Extended Functionality (P1 Priority) 📊

**Goal**: Multi-datasource support and advanced analysis

#### 2.1 PCAWG Deep Analysis Pipeline
- [ ] PCAWG correlation (o2o, o2m, m2o modes)
- [ ] PCAWG comparison (o2o, o2m, m2o modes)
- [ ] PCAWG survival (o2o, o2m, m2o modes)

#### 2.2 CCLE Deep Analysis Pipeline
- [ ] CCLE correlation (o2o, m2o modes)
- [ ] CCLE comparison (o2o, m2o modes)

#### 2.3 Immune Analysis Enhancement
- [ ] Immune correlation analysis
- [ ] TIL correlation analysis
- [ ] MSI correlation analysis
- [ ] TMB correlation analysis

#### 2.4 Mutation Analysis Enhancement
- [ ] Mutation frequency analysis
- [ ] Mutation vs expression
- [ ] Mutation survival analysis
- [ ] Co-occurrence analysis

---

### Phase 3: Advanced Features (P2 Priority) 🔬

**Goal**: Professional analysis capabilities

#### 3.1 Pharmacogenomics Module
- [ ] Drug-omic correlation analysis
- [ ] Drug-mutation association
- [ ] Drug sensitivity profiling
- [ ] Feature across cell types

#### 3.2 Dimension Reduction Module
- [ ] PCA analysis
- [ ] t-SNE analysis
- [ ] UMAP analysis

#### 3.3 Batch Analysis Module
- [ ] Multi-gene correlation batch
- [ ] Multi-gene Cox batch
- [ ] Result aggregation

#### 3.4 Cross-Omics Analysis
- [ ] Gene cross-omics (o2m)
- [ ] Pathway cross-omics (o2m)

---

### Phase 4: Testing Coverage (P1 Priority) ✅

**Goal**: 100% core functionality test coverage

#### 4.1 Unit Tests
- [ ] All exported functions have tests
- [ ] Edge case coverage
- [ ] Error handling tests

#### 4.2 Shiny Tests
- [ ] Module UI rendering tests
- [ ] User interaction flow tests
- [ ] Async operation tests

#### 4.3 Integration Tests
- [ ] End-to-end workflows
- [ ] Data flow validation
- [ ] Performance benchmarks

---

### Phase 5: CRAN Standardization (P2 Priority) 📦

**Goal**: CRAN-ready package

#### 5.1 Documentation
- [ ] All functions have complete roxygen2 docs
- [ ] All examples are runnable
- [ ] Vignettes for major features

#### 5.2 Quality Checks
- [ ] R CMD check passes with 0 errors/warnings/notes
- [ ] Code coverage > 80%
- [ ] pkgdown site complete

#### 5.3 Release Preparation
- [ ] NEWS.md updated
- [ ] Version bumped
- [ ] CRAN submission ready

---

## Module Development Workflow

### For Each Module:

1. **Analysis** (15 min)
   - Review UCSCXenaShiny source
   - Identify business logic
   - Plan function signatures

2. **Function Extraction** (30 min)
   - Create pure R functions
   - Move to appropriate R/ file
   - Add roxygen2 documentation

3. **Unit Testing** (30 min)
   - Write testthat tests
   - Test happy path and edge cases
   - Run tests: `devtools::test()`

4. **Shiny Integration** (30 min)
   - Update/create module file
   - Call extracted functions
   - Handle reactive logic

5. **Shiny Testing** (30 min)
   - Write shinytest2 tests
   - Test user interactions
   - Run tests

6. **Documentation** (15 min)
   - Update examples
   - Run `devtools::document()`

7. **Validation** (15 min)
   - Run `devtools::check()`
   - Fix any issues

8. **Commit** (5 min)
   - `git add .`
   - `git commit -m "feat(module): description"`

**Total per module**: ~2.5 hours

---

## Coding Standards

### Function Design
```r
#' Brief description
#'
#' Detailed description of what the function does
#'
#' @param data Input data frame
#' @param gene Gene symbol
#' @param cancer Cancer type code
#' @return A ggplot object or analysis result
#' @export
#' @examples
#' \donttest{
#' result <- analyze_feature(data, "TP53", "BRCA")
#' }
analyze_feature <- function(data, gene, cancer) {
  # Implementation
}
```

### Module Structure
```r
# mod_feature.R
mod_feature_ui <- function(id) {
  ns <- NS(id)
  # UI definition
}

mod_feature_server <- function(id, data_source) {
  moduleServer(id, function(input, output, session) {
    # Reactive logic
    # Call extracted functions
  })
}
```

### Testing Pattern
```r
test_that("analyze_feature works correctly", {
  result <- analyze_feature(test_data, "TP53", "BRCA")
  expect_s3_class(result, "ggplot")
  expect_error(analyze_feature(NULL, "TP53", "BRCA"))
})
```

---

## Progress Tracking

See [FUNCTION_CHECKLIST.md](FUNCTION_CHECKLIST.md) for detailed progress tracking.

---

## Notes

- Each task should be completed in isolation
- Run tests after every change
- Commit after each completed feature
- Update checklist immediately after completion
- Ask for clarification when requirements are unclear
