# ZinaSuite <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start --
[![R-CMD-check](https://github.com/WangLabCSU/ZinaSuite/workflows/R-CMD-check/badge.svg)](https://github.com/WangLabCSU/ZinaSuite/actions)
[![Codecov test coverage](https://codecov.io/gh/WangLabCSU/ZinaSuite/branch/main/graph/badge.svg)](https://app.codecov.io/gh/WangLabCSU/ZinaSuite?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/ZinaSuite)](https://CRAN.R-project.org/package=ZinaSuite)
<!-- badges: end -->

## Overview

**ZinaSuite** is a modernized R package for interactively exploring UCSC Xena datasets. It provides a comprehensive framework for downloading, analyzing, and visualizing cancer genomics data with enhanced performance through mirai-based asynchronous computing and modular Shiny application design.

### Key Features

- **🔬 Multi-Dataset Support**: TCGA, PCAWG, and CCLE datasets
- **⚡ High Performance**: Asynchronous computing with mirai for parallel processing
- **📊 Comprehensive Analysis**: Correlation, survival, mutation, immune, and dimension analysis
- **🎨 Publication-Ready Visualizations**: ggplot2-based plots with customizable themes
- **🖥️ Interactive Shiny App**: Modular design with 11 functional modules
- **💾 Smart Caching**: Memory and disk caching to reduce redundant queries
- **📑 Report Generation**: Export results to HTML, PDF, Word, Excel
- **✅ Data Validation**: Built-in validation for gene symbols, cancer types, and sample IDs

## Installation

```r
# Install from GitHub
remotes::install_github("WangLabCSU/ZinaSuite")

# Or install from CRAN (when available)
install.packages("ZinaSuite")
```

## Quick Start

### Query Gene Expression

```r
library(ZinaSuite)

# Query TP53 expression in TCGA BRCA
data <- query_gene_expression("TP53", cancer = "BRCA")
head(data)

# Query multiple genes
multi_data <- query_molecules(c("TP53", "BRCA1", "EGFR"), cancer = "BRCA")
```

### Correlation Analysis

```r
# Analyze correlation between two genes
cor_result <- analyze_correlation("TP53", "BRCA1", cancer = "BRCA")
print(cor_result)

# Visualize correlation
vis_correlation("TP53", "BRCA1", cancer = "BRCA")
```

### Survival Analysis

```r
# Survival analysis by gene expression
surv_result <- analyze_survival_by_expression("TP53", cancer = "BRCA")
vis_survival_by_gene("TP53", cancer = "BRCA")
```

### Launch Shiny App

```r
# Start the interactive Shiny application
run_zinasuite()
```

## Shiny Application

ZinaSuite includes a comprehensive Shiny application with 11 functional modules:

| Module | Description |
|--------|-------------|
| **Home** | Welcome page and quick start guide |
| **Data Query** | Query gene expression, mutation, CNV, methylation data |
| **Analysis** | Correlation and survival analysis |
| **Visualization** | Distribution plots, correlation plots, heatmaps |
| **Pan-Cancer** | Pan-cancer analysis across 33 cancer types |
| **Mutation** | Mutation frequency and co-occurrence analysis |
| **Dimension** | PCA, t-SNE, UMAP dimensionality reduction |
| **Immune** | Immune infiltration and correlation analysis |
| **PharmacoGenomics** | Drug sensitivity and pharmacogenomics analysis |
| **Batch** | Batch processing for multiple genes |
| **About** | Package information and documentation |

### Shiny App Architecture

The Shiny app follows a modular architecture with:

- **UI/Server Separation**: Clean separation of concerns
- **Business Logic Extraction**: Testable functions independent of Shiny
- **Reactive State Management**: Shared app state across modules
- **Async Computing**: Non-blocking operations with mirai
- **Progress Tracking**: Real-time progress updates

## Advanced Features

### Caching System

```r
# Use cached queries for better performance
data <- cached_query(
  function() query_gene_expression("TP53", cancer = "BRCA"),
  cache_key = "TP53_BRCA",
  ttl = 3600  # Cache for 1 hour
)

# Manage cache
cache <- get_zina_cache()
cache$stats()  # View cache statistics
clear_zina_cache()  # Clear all cache
```

### Batch Processing

```r
# Create analysis pipeline
pipeline <- create_pipeline(list(
  query = function(gene) query_gene_expression(gene, cancer = "BRCA"),
  analyze = function(data) analyze_survival_by_expression(data),
  visualize = function(result) vis_survival(result)
), name = "Gene Survival Analysis")

# Run pipeline
result <- run_pipeline(pipeline, gene = "TP53")
```

### Report Generation

```r
# Generate analysis report
result <- analyze_correlation("TP53", "BRCA1", cancer = "BRCA")
generate_report(
  result,
  output_file = "analysis_report.html",
  output_format = "html",
  title = "TP53-BRCA1 Correlation Analysis"
)

# Export data
export_data(result, "results.xlsx", format = "excel")
```

### Data Validation

```r
# Validate inputs before analysis
validate_gene_symbol("TP53")
validate_cancer_type("BRCA")
validate_sample_ids(c("TCGA-01-1234-01", "TCGA-02-5678-01"))

# Data quality report
report <- data_quality_report(expression_data)
print(report)
```

## API Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `query_gene_expression()` | Query gene expression data |
| `query_mutation()` | Query mutation data |
| `query_cnv()` | Query copy number variation |
| `analyze_correlation()` | Correlation analysis |
| `analyze_survival()` | Survival analysis |
| `vis_correlation()` | Correlation visualization |
| `vis_survival()` | Survival visualization |

### Shiny Functions

| Function | Description |
|----------|-------------|
| `run_zinasuite()` | Launch Shiny application |
| `build_query_params()` | Build query parameters |
| `validate_query_input()` | Validate query inputs |
| `execute_data_query()` | Execute data query |
| `create_batch_job()` | Create batch processing job |

## Testing

ZinaSuite includes comprehensive test suites:

```r
# Run all tests
devtools::test()

# Run specific test file
devtools::test(filter = "shiny-logic")
```

**Test Coverage:** 319 tests covering core functions, Shiny logic, data validation, and visualization.

## Comparison with UCSCXenaShiny

| Feature | UCSCXenaShiny | ZinaSuite |
|---------|---------------|-----------|
| Architecture | Monolithic | Modular |
| Async Computing | Limited | Full mirai support |
| Test Coverage | Basic | Comprehensive (319 tests) |
| Caching | None | Memory + Disk with TTL |
| Report Generation | Limited | HTML/PDF/Word/Excel |
| Data Validation | Basic | Extensive |
| Code Organization | Mixed | Clean separation |

## Documentation

- [Function Reference](https://wanglabcsu.github.io/ZinaSuite/reference/)
- [Vignettes](https://wanglabcsu.github.io/ZinaSuite/articles/)
- [Shiny App Guide](https://wanglabcsu.github.io/ZinaSuite/articles/shiny-app.html)

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## Citation

If you use ZinaSuite in your research, please cite:

```
Wang S, et al. (2025). ZinaSuite: A Modern R Package for UCSC Xena Data Analysis. 
Journal of Open Source Software, 10(XXX), XXXX.
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- UCSC Xena Team for providing the data hub
- UCSCXenaTools package for the foundation
- All contributors to this project

## Contact

- **Issues**: [GitHub Issues](https://github.com/WangLabCSU/ZinaSuite/issues)
- **Email**: w_shixiang@163.com
- **Website**: https://github.com/WangLabCSU/ZinaSuite

---

**Note**: This package is under active development. Please report any issues or feature requests on GitHub.
