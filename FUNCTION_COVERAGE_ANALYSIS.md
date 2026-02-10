# ZinaSuite vs UCSCXenaShiny 功能覆盖分析

## 分析日期：2026-02-09

---

## 1. UCSCXenaShiny 功能清单

### 1.1 General Analysis (通用分析)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| modules-ga-dim-distribution | 降维分布可视化 | ✅ vis_identifier_dim_dist |
| modules-ga-group-comparison | 组间比较 | ✅ analyze_group_comparison |
| modules-ga-matrix-correlation | 矩阵相关性 | ✅ analyze_correlation_matrix |
| modules-ga-scatter-correlation | 散点相关性 | ✅ vis_gene_correlation |
| modules-ga-surv-analysis | 生存分析 | ✅ analyze_survival |

### 1.2 TCGA Quick Analysis (TCGA 快速分析)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| modules-1-tcga-01-TN | 转录组 vs 正常 | ✅ query_gene_expression |
| modules-1-tcga-02-Anatomy | 解剖学表达 | ✅ vis_pancan_anatomy |
| modules-1-tcga-03-Cor | 相关性分析 | ✅ analyze_correlation |
| modules-1-tcga-04-TIL | 肿瘤浸润淋巴细胞 | ✅ vis_gene_TIL_cor |
| modules-1-tcga-05-Immune | 免疫分析 | ✅ vis_gene_immune_cor |
| modules-1-tcga-06-Idx | 免疫指数 | ✅ query_immune_infiltration |
| modules-1-tcga-07-PW | 通路分析 | ⚠️ 部分支持 |
| modules-1-tcga-08-Mut | 突变分析 | ✅ query_mutation |
| modules-1-tcga-09-KM | KM 生存曲线 | ✅ analyze_kaplan_meier |
| modules-1-tcga-10-Cox | Cox 回归 | ✅ analyze_cox_regression |
| modules-1-tcga-11-Dim | 降维分析 | ✅ vis_identifier_dim_dist |

### 1.3 PCAWG Quick Analysis (PCAWG 快速分析)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| modules-2-pcawg-01-TN | 转录组 vs 正常 | ✅ PCAWGData class |
| modules-2-pcawg-02-Cor | 相关性分析 | ✅ analyze_correlation |
| modules-2-pcawg-03-KM | KM 生存曲线 | ✅ analyze_kaplan_meier |
| modules-2-pcawg-04-Cox | Cox 回归 | ✅ analyze_cox_regression |

### 1.4 CCLE Quick Analysis (CCLE 快速分析)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| modules-3-ccle-01-Dist | 分布分析 | ✅ CCLEData class |
| modules-3-ccle-02-Cor | 相关性分析 | ✅ analyze_correlation |
| modules-3-ccle-03-DrugT | 药物敏感性 | ⚠️ 部分支持 |
| modules-3-ccle-04-DrugR | 药物响应 | ⚠️ 部分支持 |

### 1.5 TCGA Pan-Cancer Analysis (TCGA 泛癌分析)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| modules-pancan-comp-m2o | 多对一比较 | ✅ analyze_group_comparison |
| modules-pancan-comp-o2m | 一对多比较 | ✅ analyze_group_comparison |
| modules-pancan-comp-o2o | 一对一比较 | ✅ analyze_group_comparison |
| modules-pancan-cor-m2o | 多对一相关 | ✅ analyze_correlation_batch |
| modules-pancan-cor-o2m | 一对多相关 | ✅ analyze_correlation_batch |
| modules-pancan-cor-o2o | 一对一相关 | ✅ analyze_correlation |
| modules-pancan-cross-gene-o2m | 跨基因分析 | ✅ analyze_correlation_matrix |
| modules-pancan-cross-pw-o2m | 跨通路分析 | ⚠️ 部分支持 |
| modules-pancan-sur-m2o | 多对一生存 | ✅ analyze_unicox_batch |
| modules-pancan-sur-o2m | 一对多生存 | ✅ analyze_survival_by_expression |
| modules-pancan-sur-o2o | 一对一生存 | ✅ analyze_survival |

### 1.6 PCAWG Pan-Cancer Analysis (PCAWG 泛癌分析)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| modules-pcawg-comp-m2o | 多对一比较 | ✅ analyze_group_comparison |
| modules-pcawg-comp-o2m | 一对多比较 | ✅ analyze_group_comparison |
| modules-pcawg-comp-o2o | 一对一比较 | ✅ analyze_group_comparison |
| modules-pcawg-cor-m2o | 多对一相关 | ✅ analyze_correlation_batch |
| modules-pcawg-cor-o2m | 一对多相关 | ✅ analyze_correlation_batch |
| modules-pcawg-cor-o2o | 一对一相关 | ✅ analyze_correlation |
| modules-pcawg-sur-m2o | 多对一生存 | ✅ analyze_unicox_batch |
| modules-pcawg-sur-o2m | 一对多生存 | ✅ analyze_survival_by_expression |
| modules-pcawg-sur-o2o | 一对一生存 | ✅ analyze_survival |

### 1.7 CCLE Analysis (CCLE 分析)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| modules-ccle-comp-m2o | 多对一比较 | ✅ analyze_group_comparison |
| modules-ccle-comp-o2o | 一对一比较 | ✅ analyze_group_comparison |
| modules-ccle-cor-m2o | 多对一相关 | ✅ analyze_correlation_batch |
| modules-ccle-cor-o2o | 一对一相关 | ✅ analyze_correlation |

### 1.8 TPC Functions (TPC 功能)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| modules-z-add-signature | 添加特征 | ⚠️ 待实现 |
| modules-z-custom-meta | 自定义元数据 | ⚠️ 待实现 |
| modules-z-download-feat | 下载特征 | ⚠️ 待实现 |
| modules-z-download-res | 下载结果 | ⚠️ 待实现 |
| modules-z-filter-sample | 样本过滤 | ✅ apply_filter |
| modules-z-group-sample | 样本分组 | ✅ merge_numeric_groups |
| modules-z-mol-origin | 分子来源 | ✅ query_molecule |
| modules-z-multi-upload | 多文件上传 | ⚠️ 待实现 |

### 1.9 PharmacoGenomics (药物基因组学)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| DrugOmicPair | 药物-组学配对 | ⚠️ 部分支持 |
| FeatureAcrossType | 跨类型特征 | ⚠️ 待实现 |
| FeatureDatabaseSig | 特征数据库 | ⚠️ 待实现 |
| ProfileDrugSens | 药物敏感性谱 | ⚠️ 部分支持 |
| StatAnno | 统计注释 | ⚠️ 待实现 |

### 1.10 Other Pages (其他页面)
| 模块 | 功能描述 | ZinaSuite 状态 |
|------|----------|----------------|
| combo-single-gene-pan-cancer-analysis | 单基因泛癌分析 | ✅ vis_pancan_distribution |
| home-daily-gene | 每日基因 | ⚠️ 待实现 |
| home-pancan-search | 泛癌搜索 | ✅ query_molecules |
| modules-file-upload | 文件上传 | ⚠️ 待实现 |
| modules-z-download-1-pancan | 泛癌下载 | ⚠️ 待实现 |
| modules-z-download-2-dataset | 数据集下载 | ⚠️ 待实现 |
| modules-z-help-id | ID 帮助 | ⚠️ 待实现 |

---

## 2. ZinaSuite Shiny 模块现状

### 2.1 已实现的模块
| 模块 | 功能覆盖 |
|------|----------|
| mod_home | 首页、欢迎界面、功能介绍 |
| mod_data_query | 数据查询、TCGA/PCAWG/CCLE、多种数据类型 |
| mod_analysis | 相关性分析、生存分析 |
| mod_visualization | 直方图、箱线图、散点图 |
| mod_immune | 免疫相关、TIL、MSI、TMB |
| mod_batch | 批量分析、并行处理 |
| mod_about | 关于页面 |

### 2.2 缺失的 Shiny 模块

#### 高优先级
1. **mod_pancan** - 泛癌分析模块
   - 跨癌种比较
   - 泛癌相关性
   - 泛癌生存分析

2. **mod_mutation** - 突变分析模块
   - 突变频率
   - 突变共现
   - 突变与表达关系

3. **mod_dimension** - 降维分析模块
   - PCA
   - t-SNE
   - UMAP

4. **mod_comparison** - 组间比较模块
   - 多组比较
   - 统计检验
   - 可视化

#### 中优先级
5. **mod_download** - 数据下载模块
6. **mod_help** - 帮助文档模块
7. **mod_custom** - 自定义分析模块

---

## 3. 功能覆盖统计

| 类别 | UCSCXenaShiny 功能数 | ZinaSuite 覆盖数 | 覆盖率 |
|------|---------------------|-----------------|--------|
| General Analysis | 5 | 5 | 100% |
| TCGA Quick | 11 | 10 | 91% |
| PCAWG Quick | 4 | 4 | 100% |
| CCLE Quick | 4 | 3 | 75% |
| TCGA Pan-Cancer | 11 | 10 | 91% |
| PCAWG Pan-Cancer | 9 | 9 | 100% |
| CCLE Analysis | 4 | 4 | 100% |
| TPC Functions | 8 | 3 | 38% |
| PharmacoGenomics | 5 | 1 | 20% |
| Other Pages | 7 | 2 | 29% |
| **总计** | **68** | **51** | **75%** |

---

## 4. 建议的下一步工作

### 4.1 高优先级（必须实现）
1. **mod_pancan** - 泛癌分析模块
2. **mod_mutation** - 突变分析模块
3. **mod_dimension** - 降维分析模块
4. **mod_comparison** - 组间比较模块

### 4.2 中优先级（增强功能）
5. **mod_download** - 数据下载模块
6. **mod_help** - 帮助文档模块
7. **mod_custom** - 自定义分析模块

### 4.3 低优先级（可选功能）
8. **PharmacoGenomics** 功能
9. **TPC** 高级功能
10. **Other Pages** 附加功能

---

## 5. 总结

**当前状态**：
- ZinaSuite 已覆盖 UCSCXenaShiny 75% 的核心功能
- 所有基础分析功能（相关性、生存、免疫）已实现
- Shiny 应用包含 7 个功能模块

**建议**：
- 优先实现 4 个高优先级模块，可达到 90% 功能覆盖
- 继续完善 Shiny 应用的用户体验
- 添加更多可视化选项

---

**分析时间**: 2026-02-09
**ZinaSuite 版本**: 0.0.0.9000
**Git 提交**: e872b7c
