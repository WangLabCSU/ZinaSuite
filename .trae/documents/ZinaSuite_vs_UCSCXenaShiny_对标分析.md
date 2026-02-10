# ZinaSuite vs UCSCXenaShiny 功能对标分析

## 分析日期
2026-02-09

## 1. 功能覆盖对比总览

| 功能类别 | UCSCXenaShiny | ZinaSuite | 覆盖状态 |
|---------|--------------|-----------|---------|
| **核心数据获取** | | | |
| TCGA基因表达 | ✅ get_pancan_gene_value | ✅ query_gene_expression | 🟢 已覆盖 |
| TCGA转录本表达 | ✅ get_pancan_transcript_value | ⚠️ 部分支持 | 🟡 需完善 |
| TCGA蛋白质表达 | ✅ get_pancan_protein_value | ✅ query_protein | 🟢 已覆盖 |
| TCGA突变状态 | ✅ get_pancan_mutation_status | ✅ query_mutation | 🟢 已覆盖 |
| TCGA拷贝数变异 | ✅ get_pancan_cn_value | ✅ query_cnv | 🟢 已覆盖 |
| TCGA甲基化 | ✅ get_pancan_methylation_value | ✅ query_methylation | 🟢 已覆盖 |
| TCGA miRNA | ✅ get_pancan_miRNA_value | ✅ query_mirna | 🟢 已覆盖 |
| PCAWG数据获取 | ✅ 6个函数 | ⚠️ 基础支持 | 🟡 需扩展 |
| CCLE数据获取 | ✅ 4个函数 | ⚠️ 基础支持 | 🟡 需扩展 |
| **数据查询** | | | |
| 统一查询接口 | ✅ query_pancan_value | ✅ query_molecule | 🟢 已覆盖 |
| 批量查询 | ✅ | ✅ query_molecules + mirai | 🟢 已覆盖 |
| 样本分组查询 | ✅ query_tcga_group | ❌ 未实现 | 🔴 缺失 |
| **可视化** | | | |
| 肿瘤vs正常图 | ✅ vis_toil_TvsN | ⚠️ vis_tumor_normal (简化版) | 🟡 需完善 |
| Cox森林图 | ✅ vis_unicox_tree | ✅ vis_unicox_forest | 🟢 已覆盖 |
| 解剖位置图 | ✅ vis_pancan_anatomy | ❌ 未实现 | 🔴 缺失 |
| 免疫相关性热图 | ✅ vis_gene_immune_cor | ❌ 未实现 | 🔴 缺失 |
| TIL相关性 | ✅ vis_gene_TIL_cor | ❌ 未实现 | 🔴 缺失 |
| TMB相关性 | ✅ vis_gene_tmb_cor | ❌ 未实现 | 🔴 缺失 |
| MSI相关性 | ✅ vis_gene_msi_cor | ❌ 未实现 | 🔴 缺失 |
| 干性相关性 | ✅ vis_gene_stemness_cor | ❌ 未实现 | 🔴 缺失 |
| 基因-基因相关性 | ✅ vis_gene_cor | ✅ vis_gene_correlation | 🟢 已覆盖 |
| 泛癌相关性 | ✅ vis_gene_cor_cancer | ✅ vis_pancan_correlation | 🟢 已覆盖 |
| 突变vs野生型 | ✅ vis_toil_Mut | ❌ 未实现 | 🔴 缺失 |
| K-M生存曲线 | ✅ tcga_surv_plot | ✅ vis_survival | 🟢 已覆盖 |
| 降维可视化 | ✅ vis_dim_dist | ❌ 未实现 | 🔴 缺失 |
| 交叉组学图 | ✅ vis_gene_cross_omics | ❌ 未实现 | 🔴 缺失 |
| **分析功能** | | | |
| 相关性分析 | ✅ ezcor | ✅ analyze_correlation | 🟢 已覆盖 |
| 批量相关性 | ✅ ezcor_batch | ✅ analyze_correlation_batch | 🟢 已覆盖 |
| 偏相关分析 | ✅ ezcor_partial_cor | ✅ analyze_partial_correlation | 🟢 已覆盖 |
| 相关矩阵 | ❌ | ✅ analyze_correlation_matrix | 🟢 扩展功能 |
| 生存分析 | ✅ tcga_surv_get/plot | ✅ analyze_survival | 🟢 已覆盖 |
| 批量Cox分析 | ❌ | ✅ analyze_unicox_batch | 🟢 扩展功能 |
| 基因表达分组生存 | ❌ | ✅ analyze_survival_by_expression | 🟢 扩展功能 |
| **Shiny应用** | | | |
| 首页模块 | ✅ | ⚠️ 基础框架 | 🟡 需完善 |
| 快速分析模块 | ✅ (17个) | ❌ 未实现 | 🔴 缺失 |
| 高级分析模块 | ✅ (27个) | ❌ 未实现 | 🔴 缺失 |
| 通用分析模块 | ✅ (5个) | ❌ 未实现 | 🔴 缺失 |
| 药物基因组学 | ✅ (5个) | ❌ 未实现 | 🔴 缺失 |
| 数据下载模块 | ✅ | ❌ 未实现 | 🔴 缺失 |
| **数据基础设施** | | | |
| 内置数据集 | ✅ 16个 | ⚠️ 框架已建 | 🟡 需填充 |
| 远程数据支持 | ✅ 20+ | ⚠️ 框架已建 | 🟡 需填充 |
| 缓存管理 | ✅ 基础版 | ✅ R6 CacheManager | 🟢 已增强 |
| 异步计算 | ❌ future | ✅ mirai AsyncCompute | 🟢 已增强 |

## 2. 详细功能差距分析

### 2.1 数据获取功能差距

#### PCAWG数据支持不足
**UCSCXenaShiny功能：**
- get_pcawg_gene_value - 基因表达
- get_pcawg_fusion_value - 基因融合
- get_pcawg_promoter_value - 启动子活性
- get_pcawg_miRNA_value - miRNA表达
- get_pcawg_APOBEC_mutagenesis_value - APOBEC突变

**ZinaSuite现状：**
- 仅通过source="pcawg"参数支持基础查询
- 缺乏PCAWG特有数据类型支持

**优先级：** 中
**实现建议：** 扩展XenaData类，添加PCAWG特有方法

#### CCLE数据支持不足
**UCSCXenaShiny功能：**
- get_ccle_gene_value
- get_ccle_protein_value
- get_ccle_cn_value
- get_ccle_mutation_status

**ZinaSuite现状：**
- 仅通过source="ccle"参数支持基础查询

**优先级：** 中
**实现建议：** 类似PCAWG，添加CCLE特有方法

### 2.2 可视化功能差距

#### 缺失的关键可视化（按优先级排序）

1. **解剖位置表达图 (vis_pancan_anatomy)**
   - 功能：在人体解剖图上展示基因表达
   - 优先级：高
   - 实现复杂度：高（需要解剖图数据）

2. **免疫特征相关性热图 (vis_gene_immune_cor)**
   - 功能：基因与28种免疫特征的相关性
   - 优先级：高
   - 依赖：需要免疫特征数据集

3. **TIL相关性分析 (vis_gene_TIL_cor)**
   - 功能：基因与肿瘤浸润淋巴细胞相关性
   - 优先级：高
   - 依赖：需要TIL数据集

4. **突变vs野生型比较 (vis_toil_Mut)**
   - 功能：比较突变和野生型样本的分子特征
   - 优先级：高
   - 实现：相对简单，基于现有分组比较

5. **降维可视化 (vis_dim_dist)**
   - 功能：PCA/UMAP/tSNE降维展示
   - 优先级：中
   - 依赖：需要Rtsne/uwot包

6. **交叉组学可视化 (vis_gene_cross_omics)**
   - 功能：funkyheatmap展示多组学数据
   - 优先级：中
   - 依赖：需要funkyheatmap包

### 2.3 分析功能差距

#### 样本分组查询 (query_tcga_group)
**UCSCXenaShiny功能：**
- 按内置表型分组（癌症类型、分期等）
- 自定义分组
- 支持过滤和合并操作
- 返回分组后的样本信息

**ZinaSuite现状：**
- 缺乏分组查询功能

**优先级：** 高
**实现建议：** 新增analysis-grouping.R模块

### 2.4 Shiny应用差距

#### 模块缺失统计
- 快速分析模块：17个模块缺失
- 高级分析模块：27个模块缺失
- 通用分析模块：5个模块缺失
- 药物基因组学：5个模块缺失

**优先级：** 中（基础框架已建立）
**实现建议：** 按模块逐步移植，复用现有R函数

## 3. 数据基础设施对比

### 3.1 内置数据集

| 数据集 | UCSCXenaShiny | ZinaSuite | 状态 |
|--------|--------------|-----------|------|
| tcga_clinical | ✅ | ⚠️ 待添加 | 🔴 |
| tcga_surv | ✅ | ⚠️ 待添加 | 🔴 |
| tcga_gtex | ✅ | ⚠️ 待添加 | 🔴 |
| tcga_purity | ✅ | ⚠️ 待添加 | 🔴 |
| pcawg_info | ✅ | ⚠️ 待添加 | 🔴 |
| ccle_info | ✅ | ⚠️ 待添加 | 🔴 |

### 3.2 远程数据集

| 数据集 | UCSCXenaShiny | ZinaSuite | 状态 |
|--------|--------------|-----------|------|
| tcga_TIL | ✅ | ⚠️ 待添加 | 🔴 |
| tcga_PW | ✅ | ⚠️ 待添加 | 🔴 |
| tcga_tmb | ✅ | ⚠️ 待添加 | 🔴 |
| tcga_MSI | ✅ | ⚠️ 待添加 | 🔴 |
| tcga_stemness | ✅ | ⚠️ 待添加 | 🔴 |

## 4. 架构优势对比

### 4.1 ZinaSuite优势

1. **现代R6类架构**
   - CacheManager：LRU缓存 + 内存/磁盘双存储
   - AsyncCompute：mirai异步计算引擎
   - DataSource/XenaData：可扩展的数据源抽象

2. **高性能异步计算**
   - 使用mirai替代future/furrr
   - 支持批量并行处理
   - 内置进度跟踪

3. **代码质量**
   - 完整的roxygen2文档
   - 可运行的示例（\donttest）
   - testthat测试框架

### 4.2 UCSCXenaShiny优势

1. **功能完整性**
   - 50+可视化函数
   - 完整的Shiny应用
   - 丰富的内置数据集

2. **生态系统**
   - 教程文档完善
   - 用户社区活跃
   - 持续更新维护

## 5. 开发优先级建议

### 阶段1：核心功能补齐（高优先级）

1. **数据基础设施**
   - [ ] 添加内置数据集（tcga_gtex, tcga_surv等）
   - [ ] 实现远程数据加载（TIL, TMB, MSI等）
   - [ ] 添加样本分组查询功能

2. **关键可视化**
   - [ ] 免疫相关性热图
   - [ ] TIL相关性分析
   - [ ] 突变vs野生型比较
   - [ ] 解剖位置表达图

### 阶段2：数据扩展（中优先级）

1. **PCAWG数据支持**
   - [ ] 基因融合数据
   - [ ] 启动子活性数据
   - [ ] APOBEC突变特征

2. **CCLE数据支持**
   - [ ] 细胞系特有功能
   - [ ] 药物响应数据集成

### 阶段3：Shiny应用完善（中优先级）

1. **核心模块**
   - [ ] 首页模块
   - [ ] 快速分析模块（TCGA）
   - [ ] 数据下载模块

2. **高级功能**
   - [ ] 高级分析模块
   - [ ] 通用分析模块
   - [ ] 药物基因组学模块

### 阶段4：高级功能（低优先级）

1. **分析增强**
   - [ ] 降维分析
   - [ ] 交叉组学分析
   - [ ] 通路分析

2. **可视化增强**
   - [ ] 交互式图形
   - [ ] 3D可视化
   - [ ] 动画展示

## 6. 技术债务与改进

### 6.1 已知问题

1. **mirai worker限制**
   - 问题：worker进程无法访问包函数
   - 影响：批量查询时可能失败
   - 解决：需要在worker中显式加载包

2. **UCSCXenaTools依赖**
   - 问题：数据源依赖外部包
   - 风险：API变更可能导致功能中断
   - 建议：添加数据获取的fallback机制

### 6.2 性能优化机会

1. **缓存策略优化**
   - 当前：简单的LRU策略
   - 改进：按数据类型和访问频率优化

2. **并行计算优化**
   - 当前：基础mirai_map使用
   - 改进：动态worker数量调整

## 7. 总结

### 功能完成度评估

- **核心数据获取**: 70% ✅
- **数据查询**: 80% ✅
- **基础可视化**: 60% 🟡
- **高级可视化**: 30% 🔴
- **分析功能**: 75% ✅
- **Shiny应用**: 10% 🔴
- **数据基础设施**: 40% 🟡

**总体完成度：约 50%**

### 下一步行动建议

1. **立即行动**：补充内置数据集，这是其他功能的基础
2. **短期目标**：实现关键可视化（免疫、TIL、突变比较）
3. **中期目标**：完善Shiny应用核心模块
4. **长期目标**：达到与UCSCXenaShiny功能对等
