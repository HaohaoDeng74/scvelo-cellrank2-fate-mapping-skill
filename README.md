# scVelo + CellRank2 Fate Mapping Skill

一个面向 **Codex / Agent / 生信分析自动化流程** 的 scVelo + CellRank2 细胞命运轨迹分析 skill 包。

本项目用于从单细胞转录组 AnnData 对象出发，完成 RNA velocity 质控、scVelo 动态建模、CellRank2 命运概率推断、terminal / initial state 识别、候选命运相关驱动基因分析、基因动态趋势分析和报告生成。

本 skill 的核心目标不是“画一张 velocity stream 图”，而是建立一套更接近发表级分析逻辑的证据链：

```text
输入是否可靠
→ 是否存在可解释的动态过程
→ velocity 方向是否可信
→ CellRank fate probability 是否稳定
→ driver genes 是否只是相关候选还是具备机制解释价值
→ 最终结论应如何分级表达
```

---

## 1. 项目定位

本项目适用于以下分析需求：

* scVelo RNA velocity 分析
* CellRank / CellRank2 fate mapping
* initial state / terminal state 推断
* fate probability 分析
* latent time / pseudotime 动态分析
* lineage-associated driver candidate genes 分析
* 感染、免疫、发育、分化、激活、耗竭、再生等动态单细胞过程分析
* 面向论文或项目报告的命运轨迹图件和解释框架生成

需要强调的是：本项目输出的 driver genes 应理解为 **fate-associated candidate genes**，即与命运概率相关的候选基因，不应直接解释为因果驱动基因。因果结论需要扰动实验、流式、空间组、谱系示踪或其他正交证据支持。

---

## 2. 主要功能

### 2.1 输入验证

`00_validate_input.py` 会检查：

* `.h5ad` 文件是否存在
* `adata.X` / `adata.layers` 维度是否合理
* 是否包含 `spliced` 和 `unspliced`
* 是否包含 UMAP / PCA / neighbors
* 是否包含 cell type、sample、condition、timepoint 等元数据
* lineage subset 后细胞数是否足够
* unspliced fraction 是否异常
* 是否推荐走 velocity route 或 non-velocity route

### 2.2 scVelo 分析

`01_run_scvelo.py` 支持：

* scVelo 预处理
* moments 计算
* dynamical mode velocity
* stochastic fallback
* velocity graph
* latent time
* velocity 结果保存
* route / warning / fallback 信息写入 `run_manifest.json`

### 2.3 velocity 质控与可视化

`02_velocity_qc.py` 输出：

* velocity stream
* velocity confidence
* velocity length
* latent time
* root / end point 相关图件
* top velocity gene phase portrait
* 实际生成图件路径记录

### 2.4 CellRank2 fate mapping

`03_run_cellrank2.py` 支持：

* VelocityKernel
* ConnectivityKernel
* PseudotimeKernel
* RealTimeKernel
* hybrid kernel
* GPCCA estimator
* terminal states
* initial states
* fate probabilities
* coarse transition matrix
* n_states sensitivity
* fate assignment stability

如果缺少 `spliced/unspliced`，脚本会自动进入 non-velocity route，并在报告中降低方向性解释等级。

### 2.5 driver candidate 分析

`04_driver_gene_analysis.py` 支持两种模式：

1. 优先尝试 CellRank native `compute_lineage_drivers()`
2. 如果失败，则使用 correlation fallback，即基因表达与 fate probability 的相关性分析

输出表会标记：

* method
* lineage
* gene
* correlation / p value / q value
* risk flag
* mitochondrial / ribosomal / hemoglobin / cell-cycle / stress / ambient RNA 风险基因

### 2.6 gene trend 与 pathway 支持

`05_gene_trend_and_pathway.py` 支持：

* top driver candidates 的动态趋势汇总
* latent time / pseudotime 分箱趋势
* 可选 GMT gene set enrichment
* 使用 expressed genes 作为背景基因集，避免使用全基因组作为不合理背景

### 2.7 自动报告与图件规划

`06_generate_report.py` 输出：

* `06_report/scvelo_cellrank2_report.md`
* `06_report/report_status.json`
* `07_publication_figures/figure_plan.md`

报告包括：

* 输入 QC
* velocity route 状态
* CellRank2 route 状态
* fate probability 输出状态
* driver candidate 输出状态
* gene trend 输出状态
* directionality support
* evidence grade
* 主要 warning 和失败模式

---

## 3. 推荐目录结构

```text
scvelo-cellrank2-fate-mapping-skill/
├── SKILL.md
├── .gitignore
├── agents/
│   └── openai.yaml
├── config_template.yaml
├── scripts/
│   ├── 00_validate_input.py
│   ├── 01_run_scvelo.py
│   ├── 02_velocity_qc.py
│   ├── 03_run_cellrank2.py
│   ├── 04_driver_gene_analysis.py
│   ├── 05_gene_trend_and_pathway.py
│   ├── 06_generate_report.py
│   ├── smoke_test.py
│   ├── workflow_test_minimal.py
│   └── utils.py
├── references/
│   ├── api_compatibility.md
│   ├── cns_application_patterns.md
│   ├── failure_modes.md
│   ├── literature_guided_interpretation.md
│   ├── qc_checklist.md
│   ├── visualization_spec.md
│   └── v4_hygiene_revision_notes.md
├── templates/
│   ├── codex_prompt.md
│   ├── figure_plan_template.md
│   └── report_template.md
└── examples/
    └── clp_spleen_custom_example.md
```

---

## 4. 安装环境

推荐使用独立 conda 环境。

```bash
conda create -n scvelo_cellrank2 python=3.10 -y
conda activate scvelo_cellrank2
```

安装基础依赖：

```bash
pip install numpy pandas scipy scikit-learn pyyaml matplotlib
pip install anndata scanpy scvelo cellrank
```

如需 pathway enrichment，可根据项目需要安装：

```bash
pip install gseapy
```

注意：CellRank、scVelo、Scanpy、AnnData 的版本组合可能影响 API 兼容性。正式运行前建议先执行 `smoke_test.py` 和 `workflow_test_minimal.py`。

---

## 5. 快速开始

### 5.1 克隆仓库

```bash
git clone https://github.com/HaohaoDeng74/scvelo-cellrank2-fate-mapping-skill.git
cd scvelo-cellrank2-fate-mapping-skill
```

### 5.2 复制配置文件

```bash
cp config_template.yaml config.yaml
```

编辑 `config.yaml`，至少修改：

```yaml
project:
  input_h5ad: /path/to/your/input.h5ad
  output_dir: results

metadata:
  celltype_key: celltype
  sample_key: sample
  condition_key: condition
  timepoint_key: timepoint
```

如果只分析某个 lineage subset，例如髓系细胞：

```yaml
subset:
  enabled: true
  subset_name: myeloid_lineage
  obs_key: major_celltype
  include_values:
    - Monocyte
    - Macrophage
    - Neutrophil
```

### 5.3 运行分析流程

```bash
python scripts/00_validate_input.py --config config.yaml
python scripts/01_run_scvelo.py --config config.yaml
python scripts/02_velocity_qc.py --config config.yaml
python scripts/03_run_cellrank2.py --config config.yaml
python scripts/04_driver_gene_analysis.py --config config.yaml
python scripts/05_gene_trend_and_pathway.py --config config.yaml
python scripts/06_generate_report.py --config config.yaml
```

### 5.4 运行轻量检查

```bash
python scripts/smoke_test.py --skill-root .
```

如果环境已安装 `anndata / scanpy / cellrank`，可以运行最小 workflow 测试：

```bash
python scripts/workflow_test_minimal.py --skill-root .
```

---

## 6. 输入文件要求

理想输入是 `.h5ad` 文件，并包含：

```text
adata.X
adata.obs["celltype"] 或 adata.obs["leiden"]
adata.obs["sample"]
adata.obs["condition"]
adata.obs["timepoint"]
adata.obsm["X_umap"]
adata.layers["spliced"]
adata.layers["unspliced"]
```

如果缺少 `spliced` 和 `unspliced`，不能运行标准 RNA velocity。此时 workflow 会进入 non-velocity CellRank route，例如：

* PseudotimeKernel
* RealTimeKernel
* ConnectivityKernel

其中 connectivity-only route 不应解释为强方向性命运推断，只能作为 state-structure mapping 或探索性 fate structure 分析。

---

## 7. 输出结果

默认输出目录为：

```text
results/
├── run_manifest.json
├── 00_input_qc/
├── 01_velocity/
├── 02_velocity_qc/
├── 03_cellrank_fate/
├── 04_driver_genes/
├── 05_gene_trends/
├── 06_report/
└── 07_publication_figures/
```

关键输出文件包括：

```text
00_input_qc/input_validation.json
01_velocity/adata_scvelo.h5ad
02_velocity_qc/velocity_qc_summary.json
03_cellrank_fate/adata_cellrank2.h5ad
03_cellrank_fate/fate_probabilities.csv
03_cellrank_fate/terminal_states.csv
03_cellrank_fate/n_states_sensitivity.csv
04_driver_genes/lineage_driver_candidates.csv
04_driver_genes/top_lineage_driver_candidates.csv
05_gene_trends/gene_trends_summary.csv
06_report/scvelo_cellrank2_report.md
07_publication_figures/figure_plan.md
```

---

## 8. 结果解释原则

### 8.1 velocity stream 不能单独作为结论

必须结合：

* velocity confidence
* latent time
* phase portrait
* known markers
* sample timepoint
* condition distribution
* terminal state plausibility
* batch / sample composition

### 8.2 driver gene 只能称为候选基因

推荐表述：

```text
fate-associated candidate genes
lineage-associated driver candidates
candidate genes correlated with fate probabilities
```

避免表述：

```text
causal driver
master regulator
determines cell fate
proves fate decision
```

### 8.3 复杂免疫组织建议先 subset

对于脾脏、肺、肝脏、感染组织、肿瘤微环境等复杂组织，不建议直接用 all cells 解释全局 velocity。推荐先按 major lineage 拆分：

* erythroid / EMH
* monocyte–macrophage
* neutrophil / granulocyte
* B cell–plasma cell
* T cell activation / exhaustion / Treg
* dendritic cell lineage

---

## 9. CLP 脾脏数据建议

对于 CLP 脾脏单细胞数据，建议优先分析以下轴线：

### 9.1 Erythroid / EMH lineage

关注：

* CLP 是否诱导 extramedullary hematopoiesis
* 红系成熟或应激红系生成轨迹是否改变

参考 marker：

```text
Gata1, Klf1, Alas2, Hbb-bs, Hba-a1, Mki67, Top2a
```

### 9.2 Monocyte–macrophage lineage

关注：

* 炎症单核细胞是否向吞噬、修复或免疫抑制样 macrophage 状态偏移
* CLP 是否改变髓系状态命运概率

参考 marker：

```text
Ly6c2, S100a8, S100a9, Lcn2, Il1b, Ccl3, Ccl4, Mrc1, C1qa, C1qb, C1qc
```

### 9.3 B cell–plasma cell lineage

关注：

* B 细胞是否向 plasmablast / plasma cell 分化增强
* Ig 相关表达是否受到 CLP 后免疫重塑影响

参考 marker：

```text
Ms4a1, Cd79a, Cd74, Mzb1, Jchain, Xbp1, Prdm1
```

### 9.4 T cell activation / exhaustion / Treg axis

关注：

* T 细胞是否出现 activation、exhaustion、Treg bias 或长期免疫抑制状态

参考 marker：

```text
Cd3d, Cd4, Cd8a, Il7r, Ccr7, Mki67, Pdcd1, Lag3, Havcr2, Foxp3, Il2ra, Ctla4
```

---

## 10. 证据等级

自动报告中的 evidence grade 是计算审计等级，不等同于最终论文结论。

```text
Grade A:
velocity direction、latent time、fate probability、terminal states、marker dynamics、timepoint 和外部证据一致。

Grade B:
多数计算证据一致，但实验验证或文献支持不完整。

Grade C:
轨迹或 fate probability 有提示意义，但方向性或 driver 解释仍不稳。

Grade D:
velocity 假设失败、细胞群不相关、phase portrait 差、terminal states 不稳定或结果与已知生物学矛盾。
```

---

## 11. 仓库卫生规则

本仓库不应提交：

```text
__pycache__/
*.pyc
results/
run_manifest.json
*.h5ad
*.loom
*.h5
*.log
```

这些内容已写入 `.gitignore`。

正式发布前建议检查：

```bash
find . -name "__pycache__" -o -name "*.pyc" -o -name "run_manifest.json" -o -name "results"
```

---

## 12. 当前状态

当前版本定位为：

```text
release-candidate / internal beta
```

已完成：

* skill 结构整理
* YAML frontmatter
* 配置模板
* manifest 机制
* velocity / non-velocity route
* driver candidate 分析
* gene trend 基础分析
* report 生成
* smoke test
* 仓库卫生清理

仍建议在正式用于论文分析前完成：

* 使用真实或半真实 `.h5ad` 跑通完整 `00-06` 流程
* 确认 CellRank2 版本兼容性
* 检查 fate probabilities、terminal states、driver candidates 是否真实生成
* 结合 marker、timepoint、DEG/pathway、流式、空间组或文献进行人工解释
* 不将自动报告直接作为论文结论

---

## 13. 推荐引用与方法依据

本 workflow 的分析逻辑主要基于：

* RNA velocity：利用 spliced / unspliced RNA 推断细胞转录状态变化方向
* scVelo：通过 dynamical model 建模 splicing kinetics 和 transient cell states
* CellRank / CellRank2：通过 kernel 和 GPCCA 进行 fate probability、terminal state 和 driver candidate 分析
* 高水平文章中常见的多层证据链：方向性、时间点、marker、命运概率、调控程序和正交验证

---

## 14. 免责声明

本工具用于科研分析辅助和假设生成。自动输出的 trajectory、terminal states、fate probabilities 和 driver candidates 均需要结合生物学背景、样本设计、技术质量、批次结构和实验验证进行解释。

不得将本工具输出的候选 driver gene 直接解释为因果驱动基因。
