# AIDD8122 EGFR 活性预测：预先确定的分析方案

状态：2026-11-08 在下载数据、训练模型或查看测试结果之前制定并保存。实际执行记录见 `results/` 和 `RESULTS_SUMMARY.md`；以下分析规则保持原样。

## 研究问题

仅根据小分子的二维化学结构，能否预测其在体外实验中对 EGFR 的较强抑制活性？比较同一随机森林模型使用两种分子指纹时的表现：MACCS（基础）与 Morgan radius 2、2048 bit（改进）。这是一项分子活性预测研究，不预测肺癌患者治疗效果。

## 数据及标签

- 唯一主数据源：TeachOpenCADD T001 提供的冻结 ChEMBL 27 EGFR 文件 `EGFR_compounds_ea055ef.csv`，原始 5,568 行。原始 URL：<https://raw.githubusercontent.com/volkamerlab/teachopencadd/master/teachopencadd/talktorials/T001_query_chembl/data/EGFR_compounds_ea055ef.csv>。
- 下载后保存原文件、下载日期、SHA-256 和行数。训练数据一律由该文件重建。
- 保留可解析的 SMILES、有限且为正的 IC50（nM）及 pIC50。复核 `pIC50 = 9 - log10(IC50_nM)`；不一致的记录先核查，不能静默修改。
- 按 RDKit 规范化后的异构 SMILES 去重。同结构若出现相互冲突的活性标签，将全部排除并记录数量；其余重复记录保留一条。
- 标签：pIC50 >= 6.3 为强抑制（1），否则为另一类（0），对应 IC50 约 501 nM。此阈值取自 TeachOpenCADD T007，属于研究定义，并非临床阈值。
- 清洗前后行数、两类数量、排除原因均输出为表格。

## 训练、验证和测试

- 拆分前按 Bemis–Murcko 核心骨架分组；同一骨架不得跨组。无环结构如果核心骨架为空，则以规范化 SMILES 作为组名，并记录其数量。
- 固定随机种子 20261108，将组分配至约 70% 训练、15% 验证、15% 测试。记录实际样本比例及每组阳性比例，保存每个分子的分组清单。
- 训练集用于拟合模型。验证集用于比较预先指定的两个指纹方案，或发现明显问题；不根据测试集调整标签、特征、参数或拆分。测试集只评估一次。
- 如果某一组仅含单一类别，则在看任何模型表现之前，按预先列出的后续种子 20261109、20261110……寻找首个三组均含两类的拆分；记录最终种子与原因。

## 比较方案

1. 无分子信息参照：训练集阳性比例，给所有验证/测试分子相同概率。
2. 基础模型：MACCS 指纹 + 随机森林。
3. 改进模型：Morgan radius 2、2048 bit 指纹 + 同样的随机森林。

两种森林使用相同的固定设置：300 棵树、每个叶结点至少 2 个样本、`class_weight='balanced'`、随机种子 20261108。无额外超参数搜索。这样主要比较“结构编码方式”是否带来增益。

## 指标与图表

- 主要指标：测试集 ROC-AUC（按预测概率计算）。
- 次要指标：活性类别的 average precision / PR-AUC、balanced accuracy、灵敏度、特异度；后四项中需阈值的指标固定使用概率 0.5。
- 展示训练、验证、测试表现，避免只报训练分数；在测试集上给出两种指纹模型和无信息参照的并列表格。
- 图表：样本流程与类别数量、ROC/PR 曲线、测试集混淆矩阵。如模型差距很小或结果不佳，如实报告。

## 局限与边界

ChEMBL 活性记录来自不同文献/实验条件；冻结数据集还以每个 ChEMBL 分子的首条记录代表其活性，可能引入测量差异。即使骨架隔离测试良好，也不能据此推断细胞杀伤、药代动力学、安全性、EGFR 突变背景下的活性或患者疗效。单次拆分存在抽样波动；若时间允许，可对预先固定的第二组骨架拆分做敏感性分析，但主结果仍为上述拆分。

## 主要依据

- TeachOpenCADD T001: <https://projects.volkamerlab.org/teachopencadd/talktorials/T001_query_chembl.html>
- TeachOpenCADD T007: <https://projects.volkamerlab.org/teachopencadd/talktorials/T007_compound_activity_machine_learning.html>
- ChEMBL 数据库论文：Zdrazil et al., Nucleic Acids Res. 2024, doi:10.1093/nar/gkad1004。
- Morgan/ECFP：Rogers and Hahn, J Chem Inf Model. 2010, doi:10.1021/ci100050t。
- 骨架拆分：Wu et al., Chem Sci. 2018, doi:10.1039/C7SC02664A。
