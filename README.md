# EGFR in vitro inhibition prediction

**AIDD8122 final project — CEN FENG (P2600163)**

This project asks whether a molecule's two-dimensional structure can predict **strong inhibition of human EGFR in an in vitro binding assay**. The label is `pIC50 >= 6.3`, approximately `IC50 <= 501 nM`. It is a research definition, not a clinical treatment threshold. The model does not predict cancer patient outcomes.

## Main result

On one fixed Bemis–Murcko scaffold split, a Morgan fingerprint plus random forest achieved test ROC-AUC **0.865**, compared with **0.831** for MACCS plus the same forest and **0.500** for a no-structure prevalence baseline. The Morgan model was selected on validation ROC-AUC before the test comparison. This is a single-split result; the test set contains a large, relatively active scaffold group and has a higher positive fraction than training. See [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md) for the complete metrics and limitations.

## Files

| Path | Purpose |
|---|---|
| `PROJECT_PROTOCOL.md` | Analysis rules recorded before model fitting |
| `run_analysis.py` | Full data audit, scaffold split, model training and evaluation |
| `requirements.txt` | Python dependencies |
| `data/raw/EGFR_compounds_ea055ef.csv` | Unmodified frozen public input file |
| `DATA_ATTRIBUTION.md` | Data origin and reuse information |
| `results/` | Audits, split manifest, metrics, predictions and figures |

## Reproduce the analysis

Use Python 3.12 in a fresh environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_analysis.py
```

The script writes `results/` from the original CSV. The supplied run used Python 3.12.14, RDKit 2026.09.1, scikit-learn 1.9.1, pandas 3.1.1 and NumPy 2.5.4. The data checksum and cleaning counts are in `results/data_audit.json`; each molecule's fixed split is in `results/split_manifest.csv`.

## Scope and limitations

The dataset was prepared from ChEMBL 27 through TeachOpenCADD T001. It contains exact IC50 measurements for the human EGFR single-protein target, but the frozen compound CSV does not retain detailed assay metadata. Different experimental conditions and the first-measurement selection in T001 may introduce label noise. There is no evaluation of cells, EGFR mutations, kinase selectivity, toxicity, pharmacokinetics, or patient benefit. A single scaffold split is insufficient to establish a stable performance difference.

Generative AI assisted with code, manuscript and slide drafting. Reported numbers were produced by running the included script on public data; no experimental or patient data were fabricated.

## Sources

- [TeachOpenCADD T001: compound data acquisition](https://projects.volkamerlab.org/teachopencadd/talktorials/T001_query_chembl.html)
- [TeachOpenCADD T007: ligand-based screening](https://projects.volkamerlab.org/teachopencadd/talktorials/T007_compound_activity_machine_learning.html)
- Zdrazil B et al. *Nucleic Acids Research* 2024. doi:10.1093/nar/gkad1004.
- Rogers D, Hahn M. *Journal of Chemical Information and Modeling* 2010. doi:10.1021/ci100050t.
- Wu Z et al. *Chemical Science* 2018. doi:10.1039/C7SC02664A.

## 中文提示

这项研究预测的是**分子在体外实验中对 EGFR 的抑制活性**，不预测患者是否从 EGFR 靶向治疗获益。测试集 ROC-AUC 为 Morgan 0.865、MACCS 0.831；测试集的样本构成与训练集不同，报告中已如实说明。
