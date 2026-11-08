"""Reproduce the preregistered EGFR structure/activity comparison."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MPL_CACHE = ROOT / ".cache" / "matplotlib"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdkit
import sklearn
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    precision_recall_curve,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GroupShuffleSplit


RAW = ROOT / "data" / "raw" / "EGFR_compounds_ea055ef.csv"
OUT = ROOT / "results"
SEED = 20261108
LABEL_THRESHOLD = 6.3


def load_and_clean() -> tuple[pd.DataFrame, dict]:
    with RAW.open("rb") as file:
        digest = hashlib.sha256(file.read()).hexdigest()
    data = pd.read_csv(RAW)
    audit: dict = {
        "source": str(RAW.relative_to(ROOT)),
        "sha256": digest,
        "raw_rows": len(data),
        "excluded": Counter(),
    }
    records = []
    for row in data.itertuples(index=False):
        try:
            ic50 = float(row.IC50)
            pic50 = float(row.pIC50)
        except (ValueError, TypeError):
            audit["excluded"]["missing_or_non_numeric_activity"] += 1
            continue
        if not (math.isfinite(ic50) and ic50 > 0 and math.isfinite(pic50)):
            audit["excluded"]["non_positive_or_non_finite_activity"] += 1
            continue
        if abs(pic50 - (9 - math.log10(ic50))) > 1e-6:
            raise ValueError(f"IC50/pIC50 disagreement: {row.molecule_chembl_id}")
        mol = Chem.MolFromSmiles(str(row.smiles))
        if mol is None:
            audit["excluded"]["invalid_smiles"] += 1
            continue
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold = Chem.MolToSmiles(scaffold_mol, isomericSmiles=False)
        records.append(
            {
                "molecule_chembl_id": row.molecule_chembl_id,
                "IC50_nM": ic50,
                "pIC50": pic50,
                "smiles": canonical,
                "scaffold_group": scaffold or f"ACYCLIC:{canonical}",
                "label": int(pic50 >= LABEL_THRESHOLD),
                "is_acyclic": not bool(scaffold),
            }
        )

    by_smiles = defaultdict(list)
    for record in records:
        by_smiles[record["smiles"]].append(record)
    cleaned = []
    for rows in by_smiles.values():
        if len({row["label"] for row in rows}) > 1:
            audit["excluded"]["conflicting_duplicate_structure"] += len(rows)
        else:
            cleaned.append(rows[0])
            audit["excluded"]["same_label_duplicate_structure"] += len(rows) - 1
    frame = pd.DataFrame(cleaned).sort_values("molecule_chembl_id").reset_index(drop=True)
    audit["excluded"] = dict(audit["excluded"])
    audit["clean_rows"] = len(frame)
    audit["class_counts"] = {str(k): int(v) for k, v in frame.label.value_counts().sort_index().items()}
    audit["unique_scaffolds"] = frame.scaffold_group.nunique()
    audit["acyclic_rows"] = int(frame.is_acyclic.sum())
    return frame, audit


def split_data(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    groups = frame.scaffold_group.to_numpy()
    for seed in range(SEED, SEED + 21):
        train_idx, remaining_idx = next(
            GroupShuffleSplit(n_splits=1, train_size=0.70, random_state=seed).split(frame, frame.label, groups)
        )
        remaining = frame.iloc[remaining_idx]
        val_rel, test_rel = next(
            GroupShuffleSplit(n_splits=1, train_size=0.50, random_state=seed).split(
                remaining, remaining.label, remaining.scaffold_group
            )
        )
        val_idx, test_idx = remaining_idx[val_rel], remaining_idx[test_rel]
        if all(frame.iloc[idx].label.nunique() == 2 for idx in (train_idx, val_idx, test_idx)):
            break
    else:
        raise RuntimeError("No valid two-class split within preregistered seeds")
    result = frame.copy()
    result["split"] = ""
    result.loc[train_idx, "split"] = "train"
    result.loc[val_idx, "split"] = "validation"
    result.loc[test_idx, "split"] = "test"
    detail = {
        "split_seed": seed,
        "fallback_seeds_tried": seed - SEED,
        "splits": {},
        "scaffold_overlap": 0,
    }
    for name in ("train", "validation", "test"):
        part = result[result.split == name]
        detail["splits"][name] = {
            "n": len(part),
            "fraction": len(part) / len(result),
            "positive": int(part.label.sum()),
            "positive_fraction": float(part.label.mean()),
            "unique_scaffolds": part.scaffold_group.nunique(),
        }
    assert result.groupby("scaffold_group").split.nunique().max() == 1
    return result, detail


def fingerprints(smiles: pd.Series, kind: str) -> np.ndarray:
    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    output = []
    for value in smiles:
        mol = Chem.MolFromSmiles(value)
        fp = MACCSkeys.GenMACCSKeys(mol) if kind == "MACCS" else morgan.GetFingerprint(mol)
        arr = np.zeros(len(fp), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        output.append(arr)
    return np.stack(output)


def evaluate(y: np.ndarray, prob: np.ndarray) -> dict:
    predicted = (prob >= 0.5).astype(int)
    cm = confusion_matrix(y, predicted, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    return {
        "n": int(len(y)),
        "positive": int(y.sum()),
        "roc_auc": float(roc_auc_score(y, prob)),
        "average_precision": float(average_precision_score(y, prob)),
        "balanced_accuracy": float(balanced_accuracy_score(y, predicted)),
        "sensitivity": float(recall_score(y, predicted, pos_label=1, zero_division=0)),
        "specificity": float(tn / (tn + fp)),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def draw_curves(y: np.ndarray, predictions: dict[str, np.ndarray]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for name, probability in predictions.items():
        fpr, tpr, _ = roc_curve(y, probability)
        precision, recall, _ = precision_recall_curve(y, probability)
        axes[0].plot(fpr, tpr, label=f"{name} (AUC {roc_auc_score(y, probability):.3f})")
        axes[1].plot(recall, precision, label=f"{name} (AP {average_precision_score(y, probability):.3f})")
    axes[0].set(xlabel="False positive rate", ylabel="Sensitivity", title="Test ROC")
    axes[1].set(xlabel="Sensitivity", ylabel="Precision", title="Test precision-recall")
    for axis in axes:
        axis.legend(loc="best", fontsize=8)
        axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(OUT / "test_roc_pr.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    frame, audit = load_and_clean()
    frame, split = split_data(frame)
    frame.to_csv(OUT / "split_manifest.csv", index=False)
    (OUT / "data_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (OUT / "split_audit.json").write_text(json.dumps(split, indent=2), encoding="utf-8")

    masks = {name: (frame.split == name).to_numpy() for name in ("train", "validation", "test")}
    y = frame.label.to_numpy()
    baseline = float(y[masks["train"]].mean())
    predictions = {name: {} for name in masks}
    metrics = []
    for name in masks:
        prob = np.full(masks[name].sum(), baseline)
        predictions[name]["Prevalence baseline"] = prob
        metrics.append({"split": name, "model": "Prevalence baseline", **evaluate(y[masks[name]], prob)})

    for kind in ("MACCS", "Morgan"):
        x = fingerprints(frame.smiles, kind)
        model = RandomForestClassifier(
            n_estimators=300,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=SEED,
            n_jobs=-1,
        )
        model.fit(x[masks["train"]], y[masks["train"]])
        for name in masks:
            prob = model.predict_proba(x[masks[name]])[:, 1]
            predictions[name][kind] = prob
            metrics.append({"split": name, "model": kind, **evaluate(y[masks[name]], prob)})
        print(f"Finished {kind}", flush=True)

    scores = pd.DataFrame(metrics)
    scores.to_csv(OUT / "metrics.csv", index=False)
    val = scores[scores.split == "validation"].set_index("model")
    selected = "Morgan" if val.loc["Morgan", "roc_auc"] > val.loc["MACCS", "roc_auc"] else "MACCS"
    (OUT / "selection.json").write_text(
        json.dumps({"criterion": "validation ROC-AUC", "selected_model": selected}, indent=2), encoding="utf-8"
    )
    test_predictions = pd.DataFrame({
        "molecule_chembl_id": frame.loc[masks["test"], "molecule_chembl_id"].to_numpy(),
        "label": y[masks["test"]],
        **{name + "_probability": prob for name, prob in predictions["test"].items()},
    })
    test_predictions.to_csv(OUT / "test_predictions.csv", index=False)
    draw_curves(y[masks["test"]], predictions["test"])
    for name in ("MACCS", "Morgan"):
        cm = confusion_matrix(y[masks["test"]], predictions["test"][name] >= 0.5, labels=[0, 1])
        disp = ConfusionMatrixDisplay(cm, display_labels=["Other", "Strong"])
        disp.plot(cmap="Blues", colorbar=False)
        plt.title(f"{name}: test confusion matrix")
        plt.tight_layout()
        plt.savefig(OUT / f"test_confusion_{name.lower()}.png", dpi=180)
        plt.close()
    versions = {
        "python": platform.python_version(), "rdkit": rdkit.__version__, "sklearn": sklearn.__version__,
        "pandas": pd.__version__, "numpy": np.__version__,
    }
    (OUT / "versions.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")
    print(scores.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"Selected on validation: {selected}")


if __name__ == "__main__":
    main()
