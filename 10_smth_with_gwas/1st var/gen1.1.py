import re
import os
import numpy as np
import pandas as pd

# ------------------ helpers ------------------

def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

_float_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

def parse_first_float(x):
    """Extract first numeric token from a cell; return NaN if missing."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan
    s = str(x).strip()
    if s in ("", "-", "NA", "NaN", "nan", "None"):
        return np.nan
    m = _float_re.search(s)
    return float(m.group(0)) if m else np.nan

def parse_beta_signed(x):
    """Parse beta with 'increase/decrease' words (decrease -> negative)."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan
    s = str(x).strip().lower()
    if s in ("", "-", "na", "nan", "none"):
        return np.nan
    m = _float_re.search(s)
    if not m:
        return np.nan
    v = float(m.group(0))
    # GWAS export often encodes direction in text
    if "decrease" in s and v > 0:
        v = -v
    return v

def extract_rsid(x):
    """Return rsID (rs123) from a variant string, else None."""
    if x is None:
        return None
    m = re.search(r"(rs\d+)", str(x), flags=re.IGNORECASE)
    return m.group(1).lower() if m else None

# ------------------ main logic ------------------

def build_ncbi_query_sort_by_raf(file_path, top_n=100):
    if not os.path.exists(file_path):
        raise FileNotFoundError("Файл не найден. Проверь путь.")

    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    # Columns (GWAS Catalog export typical names)
    var_col = pick_col(df, ["riskAllele", "SNPS", "Variant and risk allele", "VARIANT_ID", "variantId"])
    or_col = pick_col(df, ["orValue", "OR", "odds_ratio"])
    beta_col = pick_col(df, ["beta", "BETA", "Beta"])
    raf_col = pick_col(df, ["riskFrequency", "RAF", "raf", "risk_allele_frequency"])

    if not var_col:
        raise ValueError(f"Не найдена колонка с вариантами. Колонки: {list(df.columns)}")
    if not raf_col:
        raise ValueError(f"Не найдена колонка RAF (riskFrequency/RAF). Колонки: {list(df.columns)}")
    if not or_col and not beta_col:
        raise ValueError(f"Не найдены колонки OR/beta. Колонки: {list(df.columns)}")

    # Parse RAF
    df["raf_num"] = df[raf_col].apply(parse_first_float)

    # Parse effects (for filtering)
    df["or_num"] = df[or_col].apply(parse_first_float) if or_col else np.nan
    df["beta_num"] = df[beta_col].apply(parse_beta_signed) if beta_col else np.nan

    # Filter: keep rows where OR or beta exists
    df_f = df[(df["or_num"].notna()) | (df["beta_num"].notna())].copy()

    # Keep rows with RAF
    df_f = df_f[df_f["raf_num"].notna()].copy()

    # Extract rsID
    df_f["rsid"] = df_f[var_col].apply(extract_rsid)

    # Sort by RAF (descending)
    df_sorted = df_f.sort_values("raf_num", ascending=False)

    # Unique rsIDs in that order
    rsids = (
        df_sorted["rsid"]
        .dropna()
        .drop_duplicates()
        .head(top_n)
        .tolist()
    )

    query = " OR ".join(rsids)

    # Some useful output table (top rows by RAF)
    out_table = df_sorted.head(max(top_n, 100)).copy()

    # For logging: which effect type is more prevalent
    or_count = int(df_f["or_num"].notna().sum())
    beta_count = int(df_f["beta_num"].notna().sum())
    effect_used = "OR" if or_count >= beta_count else "BETA"

    return effect_used, rsids, query, out_table

# ------------------ run (interactive) ------------------

if __name__ == "__main__":
    file_path = input("Введите имя/путь к TSV файлу: ").strip()
    top_n_in = input("Сколько rsID взять (по умолчанию 100): ").strip()
    top_n = int(top_n_in) if top_n_in else 100

    effect_used, rsids, query, out_table = build_ncbi_query_sort_by_raf(file_path, top_n=top_n)

    print("\nФильтр: строки с OR или beta (эффект чаще представлен как):", effect_used)
    print("Уникальных rsID найдено:", len(rsids))
    print("\nNCBI query:\n")
    print(query)

    # Save outputs
    with open("ncbi_query.txt", "w", encoding="utf-8") as f:
        f.write(query)

    out_table.to_csv("top_variants_sorted_by_RAF.csv", index=False, encoding="utf-8")

    print("\nСохранено:")
    print(" - ncbi_query.txt")
    print(" - top_variants_sorted_by_RAF.csv")