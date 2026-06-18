#!/usr/bin/env python3
"""
Golden-file parity comparator for wade-nf (SPEC.md section 14, layer 4).

WADE writes CSVs with quote=FALSE and column/row order that can differ harmlessly
between runs, so we compare as *data*, not bytes:
  * read both CSVs
  * align on a key column (default: first column / SampleNo)
  * sort rows by key and sort columns by name
  * trim surrounding whitespace, normalise NA/empty spellings
  * assert cell-by-cell equality, printing a clear per-cell diff

Exit code 0 = identical (within normalisation); 1 = mismatch; 2 = usage/structure error.

Usage:
    compare_csv.py RESULT.csv EXPECTED.csv [--key SampleNo]
"""
import argparse
import sys
import pandas as pd


def norm(df):
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]
    # stringify + trim + normalise empty/NA spellings
    for c in df.columns:
        df[c] = (df[c].astype(str)
                       .str.strip()
                       .replace({'nan': '', 'NA': '', '<NA>': '', 'None': ''}))
    df = df.reindex(sorted(df.columns), axis=1)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("result")
    ap.add_argument("expected")
    ap.add_argument("--key", default=None,
                    help="key column to sort/align rows on (default: first column)")
    args = ap.parse_args()

    try:
        a = pd.read_csv(args.result, dtype=str, keep_default_na=False)
        b = pd.read_csv(args.expected, dtype=str, keep_default_na=False)
    except Exception as e:
        print(f"ERROR reading inputs: {e}", file=sys.stderr)
        return 2

    a, b = norm(a), norm(b)

    if set(a.columns) != set(b.columns):
        only_a = sorted(set(a.columns) - set(b.columns))
        only_b = sorted(set(b.columns) - set(a.columns))
        print("COLUMN MISMATCH")
        if only_a:
            print(f"  only in result  : {only_a}")
        if only_b:
            print(f"  only in expected: {only_b}")
        return 1

    key = args.key or a.columns[0]
    if key not in a.columns:
        # fall back to first column
        key = a.columns[0]
    a = a.sort_values(key).reset_index(drop=True)
    b = b.sort_values(key).reset_index(drop=True)

    if len(a) != len(b):
        print(f"ROW COUNT MISMATCH: result={len(a)} expected={len(b)}")
        print(f"  result keys  : {sorted(a[key].tolist())}")
        print(f"  expected keys: {sorted(b[key].tolist())}")
        return 1

    diffs = []
    for i in range(len(a)):
        for c in a.columns:
            va, vb = a.at[i, c], b.at[i, c]
            if va != vb:
                diffs.append((a.at[i, key], c, va, vb))

    if diffs:
        print(f"{len(diffs)} CELL MISMATCH(es) (key, column, result, expected):")
        for k, c, va, vb in diffs[:50]:
            print(f"  [{k}] {c}: {va!r} != {vb!r}")
        if len(diffs) > 50:
            print(f"  ... and {len(diffs) - 50} more")
        return 1

    print(f"OK: {len(a)} rows x {len(a.columns)} cols identical (key={key})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
