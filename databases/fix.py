import re
import pandas as pd

def split_miRNA_target(in_path, out_path, out_target_col):
    # Read raw single-column file
    # auto-detect separator, but will fallback to one-column if none
    with open(in_path, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
    # First line is header
    header = lines[0]
    rows = lines[1:]

    # Prepare output rows
    out = []
    # miRNA ID pattern: hsa-miR-... or hsa-let-...
    # We’ll split at the boundary between miRNA and the rest
    # This regex captures miRNA (miR or let family) with optional arm -5p/-3p, then the rest
    pat = re.compile(r'^(hsa-(?:miR|let)-[0-9a-zA-Z\-]+?p?)(.+)$')

    for r in rows:
        m = pat.match(r)
        if not m:
            # If not matching, try a more permissive pattern (up to first uppercase letter in gene)
            # or skip/log
            # Attempt split at first uppercase letter after miRNA prefix
            m2 = re.match(r'^(hsa-(?:miR|let)-[0-9a-zA-Z\-]+?p?)([A-Z].+)$', r)
            if not m2:
                # If still not matched, skip row or log
                # print(f"Skipping unmapped row: {r}")
                continue
            mi, tgt = m2.group(1), m2.group(2)
        else:
            mi, tgt = m.group(1), m.group(2)

        mi = mi.strip()
        tgt = tgt.strip()
        out.append((mi, tgt))

    if not out:
        raise ValueError(f"No rows could be parsed from {in_path}. Check formatting.")

    df = pd.DataFrame(out, columns=["miRNA", out_target_col])
    # Drop duplicates, strip spaces just in case
    df["miRNA"] = df["miRNA"].str.strip()
    df[out_target_col] = df[out_target_col].str.strip()
    df = df.drop_duplicates().reset_index(drop=True)
    # Write as tab-delimited TSV
    df.to_csv(out_path, sep='\t', index=False)
    print(f"Wrote {len(df)} rows to {out_path}")

# Paths: adjust to your actual file names/locations
split_miRNA_target("LncBase.txt", "LncBase_w.txt", out_target_col="lncRNA")
split_miRNA_target("miRTarBase.txt", "miRTarBase_w.txt", out_target_col="mRNA")