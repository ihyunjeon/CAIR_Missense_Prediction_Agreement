"""
Phase 0 feasibility spike -- part 1: build one fully joined variant table.

ClinVar (genomic coords) -> AlphaMissense hg38 -> AlphaFold (via UniProt).

DESIGN NOTE (learned the hard way, 2026-09-01):
Do NOT join AlphaMissense to MANE on transcript_id. AlphaMissense was built on an
older Ensembl release, so MANE v1.5 versions do not match (BRCA1 .9 vs .8, CFTR
.11 vs .10, MLH1 .8 vs .6) and TP53 is annotated on a different transcript
entirely (ENST00000445888 vs MANE ENST00000269305). An equality join on
transcript_id silently drops ~97% of rows.

(chrom,pos,ref,alt) is the pivot. MANE is used only to RANK transcripts when AM
reports more than one at a locus.

    python scripts/spike/phase0_join.py
"""
from pathlib import Path
import duckdb

GENES = ["TP53", "BRCA1", "CFTR", "PTEN", "MLH1"]
TWO_STAR = ("criteria provided, multiple submitters, no conflicts",
            "reviewed by expert panel", "practice guideline")
PATH   = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
BENIGN = ("Benign", "Likely benign", "Benign/Likely benign")

OUT = Path("data/spike"); OUT.mkdir(parents=True, exist_ok=True)
con = duckdb.connect(); con.execute("PRAGMA threads=4")

print("[1] ClinVar: >=2-star P/LP + B/LB SNVs, GRCh38 ...")
con.execute(f"""
CREATE TABLE cv AS
SELECT "GeneSymbol" AS gene, 'chr' || "Chromosome" AS chrom,
       CAST("PositionVCF" AS BIGINT) AS pos,
       "ReferenceAlleleVCF" AS ref, "AlternateAlleleVCF" AS alt,
       "Name" AS cv_name, "ClinicalSignificance" AS cv_sig,
       "ReviewStatus" AS cv_review, CAST("VariationID" AS BIGINT) AS variation_id,
       CASE WHEN "ClinicalSignificance" IN {PATH} THEN 1 ELSE 0 END AS label_path
FROM read_csv('data/clinvar/variant_summary.txt.gz', delim='\t', header=true,
              quote='', all_varchar=true, ignore_errors=true)
WHERE "Assembly"='GRCh38' AND "Type"='single nucleotide variant'
  AND "ReviewStatus" IN {TWO_STAR}
  AND ("ClinicalSignificance" IN {PATH} OR "ClinicalSignificance" IN {BENIGN})
  AND "GeneSymbol" IN {tuple(GENES)} AND "PositionVCF" IS NOT NULL
""")
n_cv = con.execute("SELECT count(*) FROM cv").fetchone()[0]
print(f"    {n_cv:,} ClinVar rows")

print("[2] MANE Select (for transcript RANKING only) ...")
con.execute("""
CREATE TABLE mane AS
SELECT "symbol" AS gene, split_part("Ensembl_nuc",'.',1) AS enst_base
FROM read_csv('data/mane/mane_summary.txt.gz', delim='\t', header=true,
              all_varchar=true, ignore_errors=true)
WHERE "MANE_status"='MANE Select'
""")

print("[3] coordinate join ClinVar -> AlphaMissense ...")
con.execute("""
CREATE TABLE hits AS
SELECT cv.*, am.uniprot_id, am.transcript_id, am.protein_variant,
       am.am_pathogenicity, am.am_class,
       CASE WHEN m.enst_base IS NOT NULL THEN 1 ELSE 0 END AS is_mane
FROM cv
JOIN 'data/alphamissense/am_hg38.parquet' am
  ON am."CHROM"=cv.chrom AND am."POS"=cv.pos
 AND am."REF"=cv.ref     AND am."ALT"=cv.alt
LEFT JOIN mane m
  ON m.gene = cv.gene AND m.enst_base = split_part(am.transcript_id,'.',1)
""")
print(f"    {con.execute('SELECT count(*) FROM hits').fetchone()[0]:,} raw hits "
      f"(before transcript de-dup)")

# one row per ClinVar variant: prefer the MANE transcript, else lowest ENST
con.execute("""
CREATE TABLE joined AS
SELECT * EXCLUDE (rn) FROM (
  SELECT *, row_number() OVER (
      PARTITION BY variation_id, chrom, pos, ref, alt
      ORDER BY is_mane DESC, transcript_id) AS rn
  FROM hits) WHERE rn = 1
""")
n_j = con.execute("SELECT count(*) FROM joined").fetchone()[0]
print(f"    {n_j:,} rows after de-dup")
print("    per gene:", con.execute(
    "SELECT gene, count(*) FROM joined GROUP BY gene ORDER BY gene").fetchall())
print("    MANE-transcript coverage:", con.execute(
    "SELECT is_mane, count(*) FROM joined GROUP BY 1").fetchall())
print("    label balance:", con.execute(
    "SELECT label_path, count(*) FROM joined GROUP BY 1").fetchall())

con.execute(f"COPY joined TO '{OUT/'joined_pre_assert.parquet'}' (FORMAT parquet)")
print(f"\n[stat] ClinVar rows surviving join: {n_j}/{n_cv} = {100*n_j/n_cv:.1f}%")
print("       (non-missense ClinVar rows cannot match -- AM is missense-only)")
