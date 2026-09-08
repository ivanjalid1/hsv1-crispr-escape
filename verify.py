#!/usr/bin/env python3
"""
verify.py -- confirm that the headline numbers of the preprint reproduce.

Three levels, cheapest first. Run them in order.

    python verify.py                 # LEVEL 1  offline, < 1 s, no downloads
    python verify.py --rerun         # LEVEL 2  offline, ~30 s, needs the genome cache
    python verify.py --rerun --full  # LEVEL 3  offline, ~5 min, needs the genome cache

LEVEL 1 (default) reads the small, version-controlled result files that are the
pinned record of the run reported in the manuscript, and checks every headline
number in the paper against them. It touches the network never and the disk
barely. Use it to answer: *does the archived record actually say what the paper
says?* It also recomputes the Clopper-Pearson resolution floor from first
principles, which needs no files at all.

LEVEL 2 adds recomputation. It re-runs pipeline stages 2-4 (site enumeration,
conservation scoring, ranking) from the cached genomes with --force, then
compares every regenerated file against the SHA-256 recorded in
results/CHECKSUMS.sha256 at pin time. Use it to answer: *does the code still
produce, byte for byte, the files the paper was written from?* It needs the
183 cached FASTAs under data/raw/. If you do not have them yet:

    python run_pipeline.py           # ~80 s; needs NCBI_EMAIL set

LEVEL 3 additionally re-runs stage 6 (the SaCas9 head-to-head, ~50 s) and
stage 7 (the escape model, ~3 min). Both are offline once stage 1 has run.

NOT COVERED, deliberately: stage 8, the human GRCh38 off-target screen. It needs
a 4 GB genome cache (~1.0 GB downloaded) and ~6 minutes of scanning, and no
amount of it is needed to check the conservation, joint-coverage or escape
claims. Its pinned outputs ARE checked at level 1, from the version-controlled
results/offtarget_summary.tsv. To recompute it yourself:

    python src/offtarget.py --stage fetch      # once, explicit, ~4.1 GB on disk
    python src/offtarget.py                    # ~6 min

Exit status is 0 if every check passes and 1 otherwise.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
DATA = ROOT / "data"

# --------------------------------------------------------------------------
# The pinned run. Every number below is quoted in the manuscript or the README;
# the comment says where. Nothing here is computed -- these are the assertions.
# --------------------------------------------------------------------------
PINNED_RETRIEVAL_UTC = "2026-09-07T21:56:25Z"
PINNED_N_GENOMES = 183
PINNED_QUERY = (
    'txid10298[Organism:exp] AND 145000:160000[SLEN] AND biomol_genomic[PROP] '
    'AND "complete genome"[Title] NOT patent[PROP]'
)

# ANSI-free status markers: this has to be readable in a CI log and a Windows console.
OK, BAD, SKIP = "  ok  ", " FAIL ", " skip "


class Checker:
    def __init__(self) -> None:
        self.passed = 0
        self.failed = 0
        self.skipped = 0

    def section(self, title: str) -> None:
        print()
        print(title)
        print("-" * len(title))

    def check(self, label: str, got, want, *, tol: float | None = None) -> bool:
        if tol is not None and isinstance(got, (int, float)) and isinstance(want, (int, float)):
            good = math.isclose(float(got), float(want), rel_tol=tol, abs_tol=0.0)
        else:
            good = got == want
        if good:
            self.passed += 1
            print(f"[{OK}] {label}")
        else:
            self.failed += 1
            print(f"[{BAD}] {label}")
            print(f"          expected: {want!r}")
            print(f"          observed: {got!r}")
        return good

    def skip(self, label: str, why: str) -> None:
        self.skipped += 1
        print(f"[{SKIP}] {label}  ({why})")

    def report(self) -> int:
        print()
        print("=" * 72)
        verdict = "ALL CHECKS PASSED" if self.failed == 0 else f"{self.failed} CHECK(S) FAILED"
        print(f"{verdict}   ({self.passed} passed, {self.failed} failed, {self.skipped} skipped)")
        print("=" * 72)
        return 0 if self.failed == 0 else 1


def load_json(path: Path):
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def load_tsv(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return {row[key]: row for row in csv.DictReader(fh, delimiter="\t")}


def sha256(path: Path) -> str:
    """SHA-256 with CRLF (and lone CR) normalised to LF.

    The pipeline writes CRLF on Windows and LF elsewhere, and git rewrites line
    endings on checkout, so hashing raw bytes would make results/CHECKSUMS.sha256
    platform-specific and would fail spuriously for a reviewer on Linux or macOS.
    The manifest is generated with the same normalisation.
    """
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


# --------------------------------------------------------------------------
# LEVEL 1
# --------------------------------------------------------------------------
def level1(c: Checker) -> None:
    c.section("1. Corpus resolution floor (recomputed from first principles, no files)")
    sys.path.insert(0, str(ROOT))
    try:
        from src.escape import binomial_upper_bound  # noqa: PLC0415
    except Exception as exc:  # pragma: no cover - import guard
        c.skip("Clopper-Pearson upper limit on 0/183", f"cannot import src.escape: {exc}")
    else:
        ub = binomial_upper_bound(0, PINNED_N_GENOMES)
        # Manuscript 2.1: "The exact binomial (Clopper-Pearson) 95% upper limit on
        # 0 events in 183 trials is 0.0162."
        c.check("Clopper-Pearson 95% upper limit on 0/183 = 0.0162", round(ub, 4), 0.0162)
        # Closed form: 1 - alpha**(1/n). Independent of the bisection in src/escape.py.
        c.check(
            "  ... agrees with the closed form 1 - 0.05**(1/183)",
            round(1.0 - 0.05 ** (1.0 / PINNED_N_GENOMES), 6),
            round(ub, 6),
        )
        # Manuscript 2.1: "a single site cannot be certified below ~1.6e-2, a pair
        # below ~2.6e-4".
        c.check("pair floor 0.0162**2 ~ 2.6e-4", round(ub * ub, 6), 0.000262, tol=2e-2)
        # Manuscript 2.1 / Hanley & Lippman-Hand: rule of three, 3/n = 0.0164.
        c.check("rule of three 3/183 ~ 0.0164", round(3.0 / PINNED_N_GENOMES, 4), 0.0164)

    c.section("2. Provenance: what corpus the pinned run used")
    manifest = DATA / "manifest.tsv"
    if manifest.exists():
        rows = list(csv.DictReader(manifest.open(encoding="utf-8", newline=""), delimiter="\t"))
        c.check("data/manifest.tsv holds 183 genomes", len(rows), PINNED_N_GENOMES)
        c.check("manifest records the pinned Entrez query", rows[0]["entrez_query"], PINNED_QUERY)
        c.check("manifest records the pinned retrieval date", rows[0]["retrieval_date_utc"],
                PINNED_RETRIEVAL_UTC)
        lengths = sorted(int(r["length_bp"]) for r in rows)
        # Manuscript 2.1: "183 records, 147,898-159,092 bp".
        c.check("shortest genome 147,898 bp", lengths[0], 147898)
        c.check("longest genome 159,092 bp", lengths[-1], 159092)
    else:
        c.skip("data/manifest.tsv corpus pin", "file not present")

    run_log = RESULTS / "run_log.json"
    if run_log.exists():
        log = load_json(run_log)
        c.check("run_log.json records the pinned query", log["entrez_query"], PINNED_QUERY)
        c.check("run_log.json records the pinned retrieval date",
                log["retrieval_date_utc"], PINNED_RETRIEVAL_UTC)
        c.check("run_log.json lists 183 accessions", len(log["accessions_used"]), PINNED_N_GENOMES)
        c.check("run_log.json accession list is sorted and unique",
                log["accessions_used"], sorted(set(log["accessions_used"])))
        c.check("reference accession NC_001806.2", log["reference_accession"], "NC_001806.2")
    else:
        c.skip("results/run_log.json", "file not present")

    c.section("3. Stage 4: SpCas9 site space and conservation (results/summary.json)")
    summary_path = RESULTS / "summary.json"
    if summary_path.exists():
        s = load_json(summary_path)
        cum = s["conservation_cumulative_counts"]
        # Manuscript 2.1 / README "Results actually observed".
        c.check("183 genomes in the denominator", s["n_strains_total"], PINNED_N_GENOMES)
        c.check("4,777 SpCas9 candidate 23-mers", s["n_guides"], 4777)
        c.check("833 present in all 183 genomes (17.4%)", cum["== 1.00 (perfect)"], 833)
        c.check("  ... which is 17.4% of 4,777", round(100 * 833 / 4777, 1), 17.4)
        c.check("1,465 at >= 0.99", cum[">= 0.99"], 1465)
        c.check("3,431 at >= 0.95", cum[">= 0.95"], 3431)
        c.check("0 sites absent from every genome", cum["== 0.00"], 0)
        c.check("644 perfectly conserved AND filter-passing", s["n_passing_all_filters"], 644)
        c.check("119 carry a poly-T terminator", s["n_has_polyT"], 119)
        c.check("1,366 fall outside the GC band", s["n_gc_out_of_range"], 1366)
        per_gene = s["perfectly_conserved_per_gene"]
        c.check("perfect-by-gene: UL19 247, RL2 91, UL54 52",
                (per_gene["UL19"], per_gene["RL2"], per_gene["UL54"]), (247, 91, 52))
    else:
        c.skip("results/summary.json", "file not present")

    c.section("4. Stage 6: the SaCas9 head-to-head (results/sacas9_benchmark_summary.json)")
    bench_path = RESULTS / "sacas9_benchmark_summary.json"
    if bench_path.exists():
        b = load_json(bench_path)
        g = b["amrani_guides"]
        # Manuscript 2.3, Table 2.
        c.check("ICP0 NNGRRT pool = 46 sites", b["pool_sizes"]["RL2_ICP0"], 46)
        c.check("ICP27 NNGRRT pool = 35 sites", b["pool_sizes"]["UL54_ICP27"], 35)
        c.check("ICP0g2 ranks 31 / 46 in its gene pool", g["ICP0g2"]["rank_in_gene"], 31)
        c.check("ICP0g2 conservation 0.847 over 183 genomes",
                g["ICP0g2"]["conservation_complete_genomes"], 0.846995, tol=1e-5)
        c.check("ICP27g1 ranks 12 / 35", g["ICP27g1"]["rank_in_gene"], 12)
        c.check("ICP27g1 conservation 0.978 over 183 genomes",
                g["ICP27g1"]["conservation_complete_genomes"], 0.978142, tol=1e-5)
        c.check("local GC around ICP0g2 = 0.840", g["ICP0g2"]["local_gc_200bp"], 0.84, tol=1e-6)
        c.check("30 ICP0 sites better conserved than ICP0g2", b["icp0_sites_beating_icp0g2"], 30)
        c.check("  ... 12 of them pass every filter",
                b["icp0_sites_beating_icp0g2_passing_filters"], 12)
        c.check("  ... 4 of them are perfectly conserved",
                b["icp0_sites_beating_icp0g2_perfectly_conserved"], 4)
        # Manuscript 2.2: the joint-vs-marginal result.
        c.check("lead pair jointly intact in 151/183 genomes",
                b["lead_pair"]["n_genomes_both_present"], 151)
        c.check("lead pair joint intactness 0.8251",
                b["lead_pair"]["joint_conservation"], 0.825137, tol=1e-5)
        c.check("  ... 151/183 really is 0.8251", round(151 / 183, 4), 0.8251)
        c.check("  ... below both marginals 0.847 and 0.978",
                b["lead_pair"]["joint_conservation"] < 0.846995, True)
        c.check("  ... and at the lower bound a + b - 1 = 0.825",
                round(0.846995 + 0.978142 - 1, 3), 0.825)
        c.check("  ... product of marginals would predict 0.8285",
                round(0.846995 * 0.978142, 4), 0.8285)
    else:
        c.skip("results/sacas9_benchmark_summary.json", "file not present")

    c.section("5. Stage 7: the multiplex escape model (results/escape_summary.json)")
    esc_path = RESULTS / "escape_summary.json"
    if esc_path.exists():
        e = load_json(esc_path)
        kc = {row["k"]: row for row in e["k_curve"]}
        # Manuscript Table 1, column "P(escape), measured absences".
        c.check("escape model runs over 183 genomes", e["n_genomes"], PINNED_N_GENOMES)
        c.check("449 SaCas9 sites in the 20 nt / NNGRRT pool", e["n_sites"], 449)
        c.check("published lead pair P(escape) = 9.720e-04",
                e["amrani_lead_pair"]["p_escape"], 9.719879026815389e-04, tol=1e-9)
        c.check("published lead pair joint conservation = 0.8251",
                e["amrani_lead_pair"]["joint_conservation"], 0.825136612021858, tol=1e-9)
        c.check("k=1 best P(escape) = 4.754e-03", kc[1]["p_escape_best"], 4.75423268917531e-03,
                tol=1e-9)
        c.check("k=2 best P(escape) = 2.260e-05", kc[2]["p_escape_best"], 2.260272846282309e-05,
                tol=1e-9)
        c.check("k=3 best P(escape) = 1.075e-07", kc[3]["p_escape_best"], 1.074586305225068e-07,
                tol=1e-9)
        c.check("k=4 best P(escape) = 5.109e-10", kc[4]["p_escape_best"], 5.108833339641135e-10,
                tol=1e-9)
        c.check("k=5 best P(escape) = 2.429e-12", kc[5]["p_escape_best"], 2.4288582466870552e-12,
                tol=1e-9)
        # Manuscript 2.1: point estimate needs k=2 at 1e-3 and k=3 at 1e-6.
        c.check("point estimate reaches 1e-3 at k=2", e["minimum_k"]["0.001"], 2)
        c.check("point estimate reaches 1e-6 at k=3", e["minimum_k"]["1e-06"], 3)
        # Manuscript 2.4: the default spectrum's in-frame fraction.
        c.check("default indel spectrum is labelled an ASSUMPTION",
                e["indel_spectrum"]["provenance"], "ASSUMPTION")
        c.check("default in-frame fraction 0.199",
                round(e["indel_spectrum"]["inframe_fraction"], 3), 0.199)
        c.check("baseline repeat regime is 'redundant'", e["params"]["repeat_model"], "redundant")
    else:
        c.skip("results/escape_summary.json", "file not present")

    sens = RESULTS / "escape_sensitivity.tsv"
    if sens.exists():
        rows = list(csv.DictReader(sens.open(encoding="utf-8", newline=""), delimiter="\t"))
        # Abstract: "guide-set rankings hold in 27 of 37 swept settings".
        c.check("sensitivity sweep covers 37 settings", len(rows), 37)
    else:
        c.skip("results/escape_sensitivity.tsv", "file not present")

    c.section("6. Stage 5: what the denominator is worth (results/robustness_summary.json)")
    rob_path = RESULTS / "robustness_summary.json"
    if rob_path.exists():
        r = load_json(rob_path)
        # README "Robustness", sections A and B.
        c.check("effective sample size 132 at 99.9% identity",
                r["redundancy"]["effective_sample_size"], 132)
        c.check("24 multi-member clusters", r["redundancy"]["n_multi_member_clusters"], 24)
        c.check("833 perfect over all 183, 857 over the 132 representatives",
                (r["redundancy"]["perfect_all_genomes"], r["redundancy"]["perfect_deduplicated"]),
                (833, 857))
        c.check("N-tolerant matching moves 833 to 931",
                (r["ambiguity"]["perfect_strict"], r["ambiguity"]["perfect_n_tolerant"]),
                (833, 931))
    else:
        c.skip("results/robustness_summary.json", "file not present")

    c.section("7. Stage 8 + reconciliation: the one-guide substitution (Table 4)")
    rec_path = RESULTS / "recommendation_table.tsv"
    if rec_path.exists():
        rec = load_tsv(rec_path, "guide_id")
        pub, new, withdrawn = rec["RL2_4496+"], rec["RL2_5335+"], rec["RL2_3441+"]
        # Manuscript Table 4, and README "The reconciled recommendation".
        c.check("ICP0g2 is the published lead", pub["published_name"], "ICP0g2")
        c.check("published pair: conservation 0.847", float(pub["cons_183"]), 0.846995, tol=1e-5)
        c.check("published pair: joint 0.8251",
                float(pub["joint_conservation_with_icp27g1"]), 0.825136612021858, tol=1e-9)
        c.check("published pair: P(escape) 9.720e-04",
                float(pub["p_escape_with_icp27g1"]), 9.719879026815349e-04, tol=1e-9)
        c.check("published pair: NNGRRT <=4 mm / <=3 mm = 5 / 0",
                (int(pub["nngrrt_le4"]), int(pub["nngrrt_le3"])), (5, 0))
        c.check("published pair: NNGRRN <=4 mm / <=3 mm = 56 / 5",
                (int(pub["nngrrn_le4"]), int(pub["nngrrn_le3"])), (56, 5))
        c.check("published pair: local GC 0.840", float(pub["local_gc_200bp"]), 0.84, tol=1e-6)

        c.check("RL2_5335+ ranks 1 in the reconciled table", int(new["rank_in_gene"]), 1)
        c.check("RL2_5335+: conservation 1.000 on all 183", float(new["cons_183"]), 1.0)
        c.check("RL2_5335+: conservation 1.000 on the gene-level corpus (429 records)",
                (float(new["cons_gene_level"]), int(new["n_gene_level_records"])), (1.0, 429))
        c.check("RL2_5335+: joint with ICP27g1 = 0.9781 (179/183)",
                (round(float(new["joint_conservation_with_icp27g1"]), 4),
                 int(new["n_genomes_both_present"])), (0.9781, 179))
        c.check("RL2_5335+: P(escape) 2.5145e-04",
                float(new["p_escape_with_icp27g1"]), 2.514535098093381e-04, tol=1e-9)
        c.check("  ... which is 3.87x lower than the published pair",
                round(float(pub["p_escape_with_icp27g1"]) / float(new["p_escape_with_icp27g1"]), 2),
                3.87)
        c.check("RL2_5335+: NNGRRT <=4 mm / <=3 mm = 3 / 0",
                (int(new["nngrrt_le4"]), int(new["nngrrt_le3"])), (3, 0))
        c.check("RL2_5335+: NNGRRN <=4 mm / <=3 mm = 41 / 2",
                (int(new["nngrrn_le4"]), int(new["nngrrn_le3"])), (41, 2))
        c.check("RL2_5335+: 0 coding-exon hits at <=3 mm", int(new["cds_le3"]), 0)
        c.check("RL2_5335+: local GC 0.725", float(new["local_gc_200bp"]), 0.725, tol=1e-6)
        c.check("RL2_5335+: present in both ICP0 repeat copies",
                int(new["n_reference_copies"]), 2)
        # README: the withdrawal. This is the check that must not quietly pass.
        c.check("withdrawn RL2_3441+ is dirtier: NNGRRT <=4 mm = 19",
                int(withdrawn["nngrrt_le4"]), 19)
        c.check("withdrawn RL2_3441+ is dirtier: NNGRRT <=3 mm = 6",
                int(withdrawn["nngrrt_le3"]), 6)
    else:
        c.skip("results/recommendation_table.tsv", "file not present")

    ot_path = RESULTS / "offtarget_summary.tsv"
    if ot_path.exists():
        ot = load_tsv(ot_path, "guide_id")
        # Manuscript 2.5: the screen corroborates Amrani et al. at the canonical PAM.
        c.check("no published lead guide has a NNGRRT site at <=3 mismatches",
                (int(ot["RL2_4496+"]["nngrrt_le3"]), int(ot["UL54_115156+"]["nngrrt_le3"])),
                (0, 0))
        c.check("no published lead guide hits a coding exon at <=3 mismatches",
                (int(ot["RL2_4496+"]["n_le3_in_cds"]), int(ot["UL54_115156+"]["n_le3_in_cds"])),
                (0, 0))
    else:
        c.skip("results/offtarget_summary.tsv", "file not present")


# --------------------------------------------------------------------------
# LEVEL 2 / 3
# --------------------------------------------------------------------------
def run_stage(c: Checker, label: str, args: list[str]) -> bool:
    print(f"\n$ python {' '.join(args)}")
    proc = subprocess.run([sys.executable, *args], cwd=ROOT)
    if proc.returncode != 0:
        c.failed += 1
        print(f"[{BAD}] {label} exited {proc.returncode}")
        return False
    c.passed += 1
    print(f"[{OK}] {label} completed")
    return True


def check_checksums(c: Checker, prefixes: tuple[str, ...]) -> None:
    manifest = RESULTS / "CHECKSUMS.sha256"
    if not manifest.exists():
        c.skip("byte-identity against results/CHECKSUMS.sha256", "manifest not present")
        return
    n = 0
    for line in manifest.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        want, rel = line.split(None, 1)
        rel = rel.strip()
        if not rel.startswith(prefixes):
            continue
        path = ROOT / rel
        n += 1
        if not path.exists():
            c.failed += 1
            print(f"[{BAD}] {rel}: not regenerated")
            continue
        got = sha256(path)
        if got == want:
            c.passed += 1
            print(f"[{OK}] {rel}  sha256 {got[:16]}...")
        else:
            c.failed += 1
            print(f"[{BAD}] {rel}: sha256 differs")
            print(f"          pinned:   {want}")
            print(f"          recomputed: {got}")
    if n == 0:
        c.skip("byte-identity check", "no matching entries in CHECKSUMS.sha256")


def level2(c: Checker, full: bool) -> None:
    c.section("8. Recomputation: stages 2-4 from the cached genomes")
    raw = DATA / "raw"
    n_fasta = len(list(raw.glob("*.fasta"))) if raw.is_dir() else 0
    if n_fasta < PINNED_N_GENOMES:
        c.skip("re-run of stages 2-4", f"only {n_fasta} cached FASTAs under data/raw/; "
                                       "run `python run_pipeline.py` first (needs NCBI_EMAIL)")
        return
    if not run_stage(c, "stages 2-4", ["run_pipeline.py", "--skip-fetch", "--force"]):
        return
    check_checksums(c, ("results/guides_candidates.tsv", "results/conservation.tsv",
                        "results/guides_ranked.tsv", "results/summary.json"))

    if not full:
        print("\n(stages 6 and 7 not recomputed; add --full to include them, ~4 min)")
        return

    c.section("9. Recomputation: stage 6, the SaCas9 head-to-head")
    if run_stage(c, "stage 6", ["run_pipeline.py", "--skip-fetch", "--force",
                                "--benchmark-sacas9"]):
        check_checksums(c, ("results/sacas9_benchmark_pool.tsv",
                            "results/sacas9_benchmark_pairs.tsv"))

    c.section("10. Recomputation: stage 7, the escape model")
    if run_stage(c, "stage 7", ["run_pipeline.py", "--skip-fetch", "--force", "--escape"]):
        check_checksums(c, ("results/escape_k_curve.tsv",
                            "results/escape_k_curve_one_per_gene.tsv",
                            "results/escape_guide_sets.tsv",
                            "results/escape_site_parameters.tsv",
                            "results/escape_codon_tolerance.tsv",
                            "results/escape_sensitivity.tsv",
                            "results/escape_sensitivity_sets.tsv",
                            "results/escape_tolerance_qc.tsv"))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Confirm the headline numbers of the preprint reproduce.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--rerun", action="store_true",
                    help="recompute stages 2-4 and check byte-identity (offline, ~30 s)")
    ap.add_argument("--full", action="store_true",
                    help="with --rerun, also recompute stages 6 and 7 (offline, ~5 min total)")
    args = ap.parse_args()

    print("=" * 72)
    print("verify.py -- headline numbers of the HSV-1 CRISPR conservation analysis")
    print("=" * 72)
    print(f"repository : {ROOT}")
    print(f"python     : {sys.version.split()[0]}")
    for mod in ("numpy", "pandas", "Bio"):
        try:
            m = __import__(mod)
            print(f"{mod:<11}: {getattr(m, '__version__', '?')}")
        except ImportError:
            print(f"{mod:<11}: NOT INSTALLED")

    c = Checker()
    level1(c)
    if args.rerun:
        level2(c, args.full)
    elif args.full:
        print("\n(--full has no effect without --rerun)")
    else:
        print("\n(nothing recomputed; add --rerun to recompute stages 2-4 offline)")
    return c.report()


if __name__ == "__main__":
    raise SystemExit(main())
