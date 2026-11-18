#!/usr/bin/env python3
"""Run the four genomes through The RIPper website (https://theripper.hawk.rocks) and save what it returns.

The website is driven through the same addresses its own pages call (upload, RIP profile, LRAR list, window table),
with its default settings: window 1000, step 500, minimum composite 0.01, minimum product 1.1, maximum substrate 0.75,
chain of 7 windows, "change in GC content" off.

IMPORTANT: the website counts letter pairs in UPPERCASE sequence only. The NCBI assemblies (MZ5-1-6, B71, P131) have their
repeats in lowercase (soft-masked), so the website silently ignores about 7 to 14% of each of them and gives wrong index
values in every window that touches a repeat. The genomes are therefore converted to uppercase before upload.

Usage: 09_ripper_website.py --genome NAME=fasta [--genome ...] --raw scratch_dir [--out results/ripper_website]
  --raw   folder for the uppercase copies and the raw JSON; files already there are not fetched again
Needs curl on the PATH.
"""
import argparse, json, subprocess
from pathlib import Path
import pandas as pd

URL = "https://theripper.hawk.rocks"
SET = "1000/500/0.01/1.1/0.75/7"   # window/step/min composite/min product/max substrate/chain
ap = argparse.ArgumentParser()
ap.add_argument("--genome", action="append", required=True)
ap.add_argument("--raw", required=True)
ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "results" / "ripper_website"))
a = ap.parse_args()
raw, out = Path(a.raw), Path(a.out)
raw.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
jar = str(raw / "cookies.txt")   # the website keeps uploaded files per browser session


def curl(*args):
    subprocess.run(["curl", "-s", "-S", "--fail", "-b", jar, "-c", jar, *args], check=True)


def load(path):
    x = json.loads(path.read_text(encoding="utf-8"))
    return json.loads(x) if isinstance(x, str) else x   # the website wraps its JSON in a JSON string


prof, lrar = [], []
for spec in a.genome:
    g, fasta = spec.split("=", 1)
    name = f"{g}.upper.fa"
    want = {"profile": f"api/rip/profile/file/{name}/{SET}/false", "profile_gc": f"api/rip/profile/file/{name}/{SET}/true",
            "lrar": f"api/rip/lrar/file/{name}/{SET}/false", "windows": f"api/rip/file/{name}/1000/500"}
    if not all((raw / f"{g}.{k}.json").exists() for k in want):
        with open(fasta) as i, open(raw / name, "w", newline="\n") as o:
            for l in i:
                o.write(l if l.startswith(">") else l.upper())
        curl("-o", str(raw / f"{g}.upload.txt"), "-F", f"{name}=@{raw / name}", f"{URL}/api/files/upload")
        for k, u in want.items():
            curl("-o", str(raw / f"{g}.{k}.json"), f"{URL}/{u}")
    p, pg = load(raw / f"{g}.profile.json"), load(raw / f"{g}.profile_gc.json")
    w = pd.DataFrame(load(raw / f"{g}.windows.json"))
    w.insert(0, "genome", g)
    w["rip"] = (w.Product >= 1.1) & (w.Substrate <= 0.75) & (w.Composite >= 0.01)   # the rule in the website's own page code
    w.to_csv(out / f"{g}.website_windows.tsv.gz", sep="\t", index=False)
    full = w[w.End - w.Start == 1000]
    prof.append(dict(genome=g, genome_bp=p["FileBP"], gc_percent=p["TotalGCContent"], windows=int(p["WindowsInvestigated"]),
                     rip_affected_windows=p["RIPPositiveWindows"], rip_affected_percent=p["EstimatedGenomeRIP"],
                     lrar_count=p["Count"], lrar_total_bp=p["SumOfLRAR"], lrar_mean_bp=p["SumAverage"], lrar_mean_gc=p["LRARAverageGCContent"],
                     lrar_product=p["ProductAverage"], lrar_substrate=p["SubstrateAverage"], lrar_composite=p["CompositeAverage"],
                     rip_affected_windows_gc_option=pg["RIPPositiveWindows"], rip_affected_percent_gc_option=pg["EstimatedGenomeRIP"],
                     rip_affected_full_length_windows=int(full.rip.sum()), rip_affected_end_of_sequence_windows=int(w.rip.sum() - full.rip.sum())))
    for r in load(raw / f"{g}.lrar.json"):
        lrar.append(dict(genome=g, chrom=r["Name"].split()[0], start=r["Start"], end=r["End"], windows=r["Count"], length=r["Size"],
                         product=r["Product"], substrate=r["Substrate"], composite=r["Composite"], gc_percent=r["GCContent"]))
P = pd.DataFrame(prof)
P.to_csv(out / "ripper_website_profile.tsv", sep="\t", index=False)
pd.DataFrame(lrar, columns=["genome", "chrom", "start", "end", "windows", "length", "product", "substrate", "composite", "gc_percent"]) \
    .to_csv(out / "ripper_website_LRARs.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(P.set_index("genome").T.to_string())
