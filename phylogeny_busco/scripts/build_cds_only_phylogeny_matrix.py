#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def write_fasta(seqs: dict[str, str], path: Path, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in seqs.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def write_phylip(seqs: dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    nseq = len(seqs)
    nchar = len(next(iter(seqs.values()))) if seqs else 0
    with path.open("w") as handle:
        handle.write(f"{nseq} {nchar}\n")
        for name, seq in seqs.items():
            handle.write(f"{name} {seq}\n")


def attrs_to_dict(attrs: str) -> dict[str, str]:
    out = {}
    for part in attrs.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            out[key] = value
    return out


def fetch_remote_text(host: str, remote_path: str, attempts: int = 8) -> str:
    cmd = [
        "ssh",
        "-o",
        "BatchMode=yes",
        "-o",
        "ConnectTimeout=12",
        host,
        f"cat {remote_path}",
    ]
    last = ""
    for attempt in range(1, attempts + 1):
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if proc.returncode == 0:
            return proc.stdout
        last = proc.stderr.strip()
        if attempt < attempts:
            time.sleep(30)
    raise RuntimeError(f"Could not read {remote_path} from {host}: {last}")


def map_interval_to_tx(feature: dict, interval: tuple[int, int]) -> list[tuple[int, int]]:
    exons = feature["exon"]
    strand = feature["strand"]
    ordered = sorted(exons, key=lambda item: item[0], reverse=(strand == "-"))
    tx_pos = 1
    out: list[tuple[int, int]] = []
    for ex_start, ex_end in ordered:
        ex_len = ex_end - ex_start + 1
        ov_start = max(ex_start, interval[0])
        ov_end = min(ex_end, interval[1])
        if ov_start <= ov_end:
            if strand == "+":
                tx_start = tx_pos + (ov_start - ex_start)
                tx_end = tx_pos + (ov_end - ex_start)
            else:
                tx_start = tx_pos + (ex_end - ov_end)
                tx_end = tx_pos + (ex_end - ov_start)
            out.append((tx_start, tx_end))
        tx_pos += ex_len
    return out


def classify_sites(sites: list[dict[str, str]], gff_text: str) -> list[str]:
    needed_genes = {row["chrom"] for row in sites}
    needed_mrnas = {f"{gene}-mRNA-1" for gene in needed_genes}
    features = defaultdict(lambda: {"strand": None, "exon": [], "CDS": [], "five_prime_UTR": [], "three_prime_UTR": []})

    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        typ = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attrs = attrs_to_dict(fields[8])
        if typ == "mRNA":
            mid = attrs.get("ID")
            if mid in needed_mrnas:
                features[mid]["strand"] = strand
        elif typ in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
            parent = attrs.get("Parent")
            if parent in needed_mrnas:
                features[parent]["strand"] = strand
                features[parent][typ].append((start, end))

    tx_features: dict[str, list[tuple[int, int, str]]] = {}
    missing_exons = []
    for gene in needed_genes:
        mrna = f"{gene}-mRNA-1"
        feature = features[mrna]
        if not feature["exon"]:
            missing_exons.append(gene)
            tx_features[gene] = []
            continue
        intervals: list[tuple[int, int, str]] = []
        for typ in ("five_prime_UTR", "CDS", "three_prime_UTR"):
            for interval in feature[typ]:
                for tx_interval in map_interval_to_tx(feature, interval):
                    intervals.append((tx_interval[0], tx_interval[1], typ))
        tx_features[gene] = intervals
    if missing_exons:
        raise RuntimeError(f"Missing exon annotations for {len(missing_exons)} genes; first: {missing_exons[:5]}")

    classes = []
    for row in sites:
        gene = row["chrom"]
        pos = int(row["pos"])
        labels = [typ for start, end, typ in tx_features.get(gene, []) if start <= pos <= end]
        if "CDS" in labels:
            classes.append("CDS")
        elif "five_prime_UTR" in labels:
            classes.append("five_prime_UTR")
        elif "three_prime_UTR" in labels:
            classes.append("three_prime_UTR")
        else:
            classes.append("unannotated_transcript_non_CDS")
    return classes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--input-sites", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--ssh-host", default="user@cluster")
    parser.add_argument(
        "--remote-gff3",
        default="/path/to/scratch/finger_millet_blast_phylogeny/E2-annotation-Final/E2.gff3",
    )
    args = parser.parse_args()
    args.input_fasta = args.input_fasta.resolve()
    args.input_sites = args.input_sites.resolve()
    args.outdir = args.outdir.resolve()

    seqs = read_fasta(args.input_fasta)
    sites: list[dict[str, str]] = []
    with args.input_sites.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sites.append(row)
    if not seqs:
        raise SystemExit("No sequences found")
    if any(len(seq) != len(sites) for seq in seqs.values()):
        raise SystemExit("FASTA sequence lengths do not match sites table")

    gff_text = fetch_remote_text(args.ssh_host, args.remote_gff3)
    classes = classify_sites(sites, gff_text)
    counts = Counter(classes)
    keep = [i for i, cls in enumerate(classes) if cls == "CDS"]

    filtered = {name: "".join(seq[i] for i in keep) for name, seq in seqs.items()}
    kept_sites = [sites[i] for i in keep]

    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_out = args.outdir / f"{args.prefix}.fasta"
    phy_out = args.outdir / f"{args.prefix}.phy"
    sites_out = args.outdir / f"{args.prefix}.sites.tsv"
    summary_out = args.outdir / f"{args.prefix}.summary.tsv"

    write_fasta(filtered, fasta_out)
    write_phylip(filtered, phy_out)
    with sites_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(sites[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept_sites)
    with summary_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["input_samples", len(seqs)])
        writer.writerow(["input_sites", len(sites)])
        writer.writerow(["kept_CDS_sites", len(keep)])
        writer.writerow(["dropped_non_CDS_sites", len(sites) - len(keep)])
        for key in ("CDS", "five_prime_UTR", "three_prime_UTR", "unannotated_transcript_non_CDS"):
            writer.writerow([f"classified_{key}_sites", counts.get(key, 0)])
        writer.writerow(["fasta", str(fasta_out)])
        writer.writerow(["phylip", str(phy_out)])
        writer.writerow(["sites", str(sites_out)])

    print(fasta_out.resolve())
    print(phy_out.resolve())
    print(summary_out.resolve())


if __name__ == "__main__":
    main()
