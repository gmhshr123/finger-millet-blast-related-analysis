from pathlib import Path


base = Path(__file__).resolve().parent
files = [
    "Fosbury_LTR_lt99_corrected_element_sequences_TG_to_CA.fasta",
    "Fosbury_LTR_lt99_corrected_element_sequences_with_5bp_TSDs.fasta",
]


def records(path):
    header = None
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks).upper()
            header = line
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks).upper()


for name in files:
    count = 0
    ok = 0
    for _, sequence in records(base / name):
        count += 1
        if "with_5bp_TSDs" in name:
            ok += sequence[:5] == sequence[-5:] and sequence[5:7] == "TG" and sequence[-7:-5] == "CA"
        else:
            ok += sequence[:2] == "TG" and sequence[-2:] == "CA"
    print(name, count, ok)
