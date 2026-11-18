from pathlib import Path


fasta = Path(__file__).resolve().parent / "Fosbury_LTR_lt99_full_length_sequences.fasta"

record_id = None
chunks = []
for line in fasta.read_text(encoding="utf-8").splitlines():
    if line.startswith(">"):
        if record_id is not None:
            seq = "".join(chunks).upper()
            print(record_id, len(seq), seq[:5], seq[-5:], seq[:5] == seq[-5:], seq[5:7], seq[-7:-5])
        record_id = line[1:].split("|", 1)[0]
        chunks = []
    elif line.strip():
        chunks.append(line.strip())

if record_id is not None:
    seq = "".join(chunks).upper()
    print(record_id, len(seq), seq[:5], seq[-5:], seq[:5] == seq[-5:], seq[5:7], seq[-7:-5])
