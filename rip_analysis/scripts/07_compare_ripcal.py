#!/usr/bin/env python3
"""Put the RIPCaL results next to the counts of script 02, group by group (run from manuscripts_all/)."""
import numpy as np
import pandas as pd

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 100)
R = "rip_analysis/results/"
pc = pd.read_csv(R + "ripcal/ripcal_per_copy.tsv", sep="\t")
sm = pd.read_csv(R + "ripcal/ripcal_summary.tsv", sep="\t")
mine = pd.read_csv(R + "four_genomes/family_summary.tsv", sep="\t")
mine["group"] = mine.genome + "." + mine.family + "." + mine.subfamily
rep = pd.read_csv(R + "RIP_summary.tsv", sep="\t")
G = {"E2": "E2", "MZ5-1-6": "MZ5", "B71": "B71", "P131": "PoP131"}
rep["group"] = rep.Strain.map(G) + "." + rep.Family + "." + rep.Group.replace({"whole set": "all"})

num = ["G->A", "A->G", "C->T", "T->C", "Transitions", "Transversions", "Mismatches", "CA->TA", "TG->TA", "CT->TT", "AG->AA",
       "CC->TC", "GG->GA", "CG->TG", "CG->CA"]
a = pc.groupby(["group", "reference"])[num].sum().reset_index()
a["n_copies"] = pc.groupby(["group", "reference"]).size().values
a["rip"] = a["C->T"] + a["G->A"]; a["other_ti"] = a["T->C"] + a["A->G"]
a["ripcal_ratio"] = a.rip / a.other_ti.replace(0, np.nan)
a["ripcal_rip_share"] = a.rip / a.Mismatches.replace(0, np.nan)
a["CpA"] = a["CA->TA"] + a["TG->TA"]; a["CpT"] = a["CT->TT"] + a["AG->AA"]; a["CpC"] = a["CC->TC"] + a["GG->GA"]; a["CpG"] = a["CG->TG"] + a["CG->CA"]
tot = a[["CpA", "CpT", "CpC", "CpG"]].sum(axis=1)
for k in ("CpA", "CpT", "CpC", "CpG"):
    a[k + "_share"] = a[k] / tot
dom = sm[["group", "reference", "dominance ca->ta", "dominance ct->tt", "dominance cc->tc", "dominance cg->tg"]]
a = a.merge(dom, on=["group", "reference"])
m = mine[["group", "copy_C_to_T", "copy_G_to_A", "copy_T_to_C", "copy_A_to_G", "copy_transversions", "copy_substitutions",
          "copy_rip_vs_other_transitions", "event_rip_vs_other_transitions"]]
c = a[a.reference == "consensus"].merge(m, on="group")
c["our_rip"] = c.copy_C_to_T + c.copy_G_to_A
out = a.pivot(index="group", columns="reference", values=["ripcal_ratio", "CpA_share", "CpT_share", "dominance ca->ta", "dominance ct->tt"])
out.columns = [f"{v} ({r})" for v, r in out.columns]
out = out.reset_index().merge(m[["group", "copy_rip_vs_other_transitions", "event_rip_vs_other_transitions"]], on="group")
out = out[out.group.isin(rep.group)].merge(rep[["group", "Strain", "Family", "Group", "Copies"]], on="group")
out.to_csv(R + "ripcal/ripcal_vs_ours.tsv", sep="\t", index=False, float_format="%.3f")
c = c[c.group.isin(rep.group)]
print("copy-level counts, RIPCaL (consensus reference) vs ours:")
print(c[["group", "n_copies", "rip", "our_rip", "other_ti", "copy_T_to_C", "copy_A_to_G", "Transversions", "copy_transversions",
         "ripcal_ratio", "copy_rip_vs_other_transitions"]].round(2).to_string(index=False))
print("correlation of RIP-type counts:", round(np.corrcoef(c.rip, c.our_rip)[0, 1], 4),
      "| of ratios (log):", round(np.corrcoef(np.log(c.ripcal_ratio), np.log(c.copy_rip_vs_other_transitions))[0, 1], 4))
show = out[["Strain", "Family", "Group", "Copies", "copy_rip_vs_other_transitions", "ripcal_ratio (consensus)", "ripcal_ratio (gc)",
            "CpT_share (consensus)", "CpA_share (consensus)", "dominance ct->tt (consensus)", "dominance ca->ta (consensus)"]]
print(show.round(2).to_string(index=False))
print("dominant site by RIPCaL (consensus):", (out[["dominance ca->ta (consensus)", "dominance ct->tt (consensus)"]].idxmax(axis=1)).value_counts().to_dict())
print("median shares consensus: CpT", out["CpT_share (consensus)"].median().round(3), "CpA", out["CpA_share (consensus)"].median().round(3))
