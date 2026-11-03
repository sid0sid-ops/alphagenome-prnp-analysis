#!/usr/bin/env python3
"""
Convert all citations across slides_with_narration.md to author-year Harvard style,
e.g. (Avsec et al., 2026), (Mead et al., 2015), (Minikel et al., 2020).
"""

import re
from pathlib import Path

SLIDES_PATH = Path("alphagenome-prnp-analysis/docs/slides_with_narration.md")

MAP = {
    1: "Avsec et al., 2026",
    2: "Mead et al., 2015",
    3: "Minikel et al., 2020",
    4: "Avsec et al., 2021",
    5: "Jaganathan et al., 2019",
    6: "Nair et al., 2023",
    7: "Prusiner, 1998",
    8: "Forrest et al., 2014",
    9: "ENCODE Project Consortium, 2020",
    10: "Cheng et al., 2023",
    11: "Kent et al., 2002",
    12: "Vallabh et al., 2021"
}

def replace_citation_group(match):
    nums = [int(n) for n in re.findall(r"\d+", match.group(0))]
    authors = [MAP.get(n, f"Ref {n}") for n in nums]
    joined = "; ".join(authors)
    return f"({joined})"

def convert():
    text = SLIDES_PATH.read_text(encoding="utf-8")
    lines = text.splitlines()
    out = []
    
    in_bib_table = False
    for line in lines:
        if "### 📚 Master Reference Bibliography" in line:
            in_bib_table = True
            out.append(line)
            continue
        if in_bib_table and line.startswith("### 🎙️ Spoken Narration Script"):
            in_bib_table = False
            
        if in_bib_table and line.startswith("| **["):
            # Format: | **[1]** | Avsec, Ž. et al. ...
            m = re.match(r"\|\s*\*\*\[(\d+)\]\*\*\s*\|(.*)", line)
            if m:
                idx = int(m.group(1))
                rest = m.group(2)
                author_tag = MAP.get(idx, f"Ref {idx}")
                out.append(f"| **{author_tag}** |{rest}")
                continue
                
        # Replace bracketed numbers on narration lines
        if line.startswith("> \""):
            line = re.sub(r"\[\d+\](?:,\s*\[\d+\])*", replace_citation_group, line)
        
        # Also clean up visual bullet points
        line = re.sub(r"\[Avsec et al\., Nature 2026\]", "(Avsec et al., 2026)", line)
        line = re.sub(r"\[Avsec et al\., Nat Methods 2021\]", "(Avsec et al., 2021)", line)
        line = re.sub(r"\[Avsec et al\., Nature Methods 2021\]", "(Avsec et al., 2021)", line)
        line = re.sub(r"\[Jaganathan et al\., Cell 2019\]", "(Jaganathan et al., 2019)", line)
        line = re.sub(r"\[Nair et al\., Nat Genet 2023\]", "(Nair et al., 2023)", line)
        line = re.sub(r"\[Zhou, Nat Genet 2022\]", "(Zhou, 2022)", line)
        line = re.sub(r"\[ENCODE Consortium, Nature 2020\]", "(ENCODE Project Consortium, 2020)", line)
        line = re.sub(r"\[ENCODE, Nature 2020\]", "(ENCODE Project Consortium, 2020)", line)
        line = re.sub(r"\[Forrest et al\., Nature 2014\]", "(Forrest et al., 2014)", line)
        line = re.sub(r"\[Cheng et al\., Science 2023\]", "(Cheng et al., 2023)", line)
        line = re.sub(r"\[Kent et al\., Genome Res 2002\]", "(Kent et al., 2002)", line)
        line = re.sub(r"\[Kent et al\., 2002\]", "(Kent et al., 2002)", line)
        line = re.sub(r"\[Büeler et al\., Nature 1992; Prusiner, PNAS 1998\]", "(Büeler et al., 1992; Prusiner, 1998)", line)
        line = re.sub(r"\[Mead, Collinge et al\., Nature 2015; Prusiner, PNAS 1998\]", "(Mead et al., 2015; Prusiner, 1998)", line)
        line = re.sub(r"\[Mead et al\., Nature 2015; Prusiner, PNAS 1998\]", "(Mead et al., 2015; Prusiner, 1998)", line)
        line = re.sub(r"\[Mead et al\., Nature 2015\]", "(Mead et al., 2015)", line)
        line = re.sub(r"\[Mead et al\., 2015\]", "(Mead et al., 2015)", line)
        line = re.sub(r"\[Minikel et al\., Sci Transl Med 2020; Vallabh et al\., Nucleic Acid Ther 2021\]", "(Minikel et al., 2020; Vallabh et al., 2021)", line)
        line = re.sub(r"\[Prusiner, PNAS 1998\]", "(Prusiner, 1998)", line)
        line = re.sub(r"\[DeepMind, 2026\]", "(Google DeepMind, 2026)", line)
        
        out.append(line)
        
    converted_text = "\n".join(out)
    SLIDES_PATH.write_text(converted_text, encoding="utf-8")
    print(f"✓ Successfully converted {SLIDES_PATH} to author-year citations!")

if __name__ == "__main__":
    convert()
