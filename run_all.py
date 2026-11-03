#!/usr/bin/env python3
"""
Master Pipeline Runner for AlphaGenome PRNP Analysis.
Executes the full scientific pipeline:
1. Forensic UCSC GRCh38 NM_000311 coordinate and sequence verification
2. AlphaGenome Foundation Model query & continuous track retrieval
3. Dual-Set Publication Figure Generation (Clean & Annotated)
"""

import sys
import subprocess
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

def run_step(step_name, script_path):
    print("\n" + "=" * 70)
    print(f"  RUNNING: {step_name}")
    print("=" * 70)
    cmd = [sys.executable, str(script_path)]
    res = subprocess.run(cmd, cwd=BASE_DIR)
    if res.returncode != 0:
        print(f"❌ Error occurred during: {step_name}")
        sys.exit(res.returncode)
    print(f"✓ Completed: {step_name}")

def main():
    print("==========================================================================")
    print("  AlphaGenome (Google DeepMind, Nature 2026) PRNP Master Runner")
    print("==========================================================================")
    
    run_step("1. Forensic Coordinate & Sequence Verification (UCSC GRCh38 NM_000311)", BASE_DIR / "src" / "verify_prnp_coordinates.py")
    run_step("2. Canonical PRNP Query & Foundation Model Inference (AlphaGenome SDK)", BASE_DIR / "src" / "run_prnp_alphagenome_analysis.py")
    run_step("3. Dual-Set Publication Figure Generation (Clean & Annotated)", BASE_DIR / "src" / "generate_figures_dual_set.py")
    
    print("\n" + "=" * 70)
    print("  🎉 All pipeline steps completed successfully!")
    print("  • Clean Figures:      alphagenome-prnp-analysis/figures/clean/")
    print("  • Annotated Figures:  alphagenome-prnp-analysis/figures/annotated/")
    print("  • Nature Snapshots:   alphagenome-prnp-analysis/figures/paper_reference/")
    print("  • Verified Data:      alphagenome-prnp-analysis/data/")
    print("  • Slide Deck:         alphagenome-prnp-analysis/docs/slides_with_narration.md")
    print("==========================================================================")

if __name__ == "__main__":
    main()
