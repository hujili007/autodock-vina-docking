import argparse
import json
from pathlib import Path

from docking_pipeline import log_step


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 5: summarize AutoDock Vina docking results."
    )
    parser.add_argument("--summary", default="docking_summary.json", help="Summary JSON from Step 4.")
    return parser.parse_args()


def main():
    args = parse_args()
    summary_path = Path(args.summary)
    if not summary_path.exists():
        raise FileNotFoundError(f"Docking summary file does not exist: {summary_path}")

    log_step("Step 5 summarizes the docking result.")
    data = json.loads(summary_path.read_text(encoding="utf-8"))

    print("Docking result summary")
    print("======================")
    print(f"Receptor: {data.get('receptor')}")
    print(f"Ligand: {data.get('ligand')}")
    print(f"Pose file: {data.get('pose_file')}")
    print(f"Best affinity: {data.get('best_affinity_kcal_per_mol')} kcal/mol")
    print(f"Box mode: {data.get('box_mode')}")
    print(f"Box source: {data.get('box_source')}")
    print(f"Center: {data.get('center')}")
    print(f"Box size: {data.get('box_size')}")


if __name__ == "__main__":
    main()
