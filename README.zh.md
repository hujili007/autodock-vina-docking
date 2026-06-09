# AutoDock Vina 分子对接 Skill

AutoDock Vina 分子对接 Skill 是一个可复现的脚本化蛋白-小分子对接流程。它适合学习分子对接、小规模对接实验，以及在 Codex 或 Claude 中作为 agent skill 自动调用。

该流程接收蛋白 target/PDB ID 和小分子名称/PubChem CID，自动准备受体和配体文件、判断 docking box、运行 AutoDock Vina、提取 best pose、生成 PyMOL 图 a/图 b，并写入结果汇总。

## 功能特性

- 支持从命令行参数运行单个蛋白-小分子对接。
- 支持 CSV 单行或多行批量对接。
- 批量模式下，相同 `target+pdbid` 的受体和 docking box 只准备一次并复用。
- 根据 PubChem CID 下载小分子结构。
- 优先使用 PubChem 3D SDF；如果没有 3D SDF，则使用 PubChem 2D SDF 并用 RDKit 生成 3D 构象。
- 生成受体和配体 PDBQT 文件。
- 通过 AutoDock Vina Python API 运行对接。
- 将 best pose 提取为 PDBQT 和 PDB。
- 自动生成 PyMOL 风格的图 a 和图 b。
- 写入 `docking_results_summary.xlsx`、每次运行的 `manifest.json` 和批量汇总文件。

## 目录结构

```text
autodock-vina-docking/
├── SKILL.md
├── README.md
├── README.zh.md
├── references/
│   ├── docking_box.md
│   ├── pubchem_ligand.md
│   ├── pymol_figures.md
│   └── workflow.md
└── scripts/
    ├── setup_conda_env.sh
    ├── run_pipeline.py
    ├── run_batch.py
    ├── PrepareReceptor.py
    ├── IdentifyBindingSite.py
    ├── PrepareLigand.py
    ├── RunDocking.py
    ├── AnalyzeDocking.py
    ├── render_pymol_figures.py
    ├── docking_pipeline.py
    └── test_docking_pipeline.py
```

## 环境要求

请先安装 Miniconda、Anaconda 或 Mambaforge。然后在项目根目录创建推荐的 `docking` 环境：

```bash
bash scripts/setup_conda_env.sh
```

等价的手动安装命令：

```bash
conda create -n docking -c conda-forge python=3.11 meeko openbabel vina rdkit prody pymol-open-source pillow
```

如果已经有 `docking` 环境，只是缺少 PyMOL 出图能力：

```bash
conda install -n docking -c conda-forge pymol-open-source pillow
```

验证环境：

```bash
conda run -n docking python scripts/run_pipeline.py --help
conda run -n docking python scripts/run_batch.py --help
```

## 单个小分子对接示例

在 `autodock-vina-docking/` 目录下运行：

```bash
conda run -n docking python scripts/run_pipeline.py \
  --target S100A8 \
  --pdb-id 5HLO \
  --compound-name Epicatechin \
  --pubchem-cid 72276
```

运行后会创建：

```text
runs/<PDBID>_<CID>_<YYYYMMDD_HHMMSS>/
```

项目级 Excel 总表会被创建或追加：

```text
docking_results_summary.xlsx
```

## CSV 批量对接示例

创建 `ligands.csv`，表头必须严格使用：

```csv
compound_name,pubchem_cid,target,pdbid
Epicatechin,72276,S100A8,5HLO
DL-Tryptophan,1148,S100A8,5HLO
```

运行：

```bash
conda run -n docking python scripts/run_batch.py --input-csv ligands.csv
```

批量模式会创建：

```text
runs/batch_<YYYYMMDD_HHMMSS>/
```

如果 CSV 中多行使用相同的 `target+pdbid`，受体准备和 docking box 判断只执行一次并复用。每一行小分子仍然会有独立的 docking 目录。

## 输出文件

单次运行和批量模式中的每一行会输出：

- 原始受体 PDB。
- protein-only 受体 PDB。
- 受体 PDBQT。
- 配体 SDF。
- 配体 PDBQT。
- `vina_box_config.txt`。
- 多 pose 对接结果 PDBQT。
- best pose PDBQT。
- best pose PDB。
- `figure_a_overall_pose.png`。
- `figure_b_binding_site.png`。
- `pymol_docking_session.pse`。
- `manifest.json`。

批量模式额外输出：

- `batch_manifest.json`：记录输入 CSV、总行数、成功/失败数、受体缓存文件、每行 manifest 和错误信息。
- `batch_results_summary.csv`：记录每一行的 success 或 failed 状态。
- `docking_results_summary.xlsx`：只追加成功对接的结果行。

## 作为 Codex 或 Claude Skill 使用

本项目按 agent skill 的结构组织。若要在 Codex 中使用，可以将 `autodock-vina-docking/` 目录放到 Codex skills 目录，例如：

```text
~/.codex/skills/autodock-vina-docking/
```

然后开启新的 Codex 会话，用自然语言调用，例如：

```text
使用 autodock-vina-docking skill，帮我对接 S100A8 / 5HLO 和 DL-Tryptophan / 1148。
```

批量对接时，提供 CSV 文件并让 agent 运行其中所有行即可。

## 结果解释注意事项

Docking score 和 docking pose 是计算预测结果。它们适合用于提出假设、学习流程和初步排序，但不能作为小分子真实结合或生物活性的实验证据。

如果没有检测到可信的共晶配体，本流程会使用 blind docking，并在 config、日志、summary 和 manifest 中写入 warning。Blind docking 更容易产生假阳性 pose，需要谨慎解释。

## 致谢

本流程基于多个开源科学计算工具和公共数据库，包括 AutoDock Vina、RDKit、Meeko、Open Babel、ProDy、PubChem PUG REST 和 PyMOL open source。如果在科研产出中使用本流程，请引用相关工具和数据库。

## License

本项目采用 [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) 协议开源。欢迎非商业用途的分享与改编，但请注明出处。
