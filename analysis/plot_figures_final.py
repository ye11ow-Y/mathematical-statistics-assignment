"""Prepare the selected first-version figures for final user confirmation."""

from hashlib import sha256
import json
import warnings

import plot_figures as original

import matplotlib.pyplot as plt
import numpy as np


ROOT = original.ROOT
DESTINATION = ROOT / "output" / "figures_final"
MANUSCRIPT_NAME = "数理统计作业正文_配图确认版.md"
WORD_LABELS = {
    "fig01_consistency": {
        "panels": ["（a）左图　样本均值的估计轨迹", "（b）右图　样本方差的估计轨迹"],
        "caption": "虚线为真值 1 和 25，阴影为逐点 95% 抽样参考带。两条轨迹来自同一条嵌套样本序列，横轴为对数刻度。",
        "legend": ["估计轨迹", "真值", "逐点 95% 参考带"],
    },
    "fig02_standardized_distribution": {
        "panels": ["（a）左上　n=10", "（b）右上　n=100", "（c）左下　n=1000", "（d）右下　n=10000"],
        "caption": "直方图为每组 1000 次模拟的经验密度，曲线为标准正态密度。四组使用相同分箱和坐标范围。",
        "legend": ["经验密度", r"$N(0,1)$"],
    },
    "fig03_bayesian_update": {
        "panels": ["（a）左上　n=20 时的后验密度", "（b）右上　Beta(3,7) 先验下的后验更新",
                   "（c）左下　Jeffreys 先验下的后验更新", "（d）右下　平坦先验下的后验更新"],
        "caption": "后三图曲线为后验均值，阴影为 90% 等尾可信区间，点线为真值 0.3。三种先验使用同一条样本序列，样本量横轴为对数刻度。",
        "legend": [r"Beta$(3,7)$", "Jeffreys", "平坦先验", r"真值 $0.3$"],
    },
    "fig04_test_power": {
        "panels": ["（a）左图　备择假设下的功效", "（b）右图　原假设下的第一类错误率"],
        "caption": "蓝线为理论概率，红点及误差条为模拟比例与 95% 二项精确区间。取 θ₀=0、θ₁=0.5、临界值 1.65，各组模拟 10000 次，横轴为对数刻度。",
        "legend": ["理论概率", "模拟与 95% 区间"],
    },
    "fig05_submission_prediction": {
        "panels": ["（a）左图　全体学生的排序预测", "（b）右图　截止附近的局部放大"],
        "caption": "蓝线为最终预测中位数，青色虚线与红色点划线为偏早、偏晚情景的次序中位数。零线为截止时刻，正值为提前，负值为迟交，人数为 146。",
        "legend": ["预测中位数", "偏早情景", "偏晚情景"],
    },
}


def drawing_signature(fig):
    """Hash plotted data, styles and geometry independently of text annotations."""
    digest = sha256()
    for ax in fig.axes:
        digest.update(repr((ax.get_position().bounds, ax.get_xlim(), ax.get_ylim(),
                            ax.get_xscale(), ax.get_yscale())).encode())
        for line in ax.lines:
            digest.update(np.asarray(line.get_xydata(), dtype=float).tobytes())
            digest.update(repr((line.get_color(), line.get_linewidth(), line.get_linestyle(),
                                line.get_marker(), line.get_markersize())).encode())
        for patch in ax.patches:
            digest.update(patch.get_path().vertices.tobytes())
            digest.update(patch.get_transform().get_matrix().tobytes())
            digest.update(repr((patch.get_facecolor(), patch.get_edgecolor(), patch.get_linewidth())).encode())
        for collection in ax.collections:
            digest.update(np.asarray(collection.get_offsets()).tobytes())
            digest.update(collection.get_facecolors().tobytes())
            digest.update(collection.get_edgecolors().tobytes())
            digest.update(collection.get_linewidths().tobytes())
            for path in collection.get_paths():
                digest.update(path.vertices.tobytes())
    return digest.hexdigest()


def snapshot_previous_versions():
    paths = []
    for name in ["figures", "figures_v2", "figures_v3"]:
        paths.extend((ROOT / "output" / name).glob("*"))
    paths.extend(path for path in (ROOT / "output").glob("*.md") if path.name != MANUSCRIPT_NAME)
    paths.extend([ROOT / "analysis" / "plot_figures.py", ROOT / "analysis" / "plot_figures_v2.py",
                  ROOT / "analysis" / "plot_figures_v3.py"])
    return {path.relative_to(ROOT).as_posix(): sha256(path.read_bytes()).hexdigest()
            for path in paths if path.is_file()}


def make_overview():
    fig = plt.figure(figsize=(12, 12.5), facecolor="white")
    slots = [(0.025, 0.69, 0.46, 0.26), (0.515, 0.64, 0.46, 0.34),
             (0.025, 0.285, 0.46, 0.36), (0.515, 0.315, 0.46, 0.28),
             (0.255, 0.015, 0.49, 0.28)]
    for item, bounds in zip(original.FIGURES, slots):
        ax = fig.add_axes(bounds)
        ax.imshow(plt.imread(DESTINATION / f"{item['stem']}.png"))
        ax.set_axis_off()
    fig.savefig(DESTINATION / "confirmation_overview.png", dpi=150)
    plt.close(fig)


def write_documents():
    manuscript = (ROOT / "output" / "数理统计作业正文_排版修订版.md").read_text(encoding="utf-8")
    notes = ["# 配图最终确认版", "",
             "本套采用用户选定的第一版五幅图。图内保留坐标和简短图例，图名、子图名与说明集中在下文及配图确认版正文中，供录入 Word。", "",
             "[五幅图总览](confirmation_overview.png)", "",
             "各图均提供 17 cm 宽的 300 dpi PNG、嵌入字体的 PDF，以及保留文本的 SVG。", ""]
    for number, item in enumerate(original.FIGURES, 1):
        stem = item["stem"]
        labels = WORD_LABELS[stem]
        item["panel_titles_for_word"] = labels["panels"]
        placeholder = next(line for line in manuscript.splitlines() if line.startswith(f"【图 {number} 待插入"))
        marker = (f"![{item['title']}](figures_final/{stem}.png)\n\n"
                  f"【Word 图名】{item['title']}\n\n"
                  "【Word 子图名】" + "；".join(labels["panels"]) + "\n\n"
                  f"【图注建议】{item['caption']}")
        manuscript = manuscript.replace(placeholder, marker)
        notes.extend([f"## {item['title']}", "",
                      f"[PNG]({stem}.png) · [PDF]({stem}.pdf) · [SVG]({stem}.svg)", "",
                      "；".join(labels["panels"]) + "。", "", item["caption"], ""])
    notes.extend(["## 核对说明", "",
                  "从第一版绘图函数生成五幅图，并在文字整理前后核对曲线、分箱、区间、配色与坐标范围的签名。图内的总图名、子图名、图底说明、经验均值方差标注和提交预测注释统一移至正文或保留于已有结果表中。", "",
                  "第 6B 在 n=1024 时的模拟第一类错误率为 0.0449，该次 95% 二项精确区间约为 [0.040924, 0.049143]，理论水平约为 0.049471。图中呈现该次原始模拟结果。第七题继续使用人数 146 和既定情景假设。", "",
                  "## 复现", "", "在项目目录执行 `python analysis/plot_figures_final.py`。", ""])
    (ROOT / "output" / MANUSCRIPT_NAME).write_text(manuscript, encoding="utf-8")
    (DESTINATION / "README.md").write_text("\n".join(notes), encoding="utf-8")


def main():
    previous = snapshot_previous_versions()
    DESTINATION.mkdir(parents=True, exist_ok=True)
    original.DESTINATION = DESTINATION
    original.configure_style()
    normal_finish = original.finish_figure

    def finalize(fig, stem, title, caption, sources, note):
        signature = drawing_signature(fig)
        for ax in fig.axes:
            for location in ["left", "center", "right"]:
                ax.set_title("", loc=location)
            for text in list(ax.texts):
                text.remove()
        for legend in fig.legends:
            texts = legend.get_texts()
            labels = WORD_LABELS[stem]["legend"]
            assert len(texts) == len(labels)
            for text, label in zip(texts, labels):
                text.set_text(label)
        if stem == "fig02_standardized_distribution":
            for ax in fig.axes[2:]:
                ax.set_xlabel(r"标准化均值 $Z_n$")
        if stem == "fig03_bayesian_update":
            for ax in fig.axes[1:]:
                ax.set_ylabel("后验均值")
        if stem == "fig05_submission_prediction":
            for ax in fig.axes:
                ax.set_xlabel("提交排序位置")
        assert drawing_signature(fig) == signature, f"Drawing changed in {stem}"
        normal_finish(fig, stem, title, WORD_LABELS[stem]["caption"], sources, "")
        original.FIGURES[-1]["drawing_signature"] = signature
        original.CHECKS.append(f"{stem}: first-version plotted data, styles and axes match after text cleanup.")

    original.finish_figure = finalize
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            original.plot_consistency()
            original.plot_standardized_distribution()
            original.plot_bayesian_update()
            original.plot_test_power()
            original.plot_submission_prediction()
            make_overview()
        assert not caught, "\n".join(str(item.message) for item in caught)
    finally:
        original.finish_figure = normal_finish
    write_documents()
    assert snapshot_previous_versions() == previous, "Previous delivery changed"
    original.CHECKS.append("All previous figure versions, scripts and manuscripts retain their SHA-256 hashes.")
    original.CHECKS.append("Five confirmation figures exported to PNG/PDF/SVG without plotting warnings.")
    manifest = {"selected_version": "first", "sources_sha256": original.SOURCE_HASHES,
                "preserved_sha256": previous, "figures": original.FIGURES}
    (DESTINATION / "figure_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (DESTINATION / "verification.txt").write_text("\n".join(original.CHECKS) + "\n", encoding="utf-8")
    print("\n".join(original.CHECKS))


if __name__ == "__main__":
    main()
