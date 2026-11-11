"""Build the five assignment figures from the saved R calculation tables.

Run from any directory with ``python analysis/plot_figures.py``.
"""

from hashlib import sha256
import json
from pathlib import Path
import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
DESTINATION = ROOT / "output" / "figures"
WIDTH = 17 / 2.54
BLUE = "#0F4D92"
GREEN = "#42949E"
RED = "#B64342"
INK = "#272727"
GRAY = "#767676"
GRID = "#E8E8E8"
PRIOR_STYLES = [
    ("Beta(3,7)", r"Beta$(3,7)$", BLUE, "-"),
    ("Jeffreys", "Jeffreys", GREEN, "--"),
    ("Flat", "平坦先验", RED, "-."),
]
SOURCE_HASHES = {}
FIGURES = []
CHECKS = []


def read_table(relative_path):
    path = ROOT / "output" / relative_path
    SOURCE_HASHES[path.relative_to(ROOT).as_posix()] = sha256(path.read_bytes()).hexdigest()
    table = pd.read_csv(path)
    assert not table.isna().any().any(), f"Missing values in {relative_path}"
    return table


def configure_style():
    for family in ["Cambria", "DengXian"]:
        font_manager.findfont(font_manager.FontProperties(family=family), fallback_to_default=False)
    plt.rcParams.update({
        "font.family": ["Cambria", "DengXian"],
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.titlepad": 9,
        "axes.labelpad": 5,
        "axes.edgecolor": GRAY,
        "axes.linewidth": 0.65,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.axisbelow": True,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 1.45,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "legend.handlelength": 2.4,
        "legend.columnspacing": 1.6,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "statistics-assignment-20261110",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": True,
    })


def sample_axis(ax, ticks, logarithmic=True):
    if logarithmic:
        ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel(r"样本量 $n$")


def finish_figure(fig, stem, title, caption, sources, note):
    fig.text(0.02, 0.012, note, fontsize=7.3, color=GRAY, va="bottom")
    # MathText uses a single default face for the text surrounding formulas.
    # Give mixed Chinese/formula labels a face containing Chinese glyphs.
    for label in fig.findobj(Text):
        text = label.get_text()
        if "$" in text and any("\u4e00" <= character <= "\u9fff" for character in text):
            label.set_fontfamily("DengXian")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for label in fig.findobj(Text):
        if not label.get_visible() or not label.get_text().strip():
            continue
        box = label.get_window_extent(renderer)
        assert box.x0 >= -1 and box.y0 >= -1, f"Text outside canvas in {stem}: {label.get_text()}"
        assert box.x1 <= fig.bbox.width + 1 and box.y1 <= fig.bbox.height + 1, f"Text outside canvas in {stem}: {label.get_text()}"
    for extension in ["png", "pdf", "svg"]:
        fig.savefig(DESTINATION / f"{stem}.{extension}")
    FIGURES.append({
        "stem": stem,
        "title": title,
        "caption": caption,
        "sources": sources,
        "width_cm": round(fig.get_figwidth() * 2.54, 2),
        "height_cm": round(fig.get_figheight() * 2.54, 2),
    })
    plt.close(fig)


def plot_consistency():
    source = "q1-q3/q1_consistency_path.csv"
    data = read_table(source)
    n = data.n.to_numpy()
    assert len(data) == 501 and np.all(np.diff(n) > 0)
    assert n[0] == 10 and n[-1] == 10000
    np.testing.assert_allclose(data.mean_reference_lower, 1 + stats.norm.ppf(0.025) * 5 / np.sqrt(n))
    np.testing.assert_allclose(data.variance_reference_upper, 25 * stats.chi2.ppf(0.975, n - 1) / (n - 1))
    CHECKS.append("Figure 1: 501 nested-sample points; normal and chi-square reference limits verified.")

    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.15))
    fig.subplots_adjust(left=0.09, right=0.96, bottom=0.23, top=0.79, wspace=0.31)
    settings = [
        ("mean", "mean_reference_lower", "mean_reference_upper", 1, r"(a) 样本均值 $\bar X_n$", r"均值估计 $\bar X_n$"),
        ("variance", "variance_reference_lower", "variance_reference_upper", 25, r"(b) 样本方差 $S_n^2$", r"方差估计 $S_n^2$"),
    ]
    for ax, (estimate, lower, upper, truth, title, ylabel) in zip(axes, settings):
        ax.fill_between(n, data[lower], data[upper], color=BLUE, alpha=0.12, linewidth=0)
        ax.axhline(truth, color=INK, ls="--", lw=1)
        ax.plot(n, data[estimate], color=BLUE, lw=1.15)
        sample_axis(ax, [10, 100, 1000, 10000])
        ax.set_xlim(10, 10000)
        ax.set_title(title, loc="left")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color=GRID, lw=0.6)
    axes[0].set_ylim(-2.6, 4.6)
    axes[0].set_yticks([-2, 0, 1, 2, 4])
    axes[1].set_ylim(5, 55)
    axes[1].set_yticks([10, 25, 40, 55])
    handles = [Line2D([], [], color=BLUE, label="同一序列的估计轨迹"),
               Line2D([], [], color=INK, ls="--", lw=1, label="参数真值"),
               Patch(facecolor=BLUE, alpha=0.12, label="逐点 95% 抽样参考带")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.51, 0.98), ncol=3)
    finish_figure(fig, "fig01_consistency", "图 1　均值与方差估计的一致性",
                  "同一条正态样本序列的均值与方差估计随样本量变化。虚线标出真值 1 和 25，阴影为按已知模拟真值计算的逐点 95% 抽样参考带。横轴为对数刻度。",
                  [source], "参考带表示各样本量下估计量的抽样范围，使用逐点覆盖率；轨迹来自同一条嵌套样本序列。")


def plot_standardized_distribution():
    source = "q1-q3/q1_normality_samples.csv"
    data = read_table(source)
    assert data.groupby("n").size().to_dict() == {10: 1000, 100: 1000, 1000: 1000, 10000: 1000}
    np.testing.assert_allclose(data.z_true_sigma, np.sqrt(data.n) * (data["mean"] - 1) / 5, atol=1e-12)
    bins = np.linspace(-4.2, 4.2, 29)
    assert data.z_true_sigma.between(bins[0], bins[-1]).all()
    CHECKS.append("Figure 2: true-sigma Z statistic verified for all 4,000 records; shared bins contain all samples.")

    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 4.75), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.085, right=0.98, bottom=0.16, top=0.86, hspace=0.65, wspace=0.20)
    x = np.linspace(-4.2, 4.2, 600)
    for ax, n, letter in zip(axes.flat, [10, 100, 1000, 10000], "abcd"):
        sample = data.loc[data.n.eq(n), "z_true_sigma"].to_numpy()
        ax.hist(sample, bins=bins, density=True, color=BLUE, alpha=0.24,
                edgecolor="white", linewidth=0.5)
        ax.plot(x, stats.norm.pdf(x), color=BLUE)
        ax.set_title(rf"({letter}) $n={n}$", loc="left")
        ax.text(0.97, 0.94, f"均值 {sample.mean():.3f}\n方差 {sample.var():.3f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.6, color=GRAY)
        ax.set_xlim(-4.2, 4.2)
        ax.set_ylim(0, 0.50)
        ax.set_xticks([-4, -2, 0, 2, 4])
        ax.set_yticks([0, 0.2, 0.4])
        ax.grid(axis="y", color=GRID, lw=0.6)
    for ax in axes[:, 0]:
        ax.set_ylabel("概率密度")
    for ax in axes[1, :]:
        ax.set_xlabel(r"$Z_n=\sqrt{n}(\bar X_n-1)/5$")
    handles = [Patch(facecolor=BLUE, alpha=0.24, label="经验密度（每组 1000 次）"),
               Line2D([], [], color=BLUE, label=r"标准正态密度 $N(0,1)$")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.51, 0.98), ncol=2)
    finish_figure(fig, "fig02_standardized_distribution", "图 2　标准化均值的抽样分布",
                  "四种样本量下标准化均值的经验密度与标准正态密度。各组重复 1000 次，使用相同分箱及坐标范围，面板内列出经验均值和分母为 B 的经验方差。",
                  [source], "正态总体下，使用真标准差构造的标准化均值在每个样本量下都精确服从标准正态分布。")


def plot_bayesian_update():
    density_source = "q4-q6/q5_density_n20.csv"
    path_source = "q4-q6/q5_posterior_path.csv"
    priors_source = "q4-q6/q5_priors.csv"
    density = read_table(density_source)
    path = read_table(path_source)
    priors = read_table(priors_source)
    assert len(path) == 150 and len(density) == 3003
    np.testing.assert_allclose(path.posterior_mean, path.a_post / (path.a_post + path.b_post))
    np.testing.assert_allclose(path.lower90, stats.beta.ppf(0.05, path.a_post, path.b_post), atol=1e-12)
    np.testing.assert_allclose(path.upper90, stats.beta.ppf(0.95, path.a_post, path.b_post), atol=1e-12)
    for prior, _, _, _ in PRIOR_STYLES:
        initial = path.loc[path.prior.eq(prior)].iloc[0]
        curve = density.loc[density.prior.eq(prior)]
        np.testing.assert_allclose(curve.posterior_density,
                                   stats.beta.pdf(curve.theta, initial.a_post, initial.b_post), atol=1e-10)
    assert np.all(priors.prior_strength == priors.a + priors.b)
    CHECKS.append("Figure 3: all posterior means, 90% beta quantiles and n=20 posterior densities independently verified.")

    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 5.0))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.16, top=0.85, hspace=0.66, wspace=0.28)
    ax = axes[0, 0]
    for prior, label, color, linestyle in PRIOR_STYLES:
        curve = density.loc[density.prior.eq(prior)]
        ax.plot(curve.theta, curve.posterior_density, color=color, ls=linestyle, label=label)
    ax.axvline(0.3, color=INK, lw=0.9, ls=":")
    ax.set_title(r"(a) 初始后验密度，$n=20$", loc="left")
    ax.set_xlabel(r"参数 $\theta$")
    ax.set_ylabel("概率密度")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 8)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.grid(axis="y", color=GRID, lw=0.6)
    for ax, (prior, label, color, linestyle), letter in zip(
            [axes[0, 1], axes[1, 0], axes[1, 1]], PRIOR_STYLES, "bcd"):
        series = path.loc[path.prior.eq(prior)].sort_values("n")
        ax.fill_between(series.n, series.lower90, series.upper90, color=color, alpha=0.14, lw=0)
        ax.plot(series.n, series.posterior_mean, color=color, ls=linestyle)
        ax.axhline(0.3, color=INK, lw=0.9, ls=":")
        sample_axis(ax, [20, 100, 500, 1000])
        ax.set_xlim(20, 1000)
        ax.set_ylim(0, 0.46)
        ax.set_yticks([0, 0.15, 0.30, 0.45])
        ax.set_title(f"({letter}) {label}", loc="left")
        ax.set_ylabel(r"后验均值与可信区间")
        ax.grid(axis="y", color=GRID, lw=0.6)
    handles = [Line2D([], [], color=color, ls=linestyle, label=label)
               for _, label, color, linestyle in PRIOR_STYLES]
    handles.append(Line2D([], [], color=INK, ls=":", lw=0.9, label=r"真值 $\theta=0.3$"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.51, 0.98), ncol=4)
    finish_figure(fig, "fig03_bayesian_update", "图 3　三种先验下的贝叶斯更新",
                  "初始样本下的三种后验密度，以及同一条 Bernoulli 样本序列下的后验更新。后三个面板采用相同纵轴范围，曲线为后验均值，阴影为各样本量的 90% 等尾可信区间，点线为真值 0.3。样本量横轴为对数刻度。",
                  [density_source, path_source, priors_source], "(b)–(d) 的阴影为 90% 等尾可信区间，曲线为后验均值；三种先验共享同一条样本序列。")


def plot_test_power():
    summary_source = "q4-q6/q6b_summary.csv"
    curve_source = "q4-q6/q6b_theoretical_curve.csv"
    summary = read_table(summary_source)
    curve = read_table(curve_source)
    expected_size = stats.norm.sf(1.65)
    np.testing.assert_allclose(curve.power, stats.norm.sf(1.65 - 0.5 * np.sqrt(curve.n)), atol=1e-14)
    np.testing.assert_allclose(curve.type2, stats.norm.cdf(1.65 - 0.5 * np.sqrt(curve.n)), rtol=1e-11)
    np.testing.assert_allclose(summary.empirical_probability, summary["count"] / summary.total)
    for row in summary.itertuples():
        lower = 0.0 if row.count == 0 else stats.beta.ppf(0.025, row.count, row.total - row.count + 1)
        upper = 1.0 if row.count == row.total else stats.beta.ppf(0.975, row.count + 1, row.total - row.count)
        np.testing.assert_allclose([row.mc_lower95, row.mc_upper95], [lower, upper], atol=1e-12)
    # Evaluate both tails directly; power can round to one for large n.
    stable_mcse = np.sqrt(stats.norm.sf(1.65 - 0.5 * np.sqrt(curve.n)) * curve.type2 / 10000)
    assert stable_mcse.iloc[-1] > 0
    CHECKS.append("Figure 4: direct normal tails and all 10 Clopper-Pearson intervals verified; beta_1024 remains positive.")
    missed = summary.loc[(summary.theoretical_probability < summary.mc_lower95)
                         | (summary.theoretical_probability > summary.mc_upper95)]
    CHECKS.append("Figure 4: intervals missing their theoretical probability: "
                  + ", ".join(f"n={row.n}/{row.experiment}" for row in missed.itertuples()) + ".")

    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.2))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.23, top=0.79, wspace=0.33)
    for ax, experiment, title, ylabel in zip(axes, ["power", "size"],
            ["(a) 备择假设下的功效", "(b) 原假设下的第一类错误率"],
            [r"功效 $1-\beta_n$", "第一类错误率"]):
        points = summary.loc[summary.experiment.eq(experiment)].sort_values("n")
        reference = curve.power if experiment == "power" else np.repeat(expected_size, len(curve))
        ax.plot(curve.n, reference, color=BLUE, lw=1.5)
        ax.errorbar(points.n, points.empirical_probability,
                    yerr=np.vstack([points.empirical_probability - points.mc_lower95,
                                    points.mc_upper95 - points.empirical_probability]),
                    fmt="o", ms=4, color=RED, markerfacecolor="white", markeredgewidth=1,
                    elinewidth=1, capsize=3, capthick=1, clip_on=False)
        sample_axis(ax, [4, 16, 64, 256, 1024])
        ax.set_xlim(3, 1400)
        ax.set_title(title, loc="left")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color=GRID, lw=0.6)
    axes[0].set_ylim(0, 1.04)
    axes[0].set_yticks([0, 0.25, 0.5, 0.75, 1])
    axes[1].set_ylim(0.035, 0.065)
    axes[1].set_yticks([0.035, 0.045, 0.055, 0.065])
    handles = [Line2D([], [], color=BLUE, label="理论概率（临界值 1.65）"),
               Line2D([], [], color=RED, marker="o", markerfacecolor="white", ls="none",
                      label="模拟比例与 95% 二项精确区间")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.51, 0.98), ncol=2)
    finish_figure(fig, "fig04_test_power", "图 4　检验功效与第一类错误率",
                  "在 θ₀=0、θ₁=0.5、临界值为 1.65 的设置下，展示理论功效和第一类错误率，以及五种样本量各 10000 次模拟的比例与 95% Clopper–Pearson 二项精确区间。横轴为对数刻度。n=1024 时第一类错误率的该次模拟区间略低于理论水平。",
                  [summary_source, curve_source], r"$\theta_0=0,\ \theta_1=0.5$；每种样本量、每个假设各模拟 10000 次；理论检验水平约为 0.04947。")


def plot_submission_prediction():
    source = "q7/predictions.csv"
    scenario_source = "q7/scenarios.csv"
    data = read_table(source)
    scenarios = read_table(scenario_source)
    population = len(data)
    assert data.N_students.eq(population).all()
    np.testing.assert_array_equal(data["rank"], np.arange(1, population + 1))
    np.testing.assert_allclose(scenarios.weight.sum(), 1)
    assert np.all(np.diff(data.hours) <= 0)
    np.testing.assert_allclose(data.predictive_cdf_at_prediction, 0.5, atol=1e-7)
    assert np.all(data.later_order_median <= data.hours)
    assert np.all(data.hours <= data.earlier_order_median)
    CHECKS.append(f"Figure 5: {population} descending, unrounded predictions; all median CDF residuals within tolerance.")

    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.35), gridspec_kw={"width_ratios": [1.12, 1]})
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.24, top=0.79, wspace=0.32)
    for ax in axes:
        ax.plot(data["rank"], data.earlier_order_median, color=GREEN, ls="--", lw=1.05)
        ax.plot(data["rank"], data.later_order_median, color=RED, ls="-.", lw=1.05)
        ax.plot(data["rank"], data.hours, color=BLUE, lw=1.7, zorder=4)
        ax.axhline(0, color=INK, lw=0.9, ls=":")
        ax.set_xlabel("提交排序位置（由早到晚）")
        ax.set_ylabel("相对提交时间（小时）")
        ax.grid(axis="y", color=GRID, lw=0.6)
    axes[0].set_title(f"(a) 全体排序预测，{population} 人", loc="left")
    axes[0].set_xlim(1, population)
    axes[0].set_xticks([1, 30, 60, 90, 120, population])
    axes[0].set_ylim(-6, 132)
    axes[0].set_yticks([0, 24, 60, 90, 120])
    axes[0].text(0.97, 0.94, "正值为提前\n负值为迟交", transform=axes[0].transAxes,
                 ha="right", va="top", fontsize=8, color=GRAY)
    axes[1].set_title("(b) 截止附近放大", loc="left")
    axes[1].set_xlim(population - 26, population)
    axes[1].set_xticks([120, 130, 140, population])
    axes[1].set_ylim(-3.7, 5.2)
    axes[1].set_yticks([-3, -1.5, 0, 1.5, 3, 4.5])
    late = data.loc[data.hours.lt(0)]
    axes[1].scatter(late["rank"], late.hours, color=BLUE, s=15, zorder=5)
    axes[1].annotate(f"最晚预测 {data.hours.iloc[-1]:.2f} 小时", xy=(population, data.hours.iloc[-1]),
                     xytext=(population - 18, -2.7), fontsize=7.8,
                     arrowprops={"arrowstyle": "-", "color": GRAY, "lw": 0.7})
    handles = [Line2D([], [], color=BLUE, lw=1.7, label="最终情景预测中位数"),
               Line2D([], [], color=GREEN, ls="--", label="整体偏早情景"),
               Line2D([], [], color=RED, ls="-.", label="整体偏晚情景")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.51, 0.98), ncol=3)
    finish_figure(fig, "fig05_submission_prediction", "图 5　相对提交时间的排序预测",
                  f"{population} 人的相对提交时间排序预测与截止附近放大图。蓝线为最终情景预测中位数，虚线与点划线分别为整体偏早、偏晚情景下的次序中位数，水平点线为截止时刻。正值表示提前提交，负值表示迟交。情景曲线用于展示既定假设下的敏感性。",
                  [source, scenario_source], "截止时刻为北京时间 2026-11-12 17:25；偏早、偏晚曲线表示情景中位数；计算与作图使用未舍入数值。")


def build_overview():
    fig = plt.figure(figsize=(12, 12.5), facecolor="white")
    slots = [(0.03, 0.69, 0.46, 0.26), (0.51, 0.63, 0.46, 0.35),
             (0.03, 0.285, 0.46, 0.36), (0.51, 0.31, 0.46, 0.29),
             (0.25, 0.015, 0.50, 0.28)]
    for item, bounds in zip(FIGURES, slots):
        ax = fig.add_axes(bounds)
        ax.imshow(plt.imread(DESTINATION / f"{item['stem']}.png"))
        ax.set_axis_off()
        ax.set_title(item["title"], fontsize=12, loc="left", pad=7)
    fig.savefig(DESTINATION / "figure_overview.png", dpi=150)
    plt.close(fig)


def write_delivery_files():
    notes = [
        "# 作业配图", "",
        "五幅图对应最新版正文预留的图 1–5。每幅提供 300 dpi PNG、嵌入字体的 PDF 和保留文本的 SVG，宽度统一为 17 cm。中文使用等线，西文使用 Cambria，数学符号使用 STIX。",
        "", "## 图形选择", "",
        "轨迹图呈现估计随样本量的变化，分面直方图比较标准化分布，后验密度与区间图呈现先验影响及不确定性，概率曲线与二项区间图核对检验表现，排序曲线与局部放大图呈现提交预测。第 2–4 题沿用正文中的推导和数值表。",
        "", "## 文件与图注", "",
    ]
    for item in FIGURES:
        stem = item["stem"]
        notes.extend([f"### {item['title']}", "",
                      f"[PNG]({stem}.png) · [PDF]({stem}.pdf) · [SVG]({stem}.svg)", "",
                      item["caption"], "",
                      "数据来自 " + "、".join(f"`output/{source}`" for source in item["sources"]) + "。", ""])
    notes.extend([
        "## 统计口径", "",
        "图 1 使用逐点抽样参考带，按模拟真值计算。图 3 使用各样本量的后验可信区间。图 4 的误差条表示 Monte Carlo 模拟比例的二项精确区间。图 5 的三条曲线表示不同情景下的次序中位数。",
        "",
        "第 6B 的 n=1024、H0 模拟中，第一类错误率为 0.0449，95% 二项精确区间约为 [0.040924, 0.049143]，理论水平约为 0.049471。图中保留该次模拟的结果。n=256、1024 的经验功效为 1，相应 95% 区间仍有宽度。正态上下尾分别计算，避免第二类错误率在计算中消失。",
        "",
        "第七题继续使用人数 146 和正文已声明的情景权重及参数，预测误差需观测真实提交数据后计算。图中使用未舍入值，排名 143 的 0.00 小时在原始结果中为正值。",
        "", "## 复现", "",
        "在项目目录执行 `python analysis/plot_figures.py`。依赖 Python、NumPy、pandas、SciPy、Matplotlib，以及系统等线和 Cambria 字体。该脚本读取已有 R 结果表，核对理论公式与区间，并生成图片及配图版 Markdown。",
        "",
        "`figure_manifest.json` 记录源数据 SHA-256、输出尺寸、风格和图注。`verification.txt` 记录本次数据与导出检查。", "",
    ])
    (DESTINATION / "README.md").write_text("\n".join(notes), encoding="utf-8")
    original = ROOT / "output" / "数理统计作业正文_排版修订版.md"
    manuscript = original.read_text(encoding="utf-8")
    for number, item in enumerate(FIGURES, start=1):
        placeholders = [line for line in manuscript.splitlines() if line.startswith(f"【图 {number} 待插入")]
        assert len(placeholders) == 1, f"Expected one placeholder for figure {number}"
        replacement = f"![{item['title']}](figures/{item['stem']}.png)\n\n{item['title']}。{item['caption']}"
        manuscript = manuscript.replace(placeholders[0], replacement)
    (ROOT / "output" / "数理统计作业正文_配图版.md").write_text(manuscript, encoding="utf-8")
    manifest = {
        "source_calculations": "Existing R simulation outputs; no new random draws",
        "style": {"latin": "Cambria", "chinese": "DengXian", "math": "STIX",
                  "blue": BLUE, "teal": GREEN, "red": RED, "dpi": 300},
        "sources_sha256": SOURCE_HASHES,
        "figures": FIGURES,
    }
    (DESTINATION / "figure_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def main():
    DESTINATION.mkdir(parents=True, exist_ok=True)
    configure_style()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        plot_consistency()
        plot_standardized_distribution()
        plot_bayesian_update()
        plot_test_power()
        plot_submission_prediction()
        build_overview()
    assert not caught, "Plot warnings: " + "\n".join(str(item.message) for item in caught)
    write_delivery_files()
    for item in FIGURES:
        for extension in ["png", "pdf", "svg"]:
            path = DESTINATION / f"{item['stem']}.{extension}"
            assert path.is_file() and path.stat().st_size > 1000, f"Invalid export {path}"
    CHECKS.append("Export: 5 figures x 3 formats, overview and illustrated manuscript generated; no plotting warnings.")
    (DESTINATION / "verification.txt").write_text("\n".join(CHECKS) + "\n", encoding="utf-8")
    for line in CHECKS:
        print(line)


if __name__ == "__main__":
    main()
