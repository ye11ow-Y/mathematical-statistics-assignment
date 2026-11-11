"""Generate English figures using the established plots and saved R outputs."""

import json
import re
import warnings

import matplotlib.pyplot as plt
from matplotlib.text import Text

import plot_figures as original
from plot_figures_final import drawing_signature


DESTINATION = original.ROOT / 'output/figures_en'
LABELS = {
    r'样本量 $n$': r'Sample size $n$',
    r'均值估计 $\bar X_n$': r'Mean estimate $\bar X_n$',
    r'方差估计 $S_n^2$': r'Variance estimate $S_n^2$',
    '概率密度': 'Probability density',
    r'参数 $\theta$': r'Parameter $\theta$',
    r'后验均值与可信区间': 'Posterior mean',
    r'功效 $1-\beta_n$': r'Power $1-\beta_n$',
    '第一类错误率': 'Type I error probability',
    '提交排序位置（由早到晚）': 'Submission rank',
    '相对提交时间（小时）': 'Relative submission time (hours)',
}
FIGURES = {
    'fig01_consistency': {
        'title': 'Consistency of the mean and variance estimators',
        'legend': ['Estimate', 'True value', 'Pointwise 95% band'],
    },
    'fig02_standardized_distribution': {
        'title': 'Sampling distribution of the standardized mean',
        'legend': ['Empirical density', r'$N(0,1)$'],
    },
    'fig03_bayesian_update': {
        'title': 'Bayesian updating under three priors',
        'legend': [r'Beta$(3,7)$', 'Jeffreys', 'Uniform', r'True value $0.3$'],
    },
    'fig04_test_power': {
        'title': 'Power and Type I error probability',
        'legend': ['Theoretical probability', 'Simulation and 95% interval'],
    },
    'fig05_submission_prediction': {
        'title': 'Ordered predictions of relative submission times',
        'legend': ['Predictive median', 'Earlier scenario', 'Later scenario'],
    },
}


def finish_english_figure(fig, stem, title, caption, sources, note):
    signature = drawing_signature(fig)
    for ax in fig.axes:
        for location in ['left', 'center', 'right']:
            ax.set_title('', loc=location)
        for annotation in list(ax.texts):
            annotation.remove()
        ax.set_xlabel(LABELS.get(ax.get_xlabel(), ax.get_xlabel()))
        ax.set_ylabel(LABELS.get(ax.get_ylabel(), ax.get_ylabel()))
    for legend in fig.legends:
        for label, english in zip(legend.get_texts(), FIGURES[stem]['legend'], strict=True):
            label.set_text(english)
    if stem == 'fig02_standardized_distribution':
        for ax in fig.axes[2:]:
            ax.set_xlabel(r'Standardized mean $Z_n$')
    assert drawing_signature(fig) == signature
    assert all(not re.search(r'[\u4e00-\u9fff]', item.get_text()) for item in fig.findobj(Text))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for label in fig.findobj(Text):
        if not label.get_visible() or not label.get_text().strip():
            continue
        box = label.get_window_extent(renderer)
        assert box.x0 >= -1 and box.y0 >= -1, label.get_text()
        assert box.x1 <= fig.bbox.width + 1 and box.y1 <= fig.bbox.height + 1, label.get_text()
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(DESTINATION / f'{stem}.{extension}')
    FIGURES[stem]['drawing_signature'] = signature
    FIGURES[stem]['sources'] = sources
    original.CHECKS.append(f'{stem}: English text verified, plotted data and geometry preserved, all labels within canvas.')
    plt.close(fig)


def main():
    DESTINATION.mkdir(parents=True, exist_ok=True)
    original.configure_style()
    original.finish_figure = finish_english_figure
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        original.plot_consistency()
        original.plot_standardized_distribution()
        original.plot_bayesian_update()
        original.plot_test_power()
        original.plot_submission_prediction()
    assert not caught, '\n'.join(str(item.message) for item in caught)
    reference = json.loads((original.ROOT / 'output/figures_final/figure_manifest.json').read_text(encoding='utf-8'))
    expected = {item['stem']: item['drawing_signature'] for item in reference['figures']}
    assert {stem: item['drawing_signature'] for stem, item in FIGURES.items()} == expected
    manifest = {'language': 'English', 'sources_sha256': original.SOURCE_HASHES, 'figures': FIGURES}
    (DESTINATION / 'figure_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    (DESTINATION / 'verification.txt').write_text('\n'.join(original.CHECKS) + '\n', encoding='utf-8')
    print('\n'.join(original.CHECKS))


if __name__ == '__main__':
    main()
