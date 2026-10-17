"""Plot palette and typography extracted from the submitted figure generator."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


CLASS_ORDER = [
    "Class I",
    "Class II",
    "Class III",
    "Class IV",
    "Class V",
    "Class VI",
]

CLASS_COLORS = {
    "Class I": "#0F4D92",
    "Class II": "#E28E2C",
    "Class III": "#4D8A64",
    "Class IV": "#B64342",
    "Class V": "#8069A8",
    "Class VI": "#8F8F8F",
}

AXIS_COLORS = {
    "translational": "#0F4D92",
    "functional": "#E28E2C",
    "baseline": "#767676",
}

NATURE = {
    "ink": "#272727",
    "muted": "#606060",
    "line": "#D8D8D8",
    "panel": "#F7F8FA",
    "blue_soft": "#DCEBFA",
    "orange_soft": "#F0E0D0",
    "green_soft": "#DDF3DE",
    "red_soft": "#F6CFCB",
    "teal": "#42949E",
}

matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 8.2,
        "axes.titlesize": 9.0,
        "axes.labelsize": 8.2,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.4,
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.75,
        "axes.edgecolor": "#303030",
        "xtick.color": "#303030",
        "ytick.color": "#303030",
        "axes.labelcolor": "#303030",
    }
)

def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.075, 1.035, letter, transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom", color=NATURE["ink"])

def add_grid(ax: plt.Axes, axis: str = "both") -> None:
    ax.grid(axis=axis, color=NATURE["line"], linewidth=0.45, alpha=0.78)
    ax.set_axisbelow(True)
