import pandas as pd, numpy as np, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

c = pd.read_csv("候选蛋白_汇总.csv")
os.makedirs("figures", exist_ok=True)

color = {"强（H4≥0.8）": "#1f4e79", "中（0.5≤H4<0.8）": "#c55a11",
         "弱（H4<0.5，支持 LD 混杂）": "#a6a6a6", "未评估（共定位缺失）": "#d9d9d9"}
legend_label = {"强（H4≥0.8）": "Strong (PP.H4 ≥ 0.8)",
                "中（0.5≤H4<0.8）": "Moderate (0.5–0.8)",
                "弱（H4<0.5，支持 LD 混杂）": "Weak (PP.H3 dominant)",
                "未评估（共定位缺失）": "Not evaluated"}

# 图1 森林图（UACR）
d = c.dropna(subset=["UACR_beta"]).sort_values("UACR_beta")
fig, ax = plt.subplots(figsize=(6.5, 4.6))
y = np.arange(len(d))
ax.errorbar(d["UACR_beta"], y, xerr=1.96*d["UACR_se"], fmt="o", ms=5, lw=1.2,
            ecolor="#555555", mfc="white", mec="#222222", capsize=2, zorder=3)
for i, (_, r) in enumerate(d.iterrows()):
    ax.plot(r["UACR_beta"], i, "o", ms=6, color=color.get(r["证据等级"], "#888888"), zorder=4)
    if r["gene"] == "BAG3":
        ax.plot(r["UACR_beta"], i, "o", ms=9, mfc="#c00000", mec="black", mew=0.8, zorder=5)
ax.axvline(0, color="#888888", lw=0.8, ls="--")
ax.set_yticks(y); ax.set_yticklabels(d["gene"], fontsize=9)
ax.set_xlabel("MR effect on log UACR (per 1-SD higher genetically predicted protein)", fontsize=9)
ax.set_title("Cis-pQTL Mendelian randomization: UACR", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=v, markersize=7,
                  label=legend_label[k]) for k, v in color.items()]
handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#c00000",
                      markeredgecolor="black", markersize=9,
                      label="Not replicated (opposite direction)"))
ax.legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
ax.grid(axis="x", color="#eeeeee", lw=0.6)
ax.set_axisbelow(True)
fig.tight_layout(); fig.savefig("figures/Fig1_forest_UACR.png", dpi=300)
fig.savefig("figures/Fig1_forest_UACR.pdf"); plt.close(fig)

# 图2 双终点散点
d2 = c.dropna(subset=["UACR_beta", "eGFR_beta"])
fig, ax = plt.subplots(figsize=(5.2, 5.0))
for _, r in d2.iterrows():
    ax.scatter(r["UACR_beta"], r["eGFR_beta"], s=42, color=color.get(r["证据等级"], "#888888"),
               edgecolor="#333333", linewidth=0.5, zorder=3)
    ax.annotate(r["gene"], (r["UACR_beta"], r["eGFR_beta"]), textcoords="offset points",
                xytext=(5, 4), fontsize=7.5)
ax.axhline(0, color="#888888", lw=0.8, ls="--"); ax.axvline(0, color="#888888", lw=0.8, ls="--")
ax.set_xlabel("MR effect on log UACR", fontsize=9)
ax.set_ylabel("MR effect on eGFR", fontsize=9)
ax.set_title("Concordance between UACR and eGFR effects", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(color="#f0f0f0", lw=0.6); ax.set_axisbelow(True)
fig.tight_layout(); fig.savefig("figures/Fig2_concordance.png", dpi=300)
fig.savefig("figures/Fig2_concordance.pdf"); plt.close(fig)

# 图3 共定位后验
d3 = pd.read_csv("coloc_final.csv").sort_values("H4")
fig, ax = plt.subplots(figsize=(6.0, 4.2))
y = np.arange(len(d3))
ax.barh(y, d3["H4"], color="#1f4e79", height=0.55, label="PP.H4 (shared causal variant)")
ax.barh(y, d3["H3"], left=d3["H4"], color="#bfbfbf", height=0.55, label="PP.H3 (distinct variants)")
ax.axvline(0.8, color="#c00000", lw=1.0, ls="--")
ax.set_yticks(y); ax.set_yticklabels(d3["gene"], fontsize=9)
ax.set_xlim(0, 1.0); ax.set_xlabel("Posterior probability", fontsize=9)
ax.set_title("Colocalization between pQTL and UACR signals", fontsize=10)
ax.legend(fontsize=7, frameon=False, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/Fig3_coloc.png", dpi=300)
fig.savefig("figures/Fig3_coloc.pdf"); plt.close(fig)

for f in sorted(os.listdir("figures")):
    print(f, os.path.getsize(os.path.join("figures", f)))
