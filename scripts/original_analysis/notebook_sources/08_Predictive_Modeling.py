# %% [markdown]
# # 08 · Predictive Modelling — Feasibility and Exploratory Models
# **Research objective:** Objective 4 — whether the data support predictive modelling of AMR-related patterns.
#
# **Feasibility first.** No sample-level clinical, epidemiological, food-safety or business outcome is linked to any metagenome by place and time. Disease prediction, outbreak forecasting, clinical AMR forecasting and food-business forecasting are therefore **not feasible** and are not attempted. Two exploratory, internally valid targets exist:
#
# 1. **Model A — Environmental-matrix signature (ResFinder/KMA).** Can resistome class composition distinguish municipal wastewater from other urban matrices *in a study not seen during training*? Grouped leave-one-study-out cross-validation prevents samples from one study appearing in both training and test folds.
# 2. **Model B — Resistome load from taxonomy (DS02).** Can prokaryotic phylum composition predict total ARG abundance (log10 RPKM) (i) in a later period (forward temporal hold-out) and (ii) in a city not used for training (leave-one-city-out)?
#
# Both are tests of transferable signal, not forecasts of health outcomes. Baselines: majority-class / prior-probability classifier (A) and training-mean predictor (B). Feature transformations (CLR) are computed per sample and contain no information from other samples, so they cannot leak across folds; scaling and model fitting are nested inside each training fold through pipelines. Random seed fixed.

# %%
#BOOTSTRAP
NB = "NB08"
nb_header(NB)
import seaborn as sns
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, HistGradientBoostingClassifier
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.model_selection import LeaveOneGroupOut, StratifiedKFold, cross_val_predict
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, precision_score, recall_score, f1_score, roc_auc_score,
                             average_precision_score, brier_score_loss, confusion_matrix, roc_curve, precision_recall_curve,
                             r2_score, mean_absolute_error, mean_squared_error)
from sklearn.calibration import calibration_curve
from sklearn.inspection import permutation_importance

# %% [markdown]
# ## Model A — municipal wastewater vs other urban matrices from resistome class composition

# %%
us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0)
cfp = pd.read_csv(DER / "rf_class_fpkm.csv", index_col=0)
urb = us[us.Stratum == "Urban environmental"].copy()
Xc = cfp.loc[urb.index].drop(columns=[c for c in ["Unclassified"] if c in cfp.columns])
from sklearn.base import BaseEstimator, TransformerMixin
class FoldwiseZeroReplaceCLR(BaseEstimator, TransformerMixin):
    """AUD-17 fix: multiplicative zero replacement whose constant (frac x smallest non-zero share) is
    estimated on the training rows of each fold only, followed by the centred log-ratio transform."""
    def __init__(self, frac=0.65): self.frac = frac
    def fit(self, X, y=None):
        P = closure(np.asarray(X, float)); self.delta_ = float(self.frac * P[P > 0].min()); return self
    def transform(self, X):
        Pz, _ = mult_replace(np.asarray(X, float), delta=self.delta_); return clr(Pz)
XA = pd.DataFrame(closure(Xc.values), index=Xc.index, columns=Xc.columns)   # closed class shares; CLR is fitted inside each training fold
delta_all = 0.65 * XA.values[XA.values > 0].min()
print("zero-replacement constant on all units (previous, transductive):", delta_all)
yA = (urb.Sample_Type == "Municipal wastewater").astype(int).values; gA = urb.Project_Accession.values
designA = pd.DataFrame([dict(item="Unit of analysis", value="Analytic unit (biological sample)"), dict(item="Units", value=len(yA)),
    dict(item="Positive class (municipal wastewater)", value=int(yA.sum())), dict(item="Negative class (surface water, air, built environment, transit)", value=int((1 - yA).sum())),
    dict(item="Studies (groups)", value=len(set(gA))), dict(item="Features", value=f"CLR of {XA.shape[1]} class-group shares of total FPKM"),
    dict(item="Pre-processing (AUD-17 fix)", value="Multiplicative zero-replacement constant (0.65 x smallest non-zero share) and CLR estimated on the training fold only, inside each model pipeline; test folds never inform pre-processing"),
    dict(item="Validation", value="Stratified group 4-fold CV (whole studies held out), 10 repeats; per-fold metrics averaged (leave-one-study-out used for per-study accuracy, Table M03)"), dict(item="Comparison validation (leakage illustration)", value="Stratified 4-fold ignoring study, 10 repeats")])
save_table(designA, "M01_modelA_design", "Model A design: predicting municipal wastewater from resistome class composition", NB, sheet="Predictive_Models")
models = {"Baseline (prior)": make_pipeline(FoldwiseZeroReplaceCLR(), DummyClassifier(strategy="prior")),
          "Logistic regression (L2)": make_pipeline(FoldwiseZeroReplaceCLR(), StandardScaler(), LogisticRegression(C=1.0, max_iter=5000)),
          "Random forest": make_pipeline(FoldwiseZeroReplaceCLR(), RandomForestClassifier(n_estimators=500, min_samples_leaf=2, class_weight="balanced", random_state=SEED, n_jobs=-1)),
          "Gradient boosting (HistGB)": make_pipeline(FoldwiseZeroReplaceCLR(), HistGradientBoostingClassifier(max_depth=3, learning_rate=.05, max_iter=300, random_state=SEED))}
logo = LeaveOneGroupOut()
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.base import clone
def clf_metrics(y, p, thr=.5):
    yh = (p >= thr).astype(int)
    return dict(accuracy=accuracy_score(y, yh), balanced_accuracy=balanced_accuracy_score(y, yh), precision=precision_score(y, yh, zero_division=0),
                recall=recall_score(y, yh, zero_division=0), f1=f1_score(y, yh, zero_division=0), brier=brier_score_loss(y, p))
def repeated_cv(mdl, X, y, g, grouped=True, seeds=range(10), k=4):
    per, oofs = [], []
    for sd in seeds:
        cv = StratifiedGroupKFold(k, shuffle=True, random_state=sd) if grouped else StratifiedKFold(k, shuffle=True, random_state=sd)
        p = np.zeros(len(y))
        for tr, te in (cv.split(X, y, g) if grouped else cv.split(X, y)):
            m = clone(mdl).fit(X[tr], y[tr]); p[te] = m.predict_proba(X[te])[:, 1]
            if len(set(y[te])) > 1:
                per.append(dict(seed=sd, roc_auc=roc_auc_score(y[te], p[te]), pr_auc=average_precision_score(y[te], p[te]),
                                balanced_accuracy=balanced_accuracy_score(y[te], (p[te] >= .5).astype(int)), n_test=len(te), studies_in_test=len(set(g[te]))))
        oofs.append(p)
    return pd.DataFrame(per), oofs[0]
# Pooling out-of-fold predictions across leave-one-study-out folds is biased when most held-out studies contain one class
# (a constant prior predictor then scores AUC ~ 0). Primary validation therefore uses stratified *group* 4-fold CV
# (whole studies held out, both classes in every test fold), repeated with 10 seeds; metrics are averaged over folds.
oof = {}; rowsA = []; perfold = []
for nm, mdl in models.items():
    for scheme, grouped in [("Stratified group 4-fold (studies held out), 10 repeats", True), ("Stratified 4-fold ignoring study (leakage-prone), 10 repeats", False)]:
        per, p0 = repeated_cv(mdl, XA.values, yA, gA, grouped=grouped)
        per["model"] = nm; per["validation"] = scheme; perfold.append(per)
        rowsA.append(dict(model=nm, validation=scheme, n=len(yA), folds_evaluated=len(per), roc_auc_mean=per.roc_auc.mean(), roc_auc_sd=per.roc_auc.std(),
                          roc_auc_p2_5=per.roc_auc.quantile(.025), roc_auc_p97_5=per.roc_auc.quantile(.975), pr_auc_mean=per.pr_auc.mean(),
                          balanced_accuracy_mean=per.balanced_accuracy.mean(), **{f"pooled_seed0_{k}": v for k, v in clf_metrics(yA, p0).items()}))
        if grouped: oof[nm] = p0
resA = pd.DataFrame(rowsA)
save_table(resA.pipe(rnd, 4), "M02_modelA_performance", "Model A performance: per-fold metrics averaged over 10 repeats of 4-fold CV (grouped by study vs ignoring study); pooled threshold metrics from seed 0", NB, sheet="Predictive_Models")
save_table(pd.concat(perfold).pipe(rnd, 4), "M02b_modelA_per_fold", "Model A per-fold metrics (all repeats)", NB, sheet="Predictive_Models")
for _, r in resA.iterrows():
    record_test(NB, f"Model A: {r.model} [{r.validation}]", "Classification performance (mean per-fold ROC-AUC)", int(r.n), r.roc_auc_mean, None, np.nan, "ROC-AUC", r.roc_auc_mean, r.roc_auc_p2_5, r.roc_auc_p97_5, family="Predictive models", note=f"SD {r.roc_auc_sd:.3f}; balanced accuracy {r.balanced_accuracy_mean:.3f}; interval = 2.5-97.5 percentiles of fold AUCs", objective="Objective 4")
resA.round(3)

# %%
fold = []
for nm in ["Logistic regression (L2)", "Random forest"]:
    for tr, te in logo.split(XA, yA, gA):
        mdl = models[nm]; mdl.fit(XA.values[tr], yA[tr]); p = mdl.predict_proba(XA.values[te])[:, 1]
        fold.append(dict(model=nm, held_out_study=gA[te][0], sample_type=urb.Sample_Type.values[te][0], n=len(te), true_positive_share=yA[te].mean(),
                         accuracy=accuracy_score(yA[te], (p >= .5).astype(int)), mean_predicted_prob=p.mean()))
fold = pd.DataFrame(fold)
save_table(fold.pipe(rnd, 4), "M03_modelA_per_study", "Model A: accuracy on each held-out study (leave-one-study-out)", NB, sheet="Predictive_Models")
rf = models["Random forest"]; imps = []
for tr, te in logo.split(XA, yA, gA):
    if len(set(yA[te])) < 1: continue
    rf.fit(XA.values[tr], yA[tr])
    pi = permutation_importance(rf, XA.values[te], yA[te], scoring="accuracy", n_repeats=20, random_state=SEED)
    imps.append(pd.Series(pi.importances_mean, index=XA.columns) * len(te))
imp = (pd.concat(imps, axis=1).sum(axis=1) / len(yA)).sort_values(ascending=False).rename("mean_accuracy_drop")
rf.fit(XA.values, yA); gini = pd.Series(rf[-1].feature_importances_, index=XA.columns, name="impurity_importance_full_fit")
lr = models["Logistic regression (L2)"].fit(XA.values, yA); coef = pd.Series(lr[-1].coef_[0], index=XA.columns, name="logistic_coefficient_standardised")
impt = pd.concat([imp, gini, coef], axis=1).sort_values("mean_accuracy_drop", ascending=False)
save_table(impt.pipe(rnd, 5).reset_index().rename(columns={"index": "feature (CLR class share)"}), "M04_modelA_feature_importance", "Model A feature importance: held-out permutation importance (RF), impurity importance and standardised logistic coefficients", NB, sheet="Predictive_Models")
impt.round(3)

# %%
best = "Random forest"
fig, axes = plt.subplots(2, 2, figsize=(6.9, 5.8)); axes = axes.ravel()
for nm, c in [("Logistic regression (L2)", "#4C72B0"), ("Random forest", "#C44E52"), ("Gradient boosting (HistGB)", "#55A868")]:
    fpr, tpr, _ = roc_curve(yA, oof[nm]); pr, rc, _ = precision_recall_curve(yA, oof[nm])
    r = resA[(resA.model == nm) & resA.validation.str.startswith("Stratified group")].iloc[0]
    axes[0].plot(fpr, tpr, color=c, lw=1.2, label=f"{nm.split(' (')[0]} fold-mean AUC={r.roc_auc_mean:.2f}"); axes[1].plot(rc, pr, color=c, lw=1.2, label=f"{nm.split(' (')[0]} fold-mean AP={r.pr_auc_mean:.2f}")
axes[0].plot([0, 1], [0, 1], "k--", lw=.6); axes[0].set_xlabel("False positive rate"); axes[0].set_ylabel("True positive rate"); axes[0].set_title("A  ROC (grouped CV, seed 0 out-of-fold)", fontsize=7.5); axes[0].legend(frameon=False, fontsize=6.2)
axes[1].axhline(yA.mean(), color="k", ls="--", lw=.6); axes[1].set_xlabel("Recall"); axes[1].set_ylabel("Precision"); axes[1].set_title("B  Precision-recall", fontsize=7.5); axes[1].legend(frameon=False, fontsize=6.2)
cm = confusion_matrix(yA, (oof[best] >= .5).astype(int))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=False, ax=axes[2], xticklabels=["Other matrix", "Wastewater"], yticklabels=["Other matrix", "Wastewater"])
axes[2].set_xlabel("Predicted"); axes[2].set_ylabel("Observed"); axes[2].set_title(f"C  Confusion matrix ({best})", fontsize=7.5)
fp_, mp_ = calibration_curve(yA, oof[best], n_bins=8, strategy="quantile")
axes[3].plot(mp_, fp_, "o-", color="#C44E52", ms=3); axes[3].plot([0, 1], [0, 1], "k--", lw=.6)
axes[3].set_xlabel("Mean predicted probability"); axes[3].set_ylabel("Observed fraction"); axes[3].set_title("D  Calibration (RF, quantile bins)", fontsize=7.5)
fig.tight_layout()
save_fig(fig, 27, "modelA_classification_performance", "Exploratory classification of municipal wastewater versus other urban matrices from resistome class composition",
         f"Out-of-fold predictions from stratified group 4-fold cross-validation in which whole studies are held out ({len(yA)} units, {len(set(gA))} studies; curves from one repeat, legend values = mean over 10 repeats x 4 folds). (A) ROC curves; (B) precision-recall curves (dashed = prevalence); (C) confusion matrix at threshold 0.5 for the random forest; (D) calibration. Study-ignorant CV (leakage-prone) is reported in Table M02 for comparison. This is a matrix-signature test, not a prediction of health or food-safety outcomes.",
         NB, "Objective 4", "tables/M02-M04")

# %% [markdown]
# ## Model B — total normalized ARG abundance from prokaryotic phylum composition (DS02)

# %%
m2 = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID")
ph = pd.read_csv(DER / "ds02_prokaryote_phylum_counts.csv", index_col=0).loc[m2.index]
s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0).loc[m2.index]
P = closure(ph.values); keep = ((ph > 0).mean() >= .9) & (P.mean(axis=0) >= 1e-3)
XB = pd.DataFrame(clr(closure(ph.loc[:, keep].values + 0.5)), index=ph.index, columns=ph.columns[keep])
yB = s2.log10_total_RPKM.values
date = pd.to_datetime(m2.Collection_Date + "-01"); cut = pd.Timestamp("2023-09-01")
tr_t, te_t = np.where(date < cut)[0], np.where(date >= cut)[0]
regs = {"Baseline (training mean)": DummyRegressor(), "Ridge regression (CV alpha)": make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 25))),
        "Random forest": RandomForestRegressor(n_estimators=500, min_samples_leaf=3, random_state=SEED, n_jobs=-1)}
def reg_metrics(y, p): return dict(R2=r2_score(y, p), MAE=mean_absolute_error(y, p), RMSE=np.sqrt(mean_squared_error(y, p)), spearman_rho=stats.spearmanr(y, p)[0])
def boot_reg(y, p, n=2000, seed=SEED):
    rng = np.random.default_rng(seed); b = []
    for _ in range(n):
        ix = rng.integers(0, len(y), len(y)); b.append(r2_score(y[ix], p[ix]))
    return np.percentile(b, [2.5, 97.5])
rowsB, predB = [], {}
for nm, mdl in regs.items():
    mdl.fit(XB.values[tr_t], yB[tr_t]); p = mdl.predict(XB.values[te_t]); lo, hi = boot_reg(yB[te_t], p)
    rowsB.append(dict(model=nm, validation=f"Temporal hold-out (train < {cut:%b %Y}, n={len(tr_t)}; test >= {cut:%b %Y}, n={len(te_t)})", n_test=len(te_t), **reg_metrics(yB[te_t], p), R2_ci_low=lo, R2_ci_high=hi))
    predB[("temporal", nm)] = p
    p = cross_val_predict(mdl, XB.values, yB, groups=m2.City.values, cv=LeaveOneGroupOut())
    lo, hi = boot_reg(yB, p)
    rowsB.append(dict(model=nm, validation="Leave-one-city-out (pooled out-of-fold)", n_test=len(yB), **reg_metrics(yB, p), R2_ci_low=lo, R2_ci_high=hi))
    predB[("city", nm)] = p
    for c in m2.City.unique():
        ix = (m2.City == c).values
        rowsB.append(dict(model=nm, validation=f"Leave-one-city-out: held-out {c}", n_test=int(ix.sum()), **reg_metrics(yB[ix], p[ix]), R2_ci_low=np.nan, R2_ci_high=np.nan))
resB = pd.DataFrame(rowsB)
save_table(resB.pipe(rnd, 4), "M05_modelB_performance", f"Model B: predicting log10 total ARG RPKM from CLR prokaryotic phylum composition ({XB.shape[1]} phyla); R2 95% CI by bootstrap", NB, sheet="Predictive_Models")
for _, r in resB[~resB.validation.str.contains("held-out")].iterrows():
    record_test(NB, f"Model B: {r.model} [{r.validation[:40]}]", "Regression performance", int(r.n_test), r.R2, None, np.nan, "R2 (out-of-sample)", r.R2, r.R2_ci_low, r.R2_ci_high, family="Predictive models", note=f"MAE {r.MAE:.3f} log10 units", objective="Objective 4")
rfB = regs["Random forest"].fit(XB.values[tr_t], yB[tr_t])
piB = permutation_importance(rfB, XB.values[te_t], yB[te_t], scoring="r2", n_repeats=30, random_state=SEED)
impB = pd.DataFrame({"phylum": XB.columns, "permutation_importance_R2_drop": piB.importances_mean, "sd": piB.importances_std}).sort_values("permutation_importance_R2_drop", ascending=False)
save_table(impB.pipe(rnd, 5), "M06_modelB_feature_importance", "Model B: permutation importance on the temporal test set (random forest)", NB, sheet="Predictive_Models")
resB[~resB.validation.str.contains("held-out")].round(3)

# %%
fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.9))
for ax, key, ttl in [(axes[0], ("temporal", "Random forest"), "A  Temporal hold-out, random forest"), (axes[1], ("temporal", "Ridge regression (CV alpha)"), "B  Temporal hold-out, ridge")]:
    p = predB[key]; y = yB[te_t]; r = resB[(resB.model == key[1]) & resB.validation.str.startswith("Temporal")].iloc[0]
    for c in CITY_COLORS:
        ix = (m2.City.values[te_t] == c); ax.scatter(y[ix], p[ix], s=10, color=CITY_COLORS[c], label=c, alpha=.75, edgecolors="none")
    lim = [min(y.min(), p.min()) - .05, max(y.max(), p.max()) + .05]; ax.plot(lim, lim, "k--", lw=.6); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Observed log10 total ARG RPKM"); ax.set_ylabel("Predicted"); ax.set_title(f"{ttl}\nR² = {r.R2:.2f} (95% CI {r.R2_ci_low:.2f} to {r.R2_ci_high:.2f}), MAE = {r.MAE:.3f}", fontsize=7.2)
axes[0].legend(frameon=False, fontsize=6.2)
t = impB.head(12).iloc[::-1]
axes[2].barh(t.phylum, t.permutation_importance_R2_drop, xerr=t.sd, color="#8172B2"); axes[2].set_xlabel("Drop in test R² when permuted"); axes[2].set_title("C  Permutation importance (RF, temporal test)", fontsize=7.2)
axes[2].tick_params(axis="y", labelsize=7)
fig.tight_layout()
save_fig(fig, 28, "modelB_regression_performance", "Exploratory prediction of total normalized ARG abundance from prokaryotic phylum composition (DS02)",
         f"Models trained on samples collected before {cut:%B %Y} (n = {len(tr_t)}) and tested on later samples (n = {len(te_t)}); leave-one-city-out results are in Table M05. Features: CLR of {XB.shape[1]} prevalent prokaryotic phyla; target: log10 total ARG RPKM. R² CIs from 2,000 bootstrap resamples of the test set. Association-based prediction within one study; not a forecast of health outcomes.",
         NB, "Objective 4", "tables/M05, M06")

# %%
feasP = pd.DataFrame([
    dict(target="Disease incidence / outbreak prediction", required_data="Clinical or epidemiological outcomes matched to samples by place and time", available="None linked", decision="Not feasible"),
    dict(target="Clinical AMR prevalence (e.g., GLASS indicators)", required_data="Country/city-year surveillance indicators with matching denominators", available="Not downloaded; country-level only; temporal mismatch with 2016 pilot", decision="Not feasible in this analysis; future work"),
    dict(target="Food-safety or food-business impact", required_data="Food-chain contamination, recall, trade or economic outcome data", available="No such data available", decision="Not feasible"),
    dict(target="Environmental matrix from resistome (Model A)", required_data="Matrix labels, resistome profiles, >= 2 studies per class", available=f"{len(yA)} units, {len(set(gA))} studies", decision="Performed (exploratory)"),
    dict(target="normalized ARG abundance from taxonomy (Model B)", required_data="Paired taxonomy and resistome from one pipeline; temporal or spatial hold-out", available="DS02: 447 samples, 4 cities, 25 months", decision="Performed (exploratory)"),
    dict(target="Temporal forecasting of normalized ARG abundance", required_data="Long, regularly spaced series per site with site identifiers", available="25 monthly points per city; sites unidentified", decision="Not performed beyond the temporal hold-out; insufficient series length per site")])
# AUD-17 documentation: Model B uses a fixed pseudo-count (0.5) for CLR, which is not estimated from data; its prevalence filter uses no outcome information.
save_table(feasP, "M07_predictive_feasibility", "Feasibility of predictive targets with the available data", NB, sheet="Predictive_Models")
rA = resA.set_index(["model", "validation"])
G_ = "Stratified group 4-fold (studies held out), 10 repeats"; U_ = "Stratified 4-fold ignoring study (leakage-prone), 10 repeats"
la = rA.loc[("Random forest", G_)]; ka = rA.loc[("Random forest", U_)]; lra = rA.loc[("Logistic regression (L2)", G_)]; ba = rA.loc[("Baseline (prior)", G_)]
record_finding(NB, "F17", "Objective 4", f"With whole studies held out, resistome class composition discriminated municipal wastewater from other urban matrices with mean fold ROC-AUC {la.roc_auc_mean:.2f} (SD {la.roc_auc_sd:.2f}; 2.5-97.5 percentile {la.roc_auc_p2_5:.2f}-{la.roc_auc_p97_5:.2f}; random forest) and {lra.roc_auc_mean:.2f} (logistic regression), against {ba.roc_auc_mean:.2f} for the prior baseline; ignoring study in the folds raised the random-forest AUC to {ka.roc_auc_mean:.2f}.",
               len(yA), f"Urban ResFinder units ({len(set(gA))} studies)", "Repeated stratified group k-fold CV (10 x 4 folds)", "Tables M01-M04; Figure 27",
               ("A wastewater resistome signature transfers to unseen studies" if la.roc_auc_p2_5 > .5 else "Transfer of the signature to unseen studies is uncertain (fold AUCs overlap 0.5)")
               + "; the gap between grouped and ungrouped validation quantifies the optimism caused by study leakage.",
               "Few non-wastewater studies; matrix labels are study-level; not a health-outcome model.")
rb = resB.set_index(["model", "validation"])
tB = resB[(resB.model == "Random forest") & resB.validation.str.startswith("Temporal")].iloc[0]; cB = resB[(resB.model == "Random forest") & (resB.validation == "Leave-one-city-out (pooled out-of-fold)")].iloc[0]
tR = resB[(resB.model == "Ridge regression (CV alpha)") & resB.validation.str.startswith("Temporal")].iloc[0]; cR = resB[(resB.model == "Ridge regression (CV alpha)") & (resB.validation == "Leave-one-city-out (pooled out-of-fold)")].iloc[0]
pcR = resB[(resB.model == "Ridge regression (CV alpha)") & resB.validation.str.contains("held-out")]
record_finding(NB, "F18", "Objective 4", f"Prokaryotic phylum composition predicted total normalized ARG abundance in later months with R² = {tR.R2:.2f} (95% CI {tR.R2_ci_low:.2f}-{tR.R2_ci_high:.2f}; ridge) and {tB.R2:.2f} ({tB.R2_ci_low:.2f}-{tB.R2_ci_high:.2f}; random forest); in unseen cities pooled R² was {cR.R2:.2f} (ridge) and {cB.R2:.2f} (random forest), with held-out-city R² ranging {pcR.R2.min():.2f} to {pcR.R2.max():.2f} (ridge). The training-mean baseline has R² <= 0 by construction.",
               len(yB), "DS02 samples", "Temporal hold-out and leave-one-city-out; bootstrap CI", "Tables M05-M06; Figure 28",
               ("Community composition carries information about normalized resistome abundance that transfers to later months" if max(tB.R2_ci_low, tR.R2_ci_low) > 0 else "Transfer of taxonomy-based prediction to later months is not established (R2 CI includes 0)")
               + f"; transfer to unseen cities is inconsistent (negative R² for {', '.join(pcR.validation.str.replace('Leave-one-city-out: held-out ', '')[pcR.R2 < 0]) or 'no city'})"
               + ". This is not a forecast of human AMR.",
               "Single study and pipeline; phylum resolution; sites unidentified; no external validation dataset.")
