"""06a: export the tables quoted in the QC report from the built workbooks and
harmonized CSVs to harmonized/report_data.json (consumed by 06_build_report.js)."""
import json, os
import pandas as pd
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
O = lambda f: os.path.join(B, "outputs", f)
I, X = O("01_Public_Dataset_Inventory.xlsx"), O("00_Master_Data_Index.xlsx")
q = pd.read_excel(I, "Search_Query_Log")
qs = q.groupby("Database").agg(Queries=("Query_ID", "count"), Returned=("N_Results_Returned", "sum"), Inspected=("N_Records_Inspected", "sum"),
                               Relevant=("N_Potentially_Relevant", "sum"), Verified=("N_Verified", "sum"),
                               Retained=("N_Included_or_Contextual_After_Review", "sum")).reset_index()
inv = pd.read_excel(I, "Study_Inventory")
invs = inv.groupby(["Record_Type", "Decision"]).size().unstack(fill_value=0).reset_index()
sm = pd.read_csv(os.path.join(B, "harmonized", "sample_master.csv"), dtype=str)
inc = sm[sm.Inclusion_Status == "Included"]
core = [["ID", "Dataset", "Accession(s)", "Samples", "Included/contextual", "Cities", "Samples with AMR table"]]
for ds in ["DS01", "DS02", "DS03", "DS04", "DS05"]:
    x = sm[sm.Dataset_ID == ds]; i = x[x.Inclusion_Status.isin(["Included", "Separate contextual dataset"])]
    core.append([ds, x.Dataset_Name.iloc[0], "; ".join(sorted(set(x.Study_ID)))[:90], str(len(x)), str(len(i)),
                 str(i[i.Geographic_Unit == "City (source-labelled)"].City.nunique()),
                 "; ".join(f"{lab}: {n}" for lab, n in [("CARD/RGI", int(i.AMR_Data_Available.str.contains("CARD").sum())),
                                                         ("ResFinder/KMA", int(i.AMR_Data_Available.str.contains("ResFinder").sum()))] if n) or "None (raw reads only)"])
st = pd.read_excel(X, "Sample_Type_Coverage"); rd = pd.read_excel(X, "Analysis_Readiness"); mm = pd.read_excel(X, "Missing_Metadata")
qc = pd.read_excel(O("02_Sample_Metadata_Master.xlsx"), "QC_Report"); ver = pd.read_excel(I, "Source_Verification")
summ = pd.read_excel(X, "Summary")
sa = inc[inc.Country.isin(["India", "Bangladesh", "Pakistan", "Nepal", "Sri Lanka"])].groupby("Country").agg(
    Included_samples=("Sample_ID", "count"), Datasets=("Dataset_ID", "nunique")).reset_index()
ex = sm[sm.Inclusion_Status == "Excluded"].Inclusion_Reason.value_counts().reset_index()
ex.columns = ["Exclusion reason", "Samples"]
def tl(df): return [list(map(str, df.columns))] + [[str(v) for v in r] for r in df.values.tolist()]
json.dump({"qs": tl(qs), "invs": tl(invs), "core": core, "st": tl(st), "rd": tl(rd[["Proposed_Analysis", "Readiness", "Reason"]]),
           "mm": tl(mm), "qc": tl(qc[["Check_ID", "Check", "Result", "N_Records_Flagged"]]), "ver": tl(ver[["Check_ID", "Item", "Finding", "Status"]]),
           "summ": dict(zip(summ.Metric, summ.Value.astype(str))), "sa": tl(sa), "ex": tl(ex)},
          open(os.path.join(B, "harmonized", "report_data.json"), "w"), indent=1)
print("exported")
