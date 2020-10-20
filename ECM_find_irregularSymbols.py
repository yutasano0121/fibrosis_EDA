import mygene
import csv
import pandas as pd

d = pd.read_csv("ECM_interaction.csv")
genes1 = list(set(d["Gene.symbol.interactor.A"]))
genes2 = list(set(d["Gene.symbol.interactor.B"]))

mg = mygene.MyGeneInfo()
query1 = mg.querymany(
    genes1, scopes="symbol", 
    fields=["ensembl.gene"], species="human", 
    as_dataframe=True
)
query2 = mg.querymany(
    genes2, scopes="symbol", 
    fields=["ensembl.gene"], species="human", 
    as_dataframe=True
)

symbol1 = list(set(query1.loc[query1["notfound"] == True].index))
symbol2 = list(set(query2.loc[query2["notfound"] == True].index))

with open("ECM_irregularSymbolsA.csv", "w") as o:
    csv_o = csv.writer(o)
    for i in symbol1:
        csv_o.writerow([i]])
with open("ECM_irregularSymbolsB.csv", "w") as o:
    csv_o = csv.writer(o)
    for i in symbol2:
        csv_o.writerow([i]])


