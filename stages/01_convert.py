#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,os
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq
from rdflib import Graph,Literal,Namespace,RDF,RDFS,URIRef
from rdflib.namespace import DCTERMS,PROV,SKOS
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

BASE=Namespace("https://biobricks.ai/epa-tsca-inventory-rdf/")
BB=Namespace("https://biobricks.ai/ontology/")
COMPOUND=Namespace("https://biobricks.ai/compound/")
CHEMICAL=URIRef("https://schema.org/ChemicalSubstance")
SCHEMA_NAME=URIRef("https://schema.org/name")
SCHEMA_IDENTIFIER=URIRef("https://schema.org/identifier")
CHEMINF_IDENTITY=URIRef("http://semanticscience.org/resource/CHEMINF_000477")
SOURCE=URIRef("https://www.epa.gov/tsca-inventory")
FIELDS={"tsca_inventory":("casrn","chemical_name","smiles","inchi","activity_status","flag","uvcb","definition"),"tsca_pmn_accession":("id","pmnno","accno","uid","exp","genericname","flag","activity")}
PRED={"activity_status":BB.activityStatus,"activity":BB.activityStatus,"flag":BB.inventoryFlag,"uvcb":BB.uvcb,"definition":DCTERMS.description,"exp":BB.expiration,"smiles":BB.canonicalSmiles,"inchi":BB.inchi,"pmnno":SCHEMA_IDENTIFIER,"accno":SCHEMA_IDENTIFIER,"uid":SCHEMA_IDENTIFIER,"id":SCHEMA_IDENTIFIER}

def text(v):
    if v is None or pd.isna(v): return None
    v=str(v).strip(); return v or None

def identity(smiles):
    mol=Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:return None
    key=Chem.MolToInchiKey(mol)
    return (Chem.MolToSmiles(mol,canonical=True),key) if key else None

def convert_table(name,frame,offset=0):
    graph=Graph(); stats={"source_rows":len(frame),"represented_rows":0,"nonempty_fields":0,"preserved_fields":0,"triples":0}
    for pos,(_,row) in enumerate(frame.iterrows()):
        if name=="tsca_inventory":
            cas=text(row.casrn); subject=BASE["substance/cas/"+cas]
            graph.add((subject,RDF.type,CHEMICAL)); graph.add((subject,SKOS.notation,Literal(cas)))
            found=identity(text(row.smiles))
            if found:
                hub=COMPOUND[found[1]];graph.add((subject,CHEMINF_IDENTITY,hub));graph.add((hub,DCTERMS.identifier,Literal(found[1])))
            label=text(row.chemical_name)
        else:
            token=text(row.pmnno) or text(row.accno) or str(offset+pos)
            subject=BASE["substance/pmn/"+hashlib.sha256(token.encode()).hexdigest()[:20]]
            graph.add((subject,RDF.type,CHEMICAL));label=text(row.genericname)
        if label:graph.add((subject,SCHEMA_NAME,Literal(label)))
        graph.add((subject,PROV.wasDerivedFrom,SOURCE))
        for field in FIELDS[name]:
            value=text(row.get(field))
            if value is None:continue
            stats["nonempty_fields"]+=1
            predicate=SKOS.notation if field=="casrn" else SCHEMA_NAME if field in {"chemical_name","genericname"} else PRED.get(field,BB[field])
            graph.add((subject,predicate,Literal(value)));stats["preserved_fields"]+=1
        stats["represented_rows"]+=1
    stats["triples"]=len(graph);stats["row_coverage"]=stats["represented_rows"]/stats["source_rows"] if stats["source_rows"] else 1;stats["field_coverage"]=stats["preserved_fields"]/stats["nonempty_fields"] if stats["nonempty_fields"] else 1
    return graph,stats

def root():
    for p in (os.getenv("EPA_TSCA_BRICK"),".bb/dependencies/epa-tsca-inventory/brick","/mnt/raid2/biobricks/epa-tsca-inventory/brick"):
        if p and Path(p).is_dir():return Path(p)
    raise FileNotFoundError("run biobricks pull or set EPA_TSCA_BRICK")

def main():
    out=Path("brick/epa-tsca-inventory-rdf.nt");out.parent.mkdir(exist_ok=True)
    report={"tables":{},"source_rows":0,"represented_rows":0,"nonempty_fields":0,"preserved_fields":0,"triples":0}
    with out.open("w") as stream:
        for name in FIELDS:
            totals={k:0 for k in ("source_rows","represented_rows","nonempty_fields","preserved_fields","triples")};offset=0
            for batch in pq.ParquetFile(root()/f"{name}.parquet").iter_batches(batch_size=3000):
                frame=batch.to_pandas();g,part=convert_table(name,frame,offset);stream.write(g.serialize(format="nt"));offset+=len(frame)
                for k in totals:totals[k]+=part[k]
            totals["row_coverage"]=totals["represented_rows"]/totals["source_rows"];totals["field_coverage"]=totals["preserved_fields"]/totals["nonempty_fields"];report["tables"][name]=totals
            for k in totals:
                if k in report:report[k]+=totals[k]
    report["row_coverage"]=report["represented_rows"]/report["source_rows"];report["field_coverage"]=report["preserved_fields"]/report["nonempty_fields"]
    Path("brick/coverage.json").write_text(json.dumps(report,indent=2)+"\n");print(json.dumps(report,indent=2))
if __name__=="__main__":main()
