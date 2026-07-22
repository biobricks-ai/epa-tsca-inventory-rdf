import importlib.util
from pathlib import Path
import pandas as pd
from rdflib import Literal,RDF
S=importlib.util.spec_from_file_location("convert",Path("stages/01_convert.py"));M=importlib.util.module_from_spec(S);S.loader.exec_module(M)

def test_nhexane_standard_identifiers():
    frame=pd.DataFrame([{"casrn":"110-54-3","chemical_name":"Hexane","smiles":"CCCCCC","inchi":"InChI=1S/C6H14/c1-3-5-6-4-2/h3-6H2,1-2H3","activity_status":"ACTIVE"}])
    g,c=M.convert_table("tsca_inventory",frame);s=M.BASE["substance/cas/110-54-3"];hub=M.COMPOUND["VLKZOEOYAKHREP-UHFFFAOYSA-N"]
    assert (s,M.SKOS.notation,Literal("110-54-3")) in g;assert (s,M.CHEMINF_IDENTITY,hub) in g;assert c["row_coverage"]==c["field_coverage"]==1

def test_pmn_entry_is_not_given_guessed_compound_identity():
    frame=pd.DataFrame([{"id":1,"pmnno":"P000005","accno":232689,"genericname":"Confidential polymer","activity":"ACTIVE"}])
    g,c=M.convert_table("tsca_pmn_accession",frame);subjects=set(g.subjects(RDF.type,M.CHEMICAL));assert len(subjects)==1;assert not list(g.objects(next(iter(subjects)),M.CHEMINF_IDENTITY));assert c["field_coverage"]==1

def test_invalid_structure_retains_cas_entry():
    frame=pd.DataFrame([{"casrn":"123-45-6","chemical_name":"Unknown","smiles":"bad","activity_status":"ACTIVE"}])
    g,_=M.convert_table("tsca_inventory",frame);s=M.BASE["substance/cas/123-45-6"];assert (s,RDF.type,M.CHEMICAL) in g;assert not list(g.objects(s,M.CHEMINF_IDENTITY))

def test_model_has_labels_provenance_and_no_blank_nodes():
    frame=pd.DataFrame([{"casrn":"110-54-3","chemical_name":"Hexane","smiles":"CCCCCC","activity_status":"ACTIVE"}]);g,_=M.convert_table("tsca_inventory",frame)
    for s in g.subjects(RDF.type,M.CHEMICAL):assert g.value(s,M.SCHEMA_NAME);assert (s,M.PROV.wasDerivedFrom,M.SOURCE) in g
    assert not any(type(x).__name__=="BNode" for t in g for x in t)
