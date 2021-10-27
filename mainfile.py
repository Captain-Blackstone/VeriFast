from app import app, db
from app.models import Article, DoiPmcidPmid

from searcher import Searcher
import json


with open("config.json", "r") as config:
    config_dict = json.load(config)
    pubmed_email_config = config_dict["pubmed_email"]
    elsevier_apikey_config = config_dict["elsevier_apikey"]

SearcherObj = Searcher(doi_pmid_pmcid_db=db,
                       api_key="",
                       pubmed_email=pubmed_email_config,
                       elsevier_apikey_config=elsevier_apikey_config,
                       pmc_db="pmc",
                       pubmed_db="pubmed",
                       doi_pmcid_pmid_class=DoiPmcidPmid)


@app.shell_context_processor
def make_shell_context():
    return {"db": db, "Article": Article, "DoiPmcidPmid": DoiPmcidPmid}
