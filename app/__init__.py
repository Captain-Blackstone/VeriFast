from flask import Flask
from config import Config
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from searcher import Searcher
import json

app = Flask(__name__)
app.config.from_object(Config)
db = SQLAlchemy(app)
migrate = Migrate(app, db)

with open("config.json", "r") as config:
    config_dict = json.load(config)
    pubmed_email_config = config_dict["pubmed_email"]
    elsevier_apikey_config = config_dict["elsevier_apikey"]

#with open("PMID_PMCID_DOI.csv", "r") as fl:
    # for demonstration only
#    text = [next(fl) for x in range(10)]
#    text = fl.readlines()



SearcherObj = Searcher(api_key="",
                       pubmed_email=pubmed_email_config,
                       elsevier_apikey_config=elsevier_apikey_config,
 #                      pmid_pmcid_doi_db_text=text,
                       pmc_db="pmc",
                       pubmed_db="pubmed")

from app import routes, models
