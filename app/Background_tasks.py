import json
import time
import traceback
from sqlalchemy.orm.exc import NoResultFound

from app import db
from app.models import Article
from mainfile import SearcherObj
from aparser import Parser
from rq import get_current_job

async def download_article(doi_or_name, log_file):
    doi_or_name = " ".join(doi_or_name.split())  # remove multiple spaces
    if len(doi_or_name) < 32:
        search_parameter = Article.doi
    else:
        search_parameter = Article.title
    try:  # If the article is already downloaded, we just find it locally
        print(f"trying to find {doi_or_name} in the db")
        record = db.session.query(Article).filter(search_parameter == doi_or_name).one()
        with open(f"article_base/{record.filename}", "r") as fl:
            json_file = json.load(fl)
        print(f"found {doi_or_name} in the db")
        job = get_current_job()
        job.meta["json"] = json_file
        job.meta["article_id"] = record.id
        job.meta["complete"] = True
        job.save_meta()
    except NoResultFound:  # If not, we download it
        print("article is not in the database")
    try:
        # Find the article
        print("starting web search")
        xml_file_text, source = await SearcherObj(doi_or_name)
        # time.sleep(2)

        # Parse the article
        print("starting parsing")
        json_file = Parser()(xml_file_text, source, parse_citations=True)
        # time.sleep(2)

        # Add the article to db
        print("adding to db")
        doi, title = json_file.get("doi", None), json_file.get("full_title", None)
        new_article = Article(doi=doi, title=title)
        if len(list(db.session.query(Article).filter(Article.doi == doi))) == 0:
            db.session.add(new_article)
            db.session.flush()
            new_article.filename = f"{new_article.id}.json"
            with open(f"article_base/{new_article.filename}", "w") as fl:
                json.dump(json_file, fl)
            db.session.commit()
        # time.sleep(2)
        print("finished")
        job = get_current_job()
        job.meta["json"] = json_file
        job.meta["article_id"] = new_article.id
        job.meta["complete"] = True
        job.save_meta()
    except Exception:
        exception_text = "|".join(traceback.format_exc().split("\n"))
        with open(log_file, "a") as fl:
            fl.write(f"'{doi_or_name}';{exception_text}\n")

