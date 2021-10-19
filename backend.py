import json
from collections import defaultdict

from app.models import Article
from sqlalchemy.orm.exc import NoResultFound
from app import db

import time

from app import SearcherObj
from parser import Parser
from nn import SemanticSearcher


def load_article(doi_or_name=None):
    doi_or_name = " ".join(doi_or_name.split())  # remove multiple spaces
    if len(doi_or_name) < 32:
        search_parameter = Article.doi
    else:
        search_parameter = Article.title
    try:
        print("trying", doi_or_name)
        record = db.session.query(Article).filter(search_parameter == doi_or_name).one()
        with open(f"article_base/{record.filename}", "r") as fl:
            json_file = json.load(fl)
            article_id = record.id
    except NoResultFound:
        start = time.time()
        xml_file_text, source = SearcherObj(doi_or_name)
        with open("test.xml", "w") as fl:
            fl.write(xml_file_text)
        end = time.time()
        print(end-start)
        print(source)
        json_file = Parser()(xml_file_text, source, parse_citations=True)
        doi, title = json_file.get("doi", None), json_file.get("full_title", None)
        new_article = Article(doi=doi, title=title)
        if len(list(db.session.query(Article).filter(Article.doi == doi))) == 0:
            db.session.add(new_article)
            db.session.flush()
            new_article.filename = f"{new_article.id}.json"
            with open(f"article_base/{new_article.filename}", "w") as fl:
                json.dump(json_file, fl)
            db.session.commit()
        article_id = new_article.id

    citations_dct = defaultdict(lambda: [])
    for citation in json_file["citations"]["sentence_ids"]:
        key = tuple(citation[:3])
        value = citation[3]
        citations_dct[key].append(value)
    return json_file, citations_dct, article_id


def paragraphs_from_article(isection, iparagraph, isentence, citation_instance, article_id):
    print("button pressed")
    isection, iparagraph, isentence = int(isection), int(iparagraph), int(isentence)
    with open(f"article_base/{article_id}.json", "r") as fl:
        article = json.load(fl)
    sentence = article["text"][isection]["section_text"][iparagraph][isentence]

    # Only for demonstration
    with open(f"article_base/{article_id}.json", "r") as fl:
        source_article = json.load(fl)
    # source_article = article["citations"]["papers"][citation_instance]

    if source_article["text"]:
        out_dict = dict()
        out_dict["full_title"] = source_article["full_title"]
        return SemanticSearcher().dirty_call(sentence, source_article, out_dict)
    query = None
    if source_article["doi"]:
        query = source_article["doi"]
    elif source_article["full_title"]:
        query = source_article["full_title"]
    if query is not None:
        print("starting search")
        source_article, source_db = SearcherObj(query)
        print("starting parsing")
        source_article_json = Parser()(source_article, source_db)
        article["citations"]["papers"][citation_instance] = source_article_json
        # with open(f"article_base/{article_id}.json", "r") as fl:
        #    json.dump(article, fl)
        print("starting semantic search")
        out_dict = dict()
        out_dict["full_title"] = source_article_json["full_title"]
        return SemanticSearcher().dirty_call(sentence, source_article_json, out_dict)

    
