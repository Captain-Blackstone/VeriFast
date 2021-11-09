import json
from collections import defaultdict

from app import app
from app.models import Task
from app import db
from nn import SemanticSearcher


def launch_background_task(task_name, high_priority=False, *args, **kwargs):
    rq_job = app.task_queue.enqueue('app.Background_tasks.' + task_name, at_front=high_priority, *args, **kwargs)
    rq_job.meta["complete"] = False
    rq_job.save_meta()
    task = Task(id=rq_job.get_id(), *args, **kwargs)
    db.session.add(task)
    return task


def get_task_results(task):
    js, art_id = None, None
    while not task.complete:
        task_report = task.get_task_results()
        task.complete = task_report[0]
        if len(task_report) == 3:  # the task is complete without errors
            js, art_id = task_report[1:]
        elif task_report[1] is not None:  # the task failed with exception
            js, art_id = None, None
    return js, art_id


def load_article(doi_or_name=None):
    task_to_find_article = launch_background_task("download_article",
                                                  doi_or_name=doi_or_name,
                                                  log_file=f"logs/{app.session_time}_failed_tasks.log")
    json_file, article_id = get_task_results(task_to_find_article)

    if json_file is None:
        return show_error(), None, None

    # This is needed for citations rendering as buttons
    citations_dct = defaultdict(lambda: [])
    for citation in json_file["citations"]["sentence_ids"]:
        key = tuple(citation[:3])
        value = citation[3]
        citations_dct[key].append(value)

    # Send to background downloading of all the cited articles
    #for val in json_file["citations"]["papers"].values():
    #    cited_paper_name = val["full_title"]
    #    launch_background_task("download_article",
    #                           doi_or_name=cited_paper_name,
    #                           log_file=f"logs/{app.session_time}_failed_tasks.log")
    return json_file, citations_dct, article_id


async def paragraphs_from_article(isection, iparagraph, isentence, citation_instance, article_id):
    print("button pressed")
    isection, iparagraph, isentence = int(isection), int(iparagraph), int(isentence)
    with open(f"article_base/{article_id}.json", "r") as fl:
        article = json.load(fl)
    sentence = article["text"][isection]["section_text"][iparagraph][isentence]

    source_article = article["citations"]["papers"][citation_instance]

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
        task_to_find_cited_article = launch_background_task("download_article",
                                                            high_priority=True,
                                                            doi_or_name=query,
                                                            log_file=f"logs/{app.session_time}_failed_tasks.log")
        source_article_json, article_id = get_task_results(task_to_find_cited_article)
        if source_article_json is None:
            return show_error()

        print("starting semantic search")
        out_dict = dict()
        out_dict["full_title"] = source_article_json["full_title"]
        return SemanticSearcher().dirty_call(sentence, source_article_json, out_dict)
    else:
        return show_error()


def show_error():
    with open("article_base/error.json", "r") as fl:
        json_file = json.load(fl)
    return json_file
