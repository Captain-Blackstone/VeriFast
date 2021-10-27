from flask import render_template, redirect, url_for, jsonify
from app import app
from app.forms import SearchForm
from backend import load_article, paragraphs_from_article
from config import Config


@app.route("/", methods=["GET", "POST"])
@app.route("/search", methods=["GET", "POST"])
def search():
    form = SearchForm()
    if form.validate_on_submit():
        doi_or_name = form.searchfield.data
        return redirect(url_for("view", doi_or_name=doi_or_name))
    return render_template("search.html", title="Search", form=form)


@app.route("/view_article/<doi_or_name>")
def view(doi_or_name):
    article, citations, article_id = load_article(doi_or_name)
    return render_template("view_article.html", article=article,
                           citations=citations,
                           article_id=article_id)


@app.route("/search_paragraphs/<citation_id>/<article>")
def search_paragraphs(citation_id, article):
    paragraphs_from_article(*citation_id, article)
    return None


@app.route("/get_paragraphs/<isection>/<iparagraph>/<isentence>/<citation_instance>/<article_id>")
async def get_paragraphs(isection, iparagraph, isentence, citation_instance, article_id):
    res = await paragraphs_from_article(isection, iparagraph, isentence, citation_instance, article_id)
    return jsonify(res)


@app.route("/contact")
def contact():
    return render_template("contacts.html",
                           mail_for_feedback=Config.feedback_email,
                           mail_for_collaborators=Config.collaboration_email)
