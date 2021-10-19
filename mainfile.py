from app import app, db
from app.models import Article, DoiPmcidPmid


@app.shell_context_processor
def make_shell_context():
    return {"db": db, "Article": Article, "DoiPmcidPmid": DoiPmcidPmid}

