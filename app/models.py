from app import db

class Article(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    doi = db.Column(db.String(32), index=True, unique=True)
    title = db.Column(db.String(1000), index=True, unique=True)
    filename = db.Column(db.String(50), index=True, unique=True)

    def __repr__(self):
        return "<Article {}>".format(self.title)

class DoiPmcidPmid(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    doi = db.Column(db.String(64), index=True)
    pmid = db.Column(db.Integer, index=True)
    pmcid = db.Column(db.String(32), index=True)
