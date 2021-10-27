from app import app, db
import rq
from rq.exceptions import NoSuchJobError
from redis.exceptions import RedisError


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

    def __repr__(self):
        return f"DoiPmcidPmid {self.doi} {self.pmid} {self.pmcid}"

class Task(db.Model):
    id = db.Column(db.String(36), primary_key=True)
    doi_or_name = db.Column(db.String(1000), index=True)
    log_file = db.Column(db.String(1000))
    complete = db.Column(db.Boolean, default=False)

    def get_rq_job(self):
        try:
            rq_job = rq.job.Job.fetch(self.id, connection=app.redis)
        except (RedisError, NoSuchJobError):
            return None
        return rq_job

    def get_task_results(self):
        job = self.get_rq_job()
        if job.meta["complete"]:
            return [True, job.meta["json"], job.meta["article_id"]]
        elif job.__dict__["exc_info"] is not None:
            exception_text = "\t".join(job.__dict__["exc_info"].split("\n"))
            return [True, exception_text]
        else:
            return [False, None]
