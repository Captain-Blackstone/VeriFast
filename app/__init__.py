from flask import Flask
from config import Config
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from redis import Redis, StrictRedis
import redis
import rq
import datetime
import os


def start_redis():
    redtest = StrictRedis(host="localhost", port=6379)
    try:
        return redtest.ping()
    except redis.exceptions.ConnectionError:
        os.system("redis-server &")
        os.system("rq worker download-articles-tasks &")


#start_redis()
app = Flask(__name__)
app.config.from_object(Config)
db = SQLAlchemy(app)
migrate = Migrate(app, db)
app.redis = Redis.from_url(app.config['REDIS_URL'])
app.task_queue = rq.Queue('download-articles-tasks', connection=app.redis)
app.session_time = datetime.datetime.now().strftime("%Y.%m.%d_%H:%M:%S")

from app import routes
from app import models
