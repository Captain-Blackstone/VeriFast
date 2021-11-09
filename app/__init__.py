from flask import Flask
from config import Config
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from redis import Redis, StrictRedis
import redis
import rq
import datetime
import os


app = Flask(__name__)
app.config.from_object(Config)
db = SQLAlchemy(app)
migrate = Migrate(app, db)
app.redis = Redis.from_url(app.config['REDIS_URL'])
app.task_queue = rq.Queue('download-articles-tasks', connection=app.redis)
app.session_time = datetime.datetime.now().strftime("%Y.%m.%d_%H:%M:%S")

from app import routes
from app import models
