from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired

class SearchForm(FlaskForm):
    searchfield = StringField("Article name or doi: ", validators=[DataRequired()])
    submit = SubmitField("Submit")
