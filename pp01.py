from flask import Flask, render_template, url_for, redirect
from flask_bootstrap import Bootstrap
from flask_fontawesome import FontAwesome
import pandas as pd
import extract

from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired


# create a flask application object
def create_app():
    app = Flask(__name__)
    Bootstrap(app)
    FontAwesome(app)
    app.config['SECRET_KEY'] = 'change this unsecure key'

    return app
app = create_app()
# we need to set a secret key attribute for secure forms

extract.get_important_data()

class QueryForm(FlaskForm):
    tf_name = StringField('Enter a valid gene symbol for a human transcription factor: ', validators=[DataRequired(message='Click here to enter data')])
    submit = SubmitField('Submit')

# define the action for the top level route
@app.route('/', methods=['GET','POST'])
def index():
    form = QueryForm()
    tf_name = None

    if form.validate_on_submit():
        tf_name = form.tf_name.data
        return redirect(url_for('background', tf_name = str(tf_name)))
    return render_template('drugprofile.html', form=form, tf_name = tf_name)

@app.route('/drug_profile')
def drug():
    return 'drug task goes here'

@app.route('/background_info/<tf_name>', methods=['GET'])
def background(tf_name):
    data = extract.all_data_to_dict(tf_name)
    ensembl = data['Ensembl']
    symbol = data['Symbol']
    family = data['Gene_family']
    chr = data["Chr_location"]
    full_name = data['Full_name']
    uniprot = data['Uniprot']
    subcell = data['Subcellular_loc']
    func = data['Function']

    return render_template('page2.html', symbol= symbol, ensembl=ensembl , family=family , chr=chr, full_name=full_name, uniprot=uniprot,
    subcell=subcell, func=func)




@app.route('/GEO')
def geo():
    return '<h1>' + 'GEO data goes here' + '</h1>'


# start the web server
if __name__ == '__main__':
    app.run(debug=True)

