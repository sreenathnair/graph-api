#!/usr/bin/env python

"""
run.py: Initialises and starts Flask application
"""

from app import app
import os

__author__          = "Sreenath Sasidharan Nair"
__version__         = "1.0"
__email__           = "sreenath@ebi.ac.uk"

app.secret_key = os.urandom(24)
port = int(os.environ.get('PORT', 5000))
app.run(host='0.0.0.0', port=port)