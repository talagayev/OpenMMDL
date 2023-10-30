import pytest
import openmm as mm
import openmm.unit as unit
from openmm.app import PDBFile, PDBxFile
from pdbfixer.pdbfixer import PDBFixer, proteinResidues, dnaResidues, rnaResidues, _guessFileFormat
from flask import Flask, request, session, g, render_template, make_response, send_file, url_for
from werkzeug.utils import secure_filename
from multiprocessing import Process, Pipe
import datetime
import os
import shutil
import signal
import sys
import tempfile
import threading
import time
import traceback
import webbrowser
import zipfile
from openmmdl.openmmdl_setup.openmmdlsetup import *

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

def test_showSelectFileType(client):
    response = client.get('/')
    assert response.status_code == 200
    assert b"selectFileType.html" in response.data

def test_selectFiles(client):
    response = client.get('/selectFiles?type=pdb')
    assert response.status_code == 200
    assert b"configurePdbFile.html" in response.data

def test_showConfigureFiles(client):
    with client.session_transaction() as sess:
        sess['fileType'] = 'pdb'
    response = client.get('/showConfigureFiles')
    assert response.status_code == 200
    assert b"configurePdbFile.html" in response.data
