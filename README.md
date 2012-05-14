Fingerprint evaluation
----------------------

Requirements

* SQLAlchemy
* Sqlite
* To perform query alterations:
** mpg123
** ffmpeg

Optional requirements: installations of fingerprint algorithms

* Echoprint
* Landmark fungerprinting
* Ke vision FP
* Chromaprint

Echoprint notes:
----------------

To start the servers:

    cd echoprint-server/solr/solr
    java -jar solr.jar
    ttserver -host localhost -port 1978 echoprint-tt.tch

Acoustid / Chromaprint
----------------------
Install pymad or ffmpeg

Set up acoustid server
mbslave
acoustid-index: https://github.com/lalinsky/acoustid-index
      depends: libqt4-dev, libicu-dev

To see the acoustid server:
https://132.206.14.136
api base url: https://132.206.14.136/ws/v2/

Get API Key, login with google @gmail

You have an API key, and an application key
API: 0bB2HTpl
Application: PKlUB2YR

A fingerprint ingest goes into the `submission' table

run the admin/cron/import.sh script to ingest into the main database.

License
-------
BSD 2-clause

To run the stuff:

Configuration:

    cp fingerprint.conf{.dist,}

Edit `fingerprint.conf` to

Import files into the database (will keep back a percentage of files for
false negative tests - default 10%):

    python createdb.py

Import each file with all of the fingerprinters with each FP engine.

    python ingest.py

You can choose just a single fingerprint type:

    python ingest.py -f echoprint

The import script will incrementally import all files
for the engine. To delete the already imported files, use

    python ingest.py -d [-f echoprint]

Note that this does not remove any databases that the fingerprint
system might use. You need to do that yourself. Builtin fp's
do this with
    python echoprint.py -d
    python landmark.py -d
    python ke.py -d

Development:

To create a new module:
create a python file e.g. `myfp.py`. You may want to create a myfp\_support
directory to put additional dependencies in

Make a class inside myfp.py that inherits from the Fingerprint class in
fingerprint.py
set the \_\_provides member to a short name for your fingerprint
engine. This is what you refer to it with in the configuration file.

Ensure you override the fingerprint(file), ingest(data), and lookup(file)
methods.
Create a class that inherits from db.Base in order to create a reference
database to map from file->fingerprint code.
fill in the fingerprint.fp\_index dict with information on your model.

See the example file: examplefp.py to see what to do.



