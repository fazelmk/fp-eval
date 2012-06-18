#!/usr/bin/python

import fingerprint
import db
import log
import conf
import queue

import sys
import importlib


def delete(engine):
    """ Call the delete-specific method for a fingerprint engine to
        remove all traces of it from this eval test, and from
        the engine's specific storage system """
    engine_map = fingerprint.fingerprint_index.get(engine)
    engine_class = engine_map.get("instance") if engine_map else None
    if not engine_class:
        raise Exception("%s is not a valid fingerprint engine" % (engine))
    instance = engine_class()
    instance.delete_all()

def add_queue(engine):
    engine_map = fingerprint.fingerprint_index.get(engine)
    engine_table = engine_map.get("dbmodel") if engine_map else None
    cur = db.session.query(db.FPFile).filter(db.FPFile.negative == False)\
            .outerjoin(engine_table).filter(engine_table.file_id == None)
    log.info("got %d things to do stuff with" % cur.count())
    thequeue = queue.FpQueue("ingest_%s" % engine)
    for f in cur:
        d = {"id": f.id}
        thequeue.put(d)

def ingest_callback(data):
    engine_map = fingerprint.fingerprint_index.get(engine)
    engine_class = engine_map.get("instance") if engine_map else None
    if not engine_class:
        raise Exception("%s is not a valid fingerprint engine" % (engine))
    engine_table = engine_map.get("dbmodel") if engine_map else None

    instance = engine_class()

    log.info("Importing files for engine %s" % (engine))
    fp_list = []
    cur = db.session.query(db.FPFile).filter(db.FPFile.negative == False)\
            .outerjoin(engine_table).filter(engine_table.file_id == None)
    log.info("got %d things to do stuff with" % cur.count())
    count = 0
    total = cur.count()
    for f in cur:
        (trackid, fpdata) = instance.fingerprint(f.path)
        error = "error" in fpdata
        if not error:
            e = engine_table(f, trackid)
            db.session.add(e)
            fp_list.append(fpdata)
            count += 1
        else:
            log.debug("Error parsing file %s. Error was: %s" % (f, fpdata["error"]))

        # Ingest every 100 songs
        if len(fp_list) > 99:
            log.info("Ingesting 100 files at once")
            log.info("%d/%d done" % (count, total))
            db.session.commit()
            instance.ingest_many(fp_list)
            fp_list = []

def main(engine):
    thequeue = queue.FpQueue("ingest_%s" % engine)
    thequeue.consume(ingest_callback)

#XXX: Don't want to continually fingerprint files that had errors
def xmain(engine):
    engine_map = fingerprint.fingerprint_index.get(engine)
    engine_class = engine_map.get("instance") if engine_map else None
    if not engine_class:
        raise Exception("%s is not a valid fingerprint engine" % (engine))
    engine_table = engine_map.get("dbmodel") if engine_map else None

    instance = engine_class()

    log.info("Importing files for engine %s" % (engine))
    fp_list = []
    cur = db.session.query(db.FPFile).filter(db.FPFile.negative == False)\
            .outerjoin(engine_table).filter(engine_table.file_id == None)
    log.info("got %d things to do stuff with" % cur.count())
    count = 0
    total = cur.count()
    for f in cur:
        (trackid, fpdata) = instance.fingerprint(f.path)
        error = "error" in fpdata
        if not error:
            e = engine_table(f, trackid)
            db.session.add(e)
            fp_list.append(fpdata)
            count += 1
        else:
            log.debug("Error parsing file %s. Error was: %s" % (f, fpdata["error"]))

        # Ingest every 100 songs
        if len(fp_list) > 99:
            log.info("Ingesting 100 files at once")
            log.info("%d/%d done" % (count, total))
            db.session.commit()
            instance.ingest_many(fp_list)
            fp_list = []

def show_fp_engines():
    print "Available fingerprinting engines:"
    print ", ".join(fingerprint.fingerprint_index.keys())

if __name__ == "__main__":
    import argparse
    conf.import_fp_modules()

    p = argparse.ArgumentParser()
    g = p.add_argument_group()

    g.add_argument("-f", help="Fingerprint algorithm to use")
    g.add_argument("--show-engines", action="store_true", help="Show all fingerprint engines available")
    g.add_argument("-d", action="store_true", help="Delete the database for the given algorithm")
    g.add_argument("-q", action="store_true", help="add stuff to queue")
    g.add_argument("-r", action="store_true", help="run, consuming from queue")

    args = p.parse_args()

    if args.show_engines:
        show_fp_engines()
        sys.exit(1)

    if not args.f:
        p.print_help()
        sys.exit(1)

    if args.d:
        delete(args.f)
    elif args.q:
        add_queue(args.f)
    elif args.r:
        main(args.f)
