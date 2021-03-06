#!/usr/bin/python

import munge
import fingerprint
import conf
import log
import queue

import db
import sqlalchemy
from sqlalchemy.orm import relationship, backref
import sqlalchemy.dialects.mysql
import random
import operator

import sys
import argparse
import datetime
import os

conf.import_fp_modules()

class Testset(db.Base):
    """ A testset is an identifier for a collection of files that are used in an evaluation """
    __tablename__ = "testset"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String(50))
    created = sqlalchemy.Column(sqlalchemy.DateTime)

    # (auto) testfiles <- list of all Testfiles

    def __init__(self, name):
        self.name = name
        now = datetime.datetime.now()
        now = now.replace(microsecond=0)
        self.created = now

    def __repr__(self):
        return "<Testset(id=%d, name=%s, created=%s)>" % (self.id, self.name, self.created)

class Testfile(db.Base):
    """ A testfile links together a testset and a file """
    __tablename__ = "testfile"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    testset_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey('testset.id'), nullable=False)
    file_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey('file.id'), nullable=False)

    file = relationship(db.FPFile)
    testset = relationship(Testset, backref="testfiles")

    def __init__(self, testset, file):
        if isinstance(testset, Testset):
            testset = testset.id
        self.testset_id = testset
        if isinstance(file, db.FPFile):
            file = file.id
        self.file_id = file

    def __repr__(self):
        return "<Testfile(id=%d, testset=%d, file=%d)>" % (self.id, self.testset_id, self.file_id)

class Run(db.Base):
    """ A run links together a testset, a munger (or set of), and a fingerprint algorithm. """
    __tablename__ = "run"

    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    # the testset listing all the files that should be tested in this run
    testset_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey('testset.id'), nullable=False)
    # A list of munge classes
    munge = sqlalchemy.Column(sqlalchemy.String(100))
    # The fp algorithm to use
    engine = sqlalchemy.Column(sqlalchemy.String(20))
    # The date this run was created
    created = sqlalchemy.Column(sqlalchemy.DateTime)
    #started
    started = sqlalchemy.Column(sqlalchemy.DateTime)
    # finished (there's a result for each testfile)
    finished = sqlalchemy.Column(sqlalchemy.DateTime)

    testset = relationship(Testset)
    # (auto) results <- list of results

    def __init__(self, testset, munge, engine):
        if not munge:
            raise Exception("munge not set")
        if not engine:
            raise Exception("engine not set")
        if isinstance(testset, Testset):
            testset = testset.id
        munge = ",".join(map(operator.methodcaller("strip"), munge.split(",")))
        self.testset_id = testset
        self.munge = munge
        self.engine = engine
        now = datetime.datetime.now()
        now = now.replace(microsecond=0)
        self.created = now
        self.started = None
        self.finished = None

    def __repr__(self):
        return "<Run(id=%d, testset=%d, engine=%s, munge=%s, created=%s, started=%s, finished=%s)>" % \
            (self.id, self.testset_id, self.engine, self.munge, self.created, self.started, self.finished)

class Result(db.Base):
    """ A row for every testset file for a run """
    __tablename__ = "result"
    __table_args__ = {'mysql_charset': 'utf8'}

    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    run_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey('run.id'), nullable=False)
    testfile_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey('testfile.id'), nullable=False)
    result = sqlalchemy.Column(sqlalchemy.dialects.mysql.VARCHAR(length=255, charset='utf8'))
    fptime = sqlalchemy.Column(sqlalchemy.Integer)
    lookuptime = sqlalchemy.Column(sqlalchemy.Integer)

    run = relationship(Run, backref="results")
    testfile = relationship(Testfile)

    def __init__(self, run, testfile, result, fptime, lookuptime):
        """
        run -- the run id
        testfile -- a db file
        result -- what the FP returns
        fptime -- the time (in ms) taken to perform a fingerprint
        lookuptime -- the time (in ms) taken to perform a lookup
        """
        if isinstance(run, Run):
            if run.id is None:
                raise Exception("This run isn't committed: %s" % str(run))
            run = run.id
        if isinstance(testfile, db.FPFile):
            testfile = testfile.id
        self.testfile_id = testfile
        self.run_id = run
        self.result = result
        self.fptime = fptime
        self.lookuptime = lookuptime

    def __repr__(self):
        return "<Result(run=%d, testfile=%d, result=%s)>" % (self.run_id, self.testfile_id, self.result)

def create_tables():
    db.create_tables()

def create_testset(name, size, holdback):
    """
    Randomly select `size' tracks from the file list and make a testset.
    If `holdback' is set, make size-holdback positives and
    holdback negatives. Otherwise randomly select size files.
    """
    if holdback is None:
        files = db.session.query(db.FPFile).all()
        random.shuffle(files)
        todo = files[:size]
    else:
        num = size - holdback
        files = db.session.query(db.FPFile).filter(db.FPFile.negative == False).all()
        random.shuffle(files)
        todo = files[:num]
        neg = db.session.query(db.FPFile).filter(db.FPFile.negative == True).all()
        random.shuffle(neg)
        todo.extend(neg[:holdback])

    testset = Testset(name)
    db.session.add(testset)
    db.session.flush()
    for fpfile in todo:
        tfile = Testfile(testset, fpfile)
        db.session.add(tfile)
    db.session.commit()

def list_testsets():
    """
    List all the testsets that have been created, including
    their size and composition of positive/negative files
    """
    tsets = db.session.query(Testset).all()
    for t in tsets:
        handle = db.session.query(Testfile).filter(Testfile.testset == t)
        count = handle.count()
        poshandle = db.session.query(Testfile).filter(Testfile.testset == t).join(Testfile.file).filter(db.FPFile.negative == False)
        poscount = poshandle.count()
        neghandle = db.session.query(Testfile).filter(Testfile.testset == t).join(Testfile.file).filter(db.FPFile.negative == True)
        negcount = neghandle.count()
        print "%s (%d): %d files - %d positive, %d negative" % (t.name, t.id, count, poscount, negcount)

def create_run(testset, fp, mungename):
    """
    Make a run. A run is a testset with a specific fingerprinter.
    Files are put through a list of munges before looked up in
    the fingerprinter
    """
    if fp not in fingerprint.fingerprint_index.keys():
        raise Exception("Unknown fingerprint name %s" % (fp))
    munges = munge.munge_classes.keys()
    for m in mungename.split(","):
        m = m.strip()
        if m not in munges:
            raise Exception("Unknown munge %s" % (m))
    testset = int(testset)
    ts = db.session.query(Testset).get(testset)
    if ts is None:
        raise Exception("Testset %d not in database" % testset)
    run = Run(testset, mungename, fp)
    db.session.add(run)
    db.session.commit()

    # Once the run has been created, make a queue that contains
    # all the testfiles to be evaluated
    thequeue = queue.FpQueue("run_%s" % run.id)
    for tf in ts.testfiles:
        data = {"testfile_id": tf.id}
        thequeue.put(data)

def reset_run(run):
    """Reset a run.
    Clear the start/end date and all results, then clear and repopulate
    the queue.
    """
    db.session.query(Result).filter(Result.run_id==run).delete()
    r = db.session.query(Run).get(run)
    r.started = None
    r.finished = None
    db.session.add(r)
    db.session.commit()
    thequeue = queue.FpQueue("run_%s" % run)
    thequeue.clear_queue()
    thequeue = queue.FpQueue("run_%s" % run)
    for tf in r.testset.testfiles:
        data = {"testfile_id": tf.id}
        thequeue.put(data)

def list_runs():
    """
    List all the runs in the database
    """
    runs = db.session.query(Run).all()
    for r in runs:
        print r

def munge_file(file, munges):
    """ apply a list of munges to `file' and return a path
    to a file to read in
    """
    if not isinstance(munges, list):
        munges = map(operator.methodcaller("strip"), munges.split(","))

    if file is None:
        log.warning("Failed to munge file (was given None!)")
        return None

    # The first thing we do is make a copy of the input file.
    # This means we can safely delete any file that is used
    # as input or output.
    nomunge = munge.NoMunge()
    tmpfile = nomunge.perform(file)
    log.debug("input file: %s" % file.encode("utf-8"))
    log.debug("tmp copy: %s" % tmpfile)

    for m in munges:
        log.debug("performing %s" % (m))
        cls = munge.munge_classes[m]
        inst = cls()
        newfile = inst.perform(tmpfile)
        remove_file(tmpfile)
        tmpfile = newfile

    return newfile

def remove_file(file):
    if file:
        os.unlink(file)

def execute_run(run_id):
    """
    Execute a run.

    This does most of the magic. 
    """
    run = db.session.query(Run).filter(Run.id == run_id)
    if run.count() == 0:
        raise Exception("No run with this id")
    else:
        run = run.one()

    if run.started is None:
        now = datetime.datetime.now()
        now = now.replace(microsecond=0)
        run.started = now
    db.session.add(run)
    db.session.commit()

    engine = run.engine
    munges = run.munge
    fpclass = fingerprint.fingerprint_index[engine]
    fp = fpclass["instance"]()
    thequeue = queue.FpQueue("run_%s" % run.id)
    ack_handles = []
    log.info("Reading queue for run %s. Got %s files" % (run.id, thequeue.size()))
    num_lookups = fp.num_lookups()
    to_lookup = []
    count = 0

    def do_fp(to_lookup):
        try:
            res = fp.lookup(to_lookup)
            if res is None:
                return False
            for r in res:
                fptime = r["fptime"]
                lookuptime = r["lookuptime"]
                fpresult = r["result"]
                t = r["track"]
                newpath = r["file"]

                remove_file(newpath)
                result = Result(run, t.id, fpresult, int(fptime), int(lookuptime))
                db.session.add(result)
            return True
        except Exception as e:
            log.warning("Error performing fingerprint")
            log.warning(e)
            raise

    while True:
        data, handle = thequeue.get()
        if data is None:
            break
        ack_handles.append(handle)
        # Find the FpFile that this Testfile points to
        t = db.session.query(Testfile).filter(Testfile.id == data["testfile_id"]).one()
        fpfile = t.file
        metadata = fp.pre_lookup(fpfile.path)

        newpath = munge_file(fpfile.path, munges)

        if newpath is None or not os.path.exists(newpath):
            log.warning("File %s doesn't exist, not fingerprinting it" % newpath)
            continue

        to_lookup.append({"track": t, "file": newpath, "data": metadata})

        done_fp = False
        if len(to_lookup) >= num_lookups:
            done_fp = do_fp(to_lookup)
            log.debug("lookup result is: %s (bool)" % done_fp)
            to_lookup = []

        count += 1
        if (num_lookups > 1 and len(to_lookup) == 0) or (num_lookups <= 1 and count % 10 == 0):
            log.info("%s more files to evaluate" % thequeue.size())
            db.session.commit()
            if done_fp:
                for h in ack_handles:
                    thequeue.ack(h)
            ack_handles = []

    done_fp = False
    if len(to_lookup) > 0:
        # Last fingerprint
        done_fp = do_fp(to_lookup)
        print "last - done_fp is %s" % done_fp
    else:
        done_fp = True

    # Mark the run as done
    now = datetime.datetime.now().replace(microsecond=0)
    run.finished = now
    # Finish any acks that are required
    if done_fp:
        for h in ack_handles:
            thequeue.ack(h)
    db.session.add(run)
    db.session.commit()

def show_munge():
    print ", ".join(sorted(munge.munge_classes.keys()))

def show_fp_engines():
    print "Available fingerprinting engines:"
    print ", ".join(fingerprint.fingerprint_index.keys())

if __name__ == "__main__":
    """
    commands:
    testset -c -n "foo" -s S -b B - create a testset of size S
                    if B is set, then B of the S are not in the database,
                    otherwise, the amount of 'new queries' is random
    testset -l - list all testsets

    makerun -t 1 -f echoprint -m "start10,noise30"
    make a run using testset 1, echoprint, and 2 munges

    run -l -- list all available runs (also the status of them? - or -s 1)
    run 1 -- execute run 1
    """
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='subparser')
    test = sub.add_parser("testset", help="create and show test sets")
    mex = test.add_mutually_exclusive_group(required=True)
    mex.add_argument("-c", action="store_true", help="create a testset")
    test.add_argument("-n", help="name of the testset")
    test.add_argument("-s", type=int, help="size of the testset")
    test.add_argument("-b", type=int, help="if set, guarantee that at least B files in the test set aren't in the FP database")
    mex.add_argument("-l", action="store_true", help="list all testsets")

    mkrun = sub.add_parser("makerun", help="make a run")
    mkrun.add_argument("-t", type=int, help="The testset to use for this run")
    mkrun.add_argument("-f", help="Fingerprint algorithm to use")
    mkrun.add_argument("-m", help="comma separated list of munges to run (in order)")
    mkrun.add_argument("-r", help="Reset a run number")
    mkrun.add_argument("--show-engines", action="store_true", help="Show all fingerprint engines available")
    mkrun.add_argument("--show-munge", action="store_true", help="List munge methods")

    run = sub.add_parser("run", help="run an evaluation")
    run = run.add_mutually_exclusive_group(required=True)
    run.add_argument("-l", action="store_true", help="list all available runs")
    run.add_argument("r", type=int, help="execute run R", nargs="?")

    args = p.parse_args()

    if args.subparser == "testset":
        if args.c:
            name = args.n
            size = args.s
            if not size or not name:
                print "size (-s) and name (-n) required"
                sys.exit(1)
            holdback = args.b
            create_testset(name, size, holdback)
        elif args.l:
            list_testsets()
    elif args.subparser == "makerun":
        if args.show_munge:
            show_munge()
        elif args.show_engines:
            show_fp_engines()
        elif args.r:
            reset_run(int(args.r))
        else:
            testset = args.t
            fp = args.f
            themunge = args.m
            if not testset or not fp or not themunge:
                print "testset (-t), fp engine (-f) and munge (-m) needed"
                sys.exit(1)
            create_run(testset, fp, themunge)
    elif args.subparser == "run":
        if args.l:
            list_runs()
        elif args.r:
            run = args.r
            execute_run(run)
