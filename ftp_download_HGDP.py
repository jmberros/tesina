#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import ftputil
import os

from datetime import datetime
from pprint import pprint
from humanize import naturalsize
from collections import OrderedDict


def log(text, pprint_on=False):
    if type(text) != str:
        if pprint_on and len(text) < 10:
            print()
            pprint(text)
            print()
            return

    print("[{}] {}".format(timestamp(), str(text)))


def timestamp():
    return datetime.now().strftime("%H:%M")

def ftp_walk_and_download(domain, remotedir, user, passwd, basedir, localdir=None):
    original_dir = os.getcwd()
    log("Change to download directory: '{}'".format(basedir))
    os.chdir(basedir)

    if os.path.isdir(localdir):
        os.rename(localdir, remotedir)  # I'll rename back after the process

    log("Login to '{}' as: '{}', '{}'".format(domain, user, passwd))
    with ftputil.FTPHost(domain, user, passwd) as host:
        for root, dirs, files in host.walk(remotedir):
            os.makedirs(root, exist_ok=True)

            log("Walking down remote dir '{}'".format(root))
            log("{} files found remotely.".format(len(files)))
            for i, fpath in enumerate([os.path.join(root, fn) for fn in files]):
                file_exists_locally = os.path.isfile(fpath)

                filetag = "({}/{}) '{}'".format(i+1, len(files), fpath)
                if file_exists_locally:
                    if os.path.getsize(fpath) == host.path.getsize(fpath):
                        log("{}: File exists, skip download.".format(filetag))
                        continue

                remote_filesize = naturalsize(host.path.getsize(fpath))
                log("{}: Downloading {}".format(filetag, remote_filesize))
                host.download(fpath, fpath)

            if len(dirs) > 0:
                log("Continue with dirs: {}".format(str(dirs)))
                for dirname in dirs:
                    os.makedirs(dirname, exist_ok=True)

            print()

    log("Done. Renaming '{}' to '{}'.".format(remotedir, localdir))
    os.rename(remotedir, localdir)
    os.chdir(original_dir)
    print()


if __name__ == "__main__":
    dataset_dirnames = OrderedDict([
        ("hgdp_v3", "dataset_1_HGDP-CEPH_v3"),
        ("hgdp_supp1", "dataset_2_supp1_Stanford"),
        ("hgdp_supp2", "dataset_3_supp2_UMichigan"),
        ("hgdp_supp3", "dataset_4_supp3_MPlank"),
        ("hgdp_supp15", "dataset_15_supp15_UCLA"),
        ("hgdp_supp10", "dataset_11_supp10_Harvard"),
    ])

    BASEDIR = "/home/juan/tesina/HGDP_data"

    for remotedir, localdir in dataset_dirnames.items():
        ftp_walk_and_download(domain="ftp.cephb.fr", basedir=BASEDIR,
                              remotedir=remotedir, localdir=localdir,
                              user="anonymous", passwd="anonymous")
