#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import ftplib
import os

from datetime import datetime
from pprint import pprint


def log(text, pprint_on=False):
    if pprint_on:
        pprint(text)
        return
    print(timestamp() + text)


def timestamp():
    return "[{}] ".format(datetime.now().strftime("%H:%M"))


def ftp_download_all(target_dir, domain, user, passwd, dldir):
    log("FTP connecting to {}".format(domain))
    ftp = ftplib.FTP(domain)
    log("Login with credentials: '{}', '{}'".format(user, passwd))
    ftp.login(user=user, passwd=passwd)
    log("Change remote dir to: '{}'".format(dldir))
    ftp.cwd(dldir)
    all_files = ftp.nlst()
    log("{} files found:".format(len(all_files)))
    log(all_files, pprint_on=True)
    target_dir and os.chdir(target_dir)

    for f in all_files:
        if os.path.isfile(f):
            log("'{}': Existent file, skip download.".format(f))
            continue
        log("'{}': Downloading...".format(f))
        ftp.retrbinary("RETR " + f, open(f, "wb").write)

    print()
    log("Check {}".format(os.getcwd()))
    ftp.quit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("target_dir")
    args = parser.parse_args()

    ftp_download_all(target_dir=args.target_dir, url="ftp.1000genomes.ebi.ac.uk",
                     user="anonymous", passwd="anonymous",
                     dldir="vol1/ftp/release/20130502/")

