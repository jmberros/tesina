#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import ftplib
import os
from datetime import datetime


def log(text):
    print(timestamp() + text)


def timestamp():
    return "[{}] ".format(datetime.now().strftime("%H:%M"))


def download_1000genomes_data(target_dir=None):
    ftp = ftplib.FTP("ftp.1000genomes.ebi.ac.uk")
    ftp.login(user="anonymous", passwd="anonymous")
    ftp.cwd("vol1/ftp/release/20130502/")
    all_files = ftp.nlst()
    target_dir and os.chdir(target_dir)

    for f in all_files:
        print("-"*80)
        print(f)
        if os.path.isfile(f):
            log("Existent file, skip download")
            continue
        log("Downloading...")
        ftp.retrbinary("RETR " + f, open(f, "wb").write)
        log("Done")

    print()
    print("Check {}".format(os.getcwd()))
    ftp.quit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("target_dir")
    args = parser.parse_args()
    download_1000genomes_data(args.target_dir)
