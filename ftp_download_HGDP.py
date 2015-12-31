#! /usr/bin/env python3
# -*- coding: utf-8 -*-

#  import argparse
import ftplib
import os

from datetime import datetime
from pprint import pprint


def log(text, pprint_on=False):
    if pprint_on and len(text) < 10:
        pprint(text)
        return

    if type(text) != str:
        text = ", ".join(text)

    print(timestamp() + text)


def timestamp():
    return "[{}] ".format(datetime.now().strftime("%H:%M"))


def ftp_download_all(target_dir, domain, user, passwd, dldir):
    original_wd = os.getcwd()
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

    log("Check '{}'".format(os.getcwd()))
    print()
    ftp.quit()
    os.chdir(original_wd)


if __name__ == "__main__":
    #  parser = argparse.ArgumentParser()
    #  parser.add_argument("target_dir")
    #  args = parser.parse_args()

    ftp_download_all(target_dir="../HGDP_data/dataset_1_HGDP-CEPH_v3/",
                     domain="ftp.cephb.fr",
                     user="anonymous", passwd="anonymous",
                     dldir="hgdp_v3/")

    ftp_download_all(target_dir="../HGDP_data/dataset_2_supp1_Stanford/",
                     domain="ftp.cephb.fr",
                     user="anonymous", passwd="anonymous",
                     dldir="hgdp_supp1/")

    ftp_download_all(target_dir="../HGDP_data/dataset_3_supp2_UMichigan/",
                     domain="ftp.cephb.fr",
                     user="anonymous", passwd="anonymous",
                     dldir="hgdp_supp2/GENO/")

    ftp_download_all(target_dir="../HGDP_data/dataset_4_supp3_MPlank/",
                     domain="ftp.cephb.fr",
                     user="anonymous", passwd="anonymous",
                     dldir="hgdp_supp3/")

