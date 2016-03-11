from datetime import datetime


def debug(msg):
    timestamp = "{:[%H:%M:%S]} ".format(datetime.now())
    print(timestamp + msg)

