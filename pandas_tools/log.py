import logging
import time


class Log:

    def __init__(self, file_name=""):
        if file_name != "":
            logging.basicConfig(filename=file_name,
                                format='%(asctime)s %(levelname)s:%(message)s',
                                level=logging.DEBUG)
            self.to_terminal = False
        else:
            self.to_terminal = True

    def info(self, statement: str):
        if self.to_terminal:
            print("{time} ---{stat}---".format(time=time.asctime(), stat=statement))
        else:
            logging.info(statement)

    def warn(self, statement: str):
        if self.to_terminal:
            print("{time} WARN---{stat}---".format(time=time.asctime(), stat=statement))
        else:
            logging.warning(statement)

    def error(self, statement: str):
        if self.to_terminal:
            print("{time} ERROR---{stat}---".format(time=time.asctime(), stat=statement))
        else:
            logging.error(statement)
