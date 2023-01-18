import os
import pathlib
import re

def _get_files(path:pathlib.Path, regex:str):
    r = re.compile(regex)
    for file in os.listdir(path):
        print(file)
        if r.match(file):
            yield file

