import time

class _Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s] ' % self.name,end=None)
        print('Elapsed: %ss' % (time.time() - self.tstart))
