import time


class Duration:
    def __init__(self):
        self.store = {}

    def start(self, string):
        self.store[string] = time.time()
        self.last = string

    def step(self, string=None):
        if not string:
            string = self.last
        print("%s: %ss" % (string, time.time() - self.store[string]))
        self.store[string] = time.time()


duration = Duration()
