"""result.py

A class that represents the result of an OSPREY design
"""
from osprey_test_suite.design.design import HEADERS
import json

def bounds(func):
    """bounds

    Decorator function that ensures that the first two args are in
    increasing order, and that all other args and kwargs are within
    the first two args.
    """
    def inner(*args, **kwargs):
        l = args[1]
        u = args[2]

        assert l <= u, \
                "%s is not less than %s" % (l, u)

        for value in list(args[3:]) + kwargs.keys():
            assert (value is None or (value >= l and value <= u)), \
                    "%s is not in the interval [%s, %s]" % \
                    (value, l, u)
        return func(*args, **kwargs)
    return inner

class Result():
    """Result

    Class to wrap K* design results
    """
    def __init__(self,  **kwargs):
        # Get all kwargs
        self.__dict__.update(dict.fromkeys(HEADERS))
        self.__dict__.update(kwargs)

        # Specify special, required kwargs
        self.status = kwargs.get("status", None)
        self.design_name = kwargs.get("design name", None)
        self.algorithm = kwargs.get("algorithm", None)
        self.epsilon = kwargs.get("epsilon", None)
        self.runtime = kwargs.get("runtime (s)", None)

        # Make results contain Bounds
        # define function that returns None if input is none
        def converter(x):
            if x == "none":
                return None
            else:
                return float(x)

        raw_results = kwargs.get("results", None)
        if raw_results is not None:
            self.results = dict( [ (e["sequence"],
                                    Bounds(converter(e["lowerbound"]),
                                           converter(e["upperbound"]),
                                           converter(e["kscore"])
                                          )
                                   )
                                  for e in raw_results])

    def consistent_with(self, other):
        """consistent_with

        Checks to make sure that design is the same,
        epsilon is the same, and bounds are correct.
        """
        if self != other:
            return False
        for s in self.results.keys():
            if s not in other.results.keys():
                return False
            if self.results.get(s) != other.results.get(s):
                return False
        return True


    def __eq__(self, other):
        return (self.design_name == other.design_name
                and self.epsilon == other.epsilon)

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def from_file(cls, f):
        return Result(**json.load(f))

class Bounds():
    @bounds
    def __init__(self, lower, upper, point=None):
        self.lower = lower
        self.upper = upper
        self.point = point

    def __eq__(self,other):
        """True if bounds overlap
        """
        return self.upper >= other.lower and other.upper >= self.lower

    def __ne__(self, other):
        """True if bounds do not overlap
        """
        return  self.upper < other.lower or other.upper < self.lower

    def __gt__(self, other):
        """True if self.upper > other.upper
        """
        return self.upper > other.upper

    def __lt__(self, other):
        """True if self.upper < other.upper
        """
        return self.upper < other.upper
    def __str__(self):
        return str.format("[%.3f, %.3f] (log10)" % (self.lower,
                                                    self.upper))
