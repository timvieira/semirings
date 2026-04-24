import numpy as np
_float = float
from semirings.base import Semiring, Chart, Wrapped, make_semiring, \
    check_axioms_samples, check_axioms, check_metric_axioms

from semirings.mert import ConvexHull, Point
from semirings.lazysort import LazySort, flatten, make_lazysort_semiring
from semirings.maxplus import MaxPlus
from semirings.maxtimes import MaxTimes
from semirings.mintimes import MinTimes
from semirings.minplus import MinPlus
from semirings.logval import LogVal, LogValVector
from semirings.expectation import Expectation, SecondOrderExpectation, make_expectation
from semirings.entropy import Entropy
from semirings.float import Float
from semirings.boolean import Boolean
from semirings.count import Count
from semirings.bag import Bag
from semirings.util import derivation, post_process
from semirings.vector import make_vector
from semirings.free import FreeExpr, Sum, Prod, Star, backprop, toposort, weight, maxtimes, lazysort as free_lazysort, sample as free_sample
from semirings.cutsets import CutSets
from semirings.regex import Symbol
from semirings.interval import Interval, make_interval
from semirings.kleene import MatrixSemiring
from semirings.graph import WeightedGraph, scc_decomposition
from semirings.misc import *


Dual = dual(Float)
