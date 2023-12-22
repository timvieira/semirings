import numpy as np
_float = float
from semirings.base import Semiring, Chart, Wrapped, make_semiring
from semirings.mert import ConvexHull, Point
from semirings.lazysort import LazySort, flatten
from semirings.maxplus import MaxPlus
from semirings.maxtimes import MaxTimes
from semirings.mintimes import MinTimes
from semirings.minplus import MinPlus
from semirings.logval import LogVal
from semirings.float import Float
from semirings.boolean import Boolean
from semirings.lazysort import LazySort
from semirings.cutsets import CutSets
from semirings.regex import Symbol
from semirings.interval import Interval, make_interval
from semirings.misc import *


Dual = dual(Float)
