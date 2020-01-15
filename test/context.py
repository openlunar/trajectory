import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from orbit import Orbit
from patched_conic import *

import numpy as np
from scipy.linalg import norm

