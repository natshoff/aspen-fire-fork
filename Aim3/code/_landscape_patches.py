
"""
10-meter aspen patch metrics
author: maxwell.cook@colorado.edu
"""

import os, sys, time
import geopandas as gpd
import pylandstats as pls
import multiprocessing as mp
import concurrent.futures
import numpy as np
from multiprocessing import Pool, cpu_count
from tqdm.notebook import tqdm
from shapely.geometry import box

