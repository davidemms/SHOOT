import os
import sys

import pytest

d_here = os.path.dirname(__file__) + "/"
sys.path.append(d_here + "..")
d_data = d_here + "data/"

from shoot import utils


def test_og_parts_to_ogs():
    ogs = ["0000001.3", "0000014.1", "0001000", "0000001.23", "0000002.3"]
    scores = [0.0, 0.001, 0.002, 0.003, 0.004] 
    
    ogs, scores = utils.og_parts_to_ogs(ogs, scores)
    assert ogs == ["0000001", "0000014", "0001000", "0000002"]
    assert scores == [0.0, 0.001, 0.002, 0.004]
    
    
def test_og_parts_to_ogs_all_full():
    ogs = ["0000001", "0000014", "0000002"]
    scores = [0.0, 0.001, 0.002]
    
    ogs, scores = utils.og_parts_to_ogs(ogs, scores)
    assert ogs == ["0000001", "0000014", "0000002"]
    assert scores == [0.0, 0.001, 0.002]
    
    
def test_og_parts_to_ogs_empty():
    ogs, scores = utils.og_parts_to_ogs([], [])
    
    assert len(ogs) == 0
    assert len(scores) == 0
    
    

