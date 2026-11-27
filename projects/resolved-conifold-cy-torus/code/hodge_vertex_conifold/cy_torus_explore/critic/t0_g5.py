import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t0_special import run
run('fresh g5', '-5/4','19/7','2/3','17/5', gs=(5,))
