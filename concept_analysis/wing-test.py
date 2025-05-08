import sys
from pathlib import Path
ROOT_DIR = Path(__file__).parent.parent
sys.path.append(ROOT_DIR) # Replace with your own path so you can use .py files from other directories
sys.path.append(f'{ROOT_DIR}/concept_analysis') # Replace with your own path so you can use .py files from other directories

from wing import Wing
from helper import *

wing = Wing(kg_2_N(1100), kmph_2_mps(90), 12) # Example values for mtow and v_cruise
print(wing.wing_mass)
print(wing.mtom)
