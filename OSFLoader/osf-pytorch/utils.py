import os
import signal
import time
from collections import deque, defaultdict
from itertools import product
import numpy as np
import argparse
import yaml
import torch
from tqdm import tqdm
import lpips

class LPIPS:
    loss_fn_alex = None

    @staticmethod
    def calculate(img_a, img_b):
        img_a, img_b = [img.permute([2, 1, 0]).unsqueeze(0) for img in [img_a, img_b]]
        if LPIPS.loss_fn_alex == None: # lazy init
            LPIPS.loss_fn_alex = lpips.LPIPS(net='alex', version='0.1')
        return LPIPS.loss_fn_alex(img_a, img_b)
