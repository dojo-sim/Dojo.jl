"""Various camera utility functions."""

import numpy as np


def w2c_to_c2w(w2c):
    """
    Args:
        w2c: [N, 4, 4] np.float32. World-to-camera extrinsics matrix.

    Returns:
        c2w: [N, 4, 4] np.float32. Camera-to-world extrinsics matrix.
    """
    R = w2c[:3, :3]
    T = w2c[:3, 3]

    c2w = np.eye(4, dtype=np.float32)
    c2w[:3, 3] = -1 * np.dot(R.transpose(), w2c[:3, 3])
    c2w[:3, :3] = R.transpose()
    return c2w
