import os

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np


def main():
    pts_dir = "/viscam/u/mguo95/workspace/osf/logs/cartoon_s0e1/pts"
    pts_plots_dir = "/viscam/u/mguo95/workspace/osf/logs/cartoon_s0e1/pts_plots"
    
    it = 1
    pts_path = os.path.join(pts_dir, f"{it:06}.npy")

    pts = np.load(pts_path)

    x_min = np.min(pts[:, :, 0])
    x_max = np.max(pts[:, :, 0])
    y_min = np.min(pts[:, :, 1])
    y_max = np.max(pts[:, :, 1])
    z_min = np.min(pts[:, :, 2])
    z_max = np.max(pts[:, :, 2])

    end = len(pts)
    # end = 5
    # for i in list(range(end)) + [(0, end)]:
    for i in [(0, end)]:
        if type(i) == tuple:
            start, end = i
        else:
            start = i
            end = i + 1
        iter_dir = os.path.join(pts_plots_dir, f"{it:06}")
        pts_plot_path = os.path.join(iter_dir, f"{start:04}_{end:04}.png")
        os.makedirs(iter_dir, exist_ok=True)
        # if os.path.exists(pts_plot_path):
        #     continue

        xs = pts[start:end, :, 0]
        ys = pts[start:end, :, 1]
        zs = pts[start:end, :, 2]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter3D(xs, ys, zs, s=0.1)
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])
        ax.set_zlim([z_min, z_max])
        plt.savefig(pts_plot_path)
        plt.close()


if __name__ == '__main__':
    main()
