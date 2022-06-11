"""Data loader for NRF (Neural Reflectance Fields) data."""
import os
import torch
import numpy as np
import imageio 
import json




trans_t = lambda t: torch.tensor([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,1,t],
    [0,0,0,1]
], dtype=torch.float)

rot_phi = lambda phi: torch.tensor([
    [1,0,0,0],
    [0,np.cos(phi),-np.sin(phi),0],
    [0,np.sin(phi), np.cos(phi),0],
    [0,0,0,1]
], dtype=torch.float)

rot_theta = lambda th: torch.tensor([
    [np.cos(th),0,-np.sin(th),0],
    [0,1,0,0],
    [np.sin(th),0, np.cos(th),0],
    [0,0,0,1]
], dtype=torch.float)


def pose_spherical(theta, phi, radius):
    c2w = trans_t(radius)
    c2w = rot_phi(phi/180.*np.pi) @ c2w
    c2w = rot_theta(theta/180.*np.pi) @ c2w
    c2w = torch.Tensor(np.array([[-1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])) @ c2w
    return c2w



def load_nrf_data(basedir, half_res=False, testskip=1, n_examples=None, no_shuffle=False, no_split=False, train_as_val=False, start=None, end=None):
    # Load camera info from numpys.
    extrinsics = np.load(os.path.join(basedir, "in_camExtrinsics.npy"))
    focal_lengths = np.load(os.path.join(basedir, "in_camFocal.npy"))

    # Determine split indices.
    n_total_examples = len(extrinsics)
    if start is None:
        start = 0  # Inclusive
    if end is None:
        end = n_total_examples  # Non-inclusive
    idxs = np.arange(start, end)  # Use the range to define the example indices to include.
    if n_examples is None:
        n_examples = len(idxs)
    if not no_shuffle:
        np.random.shuffle(idxs)
    idxs = idxs[:n_examples]  # Take first n examples after range and shuffling are applied.
    # Make train val and test all the same.
    if no_split:
        train_idxs = idxs
        val_idxs = idxs
        test_idxs = idxs
    else:
        n_test = int(n_examples * 0.15)
        n_train = n_examples - n_test
        if train_as_val:
            n_val = n_train
            train_idxs = idxs[:n_train]
            val_idxs = train_idxs
            test_idxs = idxs[n_train:]
            assert not any(np.isin(train_idxs, test_idxs))
        else:
            n_val = int(n_train * 0.15)
            n_train -= n_val
            train_idxs = idxs[:n_train]
            val_idxs = idxs[n_train:n_train + n_val]
            test_idxs = idxs[n_train + n_val:]
            assert not any(np.isin(train_idxs, val_idxs))
            assert not any(np.isin(train_idxs, test_idxs))
        assert len(train_idxs) == n_train
        assert len(val_idxs) == n_val
        assert len(test_idxs) == n_test
    split2idxs = {
        'train': train_idxs,
        'val': val_idxs,
        'test': test_idxs,
    }

    all_imgs = []
    all_poses = []
    counts = [0]
    for s, idxs in split2idxs.items():
        # Specify the skip size depending on the split.
        skip = 1 if (s=='train' or testskip==0) else testskip

        # Extract the images for the current split.
        imgs, poses = [], []
        for idx in idxs[::skip]:
            png_path = os.path.join(basedir, 'png', f'{idx:03}.png')
            imgs.append(imageio.imread(png_path))
            poses.append(extrinsics[idx])

        # Normalize images.
        imgs = (np.array(imgs) / 255.).astype(np.float32)

        poses = np.array(poses).astype(np.float32)
        counts.append(counts[-1] + imgs.shape[0])
        all_imgs.append(imgs)
        all_poses.append(poses)

    # Create a list where each element contains example indices for each split.
    i_split = [np.arange(counts[i], counts[i+1]) for i in range(3)]

    imgs = np.concatenate(all_imgs, 0)
    poses = np.concatenate(all_poses, 0)

    # Extract the height and width from the shape of the first image example.
    H, W = imgs[0].shape[:2]

    # Compute the focal length.
    focal = focal_lengths[0][0]

    # Hard-coded render poses for generating test spiral videos.
    render_poses = torch.stack([pose_spherical(angle, -30.0, 4.0) for angle in np.linspace(-180,180,40+1)[:-1]], 0)

    if half_res:
        H = H//2
        W = W//2
        focal = focal/2.

        imgs_half_res = np.zeros((imgs.shape[0], H, W, 4))
        for i, img in enumerate(imgs):
            imgs_half_res[i] = cv2.resize(img, (H, W), interpolation=cv2.INTER_AREA)
        imgs = imgs_half_res
        # imgs = tf.image.resize_area(imgs, [400, 400]).numpy()


    return imgs, poses, render_poses, [H, W, focal], i_split
