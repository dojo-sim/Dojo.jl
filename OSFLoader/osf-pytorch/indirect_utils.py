"""Utility functions for indirect illumination."""
import math

import torch

import ray_utils


def create_ray_batch(pts, near, far, rays_i, use_viewdirs):
    """Create batch for indirect rays.

    Args:
        pts: [R, S, 3] float tensor. Primary points.
        near: Near sampling bound.
        far: Far sampling bound.
        rays_i: [R, 1] float tensor. Ray image IDs.
        use_viewdirs: bool. Whether to use view directions.

    Returns:
        ray_batch: [RS, M] float tensor. Batch of secondary rays containing one secondary
            ray for each primary ray sample. Each secondary ray originates at a primary
            point and points in the direction towards the randomly sampled (indirect) light
            source.
    """
    num_primary_rays = pts.size()[0]
    num_primary_samples = pts.size()[1]

    rays_dst = ray_utils.sample_random_lightdirs(num_primary_rays, num_primary_samples)  # [R, S, 3]

    rays_o = pts.view(-1, 3)  # [RS, 3]
    rays_dst = pts.view(-1, 3)  # [RS, 3]

    rays_near = torch.full((rays_o.size()[0], 1), near).float()  # [RS, 1]
    rays_far = torch.full((rays_o.size()[0], 1), far).float()  # [RS, 1]

    rays_i = torch.tile(rays_i[:, None, :], (1, num_primary_samples, 1))  # [R?, S, 1]
    rays_i = rays_i.view(-1, 1)  # [RS, 1]

    ray_batch = ray_utils.create_ray_batch(
        rays_o=rays_o, rays_dst=rays_dst, rays_near=rays_near, rays_far=rays_far,
        rays_i=rays_i, use_viewdirs=use_viewdirs)

    return ray_batch
