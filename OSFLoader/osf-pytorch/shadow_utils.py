"""Utility functions for shadows."""
import torch

import ray_utils
import run_osf_helpers


def create_ray_batch(ray_batch, pts, metadata, use_viewdirs, lightdirs_method):
    """Create batch for shadow rays.

    Args:
        ray_batch: [R?, M] float tensor. Batch of primary rays
        pts: [R?, S, 3] float tensor. Primary points.
        metadata: [N, 3] float tensor. Metadata about each image. Currently only light
            position is provided.
        use_viewdirs: bool. Whether to use view directions.

    Returns:
        shadow_ray_batch: [R?S, M] float tensor. Batch of shadow rays containing one shadow
            ray for each primary ray sample. Each shadow ray originates at a primary point
            and points in the direction towards the light source.
    """
    num_primary_rays = pts.size()[0]  # R
    num_primary_samples = pts.size()[1]  # S

    # Samples are shadow ray origins.
    rays_o = pts.view(-1, 3)  # [R?S, 3]

    num_shadow_rays = rays_o.size()[0]  # R?S
    num_samples = pts.size()[1]  # S

    rays_i = ray_batch[:, 11:12]  # [R, 1]

    # Get light positions for each ray as the ray destinations.
    rays_dst = ray_utils.get_lightdirs(  # [R?, S, 3]
        lightdirs_method=lightdirs_method, num_rays=num_primary_rays,
        num_samples=num_primary_samples, rays_i=rays_i, metadata=metadata,
        ray_batch=ray_batch, use_viewdirs=use_viewdirs)
    rays_dst = rays_dst.view(rays_o.size())  # [R?S, 3]

    rays_i = torch.tile(rays_i.unsqueeze(1), (1, num_primary_samples, 1))  # [R?, S, 1]
    rays_i = rays_i.view(-1, 1)  # [R?S, 1]

    shadow_ray_batch = ray_utils.create_ray_batch(rays_o, rays_dst, rays_i, use_viewdirs)
    return shadow_ray_batch


def compute_transmittance(alpha):
    """Applies shadows to outputs.
    Args:
        alpha: [R?S, S, 1] tf.float32. Alpha predictions from the model.
    Returns:
        last_trans: [R?S,] tf.float32. Shadow transmittance per ray.
    """
    trans = run_osf_helpers.compute_transmittance(alpha=alpha[..., 0])  # [R?S, S]

    # Transmittance is computed in the direction origin -> end, so we grab the
    # last transmittance.
    last_trans = trans[:, -1]  # [R?S,]
    return last_trans
