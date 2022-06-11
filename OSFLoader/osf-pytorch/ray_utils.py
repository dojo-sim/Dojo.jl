"""Utility functions for ray computation."""
import math
import numpy as np
from scipy.spatial.transform import Rotation as R
import box_utils
import torch


def apply_batched_transformations(inputs, transformations):
    """Batched transformation of inputs.

    Args:
        inputs: List of [R, S, 3]
        transformations: [R, 4, 4]

    Returns:
        transformed_inputs: List of [R, S, 3]
    """
    # if rotation_only:
    #    transformations[:, :3, 3] = torch.zeros((3,), dtype=torch.float)

    transformed_inputs = []
    for x in inputs:
        N_samples = x.size()[1]
        homog_transformations = transformations.unsqueeze(1)  # [R, 1, 4, 4]
        homog_transformations = torch.tile(homog_transformations, (1, N_samples, 1, 1))  # [R, S, 4, 4]
        homog_component = torch.ones_like(x)[..., 0:1]  # [R, S, 1]
        homog_x = torch.cat((x, homog_component), axis=-1)  # [R, S, 4]
        homog_x = homog_x.unsqueeze(2)
        transformed_x = torch.matmul(
            homog_x,
            torch.transpose(homog_transformations, 2, 3))  # [R, S, 1, 4]
        transformed_x = transformed_x[..., 0, :3]  # [R, S, 3]
        transformed_inputs.append(transformed_x)
    return transformed_inputs


def get_transformation_from_params(params):
    translation, rotation = [0, 0, 0], [0, 0, 0]
    if 'translation' in params:
        translation = params['translation']
    if 'rotation' in params:
        rotation = params['rotation']
    translation = torch.tensor(translation, dtype=torch.float)
    rotmat = torch.tensor(R.from_euler('xyz', rotation, degrees=True).as_matrix(), dtype=torch.float)
    return translation, rotmat


def rotate_dirs(dirs, rotmat):
    """
    Args:
        dirs: [R, 3] float tensor.
        rotmat: [3, 3]
    """
    if type(dirs) == np.ndarray:
        dirs = torch.tensor(dirs).float()
    rotmat = rotmat.unsqueeze(0)
    rotmat = torch.tile(rotmat, (dirs.shape[0], 1, 1))  # [R, 3, 3]
    dirs_obj = torch.matmul(dirs.unsqueeze(1), torch.transpose(rotmat, 1, 2))  # [R, 1, 3]
    dirs_obj = dirs_obj.squeeze(1)  # [R, 3]
    return dirs_obj


def transform_dirs(dirs, params, inverse=False):
    translation, rotmat = get_transformation_from_params(params)  # [3,], [3, 3]
    if inverse:
        translation *= -1
        rotmat = torch.transpose(rotmat, 0, 1)  # [3, 3]
    dirs_transformed = rotate_dirs(dirs, rotmat)
    dirs_transformed += translation[None, :]
    return dirs_transformed


def transform_rays(ray_batch, params, use_viewdirs, inverse=False):
    """Transform rays into object coordinate frame given o2w transformation params.

    Note: do not assume viewdirs is always the normalized version of rays_d (e.g., in staticcam case).

    Args:
        ray_batch: [R, M] float tensor. Batch of rays.
        params: Dictionary containing transformation parameters:
            'translation': List of 3 elements. xyz translation.
            'rotation': List of 3 euler angles in xyz.
        use_viewdirs: bool. Whether to we are using viewdirs.
        inverse: bool. Whether to apply inverse of the transformations provided in 'params'.

    Returns:
        ray_batch_obj: [R, M] float tensor. The ray batch, in object coordinate frame.
    """
    rays_o, rays_d = ray_batch[:, 0:3], ray_batch[:, 3:6]
    translation, rotmat = get_transformation_from_params(params)  # [3,], [3, 3]

    if inverse:
        translation = -1 * translation  # [3,]
        rotmat = torch.transpose(rotmat, 1, 0)  # [3, 3]

    translation_inverse = -1 * translation
    rotmat_inverse = torch.transpose(rotmat, 1, 0)

    # Transform the ray origin.
    rays_o_obj, _ = box_utils.ray_to_box_coordinate_frame_pairwise(
        box_center=translation_inverse,
        box_rotation_matrix=rotmat_inverse,
        rays_start_point=rays_o,
        rays_end_point=rays_d)

    # Only apply rotation to rays_d.
    rays_d_obj = rotate_dirs(rays_d, rotmat)

    ray_batch_obj = update_ray_batch_slice(ray_batch, rays_o_obj, 0, 3)
    ray_batch_obj = update_ray_batch_slice(ray_batch_obj, rays_d_obj, 3, 6)
    if use_viewdirs:
        # Grab viewdirs from the ray batch itself. Because it may be different from rays_d
        # (as in the staticcam case).
        viewdirs = ray_batch[:, 8:11]
        viewdirs_obj = rotate_dirs(viewdirs, rotmat)
        ray_batch_obj = update_ray_batch_slice(ray_batch_obj, viewdirs_obj, 8, 11)
    return ray_batch_obj


def transform_points_into_world_coordinate_frame(pts, params, check_numerics=False):
    translation, rotmat = get_transformation_from_params(params)  # [3,], [3, 3]

    # pts_flat = pts.view(-1, 3)  # [RS, 3]
    # num_examples = pts_flat.size()[0]  # RS

    # translation = translation.unsqueeze(0)
    # translation = torch.tile(translation, (num_examples, 1))  # [RS, 3]
    # rotmat = rotmat.unsqueeze(0)
    # rotmat = torch.tile(rotmat, (num_examples, 1, 1))

    # # pts_flat_transformed = torch.matmul(pts_flat[:, None, :], torch.transpose(rotmat, 2, 1))  # [RS, 1, 3]
    # pts_flat_transformed = pts_flat[:, None, :]  # [RS, 1, 3]
    # pts_flat_transformed += translation[:, None, :]  # [RS, 1, 3]
    # pts_transformed = pts_flat_transformed.view(pts.size())  # [R, S, 3]
    chunk = 256
    # Check batch transformations works without rotation.
    if check_numerics:
        transformations = np.eye(4)
        transformations[:3, 3] = translation
        transformations = torch.tensor(transformations, dtype=torch.float)  # [4, 4]
        transformations = torch.tile(transformations[None, ...], (pts.size()[0], 1, 1))  # [R, 4, 4]
        pts_transformed1 = []
        for i in range(0, pts.size()[0], chunk):
            pts_transformed1_chunk = apply_batched_transformations(
                inputs=[pts[i:i+chunk]], transformations=transformations[i:i+chunk])[0]
            pts_transformed1.append(pts_transformed1_chunk)
        pts_transformed1 = torch.cat(pts_transformed1, dim=0)

        pts_transformed2 = pts + translation[None, None, :]

    # Now add rotation
    transformations = np.eye(4)
    transformations = torch.tensor(transformations, dtype=torch.float)
    transformations[:3, :3] = rotmat
    transformations[:3, 3] = translation
    #transformations = torch.tensor(transformations, dtype=torch.float)  # [4, 4]
    transformations = torch.tile(transformations[None, ...], (pts.size()[0], 1, 1))  # [R, 4, 4]
    pts_transformed = []
    for i in range(0, pts.size()[0], chunk):
        pts_transformed_chunk = apply_batched_transformations(
            inputs=[pts[i:i+chunk]], transformations=transformations[i:i+chunk])[0]
        pts_transformed.append(pts_transformed_chunk)
    pts_transformed = torch.cat(pts_transformed, dim=0)
    return pts_transformed


# def transform_rays(ray_batch, translation, use_viewdirs):
#    """Apply transformation to rays.

#    Args:
#        ray_batch: [R, M] float tensor. All information necessary
#            for sampling along a ray, including: ray origin, ray direction, min
#            dist, max dist, and unit-magnitude viewing direction.
#        translation: [3,] float tensor. The (x, y, z) translation to apply.
#        use_viewdirs: Whether to use view directions.

#    Returns:
#        ray_batch: [R, M] float tensor. Transformed ray batch.
#    """
#   assert translation.size()[0] == 3, "translation.size()[0] must be 3..."

#    # Since we are only supporting translation for now, only ray origins need to be
#    # modified. Ray directions do not need to change.
#    rays_o = ray_batch[:, 0:3] + translation
#    rays_remaining = ray_batch[:, 3:]
#    ray_batch = torch.cat((rays_o, rays_remaining), dim=1)
#    return ray_batch

def compute_rays_length(rays_d):
    """Compute ray length.

    Args:
        rays_d: [R, 3] float tensor. Ray directions.

    Returns:
        rays_length: [R, 1] float tensor. Ray lengths.
    """
    rays_length = torch.norm(rays_d, dim=-1, keepdim=True)  # [N_rays, 1]
    return rays_length


def normalize_rays(rays):
    """Normalize ray directions.

    Args:
        rays: [R, 3] float tensor. Ray directions.

    Returns:
        normalized_rays: [R, 3] float tensor. Normalized ray directions.
    """
    normalized_rays = rays / compute_rays_length(rays_d=rays)
    return normalized_rays


def compute_ray_dirs_and_length(rays_o, rays_dst):
    """Compute ray directions.

    Args:
        rays_o: [R, 3] float tensor. Ray origins.
        rays_dst: [R, 3] float tensor. Ray destinations.

    Returns:
        rays_d: [R, 3] float tensor. Normalized ray directions.
    """
    # The ray directions are the difference between the ray destinations and the
    # ray origins.
    rays_d = rays_dst - rays_o  # [R, 3]  # Direction out of light source

    # Compute the length of the rays.
    rays_length = compute_rays_length(rays_d=rays_d)

    # Normalized the ray directions.
    rays_d = rays_d / rays_length  # [R, 3]  # Normalize direction
    return rays_d, rays_length


def update_ray_batch_slice(ray_batch, x, start, end):
    left = ray_batch[:, :start]  # [R, ?]
    right = ray_batch[:, end:]  # [R, ?]
    updated_ray_batch = torch.cat((left, x, right), dim=-1)
    return updated_ray_batch


def update_ray_batch_bounds(ray_batch, bounds):
    updated_ray_batch = update_ray_batch_slice(ray_batch=ray_batch, x=bounds,
                                               start=6, end=8)
    return updated_ray_batch


def create_ray_batch(
    rays_o, rays_dst, rays_i, use_viewdirs, rays_near=None, rays_far=None, epsilon=1e-10):
    # Compute the ray directions.
    rays_d = rays_dst - rays_o  # [R,3]  # Direction out of light source
    rays_length = compute_rays_length(rays_d=rays_d)  # [R, 1]
    rays_d = rays_d / rays_length  # [R, 3]  # Normalize direction
    viewdirs = rays_d  # [R, 3]

    # If bounds are not provided, set the beginning and end of ray as sampling bounds.
    if rays_near is None:
        rays_near = torch.zeros((rays_o.size()[0], 1), dtype=torch.float) + epsilon  # [R, 1]
    if rays_far is None:
        rays_far = rays_length  # [R, 1]

    ray_batch = torch.cat((rays_o, rays_d, rays_near, rays_far), dim=-1)
    if use_viewdirs:
        ray_batch = torch.cat((ray_batch, viewdirs), dim=-1)
    ray_batch = torch.cat((ray_batch, rays_i), dim=-1)
    return ray_batch


def sample_random_lightdirs(num_rays, num_samples, upper_only=False):
    """Randomly sample directions in the unit sphere.

    Args:
        num_rays: int or tensor shape dimension. Number of rays.
        num_samples: int or tensor shape dimension. Number of samples per ray.
        upper_only: bool. Whether to sample only on the upper hemisphere.

    Returns:
        lightdirs: [R, S, 3] float tensor. Random light directions sampled from the unit
            sphere for each sampled point.
    """
    if upper_only:
        min_z = 0
    else:
        min_z = -1

    phi = torch.rand(num_rays, num_samples) * (2 * math.pi)  # [R, S]
    cos_theta = torch.rand(num_rays, num_samples) * (1 - min_z) + min_z  # [R, S]
    theta = torch.acos(cos_theta)  # [R, S]

    x = torch.sin(theta) * torch.cos(phi)
    y = torch.sin(theta) * torch.sin(phi)
    z = torch.cos(theta)

    lightdirs = torch.cat((x[..., None], y[..., None], z[..., None]), dim=-1)  # [R, S, 3]
    return lightdirs


def get_light_positions(rays_i, img_light_pos):
    """Extracts light positions given scene IDs.

    Args:
        rays_i: [R, 1] float tensor. Per-ray image IDs.
        img_light_pos: [N, 3] float tensor. Per-image light positions.

    Returns:
        rays_light_pos: [R, 3] float tensor. Per-ray light positions.
    """
    img_light_pos = torch.tensor(img_light_pos)
    rays_light_pos = img_light_pos[rays_i.long()].squeeze()  # [R, 3]
    return rays_light_pos


def get_lightdirs(lightdirs_method, num_rays=None, num_samples=None, rays_i=None,
    metadata=None, ray_batch=None, use_viewdirs=False, normalize=True):
    """Compute lightdirs.

    Args:
        lightdirs_method: str. Method to use for computing lightdirs.
        num_rays: int or tensor shape dimension. Number of rays.
        num_samples: int or tensor shape dimension. Number of samples per ray.
        rays_i: [R, 1] float tensor. Ray image IDs.
        metadata: [N, 3] float tensor. Metadata about each image. Currently only light
            position is provided.
        ray_batch: [R, M] float tensor. Ray batch.
        use_viewdirs: bool. Whether to use viewdirs.
        normalize: bool. Whether to normalize lightdirs.

    Returns;
        lightdirs: [R, S, 3] float tensor. Light directions for each sample.
    """
    if lightdirs_method == 'viewdirs':
        raise NotImplementedError
        assert use_viewdirs
        lightdirs = ray_batch[:, 8:11]  # [R, 3]
        lightdirs *= 1.5
        lightdirs = torch.tile(lightdirs[:, None, :], (1, num_samples, 1))
    elif lightdirs_method == 'metadata':
        lightdirs = get_light_positions(rays_i, metadata)  # [R, 3]
        if len(lightdirs.shape) == 1:
            lightdirs = lightdirs[None, :]
        lightdirs = torch.tile(lightdirs[:, None, :], (1, num_samples, 1))  # [R, S, 3]
    elif lightdirs_method == 'random':
        lightdirs = sample_random_lightdirs(num_rays, num_samples)  # [R, S, 3]
    elif lightdirs_method == 'random_upper':
        lightdirs = sample_random_lightdirs(num_rays, num_samples, upper_only=True)  # [R, S, 3]
    else:
        raise ValueError(f'Invalid lightdirs_method: {lightdirs_method}.')
    if normalize:
        lightdirs_flat = lightdirs.view(-1, 3)  # [RS, 3]
        lightdirs_flat = normalize_rays(lightdirs_flat)  # [RS, 3]
        lightdirs = lightdirs_flat.view(lightdirs.size())  # [R, S, 3]
    return lightdirs
