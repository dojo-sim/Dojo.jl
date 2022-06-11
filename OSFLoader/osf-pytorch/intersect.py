import torch

import box_utils
import ray_utils


def apply_mask_to_tensors(mask, tensors):
    """Applies mask to a list of tensors.

    Args:
        mask: [R, ...]. Mask to apply.
        tensors: List of [R, ...]. List of tensors.

    Returns:
        intersect_tensors: List of shape [R?, ...]. Masked tensors.
    """
    intersect_tensors = []
    for t in tensors:
        #intersect_t = torch.masked_select(t, mask)  # [R?, ...]
        intersect_t = t[mask]
        intersect_tensors.append(intersect_t)
    return intersect_tensors


def get_full_intersection_tensors(ray_batch):
    """Test case that selects all rays."""
    # mask = torch.rand(ray_batch.size())[:, 0]  # [R?,]
    # mask = (mask > 0.5).bool()
    mask = torch.ones_like(ray_batch, dtype=torch.bool)[:, 0]

    n_intersect = mask.size()[0]  # R?
    indices = torch.arange(n_intersect, dtype=torch.int).unsqueeze(1)  # [R?, 1]
    indices = apply_mask_to_tensors(  # [R?, M]
        mask=mask,  # [R,]
        tensors=[indices])[0]  # [R, M]

    bounds = ray_batch[:, 6:8]  # [R?,]
    bounds = apply_mask_to_tensors(  # [R?, M]
        mask=mask,  # [R,]
        tensors=[bounds])[0]  # [R, M]
    return mask, indices, bounds


def compute_object_intersect_tensors(ray_batch, box_center, box_dims):
    """Compute rays that intersect with bounding boxes.

    Args:
        ray_batch: [R, M] float tensor. Batch of rays.
        box_center: List of 3 floats containing the (x, y, z) center of bbox.
        box_dims: List of 3 floats containing the x, y, z dimensions of the bbox.

    Returns:
        intersect_ray_batch: [R?, M] float tensor. Batch of intersecting rays.
        indices: [R?, 1] float tensor. Indices of intersecting rays.
    """
    # Check that bbox params are properly formed.
    for lst in [box_center, box_dims]:
        assert type(lst) == list
        assert len(lst) == 3
        assert all((isinstance(x, int) or isinstance(x, float)) for x in lst)

    # For now, we assume bbox has no rotation.
    num_rays = ray_batch.size()[0]  # R
    box_center = torch.tile(torch.tensor(box_center), (num_rays, 1)).float()  # [R, 3]
    box_dims = torch.tile(torch.tensor(box_dims), (num_rays, 1)).float()  # [R, 3]
    box_rotation = torch.tile(torch.eye(3).unsqueeze(0), (num_rays, 1, 1)).float()  # [R, 3, 3]

    # Compute ray-bbox intersections.
    bounds, indices, mask = box_utils.compute_ray_bbox_bounds_pairwise(  # [R', 2], [R',], [R,]
        rays_o=ray_batch[:, 0:3],  # [R, 3]
        rays_d=ray_batch[:, 3:6],  # [R, 3]
        box_length=box_dims[:, 0],  # [R,]
        box_width=box_dims[:, 1],  # [R,]
        box_height=box_dims[:, 2],  # [R,]
        box_center=box_center,  # [R, 3]
        box_rotation=box_rotation)  # [R, 3, 3]

    # Apply the intersection mask to the ray batch.
    intersect_ray_batch = apply_mask_to_tensors(  # [R?, M]
        mask=mask,  # [R,]
        tensors=[ray_batch])[0]  # [R, M]

    # Update the near and far bounds of the ray batch with the intersect bounds.
    intersect_ray_batch = ray_utils.update_ray_batch_bounds(  # [R?, M]
        ray_batch=intersect_ray_batch,  # [R?, M]
        bounds=bounds)  # [R?, 2]
    return intersect_ray_batch, indices, mask  # [R?, M], [R?, 1]
