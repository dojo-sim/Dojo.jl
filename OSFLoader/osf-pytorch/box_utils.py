"""Utility functions for bounding box computation."""
import torch
import torch.nn as nn
import torch.nn.functional as F
import ray_utils

def ray_to_box_coordinate_frame_pairwise(box_center, box_rotation_matrix,
                                         rays_start_point, rays_end_point):
    """Moves a set of rays into a box's coordinate frame.

    Args:
        box_center: A tensor of size [3] or [r, 3].
        box_rotation_matrix: A tensor of size [3, 3] or [r, 3, 3].
        rays_start_point: A tensor of size [r, 3] where r is the number of rays.
        rays_end_points: A tensor of size [r, 3] where r is the number of rays.

    Returns:
        rays_start_point_in_box_frame: A tensor of size [r, 3].
        rays_end_point_in_box_frame: A tensor if size [r, 3].
    """
    r = rays_start_point.size()[0]
    box_center = torch.broadcast_to(box_center, (r, 3))
    box_rotation_matrix = torch.broadcast_to(box_rotation_matrix, (r, 3, 3))
    rays_start_point_in_box_frame = torch.matmul(
        (rays_start_point - box_center).unsqueeze(1),
        box_rotation_matrix)
    rays_end_point_in_box_frame = torch.matmul(
        (rays_end_point - box_center).unsqueeze(1),
        box_rotation_matrix)
    return (rays_start_point_in_box_frame.view(-1, 3), 
            rays_end_point_in_box_frame.view(-1, 3))


def ray_box_intersection_pairwise(box_center,
                                  box_rotation_matrix,
                                  box_length,
                                  box_width,
                                  box_height,
                                  rays_start_point,
                                  rays_end_point,
                                  exclude_negative_t=False,
                                  exclude_enlarged_t=True,
                                  epsilon=0.000001):
    """Intersects a set of rays with a box.

    Note: The intersection points are returned in the box coordinate frame.
    Note: Make sure the start and end point of the rays are not the same.
    Note: Even though a start and end point is passed for each ray, rays are
        never ending and can intersect a box beyond their start / end points.

    Args:
        box_center: A tensor of size [3] or [r, 3].
        box_rotation_matrix: A tensor of size [3, 3] or [r, 3, 3].
        box_length: A scalar tensor or of size [r].
        box_width: A scalar tensor or of size [r].
        box_height: A scalar tensor or of size [r].
        rays_start_point: A tensor of size [r, 3] where r is the number of rays.
        rays_end_point: A tensor of size [r, 3] there r is the number of rays.
        exclude_negative_t: bool.
        exclude_enlarged_t: bool.
        epsilon: A very small number.

    Returns:
        intersection_points_in_box_frame: A tensor of size [r', 2, 3]
            that contains intersection points in box coordinate frame.
        indices_of_intersecting_rays: A tensor of size [r'].
        intersection_ts: A tensor of size [r'].
    """
    r = rays_start_point.size()[0]
    box_length = box_length.expand(r)
    box_width = box_width.expand(r)
    box_height = box_height.expand(r)
    box_center = torch.broadcast_to(box_center, (r, 3))
    box_rotation_matrix = torch.broadcast_to(box_rotation_matrix, (r, 3, 3))
    rays_start_point_in_box_frame, rays_end_point_in_box_frame = (
        ray_to_box_coordinate_frame_pairwise(
            box_center=box_center,
            box_rotation_matrix=box_rotation_matrix,
            rays_start_point=rays_start_point,
            rays_end_point=rays_end_point))
    rays_a = rays_end_point_in_box_frame - rays_start_point_in_box_frame
    intersection_masks = []
    intersection_points = []
    intersection_ts = []
    box_size = [box_length, box_width, box_height]
    for axis in range(3):
        plane_value = box_size[axis] / 2.0
        for _ in range(2):
            plane_value = -plane_value
            # Compute the scalar multiples of 'rays_a' to apply in order to intersect
            # with the plane.
            t = ((plane_value - rays_start_point_in_box_frame[:, axis]) /  # [R,]
                 rays_a[:, axis])
            # The current axis only intersects with plane if the ray is not parallel
            # with the plane. Note that this will result in 't' being +/- infinity, becasue
            # the ray component in the axis is zero, resulting in rays_a[:, axis] = 0.
            intersects_with_plane = torch.abs(rays_a[:, axis]) > epsilon
            if exclude_negative_t:  # Only allow at most one negative t
                t = torch.maximum(t, torch.tensor(0.0))  # [R,]
            if exclude_enlarged_t:
                t = torch.minimum(t, torch.tensor(1.0))  # [R,]
            intersection_ts.append(t)  # [R, 1]
            intersection_points_i = []

            # Initialize a mask which represents whether each ray intersects with the
            # current plane.
            intersection_masks_i = torch.ones_like(t, dtype=torch.int32).bool()  # [R,]
            for axis2 in range(3):
                # Compute the point of intersection for the current axis.
                intersection_points_i_axis2 = (  # [R,]
                    rays_start_point_in_box_frame[:, axis2] + t * rays_a[:, axis2])
                intersection_points_i.append(intersection_points_i_axis2)  # 3x [R,]

                # Update the intersection mask depending on whether the intersection
                # point is within bounds.
                intersection_masks_i = torch.logical_and(  # [R,]
                    torch.logical_and(intersection_masks_i, intersects_with_plane),
                    torch.logical_and(
                        intersection_points_i_axis2 <= (box_size[axis2] / 2.0 + epsilon),
                        intersection_points_i_axis2 >= (-box_size[axis2] / 2.0 - epsilon)))
            intersection_points_i = torch.stack(intersection_points_i, dim=1)  # [R, 3]
            intersection_masks.append(intersection_masks_i)  # List of [R,]
            intersection_points.append(intersection_points_i)  # List of [R, 3]
    intersection_ts = torch.stack(intersection_ts, dim=1)  # [R, 6]
    intersection_masks = torch.stack(intersection_masks, dim=1)  # [R, 6]
    intersection_points = torch.stack(intersection_points, dim=1)  # [R, 6, 3]

    # Compute a mask over rays with exactly two plane intersections out of the six
    # planes. More intersections are possible if the ray coincides with a box
    # edge or corner, but we'll ignore these cases for now.
    counts = torch.sum(intersection_masks.int(), dim=1)  # [R,]
    intersection_masks_any = torch.eq(counts, 2)  # [R,]
    indices = torch.arange(intersection_masks_any.size()[0]).int()  # [R,]
    # Apply the intersection masks over tensors.
    indices = indices[intersection_masks_any]  # [R',]
    intersection_masks = intersection_masks[intersection_masks_any]  # [R', 6]
    intersection_points = intersection_points[intersection_masks_any]  # [R', 6, 3]
    intersection_points = intersection_points[intersection_masks].view(-1, 2, 3)  # [R', 2, 3]
    # Ensure one or more positive ts.
    intersection_ts = intersection_ts[intersection_masks_any]  # [R', 6]
    intersection_ts = intersection_ts[intersection_masks]  # [R'*2]
    intersection_ts = intersection_ts.view(indices.size()[0], 2)  # [R', 2]
    positive_ts_mask = (intersection_ts >= 0)  # [R', 2]
    positive_ts_count = torch.sum(positive_ts_mask.int(), dim=1)  # [R']
    positive_ts_mask = (positive_ts_count >= 1)  # [R']
    intersection_points = intersection_points[positive_ts_mask]  # [R'', 2, 3]
    false_indices = indices[torch.logical_not(positive_ts_mask)]  # [R',]
    indices = indices[positive_ts_mask]  # [R'',]
    if len(false_indices) > 0:
        intersection_masks_any[false_indices.long()] = torch.zeros(false_indices.size(), dtype=torch.bool)
    return rays_start_point_in_box_frame, intersection_masks_any, intersection_points, indices


def compute_bounds_from_intersect_points(rays_o, intersect_indices,
                                         intersect_points):
    """Computes bounds from intersection points.

    Note: Make sure that inputs are in the same coordiante frame.

    Args:
        rays_o: [R, 3] float tensor
        intersect_indices: [R', 1] float tensor
        intersect_points: [R', 2, 3] float tensor

    Returns:
        intersect_bounds: [R', 2] float tensor

    where R is the number of rays and R' is the number of intersecting rays.
    """
    intersect_rays_o = rays_o[intersect_indices]  # [R', 1, 3]
    intersect_diff = intersect_points - intersect_rays_o  # [R', 2, 3]
    intersect_bounds = torch.norm(intersect_diff, dim=2)  # [R', 2]

    # Sort the bounds so that near comes before far for all rays.
    intersect_bounds, _ = torch.sort(intersect_bounds, dim=1)  # [R', 2]

    # For some reason the sort function returns [R', ?] instead of [R', 2], so we
    # will explicitly reshape it.
    intersect_bounds = intersect_bounds.view(-1, 2)  # [R', 2]
    return intersect_bounds


def compute_ray_bbox_bounds_pairwise(rays_o, rays_d, box_length,
                                     box_width, box_height, box_center,
                                     box_rotation, far_limit=1e10):
    """Computes near and far bounds for rays intersecting with bounding boxes.

    Note: rays and boxes are defined in world coordinate frame.

    Args:
        rays_o: [R, 3] float tensor. A set of ray origins.
        rays_d: [R, 3] float tensor. A set of ray directions.
        box_length: scalar or [R,] float tensor. Bounding box length.
        box_width: scalar or [R,] float tensor. Bounding box width.
        box_height: scalar or [R,] float tensor. Bounding box height.
        box_center: [3,] or [R, 3] float tensor. The center of the box.
        box_rotation: [3, 3] or [R, 3, 3] float tensor. The box rotation matrix.
        far_limit: float. The maximum far value to use.

    Returns:
        intersect_bounds: [R', 2] float tensor. The bounds per-ray, sorted in
            ascending order.
        intersect_indices: [R', 1] float tensor. The intersection indices.
        intersect_mask: [R,] float tensor. The mask denoting intersections.
    """
    # Compute ray destinations.
    normalized_rays_d = ray_utils.normalize_rays(rays=rays_d)
    rays_dst = rays_o + far_limit * normalized_rays_d

    # Transform the rays from world to box coordinate frame.
    rays_o_in_box_frame, intersect_mask, intersect_points_in_box_frame, intersect_indices = (  # [R,], [R', 2, 3], [R', 2]
        ray_box_intersection_pairwise(
            box_center=box_center,
            box_rotation_matrix=box_rotation,
            box_length=box_length,
            box_width=box_width,
            box_height=box_height,
            rays_start_point=rays_o,
            rays_end_point=rays_dst))
    intersect_indices = intersect_indices.unsqueeze(1).long()  # [R', 1]
    intersect_bounds = compute_bounds_from_intersect_points(
        rays_o=rays_o_in_box_frame,
        intersect_indices=intersect_indices,
        intersect_points=intersect_points_in_box_frame)
    return intersect_bounds, intersect_indices, intersect_mask
