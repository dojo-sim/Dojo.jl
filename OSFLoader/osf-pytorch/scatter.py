"""Utility functions for scattering rays."""
import torch


def create_scatter_indices_for_dim(dim, shape, indices=None):
    """Create scatter indices for a given dimension."""
    dim_size = shape[dim]
    N_dims = len(shape)
    reshape = [1] * N_dims
    reshape[dim] = -1

    if indices is None:
        indices = torch.arange(dim_size, dtype=torch.int)  # [dim_size,]

    indices = indices.view(reshape)

    indices = torch.broadcast_to(
        indices, shape)  # [Ro, S, 1] or [Ro, S, C, 1]  [0,1,1,1] vs. [512,64,1,1]

    indices = indices.int()
    return indices


def create_scatter_indices(updates, dim2known_indices):
    """Create scatter indices."""
    updates_expanded = updates.unsqueeze(-1)
    target_shape = updates_expanded.size()
    n_dims = len(updates.size())  # 2 or 3

    dim_indices_list = []
    for dim in range(n_dims):
        indices = None
        if dim in dim2known_indices:
            indices = dim2known_indices[dim]
        dim_indices = create_scatter_indices_for_dim(  # [Ro, S, C, 1]
            dim=dim,
            shape=target_shape,  # [Ro, S, 1] or [Ro, S, C, 1]
            indices=indices)  # [Ro,]
        dim_indices_list.append(dim_indices)
    scatter_indices = torch.cat((dim_indices_list), dim=-1)  # [Ro, S, C, 3]
    return scatter_indices
    

def scatter_nd(tensor, updates, dim2known_indices):
    scatter_indices = create_scatter_indices(  # [Ro, S, C, 3]
        updates=updates,  # [Ro, S]
        dim2known_indices=dim2known_indices)  # [Ro,]
    #scattered_tensor = (tensor[scatter_indices.view(-1, 3)[:, 0],
    #                           scatter_indices.view(-1, 3)[:, 1],
    #                           scatter_indices.view(-1, 3)[:, 2]] = updates.view(-1))
    #scattered_tensor = tensor[scatter_indices.view(-1, 2)[:, 0].long(),
    #                          scatter_indices.view(-1, 2)[:, 1].long()]
                              #scatter_indices.view(-1, 3)[:, 2]]
    if len(scatter_indices.shape) == 3:
        scatter_indices = scatter_indices.long().view(-1, 2)
        #tensor.scatter_(-1, scatter_indices.view(-1, 2).long(), updates)
        #tensor.index_copy_(-1, scatter_indices.view(-1, 2), updates.view(-1))
        tensor.data[scatter_indices[:, 0],
               scatter_indices[:, 1]] = updates.view(-1)
    elif len(scatter_indices.shape) == 4:
        scatter_indices = scatter_indices.long().view(-1, 3)
        tensor.data[scatter_indices[:, 0],
               scatter_indices[:, 1],
               scatter_indices[:, 2]] = updates.view(-1)
    # scattered_tensor = updates.view(-1)
    # return scattered_tensor
    #elif len(scatter_indices.shape) == 4:
    #    tensor.scatter_(-1, scatter_indices.long(), updates)
        #tensor.index_copy_(-1, scatter_indices.view(-1, 3), updates.view(-1))
    return tensor


def scatter_results(intersect, indices, N_rays, keys, N_samples, N_importance=None):
    """Scatters intersecting ray results into the original set of rays.

    Args:
        intersect: Dict. Values are tensors of intersecting rays which can be any of the
            following.
                z_vals: [R?, S]
                pts: [R?, S, 3]
                rgb: [R?, S, 3]
                raw: [R?, S, 4]
                alpha: [R?, S, 1]
        indices: [R?, 1] int tensor. Intersecting ray indices to scatter back to the full set
            of rays.
        N_rays: int or int tensor. Total number of rays.
        keys: [str]. List of keys from the 'intersect' dictionary to scatter.
        N_samples: [int]. Number of samples.
        N_importance: [int]. Number of importance (fine) samples.

    Returns:
        scattered_results: Dict. Scattered results, where each value is of shape [R, ...].
            The original intersecting ray results are padding with samples with zero density,
            so they won't have any contribution to the final render.
    """
    # We use 'None' to indicate that the intersecting set of rays is equivalent to
    # the full set if rays, so we are done.
    if indices is None:
        return {k: intersect[k] for k in keys}

    scattered_results = {}
    # N_samples = intersect['z_vals'].shape[1]
    dim2known_indices = {0: indices}  # [R', 1]
    for k in keys:
        if k == 'z_vals':
            # tensor = torch.rand((N_rays, N_samples), dtype=torch.float)
            tensor = torch.arange(N_samples)  # [S,]
            tensor = tensor.float()
            tensor = torch.stack([tensor] * N_rays)  # [R, S]
        elif k == 'z_samples':
            # tensor = torch.rand((N_rays, N_samples), dtype=torch.float)  #[R, S]
            # N_importance = intersect[k].size()[1]
            tensor = torch.arange(N_importance)  # [I,]
            tensor = tensor.float()
            tensor = torch.stack([tensor] * N_rays)  # [R, I]
        elif k == 'raw':
            tensor = torch.full((N_rays, N_samples, 4), 1000.0, dtype=torch.float, requires_grad=True)  # [R, S, 4]
        elif k == 'pts':
            tensor = torch.full((N_rays, N_samples, 3), 1000.0, dtype=torch.float, requires_grad=True)  # [R, S, 3]
        elif 'rgb' in k:
            tensor = torch.zeros((N_rays, N_samples, 3), dtype=torch.float, requires_grad=True)  # [R, S, 3]
        elif 'alpha' in k:
            tensor = torch.zeros((N_rays, N_samples, 1), dtype=torch.float, requires_grad=True)  # [R, S, 1]
        else:
            raise ValueError(f'Invalid key: {k}')
        # No intersections to scatter.
        if len(indices) == 0:
            scattered_results[k] = tensor
        else:
            scattered_v = scatter_nd(  # [R, S, K]
                tensor=tensor,
                updates=intersect[k],  # [Ro, S]
                dim2known_indices=dim2known_indices)
            # Convert the batch dimension to a known dimension.
            # For some reason 'scattered_z_vals' becomes [R, ?]. We need  to explicitly
            # reshape it with 'N_sampels'.
            if k == 'z_samples':
                scattered_v = scattered_v.view(N_rays, N_importance)  # [R, I]
            else:
                if k == 'z_vals':
                    scattered_v = scattered_v.view(N_rays, N_samples)  # [R, S]
                else:
                    # scattered_v = scattered_v.view((N_rays,) + scattered_v.size()[1:])  # [R, S, K]
                    # scattered_v = scattered_v.view((-1,) + scattered_v.size()[1:])  # [R, S, K]
                    scattered_v = scattered_v.view(N_rays, N_samples, tensor.size()[2])  # [R, S, K]
            scattered_results[k] = scattered_v
    return scattered_results


def scatter_coarse_and_fine(
    ret0, ret, indices, num_rays, N_samples, N_importance, **kwargs):
    # TODO: only process raw if retraw=True.
    ret0 = scatter_results(
        intersect=ret0,
        indices=indices,
        N_rays=num_rays,
        keys=['z_vals', 'rgb', 'alpha', 'raw'],
        N_samples=N_samples)
    ret = scatter_results(
        intersect=ret,
        indices=indices,
        N_rays=num_rays,
        keys=['z_vals', 'rgb', 'alpha', 'raw', 'z_samples'],
        N_samples=N_samples + N_importance,
        N_importance=N_importance)
    return ret0, ret
