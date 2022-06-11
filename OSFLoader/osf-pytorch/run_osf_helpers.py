import os
import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import imageio
import json

# from torchsearchsorted import searchsorted

import ray_utils


# Misc utils

def img2mse(x, y): return torch.mean((x - y) ** 2)


def mse2psnr(x): return -10. * torch.log(x) / torch.log(torch.tensor([10.]))


def to8b(x): return (255*np.clip(x,0,1)).astype(np.uint8)


# Positional encoding

class Embedder:

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.create_embedding_fn()

    def create_embedding_fn(self):
        embed_fns = []
        d = self.kwargs['input_dims']
        out_dim = 0
        if self.kwargs['include_input']:
            embed_fns.append(lambda x: x)
            out_dim += d

        max_freq = self.kwargs['max_freq_log2']
        N_freqs = self.kwargs['num_freqs']

        if self.kwargs['log_sampling']:
            freq_bands = 2.**torch.linspace(0., max_freq, steps=N_freqs)
        else:
            freq_bands = torch.linspace(2.**0., 2.**max_freq, steps=N_freqs)

        for freq in freq_bands:
            for p_fn in self.kwargs['periodic_fns']:
                embed_fns.append(lambda x, p_fn=p_fn, freq=freq : p_fn(x * freq))
                out_dim += d

        self.embed_fns = embed_fns
        self.out_dim = out_dim

    def embed(self, inputs):
        return torch.cat([fn(inputs) for fn in self.embed_fns], -1)


def get_embedder(multires, i=0):
    if i == -1:
        return nn.Identity(), 3
    
    embed_kwargs = {
                'include_input' : True,
                'input_dims' : 3,
                'max_freq_log2' : multires-1,
                'num_freqs' : multires,
                'log_sampling' : True,
                'periodic_fns' : [torch.sin, torch.cos],
    }
    
    embedder_obj = Embedder(**embed_kwargs)
    embed = lambda x, eo=embedder_obj : eo.embed(x)
    return embed, embedder_obj.out_dim


# Model architecture

class NeRF(nn.Module):
    def __init__(self, D=8, W=256, input_ch=3, input_ch_views=3, input_ch_lights=3, output_ch=4, skips=[4], use_viewdirs=False, use_lightdirs=False):
        """
        """
        super(NeRF, self).__init__()
        self.D = D
        self.W = W
        self.input_ch = input_ch
        self.input_ch_views = input_ch_views
        self.input_ch_lights = input_ch_lights
        self.skips = skips
        self.use_viewdirs = use_viewdirs
        self.use_lightdirs = use_lightdirs

        self.pts_linears = nn.ModuleList(
            [nn.Linear(input_ch, W)] + [nn.Linear(W, W) if i not in self.skips else nn.Linear(W + input_ch, W) for i in range(D-1)])

        if use_viewdirs and use_lightdirs:
            self.bottleneck_linear = nn.Linear(W, W)
            self.alpha_linear = nn.Linear(W, 1)
            self.rgb_linear = nn.Linear(W//2, 3)
            self.views_lights_linears = nn.ModuleList([nn.Linear(input_ch_views + input_ch_lights + W, W//2)])
        elif use_viewdirs:
            self.bottleneck_linear = nn.Linear(W, W)
            self.alpha_linear = nn.Linear(W, 1)
            self.rgb_linear = nn.Linear(W//2, 3)
            self.views_lights_linears = nn.ModuleList([nn.Linear(input_ch_views + W, W//2)])
        elif use_lightdirs:
            self.bottleneck_linear = nn.Linear(W, W)
            self.alpha_linear = nn.Linear(W, 1)
            self.rgb_linear = nn.Linear(W//2, 3)
            self.views_lights_linears = nn.ModuleList([nn.Linear(input_ch_lights + W, W//2)])
        else:
            self.output_linear = nn.Linear(W, output_ch)

    def forward(self, x):
        inputs_pts, inputs_views, inputs_lights = torch.split(x, [self.input_ch, self.input_ch_views, self.input_ch_lights], dim=-1)
        outputs = inputs_pts
        for i, l in enumerate(self.pts_linears):
            outputs = self.pts_linears[i](outputs)
            outputs = F.relu(outputs)
            if i in self.skips:
                outputs = torch.cat([inputs_pts, outputs], -1)

        if self.use_viewdirs or self.use_lightdirs:
            alpha = self.alpha_linear(outputs)
            bottleneck = self.bottleneck_linear(outputs)
            inputs_dirs = bottleneck
            if self.use_viewdirs:
                inputs_dirs = torch.cat([inputs_dirs, inputs_views], -1)  # concat viewdirs
            if self.use_lightdirs:
                inputs_dirs = torch.cat([inputs_dirs, inputs_lights], -1)  # concat lightdirs
            outputs = inputs_dirs
            for i, l in enumerate(self.views_lights_linears):
                outputs = self.views_lights_linears[i](outputs)
                outputs = F.relu(outputs)
            outputs = self.rgb_linear(outputs)
            outputs = torch.cat([outputs, alpha], -1)
        else:
            outputs = self.output_linear(outputs)

        return outputs


def get_rays(H, W, focal, c2w, img_id=0):
    """Get ray origins, directions from a pinhole camera.

    Args:
        H: int. Image height.
        W: int. Image width.
        focal: Focal length of the camera.
        c2w: [4, 4] np.float32. Camera-to-world extrinsics matrix.
        img_id: int. ID of the image.

    Returns:
        rays_o: [H, W, 3] float tensor. Ray origins.
        rays_d: [H, W, 3] float tensor. Ray directions.
        rays_i: [H, W, 1] float tensor. Image ID.
    """
    i, j = torch.meshgrid(torch.linspace(0, W-1, W), torch.linspace(0, H-1, H))  # pytorch's meshgrid has indexing='ij'
    i = i.t()
    j = j.t()
    dirs = torch.stack([(i-W*.5)/focal, -(j-H*.5)/focal, -torch.ones_like(i)], -1)  # [H, W, 3]
    # Rotate ray directions from camera frame to the world frame
    c2w = torch.tensor(c2w)
    rays_d = torch.sum(dirs[..., None, :] * c2w[:3,:3], -1)  # dot product, equals to: [c2w.dot(dir) for dir in dirs]  # [H, W, 3]
    # Translate camera frame's origin to the world frame. It is the origin of all rays.
    rays_o = c2w[:3,-1].expand(rays_d.shape)  # [H, W, 3]
    rays_i = torch.full(rays_d.size(), img_id).float()  # [H, W, 3]
    return rays_o, rays_d, rays_i


def get_rays_np(H, W, focal, c2w, img_id=0):
    """Get ray origins, directions from a pinhole camera.

    Args:
        H: int. Image height.
        W: int. Image width.
        focal: Focal length of the camera.
        c2w: [4, 4] np.float32. Camera-to-world extrinsics matrix.
        img_id: int. ID of the image.

    Returns:
        rays_o: [H, W, 3] np.float32. Ray origins.
        rays_d: [H, W, 3] np.float32. Ray directions.
        rays_i: [H, W, 3] np.float32. Image ID.
    """
    i, j = np.meshgrid(np.arange(W, dtype=np.float32),
                       np.arange(H, dtype=np.float32), indexing='xy')
    dirs = np.stack([(i-W*.5)/focal, -(j-H*.5)/focal, -np.ones_like(i)], -1)
    rays_d = np.sum(dirs[..., np.newaxis, :] * c2w[:3, :3], -1)
    rays_o = np.broadcast_to(c2w[:3, -1], np.shape(rays_d))
    rays_i = np.full_like(rays_d, img_id, dtype=np.float32)  # [H, W, 3]
    return rays_o, rays_d, rays_i


def ndc_rays(H, W, focal, near, rays_o, rays_d):
    """Normalized device coordinate rays.
    Space such that the canvas is a cube with sides [-1, 1] in each axis.
    Args:
      H: int. Height in pixels.
      W: int. Width in pixels.
      focal: float. Focal length of pinhole camera.
      near: float or array of shape[batch_size]. Near depth bound for the scene.
      rays_o: array of shape [batch_size, 3]. Camera origin.
      rays_d: array of shape [batch_size, 3]. Ray direction.
    Returns:
      rays_o: array of shape [batch_size, 3]. Camera origin in NDC.
      rays_d: array of shape [batch_size, 3]. Ray direction in NDC.
    """
    # Shift ray origins to near plane
    t = -(near + rays_o[...,2]) / rays_d[...,2]
    rays_o = rays_o + t[...,None] * rays_d

    # Projection
    o0 = -1./(W/(2.*focal)) * rays_o[...,0] / rays_o[...,2]
    o1 = -1./(H/(2.*focal)) * rays_o[...,1] / rays_o[...,2]
    o2 = 1. + 2. * near / rays_o[...,2]

    d0 = -1./(W/(2.*focal)) * (rays_d[...,0]/rays_d[...,2] - rays_o[...,0]/rays_o[...,2])
    d1 = -1./(H/(2.*focal)) * (rays_d[...,1]/rays_d[...,2] - rays_o[...,1]/rays_o[...,2])
    d2 = -2. * near / rays_o[...,2]

    rays_o = torch.stack([o0,o1,o2], -1)
    rays_d = torch.stack([d0,d1,d2], -1)

    return rays_o, rays_d


# Hierarchical sampling (section 5.2)
def sample_pdf(bins, weights, N_samples, det=False, pytest=False):
    # Get pdf
    weights = weights + 1e-5 # prevent nans
    pdf = weights / torch.sum(weights, -1, keepdim=True)
    cdf = torch.cumsum(pdf, -1)
    cdf = torch.cat([torch.zeros_like(cdf[...,:1]), cdf], -1)  # (batch, len(bins))

    # Take uniform samples
    if det:
        u = torch.linspace(0., 1., steps=N_samples)
        u = u.expand(list(cdf.shape[:-1]) + [N_samples])
    else:
        u = torch.rand(list(cdf.shape[:-1]) + [N_samples])

    # Pytest, overwrite u with numpy's fixed random numbers
    if pytest:
        np.random.seed(0)
        new_shape = list(cdf.shape[:-1]) + [N_samples]
        if det:
            u = np.linspace(0., 1., N_samples)
            u = np.broadcast_to(u, new_shape)
        else:
            u = np.random.rand(*new_shape)
        u = torch.Tensor(u)

    # Invert CDF
    u = u.contiguous()
    # inds = searchsorted(cdf, u, side='right')
    inds = torch.searchsorted(cdf, u, right=True)
    below = torch.max(torch.zeros_like(inds-1), inds-1)
    above = torch.min((cdf.shape[-1]-1) * torch.ones_like(inds), inds)
    inds_g = torch.stack([below, above], -1)  # (batch, N_samples, 2)

    # cdf_g = tf.gather(cdf, inds_g, axis=-1, batch_dims=len(inds_g.shape)-2)
    # bins_g = tf.gather(bins, inds_g, axis=-1, batch_dims=len(inds_g.shape)-2)
    matched_shape = [inds_g.shape[0], inds_g.shape[1], cdf.shape[-1]]
    cdf_g = torch.gather(cdf.unsqueeze(1).expand(matched_shape), 2, inds_g)
    bins_g = torch.gather(bins.unsqueeze(1).expand(matched_shape), 2, inds_g)

    denom = (cdf_g[...,1]-cdf_g[...,0])
    denom = torch.where(denom<1e-5, torch.ones_like(denom), denom)
    t = (u-cdf_g[...,0])/denom
    samples = bins_g[...,0] + t * (bins_g[...,1]-bins_g[...,0])

    return samples


# Function for computing density from model prediction. This value is
# strictly between [0, 1].
def raw2alpha(raw, dists):
    return 1.0 - torch.exp(-F.relu(raw) * dists)


def compute_alpha(
    z_vals, raw_alpha, rays_d, raw_noise_std, last_dist_method='infinity'):
    """Normalizes raw sigma predictions from the network into normalized alpha.

    Args:
        z_vals: [R, S] float tensor. Integration time.
        raw_alpha: [R, S] float tensor. Raw alpha predictions from model.
        rays_d: [R, 3] float tensor. Ray directions.
        raw_noise_std: Amount of noise to apply to raw alpha predictions.
        last_dist_method: Method to use to deal with the last dist.

    Returns:
        alpha: [R, S, 1] float tensor.
    """
    # Compute 'distance' (in time) between each integration time along a ray.
    dists = z_vals[..., 1:] - z_vals[..., :-1]  # [R, S-1]

    # The 'distance' from the last integration time is infinity.
    if last_dist_method == 'infinity':
        dists = torch.cat([dists, torch.broadcast_to(torch.tensor([1e10]), dists[..., :1].size())], dim=-1)  # [R, S]
    elif last_dist_method == 'last':
        dists = torch.cat([dists, dists[..., -1:]], dim=-1)  # [R, S]

    # Multiply each distance by the norm of its corresponding direction ray
    # to convert to real world distance (accounts for non-unit directions).
    dists = dists * torch.norm(rays_d[..., None, :], dim=-1)  # [R, S]

    # raw_alpha = raw_alpha.squeeze(-1)  # [R, S]

    # Add noise to model's predictions for density. Can be used to
    # regularize network during training (prevents floater artifacts).
    noise = 0.
    if raw_noise_std > 0.:
        noise = torch.rand(raw_alpha.size()) * raw_noise_std  # [R, S]

    # Convert from raw alpha to alpha values between [0, 1].
    # Predict density of each sample along each ray. Higher values imply
    # higher likelihood of being absorbed at this point.
    alpha = raw2alpha(raw_alpha + noise, dists)  # [R, S]
    return alpha.unsqueeze(-1)  # [R, S, 1]


def normalize_rgb(raw_rgb, scaled_sigmoid):
    """Normalizes raw RGB values from the network.

    Args:
        raw_rgb: [R, S, 3] tf.float32. Raw RGB values.
        scaled_sigmoid: Whether to apply scaled sigmoid.

    Returns:
        rgb: [R, S, 3] tf.float32. Normalized RGB values.
    """
    # Extract RGB of each sample position along each ray.
    rgb = torch.sigmoid(raw_rgb)  # [R, S, 3]
    if scaled_sigmoid:
        rgb = 1.2 * (rgb - 0.5) + 0.5  # [R, S, 3]
    return rgb


def normalize_raw(raw, z_vals, rays_d, raw_noise_std, scaled_sigmoid):
    """Normalize raw outputs of the network.

    Args:
        raw: [R, S, 4] tf.float32. Raw model predictions.
        z_vals: [R, S], tf.float32. Integration time.
        raw_noise_std: Amount of noise to apply to raw alpha predictions.
        scaled_sigmoid: Whether to apply scaled sigmoid.

    Returns:
        rgb: [R, S, 3] tf.float32. Normalized rgb.
        alpha: [R, S, 1] tf.float32. Normalized alpha.
    """
    rgb = normalize_rgb(raw_rgb=raw[..., :3], scaled_sigmoid=scaled_sigmoid)
    alpha = compute_alpha(
        z_vals=z_vals, raw_alpha=raw[..., 3], rays_d=rays_d, raw_noise_std=raw_noise_std)
    return rgb, alpha


def compute_transmittance(alpha):
    """Compute transmittance from alpha values.

    Args:
        alpha: [R, S] float tensor. Alpha predictions from the model.
    Returns:
        trans: [R, S] float tensor. Transmittance values.
    """
    trans_tmp = torch.cumprod(1.-alpha + 1e-10, dim=-1)
    trans = torch.ones_like(trans_tmp)
    trans[:, 1:] = trans_tmp[:, :-1]
    return trans


def compute_weights(alpha):
    """Computes weights from alpha.
    Args:
        alpha: [R, S] float tensor. Normalized density values.
    Returns:
        weights: [R, S] float tensor. Weights per sample.
    """
    weights = alpha * compute_transmittance(alpha)  # [R, S]
    return weights


def compose_outputs(z_vals, rgb, alpha, white_bkgd):
    """Transforms model's predictions to semantically meaningful values.
    Args:
        z_vals: [R, S] float tensor. Integration time.
        rgb: [R, S, 3] float tensor. Normalized RGB values.
        alpha: [R, S, 1] float tensor. Alpha (normalized) density.
        white_bkgd: bool. If True, assume a white background.
    Returns:
        outputs: [Dict]. Dictionary of outputs.
            rgb_map: [R, 3]. Estimated RGB color of a ray.
            disp_map: [R]. Disparity map. Inverse of depth map.
            acc_map: [R]. Sum of weights along each ray.
            weights: [R, num_samples]. Weights assigned to each sampled color.
            depth_map: [R]. Estimated distance to object.
            rgb: [R, S, 3]. Same as input `rgb`.
            alpha: [R, S, 1]. Same as input `alpha`.
    """
    # Compute weight for RGB of each sample along each ray.  A cumprod() is
    # used to express the idea of the ray not having reflected up to this
    # sample yet.
    weights = compute_weights(alpha=alpha[..., 0])  # [R, S]

    # Computed weighted color of each sample along each ray.
    # print(weights.shape)
    # print(alpha.shape)
    # print(rgb.shape)
    rgb_map = torch.sum(weights.unsqueeze(-1) * rgb, dim=-2)  # [N_rays, 3]

    # Estimated depth map is expected distance.
    depth_map = torch.sum(weights * z_vals, dim=-1)

    # Disparity map is inverse depth.
    disp_map = 1. / torch.maximum(torch.tensor(1e-10), depth_map / torch.sum(weights, dim=-1))

    # Sum of weights along each ray. This value is in [0, 1] up to numerical error.
    acc_map = torch.sum(weights, -1)

    # To composite onto a white background, use the accumulated alpha map.
    if white_bkgd:
        rgb_map = rgb_map + (1.-acc_map[..., None])

    outputs = {
        'rgb_map': rgb_map,
        'disp_map': disp_map,
        'acc_map': acc_map,
        'weights': weights,
        'depth_map': depth_map,
    }
    return outputs


def default_ray_sampling(ray_batch, N_samples, lindisp, perturb):
    """Samples coarse samples along rays.
    Args:
        ray_batch: array of shape [batch_size, ...]. All information necessary
            for sampling along a ray, including: ray origin, ray direction, min
            dist, max dist, and unit-magnitude viewing direction.
        N_samples: int. Number of different times to sample along each ray.
        lindisp: bool. If True, sample linearly in inverse depth rather than in depth.
        perturb: float, 0 or 1. If non-zero, each ray is sampled at stratified
            random points in time.

    Returns:
        z_vals: [num_rays, num_samples along ray]. Integration time.
        pts: Sampled points.
    """
    # batch size
    N_rays = ray_batch.size()[0]

    # Extract ray origin, direction.
    rays_o, rays_d = ray_batch[:, 0:3], ray_batch[:, 3:6]  # [N_rays, 3] each

    # Extract lower, upper bound for ray distance.
    bounds = ray_batch[..., 6:8].view(-1, 1, 2)
    near, far = bounds[..., 0], bounds[..., 1]  # [-1,1]

    # Decide where to sample along each ray. Under the logic, all rays will be sampled at
    # the same times.
    t_vals = torch.linspace(0., 1., N_samples)
    if not lindisp:
        # Space integration times linearly between 'near' and 'far'. Same
        # integration points will be used for all rays.
        z_vals = near * (1.-t_vals) + far * (t_vals)
    else:
        near += 1e-10
        far += 1e-10
        # Sample linearly in inverse depth (disparity).
        z_vals = 1./(1./near * (1.-t_vals) + 1./far * (t_vals))
    z_vals = torch.broadcast_to(z_vals, (N_rays, N_samples))

    # Perturb sampling time along each ray.
    if perturb > 0.:
        # get intervals between samples
        mids = .5 * (z_vals[..., 1:] + z_vals[..., :-1])
        upper = torch.cat([mids, z_vals[..., -1:]], -1)
        lower = torch.cat([z_vals[..., :1], mids], -1)
        # stratified samples in those intervals
        t_rand = torch.rand(z_vals.size())
        z_vals = lower + (upper - lower) * t_rand

    # Points in space to evaluate model at.
    pts = rays_o[..., None, :] + rays_d[..., None, :] * \
        z_vals[..., :, None]  # [N_rays, N_samples, 3]
    return z_vals, pts


def default_ray_sampling_fine(ray_batch, z_vals, weights, perturb, N_importance):
    """Samples fine samples along rays.
    Args:
      ray_batch: array of shape [batch_size, ...]. All information necessary
        for sampling along a ray, including: ray origin, ray direction, min
        dist, max dist, and unit-magnitude viewing direction.
      z_vals: [num_rays, N_samples]. Coarse integration time.
      weights: [num_rays, num_samples]. Weights assigned to each sampled color.
      perturb: float, 0 or 1. If non-zero, each ray is sampled at stratified
        random points in time.
      N_importance: int. Number of additional times to sample along each ray.
        These samples are only passed to network_fine.

    Returns:
      z_vals: [num_rays, N_samples + N_importance]. Coarse and fine integration time.
      z_samples: [num_rays, N_importance]. Fine integration time.
      pts: Coarse and fine sampled points.
    """
    # Extract ray origin, direction.
    rays_o, rays_d = ray_batch[:, 0:3], ray_batch[:, 3:6]  # [N_rays, 3] each

    # Obtain additional integration times to evaluate based on the weights
    # assigned to colors in the coarse model.
    z_vals_mid = .5 * (z_vals[..., 1:] + z_vals[..., :-1])
    z_samples = sample_pdf(
        z_vals_mid, weights[..., 1:-1], N_importance, det=(perturb == 0.))
    z_samples = z_samples.detach()

    # Obtain all points to evaluate color, density at.
    z_vals, _ = torch.sort(torch.cat((z_vals, z_samples), -1), -1)
    pts = rays_o[..., None, :] + rays_d[..., None, :] * \
        z_vals[..., :, None]  # [N_rays, N_samples + N_importance, 3]
    return z_vals, z_samples, pts


def get_dirs(ray_batch, pts, metadata, use_viewdirs, use_lightdirs, lightdirs_method):
    """Get ray directions.
    Args:
        ray_batch: [R, M] float tensor. All information necessary for sampling along a
            ray, including: ray origin, ray direction, min dist, max dist, and
            unit-magnitude viewing direction, all in object coordinate frame.
        pts: [R, S, 3] float tensor. Sampled points along rays.
        metadata: [N, 3] float tensor. Metadata about each image. Currently only light
            position is provided.
        use_viewdirs: Whether to use view directions.
        use_lightdirs: Whether to use light directions.
        lightdirs_method: Method for computing lightdirs.
    """
    viewdirs, lightdirs = None, None
    if use_viewdirs:
        assert ray_batch.size()[-1] > 8
        viewdirs = ray_batch[:, 8:11]  # [R, 3]
        viewdirs = torch.broadcast_to(viewdirs[:, None], pts.size())  # [R, S, 3]
    if use_lightdirs:
        # Use viewdirs as lightdirs.
        if lightdirs_method == 'viewdirs':
            assert viewdirs is not None, "viewdirs is None"
            lightdirs = viewdirs  # [R, S, 3]
        # Compute lightdirs based on ray metadata or randomly sample directions.
        else:
            rays_i = ray_batch[:, -1:]  # [R, 1]
            lightdirs = ray_utils.get_lightdirs(  # [R, S, 3]
                    lightdirs_method=lightdirs_method, num_rays=pts.size()[0],
                    num_samples=pts.size()[1], rays_i=rays_i, metadata=metadata,
                    normalize=False)
    return viewdirs, lightdirs


def run_single_object(
        ray_batch,
        metadata,
        network_fn,
        network_query_fn,
        use_viewdirs,
        use_lightdirs,
        lightdirs_method,
        N_samples,
        lindisp=False,
        perturb=0.,
        N_importance=0,
        network_fine=None,
        raw_noise_std=0.,
        scaled_sigmoid=False,
        retraw=True,
        verbose=False,
        **kwargs):
    """
    Args:
        ray_batch: [R, M] float tensor. All information necessary for sampling along a
            ray, including: ray origin, ray direction, min dist, max dist, and
            unit-magnitude viewing direction, all in object coordinate frame.
        metadata: [N, 3] float tensor. Metadata about each image. Currently only light
            position is provided.
        network_fn: function. Model for predicting RGB and density at each point
            in space.
        network_query_fn: function used for passing queries to network_fn.
        use_viewdirs: Whether to use view directions.
        use_lightdirs: Whether to use light directions.
        lightdirs_method: Method for computing lightdirs.
        N_samples: int. Number of different times to sample along each ray.
        lindisp: bool. If True, sample linearly in inverse depth rather than in depth.
        perturb: float, 0 or 1. If non-zero, each ray is sampled at stratified
            random points in time.
        N_importance: int. Number of additional times to sample along each ray.
            These samples are only passed to network_fine.
        network_fine: "fine" network with same spec as network_fn.
        raw_noise_std: ...
        scaled_sigmoid: ...
        retraw: bool. If True, include model's raw, unprocessed predictions.
        verbose: bool. If True, print more debugging info.

    Returns:
        ret: [Dict]. Results from the fine model, containing the following items:
            z_vals
            rgb
            alpha
            If retraw is True, we also provide:
                raw
        ret0: [Dict]. Results from the coarse model.
                z_vals
                rgb
                alpha
                z_samples
            If retraw is True, we also provide:
                raw
                pts
    """
    # Extract ray origin, direction.
    rays_o, rays_d = ray_batch[:, 0:3], ray_batch[:, 3:6]  # [R, 3] each

    # Sample coarse points along each ray.
    z_vals, pts = default_ray_sampling(ray_batch, N_samples, lindisp, perturb)

    # Extract unit-normalized viewing and lighting directions.
    viewdirs, lightdirs = get_dirs(
            ray_batch, pts, metadata, use_viewdirs, use_lightdirs, lightdirs_method)

    # Evaluate model at each point.
    raw = network_query_fn(pts, viewdirs, lightdirs, network_fn)  # [R, S, 4]
    rgb, alpha = normalize_raw(raw, z_vals, rays_d, raw_noise_std, scaled_sigmoid)
    ret = {
        'z_vals': z_vals,
        'rgb': rgb,
        'alpha': alpha,
    }
    if retraw:
        ret['raw'] = raw
        # ret['pts'] = pts

    ret0 = {}
    if N_importance > 0:
        # Store coarse model outputs.
        ret0 = ret

        # Sample fine points along each ray.
        weights = compute_weights(alpha=alpha[..., 0])  # [R, S]
        z_vals, z_samples, pts = default_ray_sampling_fine(
            ray_batch, z_vals, weights, perturb, N_importance)
        viewdirs, lightdirs = get_dirs(
                ray_batch, pts, metadata, use_viewdirs, use_lightdirs, lightdirs_method)

        # Make predictions with network_fine.
        run_fn = network_fn if network_fine is None else network_fine
        raw = network_query_fn(pts, viewdirs, lightdirs, run_fn)
        rgb, alpha = normalize_raw(raw, z_vals, rays_d, raw_noise_std, scaled_sigmoid)
        ret = {
            'z_vals': z_vals,
            'rgb': rgb,
            'alpha': alpha,
            'z_samples': z_samples,
        }
        if retraw:
            ret['raw'] = raw
            ret['pts'] = pts
    return ret0, ret


def combine_multi_object_results(name2results):
    """Combines network outputs from multiple objects.
    Args:
        name2results: Dictionary mapping from object name to object results. Object
            results must contain `z_vals` with shape [R, S]. All remaining results must
            have the shape [R, S, K].
        keys: [str]. A list of keys to combine results over.
    Returns:
        results: Dict. Combined results.
    """
    # Collect z values across all objects.
    z_vals_list = []
    for _, results in name2results.items():
        z_vals_list.append(results['z_vals'])

    # Concatenate lists of object results into a single tensor.
    z_vals = torch.cat(z_vals_list, dim=-1)  # [R, S*O]

    # Compute the argsort indices.
    z_argsort_indices = torch.argsort(z_vals, -1)  # [R, S*O]
    N_rays, N_samples = z_vals.size()[0], z_vals.size()[1]
    gather_indices = torch.arange(N_rays).unsqueeze(-1)  # [R, 1]
    gather_indices = torch.tile(gather_indices, (1, N_samples))  # [R, S]
    gather_indices = torch.cat(
        [gather_indices[..., None], z_argsort_indices[..., None]], dim=-1)

    results = {}
    for k in ['z_vals', 'rgb', 'alpha', 'raw', 'ray_batch', 'shadow_ray_batch']:
        if k == 'z_vals':
            v_combined = z_vals
        else:
            v_list = [r[k] for r in name2results.values() if k in r]
            if len(v_list) > 0:
                v_combined = torch.cat(v_list, dim=1)  # [R, S*O, K]
            else:
                v_combined = None
        if v_combined is not None:
            if k in ['ray_batch', 'shadow_ray_batch']:
                results[k] = v_combined  # [R, M*O]
            else:
                # Sort the tensors.
                v_sorted = v_combined[gather_indices.view(-1, 2)[:, 0], gather_indices.view(-1, 2)[:, 1]].view(v_combined.size())  # [R, S, K]
                results[k] = v_sorted
    return results
