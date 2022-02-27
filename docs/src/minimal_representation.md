# Minimal Coordinates

Dojo simulates systems in [maximal coordinates](maximal_representation.md). 

For a mechanism with ``M`` joints and ``N`` bodies, the maximal representation ``z`` can be efficiently converted to minimal coordinates: 

``y = (y^{(1)}, \dots, y^{(M)}) \leftarrow z = (z^{(1)}, \dots, z^{(N)})``,

where ``y^{(j)}`` depends on the degree and type of joint. **Note**: this minimal representation does not stack coordinates followed by velocities, which is a common convention; instead, **coordinates and velocities are grouped by joint**.

Each minimal state comprises:

``y = (p_{\text{translational}}, p_{\text{rotational}}, w_{\text{translational}}, w_{\text{rotational}})``

coordinates ``p`` and velocities ``w`` for both translational and rotational degrees of freedom.

In the case of a floating-base joint, the minimal-representation orientation is converted to [modified Rodrigues parameters](https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions) from a unit quaternion.