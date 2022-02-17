# Minimal Coordinates

Dojo simulates systems in [maximal coordinates](maximal_rep.md). 

For a mechanism with ``M`` joints and ``N`` bodies, the maximal representation ``z`` can be efficiently converted to minimal coordinates: 

``x = (x^{(1)}, \dots, x^{(M)}) \leftarrow z = (z^{(1)}, \dots, z^{(N)})``,

where ``x^{(j)}`` depends on the degree and type of joint. 

In the case of a floating-base joint, the minimal-representation orientation is converted to [modified Rodrigues parameters](https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions) from a unit quaternion.