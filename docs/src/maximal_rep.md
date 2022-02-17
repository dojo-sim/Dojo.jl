# Maximal Coordinates 

The ``i``-th body in a mechanism with ``N`` bodies has state:  

``z^{(i)} = (x^{(i)}, v^{(i)}, q^{(i)}, \omega^{(i)}) \in \mathbf{R}^3 \times \mathbf{R}^3 \times \mathbf{H} \times \mathbf{R}^3``,  

represented in maximal coordinates, where ``\mathbf{H}`` is the space of unit quaternions. 

- ``x``: position in world frame
- ``v``: linear velocity in the world frame
- ``q``: orientation represented as a unit quaternion
- ``\omega``: angular velocity in the body frame

The mechanism state:   

``z = (z^{(1)}, \dots, z^{(N)})``.

is the concatentation of all body states.

