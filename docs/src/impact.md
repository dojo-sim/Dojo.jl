# Impact

### Mathematical Model
We model hard contact via constraints on the system’s configuration and the applied contact forces. For a system with $P$ contact points, we define a signed-distance function,

$$\phi : \mathbf{Z} \rightarrow \mathbf{R}^P$$

subject to the following element-wise constraint:

$$ϕ(z) > 0,$$

Impact forces with magnitude $\gamma \in \mathbf{R}^P$ are applied to the bodies’ contact points in the direction of their surface normals in order to enforce (5) and prevent interpenetration. A non-negative constraint,

$$\gamma \geq 0,$$

enforces physical behavior that impulses are repulsive (e.g., the floor does not attract bodies), and the complementarity condition,

$$\gamma \circ \phi(z) = 0,$$

where $\circ$ is an element-wise product operator, enforces zero force if the body is not in contact and allows non-zero force during contact.


### Constructor
```@docs
ImpactContact
```
