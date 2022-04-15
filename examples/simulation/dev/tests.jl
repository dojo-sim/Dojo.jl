using LinearAlgebra
using Plots 

# sphere to capsule 

# capsule points
ca = [0.0; 0.0; 0.0] 
cb = [0.0; 0.0; 1.0] 

# sphere point
s = [0.0; 1.0; -1.0] 

dab = cb - ca
das = s - ca

len = dot(dab, dab)

t = dot(das, dab) ./ dot(dab, dab)
t_clamp = min(max(t, 0.0), 1.0)

# closest point 
if t <= 0.0 
    p = ca 
elseif t >= 1.0 
    p = cb 
else
    p = ca + t * dab
end


###
# https://stackoverflow.com/questions/44824512/how-to-find-the-closest-point-on-a-right-rectangular-prism-3d-rectangle/44824522#44824522
using LinearAlgebra
p = [0.0; 1.0; 1.0]
origin = [0.0; 0.0; 0.0] 
r = 0.5 
v1 = [r; 0.0; 0.0]
v2 = [0.0; r; 0.0] 
v3 = [0.0; 0.0; r] 

tx = dot(p - origin, v1) / dot(v1, v1) 
ty = dot(p - origin, v2) / dot(v2, v2) 
tz = dot(p - origin, v3) / dot(v3, v3)

tx_clamp = min(max(tx, 0.0), 1.0)
ty_clamp = min(max(ty, 0.0), 1.0)
tz_clamp = min(max(tz, 0.0), 1.0)

cp = tx_clamp * v1 + ty_clamp * v2 + tz_clamp * v3 + origin

