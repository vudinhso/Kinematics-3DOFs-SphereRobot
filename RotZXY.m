function R = RotZXY(q)
# Rotation matrix from YXY convention
    phi1 = q(1); phi2 = q(2); phi3 = q(3); 
    x = [1,0,0]';
    y = [0,1,0]';
    z = [0,0,1]';
    R = RotMatrix(z, phi1) * ...
        RotMatrix(x, phi2) * ...
        RotMatrix(y, phi3); 
end
