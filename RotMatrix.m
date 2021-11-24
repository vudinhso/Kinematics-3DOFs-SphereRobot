function R = RotMatrix(u,th)
    # Build the rotation matrix from a unit vector u and angle th
    if NearZero(th)
        R = eye(3);
    else
        w_mat = VecToso3(u);
        R = eye(3) + sin(th) * w_mat + (1 - cos(th)) * w_mat * w_mat;
    end
end
