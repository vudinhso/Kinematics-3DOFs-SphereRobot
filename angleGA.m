function q = angleGA(R)
# get parameter q1, q2, q3, q4 (bivector from GA - quaternions) from rotation matrix
    [w, th] = RotLog(R);
    q = [cos(th/2), sin(th/2)*w']';
end
