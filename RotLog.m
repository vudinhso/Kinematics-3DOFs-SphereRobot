function [u, th] = RotLog(R)
    # Output a unit vector u and angle th from the rotation matrix
    if NearZero(R-eye(3))
        th = 0;
        u  = zeros(3,1);
    else
        th   = acos((trace(R)-1)/2); 
        umat = (R-R')/(2*sin(th)); 
        u    = [umat(3,2), umat(1,3), umat(2,1)]';
    end
end
