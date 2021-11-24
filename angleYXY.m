function q = angleYXY(R)
# Angle phi1, phi2, phi3 from Rotation matrix with YXY convention
    phi1 = atan2(R(1,2) , R(3,2));
    phi2 = atan2(sqrt(1-R(2,2)^2) , R(2,2));
    phi3 = atan2(R(2,1) , -R(2,3));
    q = [phi1, phi2, phi3]';
end
