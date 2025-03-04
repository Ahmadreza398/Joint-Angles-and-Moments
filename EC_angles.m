function [a,b,g] = EC_angles(R)
    a=atan2d(-R(3,2),R(3,3));
    b=atan2d(R(3,1),norm([R(1,1) R(2,1)]));
    g=atan2d(-R(2,1),R(1,1));
end

