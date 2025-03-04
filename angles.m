function [alpha,beta,gama] = angles(R)
    a=atan2d(-R(3,2),R(3,3));
    b=atan2d(R(3,1),norm([R(1,1) R(2,1)]));
    g=atan2d(-R(2,1),R(1,1));
    Rx=[1 0 0;0 cosd(a) sind(a);0 -sind(a) cosd(a)];
    Ry=[cosd(b) 0 -sind(b);0 1 0;sind(b) 0 cosd(b)];
    ang=[a;0;0]+Rx'*[0;b;0]+Rx'*Ry'*[0;0;g];
    alpha=ang(1);
    beta=ang(2);
    gama=ang(3);
end

