function plot_axis(C,i_seg,j_seg,k_seg,p1,p2,p3)

quiver3(C(1),C(2),C(3),0.1*i_seg(1),0.1*i_seg(2),0.1*i_seg(3));
quiver3(C(1),C(2),C(3),0.1*j_seg(1),0.1*j_seg(2),0.1*j_seg(3));
quiver3(C(1),C(2),C(3),0.1*k_seg(1),0.1*k_seg(2),0.1*k_seg(3));

h=patch('Vertices',[p1;p2;p3]);
set(h,'FaceColor','r','FaceAlpha',0.5);
axis equal; 

end

