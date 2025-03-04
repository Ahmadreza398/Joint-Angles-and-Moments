function [i_seg j_seg k_seg] = local_axis(p1,p2,p3)
j_seg=(p2-p1)/norm(p2-p1);
v_seg=(0.5*(p1+p2)-p3)/norm(0.5*(p1+p2)-p3);
i_seg=cross(j_seg,v_seg);
k_seg=cross(i_seg,j_seg);
end

