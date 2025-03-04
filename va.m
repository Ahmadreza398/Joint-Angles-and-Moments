function [vel,acc] = va(t,r)
n=length(t);
dt=t(2)-t(1);
vel(1,:)=(r(2,:)-r(1,:))/dt;
vel(n,:)=(r(n,:)-r(n-1,:))/dt;
acc(1,:)=(vel(2,:)-vel(1,:))/dt;
acc(n,:)=(vel(n,:)-vel(n-1,:))/dt;
for i=2:n-1
   vel(i,:)=(r(i+1,:)-r(i-1,:))/(2*dt);
   acc(i,:)=(r(i+1,:)-2*r(i,:)+r(i-1,:))/(dt^2);
end
end

