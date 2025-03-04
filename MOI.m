function I = MOI(M,L,Rp,Rd)
delta=3*M/(L*pi*(Rp^2+Rp*Rd+Rd^2));
a1=9/(20*pi); b1=3/80;
x=Rd/Rp;
a2=(1+x+x^2+x^3+x^4)/((1+x+x^2)^2);
b2=(1+4*x+10*x^2+4*x^3+x^4)/((1+x+x^2)^2);
Ixx=(a1*a2*M^2)/(delta*L)+b1*b2*M*L^2;
Iyy=Ixx;
Izz=(2*a1*a2*M^2)/(delta*L);
I=[Ixx 0 0;0 Iyy 0;0 0 Izz];
end

