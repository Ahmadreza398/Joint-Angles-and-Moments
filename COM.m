function c = COM(x)
C=(1+2*x+3*x^2)/(4*(1+x+x^2));

if x<1
    c=C;
elseif x>1
    c=1-C;
else
    c=0.5;
end

end

