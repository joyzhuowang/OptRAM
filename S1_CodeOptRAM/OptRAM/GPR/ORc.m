function [ output ] = ORc( a , b )

output=[0 1];
output(1) = (a(1)*a(2)+b(1)*b(2))/(a(2)+b(2));
output(2) = a(2)+b(2);
end

