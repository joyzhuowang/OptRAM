function [ output ] = ANDc( a , b )
output = [0 1];
if a(1) < b(1);
    output(1) = a(1);
else
    output(1) = b(1);
end;
end
