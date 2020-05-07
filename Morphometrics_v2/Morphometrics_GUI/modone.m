%MODONE    Modulus after division starting from 1.
function r = modone(x,y)
if y==0
    r = NaN;
else
    n = floor(x/y);
    r = x - n*y;
    r(r==0)=y;
end
