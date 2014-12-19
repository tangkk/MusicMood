function y = myAutocorr(x, L)
y = zeros(1,L);
Lx = length(x);
for i = 1:1:L
    for j = i:1:Lx
        y(i) = y(i) + x(j)*x(j-i+1);
    end
end