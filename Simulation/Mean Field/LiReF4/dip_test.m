function dip_test(n)

q = [0 0 0];
a = [5.175 0 0;      
     0 5.175 0;
     0 0 10.75];
for ii = 1:size(n,2)
    d{ii} = dipole_direct(q,n(ii),a);
end

for ii = 1:size(d,2)
    D{ii} = sum(sum(d{ii},3),4)/2;
end

dzz = zeros(size(n));
for ii = 1:size(dzz,2)
    dzz(ii) = D{ii}(3,3);
end

plot(n,dzz,'-o')
xlabel('Number of unit cell')
ylabel('Dipole-dipole sum')

return