%this short script is a collection of debug tools for 'Field_scan.m'

for ii = 1:floor(length(freq_temp)/nop)
A(ii,:)=freq_temp((ii-1)*23+1:ii*23);
end