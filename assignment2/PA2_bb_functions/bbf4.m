function f = bbf4(x)
	[ps,N] = size(x);
	if N ~= 30 error('dimension must 30'); end
	load bbf4
	x = (x - abs(linspace(30,-30,N))) * M;
	f = 1;
	for i=1:N
    	f = f .* cos(x(:,i) ./ sqrt(i));
	end
	f = sum(x.^2, 2) ./ 4000 - f + 1 + 14 * pi;
end
