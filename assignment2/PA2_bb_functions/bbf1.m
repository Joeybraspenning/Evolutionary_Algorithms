function f = bbf1(x)
	[ps,N] = size(x);
	if N ~= 30 error('dimension must 30'); end
	load bbf1
	f = 11 * pi + sum(((x - linspace(-30,30,N)) * M).^2);
end
