function f = bbf2(x)
	[ps,N] = size(x);
	if N ~= 30 error('dimension must 30'); end
	load bbf2
	f = 10 * N + 12 * pi + sum(((x-20*pi) * M).^2 - 10 * cos(2 * pi * ((x-20*pi) * M)));
end
