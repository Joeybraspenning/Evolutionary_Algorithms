function f = bbf5(x)
	[ps,N] = size(x);
	if N ~= 30 error('dimension must 30'); end
	load bbf5
	f = 20 + 15 * pi - 20 * exp(-0.2 * sqrt(sum((x - linspace(60,-30,N)) .^ 2, 2) / N)) - exp(sum(cos(2 * pi * (x - linspace(60,-30,N))), 2) / N) + exp(1);
end
