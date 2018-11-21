function f = bbf3(x)
	[ps,N] = size(x);
	if N ~= 30 error('dimension must 30'); end
	load bbf3
	xn = (x - abs(linspace(30,-30,N))) + ones(1,30);
	f = sum((xn(1:end-1).^2 - xn(2:end)).^2) + 13 * pi + sum((xn(1:end-1)-1).^2);
end
