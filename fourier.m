function p = fourier(x, xj, fxj)

	a = [];
	b = [0];
	m = floor(length(xj) / 2);

	for k = 1 : m + 1
        a(k) = sum(fxj(1 : 2*m) .* cos((k-1)*xj(1 : 2*m))) / m;
    endfor

    for k = 2 : m
        b(k) = sum(fxj(1 : 2*m) .* sin((k-1)*xj(1 : 2*m))) / m;
    endfor

    p = (a(1) + a(m+1)*cos((m+1) * x)) / 2;
    for k = 2 : m
        p += (a(k) * cos((k - 1) * x) + b(k) * sin((k - 1) * x));
    endfor

endfunction