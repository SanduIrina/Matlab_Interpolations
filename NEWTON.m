function y = NEWTON(a, xs, fxs)
    n = length(xs);
    dd = difDiv(xs, fxs);
    y = [];
    for k = 1 : length(a)
        y(k) = fxs(1);
        for i = 2 : n
            y(k) += (dd(i) * prod(a(k) - xs(1 : i-1)));
        endfor
    endfor
endfunction
 
function dd = difDiv(x, y)
    n = length(x);
    dd = y;
    for i = 2 : n
        dd(i : n) = (dd(i : n) .- dd(i-1 : n-1)) ./ (x(i : n) .- x(1 : n-i+1));
    endfor
endfunction
 
function y = f(a, x)
    y = exp(a .* cos(x)) ./ (2*pi*besseli(0, a));
endfunction
