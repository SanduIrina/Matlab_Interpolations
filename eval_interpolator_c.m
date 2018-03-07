function N = eval_interpolator_c(tip, eps)
	x0 = linspace(-pi,pi,1001);
	k = 1;
	a=b=inf;
	while(a>eps)
		N = 2^k;
		x = linspace(-pi,pi,N);
		k++;
		switch(tip)
		case 1
			y = LAGRANGE(x0,x,f(x));
			err = eroare(f(x0),y);
			if(err>b) 
				N = inf
				return
			endif
			a = abs(b - err);
			b = err;

		case 2
			y = NEWTON(x0,x,f(x));
			err = eroare(f(x0),y);
			if(err>b) 
				N = inf
				return
			endif
			a = abs(b - err);
			b = err;

		case 3
			y = SPLINE_LINIAR(x0,x,f(x));
			err = eroare(f(x0),y);
			if(err>b) 
				N = inf
				return
			endif
			a = abs(b - err);
			b = err;

		case 4
			y = SPLINE1(x0,x,f(x));
			err = eroare(f(x0),y);
			if(err>b) 
				N = inf
				return
			endif
			a = abs(b - err);
			b = err;

		case 5
			y = SPLINE2(x0,x,f(x));
			err = eroare(f(x0),y);
			if(err>b) 
				N = inf
				return
			endif
			a = abs(b - err);
			b = err;

		case 6
			y = FOURIER(x0,x,f(x),length(x)/2);
			err = eroare(f(x0),y);
			if(err>b) 
				N = inf
				return
			endif
			a = abs(b - err);
			b = err;
			
		endswitch
	endwhile
endfunction

function r = eroare(y,z)
	h = (2*pi)/1001;
	s = 0;
	for i = 1:1001
		s+=abs(y(i)-z(i))^2;
	endfor
	r = (h*s)^(1/2);
endfunction

function y = f(x)
	y = (exp(3.*cos(x)))./(2*pi*besseli(0,3));
	return
endfunction

function s = LAGRANGE(x,x0,fx)
	m = length(x);
	n = length(x0);
	for k = 1:m
		s(k) = 0;
		for j = 1:n
			p = 1;
			for i = 1:n
				if(i~=j)
					p = p*(x(k)-x0(i))/(x0(j)-x0(i));
				endif
			endfor
			s(k) += fx(j).*p;
		endfor
	endfor
endfunction

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
 
function y = f(x)
    y = exp(3 .* cos(x)) ./ (2*pi*besseli(0, 3));
endfunction

function yi = SPLINE_LINIAR(x0, x, fx)
	for k = 1:length(x0)
	    for i = 2 : length(x)
	        if x0(k) <= x(i)
	            j = i - 1;
	            break;
	        endif
	    endfor
	    i = j;
	    h = x(i + 1) - x(i);
	    yi(k) =fx(i) + (x0(k) - x(i))*(fx(i+1)-fx(i))/h;
	endfor	
endfunction

function yi = SPLINE1(x0, x, fx)
	dy = derive(x);
	for ii = 1:length(x0)
	    for i = 2 : length(x)
	        if x0(ii) <= x(i)
	            j = i - 1;
	            break;
	        endif
	    endfor
	    i = j;
	    h(i) = x(i + 1) - x(i);
	    t = (x0(ii) - x(j)) / (h(j));
	    yi(ii) = fx(j) * (1 - t).^3 + (3 * fx(j) + h(j) * dy(j)) * t * (1 - t).^2 + (3 * fx(j + 1) - h(j) * dy(j + 1)) * t.^2 * (1 - t) + fx(j + 1) * t.^3;
	endfor	
endfunction

function d = derive(x0)
	d = (-3.*sin(x0).*exp(3.*cos(x0)))./(2*pi*besseli(0,3));
	return;
endfunction

function yi = SPLINE2(xi, x, y)   
	for k = 1 : length(xi)
		for j = 2 : length(x)
        	if xi(k) <= x(j)
        	    i = j - 1;
        	    break;
        	end
    	end
    	a = y;
    
   		sup(1) = 0;
    	h(1) = x(2) - x(1);
    	for w = 2 : length(x) - 1
    	    h(w) = x(w + 1) - x(w);
    	    sup(w) = h(w);
    	end
    	mij(1) = 1;
    	for w = 2 : length(x) - 1
    	    mij(w) = 2 * (h(w) + h(w - 1));
    	end
    	mij(length(x)) = 1;
    	for w = 1 : length(x) - 2
    	    inf(w) = h(w);
    	end
    	inf(length(x) - 1) = 0;
    	vec(1) = 0;
    	for w = 2 : length(x) - 1
    	    vec(w) = 3 * (a(w + 1) - a(w)) / h(w) - 3 * (a(w) - a(w - 1)) / h(w - 1);
    	end
    	vec(length(x)) = 0;
    	c = Thomas(sup, mij, inf, vec);
        
    	b(i) = (a(i + 1) - a(i)) / h(i) - (h(i) / 3) * (2 * c(i) + c(i + 1)); 
    	d(i) = (c(i + 1) - c(i)) / (3*h(i));
    
    	yi(k) = a(i);
   		yi(k) = yi(k) + b(i) * (xi(k) - x(i));
  		yi(k) = yi(k) + c(i) * (xi(k) - x(i)) * (xi(k) - x(i));
 		yi(k) = yi(k) + d(i) * (xi(k) - x(i)) * (xi(k) - x(i)) * (xi(k) - x(i));
endfor
endfunction

function s = FOURIER(x0,x,fx,m)
	for l = 1:length(x0)
		a(1,1:m+1) = 0;
		b(1,1:m+1) = 0;
		s(l) = 0;
		for k = 1:m+1
			for j = 1:2*m
				a(k) += fx(j)*cos((k-1)*x(j));
				b(k) += fx(j)*sin((k-1)*x(j));
			endfor
			a(k) = a(k)/m;
			b(k) = b(k)/m;
		endfor
		b(1) = b(m+1) = 0;
		for k = 2:m
			s(l) += a(k)*cos((k-1)*x0(l)) + b(k)*sin((k-1)*x0(l));
		endfor
		s(l) += (a(1) + a(m+1)*cos((m+1)*x0(l)))/2;
	endfor
endfunction
