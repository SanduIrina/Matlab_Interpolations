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
