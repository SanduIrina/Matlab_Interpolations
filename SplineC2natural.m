function yi = SplineC2natural(x, y, xi)
	for k = 1 : length(xi)    
		for j = 2 : length(x)
        	if xi(k) <= x(j)
        	    i = j - 1;
        	    break;
        	end
    	end
    	% construire c
    
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
