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
