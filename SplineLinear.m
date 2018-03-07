function yi = SplineLinear(xc, y, x)
	for ii = 1:length(x)
	    for i = 2 : length(xc)
	        if x(ii) <= xc(i)
	            j = i - 1;
	            break;
	        endif
	    endfor
	    i = j;
	    h = xc(i + 1) - xc(i);
	    yi(ii) =y(i) + (x(ii) - xc(i))*(y(i+1)-y(i))/h;
	endfor	
endfunction
