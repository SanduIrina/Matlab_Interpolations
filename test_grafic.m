x = linspace(-pi,pi,8);
x0 = linspace(-pi,pi,1001);
y = (exp(3.*cos(x)))./(2*pi*besseli(0,3));
fileID = fopen('sunspot.dat','r');
formatSpec = '%d %f';
size = [2,300];
m = fscanf(fileID,formatSpec,size);
fclose(fileID);
cont = fix(linspace(1,300,32));
x1 = m(cont);
y1 = m(cont);

subplot(2,1,1)
plot(x0,LAGRANGE(x0,x,y),x0,NEWTON(x0,x,y),x0,SPLINE_LINIAR(x0,x,y),x0,SPLINE1(x0,x,y),x0,SPLINE2(x0,x,y),x0,FOURIER(x0,x,y,4));
subplot(2,1,2)
plot(x0,LAGRANGE(x0,x1,y1),x0,NEWTON(x0,x1,y1),x0,SPLINE_LINIAR(x0,x1,y1),x0,SPLINE1(x0,x1,y1),x0,SPLINE2(x0,x1,y1),x0,FOURIER(x0,x1,y1,16));
