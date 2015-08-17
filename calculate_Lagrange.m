%Code Used to Create Graphs and Calculate Cost Function Limits

syms x y z
f=4-x^2-y^2+x^4+y^4
g=exp(x^2)+2*exp(y^2)  
gradf=jacobian(f,[x,y,z])

[xfcrit,yfcrit]=solve(gradf(1),gradf(2));
[xfcrit,yfcrit]  

gfun=inline(vectorize(g))
double(gfun(xfcrit,yfcrit))  
ffun=inline(vectorize(f))
ffun(xfcrit([1,2,3]),yfcrit([1,2,3]))  

gradg=jacobian(g,[x,y,z])
gradcross=cross(gradf,gradg);
lagr=gradcross(3)  
[xboundcrit,yboundcrit]=solve(g-4,lagr)  

genplot(g, -2:.1:2,-2:.1:2,'contour',[4,4],'r'); hold on;
genplot(lagr, -2:.05:2,-2:.05:2,'contour',[0,0],'b'); hold off;  

xaxiscrits=solve(subs(g-4,y,0));
[xaxiscrits,[0;0]]

yaxiscrits=solve(subs(g-4,x,0));
[[0;0],yaxiscrits]  

double([xaxiscrits,[0;0]])
double([[0;0],yaxiscrits]) 

double(ffun(xaxiscrits,[0;0]))
double(ffun([0;0],yaxiscrits))  

[xb1,yb1]=newton2d(g-4,lagr,.5,.5)  
[xb2,yb2]=newton2d(g-4,lagr,-.5,.5)  
[xb3,yb3]=newton2d(g-4,lagr,-.5,-.5)  
[xb4,yb4]=newton2d(g-4,lagr,.5,-.5)  
ffun([xb1,xb2,xb3,xb4],[yb1,yb2,yb3,yb4])  

genplot(g, -2:.1:2,-2:.1:2,'contour',[4,4],'k'); hold on;
genplot(f, -2:.05:2,-2:.05:2,'contour',[3:.1:4.2]); hold off;
colorbar 
