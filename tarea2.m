
pkg load optim;

D=load("escazu.dat");

# Construct the design matrix with the original data
Xo=[ones(rows(D),1),D(:,1)];

# The outputs vector with the original data
Yo=D(:,4);

# The slope for the data normalization
minArea = min(Xo(:,2));
maxArea = max(Xo(:,2));
mx = 2/(maxArea-minArea);
bx = 1-mx*maxArea;

NM=[1 bx; 0 mx];
X=Xo*NM; # Normalized data to interval -1 to 1




#The data matrix is modified to be used as linear.

X=[X(:,1),X(:,2).*2,X(:,2).^ 2];


# Normalize also the output
minPrice=min(Yo);
maxPrice=max(Yo);
my = 2/(maxPrice-minPrice);
by = 1-my*maxPrice;

Y = my*Yo + by;

imy = 1/my;
iby = -by/my;

function res=J(theta,X,Y)
  D=(X*theta'-Y*ones(1,rows(theta)));
  res=0.5*sum(D.*D,1)';
endfunction;


function res=gradJ(theta,X,Y)
  res=(X'*(X*theta'-Y*ones(1,rows(theta))))';
endfunction;


function resn=gradJn(theta,X,Y)
  
  
  # Centered derivatives
  delta=1e-2;
  dx=delta*600;
  dy=delta*2;

  Jtx=(J(theta+ones(rows(theta),1)*[dx,0],X,Y)-J(theta-ones(rows(theta),1)*[dx,0],X,Y))/(2*dx);
  Jty=(J(theta+ones(rows(theta),1)*[0,dy],X,Y)-J(theta-ones(rows(theta),1)*[0,dy],X,Y))/(2*dy);
  
  resn=[Jtx,Jty];
endfunction




figure(1);
hold off;

axis([-1 3 -1 3]);
xlabel("theta_0");
ylabel("theta_1");

t2=-0.5; #Initial value of thetha2

# Learning rate
alpha = 0.04;


while(1)
  hold on;


  figure(1);

  # Wait for mouse click and get the point (t0,t1) in the plot coordinate sys.
  [t0,t1,buttons] = ginput(1);
  t=[t0,t1,t2];
  gt=gradJ(t,X,Y);

  # Clean the previous plot 
  hold off;
  
   # Show the clicked point
  
  plot([t0],[t1],"*r");
  axis([-1 3 -1 3]);
  hold on;

  xlabel("theta0");
  ylabel("theta1");


  
  
  # Perform the gradient descent
  ts=t; # sequence of t's
  
  
  for i=[1:200] # max 100 iterations
    tc = ts(rows(ts),:); # Current position 
    gn = gradJ(tc,X,Y);  # Gradient at current position
    tn = tc - alpha * gn;# Next position
    ts = [ts;tn];
    if (norm(tc-tn)<0.001) break; endif;
  endfor


  # Draw the trajectory
  plot(ts(:,1),ts(:,2),"k-");
  
  plot(ts(:,1),ts(:,2),"ob");
  
  
  
  figure(2)
  hold off;
  
 
  plot3(ts(:,1),ts(:,2),ts(:,3),"ob");
  hold on;
  plot3(t0,t1,t2,"*r");
  plot3(ts(:,1),ts(:,2),ts(:,3),"k-");
  xlabel('theta0');
  ylabel('theta1');
  zlabel('theta2');
  
  
  
  
  figure(3);
  hold off;
  plot(Xo(:,2),Yo,"*b");
  hold on;
  
  # The line back in the samples
  areas=linspace(min(Xo(:,2)),max(Xo(:,2)),20);
  areasn=areas*mx+bx;
  
  # We have to de-normalize the normalized estimation
  prices=(t0+t1*2*areasn+t2*pow2(areasn))*imy+iby;
  
  plot(areas,prices,'k',"linewidth",3);
  
  xlabel('Área');
  ylabel('Precio');
  # and now with the intermediate versions
  for (i=[3:rows(ts)])
    u0=ts(i,1);
    u1=ts(i,2);
    u2=ts(i,3);
    
    prices=(u0+u1*2*areasn+u2*pow2(areasn))*imy+iby;
    plot(areas,prices,'r',"linewidth",1);
  endfor;
  plot(areas,prices,'g',"linewidth",3);

  
  axis([minArea maxArea minPrice maxPrice]); 

  
  
  #Now the error as a function of the learning rate 
  
  e=zeros(200:5);
  alp=[0.001;0.005;0.01;0.045;0.09;0.095];
  
  
  for j=[1:6]
    ta=t;
    for i=[1:200] 
       # max 100 iterations
      tc = ta(rows(ta),:); # Current position 
      ga = gradJ(tc,X,Y);  # Gradient at current position
      tn = tc - alp(j) * ga;# Next position
      e(i,j)=J(tc,X,Y); #error at current position
      ta = [ta;tn]; 

    endfor
    
  endfor
  
  
  
  figure(4);
  hold off;
  plot(e(:,1),"linewidth",3);
  hold on;
  plot(e(:,2),"linewidth",3);
  plot(e(:,3),"linewidth",3);
  plot(e(:,4),"linewidth",3);
  plot(e(:,5),"linewidth",3);
  plot(e(:,6),"linewidth",3);
  legend('alpha=0.001','alpha=0.005','alpha=0.01','alpha=0.045','alpha=0.09','alpha=0.095');
  xlabel('Iteracion');
  ylabel('error');
  
  axis([0,200,0,e(1,1)+e(1,1)*0.21]);
  
 






 #Estocástico
  
  te=t;
  
  for k=[1:19] # max ~200 iterations
    for j=[1:11]
    
      tc = te(rows(te),:); # Current position 
      gn = gradJ(tc,X(j,:),Y(j,:));  # Gradient at current position
      tn = tc - alpha * gn;# Next position
      te = [te;tn];
      ea =J(tn,X,Y);
      
      if (norm(ea)<3) break; endif;
    endfor;
    if (norm(ea)<3) break; endif;
  endfor
  

  figure(5)
  hold off;
  
 
  plot3(te(:,1),te(:,2),te(:,3),"ob");
  hold on;
  plot3(t0,t1,t2,"*r");
  plot3(te(:,1),te(:,2),te(:,3),"k-");
  xlabel('theta0');
  ylabel('theta1');
  zlabel('theta2');


  figure(6);
  hold off;
  plot(Xo(:,2),Yo,"*b");
  hold on;
  
  # The line back in the samples
  areas=linspace(min(Xo(:,2)),max(Xo(:,2)),20);
  areasn=areas*mx+bx;
  # We have to de-normalize the normalized estimation
  prices=(t0+t1*2*areasn+t2*pow2(areasn))*imy+iby;
  
  plot(areas,prices,'k',"linewidth",3);
  
  xlabel('Área');
  ylabel('Precio');
  
  # and now with the intermediate versions
  for (i=[3:rows(te)])
    u0=te(i,1);
    u1=te(i,2);
    u2=te(i,3);
    
    prices=(u0+u1*2*areasn+u2*pow2(areasn))*imy+iby;
    plot(areas,prices,'r',"linewidth",1);
  endfor;
  plot(areas,prices,'g',"linewidth",3);

  
  axis([minArea maxArea minPrice maxPrice]);

  
   #Now the error as a function of the learning rate 
 
 e=zeros(200:5);
 alp=[0.001;0.005;0.01;0.045;0.09;0.095];

  for j=[1:6]
    ta=t;
    for i=[0:18] 
       # max 200 iterations
       
      for k =[1:11]
          tc = ta(rows(ta),:); # Current position 
          ga = gradJ(tc,X(k,:),Y(k,:));  # Gradient at current position
          tn = tc - alp(j) * ga;# Next position
          e(i*11+k,j)=J(tc,X,Y); #error at current position
          ta = [ta;tn]; 
          
      endfor; 
      
    endfor
    
  endfor
   
  
  figure(7);
  hold off;
  plot(e(:,1),"linewidth",1);
  hold on;
  plot(e(:,2),"linewidth",1);
  plot(e(:,3),"linewidth",1);
  plot(e(:,4),"linewidth",1);
  plot(e(:,5),"linewidth",1);
  plot(e(:,6),"linewidth",1);
  legend('alpha=0.001','alpha=0.005','alpha=0.01','alpha=0.045','alpha=0.09','alpha=0.095');
  xlabel('Iteracion');
  ylabel('error');
  
  axis([0,200,0,e(1,1)+e(1,1)*0.21]);


 
endwhile;






