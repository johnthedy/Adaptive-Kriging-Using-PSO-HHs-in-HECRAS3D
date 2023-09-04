function [x,y,type,range_x] = pearcdf(alpha3_G,alpha4_G)
a=10*alpha4_G-12*alpha3_G^2-18;
b=alpha3_G*(alpha4_G+3);
c=4*alpha4_G-3*alpha3_G^2;
d=2*alpha4_G-3*alpha3_G^2-6;
delta=b^2-4*c*d;
r0=-b/(2*d);
r1=(-b-sqrt(delta))/(2*d);
r2=(-b+sqrt(delta))/(2*d);
k=1;tol=1e-3;
    if b~=0  && abs(b)>tol
            if delta>0 && abs(delta)>tol && d<0 && abs(d)>tol
                %TYPE I
                range_x = [r2 r1]; type=1;
            end
            if delta>0 && abs(delta)>tol && (d==0||abs(d)<tol)
                %TYPE III
                d=0;
                if b<0
                    range_x = [-inf -b/c];
                else
                    range_x = [-b/c inf];
                end
                type=3;
            end
            if delta >0 && abs(delta)>tol && d>0 && abs(d)>tol
                %TYPE VI
                if b<0
                    range_x = [-inf r1];
                else
                    range_x = [r2 inf];
                end
                type=6;
            end
            if delta <0 && abs(delta)>tol
                %TYPE IV
                range_x = [-inf inf]; type=4;
            end
            if delta==0 || abs(delta)<tol
                %TYPE V
                %delta=0;
                if b<0
                    range_x = [-inf r0];
                else
                    range_x = [r0 inf];
                end
                type=5;
            end
    else %b==0
            b=0;
            if delta>0 && abs(delta)>tol
                %TYPE II
                range_x = [sqrt(delta)/(2*d) -sqrt(delta)/(2*d)];
                type=2;
            end
            if delta==0 || abs(delta)<tol
                %TYPE Normal
                delta=0;
                range_x = [-inf inf];
                type=0;
            end
            if delta<0 && abs(delta)>tol
                %TYPE VII
                range_x = [-inf inf];
                type=7;
            end
        
    end        
K = integral(@(x)fun(x,a,b,c,d,delta,k,r0,r1,r2,range_x,type),range_x(1),range_x(2),'RelTol',1e-3,'AbsTol',1e-3,'ArrayValued',true);
K=K^-1;

if range_x(1)==-inf && range_x(2)==inf
    x=-10:0.1:10;
elseif range_x(1)~=-inf && range_x(2)==inf
    x=range_x(1):0.1:10;
elseif range_x(1)==-inf && range_x(2)~=inf
    x=-10:0.1:range_x(2);
else
    x=range_x(1):0.1:range_x(2);
end
for i=1:length(x)
    y(i)=fun(x(i),a,b,c,d,delta,K,r0,r1,r2,range_x,type);
end
end

function y = fun(x,a,b,c,d,delta,k,r0,r1,r2,range_x,type)
    zs=x;
    tol=1e-3;
    
    if min(range_x(2),max(zs,range_x(1)))~=zs
        y=[0];
    else
        switch type
            case 0
                %TYPE Normal
                y=(1/sqrt(2*pi))*exp(-zs.^2/2);
            case 1
                %TYPE I
                y=k*(zs-r2).^((-a*r2-b)/sqrt(delta)).*(r1-zs).^((a*r1+b)/sqrt(delta));
            case 2
                %TYPE II
                y=k*(-c/d-zs.^2).^(-a/(2*d));
            case 3
                %TYPE III
                y=k*(c+b*zs).^((a*c-b^2)/b^2).*exp(-a*zs/b);
            case 4
                %TYPE IV
                %ins = (a*b-2*b*d)/(d*sqrt(-delta)) * atan((b+2*d*zs)/sqrt(-delta));
                y=k*(c+b*zs+d*zs.^2).^(-a/(2*d)).*exp(((a*b-2*b*d)/(d*sqrt(-delta)))*atan((b+2*d*zs)/sqrt(-delta)));
                %y=k*(c+b*zs+d*zs.^2).^(-a/(2*d)).*exp(ins);
            case 5
                %TYPE V
                y=k*abs(zs-r0).^(-a/d).*exp((a*r0+b)./(d*(zs-r0)));
            case 6
                %TYPE VI
                if b<0
                    y=k*(r1-zs).^((a*r1+b)/sqrt(delta)).*(r2-zs).^((-a*r2-b)/sqrt(delta));
                else
                    y=k*(zs-r1).^((a*r1+b)/sqrt(delta)).*(zs-r2).^((-a*r2-b)/sqrt(delta));
                end     
            case 7
                %TYPE VII
                y=k*(c/d+zs.^2).^(-a/(2*d)); 
        end
    end
    
end