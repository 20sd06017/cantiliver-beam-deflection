clear
clc
close all
% Term Project
% Deflection at any point in cantilever beam using Castigliano's theorem
% Gokul Singh Mehta
% 20SD6017
syms P x

L=input("Length of the beam (in meter): ");
E=input("Elastic modulus (in N/m^2): ");
I=input("Area moment of inertia (m^4): ");
def_loc=input("Deflection at distance (distance from fixed end, in meter): ");
num_conc=input("number of concentrated (Point) loads: ");
conc=zeros(num_conc,2);

% Point load input.
for i=1:num_conc
    fprintf("location of point load %d from fixed end (in meter): ",i)
            conc(i,2)=input(" ");
            fprintf("load at point %d (in Newton):",i)
            conc(i,1)=input(" ");
end
   num_udl=input("number of uniformaly distributed loads: ");
   udl=zeros(num_udl,4);
   %udl load input
for j=1:num_udl
    fprintf("starting location of UDL %d from fixed end: ",j)
    udl(j,3)=input(" ");
    fprintf("length of UDL %d (in meter): ",j)
    udl(j,2)=input(" ");
    udl(j,4)=udl(j,3)+udl(j,2);
    fprintf("value of load of UDL %d (in Newton/meter): ",j)
    udl(j,1)=input(" ");
end
if num_udl==0
        j=0;
 end
num_uvl=input("number of uniformly varying loads: ");
uvl=zeros(num_uvl,5);
%uvl load input
for k=1:num_uvl
    fprintf("starting location of UVL %d from fixed end (in meter): ", k)
    uvl(k,4)=input(" ");
    fprintf("length of UVL %d (in meter): ",k) 
    uvl(k,3)=input(" ");
    uvl(k,5)=uvl(k,4)+uvl(k,3);
    fprintf("Value of UVL %d at one end (in N/m): ",k)
    uvl(k,1)=input(" ");
    fprintf("value of UVL %d at other end (in N/m): ",k)
    uvl(k,2)=input(" ");
   % If UVL and UDL combine act then they are divided into
   % UDL and UVL saperately.
   if  uvl(k,1)*uvl(k,2)>0
         j=j+1;
       if abs(uvl(k,1))>abs(uvl(k,2))
                udl(j,:)=[uvl(k,2), uvl(k,3), uvl(k,4), uvl(k,5)];
                uvl(k,1)=uvl(k,1)-uvl(k,2);
                uvl(k,2)=0;
       end
       if abs(uvl(k,2))>abs(uvl(k,1))
                udl(j,:)=[uvl(k,1), uvl(k,3), uvl(k,4), uvl(k,5)];
                uvl(k,2)=uvl(k,2)-uvl(k,1);
                uvl(k,1)=0;
       end
    end
end

num_moment=input("number of moments applied: "); % inward plane positive and vice versa.
moment=zeros(num_moment,2);
M_m=0;
%Moment input
for m=1:num_moment
    
       fprintf("location of moment %d from fixed end (in meter): ",m)
       moment(m,2)=input(" ");
       fprintf("value of moment %d (in Newton-meter): ",m);
       moment(m,1)=input(" ");
       M_m=M_m+moment(m,1);
end

udl_size=size(udl,1);
% If deflection is to be calculated at a point in such that there acts an 
% UDL or UVL, Then it UDL is divided into 2 UDLs or UVLs.
for j=1:udl_size
      r=udl_size;
      if (udl(j,3)<def_loc)&&(def_loc<udl(j,4))
                 udl(j,4)=def_loc;
                 udl(j,2)=def_loc-udl(j,3);
                 udl(r,:)=[udl(j,1),udl(j,2),udl(j,3),udl(j,4)];
                 r=r+1;
      end
end
udl_size=size(udl,1);
uvl_size=size(uvl,1);
for k=1:uvl_size
       p=uvl_size;
       q=udl_size;
    if (uvl(k,4)<def_loc)&&(def_loc<uvl(k,5))
        if uvl(k,1)==0
                 p=uvl_size;
                 w=uvl(k,2)*(def_loc-uvl(k,4))/(uvl(k,3));
                 udl(q+1,:)=[w, uvl(k,5)-def_loc, def_loc, uvl(k,5)];
                 uvl(k,:)=[uvl(k,1), w, def_loc-uvl(k,4), uvl(k,4), def_loc];
                 uvl(p+1,:)= [uvl(k,1), uvl(k,2)-w, uvl(k,5)-def_loc, def_loc, uvl(k,5)];
                 p=p+1;
                 q=q+1;
        end
        if(uvl(k,2)==0)
                 w=uvl(k,1)*(uvl(k,5)-def_loc)/uvl(k,3);
                 udl(q+1,:)=[w, def_loc-uvl(k,4), uvl(k,4),def_loc];
                 uvl(k,:)=[uvl(k,1)-w, uvl(k,2), def_loc-uvl(k,4), uvl(k,4),def_loc ];
                 uvl(p+1,:)=[ w, 0, uvl(k,5)-def_loc, def_loc, uvl(k,5)];
                 p=p+1;
                 q=q+1;
        end
    end
end
n=1;
for k=1:size(uvl(k,1))
    % If UVL load 0 is towards fixed end then it is flipped and divided
    % into UDL+UVL.
    if (uvl(k,1)==0)
                udlv(n,:)=[uvl(k,2), uvl(k,3), uvl(k,4), uvl(k,5)];
                uvlv(n,:)=[-uvl(k,2), 0, uvl(k,3), uvl(k,4), uvl(k,5)];
                uvl(k,:)=[];
                udl=[udl;udlv(n,:)];
                uvl=[uvl;uvlv(n,:)];
                n=n+1;
    end
end
flag=0;
Q=0;
for i=1:num_conc
    % If deflection calculation point coincides with one of the point
    % loads.
    if conc(i,2)==def_loc
        flag=1;
    end
    if flag==1
        Q=conc(i,1);
        conc(i,1)=0;
    end
end
d=L-def_loc;
Mp=P*(x-d);
M_conc=0;
for i=1:num_conc
    g=L-conc(i,2);
    t_conc=conc(i,1)*(x-g);
    M_conc=M_conc+t_conc;
end
udl_size_new=size(udl,1);
M_udl=0;
M_udl_0=0;
%Total moment due to UDLs
for j=1:udl_size_new
            if (udl(j,3))==0
               t_udl_0=udl(j,1)*(x-(L-udl(j,4)))^2*0.5;
               M_udl_0=M_udl_0+t_udl_0;
            end
            if (udl(j,3)~=0)
                a=udl(j,2);
                c=L-udl(j,4)+a/2;
                t_udl=udl(j,1)*a*(x-c);
                M_udl=M_udl+t_udl;     
            end   
end
M_udl=M_udl_0+M_udl;
uvl_size_new=size(uvl,1);
M_uvl_R=0;
M_uvl_R_0=0;
% Total moment due to UVL
for k=1:uvl_size_new
   if (uvl(k,4)==0)
        if (uvl(k,2)==0)
            t_uvl_R_0=uvl(k,1)*(x-uvl(k,5))^3/(6*uvl(k,3));
            M_uvl_R_0=M_uvl_R_0+t_uvl_R_0;
        end
    end
    if (uvl(k,4)~=0)
        if (uvl(k,2)==0)
            e=uvl(k,5)+(2/3)*(uvl(k,3));
            t_uvl_R=uvl(k,1)*uvl(k,3)*0.5*(x-e);
            M_uvl_R=M_uvl_R+t_uvl_R;
        end
    end   
end
    M_uvl=M_uvl_R_0+M_uvl_R; 
    def=(M_m+M_conc+M_udl+M_uvl+Mp);               % Total moment
    def_int=int(def^2, x, L-def_loc, L);           % Integrate Mx^2 wrt x 
    def_int_diff=diff(def_int, P);                 % differentiation wrt P
    P=Q;
    def_sub=subs(def_int_diff);
    deflection=(double(def_sub))*1000/(2*E*I);
    fprintf("deflection at distance %d meter from fixed end is (in mm): ",def_loc)
    disp(deflection)