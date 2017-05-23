% Computes the brightness value for a facet
% using Hapke's scattering law
%
% Code by Matti Viikinkoski and Mikko Kaasalainen

function ss=hapke_bright(E,E0,mu,mu0,p,th)
% E and EO should be row vectors (1x3)
tth=tan(pi/180*th);
cth=1/sqrt(1+pi*tth);
cal=E*E0';
alpha=acos(cal);
%mu0=normal*E0';
%mu=normal*E';
%keyboard

%fi=acos((cal-mu*mu0)/(sin(ai)*sin(ae)));
[sh,mueiefi,mu0eiefi]=shadow(mu,mu0,tth,cal);
mu0new=cth*mu0eiefi;
munew=cth*mueiefi;
dnom=munew+mu0new;
sls=mu*mu0new/dnom;
fh=hapke(p,munew,mu0new,alpha);
ss=fh*sls*sh;
end

function fh=hapke(p,mu,mu0,alpha)
alb=p(1);
h=p(2);
S0=p(3);
g=p(4);
ta=tan(0.5*alpha);
ca=cos(alpha);
B0=S0*(1+g)^2/(alb*(1-g));
bdnom=1+ta/h;
B=B0/bdnom; %Opposition surge
fhgdnom=1+g^2+2*g*ca;
fhg=(1-g^2)/(fhgdnom^1.5); %Particle Phase function
sqalb=sqrt(1-alb);
dnom=1+2*mu*sqalb;
dnom0=1+2*mu0*sqalb;
%Multiple scattering functions
chmu=(1+2*mu)/dnom; %H(alb,mu)
chmu0=(1+2*mu0)/dnom0; %H(alb,mu0)
fm=chmu*chmu0-1;
fh=(1+B)*fhg+fm;
end


function [sh,mueiefi,mu0eiefi]=shadow(mu,mu0,tth,cal)
i=acos(mu0);
e=acos(mu);
%ae and ae cannot be zero
if abs(i)<1e-6
    i=1e-6;
end
if abs(e)<1e-6
    e=1e-6;
end
innerfi=(cal-mu*mu0)/(sin(i)*sin(e));
fi=acos(innerfi);
f=exp(-2*tan(0.5*fi));
E1e=exp(-2/(tth*tan(e)*pi));
E1i=exp(-2/(tth*tan(i)*pi));
E2e=exp(-1/(tth^2*tan(e)^2*pi));
E2i=exp(-1/(tth^2*tan(i)^2*pi));
if(i<e)
    top=sin(0.5*fi)^2*E2i;
    bot=2-E1e-(fi/pi)*E1i;
    mu0ei0pi=mu0+sin(i)*tth*E2i/(2-E1i);
    mue0e0=mu+sin(e)*tth*E2e/(2-E1e);
    mueiefi=mu+sin(e)*tth*(E2e-top)/bot;
    coef=mu0/mu0ei0pi;
    sh=mueiefi*coef/(mue0e0*(1-f*(1-coef)));
    mu0eiefi=mu0+sin(i)*tth*(cos(fi)*E2e+top)/bot;
else
    top=sin(0.5*fi)^2*E2e;
    bot=2-E1i-(fi/pi)*E1e;
    mu0ei00=mu0+sin(i)*tth*E2i/(2-E1i);
    mue0epi=mu+sin(e)*tth*E2e/(2-E1e);
    mueiefi=mu+sin(e)*tth*(cos(fi)*E2i+top)/bot;
    coef=mu0/mu0ei00;
    sh=mueiefi*coef/(mue0epi*(1-f*(1-coef)));
    mu0eiefi=mu0+sin(i)*tth*(E2i-top)/bot;
end
end