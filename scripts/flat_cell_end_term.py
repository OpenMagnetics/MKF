# Derivation of the end term E(nu) of MKF's screened width-cell model (MARTINEZ, a flat turn resolved into width cells,
# WindingProximityEffectLossesMartinezModel::calculate_flat_cell_end_term).
#
# An isolated rectangular conductor w x t, zero net current, in a uniform field H0 along its wide faces (A_ext = -mu0 H0 y).
#   Exact: a fine PEEC solve -- uniform-current sub-rectangles coupled by their exact mean logarithmic distance, refined
#   to delta/4 near the ends and the faces.
#   Model: full-thickness cells, each the 1-D slab with H0 on both faces, w * 2 H0^2 (F - 2G) / (2 sigma delta).
# Solving two widths W1 < W2 and taking the per-width slope as the bulk, the remainder is the two ends' excess over the
# slab; per end, in units of t * 2 H0^2 / (2 sigma delta), that is E(nu), nu = t / delta.
#
# Usage: python3 flat_cell_end_term.py <unused width> nu1 nu2 ...   (env W1, W2: the two strip widths, m; REF: refine)
# The table in MKF came from W1 = 1.5 mm, W2 = 3 mm, t = 0.20879 mm, REF = 1 and 2 (agreeing within 3 %); nu = 8 is
# width-independent within 1 %. Beyond nu = 10 MKF uses 2.728 - 3.413 nu^(-1/3) + 1.655 nu^(-2/3): the perfect-conductor
# limit and the delta^(1/3) series of a right-angle corner, fitted to the table's tail.
import os
import sys

import numpy as np

mu0 = 4e-7 * np.pi
sigma = 5.9117e7
c = mu0 / (2 * np.pi)


def ad_rect(u,v):
    u2=u*u; v2=v*v; r2=u2+v2
    with np.errstate(all='ignore'):
        val=(u2*v2/4-u2*u2/24-v2*v2/24)*np.log(np.where(r2>0,r2,1))-25/24*u2*v2
        val=val+np.where(u!=0,u2*u*v/3*np.arctan(v/np.where(u!=0,u,1)),0)+np.where(v!=0,u*v2*v/3*np.arctan(u/np.where(v!=0,v,1)),0)
    return np.where(r2>0,val/2,0)

def ml_rr(tx0,tx1,ty0,ty1,sx0,sx1,sy0,sy1):
    s=0
    for (u,su) in [(sx1-tx0,1),(sx1-tx1,-1),(sx0-tx0,-1),(sx0-tx1,1)]:
        for (v,sv) in [(sy1-ty0,1),(sy1-ty1,-1),(sy0-ty0,-1),(sy0-ty1,1)]:
            s=s+su*sv*ad_rect(u,v)
    return s/((sx1-sx0)*(sy1-sy0)*(tx1-tx0)*(ty1-ty0))


t = 0.20879e-3
H0 = 1.0


def run(nu,w):
    d=t/nu; f=1/(np.pi*mu0*sigma*d*d); w_=2*np.pi*f
    h=min(d/4,t/8)/float(os.environ.get("REF","1"))
    ny=max(8,int(np.ceil(t/h))); 
    # y: graded toward both faces (tanh-ish): use cosine spacing
    ye=-t/2+t*(1-np.cos(np.linspace(0,np.pi,ny+1)))/2
    L=w/2; xe_end=np.arange(0,min(4*t+6*d,L)+1e-15,h)   # fine near the end
    xr=list(xe_end); 
    while xr[-1]<L-1e-12: xr.append(min(L,xr[-1]+max(h,(xr[-1]-xr[-2])*1.2)))
    xr=np.array(xr); xe=np.unique(np.concatenate([-L+xr, L-xr]))
    X0,Y0=np.meshgrid(xe[:-1],ye[:-1],indexing='ij'); X1,Y1=np.meshgrid(xe[1:],ye[1:],indexing='ij')
    x0,x1,y0,y1=X0.ravel(),X1.ravel(),Y0.ravel(),Y1.ravel(); n=len(x0); a=(x1-x0)*(y1-y0)
    Lm=np.empty((n,n))
    for j in range(n): Lm[:,j]=-c*ml_rr(x0,x1,y0,y1,x0[j],x1[j],y0[j],y1[j])
    yc=(y0+y1)/2; Aext=-mu0*H0*yc
    M=np.zeros((n+1,n+1),complex); M[:n,:n]=np.diag(1/(sigma*a))+1j*w_*Lm; M[:n,n]=-1; M[n,:n]=1
    r=np.zeros(n+1,complex); r[:n]=-1j*w_*Aext
    I=np.linalg.solve(M,r)[:n]; P=np.sum(abs(I)**2/(sigma*a))/2
    F=(np.sinh(2*nu)+np.sin(2*nu))/(np.cosh(2*nu)-np.cos(2*nu)); G=(np.sinh(nu)*np.cos(nu)+np.cosh(nu)*np.sin(nu))/(np.cosh(2*nu)-np.cos(2*nu))
    Pmodel=w/(2*sigma*d)*2*H0**2*(F-2*G)
    unit=t*2*H0**2/(2*sigma*d)
    return n,P,unit
W1,W2=float(os.environ.get("W1","1.5e-3")),float(os.environ.get("W2","3e-3"))
for nu in [float(v) for v in sys.argv[2:]]:
    n1,P1,u=run(nu,W1); n2,P2,_=run(nu,W2)
    bulk=(P2-P1)/(W2-W1); ex=(P1-bulk*W1)/2/u
    print(f"t/delta={nu:5.2f} n={n1},{n2}  end excess (two-width) = {ex:7.4f} t",flush=True)
