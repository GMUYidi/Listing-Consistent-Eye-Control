"""Symbolic and synthetic numerical checks, not experimental validation data."""
from pathlib import Path
import json
import sympy as s
import numpy as np

out = Path(__file__).resolve().parents[2] / 'outputs/theory'
out.mkdir(parents=True,exist_ok=True)
p, t, r = s.symbols('psi theta phi', real=True)
pd, td, rd = s.symbols('psi_dot theta_dot phi_dot', real=True)
pa, ta = s.symbols('psi_ddot theta_ddot', real=True)
I = s.symbols('I_0', positive=True)
tx, ty, tz = s.symbols('tau_phi tau_theta tau_psi', real=True)
v1, v2 = s.symbols('v_1 v_2', real=True)
cp, ct, sp, st = s.cos(p), s.cos(t), s.sin(p), s.sin(t)
D, a, b = 1+cp*ct, cp+ct, sp*st
T = s.Matrix([[0,-sp,ct*cp],[0,cp,ct*sp],[1,0,-st]])
qd = s.Matrix([pd,td,rd])
Tdot = T.diff(p)*pd + T.diff(t)*td
torque = s.Matrix([tx,ty,tz])
checks = {}
def simp(e): return s.trigsimp(s.cancel(s.trigsimp(e)))
def assert_zero(name, expr):
    values = list(expr) if isinstance(expr,s.MatrixBase) else [expr]
    simplified = [simp(e) for e in values]
    checks[name] = all(e == 0 for e in simplified)
    if not checks[name]: raise AssertionError((name,simplified))

Rz=s.Matrix([[cp,-sp,0],[sp,cp,0],[0,0,1]])
Ry=s.Matrix([[ct,0,st],[0,1,0],[-st,0,ct]])
Rx=s.Matrix([[1,0,0],[0,s.cos(r),-s.sin(r)],[0,s.sin(r),s.cos(r)]])
R=Rz*Ry*Rx
Rdot=R.diff(p)*pd+R.diff(t)*td+R.diff(r)*rd
W=Rdot*R.T
assert_zero('T_is_space_frame_rotation_derivative',s.Matrix([W[2,1],W[0,2],W[1,0]])-T*qd)

full=s.Matrix([
 (tx*cp*s.tan(t)+ty*sp*s.tan(t)+tz)/I+(td*rd+st*td*pd)/ct,
 (ty*cp-tx*sp)/I-ct*pd*rd,
 (tx*cp+ty*sp)/(I*ct)+(st*td*rd+pd*td)/ct])
assert_zero('six_state_dynamics_match_T',T*full+Tdot*qd-torque/I)

# The principal Listing branch has L_psi=sin(theta)/D and L_theta=sin(psi)/D.
L=2*s.atan(s.tan(p/2)*s.tan(t/2))
assert_zero('Listing_first_derivative_psi',s.diff(L,p)-st/D)
assert_zero('Listing_first_derivative_theta',s.diff(L,t)-sp/D)
Ldot=(pd*st+td*sp)/D
Ldd=s.diff(Ldot,p)*pd+s.diff(Ldot,t)*td+s.diff(Ldot,pd)*pa+s.diff(Ldot,td)*ta
Ldd_claim=(pa*st+ta*sp)/D+(2*pd*td*(ct+cp)+sp*st*(pd**2*ct+td**2*cp))/D**2
assert_zero('Listing_second_derivative',Ldd-Ldd_claim)
full_L=full.subs(rd,Ldot)
Ldd_dynamic=Ldd.subs({pa:full_L[0],ta:full_L[1]})
residual_expected=(-D*tx-sp*ct*ty+st*tz)/(I*D)
assert_zero('acceleration_compatibility_residual',Ldd_dynamic-full_L[2]-residual_expected)
tx_constraint=(-ty*sp*ct+tz*st)/D
assert_zero('torque_compatibility',residual_expected.subs(tx,tx_constraint))

F=td**2*sp+td*pd*st*(2+cp*ct)
V=pd*td*sp+pd**2*st
f=s.Matrix([F/(D*ct),-ct*V/D])
G=s.Matrix([[a/ct,b/ct],[-b,a]])/(I*D)
red=s.Matrix([full_L[0],full_L[1]]).subs(tx,tx_constraint)
assert_zero('reduced_f_and_G',red-f-G*s.Matrix([tz,ty]))
assert_zero('determinant_identity',a**2+b**2-D**2)
assert_zero('G_determinant',G.det()-1/(I**2*ct))
Ginv=I*s.Matrix([[a*ct,-b],[b*ct,a]])/D
assert_zero('G_inverse',G*Ginv-s.eye(2))
u=Ginv*(s.Matrix([v1,v2])-f)
u_expanded=s.Matrix([
 I*(a*ct*v1/D-b*v2/D-a*F/D**2-b*ct*V/D**2),
 I*(b*ct*v1/D+a*v2/D-b*F/D**2+a*ct*V/D**2)])
assert_zero('expanded_inverse_torques',u-u_expanded)
assert_zero('closed_loop_accelerations',f+G*u-s.Matrix([v1,v2]))
Q=T.T*torque
assert_zero('Euler_generalized_torque',Q-s.Matrix([tz,-sp*tx+cp*ty,ct*cp*tx+ct*sp*ty-st*tz]))

rng=np.random.default_rng(20261006)
num=s.lambdify((p,t,pd,td,I,v1,v2),[G,Ginv,f,u,red.subs({tz:u[0],ty:u[1]},simultaneous=True)],'numpy',cse=True)
max_inverse_error=0.; max_acceleration_error=0.; n=200
for _ in range(n):
 vals=(rng.uniform(-np.pi/4,np.pi/4),rng.uniform(-np.pi/6,np.pi/6),rng.uniform(-10,10),rng.uniform(-10,10),6.75e-7,rng.uniform(-10000,10000),rng.uniform(-10000,10000))
 ng,ni,nf,nu,nred=num(*vals)
 max_inverse_error=max(max_inverse_error,float(np.max(np.abs(ng@ni-np.eye(2)))))
 max_acceleration_error=max(max_acceleration_error,float(np.max(np.abs(nred.ravel()-np.array(vals[-2:])))))
assert max_inverse_error<1e-12 and max_acceleration_error<1e-8
report={'purpose':'Symbolic identity checks and synthetic numerical spot checks only; no experimental data generated.',
        'sympy_version':s.__version__,'checks':checks,'synthetic_samples':n,
        'maximum_G_inverse_residual':max_inverse_error,
        'maximum_closed_loop_acceleration_residual_rad_s2':max_acceleration_error,
        'compatibility_residual':str(residual_expected),
        'G_determinant':str(1/(I**2*ct)),
        'G_inverse_latex':s.latex(Ginv),'expanded_torques_latex':s.latex(u_expanded)}
(out/'math_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
