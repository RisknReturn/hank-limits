"""C03 follow-up diagnostic, 2026-10-22.

Run with sequence-jacobian==1.0.0, numpy==2.5.3, scipy==1.18.1,
numba==0.68.0. Pass the extracted repository directory as the only argument.
This uses the repository without editing it. The differential controls patch
the loaded SSJ class in this diagnostic process only; restoration occurs at exit.
The original 30-run mutation outputs are not changed.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[k]='1'
import sys, pathlib, json, platform
from importlib import metadata
if len(sys.argv) != 2:
    raise SystemExit('Usage: python C03_probe.py PATH_TO_HANK_LIMITS_REPOSITORY')
REPO=pathlib.Path(sys.argv[1]).resolve()
if not (REPO/'hank_limits/models.py').is_file():
    raise SystemExit('Repository must contain hank_limits/models.py')
OUT=pathlib.Path(__file__).resolve().parent.parent/'out'
sys.path.insert(0,str(REPO))
import numpy as np
from hank_limits import models as M
from sequence_jacobian.utilities.function import DifferentiableCombinedExtendedFunction as DCEF

original_diff=DCEF.diff
def centered_diff(self, shock_dict, h=None, outputs=None, hide_zeros=False, twosided=None):
    if twosided is None: twosided=self.default_twosided
    return original_diff(self,shock_dict,h,outputs,hide_zeros,twosided)

def analytic_income_diff(self,shock_dict,h=None,outputs=None,hide_zeros=False,twosided=None):
    if self.name=='income_cyclical' and set(shock_dict)=={'Y'}:
        state=self.diff_functions['income_cyclical'].input_dict
        assert state['Y']==1 and state['Yss']==1, 'Analytic control assumes Y=Yss=1'
        return {'y':state['e_grid']*(1+state['zeta'])*shock_dict['Y']}
    return centered_diff(self,shock_dict,h,outputs,hide_zeros,twosided)

T=150
dr=-.01*.8**np.arange(T)
I=np.eye(T)
L=np.eye(T,k=-1)
rows=[]
for cyc,b in [(False,.1),(True,.1),(True,.01),(True,.0001)]:
    m=M.werning_model(b,cyc)
    ss=m['ss']; block=m['block']; prim=m['prim']
    e,zeta,pi=prim['e'],prim['zeta'],prim['pi']
    ypp=float(pi @ (e*zeta*(1+zeta)))
    for factor,mode in [(1.,'original'),(.5,'original'),(.25,'original'),(1.,'centered'),(.5,'centered'),(1.,'analytic')]:
        centered=mode!='original'
        DCEF.diff={'original':original_diff,'centered':centered_diff,'analytic':analytic_income_diff}[mode]
        h=1e-3*b*factor
        J=block.jacobian(ss,['Y','r'],['C','A'],T=T,twosided=True,h=h)
        CY,AY,Cr,Ar=[np.asarray(J[o][i]) for o,i in [('C','Y'),('A','Y'),('C','r'),('A','r')]]
        # Beginning-of-period assets are zero; allow the actual tiny SS calibration residual.
        RY=CY+AY-(1+ss['r'])*L@AY-I
        Rr=Cr+Ar-(1+ss['r'])*L@Ar-ss['A']*I
        da=np.linalg.solve(AY,-Ar@dr)
        dg=np.linalg.solve(CY-I,-Cr@dr)
        yy=lambda Y:e*Y**(1+zeta)
        fd_y=(yy(1+h)-yy(1))/h
        cd_y=(yy(1+h)-yy(1-h))/(2*h)
        fd_agg=float(pi@fd_y)
        cd_agg=float(pi@cd_y)
        # Directly measure the hetinput derivative chosen by SSJ.
        ht=block.hetinputs.differentiable(block.extract_ss_dict(ss),h,True)
        actual_agg=float(pi@ht.diff({'Y':1})['y'])
        # Feed the forward derivative's aggregate income into the budget identity.
        corrected=RY-(actual_agg-1)*I
        row=dict(cyclical=cyc,b=b,h=h,centered_hetinput=centered,mode=mode,
                 budget_Y_max=float(np.abs(RY).max()),budget_r_max=float(np.abs(Rr).max()),
                 residual_after_actual_income_derivative=float(np.abs(corrected).max()),
                 diag00=float(RY[0,0]),offdiag_max=float(np.abs(RY-np.diag(np.diag(RY))).max()),
                 aggregate_income_derivative=actual_agg,forward_derivative=fd_agg,
                 central_derivative=cd_agg,half_income_second_derivative=ypp/2,
                 dY_relative_difference=float(np.abs(dg-da).max()/np.abs(da).max()),
                 cond_goods=float(np.linalg.cond(CY-I)),cond_asset=float(np.linalg.cond(AY)),
                 cross_goods_on_asset=float(np.abs((CY-I)@da+Cr@dr).max()),
                 cross_asset_on_goods=float(np.abs(AY@dg+Ar@dr).max()),
                 asset_ss=float(ss['A']))
        rows.append(row)
        print(f"b={b:g} cyclical={cyc} mode={mode} h={h:g} RY={row['budget_Y_max']:.6g} dY_gap={row['dY_relative_difference']:.6g}",flush=True)
DCEF.diff=original_diff
output=dict(environment=dict(python=sys.version,platform=platform.platform(),packages={p:metadata.version(p) for p in ('sequence-jacobian','numpy','scipy','numba','llvmlite','pytest')}),T=T,rows=rows)
(OUT/'c03_probe_v0_1.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
