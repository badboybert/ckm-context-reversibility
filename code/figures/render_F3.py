"""render_F3.py — Figure 3 (HEADLINE): no mark-intrinsic property predicts reversibility
(determinant null) + the constructive predictive null. 6 panels (a-f). Renders in matplotlib with the
CKM theme mirror (R absent). Reads only on-disk verified results; writes PNG+PDF + per-panel source CSV.
"""
import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
import theme_ckm as T
SIG=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # signature-pivot/
R=lambda *a: os.path.join(SIG,'results',*a)
OUT=os.path.dirname(os.path.abspath(__file__)); SD=os.path.join(OUT,'source_data'); os.makedirs(SD,exist_ok=True)

dm=pd.read_csv(R('determinant_meta.tsv'),sep='\t')
dt=dm[dm.layer=='transcriptome'].copy(); dp=dm[dm.layer=='proteome'].copy()
pn=pd.read_csv(R('predictive_null.tsv'),sep='\t')
ab=pd.read_csv(R('predictive_null_ablation.tsv'),sep='\t')
tau=pd.read_csv(R('orphan','tau_dissociation.tsv'),sep='\t')
NAME={'loeuf':'genetic constraint (LOEUF)','n_drug_log':'druggability','n_gwas_log':'GWAS burden',
 'tau':'tissue-specificity (τ)','arch_nsig':'# cis-eQTL signals','arch_str':'cis instrument F',
 'has_arch':'has cis-QTL','causal_nonEGFR':'causal status','is_enzyme':'enzyme','is_membrane':'membrane',
 'is_secreted':'secreted','loeuf_miss':'LOEUF missing','tau_miss':'τ missing'}
# export source data
for nm,df in [('F3a_determinant_transcriptome',dt),('F3c_determinant_proteome',dp),
              ('F3d_predictive_null',pn),('F3e_ablation',ab),('F3f_tau',tau)]:
    df.to_csv(os.path.join(SD,nm+'.csv'),index=False)

fig,axs=plt.subplots(2,3,figsize=T.figsize(183,168))
(a,b,c),(d,e,f)=axs
GRN=T.LAYER['transcriptome']; BLU=T.LAYER['proteome']; NEU=T.LAYER['neutral']

# --- a) determinant forest (transcriptome k=9) ---
da=dt.sort_values('pooled_beta').reset_index(drop=True)
y=np.arange(len(da)); mde=da.mde.median()
a.axvspan(-mde,mde,color=T.SEM['equiv_band'],zorder=0)
a.axvline(0,color='k',lw=0.6,zorder=1)
a.errorbar(da.pooled_beta,y,xerr=da.ci_halfwidth,fmt='none',ecolor=NEU,elinewidth=0.8,zorder=2,capsize=1.5)
a.scatter(da.pooled_beta,y,s=14,c=GRN,zorder=3,edgecolor='k',linewidth=0.3)  # uniform: no feature is a determinant
a.set_yticks(y); a.set_yticklabels([NAME.get(x,x) for x in da.feature],fontsize=7)
a.set_xlabel('pooled standardized β'); a.set_xlim(-0.075,0.075)
a.set_title('Determinant model (k=9, 4 tissues)',fontsize=T.BODY)
a.text(0.97,0.10,'all |β| < 0.029',transform=a.transAxes,ha='right',va='bottom',fontsize=7,style='italic')
a.text(0.97,0.03,f'grey band = ±MDE ({mde:.3f})',transform=a.transAxes,ha='right',va='bottom',fontsize=6,color='#666666')
T.panel_label(a,'a',dx=-0.55)

# --- b) causal-status null (focus) ---
cs=dt[dt.feature=='causal_nonEGFR'].iloc[0]
b.axvline(0,color='k',lw=0.6); b.axvspan(-mde,mde,color=T.SEM['equiv_band'],zorder=0)
b.errorbar([cs.pooled_beta],[0],xerr=[cs.ci_halfwidth],fmt='o',color=GRN,ecolor=NEU,
           elinewidth=1.0,capsize=3,markersize=7,markeredgecolor='k',markeredgewidth=0.4)
b.set_yticks([]); b.set_xlim(-0.05,0.05); b.set_ylim(-1,1)
b.set_xlabel('pooled standardized β'); b.set_title('Causal status',fontsize=T.BODY)
_lo=cs.pooled_beta-cs.ci_halfwidth; _hi=cs.pooled_beta+cs.ci_halfwidth
_hitxt='0.000' if abs(_hi)<5e-4 else f'{_hi:+.3f}'
b.text(0.5,0.74,f'β = {cs.pooled_beta:+.3f}\n(95% CI {_lo:+.3f}, {_hitxt})\n|β| < MDE (0.013)',
       transform=b.transAxes,ha='center',va='center',fontsize=7.5)
b.text(0.5,0.18,'causally-nominated marks\nreverse no differently',transform=b.transAxes,
       ha='center',va='center',fontsize=7,style='italic')
T.panel_label(b,'b')

# --- c) layer power: transcriptome powered vs proteome inconclusive (CI halfwidth) ---
c.scatter(np.full(len(dt),0)+np.random.default_rng(0).normal(0,0.04,len(dt)),dt.ci_halfwidth,
          s=12,c=GRN,edgecolor='k',linewidth=0.3,label=f'transcriptome (k=9)')
c.scatter(np.full(len(dp),1)+np.random.default_rng(1).normal(0,0.04,len(dp)),dp.ci_halfwidth,
          s=12,c=BLU,edgecolor='k',linewidth=0.3,label=f'proteome (k=3)')
c.set_xticks([0,1]); c.set_xticklabels(['transcriptome\n(k=9)','proteome\n(k=3)'],fontsize=7.5)
c.set_ylabel('95% CI half-width (SD)'); c.set_yscale('log')
c.set_title('Layer power',fontsize=T.BODY)
c.text(0.5,0.93,'powered          inconclusive',transform=c.transAxes,ha='center',fontsize=7,style='italic')
T.panel_label(c,'c')

# --- d) predictive null: within vs LOCO Spearman per context ---
ctx=pn.context.tolist(); x=np.arange(len(ctx))
for i in x: d.plot([i,i],[pn.within_spearman[i],pn.loco_spearman[i]],color='#CCCCCC',lw=1,zorder=1)
d.scatter(x,pn.within_spearman,s=22,c=GRN,edgecolor='k',linewidth=0.3,zorder=3,label='within-context (5-fold CV)')
d.scatter(x,pn.loco_spearman,s=22,facecolor='white',edgecolor=GRN,linewidth=1.0,zorder=3,label='leave-one-context-out')
d.axhline(0,color='k',lw=0.6)
d.axhline(pn.loco_spearman.mean(),color=NEU,lw=0.8,ls=':')
d.set_xticks(x); d.set_xticklabels(ctx,rotation=30,ha='right',fontsize=7)
d.set_ylabel('Spearman (predicted vs actual)'); d.set_ylim(-0.08,0.28)
d.set_title('Predictive null',fontsize=T.BODY)
d.text(0.97,0.55,f'LOCO mean +{pn.loco_spearman.mean():.3f}',transform=d.transAxes,ha='right',fontsize=7,style='italic')
d.legend(fontsize=6.3,loc='lower right',frameon=True,framealpha=0.9,edgecolor='none',handletextpad=0.3)
T.panel_label(d,'d',dx=-0.28)

# --- e) ablation: binary LOCO AUC by feature subset ---
x=np.arange(len(ab)); w=0.26
e.bar(x-w,ab.loco_auc_full,w,color=NEU,label='all features',edgecolor='k',linewidth=0.3)
e.bar(x,ab.loco_auc_intrinsic_only,w,color=GRN,label='intrinsic only',edgecolor='k',linewidth=0.3)
e.bar(x+w,ab.loco_auc_measurability_only,w,color=BLU,label='measurability only',edgecolor='k',linewidth=0.3)
e.axhline(0.5,color='k',lw=0.6,ls='--')
e.set_xticks(x); e.set_xticklabels(ab.context,rotation=30,ha='right',fontsize=7)
e.set_ylabel('cross-context AUC'); e.set_ylim(0.45,0.85)
e.set_title('Detectability ablation (binary, secondary)',fontsize=8)
e.legend(fontsize=6.5,loc='upper left',frameon=False,handletextpad=0.3,ncol=1)
T.panel_label(e,'e',dx=-0.28)

# --- f) tau magnitude vs direction dissociation ---
lab=['magnitude','restoration','restoration\n(adj.)']
val=[tau[tau.col=='y_mag'].tau_beta.iloc[0],tau[tau.col=='y_dir'].tau_beta.iloc[0],tau[tau.col=='y_dir_ca'].tau_beta.iloc[0]]
cols=[T.SEM['up'],T.SEM['down'],'#9ecae1']
f.bar(np.arange(3),val,color=cols,edgecolor='k',linewidth=0.4,width=0.6)
f.axhline(0,color='k',lw=0.6)
f.set_xticks(np.arange(3)); f.set_xticklabels(lab,fontsize=7); f.set_xlim(-0.6,2.6); f.set_ylim(-0.075,0.125)
f.set_ylabel('τ → outcome (β)'); f.set_title('Tissue-specificity: magnitude, not direction',fontsize=7.5)
for i,v in enumerate(val):
    f.annotate(f'{v:+.3f}',(i,v),ha='center',va='bottom' if v>=0 else 'top',
               xytext=(0,3 if v>=0 else -3),textcoords='offset points',fontsize=7,clip_on=False)
f.text(0.5,-0.30,'magnitude = |Δ|; restoration = toward lean',transform=f.transAxes,ha='center',fontsize=6,style='italic')
T.panel_label(f,'f',dx=-0.28)

fig.tight_layout(w_pad=1.6,h_pad=2.2)
png,pdf=T.save_ckm(fig,os.path.join(OUT,'Figure_3'))
print('wrote',png,'and',pdf)
print('source data ->',SD)
