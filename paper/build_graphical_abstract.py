
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import data, out
import matplotlib as mpl, matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Polygon, FancyBboxPatch, FancyArrowPatch, Circle

mpl.rcParams.update({"font.family":"Liberation Sans","font.size":12,"pdf.fonttype":42,
                     "svg.fonttype":"none","savefig.bbox":"standard","savefig.pad_inches":0})
BL="#1F5C8B"; OR_="#C2622A"; TE="#2E7D6E"; GY="#9AA5B1"; INK="#1F2933"; RD="#C1483A"
PU="#7B6AA0"; PB="#EDF2F7"; SOFT="#5A6572"
W_,H_=254.0,101.6
fig=plt.figure(figsize=(W_/25.4,H_/25.4))
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,W_); ax.set_ylim(0,H_); ax.axis("off")
lab=lambda x,y,s,sz=12,c=INK,w="normal": ax.text(x,y,s,fontsize=sz,color=c,fontweight=w,ha="center",va="bottom")
def arrow(x0,y0,x1,y1,lw=1.1,col="#4A5561"):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle="-|>",mutation_scale=9,lw=lw,
                                 color=col,shrinkA=0,shrinkB=0,zorder=3))
def bracket(x,y0,y1,stubs,side=1,col="#8A939E"):
    ax.plot([x,x],[y0,y1],lw=1.0,color=col,zorder=2)
    for y in stubs: ax.plot([x,x-3.0*side],[y,y],lw=1.0,color=col,zorder=2)

# ── inputs ────────────────────────────────────────────────────────────────
lab(34,88,"Somatic 3\u2032 UTR mutations")
for i in range(4):
    y=82.0-i*3.4
    ax.plot([13,55],[y]*2,lw=2.0,color=BL,solid_capstyle="butt")
    ax.add_patch(Polygon([[22+i*8.4,y+2.5],[26+i*8.4,y+2.5],[24+i*8.4,y+0.7]],closed=True,color=RD,zorder=3))
lab(34,52,"Transcript isoforms")
for i,wid in enumerate((42,34,26,18)):
    y=45.5-i*5.6
    ax.plot([13,13+wid],[y]*2,lw=.7,color=GY)
    for ex in range(3): ax.add_patch(Rectangle((13+ex*6.6,y-1.5),4.6,3.0,facecolor=BL,edgecolor="none"))
    ax.add_patch(Rectangle((13+19.8,y-1.5),wid-19.8,3.0,facecolor=PU,edgecolor="none"))
bracket(62,34,78,[76.0,38.0],side=1)
arrow(62,56,77,56)

# ── engine ────────────────────────────────────────────────────────────────
ax.add_patch(FancyBboxPatch((79,22),66,64,boxstyle="round,pad=0,rounding_size=4",
                            facecolor="white",edgecolor=BL,lw=1.4,zorder=1))
lab(112,78,"TxMod",16,BL,"bold")
ax.add_patch(Rectangle((85,29),54,45,facecolor=PB,edgecolor="none",zorder=1))
MX=112.0
ax.plot([MX,MX],[33,70],ls=(0,(2.0,1.6)),lw=.8,color=RD,zorder=2)
for i,(l_,r_) in enumerate(((23,20),(15,25),(6,10))):
    y=66.0-i*8.4
    ax.add_patch(FancyBboxPatch((MX-l_,y),l_+r_,3.8,boxstyle="round,pad=0,rounding_size=1.0",
                                facecolor="#DCE4EC",edgecolor=BL,lw=.7,zorder=2))
    ax.add_patch(Polygon([[MX-1.1,y+4.4],[MX+1.1,y+4.4],[MX,y+1.0]],closed=True,color=RD,zorder=4))
xs=np.linspace(87,137,220)
tr=(np.exp(-((xs-104)/3.0)**2)*0.72+np.exp(-((xs-121)/2.2)**2)*1.0+np.exp(-((xs-131)/2.6)**2)*0.34)
ax.plot(xs,31.5+tr*5.0,lw=1.2,color=TE,zorder=2)
lab(112,24.0,"m\u2076A probability per base",12,SOFT)

# ── outputs ───────────────────────────────────────────────────────────────
arrow(145,54,153,54)
OX=166.0; YS=[92.0,71.0,50.0,29.0]
bracket(156,18,86,[y-6.5 for y in YS],side=-1)
for k,(y,title) in enumerate(zip(YS,("Isoform-resolved m\u2076A change","Position inside the 3\u2032 UTR",
                                     "Support in measured m\u2076A","No signature of selection"))):
    lab(OX+24,y,title)
    arrow(156,y-6.5,OX-2,y-6.5)
    gy=y-10.8
    if k==0:
        for j,(wid,col,up) in enumerate(((28,OR_,True),(19,GY,False))):
            yy=gy+5.0-j*5.6
            ax.add_patch(FancyBboxPatch((OX,yy),wid,3.4,boxstyle="round,pad=0,rounding_size=0.9",
                                        facecolor="#DCE4EC",edgecolor=BL,lw=.7))
            xm=OX+wid*0.55
            ax.add_patch(Polygon([[xm-1.2,yy+4.2],[xm+1.2,yy+4.2],[xm,yy+0.9]],closed=True,color=RD,zorder=3))
            if up: ax.add_patch(FancyArrowPatch((OX+wid+2.5,yy+0.6),(OX+wid+8.5,yy+4.4),arrowstyle="-|>",
                                                mutation_scale=9,lw=1.5,color=col,shrinkA=0,shrinkB=0))
            else:  ax.plot([OX+wid+2.5,OX+wid+8.5],[yy+1.7]*2,lw=1.5,color=col,solid_capstyle="butt")
    elif k==1:
        xx=np.linspace(0,1,140); yy=gy+1.0+7.6*np.exp(-((xx-0.5)/0.24)**2)
        ax.plot(OX+xx*48,yy,lw=1.7,color=TE); ax.plot([OX,OX+48],[gy+0.8]*2,lw=.7,color=GY)
    elif k==2:
        for j,(h,col) in enumerate(((8.6,BL),(0.7,GY))):
            ax.add_patch(Rectangle((OX+7+j*24,gy+0.8),13,h,facecolor=col,edgecolor="none"))
        ax.plot([OX,OX+48],[gy+0.8]*2,lw=.7,color=GY)
    else:
        ax.plot([OX,OX+48],[gy+4.4]*2,lw=.8,color=GY)
        ax.plot([OX+24,OX+24],[gy+1.2,gy+7.6],ls=(0,(1.8,1.5)),lw=.9,color=GY)
        for dx in (-13,-4,5,14):
            ax.add_patch(Circle((OX+24+dx,gy+4.4),1.35,facecolor=SOFT,edgecolor="white",lw=.5,zorder=3))
# Centre the composition vertically in the 5:2 frame without changing the scale.
ax.set_ylim(8,109.6)
fig.savefig(out("graphical_abstract.pdf")); fig.savefig(out("graphical_abstract.png"),dpi=600)
