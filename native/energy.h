/* Derived from ViennaRNA 2.7.2 src/ViennaRNA/eval/eval_internal.c.
 * Authors: ViennaRNA contributors, Institute for Theoretical Chemistry,
 * University of Vienna. See COPYING.ViennaRNA.
 * Modification: inline exact specialization for <= MAXLOOP total unpaired
 * bases, default salt and permitted GU closure. All other domains delegate
 * to the unchanged RNAlib evaluator. Integer operations retain upstream order.
 */
#ifndef RDO_ENERGY_H
#define RDO_ENERGY_H
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/params/default.h>
#include <ViennaRNA/eval/internal.h>

static inline int rdo_E_internal_bounded(unsigned n1,unsigned n2,unsigned type,unsigned type_2,
                                 int si1,int sj1,int sp1,int sq1,vrna_param_t *P){
  unsigned nl=n1>n2?n1:n2,ns=n1>n2?n2:n1;
  int energy=0;
  if(nl==0)return P->stack[type][type_2]+P->SaltStack;
  switch(ns){
    case 0:
      energy=P->bulge[nl];
      if(nl==1)energy+=P->stack[type][type_2];
      else{
        if(type>2)energy+=P->TerminalAU;
        if(type_2>2)energy+=P->TerminalAU;
      }
      break;
    case 1:
      if(nl==1)energy=P->int11[type][type_2][si1][sj1];
      else if(nl==2){
        if(n1==1)energy=P->int21[type][type_2][si1][sq1][sj1];
        else energy=P->int21[type_2][type][sq1][si1][sp1];
      }else{
        energy=P->internal_loop[nl+1];
        energy+=MIN2(MAX_NINIO,(int)(nl-ns)*P->ninio[2]);
        energy+=P->mismatch1nI[type][si1][sj1]+P->mismatch1nI[type_2][sq1][sp1];
      }
      break;
    case 2:
      if(nl==2){energy=P->int22[type][type_2][si1][sp1][sq1][sj1];break;}
      else if(nl==3){
        energy=P->internal_loop[5]+P->ninio[2];
        energy+=P->mismatch23I[type][si1][sj1]+P->mismatch23I[type_2][sq1][sp1];break;
      }
      /* fall through */
    default:
      energy+=P->internal_loop[nl+ns];
      energy+=MIN2(MAX_NINIO,(int)(nl-ns)*P->ninio[2]);
      energy+=P->mismatchI[type][si1][sj1]+P->mismatchI[type_2][sq1][sp1];
      break;
  }
  return energy;
}

static inline int rdo_E_internal(unsigned n1,unsigned n2,unsigned type,unsigned type_2,
                                 int si1,int sj1,int sp1,int sq1,vrna_param_t *P){
  if(!P || n1>MAXLOOP || n2>MAXLOOP || n1+n2>MAXLOOP ||
     P->model_details.salt!=VRNA_MODEL_DEFAULT_SALT || P->model_details.noGUclosure)
    return vrna_E_internal(n1,n2,type,type_2,si1,sj1,sp1,sq1,P);
  return rdo_E_internal_bounded(n1,n2,type,type_2,si1,sj1,sp1,sq1,P);
}


static inline int rdo_E_exterior_stem(unsigned type,int n5d,int n3d,vrna_param_t*p){
 int e=0;
 if(n5d>=0&&n3d>=0)e+=p->mismatchExt[type][n5d][n3d];
 else if(n5d>=0)e+=p->dangle5[type][n5d];
 else if(n3d>=0)e+=p->dangle3[type][n3d];
 if(type>2)e+=p->TerminalAU;
 return e;
}


/* Conservative integer lower bound for every internal-loop family reached by
 * rdo_E_internal_bounded: total unpaired <=30, default salt, no GU-closure
 * prohibition, pair types1..7 and base codes0..4. Independent minima can only
 * lower a sum; no sequence identity or predicted result enters this bound. */
static inline int rdo_E_internal_legal_floors_default(const vrna_param_t *p, int bounds[8][5][5][31]) {
  int stack[8],int11[8][5][5],int21_second[8][5][5],int21_first[8][5][5],int22[8][5][5];
  int mismatch1n_first[5],mismatch1n_second[5];
  for(int x=0;x<=4;x++) mismatch1n_first[x]=mismatch1n_second[x]=INF;
  int mismatch=INF,mismatch1n=INF,mismatch23=INF;
  for(int u=0;u<=7;u++) {
    stack[u]=INF;
    for(int x=0;x<=4;x++) {
      for(int y=0;y<=4;y++) {
        int11[u][x][y]=int21_second[u][x][y]=int21_first[u][x][y]=int22[u][x][y]=INF;
        for(int a=0;a<=MAXLOOP;a++) bounds[u][x][y][a]=INF;
      }
    }
  }
  for(int t=1;t<=7;t++) {
    for(int x=0;x<=4;x++) for(int y=0;y<=4;y++) {
      mismatch=MIN2(mismatch,p->mismatchI[t][x][y]);
      mismatch1n=MIN2(mismatch1n,p->mismatch1nI[t][x][y]);
      mismatch1n_first[x]=MIN2(mismatch1n_first[x],p->mismatch1nI[t][x][y]);
      mismatch1n_second[y]=MIN2(mismatch1n_second[y],p->mismatch1nI[t][x][y]);
      mismatch23=MIN2(mismatch23,p->mismatch23I[t][x][y]);
    }
    for(int u=1;u<=7;u++) {
      stack[u]=MIN2(stack[u],p->stack[t][u]);
      for(int x=0;x<=4;x++) for(int y=0;y<=4;y++) {
        int11[u][y][x]=MIN2(int11[u][y][x],p->int11[t][u][x][y]);
        for(int z=0;z<=4;z++) {
          int value=p->int21[t][u][x][y][z];
          int21_second[u][y][x]=MIN2(int21_second[u][y][x],value);
          int21_first[t][x][z]=MIN2(int21_first[t][x][z],value);
          for(int w=0;w<=4;w++) int22[u][z][y]=MIN2(int22[u][z][y],p->int22[t][u][x][y][z][w]);
        }
      }
    }
  }
  int global=INF;
  for(int u=1;u<=7;u++) for(int x=0;x<=4;x++) for(int y=0;y<=4;y++) for(int a=0;a<=MAXLOOP;a++) {
    int minimum=INF;
    for(int b=0;a+b<=MAXLOOP;b++) {
      int nl=MAX2(a,b),ns=MIN2(a,b),E;
      if(nl==0) E=stack[u]+p->SaltStack;
      else if(ns==0) E=p->bulge[nl]+((nl==1)?stack[u]:((u>2)?p->TerminalAU:0)+MIN2(0,p->TerminalAU));
      else if(ns==1) {
        if(nl==1) E=int11[u][x][y];
        else if(nl==2) E=(a==1)?int21_second[u][x][y]:int21_first[u][x][y];
        else E=p->internal_loop[a+b]+MIN2(MAX_NINIO,abs(a-b)*p->ninio[2])+((a==1)?mismatch1n_first[y]:mismatch1n_second[x])+p->mismatch1nI[u][x][y];
      } else if(ns==2 && nl==2) E=int22[u][x][y];
      else if(ns==2 && nl==3) E=p->internal_loop[5]+p->ninio[2]+mismatch23+p->mismatch23I[u][x][y];
      else E=p->internal_loop[a+b]+MIN2(MAX_NINIO,abs(a-b)*p->ninio[2])+mismatch+p->mismatchI[u][x][y];
      minimum=MIN2(minimum,E);
    }
    bounds[u][x][y][a] = minimum;
    global=MIN2(global,minimum);
  }
  return global;
}

static inline int rdo_E_internal_base_floors_default(const vrna_param_t *p, int bounds[8][5][5][31]) {
  int stack[8],int11[8],int21_second[8][5],int21_first[8][5][5],int22[8][5][5];
  int mismatch=INF,mismatch1n=INF,mismatch23=INF;
  for(int u=0;u<=7;u++) {
    stack[u]=int11[u]=INF;
    for(int x=0;x<=4;x++) {
      int21_second[u][x]=INF;
      for(int y=0;y<=4;y++) {
        int21_first[u][x][y]=int22[u][x][y]=INF;
        for(int a=0;a<=MAXLOOP;a++) bounds[u][x][y][a]=INF;
      }
    }
  }
  for(int t=1;t<=7;t++) {
    for(int x=0;x<=4;x++) for(int y=0;y<=4;y++) {
      mismatch=MIN2(mismatch,p->mismatchI[t][x][y]);
      mismatch1n=MIN2(mismatch1n,p->mismatch1nI[t][x][y]);
      mismatch23=MIN2(mismatch23,p->mismatch23I[t][x][y]);
    }
    for(int u=1;u<=7;u++) {
      stack[u]=MIN2(stack[u],p->stack[t][u]);
      for(int x=0;x<=4;x++) for(int y=0;y<=4;y++) {
        int11[u]=MIN2(int11[u],p->int11[t][u][x][y]);
        for(int z=0;z<=4;z++) {
          int value=p->int21[t][u][x][y][z];
          int21_second[u][y]=MIN2(int21_second[u][y],value);
          int21_first[t][x][z]=MIN2(int21_first[t][x][z],value);
          for(int w=0;w<=4;w++) int22[u][z][y]=MIN2(int22[u][z][y],p->int22[t][u][x][y][z][w]);
        }
      }
    }
  }
  int global=INF;
  for(int u=1;u<=7;u++) for(int x=0;x<=4;x++) for(int y=0;y<=4;y++) for(int a=0;a<=MAXLOOP;a++) {
    int minimum=INF;
    for(int b=0;a+b<=MAXLOOP;b++) {
      int nl=MAX2(a,b),ns=MIN2(a,b),E;
      if(nl==0) E=stack[u]+p->SaltStack;
      else if(ns==0) E=p->bulge[nl]+((nl==1)?stack[u]:((u>2)?p->TerminalAU:0)+MIN2(0,p->TerminalAU));
      else if(ns==1) {
        if(nl==1) E=int11[u];
        else if(nl==2) E=(a==1)?int21_second[u][x]:int21_first[u][x][y];
        else E=p->internal_loop[a+b]+MIN2(MAX_NINIO,abs(a-b)*p->ninio[2])+mismatch1n+p->mismatch1nI[u][x][y];
      } else if(ns==2 && nl==2) E=int22[u][x][y];
      else if(ns==2 && nl==3) E=p->internal_loop[5]+p->ninio[2]+mismatch23+p->mismatch23I[u][x][y];
      else E=p->internal_loop[a+b]+MIN2(MAX_NINIO,abs(a-b)*p->ninio[2])+mismatch+p->mismatchI[u][x][y];
      minimum=MIN2(minimum,E);
    }
    bounds[u][x][y][a] = minimum;
    global=MIN2(global,minimum);
  }
  return global;
}

static inline int rdo_E_internal_type_floors_default(const vrna_param_t *p, int bounds[8][31]) {
  int stack[8],int11[8],int21_second[8],int21_first[8],int22[8];
  int mismatch[8],mismatch1n[8],mismatch23[8];
  int pred_mismatch=INF,pred_mismatch1n=INF,pred_mismatch23=INF;
  for(int u=0;u<=7;u++) {
    stack[u]=int11[u]=int21_second[u]=int21_first[u]=int22[u]=INF;
    mismatch[u]=mismatch1n[u]=mismatch23[u]=INF;
    for(int a=0;a<=MAXLOOP;a++) bounds[u][a]=INF;
  }
  for(int t=1;t<=7;t++) {
    for(int x=0;x<=4;x++) for(int y=0;y<=4;y++) {
      mismatch[t]=MIN2(mismatch[t],p->mismatchI[t][x][y]);
      mismatch1n[t]=MIN2(mismatch1n[t],p->mismatch1nI[t][x][y]);
      mismatch23[t]=MIN2(mismatch23[t],p->mismatch23I[t][x][y]);
    }
    pred_mismatch=MIN2(pred_mismatch,mismatch[t]);
    pred_mismatch1n=MIN2(pred_mismatch1n,mismatch1n[t]);
    pred_mismatch23=MIN2(pred_mismatch23,mismatch23[t]);
    for(int u=1;u<=7;u++) {
      stack[u]=MIN2(stack[u],p->stack[t][u]);
      for(int x=0;x<=4;x++) for(int y=0;y<=4;y++) {
        int11[u]=MIN2(int11[u],p->int11[t][u][x][y]);
        for(int z=0;z<=4;z++) {
          int value=p->int21[t][u][x][y][z];
          int21_second[u]=MIN2(int21_second[u],value);
          int21_first[t]=MIN2(int21_first[t],value);
          for(int w=0;w<=4;w++) int22[u]=MIN2(int22[u],p->int22[t][u][x][y][z][w]);
        }
      }
    }
  }
  int global=INF;
  for(int u=1;u<=7;u++) for(int a=0;a<=MAXLOOP;a++) {
    int minimum=INF;
    for(int b=0;a+b<=MAXLOOP;b++) {
      int nl=MAX2(a,b),ns=MIN2(a,b),E;
      if(nl==0) E=stack[u]+p->SaltStack;
      else if(ns==0) E=p->bulge[nl]+((nl==1)?stack[u]:((u>2)?p->TerminalAU:0)+MIN2(0,p->TerminalAU));
      else if(ns==1) {
        if(nl==1) E=int11[u];
        else if(nl==2) E=(a==1)?int21_second[u]:int21_first[u];
        else E=p->internal_loop[a+b]+MIN2(MAX_NINIO,abs(a-b)*p->ninio[2])+pred_mismatch1n+mismatch1n[u];
      } else if(ns==2 && nl==2) E=int22[u];
      else if(ns==2 && nl==3) E=p->internal_loop[5]+p->ninio[2]+pred_mismatch23+mismatch23[u];
      else E=p->internal_loop[a+b]+MIN2(MAX_NINIO,abs(a-b)*p->ninio[2])+pred_mismatch+mismatch[u];
      minimum=MIN2(minimum,E);
    }
    bounds[u][a] = minimum;
    global=MIN2(global,minimum);
  }
  return global;
}

static inline int rdo_E_internal_floors_default(const vrna_param_t *p, int bounds[31]) {
  int stack=INF, int11=INF, int21=INF, int22=INF;
  int mismatch=INF, mismatch1n=INF, mismatch23=INF;
  for (int t=1;t<=7;t++) {
    for (int x=0;x<=4;x++) for (int y=0;y<=4;y++) {
      mismatch=MIN2(mismatch,p->mismatchI[t][x][y]);
      mismatch1n=MIN2(mismatch1n,p->mismatch1nI[t][x][y]);
      mismatch23=MIN2(mismatch23,p->mismatch23I[t][x][y]);
    }
    for (int u=1;u<=7;u++) {
      stack=MIN2(stack,p->stack[t][u]);
      for (int x=0;x<=4;x++) for (int y=0;y<=4;y++) {
        int11=MIN2(int11,p->int11[t][u][x][y]);
        for (int z=0;z<=4;z++) {
          int21=MIN2(int21,p->int21[t][u][x][y][z]);
          for (int w=0;w<=4;w++)
            int22=MIN2(int22,p->int22[t][u][x][y][z][w]);
        }
      }
    }
  }
  int floor=INF;
  for (int a=0;a<=MAXLOOP;a++) {
    int minimum=INF;
    for (int b=0;a+b<=MAXLOOP;b++) {
      int nl=MAX2(a,b),ns=MIN2(a,b),E;
      if (nl==0) E=stack+p->SaltStack;
      else if (ns==0) E=p->bulge[nl]+((nl==1)?stack:2*MIN2(0,p->TerminalAU));
      else if (ns==1) {
        if (nl==1) E=int11;
        else if (nl==2) E=int21;
        else E=p->internal_loop[a+b]+MIN2(MAX_NINIO,abs(a-b)*p->ninio[2])+2*mismatch1n;
      } else if (ns==2 && nl==2) E=int22;
      else if (ns==2 && nl==3) E=p->internal_loop[5]+p->ninio[2]+2*mismatch23;
      else E=p->internal_loop[a+b]+MIN2(MAX_NINIO,abs(a-b)*p->ninio[2])+2*mismatch;
      minimum=MIN2(minimum,E);
    }
    bounds[a] = minimum;
    floor=MIN2(floor,minimum);
  }
  return floor;
}

static inline int rdo_E_internal_floor_default(const vrna_param_t *p) {
  int bounds[31];
  return rdo_E_internal_floors_default(p,bounds);
}

#endif
