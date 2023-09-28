
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

double static **Marginals(double* theta, int maxsize)
{
 double **res;
 int r, n;
 
 res = malloc((maxsize+1)*sizeof(double *));
 for (n=0; n<=maxsize; n++){
         res[n] = calloc(n+1, sizeof(double));
 }
 
 for (r=0; r<=maxsize; r++){
         res[maxsize][r] = theta[r];
 } 
 
 for (n=maxsize-1; n>=1; n--){
         for (r=0; r<=n; r++){
                 res[n][r] = (n+1.0-r)/(n+1.0)*res[n+1][r] + (r+1.0)/(n+1.0)*res[n+1][r+1];
         }
 }
 
 return res;    
}

SEXP ReprodEstimates(SEXP nvec, SEXP rvec, SEXP freqvec)
{
 double *theta, *thetanew, abserror, **marg;
 int i, maxsize, nr, ntot, r, ri, ni, n, fri, start;
 SEXP res;
 const double eps=1e-16;
 
 nr = LENGTH(nvec);
 maxsize = 0;
 ntot = 0;
 for (i=0; i<nr; i++){
         if (INTEGER(nvec)[i]>maxsize){
                 maxsize = INTEGER(nvec)[i];
         }
         ntot += INTEGER(freqvec)[i];
 }
                 
 theta = malloc((maxsize+1) * sizeof(double));
 thetanew = malloc((maxsize+1) * sizeof(double));
 //starting values
 for (r=0; r<=maxsize; r++){
         theta[r] = 1.0/(maxsize+1);
 }
 abserror = 1;
 //EM update
 while (abserror>eps){
         abserror = 0;
         marg = Marginals(theta, maxsize);
         for (r=0; r<=maxsize; r++) thetanew[r] = 0;
         for (i=0; i<nr; i++){
                 ri = INTEGER(rvec)[i];
                 ni = INTEGER(nvec)[i];
                 fri = INTEGER(freqvec)[i];
                 for (r=ri; r<=maxsize-ni+ri; r++){
                         thetanew[r] += choose(ni,ri)*choose(maxsize-ni,r-ri)*theta[r]*fri*1.0/
                                        marg[ni][ri] ;
                 }
         }
         for (r=0; r<=maxsize; r++){
                 thetanew[r] = thetanew[r]/(ntot*choose(maxsize,r)*1.0);
                 abserror += fabs(thetanew[r]-theta[r]);
                 theta[r] = thetanew[r];
         }
         for(n = 0; n <= maxsize; n++) free(marg[n]);
     free(marg);
 }

 PROTECT(res = allocMatrix(REALSXP, maxsize+1, maxsize));
 marg = Marginals(theta, maxsize);
 for (n=1, start=0; n<=maxsize; n++, start+=(maxsize+1)){
         for (r=0; r<=n; r++){
                 REAL(res)[start+r] = marg[n][r];
         }
         for (r=n+1; r<=maxsize; r++){
                 REAL(res)[start+r] = NA_REAL;
         }
 }
 
 for(n = 0; n <= maxsize; n++) free(marg[n]);
 free(marg);
 free(theta);
 free(thetanew);
 
 UNPROTECT(1);
 return res;
        
}

void Comb(int j, int m, int nn, int kk, int nS, int* a, int* pos, SEXP res) {
        int i, val, step;
        if (j > nn) {  
        
                        val = 0;
                        step = 0;
                        for (i=1; i<=nn; i++) {
                                if (a[i]==1){
                                        INTEGER(res)[*pos+step]=val;
                                        step += nS;
                                }
                                else { //a[i]=0;
                                        val++;
                                }
                        }
                    *pos = *pos+1;
        
    }           
        else {
                if (kk-m < nn-j+1) {
                        a[j] = 0; Comb(j+1, m, nn, kk, nS, a, pos, res);
                }
                if (m<kk) {
                        a[j] = 1; Comb(j+1, m+1, nn, kk, nS, a, pos, res);
                }
        }
}

SEXP makeSmatrix(SEXP size, SEXP ntrt){
int *a, i, pos, nn, kk, nS;
SEXP res;
   nn = asInteger(size) + asInteger(ntrt);
   kk = asInteger(ntrt);
   a = calloc(nn+1, sizeof(int));
   for (i=1; i<=kk; i++) a[i]=1;

   nS = choose(nn,kk);
   PROTECT(res = allocMatrix(INTSXP, nS, kk));
   pos = 0;
   
   Comb(1, 0, nn, kk, nS, a, &pos, res);
  
   UNPROTECT(1);
   free(a);
 
  return res;
}

double GetTabElem(SEXP tab, int size, int n, int r, int j){
         return REAL(tab)[(n-1)+size*(r+(size+1)*j)];
 }

double ***HyperTable(int size){
// dhyper(i, j, size-j, k), i=0:size; j=0:size; k=0:size        
 double ***res;
 int i, j, k;
 
 res = malloc((size+1)*sizeof(double*));
 for (i=0; i<=size; i++){
         res[i] = malloc((size+1)*sizeof(double*));
         for (j=0; j<=size; j++) res[i][j] = calloc(size+1, sizeof(double));
 }
 
 for (i=0; i<=size; i++){
         for (j=i; j<=size; j++){
                 for (k=i; k<=size-j+i; k++)
                         res[i][j][k] = dhyper(i, j, size-j, k, 0); 
                 }
         }

 return res;    
}

int static *IndexVectorC(int **S, int N, int G, int nrowS){
 int *idx, i, j;
 
 idx = calloc(nrowS, sizeof(int));
 
 for (j=G-1; j>=0; j--){
   for (i=0; i<nrowS; i++){
           idx[i] = N*idx[i] + S[i][j];
   }
 }
 
 return idx;
}
int static *IndexVector(SEXP S, int N, int G, int nrowS){
 int *idx, i, j,  start;
 
 idx = calloc(nrowS, sizeof(int));
 
 for (j=G-1, start=(G-1)* nrowS; j>=0; j--, start -= nrowS){
   for (i=0; i<nrowS; i++){
           idx[i] = N*idx[i] + INTEGER(S)[start+i];
   }
 }
 
 return idx;
}

double ***CalcMarginals(SEXP S, SEXP Q, double ***ht, int *idx, int ntrt, int size, int nS){
   int j, i, n, x, start, sj;
         double ***marg;
         //ht is from HyperTable, it is passed to avoid recalculation
         
         marg = malloc(ntrt*sizeof(double*));
         for (j=0; j<ntrt; j++){
                 marg[j] = malloc((size+1)*sizeof(double*));
                 for (n=1; n<=size; n++) marg[j][n] = calloc(n+1, sizeof(double));
         }
         
         for (i=0; i<nS; i++){
                 for (j=0, start=0; j<ntrt; j++, start+=nS){
                         sj = INTEGER(S)[start+i];
                         for (n=1; n<=size; n++){
                                 for (x=imax2(0,sj-size+n); x<=imin2(sj,n); x++){
                                         marg[j][n][x] += REAL(Q)[idx[i]]*ht[x][n][sj];
                                 }
                         }
                 }
         }
  return marg;
 }

void CalcD(SEXP D, SEXP S, SEXP tab, int *idx, double ***ht, double ***marg, int ntrt, 
           int nS, int size, int ntot){
   int j, i, n, x, sj, start, t;
         //ntot = sum(tab) -- passing it avoids having to recalculate it
         for (i=0; i<nS; i++){
                 REAL(D)[idx[i]] = -ntot;
                 for (j=0, start=0; j<ntrt; j++, start+=nS){
                         sj = INTEGER(S)[start+i];
                         for (n=1; n<=size; n++){
                                 for (x=imax2(0,sj-size+n); x<=imin2(sj,n); x++){
                                         t = GetTabElem(tab,size,n,x,j);
                                         if (t>0) REAL(D)[idx[i]] += t*ht[x][n][sj]/marg[j][n][x];
                                 }
                         }
                 }
         }
         
 }

double maxD(SEXP D, int *idx, int nS){
   int i;
   double currmax, val;
   
   currmax = 0;
   for (i=0; i<nS; i++){
           val = REAL(D)[idx[i]];
           if (val > currmax) currmax = val; 
   }
   return currmax;
}

int **CalcTopD(SEXP D, SEXP S, int *idx, int limit, int *nselect, int ntrt, int nS){
         int **res, nmx, i, j, pos, start;
         double *posDvec,  dcut;
         
           // count the number of non-negative elements in D
     nmx = 0;
     for (i=0; i<nS; i++){
             if (REAL(D)[idx[i]] >= 0){
                     nmx++;
             }
     }
     if (nmx == 0){
             res = 0;
             *nselect = 0;
             return res;
     }
     
     if (nmx > limit){ //find the limit-th largest D
                 posDvec = malloc(nmx*sizeof(double));
                 pos = 0;
             for (i=0; i<nS; i++){
                     if (REAL(D)[idx[i]] >= 0){
                             posDvec[pos] = -REAL(D)[idx[i]];  //negation is needed, because rPsort uses increasing order
                             pos++;
                     }
             }
             rPsort(posDvec, nmx, limit);
             dcut = -posDvec[limit];  //the cutoff for determining the limit-th largest values
             free(posDvec);
     }
     else dcut = 0;
         
     nmx = imin2(limit, nmx);

         res = (int**) Calloc(nmx, int*);
         pos = 0;
     for (i=0; i<nS; i++){
             if (pos >= nmx) break;
             if (REAL(D)[idx[i]] < dcut) continue;
             res[pos] =(int*) Calloc(ntrt, int);  //copy the ith row of S
             for (j=0, start=0; j<ntrt; j++, start+=nS) res[pos][j] = INTEGER(S)[start+i];
             pos++;
     }

     *nselect = nmx;
     return res;         
 } 


int ntrt, size, **lmS;
double ntot, ***ht, ***marg;

double NegLogLik(int npar, double *par, void *ex){
         //par[j] = (alpha_(j+1)/alpha_0), j=0,...,nmax-1
   int j, n, r, i, sj, x;
   double res, sum;
   SEXP tab;
   
   tab = (SEXP)ex;
   res = 0;
   
   for (j=0; j<ntrt; j++){
           for (n=1; n<=size; n++){
                   for (r=0; r<=n; r++){
                           x = GetTabElem(tab,size,n,r,j);
                           if (x>0){
                                   sum = marg[j][n][r];
                                   for (i=0; i<npar; i++){
                                           sj = lmS[i][j];
                                           sum += par[i]*ht[r][n][sj];
                                   }
                                   res += x*log(sum);    
                           }
                   }
           }
   }
   sum = 0;
   for (i=0; i<npar; i++) sum += par[i];
   res -= ntot*log1p(sum);  //log1p(sum)=log(1+sum)
   
   if (!R_FINITE(res)){
     res = 1e60;
   }
   
   return (-res); }

void NegLogLikDeriv(int npar, double *par, double *gr, void *ex){
   int j, n, r, i, sj, x;
   double alpha0, sum, ***denom;
   SEXP tab;
      
   tab = (SEXP)ex;
   //prepare the shared denominators
   denom = malloc(ntrt*sizeof(double*));
   for (j=0; j<ntrt; j++){
           denom[j] = malloc((size+1)*sizeof(double*));
           for (n=1; n<=size; n++) denom[j][n] = calloc(n+1, sizeof(double));
   }
   for (j=0; j<ntrt; j++){
           for (n=1; n<=size; n++){
                   for (r=0; r<=n; r++){
                           sum = marg[j][n][r];
                           for (i=0; i<npar; i++){
                                   sj = lmS[i][j];
                                   sum += par[i]*ht[r][n][sj];
                           }
                           denom[j][n][r] = sum;         
                   }
           }
   }
   
   alpha0 = 1;
   for (i=0; i<npar; i++)  alpha0 += par[i];
   alpha0 = 1.0/alpha0;

   //calc the gradients 
   for (i=0; i<npar; i++){
           sum = -ntot*alpha0;
           for (j=0; j<ntrt; j++){
                   for (n=1; n<=size; n++){
                           for (r=0; r<=n; r++){
                                   x = GetTabElem(tab,size,n,r,j);
                                   if (x>0){
                                           sj = lmS[i][j];
                                           sum += x*ht[r][n][sj]/denom[j][n][r];
                                   }
                           }
                   }
           }
           gr[i] = -sum;
   }
   
   for (j=0; j<ntrt; j++){
           for (n=1; n<=size; n++) free(denom[j][n]);
           free(denom[j]);
   }
   free(denom); 
}

void UpdateQ(SEXP Q, double *g, int nS, int nmax, int *idx, int *lmS_idx){
   double alpha0;
   int i;
   
   alpha0 = 1;
   for (i=0; i<nmax; i++) alpha0 += g[i];
   alpha0 = 1/alpha0;
   
   for (i=0; i<nS; i++) REAL(Q)[idx[i]] *= alpha0;
   for (i=0; i<nmax; i++) REAL(Q)[lmS_idx[i]] += alpha0 * g[i];
}

 void UpdateMarginals(double ***marg, double *g, double ***ht, int **lmS,
                      int ntrt, int size, int nmax){
   double alpha0;
   int i, j, n, r, sj;
    
    alpha0 = 1;
        for (i=0; i<nmax; i++){
                alpha0 += g[i];
        }
        alpha0 = 1/alpha0;
        
        for (j=0; j<ntrt; j++){
                for (n=1; n<=size; n++){
                        for (r=0; r<=n; r++){
                                for (i=0; i<nmax; i++){
                                        sj = lmS[i][j];
                                        marg[j][n][r] += g[i] * ht[r][n][sj];
                                }
                                marg[j][n][r] *= alpha0;
                        }
                }
        }        
}

SEXP ReprodISDM(SEXP Q, SEXP S, SEXP tab, SEXP MaxIter, SEXP MaxDirections, 
                SEXP eps, SEXP verbose){
 SEXP dims, D,  res, margSXP, tmp;
 int i, j, n, r, *idx, nS, niter, nmax, fncount, grcount, fail, 
     *boundtype,  limit, *lmS_idx, nenforced;
 double rel_error, *gamma, *lower, *upper, NLLmin;
 char msg[60];


 PROTECT(dims = GET_DIM(tab));
   size = INTEGER(dims)[0];
   ntrt = INTEGER(dims)[2];
 UNPROTECT(1);
 
 PROTECT(dims = GET_DIM(S));
   nS = INTEGER(dims)[0];
 UNPROTECT(1);

 
 PROTECT(D = duplicate(Q));
 for (i=0; i<length(Q); i++) REAL(D)[i] = 0;
 
 idx = IndexVector(S, size+1, ntrt, nS);
 ht = HyperTable(size);

 ntot=0;
 for (i=0; i<length(tab); i++) ntot += REAL(tab)[i];

 limit = INTEGER(coerceVector(MaxDirections, INTSXP))[0];
 if (limit <= 0){  //set it to the number of non-empty cells
   limit=0;
   for (i=0; i<length(tab); i++){
           if (REAL(tab)[i]>0) limit++;
    }
 }

 marg = CalcMarginals(S, Q, ht, idx, ntrt, size, nS);
 CalcD(D, S, tab, idx, ht, marg, ntrt, nS, size, ntot);
 
 rel_error = maxD(D, idx, nS);
 niter = 0;
 nenforced = 0;
 
 while (niter < asInteger(MaxIter) && rel_error > asReal(eps)){
         R_CheckUserInterrupt();
         niter++;
     
         lmS = CalcTopD(D, S, idx, limit, &nmax, ntrt, nS);
         lmS_idx = IndexVectorC(lmS, size+1, ntrt, nmax);
         
         if (nmax == limit) nenforced++;
         
         gamma = (double*) Calloc(nmax, double); 
         lower = (double*) Calloc(nmax, double); 
         upper = (double*) Calloc(nmax, double); 
         boundtype = (int*) Calloc(nmax, int); 
         
         for (i=0; i<nmax; i++){
                 gamma[i] = 0;
                 lower[i] = 0;
                 upper[i] = imin2(1e6/nmax, 100); // => alpha0>1e-6
                 boundtype[i] = 1; //lower  bound only
         }
         
         
         lbfgsb(nmax, 5, gamma, lower, upper, boundtype, &NLLmin, NegLogLik, NegLogLikDeriv,
                &fail, tab, 1e5, 0, &fncount, &grcount, 1000, msg, asInteger(verbose), 10);
                
         UpdateMarginals(marg, gamma, ht, lmS, ntrt, size, nmax);
   CalcD(D, S, tab, idx, ht, marg, ntrt, nS, size, ntot);
         UpdateQ(Q, gamma, nS, nmax, idx, lmS_idx); //only needed to be able to return Q
         rel_error = maxD(D, idx, nS);
         if (asInteger(verbose)==1)
           Rprintf("Step %d, rel.error=%f, NLL=%f\n", niter, rel_error, NLLmin);

                 
         Free(gamma);
         Free(lower);
         Free(upper);
         Free(boundtype);
         for (i=0; i<nmax; i++){
                 Free(lmS[i]);
         }
         Free(lmS);
         free(lmS_idx);
 }
 
 free(idx);
 for(j=0; j<=size; j++){
          for (n=0; n<=size; n++) free(ht[j][n]);
          free(ht[j]);
  }
 free(ht);
 
 PROTECT(margSXP = allocVector(REALSXP, ntrt*size*(size+1)));
 PROTECT(dims = allocVector(INTSXP,3));
    INTEGER(dims)[0] = size+1;     INTEGER(dims)[1] = size;     INTEGER(dims)[2] = ntrt;   
    i = 0;
    for(j=0; j<ntrt; j++){
          for (n=1; n<=size; n++){ 
                  for (r=0; r<=n; r++){
                          REAL(margSXP)[i] = marg[j][n][r];
                          i++;
                  }
                  for (r=n+1; r<=size; r++){
                          REAL(margSXP)[i] = NA_REAL;
                          i++;
                  }
          }
   }
 dimgets(margSXP, dims);
 UNPROTECT(1);    //dims
 for(j=0; j<ntrt; j++){
          for (n=1; n<=size; n++) free(marg[j][n]);
          free(marg[j]);
  }
 free(marg);

 PROTECT(res = allocVector(VECSXP,5));
  if (asInteger(verbose)==1) 
   Rprintf("Limit=%d; %d Iterations; Limit enforced %d times (%4.2f percent)\n", 
           limit, niter, nenforced, nenforced*100.0/niter);

  SET_VECTOR_ELT(res, 0, margSXP);
  SET_VECTOR_ELT(res, 1, Q);
  SET_VECTOR_ELT(res, 2, D);
  SET_VECTOR_ELT(res, 3, ScalarReal(-NLLmin));
  PROTECT(tmp = allocVector(REALSXP,2));
   REAL(tmp)[0] = rel_error;
   REAL(tmp)[1] = niter;
  SET_VECTOR_ELT(res, 4, tmp);
 UNPROTECT(4);   //tmp, res, margSXP, D

 return res;

} 


SEXP UpdateReprodQ(SEXP Q, SEXP S, SEXP tab, int size, int ntrt, int nS, 
                   double*** ht, int* idx){
 int ntot, i, j, n;
 double ***marg;
 SEXP resobj, D;
 
 PROTECT(resobj = duplicate(Q));
 for (i=0; i<length(Q); i++) REAL(resobj)[i] = 0;
 
 ntot=0;
 for (i=0; i<length(tab); i++) ntot += REAL(tab)[i];
 
 marg = CalcMarginals(S, Q, ht, idx, ntrt, size, nS);

 PROTECT(D = duplicate(Q));
 for (i=0; i<length(Q); i++) REAL(D)[i] = 0;
 CalcD(D, S, tab, idx, ht, marg, ntrt, nS, size, ntot);
  
// update Q values  
 for (i=0; i<length(Q); i++){
   REAL(resobj)[i] = REAL(Q)[i] * (1+REAL(D)[i]/ntot);
 }

 //cleanup
 for(j=0; j<ntrt; j++){
          for (n=1; n<=size; n++) free(marg[j][n]);
          free(marg[j]);
 }
 free(marg);

 
 UNPROTECT(2);
 return resobj;
}

SEXP MixReprodQ(SEXP Q, SEXP S, SEXP tab, SEXP MaxIter, SEXP eps, SEXP verbose){
 double rel_error, ***ht, ntot, re, ***marg, loglik;
 int niter, size, ntrt, nS, i, j, n, r, *idx, Qlength;
 SEXP D, Qnew, resQ, dims, resobj, tmp, margSXP, Qdims;

 PROTECT(Qdims = GET_DIM(Q));
 Qlength = 1;
 for (i=0; i<length(Qdims); i++) Qlength *= INTEGER(Qdims)[i];
 
 PROTECT(resQ = duplicate(Q));
 
 PROTECT(dims = GET_DIM(tab));
   size = INTEGER(dims)[0];
   ntrt = INTEGER(dims)[2];
 UNPROTECT(1);  //dims
 
 PROTECT(dims = GET_DIM(S));
   nS = INTEGER(dims)[0];
 UNPROTECT(1); //dims
 
 ntot=0;
 for (i=0; i<length(tab); i++) ntot += REAL(tab)[i];

 idx = IndexVector(S, size+1, ntrt, nS);
 ht = HyperTable(size);

 PROTECT(tmp = allocVector(REALSXP, 2));
 rel_error = 1;
 niter = 0;
 while ((niter<asInteger(MaxIter))&&(rel_error>asReal(eps))){
          R_CheckUserInterrupt();
                niter++;
    PROTECT(Qnew = UpdateReprodQ(resQ, S, tab, size, ntrt, nS, ht, idx));
    rel_error = 0;
    for (i=0; i<length(Qnew); i++){
            if (REAL(resQ)[i]>0){
                    re = ntot*(REAL(Qnew)[i]-REAL(resQ)[i])/REAL(resQ)[i];
                    if (rel_error < re) rel_error = re;
            }
            REAL(resQ)[i] = REAL(Qnew)[i];
    }
    UNPROTECT(1); //Qnew
    if ((asInteger(verbose) == 1)&&(niter%10 == 1)){
            REAL(tmp)[1] = rel_error;
            REAL(tmp)[0] = niter;
        PrintValue(tmp);
     } 

  }     
 UNPROTECT(1);  //tmp   
 
 //calculate ML estimates for the output
 marg = CalcMarginals(S, resQ, ht, idx, ntrt, size, nS);
 PROTECT(margSXP = allocVector(REALSXP, ntrt*size*(size+1)));
 PROTECT(dims = allocVector(INTSXP,3));
    INTEGER(dims)[0] = size+1;     INTEGER(dims)[1] = size;     INTEGER(dims)[2] = ntrt;   
    i = 0;
    for(j=0; j<ntrt; j++){
          for (n=1; n<=size; n++){ 
                  for (r=0; r<=n; r++){
                          REAL(margSXP)[i] = marg[j][n][r];
                          i++;
                  }
                  for (r=n+1; r<=size; r++){
                          REAL(margSXP)[i] = NA_REAL;
                          i++;
                  }
          }
   }            
 dimgets(margSXP, dims);   
 UNPROTECT(1); //dims

 PROTECT(D = allocVector(REALSXP, Qlength));
 dimgets(D, Qdims);
 for (i=0; i<Qlength; i++) REAL(D)[i] = 0;
 CalcD(D, S, tab, idx, ht, marg, ntrt, nS, size, ntot); 

 //calculate log-likelihood for the output
 loglik = 0;
 for(j=0; j<ntrt; j++){
  for (n=1; n<=size; n++){ 
          for (r=0; r<=n; r++){
            loglik += GetTabElem(tab, size, n, r, j)*log(marg[j][n][r]);
     }
   }
  }

 for(j=0; j<ntrt; j++){
          for (n=1; n<=size; n++) free(marg[j][n]);
          free(marg[j]);
  }
 free(marg);


 for(j=0; j<=size; j++){
          for (n=0; n<=size; n++) free(ht[j][n]);
          free(ht[j]);
  }
 free(ht);
 free(idx);
 
 PROTECT(resobj = allocVector(VECSXP, 5));
   SET_VECTOR_ELT(resobj, 0, margSXP);
   SET_VECTOR_ELT(resobj, 1, resQ);
   SET_VECTOR_ELT(resobj, 2, D);
   SET_VECTOR_ELT(resobj, 3, ScalarReal(loglik));
         PROTECT(tmp = allocVector(REALSXP,2));
           REAL(tmp)[0] = rel_error;
           REAL(tmp)[1] = niter;
         SET_VECTOR_ELT(resobj, 4, tmp);
 UNPROTECT(6); //tmp, resobj,D, margSXP, resQ, Qdims
 
 return resobj;

}

