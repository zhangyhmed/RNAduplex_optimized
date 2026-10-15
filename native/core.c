/* Derived from ViennaRNA 2.7.2 duplex.c. Original notice follows.
 * Modifications: exact native MFE optimizations and reusable process workspace.
 * See NOTICE, energy.h and the engineering report for provenance and domains. */
#define duplexfold rdo_duplexfold
#define duplex_subopt rdo_duplex_subopt
#define aliduplexfold rdo_aliduplexfold
#define aliduplex_subopt rdo_aliduplex_subopt
int rdo_last_energy;
/*
 *         compute the duplex structure of two RNA strands,
 *              allowing only inter-strand base pairs.
 *       see cofold() for computing hybrid structures without
 *                           restriction.
 *
 *                           Ivo Hofacker
 *                        Vienna RNA package
 */

#ifdef _OPENMP
#error "RNAduplex_optimized core uses process workers; compile this unit without -fopenmp"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/subopt/wuchty.h"
#include "ViennaRNA/duplex.h"
#include "energy.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define STACK_BULGE1  1     /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1     /* new asymetry penalty */
#define MAXSECTORS    500   /* dimension for a backtrack array */
#define LOCALITY      0.    /* locality parameter for base-pairs */
#define UNIT 100
#define MINPSCORE -2 * UNIT
#define NONE -10000         /* score for forbidden pairs */


/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */
/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */
PRIVATE vrna_param_t  *P = NULL;
PRIVATE int           **c = NULL;     /* energy array, given that i-j pair */
PRIVATE short         *S1 = NULL, *SS1 = NULL, *S2 = NULL, *SS2 = NULL;
PRIVATE unsigned char *rdo_types;
PRIVATE uint16_t *rdo_next;
PRIVATE int rdo_stride;
PRIVATE int *rdo_cells;
PRIVATE size_t rdo_capacity, rdo_row_capacity;
PRIVATE int *rdo_generic;
PRIVATE int *rdo_suffix_min;
PRIVATE int rdo_energy_floor;
PRIVATE int rdo_floor_by_base[8][5][5][31];
PRIVATE int rdo_geometry[31][31], rdo_factor_ready, rdo_factor_valid;
PRIVATE int           n1, n2;         /* sequence lengths */

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
 */
#pragma omp threadprivate(P, c, S1, SS1, S2, SS2, n1, n2)

#endif

PRIVATE size_t rdo_enc1_capacity, rdo_enc2_capacity;
PRIVATE size_t rdo_st1_capacity, rdo_st2_capacity, rdo_result_capacity;
PRIVATE char *rdo_st1, *rdo_st2, *rdo_result;

PRIVATE unsigned rdo_encode_reuse(const char *sequence, int n, short **S, short **SS, size_t *capacity) {
  if ((size_t)(n+2) > *capacity) {
    *S = vrna_realloc(*S, (size_t)(n+2)*sizeof(short));
    *SS = vrna_realloc(*SS, (size_t)(n+2)*sizeof(short));
    *capacity = n+2;
  }
  unsigned observed=0;
  for (int i=1;i<=n;i++) {
    int code=encode_char(sequence[i-1]);
    observed |= 1u << code;
    (*S)[i]=code;
    (*SS)[i]=alias[code];
  }
  (*S)[n+1]=(*S)[1]; (*S)[0]=(short)n;
  (*SS)[n+1]=(*SS)[1]; (*SS)[0]=(*SS)[n];
  return observed;
}

/* Namespaced private API: result belongs to this per-process workspace and
 * remains valid until the next fold; callers copy/consume it before release. */
void rdo_release_result(char *result) {(void)result;}

PRIVATE char *backtrack(int i, int j);
PRIVATE duplexT duplexfold_cu(const char *, const char *, int);

PUBLIC duplexT
duplexfold(const char *s1,
           const char *s2)
{
  return duplexfold_cu(s1, s2, 1);
}


PRIVATE duplexT
duplexfold_cu(const char  *s1,
              const char  *s2,
              int         clean_up)
{
  int       i, j, Emin = INF, i_min = 0, j_min = 0;
  char      *struc;
  duplexT   mfe;
  vrna_md_t md;

  n1  = (int)strlen(s1);
  n2  = (int)strlen(s2);

  set_model_details(&md);
  if (!P) {
    P = vrna_params(&md);
    rdo_factor_ready = 0;
    make_pair_matrix();
  } else {
    md.window_size = P->model_details.window_size;
    md.max_bp_span = P->model_details.max_bp_span;
    /* check if model_details are the same as before */
    if (memcmp(&md, &(P->model_details), sizeof(vrna_md_t))) {
      free(P);
      P = vrna_params(&md);
      rdo_factor_ready = 0;
      make_pair_matrix();
    }
  }

  /* This optimized module exposes immutable default-model MFE only.
   * Configurable modes are handled by the independently built upstream CLI. */
  if (n1 > 32760 || n2 > 32760 || P->model_details.salt != VRNA_MODEL_DEFAULT_SALT || P->model_details.noGUclosure) {
    duplexT failed = {0};
    vrna_log_error("unsupported native model: use unchanged 2.7.2 backend");
    return failed;
  }

  unsigned observed1=rdo_encode_reuse(s1,n1,&S1,&SS1,&rdo_enc1_capacity);
  unsigned observed2=rdo_encode_reuse(s2,n2,&S2,&SS2,&rdo_enc2_capacity);
  int any_pair=0, all_pairs=1;
  for (unsigned first=observed1;first;first &= first-1) {
    int a=__builtin_ctz(first);
    for (unsigned second=observed2;second;second &= second-1) {
      int type=pair[a][__builtin_ctz(second)];
      any_pair |= type;
      all_pairs &= type != 0;
    }
  }
  /* Exact impossibility proof: the Cartesian product of observed encoded
   * symbols has no allowed pair. Every original DP cell would be INF;
   * invoke the same traceback below with the unchanged zero endpoint. */
  if (!any_pair) goto prediction_ready;

  rdo_stride = n2 + 2;
  size_t cells = (size_t)(n1 + 1) * rdo_stride;
  if (cells > rdo_capacity) {
    rdo_cells = vrna_realloc(rdo_cells, cells*sizeof(int));
    rdo_types = vrna_realloc(rdo_types, cells);
    rdo_next = vrna_realloc(rdo_next, cells*sizeof(uint16_t));
    rdo_generic = vrna_realloc(rdo_generic, cells*sizeof(int));
    rdo_suffix_min = vrna_realloc(rdo_suffix_min, cells*sizeof(int));
    rdo_capacity = cells;
  }
  if ((size_t)(n1+1) > rdo_row_capacity) {
    c = vrna_realloc(c, (size_t)(n1+1)*sizeof(int*));
    rdo_row_capacity = n1+1;
  }
  for (i=1;i<=n1;i++)c[i]=rdo_cells+(size_t)i*rdo_stride;


  if (!rdo_factor_ready) {
    rdo_factor_valid = P->model_details.salt == VRNA_MODEL_DEFAULT_SALT && !P->model_details.noGUclosure;
    for (int a = 0; a <= MAXLOOP; a++)
      for (int b = 0; a + b <= MAXLOOP; b++)
        rdo_geometry[a][b] = P->internal_loop[a+b] + MIN2(MAX_NINIO, abs(a-b)*P->ninio[2]);
    rdo_energy_floor = rdo_E_internal_legal_floors_default(P,rdo_floor_by_base);
    rdo_factor_ready = 1;
  }
  rdo_stride = n2 + 2;
  for (i = 1; i <= n1; i++) {
    int next = n2 + 1;
    rdo_next[i * rdo_stride + n2 + 1] = next;
    for (j = n2; j >= 1; j--) {
      int type = pair[S1[i]][S2[j]];
      rdo_types[i * rdo_stride + j] = type;
      rdo_generic[i * rdo_stride + j] = P->mismatchI[type][SS1[i+1]][SS2[j-1]];
      if (type) next = j;
      rdo_next[i * rdo_stride + j] = next;
    }
  }

  /* Dispatch once per fold. A complete observed-symbol pairing relation
   * means every cell is compatible, so contiguous l ascending visits exactly
   * the same transitions as the successor index. Sparse and traceback order
   * stay unchanged; there is no measured threshold or data-ID lookup. */
  if (all_pairs) {
  for (i = 1; i <= n1; i++) {
    int *suffix = rdo_suffix_min + i*rdo_stride;
    suffix[n2+1] = INF;
    for (j = n2; j > 0; j--) {
      int type, type2, E, k, l;
      type    = rdo_types[i * rdo_stride + j];
      c[i][j] = type ? P->DuplexInit : INF;
      suffix[j] = suffix[j+1];
      if (!type)
        continue;

      const int *floor_by_first = rdo_floor_by_base[rtype[type]][SS2[j+1]][SS1[i-1]];
      int current_generic = P->mismatchI[rtype[type]][SS2[j+1]][SS1[i-1]];
      c[i][j] += rdo_E_exterior_stem(type, (i > 1) ? SS1[i - 1] : -1, (j < n2) ? SS2[j + 1] : -1, P);
      int best = c[i][j];
      for (k = i - 1; k > 0 && k > i - MAXLOOP - 2; k--) {
        /* Suffix minimum includes every admissible l>j (and possibly more).
         * Its sum with the proven loop floor bounds this entire k block.
         * DP stores only energies; full original traceback preserves ties. */
        if (rdo_suffix_min[k*rdo_stride+j+1] + floor_by_first[i-k-1] >= best) continue;
        const int *previous = c[k];
        const uint16_t *nextrow = rdo_next + k*rdo_stride;
        const unsigned char *typerow = rdo_types + k*rdo_stride;
        const int *mismatchrow = rdo_generic + k*rdo_stride;
        for (l = j + 1; l <= MIN2(n2, MAXLOOP - i + k + j + 2); l++) {
          type2 = typerow[l];
          int a = i-k-1, b = l-j-1;
          if (((a >= 3 && b >= 3) || (a == 2 && b >= 4) || (b == 2 && a >= 4)))
            E = rdo_geometry[a][b] + (mismatchrow[l] + current_generic);
          else
            E = rdo_E_internal_bounded(a, b, type2, rtype[type],
                        SS1[k + 1], SS2[l - 1], SS1[i - 1], SS2[j + 1], P);
          best = MIN2(best, previous[l] + E);
        }
      }
      c[i][j] = best;
      suffix[j] = MIN2(best,suffix[j+1]);
      E = best;
      E += rdo_E_exterior_stem(rtype[type], (j > 1) ? SS2[j - 1] : -1, (i < n1) ? SS1[i + 1] : -1, P);
      if (E < Emin) {
        Emin  = E;
        i_min = i;
        j_min = j;
      }
    }
  }
  } else {
  for (i = 1; i <= n1; i++) {
    int *suffix = rdo_suffix_min + i*rdo_stride;
    suffix[n2+1] = INF;
    for (j = n2; j > 0; j--) {
      int type, type2, E, k, l;
      type    = rdo_types[i * rdo_stride + j];
      c[i][j] = type ? P->DuplexInit : INF;
      suffix[j] = suffix[j+1];
      if (!type)
        continue;

      const int *floor_by_first = rdo_floor_by_base[rtype[type]][SS2[j+1]][SS1[i-1]];
      int current_generic = P->mismatchI[rtype[type]][SS2[j+1]][SS1[i-1]];
      c[i][j] += rdo_E_exterior_stem(type, (i > 1) ? SS1[i - 1] : -1, (j < n2) ? SS2[j + 1] : -1, P);
      int best = c[i][j];
      for (k = i - 1; k > 0 && k > i - MAXLOOP - 2; k--) {
        /* Suffix minimum includes every admissible l>j (and possibly more).
         * Its sum with the proven loop floor bounds this entire k block.
         * DP stores only energies; full original traceback preserves ties. */
        if (rdo_suffix_min[k*rdo_stride+j+1] + floor_by_first[i-k-1] >= best) continue;
        const int *previous = c[k];
        const uint16_t *nextrow = rdo_next + k*rdo_stride;
        const unsigned char *typerow = rdo_types + k*rdo_stride;
        const int *mismatchrow = rdo_generic + k*rdo_stride;
        for (l = nextrow[j + 1]; l <= MIN2(n2, MAXLOOP - i + k + j + 2); l = nextrow[l + 1]) {
          type2 = typerow[l];
          int a = i-k-1, b = l-j-1;
          if (((a >= 3 && b >= 3) || (a == 2 && b >= 4) || (b == 2 && a >= 4)))
            E = rdo_geometry[a][b] + (mismatchrow[l] + current_generic);
          else
            E = rdo_E_internal_bounded(a, b, type2, rtype[type],
                        SS1[k + 1], SS2[l - 1], SS1[i - 1], SS2[j + 1], P);
          best = MIN2(best, previous[l] + E);
        }
      }
      c[i][j] = best;
      suffix[j] = MIN2(best,suffix[j+1]);
      E = best;
      E += rdo_E_exterior_stem(rtype[type], (j > 1) ? SS2[j - 1] : -1, (i < n1) ? SS1[i + 1] : -1, P);
      if (E < Emin) {
        Emin  = E;
        i_min = i;
        j_min = j;
      }
    }
  }
  }

prediction_ready:
  rdo_last_energy = Emin;
  struc = backtrack(i_min, j_min);
  if (i_min < n1)
    i_min++;

  if (j_min > 1)
    j_min--;

  mfe.i         = i_min;
  mfe.j         = j_min;
  mfe.energy    = (float)Emin / 100.;
  mfe.structure = struc;
  if (clean_up) {
  }

  return mfe;
}


PRIVATE char *
backtrack(int i,
          int j)
{
  /* backtrack structure going backwards from i, and forwards from j
   * return structure in bracket notation with & as separator */
  int   k, l, type, type2, E, traced, i0, j0;
  char  *st1, *st2, *struc;

  if ((size_t)(n1+1)>rdo_st1_capacity) {
    rdo_st1=vrna_realloc(rdo_st1,(size_t)n1+1);rdo_st1_capacity=n1+1;
  }
  if ((size_t)(n2+1)>rdo_st2_capacity) {
    rdo_st2=vrna_realloc(rdo_st2,(size_t)n2+1);rdo_st2_capacity=n2+1;
  }
  st1=rdo_st1;st2=rdo_st2;
  memset(st1,0,(size_t)n1+1);memset(st2,0,(size_t)n2+1);

  i0  = MIN2(i + 1, n1);
  j0  = MAX2(j - 1, 1);

  while (i > 0 && j <= n2) {
    E           = c[i][j];
    traced      = 0;
    st1[i - 1]  = '(';
    st2[j - 1]  = ')';
    type        = rdo_types[i * rdo_stride + j];
    if (!type) {
      vrna_log_error("backtrack failed in fold duplex");
      return NULL;
    }

    for (k = i - 1; k > 0 && k > i - MAXLOOP - 2; k--) {
      for (l = rdo_next[k * rdo_stride + j + 1]; l <= MIN2(n2, MAXLOOP - i + k + j + 2); l = rdo_next[k * rdo_stride + l + 1]) {
        int LE;
        if (i - k + l - j - 2 > MAXLOOP)
          break;

        type2 = rdo_types[k * rdo_stride + l];
        if (!type2)
          continue;

        LE = rdo_E_internal_bounded(i - k - 1, l - j - 1, type2, rtype[type],
                       SS1[k + 1], SS2[l - 1], SS1[i - 1], SS2[j + 1], P);
        if (E == c[k][l] + LE) {
          traced  = 1;
          i       = k;
          j       = l;
          break;
        }
      }
      if (traced)
        break;
    }
    if (!traced) {
      E -= rdo_E_exterior_stem(type, (i > 1) ? SS1[i - 1] : -1, (j < n2) ? SS2[j + 1] : -1, P);
      if (E != P->DuplexInit) {
        vrna_log_error("backtrack failed in fold duplex");
            return NULL;
      } else {
        break;
      }
    }
  }
  if (i > 1)
    i--;

  if (j < n2)
    j++;

  size_t result_size=i0-i+1+j-j0+1+2;
  if (result_size>rdo_result_capacity) {
    rdo_result=vrna_realloc(rdo_result,result_size);rdo_result_capacity=result_size;
  }
  struc=rdo_result;
  for (k = MAX2(i, 1); k <= i0; k++)
    if (!st1[k - 1])
      st1[k - 1] = '.';

  for (k = j0; k <= j; k++)
    if (!st2[k - 1])
      st2[k - 1] = '.';

  strcpy(struc, st1 + MAX2(i - 1, 0));
  strcat(struc, "&");
  strcat(struc, st2 + j0 - 1);

  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */

  return struc;
}


/*------------------------------------------------------------------------*/

