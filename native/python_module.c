/* RNAduplex_optimized CPython interface. See COPYING.ViennaRNA and NOTICE. */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <ctype.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <ViennaRNA/duplex.h>
#include <ViennaRNA/vrna_config.h>
#ifndef RDO_REVISION
#define RDO_REVISION "unbuilt"
#endif
#ifndef RDO_BUILD_FLAGS
#define RDO_BUILD_FLAGS "unspecified"
#endif

extern duplexT rdo_duplexfold(const char *, const char *);
extern int rdo_last_energy;
extern void rdo_release_result(char *);

static PyObject *prediction(PyObject *self, PyObject *args) {
  const char *a, *b;
  Py_ssize_t na, nb;
  if (!PyArg_ParseTuple(args,"s#s#",&a,&na,&b,&nb)) return NULL;
  if (na == 0 || nb == 0 || na > 32760 || nb > 32760) {
    PyErr_SetString(PyExc_ValueError,"sequences must contain 1..32760 letters");
    return NULL;
  }
  for (Py_ssize_t i=0;i<na;i++) if (!isalpha((unsigned char)a[i]) || (unsigned char)a[i]>127) {
    PyErr_SetString(PyExc_ValueError,"sequences must contain ASCII letters only"); return NULL;
  }
  for (Py_ssize_t i=0;i<nb;i++) if (!isalpha((unsigned char)b[i]) || (unsigned char)b[i]>127) {
    PyErr_SetString(PyExc_ValueError,"sequences must contain ASCII letters only"); return NULL;
  }
  duplexT r = rdo_duplexfold(a,b);
  if (!r.structure) {PyErr_SetString(PyExc_RuntimeError,"native traceback failed");return NULL;}
  int native_energy=rdo_last_energy;
  int first = (int)(strchr(r.structure,'&')-r.structure);
  int ls=r.i+1-first, re=r.j+(int)strlen(r.structure)-first-2;
  size_t cap=strlen(r.structure)+160;
  char *text=PyMem_Malloc(cap);
  if (!text) {rdo_release_result(r.structure);return PyErr_NoMemory();}
  char energy[64];
  snprintf(energy,sizeof(energy),"%.2f",r.energy);
  snprintf(text,cap,"%s %3d,%-3d : %3d,%-3d (%5.2f)\n",r.structure,ls,r.i,r.j,re,r.energy);
  /* Copy the borrowed result before any GC-tracked container allocation can
   * run a Python callback. Scalars and formatted text are already local. */
  PyObject *structure=PyUnicode_FromString(r.structure);
  rdo_release_result(r.structure);
  if (!structure) {PyMem_Free(text);return NULL;}
  PyObject *result=Py_BuildValue("{s:O,s:i,s:i,s:i,s:i,s:s,s:i,s:s}",
    "structure",structure,"left_start",ls,"left_end",r.i,"right_start",r.j,
    "right_end",re,"energy",energy,"native_energy",native_energy,"text",text);
  Py_DECREF(structure);PyMem_Free(text);return result;
}

static PyObject *build_info(PyObject *self, PyObject *args) {
  return Py_BuildValue("{s:s,s:s,s:s,s:s}","revision",RDO_REVISION,
                      "viennarna",VRNA_VERSION,"compiler",__VERSION__,"flags",RDO_BUILD_FLAGS);
}

/* Keep the GIL: the legacy core owns mutable process-local state. Batch
 * concurrency uses independent processes, never concurrent native calls. */
static PyMethodDef methods[]={
 {"predict",prediction,METH_VARARGS,"Full default-model prediction, including native integer energy and exact text."},
 {"build_info",build_info,METH_NOARGS,"Compiled source identity and ViennaRNA header version."},
 {NULL,NULL,0,NULL}
};
static struct PyModuleDef module={PyModuleDef_HEAD_INIT,"_native",NULL,-1,methods};
PyMODINIT_FUNC PyInit__native(void) {return PyModule_Create(&module);}
