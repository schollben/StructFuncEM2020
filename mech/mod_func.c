#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _CaT_reg();
extern void _ampanmda_reg();
extern void _ca_reg();
extern void _exp2synNMDA_reg();
extern void _h_reg();
extern void _kca_reg();
extern void _km_reg();
extern void _kv_reg();
extern void _na_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," CaT.mod");
fprintf(stderr," ampanmda.mod");
fprintf(stderr," ca.mod");
fprintf(stderr," exp2synNMDA.mod");
fprintf(stderr," h.mod");
fprintf(stderr," kca.mod");
fprintf(stderr," km.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," na.mod");
fprintf(stderr, "\n");
    }
_CaT_reg();
_ampanmda_reg();
_ca_reg();
_exp2synNMDA_reg();
_h_reg();
_kca_reg();
_km_reg();
_kv_reg();
_na_reg();
}
