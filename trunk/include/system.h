#ifndef _tpm_system_h
#define _tpm_system_h
#include "types.h"

void 
init_system (System *tpm_system, char *input_file);

void 
relax_system(System *tpm_system, int graphic_enabled, GLWindow *gl_window);

void
release_system (System *tpm_system);
void
run_system(System *tpm_system, int verbose_enabled, int graphic_enabled, GLWindow *gl_window);

#endif //_tpm_system_h
