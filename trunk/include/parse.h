#ifndef _tpm_parse_h
#define _tpm_parse_h
void
parse_option (char **input_file, int *graphic_enabled, 
                    int *verbose_enabled, int argc, char** argv);
void 
parse_input_file (FILE *fp, System *tpm_system);

#endif //_tpm_parse_h
