// #include "pch.h"
#include "suitesparse_vswrap.h"

namespace {

HINSTANCE dllHandle = LoadLibrary(TEXT("suitesparse_wrap.dll"));

}

cholmodw_initTp cholmodw_init = (cholmodw_initTp)GetProcAddress(dllHandle, "cholmodw_init");
cholmodw_setmat_tripletsTp cholmodw_setmat_triplets = (cholmodw_setmat_tripletsTp)GetProcAddress(dllHandle, "cholmodw_setmat_triplets");
cholmodw_setmat_rowmajTp cholmodw_setmat_rowmaj = (cholmodw_setmat_rowmajTp)GetProcAddress(dllHandle, "cholmodw_setmat_rowmaj");
cholmodw_setmat_csrTp cholmodw_setmat_csr = (cholmodw_setmat_csrTp)GetProcAddress(dllHandle, "cholmodw_setmat_csr");
cholmodw_solveTp cholmodw_solve = (cholmodw_solveTp)GetProcAddress(dllHandle, "cholmodw_solve");
cholmodw_print_infoTp cholmodw_print_info = (cholmodw_print_infoTp)GetProcAddress(dllHandle, "cholmodw_print_info");
cholmodw_freeTp cholmodw_free = (cholmodw_freeTp)GetProcAddress(dllHandle, "cholmodw_free");