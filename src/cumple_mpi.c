#define _POSIX_C_SOURCE 200809L
#include <mpi.h>

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rng.h"



enum { DIAS_EN_ANIO = 365, PALABRAS_BITSET = 6 };


typedef enum { MODO_RNG_RANK, MODO_RNG_TRIAL } modo_rng;

typedef struct {
    uint64_t tam_grupo;    
    uint64_t total_ensayos; 
    uint64_t semilla;      
    bool formato_csv;      
    bool ayuda;             
    modo_rng modo;          
} opciones_programa;

static __attribute__((unused)) bool parsear_u64(const char *texto, uint64_t *valor)
{
    char *fin = NULL;
    unsigned long long parseado;

    if (texto == NULL || *texto == '\0' || *texto == '-') return false;
    errno = 0;
    parseado = strtoull(texto, &fin, 10);
    if (errno == ERANGE || fin == texto || *fin != '\0') return false;
    if (parseado > UINT64_MAX) return false;
    *valor = (uint64_t)parseado;
    return true;
}


static __attribute__((unused)) void mostrar_ayuda(FILE *flujo, const char *nombre_programa)
{
    fprintf(flujo,
            "Uso: %s <K> <N> <SEED> [--csv] [--rng rank|trial]\n"
            "     %s --help\n\n"
            "Argumentos:\n"
            "  K    : Tamaño del grupo de personas (ej. 23)\n"
            "  N    : Total de simulaciones Monte Carlo (ej. 10000000)\n"
            "  SEED : Semilla numérica inicial (uint64)\n\n"
            "Opciones:\n"
            "  --csv         : Salida en formato CSV para scripts\n"
            "  --rng rank    : Generador único por proceso MPI (por defecto)\n"
            "  --rng trial   : Generador independiente del número de procesos\n",
            nombre_programa, nombre_programa);
}



int main(int argc, char **argv)
{
   
    int id_proceso, cantidad_procesos;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id_proceso);
    MPI_Comm_size(MPI_COMM_WORLD, &cantidad_procesos);

    if (id_proceso == 0) {
        printf("MPIinicializado correctamente. Procesos activos: %d\n", cantidad_procesos);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}