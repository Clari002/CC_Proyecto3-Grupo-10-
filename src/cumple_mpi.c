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
static __attribute__((unused)) bool parsear_opciones(int argc, char **argv, opciones_programa *opts, const char **error) {
    int posicionales = 0;

    opts->formato_csv = false;
    opts->ayuda = false;
    opts->modo = MODO_RNG_RANK;

    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        opts->ayuda = true;
        return true;
    }

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--csv") == 0) {
            opts->formato_csv = true;
        } else if (strcmp(argv[i], "--rng") == 0) {
            if (i + 1 >= argc) {
                *error = "La opción --rng requiere un argumento ('rank' o 'trial').";
                return false;
            }
            const char *val_rng = argv[++i];
            if (strcmp(val_rng, "rank") == 0) {
                opts->modo = MODO_RNG_RANK;
            } else if (strcmp(val_rng, "trial") == 0) {
                opts->modo = MODO_RNG_TRIAL;
            } else {
                *error = "Modo RNG inválido. Debe ser 'rank' o 'trial'.";
                return false;
            }
        } else if (argv[i][0] == '-') {
            *error = "Opción desconocida detectada.";
            return false;
        } else {
            if (posicionales == 0) {
                if (!parsear_u64(argv[i], &opts->tam_grupo) || opts->tam_grupo < 2) {
                    *error = "El tamaño del grupo K debe ser un entero >= 2.";
                    return false;
                }
            } else if (posicionales == 1) {
                if (!parsear_u64(argv[i], &opts->total_ensayos) || opts->total_ensayos == 0) {
                    *error = "El total de ensayos N debe ser un entero > 0.";
                    return false;
                }
            } else if (posicionales == 2) {
                if (!parsear_u64(argv[i], &opts->semilla)) {
                    *error = "La semilla SEED debe ser un entero de 64 bits válido.";
                    return false;
                }
            } else {
                *error = "Demasiados argumentos posicionales proporcionados.";
                return false;
            }
            posicionales++;
        }
    }

    if (posicionales != 3) {
        *error = "Faltan argumentos obligatorios: se requiere <K> <N> <SEED>.";
        return false;
    }

    return true;
}


static long double calcular_probabilidad_teorica(uint64_t k)
{
    long double sin_colision = 1.0L;
    if (k <= 1)
    {
        return 0.0L;
    }
    if (k > DIAS_EN_ANIO)
    {
        return 1.0L;
    }
    for (uint64_t i = 0; i < k; i++)
    {
        sin_colision *= ((long double)DIAS_EN_ANIO - (long double)i) /
                        (long double)DIAS_EN_ANIO;
    }
    return 1.0L - sin_colision;
    
}

static __attribute__((unused)) uint64_t contar_colisiones_locales(uint64_t k, uint64_t ensayos_locales, uint64_t primer_ensayo, uint64_t semilla, modo_rng modo, int rango)
{
    uint64_t colisiones = 0;
    estado_pcg32 rng_rango;

    if (k <= 1)
    {
        return 0;
    }
    if (k > DIAS_EN_ANIO)
    {
        return ensayos_locales;
    }
    if (modo == MODO_RNG_RANK)
    {
        pcg32_inicializar_rango(&rng_rango, semilla, (uint64_t)rango);
    }

    for (uint64_t id_local = 0; id_local < ensayos_locales; ++id_local)
    {
        uint64_t bitset[PALABRAS_BITSET] = { 0, 0, 0, 0, 0, 0};
        estado_pcg32 rng_ensayo;
        estado_pcg32 *rng = &rng_rango;

        if (modo == MODO_RNG_TRIAL)
        {
            pcg32_inicializar_ensayo(&rng_ensayo, semilla, primer_ensayo + id_local);
            rng = &rng_ensayo;
        }

        for (uint64_t persona = 0; persona < k; ++persona)
        {
            uint64_t dia = pcg32_acotado(rng, DIAS_EN_ANIO);
            uint64_t mascara = UINT64_C(1) << (dia & 63U);
            uint64_t *palabra = &bitset[dia >> 6U];

            if ((*palabra & mascara) != 0)
            {
                ++colisiones;
                break;
            }
            *palabra |= mascara;
        }
        
        
    }
    return colisiones;
    
    
    
}

int main(int argc, char **argv)
{
   
    int id_proceso, cantidad_procesos;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id_proceso);
    MPI_Comm_size(MPI_COMM_WORLD, &cantidad_procesos);

    opciones_programa opts = {0};
    const char *error_mensaje = NULL;

    if (!parsear_opciones(argc, argv, &opts, &error_mensaje))
    {
        if (id_proceso == 0)
        {
            fprintf(stderr, "Error de argumentos: %s\n\n", error_mensaje);
            mostrar_ayuda(stderr, argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
        
    }

    if (opts.ayuda)
    {
        if (id_proceso == 0)
        {
            mostrar_ayuda(stdout, argv[0]);
        }
        MPI_Finalize();
        return EXIT_SUCCESS;
    }

    uint64_t ensayos_base = opts.total_ensayos / cantidad_procesos;
    uint64_t residuo = opts.total_ensayos % cantidad_procesos;
    uint64_t ensayos_locales = ensayos_base + (id_proceso < (int)residuo ? 1 : 0);

    uint64_t primer_ensayo = 0;
    for (int p = 0; p < id_proceso; p++)
    {
        primer_ensayo += ensayos_base + (p < (int)residuo ? 1 : 0);
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    double tiempo_inicio = MPI_Wtime();

    uint64_t colisiones_locales = contar_colisiones_locales(opts.tam_grupo, ensayos_locales, primer_ensayo, opts.semilla, opts.modo, id_proceso );

    
 

    uint64_t colisiones_globales = 0;
    MPI_Reduce(&colisiones_locales, &colisiones_globales, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    double tiempo_fin = MPI_Wtime();
    double tiempo_paralelo = tiempo_fin - tiempo_inicio;

    if (id_proceso == 0)
    {
        long double simulado = (long double)colisiones_globales / (long double)opts.total_ensayos;
        long double teorico = calcular_probabilidad_teorica(opts.tam_grupo);
        long double error_abs = fabsl(simulado - teorico);

        if (opts.formato_csv)
        {
            printf("%" PRIu64 ",%" PRIu64 ",%d,%" PRIu64 ",%" PRIu64 ",%.12Lf,%.12Lf,%.12Lf,%.9f\n",
            opts.tam_grupo, opts.total_ensayos, cantidad_procesos, opts.semilla, colisiones_globales,
            simulado, teorico, error_abs, tiempo_paralelo);
        }else{
            printf("K=%" PRIu64 "\nN=%" PRIu64 "\nprocesos=%d\nsemilla=%" PRIu64
            "\nrng=%s\ncolisiones=%" PRIu64 "\np_sim=%.12Lf\n"
            "p_teorico=%.12Lf\nerror_abs=%.12Lf\ntiempo=%.9f\n",
        opts.tam_grupo, opts.total_ensayos, cantidad_procesos, opts.semilla,
            opts.modo == MODO_RNG_TRIAL ? "trial" : "rank", colisiones_globales,
        simulado, teorico, error_abs, tiempo_paralelo);
        }
        
    }
    
    
    
    
    

    MPI_Finalize();
    return EXIT_SUCCESS;
}