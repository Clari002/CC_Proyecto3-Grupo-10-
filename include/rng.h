#ifndef CUMPLE_RNG_H
#define CUMPLE_RNG_H

#include <stdint.h>

typedef struct {
    uint64_t estado;    
    uint64_t incremento; 
} estado_pcg32;


static inline uint64_t splitmix64_siguiente(uint64_t *estado_ptr)
{
    uint64_t z = (*estado_ptr += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31U);
}

static inline uint32_t pcg32_siguiente(estado_pcg32 *rng)
{
    uint64_t estado_anterior = rng->estado;
    uint32_t desplazamiento_xor = (uint32_t)(((estado_anterior >> 18U) ^ estado_anterior) >> 27U);
    uint32_t rotacion = (uint32_t)(estado_anterior >> 59U);
    

    rng->estado = estado_anterior * UINT64_C(6364136223846793005) + rng->incremento;
    
    return (desplazamiento_xor >> rotacion) | (desplazamiento_xor << ((-rotacion) & 31U));
}


static inline void pcg32_inicializar(estado_pcg32 *rng, uint64_t estado_inicial, uint64_t flujo)
{
    rng->estado = 0;
    rng->incremento = (flujo << 1U) | 1U;
    (void)pcg32_siguiente(rng);
    rng->estado += estado_inicial;
    (void)pcg32_siguiente(rng);
}

static inline void pcg32_inicializar_rango(estado_pcg32 *rng, uint64_t semilla, uint64_t id_proceso)
{
    uint64_t mezclador = semilla ^ (id_proceso * UINT64_C(0xd2b74407b1ce6e93));
    uint64_t estado_inicial = splitmix64_siguiente(&mezclador);
    uint64_t flujo = splitmix64_siguiente(&mezclador);
    pcg32_inicializar(rng, estado_inicial, flujo);
}


static inline void pcg32_inicializar_ensayo(estado_pcg32 *rng, uint64_t semilla, uint64_t id_ensayo)
{
    uint64_t mezclador = semilla ^ (id_ensayo * UINT64_C(0x9e3779b97f4a7c15));
    uint64_t estado_inicial = splitmix64_siguiente(&mezclador);
    uint64_t flujo = splitmix64_siguiente(&mezclador);
    pcg32_inicializar(rng, estado_inicial, flujo);
}


static inline uint32_t pcg32_acotado(estado_pcg32 *rng, uint32_t cota)
{
    uint32_t umbral = (uint32_t)(-cota) % cota;
    for (;;) {
        uint32_t valor = pcg32_siguiente(rng);
        if (valor >= umbral)
        {
            return valor % cota;
        } 
    }
}

#endif