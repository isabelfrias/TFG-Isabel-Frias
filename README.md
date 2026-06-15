# TFG "Implementación eficiente de la FFT para una CGRA integrada de bajo consumo en la plataforma X-HEEP"
Repositorio correspondiente al Trabajo de Fin de Grado de Isabel L. Frías Castillo, realizado para la asignatura "Trabajo de Fin de Grado" (Curso 2025/2026, Convocatoria de junio 2026).

## Descripción del proyecto
Este proyecto es centra en la optimización e implementación del algoritmo FFT (*Fast Fourier Transform*) sobre una arquitectura CGRA (*Coarse-Grained Reconfigurable Architechture*), integrada en el ecosistema X-HEEP. El objetivo principal es lograr una ejecución altamente eficiente y de ultrabajo consumo.

## Requisitos Previos
Para el correcto funcionamiento del código, es necesario contar con:
- **Docker Desktop** (con el entorno *HEEPsilon* configurado)
- **Repositorio base HEEPsilon** (Fork del grupo UMALabRISCV)
- **WSL** o **Linux** como sistema operativo

## Ejecución de las implementaciones 
1. Clonar el repositorio **HEEPsilon** del fork del grupo *UMALabRISCV*: https://github.com/UMALabRISCV/HEEPsilon.git. y realizar la instalación de Docker siguiendo los pasos del README.
2. Incluir las carpetas de este repositorio dentro de ```sw/applications```.
### Guía de ejecución
Los comandos para ejecutar son los siguientes comenzando desde la raíz del repositorio:
```bash
make clean-app
make app PROJECT=NOMBRE_DEL_PROYECTO/CARPETA TARGET=sim
make verilator-sim
cd build/eslepfl_systems_heepsilon_0/sim-verilator
./Vtestharness +firmware=../../../sw/build/main.hex
cat uart0.log
```

## Estructura de las carpetas
###  `/CGRA` - Implementación Base en CGRA
Versión sin optimizaciones del algoritmo FFT ejecutado en la CGRA.

###  `/CGRA_O1` - Primera Optimización
Primera versión optimizada que mejora el rendimiento del CGRA mediante cambios en la configuración del bus.
- **Mejora principal**: El CGRA computa todos los butterflies en una etapa en un único lanzamiento de kernel, incrementando la granularidad y reduciendo el overhead de reconfiguración.
- **Optimización adicional**: Bus configurado en formato `ntoM`, permitiendo accesos concurrentes a memoria desde el CGRA, reduciendo bloqueos.

###  `/CGRA_O2` - Segunda Optimización
Segunda versión optimizada con desenrollado de bucles en el kernel del CGRA.
- **Mejora principal**: Loop unrolling de las etapas FFT en el kernel del CGRA, permitiendo computar múltiples etapas en un único lanzamiento de kernel.
- **Parámetro configurable**: `CGRA_STAGES` define el factor de despliegue (máximo: log₂(FFT_SIZE)).
- **Impacto**: Reducción significativa del overhead por lanzamientos de kernel y reconfiguración de punteros.
- 
###  `/CPU` - Implementación en CPU
Implementación de referencia del algoritmo FFT ejecutado completamente en la CPU para comparación de rendimiento. Incluye el archivo modificado ```CMakeLists.txt``` con la opción añadida para la optimización con flags (CPU_O1)

## Archivos Comunes
Los siguientes archivos aparecen en las carpetas CGRA, CGRA_O1, CGRA_O2 y CPU.

### **main.c**
Archivo principal de la aplicación. Contiene:
- Inicialización del CGRA/CPU
- Configuración de interrupciones y contadores de rendimiento
- Bit-reversal del input para FFT Decimation-In-Time (DIT)
- Bucle principal de cálculo de FFT (varía según optimización)
- Verificación de resultados contra valores esperados
- Medición y reporte de métricas de rendimiento (ciclos, latencia, stalls)

### **fft_data.h**
Archivo de datos que contiene:
- `input_signal[]`: Vector de entrada con valores de señal complejos (parte real e imaginaria intercaladas)
- `exp_output_real[]` y `exp_output_imag[]`: Valores esperados de salida para verificación

### **fft_factors_*_32b_int.h**
Archivos de precálculo de factores twiddle para diferentes tamaños de FFT:
- `fft_factors_8_32b_int.h`: Factores para FFT de 8 puntos
- `fft_factors_256_32b_int.h`: Factores para FFT de 256 puntos
- `fft_factors_512_32b_int.h`: Factores para FFT de 512 puntos
- `fft_factors_1024_32b_int.h`: Factores para FFT de 1024 puntos
- `fft_factors_2048_32b_int.h`: Factores para FFT de 2048 puntos

Estos archivos contienen arrays precalculados con valores reales e imaginarios (`f_real[]` y `f_imag[]`) usados durante el cálculo de cada etapa de la FFT.

### **fxp.h / fxp.c**
Implementación de aritmética de punto fijo (Fixed-Point eXtended):
- Define el tipo `fxp` utilizado para los cálculos numéricos
- Implementa operaciones aritméticas especializadas para punto fijo
- Crítico para la precisión de los cálculos FFT en aritmética de punto fijo

### **defines.h**
Archivo de definiciones de configuración:
- Define macros de configuración del tamaño FFT, valores de factores de escala y otras constantes del algoritmo
- Permite compilación flexible para diferentes tamaños de entrada

### **cgra_bitstream.h**
Archivo específico del CGRA que contiene:
- Definición del bitstream de configuración del CGRA
- Inicialización de memoria de contexto (context memory) y kernel memory
- Configuración necesaria para programar el CGRA

**Archivos específicos de CGRA/CGRA_O1/CGRA_O2**:
- `cgra_bitstream.h`: Contiene bitstreams diferentes para cada versión según su configuración de kernel optimizada

### **instructions.csv**
Secuencia de instrucciones CGRA en formato CSV:
- Describe el kernel ejecutado por el CGRA
- Formato: columnas representan los 4 elementos procesadores del CGRA
- Estructura: Etapas de instrucción (ciclos de reloj) en filas

