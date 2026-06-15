# TFG "Implementación eficiente de la FFT para una CGRAintegrada de bajo consumo en la plataforma X-HEEP"
Repositorio correspondiente al Trabajo de Fin de Grado (TFG) de Isabel L. Frías Castillo, realizado para la asignatura "Trabajo de Fin de Grado" (Curso 2025/2026, Convocatoria de junio 2026).

## Descripción del proyecto
Este proyecto es centra en la optimización e implementación del algoritmo FFT (*Fast Fourier Transform*) sobre una arquitectura CGRA (*Coarse-Grained Reconfigurable Architechture*), integrada en el ecosistema X-HEEP. El objetivo principal es lograr una ejecución altamente eficiente y de ultrabajo consumo.

## Requisitos Previos
Para el correcto funcionamiento del código, es necesario contar con:
- **Docker Desktop** (con el entorno *HEEPsilon* configurado)
- **Repositorio base HEEPsilon** (Fork del grupo UMALabRISCV)
- **WSL** o **Linux** como sistema operativo

## Ejecución de las implementaciones 
1. Clonar el repositorio **HEEPsilon** del fork del grupo *UMALabRISCV*: https://github.com/UMALabRISCV/HEEPsilon.git. y realizar la instalación de Docker siguiendo los pasos del README.
2. Incluir las carpetas de este repositorio dentro de ```bash sw/applications ```.
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
