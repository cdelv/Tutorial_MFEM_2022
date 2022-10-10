# Tutorial_MFEM_2022

## Instalar MFEM

Hay 2 formas de instalar MFEM. La primera es con Spack (https://spack.io/) y la segunda es desde la fuente. Para la instalación deben tener MPI instalado o seguir las instrucciones en https://mfem.org/building/ para instalar la versión serial. Nosotros vamos a usar la versión paralela. Para usar spack deben tener python.

### Spack:

escriban en una terminal 

```
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
```

para instalar MFEM pueden hacer

```
cd spack/bin
./spack install mfem
```

o alternativamente

```
source spack/share/spack/setup-env.sh
spack install mfem
```

En MAC hay que cambiar la extensión por .csh y para arch linux (fish shell) .fish en vez de .sh

Desde la fuente pueden usar el script de bash adjunto o usar los comandos del script como guía para hacer su propia instalación.  Deben tener MPI, Cmake y Make instalados. También pueden seguir las intstrucciones de la página de MFEM.

```
bash install_mfem.sh
```

Todo esto es para linux, sé que se puede en Windows. La instalación en MAC es parecida a la de linux, para procesadores M1 y M2 hay que hacer algunas modificaciones.

## Correr un programa 

Linkear librerias en C++ puede ser complicado, así que pueden usar el Makefile que preparamos para eso. Deben entrar al Makefile y que colocar el Path de la carpeta de la instalación de mfem y el Makefile va a hacer todo el trabajo. En la misma carpeta debe estar el programa que quieren correr. Dentro del Makefile deben cambiar el nombre del programa que desean correr. 

Escribir make va a compilar el programa, para correrlo pueden hacer 

```
./a.out 
```

o si quieren correr en paralelo 

```
mpirun -np 4 ./a.out
```

El 4 significa que quiero usar 4 procesadores. 

Despues de correr el programa se va a crear una carpeta y dentro está la solución de la ecuación. La pueden ver con Paraview.