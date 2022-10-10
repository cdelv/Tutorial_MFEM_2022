#Para instalación de la fuente:
#MY_MFEM_INSTALL_DIR = mfem/build
#HYPRE_INC = -I$(MY_MFEM_INSTALL_DIR)/../../hypre/src/hypre/include

#Para instalación con Spack:
MY_MFEM_INSTALL_DIR = spack/opt/spack/linux-pop22-icelake/gcc-11.2.0/mfem-4.4.0-7asre5epqp3dfcqapvavyggiqiuhdb55

CONFIG_MK=$(MY_MFEM_INSTALL_DIR)/share/mfem/config.mk 

include $(CONFIG_MK)

###

CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS) $(HYPRE_INC)
FILE = Laplace.cpp

all:
	@$(CXX) $(FLAGS) $(FILE) $(MFEM_LIBS)
