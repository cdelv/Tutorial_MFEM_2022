#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
//Ecuación de Laplace
// -∇²u = f

// Declaración de Funciones
double Condicion_Frontera(const Vector &X);
double forzante(const Vector &X);
double difusividad(const Vector &X);

int main(int argc, char *argv[])
{
    // 1. Inicializar MPI y HYPRE.
    Mpi::Init(argc, argv);
    Hypre::Init();
    tic_toc.Clear(); // Para medir tiempo
    tic_toc.Start();

    // 2. Opciones.
    double Lx = 1.0;
    int orden = 2;
    int refinamientos_seriales = 2;
    int refinamientos_paralelos = 2;

    // 3. Crear la malla (mesh). 
    // La creo dentro de un scope para eliminar el mesh serial de la memoria 
    // porque después no lo voy a necesitar. Un mesh bien refinado puede usar mucho espacio.
    ParMesh pmesh;
    {
        // Crear Mesh Serial.
        Mesh mesh = Mesh::MakeCartesian2D(10, 10, Element::TRIANGLE, false, Lx, Lx, true);

        // Refinar Mesh Serial.
        for (int i = 0; i < refinamientos_seriales; ++i)
            mesh.UniformRefinement();

        // Crear Mesh Paralelo.
        pmesh = ParMesh(MPI_COMM_WORLD, mesh);
    }

    // 4. Refinar Mesh Paralelo.
    for (int ii = 0; ii < refinamientos_paralelos; ii++)
        pmesh.UniformRefinement();

    // 5. Definir el espacio de elementos finitos. Usamos H1.
    H1_FECollection fec(orden, pmesh.Dimension());
    ParFiniteElementSpace fespace(&pmesh, &fec);
    HYPRE_BigInt n_elementos = fespace.GlobalTrueVSize();
    if (Mpi::Root())
       std::cout << "Número de Elementos: " << n_elementos << std::endl;

    // 6. Extraer lista de DOFs de frontera.
    // Definir condiciones de fronteras.
    /****
     *                         2
     *         /-------------------------------\
     *         |                               |
     *         |                               |
     *         |                               | 
     *         |                               |
     *       3 |                               | 1
     *         |                               |
     *         |                               |
     *         |                               | 
     *         |                               |
     *         \-------------------------------/
     *                          0
     *
     ****/
    Array<int> ess_bdr(pmesh.bdr_attributes.Max());
    Array<int> ess_tdof;
    ess_bdr = 0;
    ess_bdr [0] = 1;   ess_bdr [1] = 1;   
    ess_bdr [2] = 1;   ess_bdr [3] = 1;   
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof);

    // 7. Definir la solución u como un grid function del espacio de elementos finitos. 
    ParGridFunction u(&fespace);
    ConstantCoefficient zero(0.);
    FunctionCoefficient bdr(Condicion_Frontera);
    u.ProjectCoefficient(zero);
    u.ProjectBdrCoefficient(bdr, ess_bdr);

    // 8. Crear la forma lineal b(.) que corresponde al lado derecho de la ecuación.
    FunctionCoefficient f(forzante);
    ParLinearForm b(&fespace);
    b.AddDomainIntegrator(new DomainLFIntegrator(f));
    b.Assemble();

    // 9. Crear la forma bilineal a(.,.) que corresponde al operador -∇².
    FunctionCoefficient k(difusividad);
    ParBilinearForm a(&fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(k));
    a.Assemble();

    // 10. Crear el sistema lineal A X = B.
    HypreParMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof, u, b, A, X, B);

    // 11. Resolver el sistema usando PCG con el precondicionador de hypre, BoomerAMG.
    HypreBoomerAMG prec(A);
    HyprePCG cg(MPI_COMM_WORLD);    
    cg.SetTol(1e-16);
    cg.SetAbsTol(0.0);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(0);
    prec.SetPrintLevel(0);
    cg.SetPreconditioner(prec);
    cg.SetOperator(A);
    cg.Mult(B, X);

    // 12. Recuperar la solución x.
    a.RecoverFEMSolution(X, b, u);

    // 13. Visualizar la solución con Paraview. 
    ParaViewDataCollection paraview_out = ParaViewDataCollection("Laplace", &pmesh);
    paraview_out.SetLevelsOfDetail(orden);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Solución", &u);
    paraview_out.Save();

    //14. Reportar tiempo de ejecución del programa
    tic_toc.Stop();
    if (Mpi::Root())
       std::cout << "Tiempo de ejecución: " << tic_toc.RealTime() << "s." << std::endl;

    return 0;
}

// Implementación de funciones
double Condicion_Frontera(const Vector &X)
{
    return 0.;
}
double forzante(const Vector &X)
{
    return 1.;
}
double difusividad(const Vector &X)
{
    return 1.;
}
