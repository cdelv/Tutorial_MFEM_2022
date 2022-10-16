#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
//Ecuación de Difusión
// ∂ₜu - ∇²u = f

// Declaración de Funciones
double Condicion_Inicial(const Vector &X);
double Condicion_Frontera(const Vector &X, double t);
double forzante(const Vector &X, double t);
double difusividad(const Vector &X);
void progress_bar(double x, int width);

/** Después de la discretización la ecuación se puede escribir como
 *
 *     M*∂ₜu + K*u = f
 *
 *  u es la solución de la ecuación, M es la matriz de masa
 *  y K = ∇² es el operador de difusión.
 *
 *  La clase ConductionOperator representa el lado derecho de la ODE.
 */
class ConductionOperator : public TimeDependentOperator
{
protected:
    ParFiniteElementSpace &fespace;
    Array<int> ess_bdr;          // Lista con las condiciones de frontera
    Array<int> ess_tdof;         // Lista con las condiciones de frontera
    HypreParMatrix *Mmat = nullptr;
    HypreParMatrix *Kmat = nullptr;
    HypreParMatrix *Kelim = nullptr;
    HypreParMatrix *Tmat = nullptr;    // T = M + dt * K
    HyprePCG M_solver;      // Krylov solver para invertir la matriz de masa M
    HypreBoomerAMG M_prec;  // Precondicionador para la matriz de masa M
    HyprePCG T_solver;      // Solver implícito para T = M + dt K
    HypreBoomerAMG T_prec;  // Precondicionador para el Solver implícito
    Vector B;               // Lado derecho de la ecuación
    mutable Vector z;       // vector auxiliar
public:
    ConductionOperator(ParFiniteElementSpace &f, const Vector &U, Array<int> ess_bdr, Array<int> ess_tdof);
    virtual void Mult(const Vector &U, Vector &du_dt) const;
    virtual void ImplicitSolve(const double dt, const Vector &U, Vector &k);
    virtual ~ConductionOperator();    
};

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
    int refinamientos_paralelos = 1;

    double dt = 0.005;
    double t, tmax = 1;
    int vis_print = 1, NFrames = 500;
    double tprint, print_every=tmax/NFrames;

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
       std::cout << "Numero de Elementos: " << n_elementos << std::endl;

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
    ess_bdr [0] = 0;   ess_bdr [1] = 0;   
    ess_bdr [2] = 0;   ess_bdr [3] = 0;   
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof);

    // 7. Definir la solución u como un grid function del espacio. 
    ParGridFunction u(&fespace);
    FunctionCoefficient u_0(Condicion_Inicial);
    FunctionCoefficient bdr(Condicion_Frontera);
    u.ProjectCoefficient(u_0);
    u.ProjectBdrCoefficient(bdr, ess_bdr);

    // 8. Definir el integrador temporal. Dentro del integrador creamos 
    // las formas lineales y bilineales del sistema. 
    Vector U;
    u.GetTrueDofs(U);
    ConductionOperator oper(fespace, U, ess_bdr, ess_tdof);
    ODESolver *ode_solver = new BackwardEulerSolver;
    ode_solver->Init(oper);

    // 9. Preparar la Visualización de la solución con Paraview. 
    ParaViewDataCollection paraview_out = ParaViewDataCollection("Difusion", &pmesh);
    paraview_out.SetLevelsOfDetail(orden);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Solución", &u);

    // 10. Imprimir estado inicial.
    paraview_out.SetCycle(vis_print);
    paraview_out.SetTime(t);
    paraview_out.Save();
    if(Mpi::Root())
        std::cout << "\n"; // Para que la barra se vea bien.

    // 11. Evolución temporal.
    for (tprint=0, t=0; t<tmax; tprint+=dt)
    {   
        ode_solver->Step(U, t, dt);

        // Imprimir Frame.
        if(tprint > print_every){
            vis_print+=1;
            u.Distribute(U);
            paraview_out.SetCycle(vis_print);
            paraview_out.SetTime(t);
            paraview_out.Save();
        }

        // Barra de progreso.
        if(Mpi::Root())
            progress_bar(t,60);
    }
    if(Mpi::Root())
        std::cout << "\n\n"; // Para que la barra se vea bien.

    // 12. Liberar la memoria.
    delete ode_solver;

    //13. Reportar tiempo de ejecución del programa
    tic_toc.Stop();
    if (Mpi::Root())
       std::cout << "Tiempo de ejecución: " << tic_toc.RealTime() << "s." << std::endl;

    return 0;
}

// Constructor de ConductionOperator
ConductionOperator::ConductionOperator(ParFiniteElementSpace &f, const Vector &U, 
                                       Array<int> ess_bdr, Array<int> ess_tdof) : 
    TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), 
    M_solver(f.GetComm()), T_solver(f.GetComm()), 
    z(height), B(height), 
    ess_bdr(ess_bdr), ess_tdof(ess_tdof)
{   
    // Crear la matriz de masa.
    ParBilinearForm M(&fespace);
    M.AddDomainIntegrator(new MassIntegrator());
    M.Assemble();
    M.EliminateEssentialBC(ess_bdr, Operator::DIAG_ONE);
    M.Finalize();
    Mmat = M.ParallelAssemble();

    // Crear la matriz asociada al operador de difusión.
    ParBilinearForm K(&fespace);
    FunctionCoefficient D(difusividad);
    K.AddDomainIntegrator(new DiffusionIntegrator(D));
    K.Assemble();
    K.Finalize();
    Kmat = K.ParallelAssemble();
    Kelim = Kmat->EliminateRowsCols(ess_tdof); // Para eliminar los tdof después.

    // Crear lado derecho de la ecuación.
    FunctionCoefficient RHS(forzante);
    ParLinearForm b(&fespace);
    b.AddDomainIntegrator(new DomainLFIntegrator(RHS));
    b.Assemble();
    b.ParallelAssemble(B);

    // Imponer condiciones de frontera de Dirichlet.
    Kmat->EliminateBC(*Kelim, ess_tdof, U, B);

    //Configure M solver.
    M_prec.SetPrintLevel(0);
    M_prec.SetOperator(*Mmat);
    M_solver.SetTol(1e-15);
    M_solver.SetAbsTol(0.0);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(*Mmat);

    //Configure T solver.
    T_prec.SetPrintLevel(0);
    T_solver.SetTol(1e-15);
    T_solver.SetAbsTol(0.0);
    T_solver.SetMaxIter(100);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);
}
void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // ∂ₜu = M^{-1}*(f-Ku)
    // Para ∂ₜu, donde K es linealizado usando u del paso anterior.
    Kmat->Mult(u, z);
    z.Neg(); // z = -z
    z+=B;
    M_solver.Mult(z, du_dt);
}
void ConductionOperator::ImplicitSolve(const double dt, const Vector &U, Vector &du_dt)
{
    // Solucionar la ecuación
    // ∂ₜu = M^{-1}*[f - K(u + dt*∂ₜu)]
    // Para ∂ₜu, donde K es linealizado usando u del paso anterior.
    if (Tmat)
        delete Tmat;
    
    Tmat = Add(1.0, *Mmat, dt, *Kmat);
    T_prec.SetOperator(*Tmat);
    T_solver.SetOperator(*Tmat);

    Kmat->Mult(U, z);
    z.Neg();
    z+=B;
    T_solver.Mult(z, du_dt);
}
ConductionOperator::~ConductionOperator(){
    delete Tmat;
    delete Mmat;
    delete Kmat;
    delete Kelim;
}

// Implementación de funciones.
double Condicion_Frontera(const Vector &X, double t){
    return 0.;
}
double forzante(const Vector &X, double t){
    return 1.;
}   
double difusividad(const Vector &X){
    return 1.;
}
double Condicion_Inicial(const Vector &X){
    return 0.;
}
void progress_bar(double x, int width){
    std::cout << "[";
    int pos = width * x;
    for (int i = 0; i < width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(x * 100.0) << " %\r";
    std::cout.flush();
}
