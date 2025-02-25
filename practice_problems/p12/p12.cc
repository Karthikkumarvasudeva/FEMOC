#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

using namespace dealii;

class BiharmonicProblem
{
public:
    BiharmonicProblem();
    void run();

private:
    void setup_system();
    void assemble_system();
    void solve();
    void output_results() const;

    Triangulation<2> triangulation;
    FESystem<2> fe;
    DoFHandler<2> dof_handler;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> solution;
    Vector<double> system_rhs;
};

BiharmonicProblem::BiharmonicProblem()
    : fe(FE_Q<2>(2), 1), dof_handler(triangulation)
{}

void BiharmonicProblem::setup_system()
{
    dof_handler.distribute_dofs(fe);
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}

void BiharmonicProblem::assemble_system()
{
    QGauss<2> quadrature_formula(fe.degree + 1);
    FEValues<2> fe_values(fe, quadrature_formula,
                          update_values | update_gradients |
                              update_hessians | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        cell_matrix = 0;
        cell_rhs = 0;
        fe_values.reinit(cell);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                for (unsigned int q_index = 0; q_index < quadrature_formula.size(); ++q_index)
                    cell_matrix(i, j) += fe_values.shape_hessian(i, q_index) *
                                         fe_values.shape_hessian(j, q_index) *
                                         fe_values.JxW(q_index);

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
    }
}

void BiharmonicProblem::solve()
{
    SolverControl solver_control(1000, 1e-12);
    SolverCG<> solver(solver_control);
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
    solver.solve(system_matrix, solution, system_rhs, preconditioner);
}

void BiharmonicProblem::output_results() const
{
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();
    std::ofstream output("solution.vtk");
    data_out.write_vtk(output);
}

void BiharmonicProblem::run()
{
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(5);
    setup_system();
    assemble_system();
    solve();
    output_results();
}

int main()
{
    BiharmonicProblem biharmonic_problem;
    biharmonic_problem.run();
    return 0;
}
