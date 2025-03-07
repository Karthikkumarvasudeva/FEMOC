#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/lac/vector.h>

#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>


#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

class CoupledSystem {
public:
    //CoupledSystem(const unsigned int degree, const unsigned int refinement_steps, const double nu)
    //    : degree(degree), refinement_steps(refinement_steps), nu(nu),
    //      fe(degree), dof_handler(triangulation) {}


    //CoupledSystem(const unsigned int degree, const unsigned int refinement_steps, const double nu)
    //: degree(degree),
    //  nu(nu),
    //  refinement_steps(refinement_steps),
    //  fe(degree),
    //  dof_handler(triangulation),
    //  data_out() {}


    CoupledSystem(const unsigned int degree, const unsigned int refinement_steps, const double nu)
    : triangulation(),        // Initialize first
      fe(degree),             // FE_Q depends on degree
      dof_handler(triangulation), // DoFHandler depends on triangulation
      data_out(),
      mass_matrix(),
      b_matrix(),
      c1_matrix(),
      c2_matrix(),
      s1_n(),
      s2_n(),
      w_n1(),
      u_n1(),
      nu(nu),                // nu comes before degree
      degree(degree),        // Match declaration order
      refinement_steps(refinement_steps),
      preconditioner() {}    // Last member



    void run();

private:
    void setup_system();
    void assemble_m_b_matrices();
    void solve_wu_system(Vector<double>& w, Vector<double>& u, const Vector<double>& s1, const Vector<double>& s2);
    void solve_si_system(Vector<double>& si, const SparseMatrix<double>& ci, const Vector<double>& u);
    void output_results(const unsigned int cycle) const;

    Triangulation<2> triangulation;
    FE_Q<1> fe;  //quadratic elements
    DoFHandler<2> dof_handler;
    DataOut<2> data_out;

    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> b_matrix;
    SparseMatrix<double> c1_matrix;
    SparseMatrix<double> c2_matrix;

    Vector<double> s1_n, s2_n, w_n1, u_n1;
    double nu;
    const unsigned int degree;
    const unsigned int refinement_steps;
        // Reusable preconditioner
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
};


/*void CoupledSystem::setup_system() {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern mass_dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, mass_dsp);
    mass_matrix.reinit(mass_dsp);
    b_matrix.reinit(dof_handler.n_dofs(), dof_handler.n_dofs());

    dealii::SparsityPattern sparsity_pattern;

// Generate SparsityPattern using DynamicSparsityPattern
    dealii::DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

// Reinitialize the matrix with the sparsity pattern
    c1_matrix.reinit(sparsity_pattern);

   // // // c1_matrix.reinit(dof_handler.n_dofs(), dof_handler.n_dofs());

// Reinitialize the matrix with the sparsity pattern
    c2_matrix.reinit(sparsity_pattern);
    // // // c2_matrix.reinit(dof_handler.n_dofs(), dof_handler.n_dofs());
}*/

/*void CoupledSystem::setup_system() {
    dof_handler.distribute_dofs(fe);
   dealii::
    // Create DynamicSparsityPattern
    DynamicSparsityPattern mass_dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, mass_dsp);

    // Reinitialize SparsityPattern using the DynamicSparsityPattern
    SparsityPattern mass_sp;
    mass_sp.copy_from(mass_dsp);
    mass_matrix.reinit(mass_sp);

    // Same for b_matrix, c1_matrix, c2_matrix
    DynamicSparsityPattern b_dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, b_dsp);
    SparsityPattern b_sp;
    b_sp.copy_from(b_dsp);
    b_matrix.reinit(b_sp);

    DynamicSparsityPattern c1_dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, c1_dsp);
    SparsityPattern c1_sp;
    c1_sp.copy_from(c1_dsp);
    c1_matrix.reinit(c1_sp);

    DynamicSparsityPattern c2_dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, c2_dsp);
    SparsityPattern c2_sp;
    c2_sp.copy_from(c2_dsp);
    c2_matrix.reinit(c2_sp);
}*/


void CoupledSystem::setup_system() {
    dof_handler.distribute_dofs(fe);

    // Create DynamicSparsityPattern
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    // Reinitialize SparsityPattern using the DynamicSparsityPattern
    SparsityPattern sp;
    sp.copy_from(dsp);
    mass_matrix.reinit(sp);

    // Same for b_matrix, c1_matrix, c2_matrix
    b_matrix.reinit(sp);
    c1_matrix.reinit(sp);
    c2_matrix.reinit(sp);
}

void CoupledSystem::assemble_m_b_matrices() {
        const unsigned int dofs_per_cell = fe.dofs_per_cell;

    QGauss<2> quadrature_formula(degree + 1);
    FEValues<2> fe_values(fe, quadrature_formula,
                                  update_values | update_gradients |
                                      update_quadrature_points | update_JxW_values);


    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators()) {
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                for (unsigned int q_index = 0; q_index < quadrature_formula.size(); ++q_index) {
                    mass_matrix.add(local_dof_indices[i], local_dof_indices[j],
                                    fe_values.shape_value(i, q_index) *
                                        fe_values.shape_value(j, q_index) *
                                        fe_values.JxW(q_index));
                    b_matrix.add(local_dof_indices[i], local_dof_indices[j],
                                   (fe_values.shape_grad(i, q_index)[0] *
                                        fe_values.shape_grad(j, q_index)[0] +
                                    fe_values.shape_grad(i, q_index)[1] *
                                        fe_values.shape_grad(j, q_index)[1]) *
                                        fe_values.JxW(q_index));
                    c1_matrix.add(local_dof_indices[i], local_dof_indices[j],
                                   0.5 * (fe_values.shape_grad(i, q_index)[0] *
                                              fe_values.shape_grad(j, q_index)[0] -
                                          fe_values.shape_grad(i, q_index)[1] *
                                              fe_values.shape_grad(j, q_index)[1]) *
                                          fe_values.JxW(q_index));
                    c2_matrix.add(local_dof_indices[i], local_dof_indices[j],
                                   0.5 * (fe_values.shape_grad(i, q_index)[1] *
                                              fe_values.shape_grad(j, q_index)[0] -fe_values.shape_grad(i, q_index)[0] *
                                              fe_values.shape_grad(j, q_index)[1]) *
                                          fe_values.JxW(q_index));
                }
            }
        }
    }
        // Now initialize the preconditioner with the mass matrix.
            // Ensure mass_matrix is non-empty before initializing preconditioner
    if (mass_matrix.m() > 0) {
        preconditioner.initialize(mass_matrix, 1.2);
    }
  // // //  preconditioner.initialize(mass_matrix, 1.2);
}

void CoupledSystem::run() {
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(refinement_steps);

    setup_system();
    assemble_m_b_matrices();


    s1_n.reinit(dof_handler.n_dofs());
    s2_n.reinit(dof_handler.n_dofs());
    w_n1.reinit(dof_handler.n_dofs());
    u_n1.reinit(dof_handler.n_dofs());

    if (std::abs(nu - 1.0) < 1e-10) {
        // ν = 1: Block-triangular solve
        solve_wu_system(w_n1, u_n1, s1_n, s2_n);
        solve_si_system(s1_n, c1_matrix, u_n1);
        solve_si_system(s2_n, c2_matrix, u_n1);
        output_results(0);
    } else {
        // ν ≠ 1: Relaxation algorithm
        for (unsigned int cycle = 0; cycle < 10; ++cycle) {
            solve_wu_system(w_n1, u_n1, s1_n, s2_n);
            solve_si_system(s1_n, c1_matrix, u_n1);
            solve_si_system(s2_n, c2_matrix, u_n1);
            output_results(cycle);
        }
    }
}



void CoupledSystem::solve_wu_system(Vector<double>& w, Vector<double>& u, const Vector<double>& s1, const Vector<double>& s2) {

    // Define the sparsity pattern
    DynamicSparsityPattern dsp(2 * dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    SparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);

    // Create the system matrix with the sparsity pattern
   // // //  SparseMatrix<double> system_matrix(sparsity_pattern);

     SparseMatrix<double> system_matrix;
    system_matrix.reinit(sparsity_pattern);  // Now properly initialized


    // Create system_rhs and solution vectors with the appropriate size
    Vector<double> system_rhs(2 * dof_handler.n_dofs());
    system_rhs = 0.0;  // Initialize to zero
    Vector<double> solution(2 * dof_handler.n_dofs());
    solution = 0.0;  // Initialize to zero



    // Fill the system_matrix
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) {
        for (unsigned int j = 0; j < dof_handler.n_dofs(); ++j) {
            system_matrix.set(i, j, mass_matrix.el(i, j));
            system_matrix.set(i, j + dof_handler.n_dofs(), b_matrix.el(j, i));
            system_matrix.set(i + dof_handler.n_dofs(), j, b_matrix.el(i, j));
            system_matrix.set(i + dof_handler.n_dofs(), j + dof_handler.n_dofs(), 0.0);
        }
    }

    // Initialize FEValues
    QGauss<2> quadrature_formula(degree + 1);
    FEValues<2> fe_values(fe, quadrature_formula,
                          update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Assemble the system_rhs
    for (const auto &cell : dof_handler.active_cell_iterators()) {
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int q_index = 0; q_index < quadrature_formula.size(); ++q_index) {
                system_rhs(local_dof_indices[i] + dof_handler.n_dofs()) +=
                    -nu * std::sin(M_PI * fe_values.quadrature_point(q_index)[0]) *
                    std::sin(M_PI * fe_values.quadrature_point(q_index)[1]) *
                    fe_values.shape_value(i, q_index) * fe_values.JxW(q_index);

                system_rhs(local_dof_indices[i] + dof_handler.n_dofs()) -=
                    nu * (c1_matrix.el(local_dof_indices[i], 0) * s1(0) +
                          c1_matrix.el(local_dof_indices[i], 1) * s1(1) +
                          c1_matrix.el(local_dof_indices[i], 2) * s1(2) +
                          c2_matrix.el(local_dof_indices[i], 0) * s2(0) +
                          c2_matrix.el(local_dof_indices[i], 1) * s2(1) +
                          c2_matrix.el(local_dof_indices[i], 2) * s2(2));
            }
        }
    }

    // Set up the solver and solve the system
    SolverControl solver_control(1000, 1e-10);
    SolverCG<> solver_cg(solver_control);

    solver_cg.solve(system_matrix, solution, system_rhs, preconditioner);

 // //       solver_cg.solve(system_matrix, solution, system_rhs, preconditionIdentity());


    // Post-process the solution
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) {
        w(i) = solution(i);
        u(i) = solution(i + dof_handler.n_dofs());
    }
}


void CoupledSystem::solve_si_system(Vector<double>& si, const SparseMatrix<double>& ci, const Vector<double>& u) {
    Vector<double> rhs(dof_handler.n_dofs());
    ci.vmult(rhs, u); ////////////////////////////////////////////
    rhs *= -1.0;

    SolverControl solver_control(1000, 1e-10);
    SolverCG<> solver_cg(solver_control);

    solver_cg.solve(mass_matrix, si, rhs, preconditioner);
}

void CoupledSystem::output_results(const unsigned int cycle) const {
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(w_n1, "W");
    data_out.add_data_vector(u_n1, "U");
    data_out.add_data_vector(s1_n, "S1");
    data_out.add_data_vector(s2_n, "S2");
    data_out.build_patches();

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
    data_out.write_vtk(output);
}

int main() {
    try {
        CoupledSystem coupled_system(1, 4, 0.5); // degree, refinement, nu
        coupled_system.run();
    } catch (std::exception &exc) {
        std::cerr << std::endl
                  << "----------------------------------------------------"
                  << std::endl
                  << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    } catch (...) {
        std::cerr << std::endl
                  << "----------------------------------------------------"
                  << std::endl
                  << "Exception of unknown type!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }
    return 0;
}
