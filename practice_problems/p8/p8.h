#ifndef FEM_H
#define FEM_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <map>
#include <vector>
#include <iostream>
#include <cmath>

using namespace dealii;

template <int dim>
class FEM {
public:
    FEM(unsigned int order, unsigned int problem);
    ~FEM();

    void generate_mesh(unsigned int numberOfElements);
    void define_boundary_conds();
    void setup_system();
    void assemble_system();
    void solve();
    void output_results();
    double l2norm_of_error();

private:
    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> K;
    Vector<double> D, F;

    std::vector<double> nodeLocation;
    std::map<unsigned int, double> boundary_values;

    unsigned int quadRule;
    std::vector<double> quad_points;
    std::vector<double> quad_weight;

    unsigned int basisFunctionOrder;
    unsigned int prob;
    double L, g1, g2;

    double xi_at_node(unsigned int dealNode);
    double basis_function(unsigned int node, double xi);
    double basis_gradient(unsigned int node, double xi);
};

#endif  // FEM_H
