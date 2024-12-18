#include "p8.h"

#include "p8.h"

int main() {
    // You can call functions from your FEM class here
    FEM<2> fem(2, 1);  // Example: create FEM object with 2D and problem 1
    fem.generate_mesh(10);
    fem.setup_system();
    fem.assemble_system();
    fem.solve();
    fem.output_results();

    double l2norm = fem.l2norm_of_error();  // Declare L2 norm function
    std::cout << "L2 norm of the error: " << l2norm << std::endl;



    return 0;  // Return success
}


template <int dim>
FEM<dim>::FEM(unsigned int order, unsigned int problem)
    : fe(order), dof_handler(triangulation), basisFunctionOrder(order), prob(problem) {
    if (prob != 1 && prob != 2) {
        std::cerr << "Error: Problem number should be 1 or 2." << std::endl;
        exit(EXIT_FAILURE);
    }
}

template <int dim>
FEM<dim>::~FEM() {
    dof_handler.clear();
}
template <int dim>
void FEM<dim>::generate_mesh(unsigned int numberOfElements) {
    L = 1.0;
    Point<dim> min(0.0), max(L);
    // Ensure two subdivisions for 2D mesh
    std::vector<unsigned int> n_subdivisions(2, numberOfElements);

    // Debug print
    std::cout << "Number of elements: " << numberOfElements << std::endl;
    std::cout << "Subdivisions: " << n_subdivisions[0] << ", " << n_subdivisions[1] << std::endl;

    GridGenerator::subdivided_hyper_rectangle(triangulation, n_subdivisions, min, max);
    std::cout << "Mesh generated with " << triangulation.n_active_cells() << " cells." << std::endl;
}


template <int dim>
void FEM<dim>::define_boundary_conds() {
    for (unsigned int globalNode = 0; globalNode < nodeLocation.size(); globalNode++) {
        if (nodeLocation[globalNode] == 0) {
            boundary_values[globalNode] = g1;
        }
        if (nodeLocation[globalNode] == L && prob == 1) {
            boundary_values[globalNode] = g2;
        }
    }
}

template <int dim>
void FEM<dim>::setup_system() {
    g1 = 0;
    g2 = 0;

    dof_handler.distribute_dofs(fe);

    MappingQ1<dim, dim> mapping;
    std::vector<Point<dim, double>> dof_coords(dof_handler.n_dofs());
    nodeLocation.resize(dof_handler.n_dofs());
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, dof_coords);

    for (unsigned int i = 0; i < dof_coords.size(); i++) {
        nodeLocation[i] = dof_coords[i][0];
    }

    define_boundary_conds();

    sparsity_pattern.reinit(dof_handler.n_dofs(), dof_handler.n_dofs(), dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
    sparsity_pattern.compress();

    K.reinit(sparsity_pattern);
    F.reinit(dof_handler.n_dofs());
    D.reinit(dof_handler.n_dofs());

    quadRule = 2;
    quad_points = {-std::sqrt(1.0 / 3.0), std::sqrt(1.0 / 3.0)};
    quad_weight = {1.0, 1.0};

    std::cout << "Number of active elements: " << triangulation.n_active_cells() << std::endl;
    std::cout << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;
}

template <int dim>
void FEM<dim>::assemble_system() {
    K = 0;
    F = 0;

    const unsigned int dofs_per_elem = fe.dofs_per_cell;
    FullMatrix<double> Klocal(dofs_per_elem, dofs_per_elem);
    Vector<double> Flocal(dofs_per_elem);
    std::vector<unsigned int> local_dof_indices(dofs_per_elem);

    typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(), endc = dof_handler.end();

    for (; elem != endc; ++elem) {
        elem->get_dof_indices(local_dof_indices);

        // Remove or comment this line if unused
        // double h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]];

        Flocal = 0.0;
        Klocal = 0.0;
        for (unsigned int A = 0; A < dofs_per_elem; A++) {
            for (unsigned int B = 0; B < dofs_per_elem; B++) {
                for (unsigned int q = 0; q < quadRule; q++) {
                    // Add stiffness matrix assembly logic here
                }
            }
        }

        for (unsigned int A = 0; A < dofs_per_elem; A++) {
            for (unsigned int B = 0; B < dofs_per_elem; B++) {
                K.add(local_dof_indices[A], local_dof_indices[B], Klocal(A, B));
            }
            F(local_dof_indices[A]) += Flocal(A);
        }
    }

    MatrixTools::apply_boundary_values(boundary_values, K, D, F, false);
}

template <int dim>
void FEM<dim>::solve() {
    SparseDirectUMFPACK A;
    A.initialize(K);
    A.vmult(D, F);
}

template <int dim>
void FEM<dim>::output_results() {
    std::ofstream output("solution.vtk");
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(D, "solution");
    data_out.build_patches();
    data_out.write_vtk(output);
    std::cout << "Solution written to solution.vtk" << std::endl;
}


template <int dim>
class ExactSolution {
public:
    // Example exact solution: u(x) = x^2 (for 1D or 2D problems)
    double value(const Point<dim> &p) const {
        if (dim == 1) {
            return p[0] * p[0];  // u(x) = x^2 in 1D
        } else if (dim == 2) {
            return p[0] * p[0] + p[1] * p[1];  // u(x, y) = x^2 + y^2 in 2D
        }
        return 0.0;  // Default case (you can modify for higher dimensions)
    }
};




template <int dim>
double FEM<dim>::l2norm_of_error() {
    ExactSolution<dim> exact_solution;  // Assuming ExactSolution class is defined elsewhere

    double l2norm = 0.0;
    typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(), endc = dof_handler.end();

    for (; elem != endc; ++elem) {
        std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);
        elem->get_dof_indices(local_dof_indices);

        const QGauss<dim> quadrature_formula(fe.degree + 1);
        for (unsigned int q = 0; q < quadrature_formula.size(); ++q) {
            Point<dim> qp = quadrature_formula.point(q);

            double u_exact = exact_solution.value(qp);  // Get exact solution at quadrature point

            double u_numerical = 0.0;
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) {
                u_numerical += D(local_dof_indices[i]) * fe.shape_value(i, qp);
            }

            double diff = u_numerical - u_exact;
            l2norm += diff * diff * elem->measure();
        }
    }

    return std::sqrt(l2norm);  // Return square root of accumulated L2 norm
}







//template <int dim>
//double FEM<dim>::l2norm_of_error() {
//    double l2norm = 0.0;

  //  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(), endc = dof_handler.end();
   // for (; elem != endc; ++elem) {
        // Add L2 norm computation logic here
  //  }

  //  return std::sqrt(l2norm);
//}

template class FEM<1>; // Explicit instantiation for 1D case
template class FEM<2>; // Explicit instantiation for 2D case
