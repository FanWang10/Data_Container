#ifndef FIELD_HPP
#define FIELD_HPP

#include "Solution.hpp"
#include "Solution_Vector_Field.hpp"
#include "Grid.hpp"
#include <vector>

class Field
{
    private:
        Solution solution;
        Solution_Vector_Field solution_vector_field;
        Grid grid;
        std::vector<double> get_vorticity(double x, double y, double z);
    public:
        Field(Solution s, Grid g);
        // new -- implemented
        
        Field(Solution_Vector_Field s, Grid g);
        void assign_grid(Grid g);
        void assign_solution(Solution s);
        // new -- implemented
        void assign_solution_vector_field(Solution_Vector_Field s);
        double get_value(double x, double y, double z);
        // new -- implemented
        std::vector<double> get_vector_value(double x, double y, double z);
        std::vector<double> get_gradient(double x, double y, double z);
        std::vector<std::vector<double>> get_gradient_vector_field(double x, double y, double z);
        bool is_point_in_bound(double x, double y, double z);
        std::vector<int> get_length_per_dim_grid();
        std::vector<int> get_length_per_dim_solution();
        double get_vorticity_magnitude(std::vector<double> vorticity);
        std::vector<double> project_on_solution_coord(std::vector<double> v);
        std::vector<double> get_vorticity(double x, double y, double z, std::vector<std::vector<double>> gradient);
        void write_to_raw(char * filename);
};


#endif