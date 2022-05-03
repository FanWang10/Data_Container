#ifndef SOLUTION_VECTOR_FIELD_HPP
#define SOLUTION_VECTOR_FIELD_HPP

#include <vector>
#include <fstream>
#include <string>
class Solution_Vector_Field
{
    private:
        // The Vector That Have All Vector Field Info
        std::vector<std::vector<double>> solution_vector_field;
        // The Size of Vector Field
        int size_x, size_y, size_z;
        // Interpolate Vector
        std::vector<double> interpolation_solution(double x, double y, double z, int x_cell, int y_cell, int z_cell);
        std::vector<double> interpolation_gradient_x(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);
        std::vector<double> interpolation_gradient_y(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);
        std::vector<double> interpolation_gradient_z(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);

    public:
        Solution_Vector_Field();
        Solution_Vector_Field(std::vector<std::vector<double>> solution, int x, int y, int z);
        void resize(int x, int y, int z);
        void assign(int x, int y, int z, std::vector<double> value_vector);
        std::vector<double> read_at_vertex(int x, int y, int z);
        std::vector<double> read_at(double x, double y, double z, int x_cell, int y_cell, int z_cell);
        std::vector<std::vector<double>> get_solution_vector();
        std::vector<std::vector<double>> get_gradient_at_vertex(int x, int y, int z, std::vector<double> & delta);
        std::vector<std::vector<double>> get_gradient(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);
        std::vector<int> get_size();
};

#endif