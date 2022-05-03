#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <vector>
#include <fstream>
#include <string>
class Solution
{
    private:
        std::vector<double> solution;
        int size_x, size_y, size_z;
        double interpolation_solution(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);
        double interpolation_gradient_x(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);
        double interpolation_gradient_y(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);
        double interpolation_gradient_z(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);
    public:
        Solution();
        Solution(std::vector<double> solution, int x, int y, int z);
        void resize(int x, int y, int z);
        void assign(int x, int y, int z, double value);
        double read_at_vertex(int x, int y, int z);
        std::vector<double> gradient_at_vertex(int x, int y, int z, std::vector<double> & delta);
        double read_at(double x, double y, double z, int x_cell, int y_cell, int z_cell,  std::vector<double> delta);
        std::vector<double> read_gradient_at(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta);
        std::vector<double> get_solution_vector();
        void read_from_raw(std::string raw_name);
};


#endif