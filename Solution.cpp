#include "Solution.hpp"
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <string>
Solution::Solution()
{

}

Solution::Solution(std::vector<double> newSolution, int x, int y, int z)
{
    solution = newSolution;
    size_x = x;
    size_y = y;
    size_z = z;
}

// NOT WORKING
void Solution::resize(int x, int y, int z)
{
    std::vector<double> newSolution (x*y*z, 0);
    size_x = x;
    size_y = y;
    size_z = z;
}

//  x + (y * size_x) + (z * size_x * size_y)
void Solution::assign(int x, int y, int z, double value)
{
    
    if((x + (y * size_x) + (z * size_x * size_y)) <= size_x * size_y * size_z)
    {
        solution[x + (y * size_x) + (z * size_x * size_y)] = value;
    }else{
        std::cout << "Given Point Out Of Field \n" << x << " " << y << " " << z << "\n";
    }
   
}

double Solution::read_at_vertex(int x, int y, int z)
{
    return solution[x + (y * size_x) + (z * size_x * size_y)];
}


std::vector<double> Solution::gradient_at_vertex(int x, int y, int z, std::vector<double> & delta)
{
    // The result.
    std::vector<double> gradient(3, 0);
    double x_delta = delta[0];
    double y_delta = delta[1];
    double z_delta = delta[2];

    //std::cout << "Delta is " << x_delta << " , " << y_delta << " , " << z_delta << " \n";
    if(x == (size_x - 1))
    {
        // gradient(x) = value(x) - value(x-1) / delta
        gradient[0] = (solution[x + (y * size_x) + (z * size_x * size_y)] - solution[x - 1 + (y * size_x) + (z * size_x * size_y)]) / delta[0];
    }else if(x == 0)
    {
        // gradient(x) = value(x + 1) - value(x) / delta
        gradient[0] = (solution[x + 1 + (y * size_x) + (z * size_x * size_y)] - solution[x + (y * size_x) + (z * size_x * size_y)]) / delta[0];
    }else{

        // gradient(x) = (value(x + 1) - value(x - 1)) / 2 * delta
        gradient[0] = (solution[x + 1 + (y * size_x) + (z * size_x * size_y)] - solution[x - 1 + (y * size_x) + (z * size_x * size_y)]) / (2 * delta[0]);
    }

    if(y == (size_y - 1))
    {
        // gradient(y) = value(y) - value(y-1)
        gradient[1] = (solution[x + (y * size_x) + (z * size_x * size_y)] - solution[x + ((y - 1) * size_x) + (z * size_x * size_y)]) / delta[1];
    }else if(y == 0)
    {
        // gradient(y) = value(y + 1) - value(y)
        gradient[1] = (solution[x + ((y + 1) * size_x) + (z * size_x * size_y)] - solution[x + (y * size_x) + (z * size_x * size_y)]) / delta[1];
    }else{

        // gradient(y) = (value(y + 1) - value(y - 1)) / 2
        gradient[1] = (solution[x + ((y + 1) * size_x) + (z * size_x * size_y)] - solution[x + ((y - 1) * size_x) + (z * size_x * size_y)]) / (2 * delta[1]);
    }

    if(z == (size_z - 1))
    {
        // gradient(z) = value(z) - value(z-1)
        gradient[2] = (solution[x + (y * size_x) + (z * size_x * size_y)] - solution[x + (y * size_x) + ((z - 1) * size_x * size_y)]) / delta[2];
    }else if(z == 0)
    {
        // gradient(z) = value(z + 1) - value(z)
        gradient[2] = (solution[x + (y * size_x) + ((z + 1) * size_x * size_y)] - solution[x + (y * size_x) + (z * size_x * size_y)]) / delta[2];
    }else{

        // gradient(z) = (value(z + 1) - value(z - 1)) / 2
        gradient[2] = (solution[x + (y * size_x) + ((z + 1) * size_x * size_y)] - solution[x + (y * size_x) + ((z - 1) * size_x * size_y)]) / (2 * delta[2]);
    }


    return gradient;
}

double Solution::interpolation_solution(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    double solution_000 = read_at_vertex(x_cell, y_cell, z_cell);
    double solution_100 = read_at_vertex(x_cell + 1, y_cell, z_cell);
    double solution_001 = read_at_vertex(x_cell, y_cell, z_cell + 1);
    double solution_101 = read_at_vertex(x_cell + 1, y_cell, z_cell + 1);
    double solution_010 = read_at_vertex(x_cell, y_cell + 1, z_cell);
    double solution_110 = read_at_vertex(x_cell + 1, y_cell + 1, z_cell);
    double solution_011 = read_at_vertex(x_cell, y_cell + 1, z_cell + 1);
    double solution_111 = read_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1);

    double k_x = (x - x_cell);
    double k_y = (y - y_cell);
    double k_z = (z - z_cell);

    double itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
    double itpl_b = solution_000 * (1.0 - k_x) + solution_100 * k_x;
    double itpl_c = solution_011 * (1.0 - k_x) + solution_111 * k_x;
    double itpl_d = solution_010 * (1.0 - k_x) + solution_110 * k_x;

    double itpl_e = itpl_b * (1.0 - k_z) + itpl_a * k_z;
    double itpl_f = itpl_d * (1.0 - k_z) + itpl_c * k_z;

    double result = itpl_e * (1.0 - k_y) + itpl_f * k_y;

    return result;
}

double Solution::interpolation_gradient_x(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    double solution_000 = gradient_at_vertex(x_cell, y_cell, z_cell, delta).at(0);
    double solution_100 = gradient_at_vertex(x_cell + 1, y_cell, z_cell, delta).at(0);
    double solution_001 = gradient_at_vertex(x_cell, y_cell, z_cell + 1, delta).at(0);
    double solution_101 = gradient_at_vertex(x_cell + 1, y_cell, z_cell + 1, delta).at(0);
    double solution_010 = gradient_at_vertex(x_cell, y_cell + 1, z_cell, delta).at(0);
    double solution_110 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell, delta).at(0);
    double solution_011 = gradient_at_vertex(x_cell, y_cell + 1, z_cell + 1, delta).at(0);
    double solution_111 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1, delta).at(0);

    double k_x = (x - x_cell);  
    double k_y = (y - y_cell);
    double k_z = (z - z_cell);

    double itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
    double itpl_b = solution_000 * (1.0 - k_x) + solution_100 * k_x;
    double itpl_c = solution_011 * (1.0 - k_x) + solution_111 * k_x;
    double itpl_d = solution_010 * (1.0 - k_x) + solution_110 * k_x;

    double itpl_e = itpl_b * (1.0 - k_z) + itpl_a * k_z;
    double itpl_f = itpl_d * (1.0 - k_z) + itpl_c * k_z;

    double result = itpl_e * (1.0 - k_y) + itpl_f * k_y;

    return result;
}

double Solution::interpolation_gradient_y(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    double solution_000 = gradient_at_vertex(x_cell, y_cell, z_cell, delta).at(1);
    double solution_100 = gradient_at_vertex(x_cell + 1, y_cell, z_cell, delta).at(1);
    double solution_001 = gradient_at_vertex(x_cell, y_cell, z_cell + 1, delta).at(1);
    double solution_101 = gradient_at_vertex(x_cell + 1, y_cell, z_cell + 1, delta).at(1);
    double solution_010 = gradient_at_vertex(x_cell, y_cell + 1, z_cell, delta).at(1);
    double solution_110 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell, delta).at(1);
    double solution_011 = gradient_at_vertex(x_cell, y_cell + 1, z_cell + 1, delta).at(1);
    double solution_111 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1, delta).at(1);
   

    double k_x = (x - x_cell);  
    double k_y = (y - y_cell);
    double k_z = (z - z_cell);

    double itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
    double itpl_b = solution_000 * (1.0 - k_x) + solution_100 * k_x;
    double itpl_c = solution_011 * (1.0 - k_x) + solution_111 * k_x;
    double itpl_d = solution_010 * (1.0 - k_x) + solution_110 * k_x;

    double itpl_e = itpl_b * (1.0 - k_z) + itpl_a * k_z;
    double itpl_f = itpl_d * (1.0 - k_z) + itpl_c * k_z;

    double result = itpl_e * (1.0 - k_y) + itpl_f * k_y;

    return result;
}

double Solution::interpolation_gradient_z(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    double solution_000 = gradient_at_vertex(x_cell, y_cell, z_cell, delta).at(2);
    double solution_100 = gradient_at_vertex(x_cell + 1, y_cell, z_cell, delta).at(2);
    double solution_001 = gradient_at_vertex(x_cell, y_cell, z_cell + 1, delta).at(2);
    double solution_101 = gradient_at_vertex(x_cell + 1, y_cell, z_cell + 1, delta).at(2);
    double solution_010 = gradient_at_vertex(x_cell, y_cell + 1, z_cell, delta).at(2);
    double solution_110 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell, delta).at(2);
    double solution_011 = gradient_at_vertex(x_cell, y_cell + 1, z_cell + 1, delta).at(2);
    double solution_111 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1, delta).at(2);

    double k_x = (x - x_cell);  
    double k_y = (y - y_cell);
    double k_z = (z - z_cell);

    double itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
    double itpl_b = solution_000 * (1.0 - k_x) + solution_100 * k_x;
    double itpl_c = solution_011 * (1.0 - k_x) + solution_111 * k_x;
    double itpl_d = solution_010 * (1.0 - k_x) + solution_110 * k_x;

    double itpl_e = itpl_b * (1.0 - k_z) + itpl_a * k_z;
    double itpl_f = itpl_d * (1.0 - k_z) + itpl_c * k_z;

    double result = itpl_e * (1.0 - k_y) + itpl_f * k_y;

    return result;
}

double Solution::read_at(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> delta)
{
    return interpolation_solution(x, y, z, x_cell, y_cell, z_cell, delta);
}

std::vector<double> Solution::read_gradient_at(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    double x_g = interpolation_gradient_x(x, y, z, x_cell, y_cell, z_cell, delta);
    double y_g = interpolation_gradient_y(x, y, z, x_cell, y_cell, z_cell, delta);
    double z_g = interpolation_gradient_z(x, y, z, x_cell, y_cell, z_cell, delta);
    std::vector<double> result {x_g, y_g, z_g};
    return result;
}

std::vector<double> Solution::get_solution_vector()
{
    return solution;
}

void Solution::read_from_raw(std::string raw_name)
{
    std::ifstream newRaw(raw_name, std::ios::binary);
    std::vector<double> rawValues(0);
    float rawValue;
    char buf[sizeof(float)];

    while(newRaw.read(buf, sizeof(buf)))
    {
        memcpy(&rawValue, buf, sizeof(rawValue));
        rawValues.push_back(rawValue);
    }
    solution = rawValues;
}