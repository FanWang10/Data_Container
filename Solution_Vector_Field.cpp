#include "Solution_Vector_Field.hpp"
#include <iostream>

Solution_Vector_Field::Solution_Vector_Field()
{
    
}

Solution_Vector_Field::Solution_Vector_Field(std::vector<std::vector<double>> solution, int x, int y, int z)
{
    solution_vector_field = solution;
    size_x = x;
    size_y = y;
    size_z = z;

}

std::vector<double> scalar_vector_multiplication(std::vector<double> vec, double k)
{
    std::vector<double> result(0, 0);
    for(int i = 0; i < vec.size(); i++)
    {
        result.push_back(vec[i] * k);
    }
    return result;
}

std::vector<double> vector_vector_addition(std::vector<double> vec1, std::vector<double> vec2)
{
    std::vector<double> result(0,0);
    for(int i = 0; i < vec1.size(); i++)
    {
        result.push_back(vec1[i] + vec2[i]);
    }
    return result;
}

// itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
std::vector<double> calculate_interpolation_vector(std::vector<double> & vec1, std::vector<double> & vec2, double k)
{
    std::vector<double> result;

    result = vector_vector_addition(
    scalar_vector_multiplication(vec1, (1.0 - k)),
    scalar_vector_multiplication(vec2, k)
    );

    return result;
}


std::vector<double> Solution_Vector_Field::interpolation_solution(double x, double y, double z, int x_cell, int y_cell, int z_cell)
{

    std::vector<double> solution_000 = read_at_vertex(x_cell, y_cell, z_cell);
    std::vector<double> solution_100 = read_at_vertex(x_cell + 1, y_cell, z_cell);
    std::vector<double> solution_001 = read_at_vertex(x_cell, y_cell, z_cell + 1);
    std::vector<double> solution_101 = read_at_vertex(x_cell + 1, y_cell, z_cell + 1);
    std::vector<double> solution_010 = read_at_vertex(x_cell, y_cell + 1, z_cell);
    std::vector<double> solution_110 = read_at_vertex(x_cell + 1, y_cell + 1, z_cell);
    std::vector<double> solution_011 = read_at_vertex(x_cell, y_cell + 1, z_cell + 1);
    std::vector<double> solution_111 = read_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1);

    double k_x = (x - x_cell);
    double k_y = (y - y_cell);
    double k_z = (z - z_cell);

    // double itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
    std::vector<double> itpl_a = calculate_interpolation_vector(solution_001, solution_101, k_x);
    // double itpl_b = solution_000 * (1.0 - k_x) + solution_100 * k_x;
    std::vector<double> itpl_b = calculate_interpolation_vector(solution_000, solution_100, k_x);
    // double itpl_c = solution_011 * (1.0 - k_x) + solution_111 * k_x;
    std::vector<double> itpl_c = calculate_interpolation_vector(solution_011, solution_111, k_x);
    // double itpl_d = solution_010 * (1.0 - k_x) + solution_110 * k_x;
    std::vector<double> itpl_d = calculate_interpolation_vector(solution_010, solution_110, k_x);

    // double itpl_e = itpl_b * (1.0 - k_z) + itpl_a * k_z;
    std::vector<double> itpl_e = calculate_interpolation_vector(itpl_b, itpl_a, k_z);
    // double itpl_f = itpl_d * (1.0 - k_z) + itpl_c * k_z;
    std::vector<double> itpl_f = calculate_interpolation_vector(itpl_d, itpl_c, k_z);

    // double result = itpl_e * (1.0 - k_y) + itpl_f * k_y;
    std::vector<double> result = calculate_interpolation_vector(itpl_e, itpl_f, k_y);

    return result;
}

void Solution_Vector_Field::resize(int x, int y, int z)
{
    std::vector<double> offset{0, 0, 0};
    std::vector<std::vector<double>> newSolution (x*y*z, offset);
    size_x = x;
    size_y = y;
    size_z = z;
}

void Solution_Vector_Field::assign(int x, int y, int z, std::vector<double> value_vector)
{
    if((x + (y * size_x) + (z * size_x * size_y)) <= size_x * size_y * size_z)
    {
        solution_vector_field[x + (y * size_x) + (z * size_x * size_y)] = value_vector;
    }else{
        std::cout << "Given Point Out Of Field \n" << x << " " << y << " " << z << "\n";
    }
}

std::vector<double> Solution_Vector_Field::read_at_vertex(int x, int y, int z)
{
    return  solution_vector_field[x + (y * size_x) + (z * size_x * size_y)];
}

std::vector<double> Solution_Vector_Field::read_at(double x, double y, double z, int x_cell, int y_cell, int z_cell)
{
    return interpolation_solution(x, y, z, x_cell, y_cell, z_cell);
}

std::vector<std::vector<double>> Solution_Vector_Field::get_solution_vector()
{
    return solution_vector_field;
}

std::vector<double> Solution_Vector_Field::interpolation_gradient_x(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    std::vector<double> gradient_x (3, 0);

    std::vector<double> solution_000 = get_gradient_at_vertex(x_cell, y_cell, z_cell, delta).at(0);
    std::vector<double> solution_100 = get_gradient_at_vertex(x_cell + 1, y_cell, z_cell, delta).at(0);
    std::vector<double> solution_001 = get_gradient_at_vertex(x_cell, y_cell, z_cell + 1, delta).at(0);
    std::vector<double> solution_101 = get_gradient_at_vertex(x_cell + 1, y_cell, z_cell + 1, delta).at(0);
    std::vector<double> solution_010 = get_gradient_at_vertex(x_cell, y_cell + 1, z_cell, delta).at(0);
    std::vector<double> solution_110 = get_gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell, delta).at(0);
    std::vector<double> solution_011 = get_gradient_at_vertex(x_cell, y_cell + 1, z_cell + 1, delta).at(0);
    std::vector<double> solution_111 = get_gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1, delta).at(0);

    double k_x = (x - x_cell);  
    double k_y = (y - y_cell);
    double k_z = (z - z_cell);


    
    // double itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
    std::vector<double> itpl_a (3, 0);
    itpl_a[0] = solution_000[0] * (1.0 - k_x) + solution_100[0] * k_x;
    itpl_a[1] = solution_000[1] * (1.0 - k_x) + solution_100[1] * k_x;
    itpl_a[2] = solution_000[2] * (1.0 - k_x) + solution_100[2] * k_x;
    // double itpl_b = solution_000 * (1.0 - k_x) + solution_100 * k_x;
    std::vector<double> itpl_b (3, 0);
    itpl_b[0] = solution_000[0] * (1.0 - k_x) + solution_100[0] * k_x;
    itpl_b[1] = solution_000[1] * (1.0 - k_x) + solution_100[1] * k_x;
    itpl_b[2] = solution_000[2] * (1.0 - k_x) + solution_100[2] * k_x;
    // double itpl_c = solution_011 * (1.0 - k_x) + solution_111 * k_x;
    std::vector<double> itpl_c (3, 0);
    itpl_c[0] = solution_011[0] * (1.0 - k_x) + solution_111[0] * k_x;
    itpl_c[1] = solution_011[1] * (1.0 - k_x) + solution_111[1] * k_x;
    itpl_c[2] = solution_011[2] * (1.0 - k_x) + solution_111[2] * k_x;
    // double itpl_d = solution_010 * (1.0 - k_x) + solution_110 * k_x;
    std::vector<double> itpl_d (3, 0);
    itpl_d[0] = solution_010[0] * (1.0 - k_x) + solution_110[0] * k_x;
    itpl_d[1] = solution_010[1] * (1.0 - k_x) + solution_110[1] * k_x;
    itpl_d[2] = solution_010[2] * (1.0 - k_x) + solution_110[2] * k_x;
    // double itpl_e = itpl_b * (1.0 - k_z) + itpl_a * k_z;
    std::vector<double> itpl_e (3, 0);
    itpl_e[0] = itpl_b[0] * (1.0 - k_z) + itpl_a[0] * k_z;
    itpl_e[1] = itpl_b[1] * (1.0 - k_z) + itpl_a[1] * k_z;
    itpl_e[2] = itpl_b[2] * (1.0 - k_z) + itpl_a[2] * k_z;
    // double itpl_f = itpl_d * (1.0 - k_z) + itpl_c * k_z;
    std::vector<double> itpl_f (3, 0);
    itpl_f[0] = itpl_d[0] * (1.0 - k_z) + itpl_c[0] * k_z;
    itpl_f[1] = itpl_d[1] * (1.0 - k_z) + itpl_c[1] * k_z;
    itpl_f[2] = itpl_d[2] * (1.0 - k_z) + itpl_c[2] * k_z;
    // double result = itpl_e * (1.0 - k_y) + itpl_f * k_y;
    std::vector<double> result (3, 0);
    result[0] = itpl_e[0] * (1.0 - k_y) + itpl_f[0] * k_y;
    result[1] = itpl_e[1] * (1.0 - k_y) + itpl_f[1] * k_y;
    result[2] = itpl_e[2] * (1.0 - k_y) + itpl_f[2] * k_y;

    return result;
}

std::vector<double> Solution_Vector_Field::interpolation_gradient_y(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    std::vector<double> gradient_y (3, 0);

    // double solution_000 = gradient_at_vertex(x_cell, y_cell, z_cell, delta).at(1);
    std::vector<double> solution_000 = get_gradient_at_vertex(x_cell, y_cell, z_cell, delta).at(1);
    // double solution_100 = gradient_at_vertex(x_cell + 1, y_cell, z_cell, delta).at(1);
    std::vector<double> solution_100 = get_gradient_at_vertex(x_cell + 1, y_cell, z_cell, delta).at(1);
    // double solution_001 = gradient_at_vertex(x_cell, y_cell, z_cell + 1, delta).at(1);
    std::vector<double> solution_001 = get_gradient_at_vertex(x_cell, y_cell, z_cell + 1, delta).at(1);
    // double solution_101 = gradient_at_vertex(x_cell + 1, y_cell, z_cell + 1, delta).at(1);
    std::vector<double> solution_101 = get_gradient_at_vertex(x_cell + 1, y_cell, z_cell + 1, delta).at(1);
    // double solution_010 = gradient_at_vertex(x_cell, y_cell + 1, z_cell, delta).at(1);
    std::vector<double> solution_010 = get_gradient_at_vertex(x_cell, y_cell + 1, z_cell, delta).at(1);
    // double solution_110 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell, delta).at(1);
    std::vector<double> solution_110 = get_gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell, delta).at(1);
    // double solution_011 = gradient_at_vertex(x_cell, y_cell + 1, z_cell + 1, delta).at(1);
    std::vector<double> solution_011 = get_gradient_at_vertex(x_cell, y_cell + 1, z_cell + 1, delta).at(1);
    // double solution_111 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1, delta).at(1);
    std::vector<double> solution_111 = get_gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1, delta).at(1);

    double k_x = (x - x_cell);  
    double k_y = (y - y_cell);
    double k_z = (z - z_cell);

    // double itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
    std::vector<double> itpl_a (3, 0);
    itpl_a[0] = solution_001[0] * (1.0 - k_x) + solution_101[0] * k_x;
    itpl_a[1] = solution_001[1] * (1.0 - k_x) + solution_101[1] * k_x;
    itpl_a[2] = solution_001[2] * (1.0 - k_x) + solution_101[2] * k_x;
    // double itpl_b = solution_000 * (1.0 - k_x) + solution_100 * k_x;
    std::vector<double> itpl_b (3, 0);
    itpl_b[0] = solution_000[0] * (1.0 - k_x) + solution_100[0] * k_x;
    itpl_b[1] = solution_000[1] * (1.0 - k_x) + solution_100[1] * k_x;
    itpl_b[2] = solution_000[2] * (1.0 - k_x) + solution_100[2] * k_x;
    // double itpl_c = solution_011 * (1.0 - k_x) + solution_111 * k_x;
    std::vector<double> itpl_c (3, 0);
    itpl_c[0] = solution_011[0] * (1.0 - k_x) + solution_111[0] * k_x;
    itpl_c[1] = solution_011[1] * (1.0 - k_x) + solution_111[1] * k_x;
    itpl_c[2] = solution_011[2] * (1.0 - k_x) + solution_111[2] * k_x;
    // double itpl_d = solution_010 * (1.0 - k_x) + solution_110 * k_x;
    std::vector<double> itpl_d (3, 0);
    itpl_d[0] = solution_010[0] * (1.0 - k_x) + solution_110[0] * k_x;
    itpl_d[1] = solution_010[1] * (1.0 - k_x) + solution_110[1] * k_x;
    itpl_d[2] = solution_010[2] * (1.0 - k_x) + solution_110[2] * k_x;
    // double itpl_e = itpl_b * (1.0 - k_z) + itpl_a * k_z;
    std::vector<double> itpl_e (3, 0);
    itpl_e[0] = itpl_b[0] * (1.0 - k_z) + itpl_a[0] * k_z;
    itpl_e[1] = itpl_b[1] * (1.0 - k_z) + itpl_a[1] * k_z;
    itpl_e[2] = itpl_b[2] * (1.0 - k_z) + itpl_a[2] * k_z;
    // double itpl_f = itpl_d * (1.0 - k_z) + itpl_c * k_z;
    std::vector<double> itpl_f (3, 0);
    itpl_f[0] = itpl_d[0] * (1.0 - k_z) + itpl_c[0] * k_z;
    itpl_f[1] = itpl_d[1] * (1.0 - k_z) + itpl_c[1] * k_z;
    itpl_f[2] = itpl_d[2] * (1.0 - k_z) + itpl_c[2] * k_z;
    // double result = itpl_e * (1.0 - k_y) + itpl_f * k_y;
    std::vector<double> result (3, 0);
    result[0] = itpl_e[0] * (1.0 - k_y) + itpl_f[0] * k_y;
    result[1] = itpl_e[1] * (1.0 - k_y) + itpl_f[1] * k_y;
    result[2] = itpl_e[2] * (1.0 - k_y) + itpl_f[2] * k_y;

    return result;
}

std::vector<double> Solution_Vector_Field::interpolation_gradient_z(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    std::vector<double> gradient_z (3, 0);

    // double solution_000 = gradient_at_vertex(x_cell, y_cell, z_cell, delta).at(2);
    std::vector<double> solution_000 = get_gradient_at_vertex(x_cell, y_cell, z_cell, delta).at(2);
    // double solution_100 = gradient_at_vertex(x_cell + 1, y_cell, z_cell, delta).at(2);
    std::vector<double> solution_100 = get_gradient_at_vertex(x_cell + 1, y_cell, z_cell, delta).at(2);
    // double solution_001 = gradient_at_vertex(x_cell, y_cell, z_cell + 1, delta).at(2);
    std::vector<double> solution_001 = get_gradient_at_vertex(x_cell, y_cell, z_cell + 1, delta).at(2);
    // double solution_101 = gradient_at_vertex(x_cell + 1, y_cell, z_cell + 1, delta).at(2);
    std::vector<double> solution_101 = get_gradient_at_vertex(x_cell + 1, y_cell, z_cell + 1, delta).at(2);
    // double solution_010 = gradient_at_vertex(x_cell, y_cell + 1, z_cell, delta).at(2);
    std::vector<double> solution_010 = get_gradient_at_vertex(x_cell, y_cell + 1, z_cell, delta).at(2);
    // double solution_110 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell, delta).at(2);
    std::vector<double> solution_110 = get_gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell, delta).at(2);
    // double solution_011 = gradient_at_vertex(x_cell, y_cell + 1, z_cell + 1, delta).at(2);
    std::vector<double> solution_011 = get_gradient_at_vertex(x_cell, y_cell + 1, z_cell + 1, delta).at(2);
    // double solution_111 = gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1, delta).at(2);
    std::vector<double> solution_111 = get_gradient_at_vertex(x_cell + 1, y_cell + 1, z_cell + 1, delta).at(2);

    double k_x = (x - x_cell);  
    double k_y = (y - y_cell);
    double k_z = (z - z_cell);

    // double itpl_a = solution_001 * (1.0 - k_x) + solution_101 * k_x;
    std::vector<double> itpl_a (3, 0);
    itpl_a[0] = solution_001[0] * (1.0 - k_x) + solution_101[0] * k_x;
    itpl_a[1] = solution_001[1] * (1.0 - k_x) + solution_101[1] * k_x;
    itpl_a[2] = solution_001[2] * (1.0 - k_x) + solution_101[2] * k_x;
    // double itpl_b = solution_000 * (1.0 - k_x) + solution_100 * k_x;
    std::vector<double> itpl_b (3, 0);
    itpl_b[0] = solution_000[0] * (1.0 - k_x) + solution_100[0] * k_x;
    itpl_b[1] = solution_000[1] * (1.0 - k_x) + solution_100[1] * k_x;
    itpl_b[2] = solution_000[2] * (1.0 - k_x) + solution_100[2] * k_x;
    // double itpl_c = solution_011 * (1.0 - k_x) + solution_111 * k_x;
    std::vector<double> itpl_c(3, 0);
    itpl_c[0] = solution_011[0] * (1.0 - k_x) + solution_111[0] * k_x;
    itpl_c[1] = solution_011[1] * (1.0 - k_x) + solution_111[1] * k_x;
    itpl_c[2] = solution_011[2] * (1.0 - k_x) + solution_111[2] * k_x;
    // double itpl_d = solution_010 * (1.0 - k_x) + solution_110 * k_x;
    std::vector<double> itpl_d (3, 0);
    itpl_d[0] = solution_010[0] * (1.0 - k_x) + solution_110[0] * k_x;
    itpl_d[1] = solution_010[1] * (1.0 - k_x) + solution_110[1] * k_x;
    itpl_d[2] = solution_010[2] * (1.0 - k_x) + solution_110[2] * k_x;
    // double itpl_e = itpl_b * (1.0 - k_z) + itpl_a * k_z;
    std::vector<double> itpl_e (3, 0);
    itpl_e[0] = itpl_b[0] * (1.0 - k_z) + itpl_a[0] * k_z;
    itpl_e[1] = itpl_b[1] * (1.0 - k_z) + itpl_a[1] * k_z;
    itpl_e[2] = itpl_b[2] * (1.0 - k_z) + itpl_a[2] * k_z;
    // double itpl_f = itpl_d * (1.0 - k_z) + itpl_c * k_z;
    std::vector<double> itpl_f (3, 0);
    itpl_f[0] = itpl_d[0] * (1.0 - k_z) + itpl_c[0] * k_z;
    itpl_f[1] = itpl_d[1] * (1.0 - k_z) + itpl_c[1] * k_z;
    itpl_f[2] = itpl_d[2] * (1.0 - k_z) + itpl_c[2] * k_z;
    // double result = itpl_e * (1.0 - k_y) + itpl_f * k_y;
    std::vector<double> result (3, 0);
    result[0] = itpl_e[0] * (1.0 - k_y) + itpl_f[0] * k_y;
    result[1] = itpl_e[1] * (1.0 - k_y) + itpl_f[1] * k_y;
    result[2] = itpl_e[2] * (1.0 - k_y) + itpl_f[2] * k_y;

    return gradient_z;
}

std::vector<std::vector<double>> Solution_Vector_Field::get_gradient_at_vertex(int x, int y, int z, std::vector<double> & delta)
{   
    std::vector<double> g_x (3, 0);
    std::vector<double> g_y (3, 0);
    std::vector<double> g_z (3, 0);
    double d_x = delta[0];
    double d_y = delta[1];
    double d_z = delta[2];

    if(x == (size_x - 1))
    {
        // gradient(x) = value(x) - value(x-1) / delta
        //gradient[0] = (solution[x + (y * size_x) + (z * size_x * size_y)] - solution[x - 1 + (y * size_x) + (z * size_x * size_y)]) / delta[0];
        g_x[0] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][0] - solution_vector_field[x - 1 + (y * size_x) + (z * size_x * size_y)][0]) / delta[0];
        g_x[1] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][1] - solution_vector_field[x - 1 + (y * size_x) + (z * size_x * size_y)][1]) / delta[0];
        g_x[2] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][2] - solution_vector_field[x - 1 + (y * size_x) + (z * size_x * size_y)][2]) / delta[0];
    }else if(x == 0)
    {
        // gradient(x) = value(x + 1) - value(x) / delta
        //gradient[0] = (solution[x + 1 + (y * size_x) + (z * size_x * size_y)] - solution[x + (y * size_x) + (z * size_x * size_y)]) / delta[0];
        g_x[0] = (solution_vector_field[x + 1 + (y * size_x) + (z * size_x * size_y)][0] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][0]) / delta[0];
        g_x[1] = (solution_vector_field[x + 1 + (y * size_x) + (z * size_x * size_y)][1] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][1]) / delta[0];
        g_x[2] = (solution_vector_field[x + 1 + (y * size_x) + (z * size_x * size_y)][2] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][2]) / delta[0];
    }else{

        // gradient(x) = (value(x + 1) - value(x - 1)) / 2 * delta
        //gradient[0] = (solution[x + 1 + (y * size_x) + (z * size_x * size_y)] - solution[x - 1 + (y * size_x) + (z * size_x * size_y)]) / (2 * delta[0]);
        g_x[0] = (solution_vector_field[x + 1 + (y * size_x) + (z * size_x * size_y)][0] - solution_vector_field[x - 1 + (y * size_x) + (z * size_x * size_y)][0]) / (2 * delta[0]);
        g_x[1] = (solution_vector_field[x + 1 + (y * size_x) + (z * size_x * size_y)][1] - solution_vector_field[x - 1 + (y * size_x) + (z * size_x * size_y)][1]) / (2 * delta[0]);
        g_x[2] = (solution_vector_field[x + 1 + (y * size_x) + (z * size_x * size_y)][2] - solution_vector_field[x - 1 + (y * size_x) + (z * size_x * size_y)][2]) / (2 * delta[0]);
    }

    if(y == (size_y - 1))
    {
        // gradient(y) = value(y) - value(y-1)
        //gradient[1] = (solution[x + (y * size_x) + (z * size_x * size_y)] - solution[x + ((y - 1) * size_x) + (z * size_x * size_y)]) / delta[1];
        g_y[0] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][0] - solution_vector_field[x + ((y - 1) * size_x) + (z * size_x * size_y)][0]) / delta[1];
        g_y[1] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][1] - solution_vector_field[x + ((y - 1) * size_x) + (z * size_x * size_y)][1]) / delta[1];
        g_y[2] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][2] - solution_vector_field[x + ((y - 1) * size_x) + (z * size_x * size_y)][2]) / delta[1];
    }else if(y == 0)
    {
        // gradient(y) = value(y + 1) - value(y)
        //gradient[1] = (solution[x + ((y + 1) * size_x) + (z * size_x * size_y)] - solution[x + (y * size_x) + (z * size_x * size_y)]) / delta[1];
        g_y[0] = (solution_vector_field[x + ((y + 1) * size_x) + (z * size_x * size_y)][0] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][0]) / delta[1];
        g_y[1] = (solution_vector_field[x + ((y + 1) * size_x) + (z * size_x * size_y)][1] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][1]) / delta[1];
        g_y[2] = (solution_vector_field[x + ((y + 1) * size_x) + (z * size_x * size_y)][2] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][2]) / delta[1];
    }else{

        // gradient(y) = (value(y + 1) - value(y - 1)) / 2
        //gradient[1] = (solution[x + ((y + 1) * size_x) + (z * size_x * size_y)] - solution[x + ((y - 1) * size_x) + (z * size_x * size_y)]) / (2 * delta[1]);
        g_y[0] = (solution_vector_field[x + ((y + 1) * size_x) + (z * size_x * size_y)][0] - solution_vector_field[x + ((y - 1) * size_x) + (z * size_x * size_y)][0]) / (2 * delta[1]);
        g_y[1] = (solution_vector_field[x + ((y + 1) * size_x) + (z * size_x * size_y)][1] - solution_vector_field[x + ((y - 1) * size_x) + (z * size_x * size_y)][1]) / (2 * delta[1]);
        g_y[2] = (solution_vector_field[x + ((y + 1) * size_x) + (z * size_x * size_y)][2] - solution_vector_field[x + ((y - 1) * size_x) + (z * size_x * size_y)][2]) / (2 * delta[1]);
    }

    if(z == (size_z - 1))
    {
        // gradient(z) = value(z) - value(z-1)
        //gradient[2] = (solution[x + (y * size_x) + (z * size_x * size_y)] - solution[x + (y * size_x) + ((z - 1) * size_x * size_y)]) / delta[2];
        g_z[0] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][0] - solution_vector_field[x + (y * size_x) + ((z - 1) * size_x * size_y)][0]) / delta[2];
        g_z[1] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][1] - solution_vector_field[x + (y * size_x) + ((z - 1) * size_x * size_y)][1]) / delta[2];
        g_z[2] = (solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][2] - solution_vector_field[x + (y * size_x) + ((z - 1) * size_x * size_y)][2]) / delta[2];
    }else if(z == 0)
    {
        // gradient(z) = value(z + 1) - value(z)
        //gradient[2] = (solution[x + (y * size_x) + ((z + 1) * size_x * size_y)] - solution[x + (y * size_x) + (z * size_x * size_y)]) / delta[2];
        g_z[0] = (solution_vector_field[x + (y * size_x) + ((z + 1) * size_x * size_y)][0] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][0]) / delta[2];
        g_z[1] = (solution_vector_field[x + (y * size_x) + ((z + 1) * size_x * size_y)][1] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][1]) / delta[2];
        g_z[2] = (solution_vector_field[x + (y * size_x) + ((z + 1) * size_x * size_y)][2] - solution_vector_field[x + (y * size_x) + (z * size_x * size_y)][2]) / delta[2];
    }else{

        // gradient(z) = (value(z + 1) - value(z - 1)) / 2
        //gradient[2] = (solution[x + (y * size_x) + ((z + 1) * size_x * size_y)] - solution[x + (y * size_x) + ((z - 1) * size_x * size_y)]) / (2 * delta[2]);
        g_z[0] = (solution_vector_field[x + (y * size_x) + ((z + 1) * size_x * size_y)][0] - solution_vector_field[x + (y * size_x) + ((z - 1) * size_x * size_y)][0]) / (2 * delta[2]);
        g_z[1] = (solution_vector_field[x + (y * size_x) + ((z + 1) * size_x * size_y)][1] - solution_vector_field[x + (y * size_x) + ((z - 1) * size_x * size_y)][1]) / (2 * delta[2]);
        g_z[2] = (solution_vector_field[x + (y * size_x) + ((z + 1) * size_x * size_y)][2] - solution_vector_field[x + (y * size_x) + ((z - 1) * size_x * size_y)][2]) / (2 * delta[2]);
    }

    std::vector<std::vector<double>> gradient {g_x, g_y, g_z};

    return gradient;
}

std::vector<std::vector<double>> Solution_Vector_Field::get_gradient(double x, double y, double z, int x_cell, int y_cell, int z_cell, std::vector<double> & delta)
{
    std::vector<double> g_x = interpolation_gradient_x(x, y, z, x_cell, y_cell, z_cell, delta);
    std::vector<double> g_y = interpolation_gradient_y(x, y, z, x_cell, y_cell, z_cell, delta);
    std::vector<double> g_z = interpolation_gradient_z(x, y, z, x_cell, y_cell, z_cell, delta);
    std::vector<std::vector<double>> gradient {g_x, g_y, g_z};

    return gradient;
}

std::vector<int> Solution_Vector_Field::get_size()
{
    std::vector<int> size {size_x, size_y, size_z};
    return size;
}