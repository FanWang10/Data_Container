#include <vector>
#include "Field.hpp"
#include "Solution.hpp"
#include "Solution_Vector_Field.hpp"
#include "Grid.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <cstring>

Field::Field(Solution s, Grid g)
{
    solution = s;
    grid = g;
}

Field::Field(Solution_Vector_Field s, Grid g)
{
    solution_vector_field = s;
    grid = g;
}

void Field::assign_grid(Grid g)
{
    grid = g;
}

void Field::assign_solution(Solution s)
{
    solution = s;
}

void Field::assign_solution_vector_field(Solution_Vector_Field s)
{
    solution_vector_field = s;
}

std::vector<double> Field::project_on_solution_coord(std::vector<double> v)
{
    std::vector<double> result = v;
    std::vector<double> min = grid.get_min_coord();
    std::vector<double> max = grid.get_max_coord();
    std::vector<double> delta = grid.get_delta();
    for(int i = 0; i < v.size(); i++)
    {
        v.at(i) = (v.at(i) - min.at(i)) / delta.at(i);
    }
    return v;
}

double Field::get_value(double x, double y, double z)
{
    /* Do Some Thing Here So That x, y, z: [-1, 1] -> [0, 255] */
    std::vector<double> pt_old {x, y, z};
    //std::cout << "pt_old " << x << " " << y << " " << z << " \n";
    std::vector<double> pt = project_on_solution_coord(pt_old);
    //std::cout << "pt " << pt.at(0) << " " << pt.at(1) << " " << pt.at(2) << "\n";
    std::vector<int> cell_index = grid.get_cell(pt);
    //std::cout << "cell_index " << cell_index.at(0) << " " << cell_index.at(1) << " " << cell_index.at(2) << " \n";
    std::vector<double> delta = grid.get_delta();
    return solution.read_at(pt.at(0), pt.at(1), pt.at(2), cell_index[0], cell_index[1], cell_index[2], delta);
}

std::vector<double> Field::get_vector_value(double x, double y, double z)
{
    std::vector<double> pt_old {x, y, z};
    std::vector<double> pt = project_on_solution_coord(pt_old);
    //std::cout << "pt " << pt.at(0) << " " << pt.at(1) << " " << pt.at(2) << "\n";
    std::vector<int> cell_index = grid.get_cell(pt);
    //std::cout << "cell_index " << cell_index.at(0) << " " << cell_index.at(1) << " " << cell_index.at(2) << " \n";
    return solution_vector_field.read_at(pt.at(0), pt.at(1), pt.at(2), cell_index[0], cell_index[1], cell_index[2]);

}

std::vector<double> Field::get_gradient(double x, double y, double z)
{
    std::vector<double> pt {x, y, z};
    std::vector<double> pt_new = project_on_solution_coord(pt);
    std::vector<double> delta = grid.get_delta();
    std::vector<int> cell_index = grid.get_cell(pt_new);
    return solution.read_gradient_at(pt_new.at(0), pt_new.at(1), pt_new.at(2), cell_index[0], cell_index[1], cell_index[2], delta);
}

std::vector<std::vector<double>> Field::get_gradient_vector_field(double x, double y, double z)
{
    std::vector<double> pt {x, y, z};
    std::vector<double> pt_new = project_on_solution_coord(pt);
    std::vector<double> delta = grid.get_delta();
    std::vector<int> cell_index = grid.get_cell(pt_new);
    return solution_vector_field.get_gradient(pt_new.at(0), pt_new.at(1), pt_new.at(2), cell_index[0], cell_index[1], cell_index[2], delta);
}


bool Field::is_point_in_bound(double x, double y, double z)
{
    return false;
}

std::vector<double> Field::get_vorticity(double x, double y, double z, std::vector<std::vector<double>> gradient)
{
    std::vector<double> vorticity (3, 0);
    
    vorticity[0] = gradient[1][2] - gradient[2][1];
    vorticity[1] = gradient[2][0] - gradient[0][2];
    vorticity[2] = gradient[0][1] - gradient[1][0];

    return vorticity;
}

double Field::get_vorticity_magnitude(std::vector<double> vorticity)
{
    double vorticity_magnitude = sqrt(vorticity[0] * vorticity[0] + vorticity[1] * vorticity[1] + vorticity[2] * vorticity[2]);

    std::cout << "Vorticity_magnitude is: " << vorticity_magnitude << "\n";

    return vorticity_magnitude;
}

void Field::write_to_raw(char * filename)
{
    //std::ofstream file(filename);
    //std::vector<std::vector<double>> solution_tmp = solution_vector_field;
    FILE* out = fopen(filename, "w");
    std::vector<int> size = solution_vector_field.get_size();
    for(int i = 0; i < size[2]; i++)
    {
        for(int j = 0; j < size[1]; j++)
        {
            for(int k = 0; k < size[0]; k++)
            {
                std::vector<double> pt = solution_vector_field.read_at_vertex(k, j, i);
                float pt_float[] = {(float)pt[0], (float)pt[1], (float)pt[2]};
                fwrite(pt_float, sizeof(float), 3, out);
            }
        }
    }
    fclose(out);
}

