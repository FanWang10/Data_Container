#include "Grid.hpp"
#include <iostream>

Grid::Grid()
{
    
}

Grid::Grid(std::vector<double> min, std::vector<double> max, std::vector<int> num_of_pts, std::vector<double> delta)
{
    int size_min_vector = min.size();
    int size_max_vector = max.size();
    int size_delta_vector = delta.size();
    int size_num_of_pt = num_of_pts.size();
    if(size_min_vector == size_max_vector && size_max_vector == size_delta_vector && size_delta_vector == size_num_of_pt)
    {
        //std::cout << "Each Size equals.\n";
        min_coord_each_dimension = min;
        max_coord_each_dimension = max;
        delta_coord_each_dimension = delta;
        num_points_each_dimension = num_of_pts;
        // num_points_each_dimension.resize(size_min_vector, 0);
        // for(int i = 0; i < size_min_vector; i++)
        // {
        //     int num = ((max_coord_each_dimension[i] - min_coord_each_dimension[i]) / delta_coord_each_dimension[i]) + 1;
        //     std::cout << "num_pts_each_dimension " << i << " is equal to " << num << "\n";
        //     //std::cout << "delta " << i << " is equal to " << "\n";
        // }
    }else{
            std::cout << "Size of each vector should be same.";
    }
}

bool Grid::is_point_in_bound(std::vector<double> point)
{
    // std::cout << "The point you give is: ";
    // for(const int& i : point)
    // {
    //     std::cout << i << " ";
    // }
    // std::cout << "\n";
    for(int i = 0; i < point.size(); i++)
    {
        if(!(point[i] <= max_coord_each_dimension[i] && point[i] >= min_coord_each_dimension[i]))
        {
            return false;
        }
    }

    return true;
}

std::vector<int> Grid::get_cell(std::vector<double> point)
{
    std::vector<int> cell(point.size(), 0);
    
    for(int i = 0; i < point.size(); i++)
    {
        cell[i] = point[i];
    }
    
    return cell;
}

std::vector<double> Grid::get_min_coord()
{
    return min_coord_each_dimension;
}

std::vector<double> Grid::get_max_coord()
{
    return max_coord_each_dimension;
}

std::vector<double> Grid::get_delta()
{
    return delta_coord_each_dimension;
}