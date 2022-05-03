#ifndef GRID_HPP
#define GRID_HPP

#include <vector>

class Grid
{
private:
    std::vector<int> num_points_each_dimension;
    std::vector<double> min_coord_each_dimension;
    std::vector<double> max_coord_each_dimension;
    std::vector<double> delta_coord_each_dimension;

public:
    Grid();
    Grid(std::vector<double> min, std::vector<double> max, std::vector<int> num_of_pts, std::vector<double> delta);
    bool is_point_in_bound(std::vector<double> point); 
    std::vector<int> get_cell(std::vector<double> point);
    std::vector<double> get_min_coord();
    std::vector<double> get_max_coord();
    std::vector<double> get_delta();
};


#endif