#include <iostream>
#include "Grid.hpp"
#include "Solution.hpp"
#include "Solution_Vector_Field.hpp"
#include "Field.hpp"
#include "png.h"
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include "matplotlibcpp.h"
#include <stdlib.h>
#include <netcdf>



#define STEP_SIZE 1
using namespace std;
namespace plt = matplotlibcpp;


void print_vector(std::vector<double> vector)
{
    cout << "Vector Print: ";
    for(int i = 0; i < vector.size(); i++)
    {
        cout << vector[i] << "  ";
    }
    cout << "\n";
}

std::vector<int> read_dimension(ifstream & istrm)
{
    int dim_x, dim_y, dim_z;
    std::vector<int> result(3, 0);
    istrm.read(reinterpret_cast<char*>(&dim_x), sizeof dim_x);
    istrm.read(reinterpret_cast<char*>(&dim_y), sizeof dim_y);
    istrm.read(reinterpret_cast<char*>(&dim_z), sizeof dim_z);
    result[0] = dim_x;
    result[1] = dim_y;
    result[2] = dim_z;
    return result;
}

void read_vector(ifstream & istrm, std::vector<double> & vec)
{
    float y, x, z;
    istrm.read(reinterpret_cast<char*>(&y), sizeof y);
    istrm.read(reinterpret_cast<char*>(&x), sizeof x);
    istrm.read(reinterpret_cast<char*>(&z), sizeof z);
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;

}

std::vector<std::vector<double>> euler_method(double x, double y, double z, int steps, 
	Field & field)
{
	int i = 0;
	std::vector<double> init_position{x, y, z};
	std::vector<std::vector<double>> trace;
	trace.push_back(init_position);
	std::vector<double> v = field.get_vector_value(x, y, z);
	for(i; i < steps; i++)
	{
		// next point p = init_p + v * h;
		x = x + v[0] * STEP_SIZE;
		y = y + v[1] * STEP_SIZE;
		z = z + v[2] * STEP_SIZE;

		// create new point
		std::vector<double> new_pt {x, y, z};
		
		// push it to trace
		trace.push_back(new_pt); 

		// Update the velocity
		v = field.get_vector_value(x, y, z);
	}
	return trace;	
}

std::vector<std::vector<double>> runge_kutta_four(double x, double y, double z, int steps,
        Field & field)
{
	int i = 0;
        std::vector<double> init_position{x, y, z};
	std::vector<std::vector<double>> trace;
	trace.push_back(init_position);
	std::vector<double> v1 = field.get_vector_value(x, y, z);
	std::vector<double> v2 = field.get_vector_value(x + (STEP_SIZE * v1[0]) / 2.0, 
							y + (STEP_SIZE * v1[1]) / 2.0, 
							z + (STEP_SIZE * v1[2]) / 2.0);
	std::vector<double> v3 = field.get_vector_value(x + (STEP_SIZE * v2[0]) / 2.0,
							y + (STEP_SIZE * v2[1]) / 2.0,
							z + (STEP_SIZE * v2[2]) / 2.0);
	std::vector<double> v4 = field.get_vector_value(x + STEP_SIZE * v3[0],
							y + STEP_SIZE * v3[1],
							z + STEP_SIZE * v3[2]);
	for(i; i < steps; i ++)
	{
		x = x + (STEP_SIZE * (v1[0] + 2 * v2[0] + 2 * v3[0] + v4[0])) / 6.0;
		y = y + (STEP_SIZE * (v1[1] + 2 * v2[1] + 2 * v3[1] + v4[1])) / 6.0;
		z = z + (STEP_SIZE * (v1[2] + 2 * v2[2] + 2 * v3[2] + v4[2])) / 6.0; 
		
		//Create New Point
		std::vector<double> new_pt {x, y, z};
		
		// Push it to trace
		trace.push_back(new_pt);

		// Update all the velocity
        	v1 = field.get_vector_value(x, y, z);
	        v2 = field.get_vector_value(x + (STEP_SIZE * v1[0]) / 2.0,
                                                        y + (STEP_SIZE * v1[1]) / 2.0,
                                                        z + (STEP_SIZE * v1[2]) / 2.0);
        	v3 = field.get_vector_value(x + (STEP_SIZE * v2[0]) / 2.0,
                                                        y + (STEP_SIZE * v2[1]) / 2.0,
                                                        z + (STEP_SIZE * v2[2]) / 2.0);
        	v4 = field.get_vector_value(x + STEP_SIZE * v3[0],
                                                        y + STEP_SIZE * v3[1],
                                                        z + STEP_SIZE * v3[2]);
	}

	return trace;
}

void segment_data(std::vector<std::vector<double>> & data, std::vector<double> & x, 
			std::vector<double> & y, std::vector<double> & z)
{
	std::vector<double> curr_pt;
	for(int i = 0; i < data.size(); i++)
	{
		curr_pt = data[i];
		x.push_back(curr_pt[0]);
		y.push_back(curr_pt[1]);
		z.push_back(curr_pt[2]); 
	}
}

void draw(std::vector<std::vector<double>> & data)
{
	std::vector<double> x, y, z;
	segment_data(data, x, y, z);
	plt::plot3(x, y, z);
}

int read_max_steps(ifstream & file)
{
	int max_steps = 0;
	string line; 
	std::getline(file, line);
	max_steps = atoi(&line[9]);
	return max_steps;
}

void read_seed(ifstream & file, std::vector<double> & seed)
{
	file >> seed[0];
	file.ignore();
	file >> seed[1];
	file.ignore();
	file >> seed[2];
}

void add_two_vector(vector<double> & v1, vector<double> & v2)
{
    v1.at(0) = v1.at(0) + v2.at(0);
    v1.at(1) = v1.at(1) + v2.at(1);
    v1.at(2) = v1.at(2) + v2.at(2);
}

void set_color(vector<double> & color, double red, double green, double blue)
{
    color[0] = red;
    color[1] = green;
    color[2] = blue;
}

vector<double> color_lookup(double scalar_value)
{
    vector<double> color {0, 0, 0};
    double k, r, g, b;

    // Color LookUp For Heart Data
    if(scalar_value >= 0)
    {
        set_color(color, (100 * scalar_value)  * 0.9, (100 * scalar_value) * 0.5, (100 * scalar_value) / 10 * 0.4);
    }


    return color;
} 

double alpha_lookup(double scalar_value)
{
    double alpha = 1;
    double k = 0;

    if(scalar_value >= 0)
    {
        alpha = 100 * scalar_value;
    }


    return alpha;
}

vector<double> blend_cout(vector<double> & cin, vector<double> & color, double alpha, double alpha_in)
{
    vector<double> cout {0, 0, 0};
    cout[0] = cin[0] + color[0] * alpha * (1 - alpha_in);
    cout[1] = cin[1] + color[1] * alpha * (1 - alpha_in);
    cout[2] = cin[2] + color[2] * alpha * (1 - alpha_in);
    return cout;
}

double blend_aout(double alpha_in, double alpha)
{
    return alpha_in + alpha * (1 - alpha_in);
}


void generate_image(Field field, int num_of_pixel_x, int num_of_pixel_y, int max_depth, string fileName, std::vector<double> min, std::vector<double> max)
{
	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;

    char * name = new char[fileName.length() + 1];
    strcpy(name, fileName.c_str());

	fp = fopen(name, "wb");

	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, num_of_pixel_x, num_of_pixel_y,
			8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);
    row = (png_bytep) malloc(3 * num_of_pixel_x * sizeof(png_byte));

    vector<double> step_ray_direction {0, 0, -max[3] / max_depth};

	for(int j = 0; j < num_of_pixel_y; j++)
    {
        for(int i = 0; i < num_of_pixel_x; i++)
        {
            double length_of_each_pixel = (max[0] - min[0]) / num_of_pixel_x;
            double start_position0 = min[0];
            double start_position1 = min[1];
            double start_position2 = min[2];

            std::vector<double> location_in_volume {(i + 0.5) * (max[0] / num_of_pixel_x), (j + 0.5) * (max[1]/ num_of_pixel_y), max[2] - 0.5};

            vector<double> color {0, 0, 0};
            double alpha = 0;

            // Initial Color Black
            vector<double> cout{0, 0, 0};
            vector<double> cin = cout;

            // Initial Alpha 0
            double alpha_out = 0;
            double alpha_in = 0;

            for(int k = 0; k < max_depth; k++)
            {
                // p = p + step_ray_direction
                add_two_vector(location_in_volume, step_ray_direction);
                
                // Get Value 
                std::vector<std::vector<double>> gradient = field.get_gradient_vector_field(location_in_volume[0], location_in_volume[1], location_in_volume[2]);
                std::vector<double> vorticity = field.get_vorticity(location_in_volume[0], location_in_volume[1], location_in_volume[2], gradient);
                double scalar_value = field.get_vorticity_magnitude(vorticity);

                
                // Color Look Ups
                color = color_lookup(scalar_value);

                // Alpha Look Ups
                alpha = alpha_lookup(scalar_value);

                cout = blend_cout (cin, color, alpha, alpha_in);
                alpha_out = blend_aout(alpha_in, alpha);
                cin = cout;
                alpha_in = alpha_out;
            }
            row[0 + i * 3] = 255.0 * cout[0];
            row[1 + i * 3] = 255.0 * cout[1];
            row[2 + i * 3] = 255.0 * cout[2];

        }
        png_write_row(png_ptr, row);
    }

    png_write_end(png_ptr, NULL);
    if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);
}



double calculate_magnitude(std::vector<double> vec)
{
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}



int main(int argc, char *argv[])
{
    // Read File
    string filename = "tornadoPC_96.vec";
    ifstream istrm(filename, std::ios::binary);

    int xdim, ydim, zdim;
    std::vector<double> vec(3, 0);
    if(!istrm.is_open())
    {
        cout << "Failed to open: " << filename << '\n';
    } else {
        /*
        **  Initialization
        **  
        */
        // read dimension
        std::vector<int> dim = read_dimension(istrm);
        xdim = dim[0];
        ydim = dim[1];
        zdim = dim[2];

        // read all vec in
        std::vector<std::vector<double>> vector_fields;
        for(int i = 0; i < xdim * ydim * zdim; i++)
        {   
            std::vector<double> vec (3,0);
            read_vector(istrm, vec);
            vector_fields.push_back(vec);
        }

        // Create Solution Vector
        Solution_Vector_Field solution(vector_fields, xdim, ydim, zdim);

        // Create Grid
        std::vector<double> min{0, 0, 0};
        std::vector<double> max{(double)(xdim -1), (double)(ydim -1), (double)(zdim -1)};
        std::vector<int> num_of_pt{xdim, ydim, zdim};
        std::vector<double> delta{1, 1, 1};
        Grid grid(min, max, num_of_pt, delta);

        // Create Field
        Field field(solution, grid);


        /*
        ** Initialization Finished
        */

	/*
	** Read Seeds
	*/
	std::cout << "Enter Seeds File Name: " ;
	string seeds_file_name;
	std::cin >> seeds_file_name;
	ifstream seeds_file(seeds_file_name);
	if(!seeds_file.is_open())
	{
		std::cout << "File Not Open! \n";
	}
	int maxSteps = read_max_steps(seeds_file);
	std::vector<double> seed {0, 0, 0};

	bool perform_euler = false;
	std::cout << "Perform Euler Method? (Y/y for yes, N/n for no) :";
	char input;
	std::cin >> input;

	if(input == 'y' || input == 'Y')
	{
		perform_euler = true;
	}
    	/*
	**  Task 1: Euler Method
    	*/
	if(perform_euler){
	std::vector<double> x_euler, y_euler, z_euler;
	for(int j = 0; j < 100; j ++){
		read_seed(seeds_file, seed);
		std::vector<std::vector<double>> t = euler_method(seed[0], seed[1], seed[2], maxSteps, field);
		//std::vector<double> x_euler, y_euler, z_euler;
		for(int i = 0; i < t.size(); i++)
		{
		 	std::vector<double> p = t[i];
		 	x_euler.push_back(p[0]);
		 	y_euler.push_back(p[1]);
		 	z_euler.push_back(p[2]);
		}
	}
	plt::plot3(x_euler, y_euler, z_euler);
	plt::title("Euler Method");
	plt::legend();
	plt::show();
	}
	
	PyObject * test;

	/*
	** Task 2: Runge-Kutta 
	*/
	if(!perform_euler){
    bool perform_rk4 = false;
	std::cout << "Perform Runge-Kutta4? (Y/y for yes, N/n for no) :";
	char input_rk4;
	std::cin >> input_rk4;

	if(input_rk4 == 'y' || input_rk4 == 'Y')
	{
		perform_rk4 = true;
	}
    if(perform_rk4){
	std::vector<double> x_runge, y_runge, z_runge;
	for(int j = 0; j < 100; j ++){
		read_seed(seeds_file, seed);
		std::vector<std::vector<double>> t = runge_kutta_four(seed[0], seed[1], seed[2], maxSteps, field);
		for(int i = 0; i < t.size(); i++)
		{
		 	std::vector<double> p = t[i];
		 	x_runge.push_back(p[0]);
		 	y_runge.push_back(p[1]);
		 	z_runge.push_back(p[2]);
		}
	}
	plt::plot3(x_runge, y_runge, z_runge);
	plt::title("Runge-Kutta-4");
	plt::legend();
	plt::show();
    }
	}
	/*
	** Task 3: Calculate vorticiy and generate image
	*/
    bool perform_vorticity = false;
	std::cout << "Generate Vorticity Magnitude Volume Render Image? (Y/y for yes, N/n for no) :";
	char input_vorticity;
	std::cin >> input_vorticity;

	if(input_vorticity == 'y' || input_vorticity == 'Y')
	{
		perform_vorticity = true;
	}
    if(perform_vorticity)
    {
        string filename_image;
        std::cout << "What's the name of your image?: ";
        std::cin >> filename_image;
        int pixel_x, pixel_y, depth;
        std::cout << "What's the number of pixels in x direction? ";
        std::cin >> pixel_x;
        std::cout << "What's the number of pixels in y direction? ";
        std::cin >> pixel_y;
        std::cout << "What's the depth of volume rendering? ";
        std::cin >> depth;
        generate_image(field, pixel_x, pixel_y, depth, filename_image, min, max);

    }
    

    /*
    ** Task 4: Dump Into paraview
    */
    bool perform_paraview = false;
	std::cout << "Generate Paraview Readable File? (Y/y for yes, N/n for no) :";
	char input_paraview;
	std::cin >> input_paraview;

	if(input_paraview == 'y' || input_paraview == 'Y')
	{
		perform_paraview = true;
	}
    if(perform_paraview){
    std::cout << "Enter paraview filename: (have to end with .raw): ";
    string name;
    std::cin >> name;
    char * namecstr = new char[name.length() + 1];
    strcpy(namecstr, name.c_str());
    field.write_to_raw(namecstr);
    }
	}



}


    //
    //  How to read from binary
    //
    // if(!istrm.is_open())
    // {
    //     cout << "failed to open":  << filename << "\n";
    // } else {
    //     int i1, i2, i3;
    //     istrm.read(reinterpret_cast<char*>(&i1), sizeof i1);
    //     istrm.read(reinterpret_cast<char*>(&i2), sizeof i2);
    //     istrm.read(reinterpret_cast<char*>(&i3), sizeof i3);
    //     std::cout << "read from file: " << i1 << " " << i2 << " " << i3 << '\n';

    //     float y, x, z;
    //     istrm.read(reinterpret_cast<char*>(&y), sizeof y);
    //     istrm.read(reinterpret_cast<char*>(&x), sizeof x);
    //     istrm.read(reinterpret_cast<char*>(&z), sizeof z);
    //     std::cout << "read from file: " << y << " " << x << " " << z << '\n';
    // }

    // std::vector<double> fist_vector = field.get_vector_value(1.5, 0, 0);
    // cout << fist_vector[0] << " , " << fist_vector[1] << " , " << fist_vector[2] << '\n';
