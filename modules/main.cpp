
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <chrono>

#include "effolkronium/random.hpp"

using namespace std::chrono;

using Random = effolkronium::random_static;

namespace py = pybind11;

struct state // Declare a structure to define a state
{        
float x = 0, y = 0, theta = 0;
};

class Node // Declare a structure to define a node
{
    public:

    bool is_init = false;
    state q;
    Node* parent;
    float time_cost = 0.0;
    int rows = 100, cols = 2;
    float trayectory[200];

    Node(){}

    Node(float x, float y, float theta)
    {
        this->q.x = x;
        this->q.y = y;
        this->q.theta = theta;
    }
    ~Node()
    {
        //std::cout << "node detroyed" << std::endl;
    }

};

float wrap_angle(float angle){return (angle + (2.0*M_PI*std::floor((M_PI - angle)/(2.0*M_PI))));}

class rrt_star_dubins
{
    private:
    float step;
    float dominion_x;
    float dominion_y;
    int iterations;
    float resolution;


    int index = 0;
    int parent_id;
    

    Node start; // start node
    
    Node end; // end node

    // pointers to store the addresses of incoming information
    Node* x_rand = new Node;    
    Node* x_nearest = new Node;
    Node* x_new = new Node;     

    Node* V = new Node[10001]; // Planning Tree

    //RRT* utils
    int near[10001];                //list to store the indexes of the near nodes maximum is the same number of nods to store
    Node* x_min = new Node;       
    Node* x_near = new Node;
    Node* temp_x_new = new Node;
    Node* temp_x_near = new Node;

    float c_min;                    // minimal total cost
    float time_to_drive;            // time to drive to a desired goal
    int neighbors;                  // Number of neighbors close to a current node

    //grid mapping checking utils

    std::vector<float> origin  = {0, 0};



    public:
    rrt_star_dubins(float step,
                    float dominion_x,
                    float dominion_y,
                    float resolution,
                    int iterations):
                    step(step),
                    dominion_x(dominion_x),
                    dominion_y(dominion_y),
                    resolution(resolution),
                    iterations(iterations){}



    std::vector<std::vector<float>>  plan(float start_x, float start_y, float start_theta, float end_x, float end_y, py::array_t<int, py::array::c_style | py::array::forcecast> map, float xo, float yo){


        index = 0;       // reset index counter the index counter 

        // Initialization stuff ///////////////////////////////////////
        std::cout << "COMPUTING PATH... " << std::endl;


        // origin initialization
        origin[0] = xo; origin[1] = yo;

        // map parameters
        py::buffer_info map_buff = map.request();
        int *map_ptr = static_cast<int *>(map_buff.ptr);
        int map_Rows = map_buff.shape[0];int map_Cols = map_buff.shape[1];

        // Initializing the start node
        start.is_init = true;
        start.q.x = start_x;
        start.q.y = start_y;
        start.q.theta = start_theta;

        // inserting the start node into the vertex list
        V[0] = start;

        // Initializing the end goal node
        end.q.x = end_x;
        end.q.y = end_y;

        ///////////////////////////////////////////////////////////////

        auto start = high_resolution_clock::now();

        for(int i = 0; i < iterations; i++)
        {
            x_rand->q = Samplefree();

            x_nearest = Nearest(*x_rand); //parent_id = x_nearest - V;

            Steer(*x_nearest, *x_rand, x_new, 1.0); // Dubins Path with time cost 1
            

            if (is_valid(*x_new, map_ptr, map_Rows, map_Cols))
            {

                neighbors = Near(*x_new, 20); //Obtaining the number of neighbor indexes stored in the near array

                *x_min = *x_new; c_min = Cost(*x_nearest) + x_new->time_cost;

                //choose parent

                for(int n = 0; n < neighbors; n++)
                {
                    time_to_drive = std::hypot(V[near[n]].q.x - x_new->q.x,V[near[n]].q.y - x_new->q.y );
                    
                    Steer(V[near[n]], *x_new, temp_x_new,  time_to_drive); // Dubins Path with diffrent time cost

                    if (is_valid(*temp_x_new, map_ptr, map_Rows, map_Cols) && (Cost(V[near[n]]) + time_to_drive) < c_min)
                    {
                        *x_min = *temp_x_new; c_min = (Cost(V[near[n]]) + time_to_drive); 
                    }
                }

                index++; V[index] = *x_min; //Adding a new node

                //rewire

                for(int n = 0; n < neighbors; n++)
                {
                    time_to_drive = std::hypot(V[near[n]].q.x - V[index].q.x,V[near[n]].q.y - V[index].q.y );

                    Steer(V[index], V[near[n]], temp_x_near, time_to_drive+3); // Dubins Path with diffrent time cost

                    if (is_valid(*temp_x_near, map_ptr, map_Rows, map_Cols) && Cost(V[index]) + temp_x_near->time_cost < Cost(V[near[n]]))
                    {
                            V[near[n]] = *temp_x_near;
                    }
                }
            }
        }

        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(stop - start);

        std::cout << "Path computed in: " << duration.count() << " microseconds" << std::endl;

        return return_points();

    }

bool is_valid(Node node_to_check, int* map, int rows, int cols)
{
    float x, y;
    int path_xi, path_yi;

    std::vector<int> cell(2, 0);

    std::vector<int> checker(2, 0);

    for(int j = 0; j < 100; j+=5){

        path_xi = (node_to_check.rows * node_to_check.cols) - 2 - j*node_to_check.cols;
        path_yi = (node_to_check.rows * node_to_check.cols) - 2 - j*node_to_check.cols + 1;
        x = node_to_check.trayectory[path_xi];
        y = node_to_check.trayectory[path_yi];

        cell = __position_to_cell__(x, y);

        if ((cell[0] > 0) && (cell[1] > 0) && (cell[0] < rows) && (cell[1] < cols) && ( (map[(cell[0])*cols + (cell[1])] > 50) || (map[(cell[0])*cols + (cell[1])] < 0) ))
        {
            return false;
        }
        // if (((cell[0] < 0) &&  (cell[0] > rows) && (cell[1] < 0) && (cell[1] > cols)) || (map[(cell[0])*cols + (cell[1])] > 50) || (map[(cell[0])*cols + (cell[1])] < 0))
        // {
        //     return false;
        // }

        // for(int obs = 0; obs < 9; obs++)
        // {
        //     checker = window_checker(obs);
        //     if ((cell[0] > 0) && (cell[1] > 0) && (cell[0] < rows) && (cell[1] < cols) && (map[(cell[0] + checker[0])*cols + (cell[1] + checker[1])] > 50))
        //     {
        //         return false;
        //     }
        // }
    }

    return true;
}

std::vector<int> __position_to_cell__(float x, float y)
{   
    int x_cell, y_cell;

    x_cell = int((x - origin[0]) / resolution);

    y_cell = int((y - origin[1]) / resolution);

    return {x_cell, y_cell};
}

std::vector<int>  window_checker(int index)
{   
    switch(index) {

    case 0 :
        return {0, -1};
    case 1 :
        return {-1, -1};
    case 2 :
        return {-1, 0};
    case 3 :
        return {-1, 1};
    case 4 :
        return {0, 1};
    case 5 :
        return {1, 1};
    case 6:
        return {1, 0};
    case 7:
        return {1, -1};
    case 8:
        return {0, 0};
    default:
        return {0, 0};

    }

}

bool is_in_tree(Node n)
{
    int node_count = index  + 1;

    for (int i = 0; i < node_count; i++)
    {
        if (n.q.x == V[i].q.x && n.q.y == V[i].q.y )
        {
            return true;
        }
    }
    return false;

}

int Near(Node x_new,float r)
{

    int node_count = index  + 1;
    float n = node_count;
    float search_r = std::min(r*std::sqrt(std::log(n)/(n)),step);
    int j = 0;
    float dist;
    for (int i = 0; i < node_count; i++)
    {
        dist = std::hypot(x_new.q.x - V[i].q.x, x_new.q.y - V[i].q.y);
        if (dist < search_r && dist > 0.1)
        {
            near[j] = i;
            j++;
        }
    }
    return j;

}

float Cost(Node current)
{
    float cost = 0;

    while (!current.is_init )
    {   

        cost += current.time_cost;
        current = *current.parent;
    }

    return cost;
}

state Samplefree(float prob = 0.2)
{   
    struct state q;
    if (Random::get(0, 1) < prob)
    {
        q.x = Random::get(-dominion_x, dominion_x);
        q.y = Random::get(-dominion_y, dominion_y);
        q.theta = Random::get(-M_PI, M_PI);
    }
    else
    {
        q = end.q;
    }

    return q;
}

Node* Nearest(Node x_rand)
{
    int node_count = (V + index) - V + 1 ;
    float distances[node_count];
    float* distances_ptr = distances;

    for (int i = 0; i < node_count; i++)
    {
        distances[i] = std::hypot(x_rand.q.x - V[i].q.x, x_rand.q.y - V[i].q.y);
    }

    float* min_element = std::min_element( distances_ptr, distances_ptr + index);

    int x_nearest_id = min_element - distances_ptr;


    return &V[x_nearest_id];
    
}

Node* find_goal(Node goal) //Function that returns the closest node to the goal
{
    int node_count = (V + index) - V + 1 ;

    float distances[node_count];
    float* distances_ptr = distances;

    for (int i = 0; i < node_count; i++)
    {
        distances[i] = std::hypot(goal.q.x - V[i].q.x, goal.q.y - V[i].q.y);
    }

    float* min_element = std::min_element( distances_ptr, distances_ptr + index);

    int x_nearest_id = min_element - distances_ptr;

    //std::cout << "Goal x: " << V[x_nearest_id].q.x << "Goal y: " << V[x_nearest_id].q.y << std::endl;

    return &V[x_nearest_id];
}

void Steer(Node& x_nearest, Node& x_rand, Node* x_new, float T)
{

    // Dubins car steering module

    float angle_lim = 1.2;

    int u = 1;

    float x = x_nearest.q.x, y = x_nearest.q.y;
    float x_d = x_rand.q.x, y_d = x_rand.q.y;
    float dx , dy;
    float theta = x_nearest.q.theta, theta_d, gamma;

    float dt = T/x_new->rows;

    float L = 0.3;


    x_new->time_cost = T;
    x_new->parent = &x_nearest;


    int path_x_i, path_y_i;

    

    for(int i = 0; i < x_new->rows; i++)
    {
        path_x_i = i*x_new->cols;        // x coordinate stored
        path_y_i = (i*x_new->cols) + 1;    // y coordinate stored
        
        x_new->trayectory[path_x_i] = x; // x coordinate stored
        x_new->trayectory[path_y_i] = y; // y coordinate stored

        //std::cout << "x: " << x << " y: " << y <<std::endl;

        //add point to trayectory
        
        dx = x_d - x;
        dy = y_d - y;
        
        theta_d = atan2(dy,dx);
        
        gamma = wrap_angle(theta_d - theta);
        gamma = (gamma < -angle_lim)?-angle_lim:gamma;
        gamma = (gamma > angle_lim)?angle_lim:gamma;

        x += (u*dt)*cos(theta);
        y += (u*dt)*sin(theta);
        
        theta += (u*dt/L) * tan(gamma); theta = wrap_angle(theta);

    }

    x_new->q.x = x;
    x_new->q.y = y;
    x_new->q.theta = theta;

}


// testing and visualization code

std::vector<std::vector<float>> points()
{
    std::vector<std::vector<float>> points(index, std::vector<float>(2));

    for(int i = 0; i < index; i++)
    {
        points[i][0] = V[i].q.x;
        points[i][1] = V[i].q.y;
    }
    
    return points;
} 

std::vector<std::vector<float>> test_dubins_trayectory(int T)
{

    Node start {0.0, 0.0, 0};

    Node end {0.0, 5.0, 0};

    Node* dubins = new Node;

    std::vector<std::vector<float>> points(dubins->rows, std::vector<float>(2));


    Steer(start, end, dubins, T);


    for(int i = 0; i < dubins->rows; i++)
    {
        std::cout << "x: " << dubins->trayectory[i*dubins->cols]<< " y: " << dubins->trayectory[i*dubins->cols + 1] << std::endl; 
        points[i][0] = dubins->trayectory[i*dubins->cols];
        points[i][1] = dubins->trayectory[i*dubins->cols + 1];
    }


    return points;

} 

std::vector<std::vector<float>>  return_points() //std::vector<std::vector<float>>

{
    int points_size = index + 1;
    std::vector<std::vector<float>> points(points_size, std::vector<float>(2));
    for(int i = 0; i < points_size; i++)
    {
            points[i][0] = V[i].q.x;
            points[i][1] = V[i].q.y;
    }

    return points;
}

std::vector<std::vector<float>>  return_path_points() //std::vector<std::vector<float>>
{
    int path_size = 1;
    Node goal = *find_goal(end);
    Node current = goal;

    while (!current.is_init)
    {
        current = *current.parent;
        path_size++;
    }



    int resulting_path_size = path_size;
    std::vector<std::vector<float>> points(resulting_path_size, std::vector<float>(2));

    current = goal;

    int i = 0;
    int index; //index for indexing the dubins path trayectory
    while (!current.is_init)
    {
        points[i][0] = current.q.x;
        points[i][1] = current.q.y;
        i++;
        current = *current.parent;
    }
    return points;
}

std::vector<std::vector<float>> return_path() //std::vector<std::vector<float>>

{
    int path_size = 0;
    Node goal = *find_goal(end);
    Node current = goal;

    while (!current.is_init)
    {
        current = *current.parent;
        path_size++;
    }

    int resulting_path_size = (path_size)*100;
    std::vector<std::vector<float>> points(resulting_path_size, std::vector<float>(2));

    current = goal;

    int path_xi, path_yi;
    int i = 0;

    while (!current.is_init)
    {
        for(int j = 0; j < 100; j++){

            path_xi = (goal.rows * goal.cols) - 2 - j*goal.cols;
            path_yi = (goal.rows * goal.cols) - 2 - j*goal.cols +1;

            points[i][0] = current.trayectory[path_xi];
            points[i][1] = current.trayectory[path_yi];


        i++;
        }

        current = *current.parent;
    }

    return points;
}

std::vector<std::vector<float>> return_tree() //std::vector<std::vector<float>>
{
    int path_size = index+1;

    int tree_size = (path_size)*100;
    std::vector<std::vector<float>> points(tree_size, std::vector<float>(2));

    int path_xi, path_yi;
    int i = 0;
    int w = 0;

    while (i < tree_size)
    {
        for(int j = 0; j < 100; j++){

            path_xi = (200) - 2 - j*2;
            path_yi = (200) - 2 - j*2 +1;

            points[i][0] = V[w].trayectory[path_xi];
            points[i][1] = V[w].trayectory[path_yi];
        i++;
        }
        w ++;
    }
    return points;
}

};




PYBIND11_MODULE(rrt_star_planner, handle) {
    handle.doc() = "There is no wrong time to rock!";
    py::bind_vector<std::vector<std::vector<float>>>(handle, "FloatVector2D");

    py::class_<rrt_star_dubins>(
        handle, "rrt_star_dubins"
    )
    .def(py::init<float, float, float, float, int>())
    .def("plan", &rrt_star_dubins::plan)
    .def("steer", &rrt_star_dubins::Steer)
    .def("test_dubins_trayectory", &rrt_star_dubins::test_dubins_trayectory)
    .def("path", &rrt_star_dubins::return_path)
    .def("path_points", &rrt_star_dubins::return_path_points)
    .def("tree",  &rrt_star_dubins::return_tree)
    ;
}
