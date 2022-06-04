#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm> 
#include <string>
#include <fstream>
#include <map> 
#include <iterator>
#include <random>
#include <sstream>
#include <filesystem>

namespace fs = std::filesystem;


double tsp(std::vector<std::vector<double>> adj, int n)
{
    double sum = 0;
    int counter = 0;
    int j = 0, i = 0;
    double min = INT_MAX;
    std::map<double, double> mark;

    mark[0] = 1;
    std::vector <double> track(n);
    for (j = 0; j < n; j++)
    {
        if (counter >= n - 1)
        {
            break;
        }
        if (j != 0 && (mark[j] == 0))
        {
            if (adj[0][j] < min)
            {
                min = adj[0][j];
                track[counter] = j + 1;
            }
        }
    }
    for (i = 1; i < n; i++)
    {
        sum += min;
        min = INT_MAX;
        mark[i] = 1;
        counter++;
        for (j = 0; j < n; j++)
        {
            if (counter >= n - 1)
            {
                break;
            }
            if (j != i && (mark[j] == 0))
            {
                if (adj[i][j] < min)
                {
                    min = adj[i][j];
                    track[counter] = j + 1;
                }
            }
        }
    }
    i = track[counter - 1] - 1;
    for (j = 0; j < n; j++)
    {

        if ((i != j) && adj[i][j] < min)
        {
            min = adj[i][j];
            track[counter] = j + 1;
        }
    }
    return sum + min;
}

#include "google/protobuf/duration.pb.h"
#include "ortools/constraint_solver/routing.h"
#include "ortools/constraint_solver/routing_enums.pb.h"
#include "ortools/constraint_solver/routing_index_manager.h"
#include "ortools/constraint_solver/routing_parameters.h"

struct Point 
{
    double x;
    double y;
};

double dist(Point p_1, Point p_2) 
{
    double d = hypot((p_1.x - p_2.x), (p_1.y - p_2.y));
}

void readInput(std::vector<std::vector<double>>& Adj, int N, std::vector<Point>& points) 
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) 
        {
            if (i == j) 
            {
                Adj[i][j] = 0;
            }
            else 
            {
                Adj[i][j] = dist(points[i], points[j]);
            }
        }
    }
}

namespace operations_research 
{
    struct DataModel 
    {
        std::vector<std::vector<double>> distance_Adj;
        std::vector<int64_t> demands;
        std::vector<int64_t> vehicle_capacities;
        int64_t num_vehicles = 0;
        RoutingIndexManager::NodeIndex depot{ 0 };
    };

    void PrintSolution(const DataModel& data, const RoutingIndexManager& manager,
        const RoutingModel& routing, const Assignment& solution) 
    {
        int64_t total_distance{ 0 };
        int64_t total_load{ 0 };
        for (int vehicle_id = 0; vehicle_id < data.num_vehicles; ++vehicle_id) 
        {
            int64_t index = routing.Start(vehicle_id);
            int64_t route_distance{ 0 };
            int64_t route_load{ 0 };
            stringstream route;
            while (routing.IsEnd(index) == false) 
            {
                int64_t node_index = manager.IndexToNode(index).value();
                route_load += data.demands[node_index];
                int64_t previous_index = index;
                index = solution.Value(routing.NextVar(index));
                route_distance += routing.GetArcCostForVehicle(previous_index, index,
                    int64_t{ vehicle_id });
            }
            total_distance += route_distance;
            total_load += route_load;
        }
        std::cout << total_distance << std::endl; 
    }

    void VrpCapacity(DataModel init_data) 
    {
        DataModel data = init_data;
        RoutingIndexManager manager(data.distance_Adj.size(), data.num_vehicles,
            data.depot);
        RoutingModel routing(manager);
        const int transit_callback_index = routing.RegisterTransitCallback(
            [&data, &manager](int64_t from_index, int64_t to_index) -> int64_t 
            {
                int from_node = manager.IndexToNode(from_index).value();
                int to_node = manager.IndexToNode(to_index).value();
                return data.distance_Adj[from_node][to_node];
            });
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);
        const int demand_callback_index = routing.RegisterUnaryTransitCallback(
            [&data, &manager](int64_t from_index) -> int64_t 
            {
                int from_node = manager.IndexToNode(from_index).value();
                return data.demands[from_node];
            });
        routing.AddDimensionWithVehicleCapacity(
            demand_callback_index,   
            int64_t{ 0 },             
            data.vehicle_capacities,
            true,                   
            "Capacity");
        RoutingSearchParameters search_parameters = DefaultRoutingSearchParameters();
        search_parameters.set_first_solution_strategy(
            FirstSolutionStrategy::PATH_CHEAPEST_ARC);
        search_parameters.set_local_search_metaheuristic(
            LocalSearchMetaheuristic::GUIDED_LOCAL_SEARCH);
        search_parameters.mutable_time_limit()->set_seconds(300);
        const Assignment* solution = routing.SolveWithParameters(search_parameters);
        PrintSolution(data, manager, routing, *solution);
    }
}

int main(int argc, char** argv) {
    std::string path = "/Users/annas/vrp_data";
    auto it = fs::directory_iterator(path);
    std::vector<fs::path> array_path;
    copy_if(fs::begin(it), fs::end(it), std::back_inserter(array_path),
        [](const auto& entry) {
            return fs::is_regular_file(entry);
        });
    for (auto& p : array_path) {
        std::ifstream fin;
        fin.open(p.string());
        std::cout << p.string() << std::endl;
        int64_t N, num_vehicles_one, vehicle_capacities_one;
        fin >> N >> num_vehicles_one >> vehicle_capacities_one;
        std::vector<int64_t> buffer(num_vehicles_one, vehicle_capacities_one);
        operations_research::DataModel data;
        data.vehicle_capacities = buffer;
        data.num_vehicles = num_vehicles_one;
        std::vector<Point> points(N);
        std::vector<int64_t> temp(N);
        for (int i = 0; i < N; i++) {
            int64_t demand;
            Point p;
            fin >> demand >> p.x >> p.y;
            points[i] = p;
            temp[i] = demand;
        }
        std::vector<std::vector<double> > Adj(N, std::vector<double>(N));
        readInput(Adj, N, points);
        data.demands = temp;
        data.distance_Adj = Adj;
        operations_research::VrpCapacity(data);
    }
    return EXIT_SUCCESS;
}
