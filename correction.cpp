#include <bits/stdc++.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <python3.10/Python.h>

using namespace std;

double toRadians(double degrees) {
    return degrees * M_PI / 180.0;
}

// Function to calculate the distance between two coordinates using Haversine formula
double haversineDistance(pair<double, double> &p1, double lat2, double lon2){

    double lat1 = p1.first, lon1 = p1.second;

    // Earth's radius in kilometers
    const double earthRadiusKm = 6371.0;

    // Convert latitude and longitude from degrees to radians
    lat1 = toRadians(lat1);
    lon1 = toRadians(lon1);
    lat2 = toRadians(lat2);
    lon2 = toRadians(lon2);

    // Haversine formula
    double dLat = lat2 - lat1;
    double dLon = lon2 - lon1;
    double a = sin(dLat / 2) * sin(dLat / 2) + cos(lat1) * cos(lat2) * sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double distance = earthRadiusKm * c;

    return distance;
}

int timeStringToSeconds(const std::string& timeStr) {
    int hours, minutes, seconds;
    std::istringstream iss(timeStr);
    char delimiter;
    iss >> hours >> delimiter >> minutes >> delimiter >> seconds;
    return hours * 3600 + minutes * 60 + seconds;
}

int timeDifferenceInSeconds(const std::string& time1, const std::string& time2) {
    int seconds1 = timeStringToSeconds(time1);
    int seconds2 = timeStringToSeconds(time2);
    return seconds2 - seconds1;
}

// Custom comparison function for pairs of strings
bool comparePairs(const pair<string, string>& a, const pair<string, string>& b) {
    return a.first <= b.first;
}

std::string executePythonCode(const std::string& pythonCode) {
    Py_Initialize();

    // Import the 'sys' module
    PyObject* sys = PyImport_ImportModule("sys");
    if (!sys) {
        PyErr_Print();
        std::cerr << "Failed to import the 'sys' module\n";
        Py_Finalize();
        return "";
    }

    // Redirect Python's stdout to a string
    PyObject* sys_dict = PyModule_GetDict(sys);
    PyObject* io_module = PyImport_ImportModule("io");
    if (!io_module) {
        PyErr_Print();
        std::cerr << "Failed to import the 'io' module\n";
        Py_XDECREF(sys);
        Py_Finalize();
        return "";
    }
    PyObject* io_dict = PyModule_GetDict(io_module);
    PyObject* stdout = PyObject_GetAttrString(io_module, "StringIO");
    if (!stdout) {
        PyErr_Print();
        std::cerr << "Failed to access StringIO class\n";
        Py_XDECREF(sys);
        Py_XDECREF(io_module);
        Py_Finalize();
        return "";
    }
    PyObject* string_io = PyObject_CallObject(stdout, nullptr);
    Py_XDECREF(stdout);
    if (!string_io) {
        PyErr_Print();
        std::cerr << "Failed to create StringIO object\n";
        Py_XDECREF(sys);
        Py_XDECREF(io_module);
        Py_Finalize();
        return "";
    }
    PyObject_SetAttrString(sys, "stdout", string_io);
    Py_XDECREF(string_io);
    Py_XDECREF(io_module);

    // Execute the Python code
    PyObject* mainModule = PyImport_AddModule("__main__");
    PyObject* globalDict = PyModule_GetDict(mainModule);
    PyObject* localDict = PyDict_New();
    PyObject* result = PyRun_String(pythonCode.c_str(), Py_file_input, globalDict, localDict);
    if (PyErr_Occurred()) {
        PyErr_Print();
        std::cerr << "Failed to execute Python code\n";
        Py_XDECREF(localDict);
        Py_XDECREF(sys);
        Py_Finalize();
        return "";
    }
    Py_XDECREF(result);

    // Retrieve the captured stdout
    PyObject* output = PyObject_GetAttrString(sys, "stdout");
    PyObject* output_value = PyObject_CallMethod(output, "getvalue", nullptr);
    std::string resultStr = "";
    if (output_value != nullptr) {
        resultStr = PyUnicode_AsUTF8(output_value);
        Py_DECREF(output_value);
    }
    Py_XDECREF(output);

    // Clean up and finalize Python
    Py_XDECREF(localDict);
    Py_XDECREF(sys);
    Py_Finalize();

    return resultStr;
}

// Define the property of the graph edges (using "length" as the weight)
struct EdgeProperty {
    double weight;
};

// Define the property of the graph vertices
typedef boost::property<boost::vertex_name_t, std::string> VertexProperty;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexProperty, EdgeProperty> DirectedGraph;
typedef boost::graph_traits<DirectedGraph>::vertex_descriptor VertexDirected;
typedef boost::graph_traits<DirectedGraph>::edge_descriptor EdgeDirected;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty> UndirectedGraph;
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor VertexUndirected;
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor EdgeUndirected;

// Custom function to print the vertex properties
ostream& operator<<(ostream& os, const VertexProperty& vp) {
    os << vp;
    return os;
}

unordered_map<string, string> stopidToVertex(const std::string& filename) {
    
    unordered_map<std::string, string> stopMap;

    // Read data from CSV file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Error: Unable to open the CSV file: " << filename << std::endl;
    }

    // Assuming the CSV format is: stop_id,nearest_node_id
    std::string line;
    std::getline(file, line); // Skip header line

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string stop_id, nearest_node_id;

        std::getline(ss, stop_id, ',');
        ss >> nearest_node_id;
        // std::getline(ss, nearest_node_id, ',');

        // Add stop_ids if they don't exist already
        if (stopMap.find(stop_id) == stopMap.end()) {
            // Assuming VertexDirected has an attribute to store nearest_node_id
            stopMap[stop_id] = nearest_node_id;
        }
    }

    // Close the file after reading
    file.close();

    return stopMap;
}

void createDirectedGraph(const string& filename, DirectedGraph& g) {
    unordered_map<string, VertexDirected> vertexMap;

    // Read data from CSV file
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Error: Unable to open the CSV file: " << filename << endl;
        return;
    }

    // Assuming the CSV format is: Source, Target, OSM ID, Highway, Oneway, Reversed, Length
    string line;
    getline(file, line); // Skip header line

    while (getline(file, line)) {
        stringstream ss(line);
        string source, target, osmId, highway, oneway, reversed;
        double length;

        getline(ss, source, ',');
        getline(ss, target, ',');
        getline(ss, osmId, ',');
        getline(ss, highway, ',');
        getline(ss, oneway, ',');
        getline(ss, reversed, ',');
        ss >> length;

        // Add vertices if they don't exist already
        VertexDirected v_source, v_target;
        auto it_source = vertexMap.find(source);
        auto it_target = vertexMap.find(target);

        if (it_source == vertexMap.end()) {
            v_source = boost::add_vertex(VertexProperty(source), g);
            vertexMap[source] = v_source;
        } else {
            v_source = it_source->second;
        }

        if (it_target == vertexMap.end()) {
            v_target = boost::add_vertex(VertexProperty(target), g);
            vertexMap[target] = v_target;
        } else {
            v_target = it_target->second;
        }

        // Add edge with weight (length)
        EdgeDirected e;
        bool added;
        tie(e, added) = boost::add_edge(v_source, v_target, g);

        if (added) {
            g[e].weight = length;
        }
    }

    // Close the file after reading
    file.close();

    return;
}

void createUndirectedGraph(const string& filename, UndirectedGraph& g) {
    unordered_map<string, VertexUndirected> vertexMap;

    // Read data from CSV file
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Error: Unable to open the CSV file: " << filename << endl;
        return;
    }

    // Assuming the CSV format is: Source, Target, OSM ID, Highway, Oneway, Reversed, Length
    string line;
    getline(file, line); // Skip header line

    while (getline(file, line)) {
        stringstream ss(line);
        string source, target, osmId, highway, oneway, reversed;
        double length;

        getline(ss, source, ',');
        getline(ss, target, ',');
        getline(ss, osmId, ',');
        getline(ss, highway, ',');
        getline(ss, oneway, ',');
        getline(ss, reversed, ',');
        ss >> length;

        // Add vertices if they don't exist already
        VertexUndirected v_source, v_target;
        auto it_source = vertexMap.find(source);
        auto it_target = vertexMap.find(target);

        if (it_source == vertexMap.end()) {
            v_source = boost::add_vertex(VertexProperty(source), g);
            vertexMap[source] = v_source;
        } else {
            v_source = it_source->second;
        }

        if (it_target == vertexMap.end()) {
            v_target = boost::add_vertex(VertexProperty(target), g);
            vertexMap[target] = v_target;
        } else {
            v_target = it_target->second;
        }

        // Add edge with weight (length) for an undirected graph
        EdgeUndirected e1, e2;
        bool added1, added2;
        tie(e1, added1) = boost::add_edge(v_source, v_target, g);
        tie(e2, added2) = boost::add_edge(v_target, v_source, g);

        if (added1 && added2) {
            g[e1].weight = length;
            g[e2].weight = length;
        }
    }

    // Close the file after reading
    file.close();
}

double printShortestPath(const string& graphType, const string& source, const string& target, vector<string> &optimalPath) {
    // cout << "Calculating shortest path for " << graphType << " graph from " << source << " to " << target << ":\n";

    // Use different graph types and vertex descriptors based on directed/undirected graph
    if (graphType == "directed") {
        DirectedGraph graph;
        createDirectedGraph("digraph_details.csv", graph);

        boost::property_map<DirectedGraph, boost::vertex_name_t>::type vertexNameMapDirected = get(boost::vertex_name, graph);
        std::unordered_map<std::string, VertexDirected> vertexDescriptorMapDirected;

        // Populate the map to access vertex descriptor by name
        for (auto vd : boost::make_iterator_range(boost::vertices(graph))) {
            vertexDescriptorMapDirected[vertexNameMapDirected[vd]] = vd;
        }

        VertexDirected startVertex = vertexDescriptorMapDirected[source];
        VertexDirected endVertex = vertexDescriptorMapDirected[target];

        // Create a container to store the predecessor vertices
        std::vector<VertexDirected> predecessors(boost::num_vertices(graph));

        // Create a map to store the distance from the start vertex
        std::unordered_map<VertexDirected, double> distances;

        // Run Dijkstra's algorithm
        boost::dijkstra_shortest_paths(graph, startVertex,
            boost::predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(boost::vertex_index, graph)))
            .distance_map(boost::make_assoc_property_map(distances))
            .weight_map(boost::get(&EdgeProperty::weight, graph)) // Specify edge weight property map
        );

        // Check if there is a valid path from startVertex to endVertex
        if (distances[endVertex] < std::numeric_limits<double>::max()) {
            // Retrieve the shortest path from start to end
            std::vector<VertexDirected> path;
            for (VertexDirected v = endVertex; v != startVertex; v = predecessors[v]) {
                path.push_back(v);
            }
            path.push_back(startVertex);

            // Print the shortest path and distance
            // cout << "Shortest Path: ";
            for (auto it = path.rbegin(); it != path.rend(); ++it) {
                // cout << vertexNameMapDirected[*it];
                // if (it + 1 != path.rend()) {
                //     cout << ", ";
                // }
                optimalPath.push_back(vertexNameMapDirected[*it]);
                // cout<<vertexNameMapDirected[*it]<<" ";
            }
            // cout << "\nShortest Distance: " << distances[endVertex] << " units\n";
        } else {
            // Print a message indicating there is no path between the two points
            cout << "There is no valid path from " << source << " to " << target << ".\n";
        }

        return distances[endVertex];
    }
    else {
        UndirectedGraph graph;
        createUndirectedGraph("graph_details.csv", graph);

        boost::property_map<UndirectedGraph, boost::vertex_name_t>::type vertexNameMapUndirected = get(boost::vertex_name, graph);
        std::unordered_map<std::string, VertexUndirected> vertexDescriptorMapUndirected;

        // Populate the map to access vertex descriptor by name
        for (auto vd : boost::make_iterator_range(boost::vertices(graph))) {
            vertexDescriptorMapUndirected[vertexNameMapUndirected[vd]] = vd;
        }

        VertexUndirected startVertex = vertexDescriptorMapUndirected[source];
        VertexUndirected endVertex = vertexDescriptorMapUndirected[target];

        // Create a container to store the predecessor vertices
        std::vector<VertexUndirected> predecessors(boost::num_vertices(graph));

        // Create a map to store the distance from the start vertex
        std::unordered_map<VertexUndirected, double> distances;

        // Run Dijkstra's algorithm
        boost::dijkstra_shortest_paths(graph, startVertex,
            boost::predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(boost::vertex_index, graph)))
            .distance_map(boost::make_assoc_property_map(distances))
            .weight_map(boost::get(&EdgeProperty::weight, graph)));

        // Retrieve the shortest path from start to end
        std::vector<VertexUndirected> path;
        for (VertexUndirected v = endVertex; v != startVertex; v = predecessors[v]) {
            path.push_back(v);
        }
        path.push_back(startVertex);

        // If the endVertex is unreachable, distances[endVertex] will be the default value (infinity)
        if (distances[endVertex] != std::numeric_limits<double>::infinity()) {
            // Print the shortest path and distance
            for (auto it = path.rbegin(); it != path.rend(); ++it) {
                optimalPath.push_back(vertexNameMapUndirected[*it]);
            }
        }

        // // Print the shortest path and distance
        // cout << "Shortest Path: ";
        // for (auto it = path.rbegin(); it != path.rend(); ++it) {
        //     cout << vertexNameMapUndirected[*it];
        //     if (it + 1 != path.rend()) {
        //         cout << ", ";
        //     }
        // }
        // cout << "\nShortest Distance: " << distances[endVertex] << " units\n";

        return distances[endVertex];
    }
}

vector<string> nearestNodes(double max_distance_meters, double source_lat, double source_lon, double dest_lat, double dest_lon) {

    // Prepare the Python code with formatted arguments
    std::string pythonCodeTemplate = R"(
import osmnx as ox

# Create the graph for the specified bounding box
G = ox.graph.graph_from_bbox(45.5097368, 45.4997368, -73.5639306, -73.5739306, network_type='all', simplify=True, retain_all=False, truncate_by_edge=False, clean_periphery=True, custom_filter=None)

# Maximum distance (in meters) within which to find the nearest nodes
max_distance_meters = {}

# Find all the nodes in the graph
all_nodes = list(G.nodes())

#Find nearest node of the source_lat and source_lon AND dest_lat and dest_lon
source_lat = {}
source_lon = {}
dest_lat = {}
dest_lon = {}
source_nearest_node = ox.distance.nearest_nodes(G, source_lon, source_lat, return_dist=False)
dest_nearest_node = ox.distance.nearest_nodes(G, dest_lon, dest_lat, return_dist=False)

# Find the latitude and longitude of the given node
given_node_latitude = source_lat
given_node_longitude = source_lon

# Find the nearest nodes to the given node within the specified distance
nearest_nodes = []
for node in all_nodes:
    if node == source_nearest_node:
        continue  # Skip the given node itself
    distance = ox.distance.great_circle_vec(given_node_latitude, given_node_longitude, G.nodes[node]['y'], G.nodes[node]['x'])
    if distance <= max_distance_meters:
        nearest_nodes.append(node)
        
# Print the result
nearest_nodes.append(source_nearest_node)
nearest_nodes.append(dest_nearest_node)
print(nearest_nodes)
)";

    // Format the Python code with the provided arguments
    std::string pythonCode = pythonCodeTemplate;
    size_t distancePos = pythonCode.find("{}");
    if (distancePos != std::string::npos) {
        pythonCode.replace(distancePos, 2, std::to_string(max_distance_meters));
    }
    size_t source_lat_pos = pythonCode.find("{}");
    if (source_lat_pos != std::string::npos) {
        pythonCode.replace(source_lat_pos, 2, std::to_string(source_lat));
    }
    size_t source_lon_pos = pythonCode.find("{}");
    if (source_lon_pos != std::string::npos) {
        pythonCode.replace(source_lon_pos, 2, std::to_string(source_lon));
    }
    size_t dest_lat_pos = pythonCode.find("{}");
    if (dest_lat_pos != std::string::npos) {
        pythonCode.replace(dest_lat_pos, 2, std::to_string(dest_lat));
    }
    size_t dest_lon_pos = pythonCode.find("{}");
    if (dest_lon_pos != std::string::npos) {
        pythonCode.replace(dest_lon_pos, 2, std::to_string(dest_lon));
    }

    // Execute the Python code and get the result
    std::string res = executePythonCode(pythonCode);
    vector<string> result;
    string temp = "";

    for(int i = 1; i < res.size() - 2; i++){
        if(res[i] == ',') result.push_back(temp), temp = "", i++;
        else temp += res[i];
    }
    result.push_back(temp);

    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " max_distance_meters source_node_id dest_node_id dep_time\n";
        return 1;
    }

    // Get the command-line arguments for max_distance_meters and source_node_id
    double max_distance_meters = std::stod(argv[1]);
    double source_lat = std::stod(argv[2]);
    double source_lon = std::stod(argv[3]);
    double dest_lat = std::stod(argv[4]);
    double dest_lon = std::stod(argv[5]);
    string dep_time = argv[6];

    // Unordered map to store stop_id to stop_lat and stop_lon pair
    unordered_map<string, pair<double, double>> stopLatLonMap;

    // Unordered map to store stop_id to sorted list of pair of arrival_time and trip_id
    unordered_map<string, vector<pair<string, string>>> stopArrivalMap;

    // Unordered_map to store trip_id to vector of pair of stop_id and arrival_time
    unordered_map<string, vector<pair<string, string>>> tripDetails;

    ifstream inputFile("filtered_data.txt");
    if (!inputFile) {
        cerr << "Error: Unable to open the input file." << endl;
        return 1;
    }

    string line;
    while (getline(inputFile, line)) {
        stringstream ss(line);
        string trip_id, arrival_time, departure_time, stop_id, stop_sequence, stop_code, stop_name;
        double stop_lat, stop_lon;
        // Read CSV fields from the line
        if (getline(ss, trip_id, ',') &&
            getline(ss, arrival_time, ',') &&
            getline(ss, departure_time, ',') &&
            getline(ss, stop_id, ',') &&
            getline(ss, stop_sequence, ',') &&
            getline(ss, stop_code, ',') &&
            getline(ss, stop_name, ',') &&
            ss >> stop_lat &&
            ss.ignore() &&
            ss >> stop_lon) {

            // Store stop_id to stop_lat and stop_lon pair in the map
            stopLatLonMap[stop_id] = make_pair(stop_lat, stop_lon);
            stopArrivalMap[stop_id].push_back({arrival_time, trip_id});
            tripDetails[trip_id].push_back({stop_id, arrival_time});
        }
    }

    // Sort the vectors in the stopArrivalMap by arrival_time
    for (auto& entry : stopArrivalMap) {
        // First, sort the vector by arrival_time
        sort(entry.second.begin(), entry.second.end());
    }

    // // Print the contents of the maps (for demonstration purposes)
    // for (const auto& entry : stopLatLonMap) {
    //     cout << "Stop ID: " << entry.first << ", Lat: " << entry.second.first << ", Lon: " << entry.second.second << endl;
    // }

    // for (const auto& entry : stopArrivalMap) {
    //     cout << "Stop ID: " << entry.first << ", Arrivals: ";
    //     for (const auto& arrival : entry.second) {
    //         cout << "(" << arrival.first << ", " << arrival.second <<") ";
    //     }
    //     cout << endl;
    // }

    vector<string> nodes = nearestNodes(max_distance_meters, source_lat, source_lon, dest_lat, dest_lon);
    string source_nearest_node = nodes[nodes.size() - 2];
    string dest_nearest_node = nodes[nodes.size() - 1];
    cout<<source_nearest_node<<" "<<dest_nearest_node<<"\n";
    nodes.pop_back(), nodes.pop_back();

    DirectedGraph directedGraph;
    UndirectedGraph undirectedGraph;

    // Create directed graph from "digraph_details.csv"
    createDirectedGraph("digraph_details.csv", directedGraph);

    // Create undirected graph from "graph_details.csv"
    createUndirectedGraph("graph_details.csv", undirectedGraph);

    unordered_map<string, string> stopMap = stopidToVertex("nearest_nodes.csv");
    unordered_map<string, string> nodeMap;
    
    for (const auto& pair : stopMap) {
        const string& stop_id = pair.first;
        const string& nearest_node_id = pair.second;
        nodeMap[nearest_node_id] = stop_id;
        // // Print the stop_id and its corresponding nearest_node_id (assuming it's an attribute in VertexDirected)
        // std::cout << "Stop ID: " << stop_id << ", Nearest Node ID: " << nearest_node_id << std::endl;
    }

    // // Print the directed graph
    // cout << "Directed Graph:" << endl;
    // boost::property_map<DirectedGraph, boost::vertex_name_t>::type vertexNameMapDirected = get(boost::vertex_name, directedGraph);
    // boost::graph_traits<DirectedGraph>::vertex_iterator viDirected, vi_endDirected;
    // for (tie(viDirected, vi_endDirected) = boost::vertices(directedGraph); viDirected != vi_endDirected; ++viDirected) {
    //     cout << "Vertex: " << vertexNameMapDirected[*viDirected] << endl;
    // }

    // boost::graph_traits<DirectedGraph>::edge_iterator eiDirected, ei_endDirected;
    // for (tie(eiDirected, ei_endDirected) = boost::edges(directedGraph); eiDirected != ei_endDirected; ++eiDirected) {
    //     VertexDirected sourceDirected = boost::source(*eiDirected, directedGraph);
    //     VertexDirected targetDirected = boost::target(*eiDirected, directedGraph);
    //     cout << "Edge: " << vertexNameMapDirected[sourceDirected] << " -> " << vertexNameMapDirected[targetDirected] << ", Weight: " << directedGraph[*eiDirected].weight << endl;
    // }

    // // Print the undirected graph
    // cout << "\nUndirected Graph:" << endl;
    // boost::property_map<UndirectedGraph, boost::vertex_name_t>::type vertexNameMapUndirected = get(boost::vertex_name, undirectedGraph);
    // boost::graph_traits<UndirectedGraph>::vertex_iterator viUndirected, vi_endUndirected;
    // for (tie(viUndirected, vi_endUndirected) = boost::vertices(undirectedGraph); viUndirected != vi_endUndirected; ++viUndirected) {
    //     cout << "Vertex: " << vertexNameMapUndirected[*viUndirected] << endl;
    // }

    // boost::graph_traits<UndirectedGraph>::edge_iterator eiUndirected, ei_endUndirected;
    // for (tie(eiUndirected, ei_endUndirected) = boost::edges(undirectedGraph); eiUndirected != ei_endUndirected; ++eiUndirected) {
    //     VertexUndirected sourceUndirected = boost::source(*eiUndirected, undirectedGraph);
    //     VertexUndirected targetUndirected = boost::target(*eiUndirected, undirectedGraph);
    //     cout << "Edge: " << vertexNameMapUndirected[sourceUndirected] << " <-> " << vertexNameMapUndirected[targetUndirected] << ", Weight: " << undirectedGraph[*eiUndirected].weight << endl;
    // }

    // Calculate and print nodes within specified distance from a source vertex for undirected graph
    // string sourceVertexUndir = "1768772432"; // Replace with any source vertex
    // double distance = 200; // Replace with the desired distance

    // // Calculate and print shortest path for directed graph
    // string sourceVertexDir = "1768772432"; // Replace with any source vertex
    // string targetVertexDir = "32125080"; // Replace with any target vertex
    // printShortestPath("directed", sourceVertexDir, targetVertexDir);

    // // Calculate and print shortest path for undirected graph
    // string sourceVertexUndir = "1768772432"; // Replace with any source vertex
    // string targetVertexUndir = "32125080"; // Replace with any target vertex
    // printShortestPath("undirected", sourceVertexUndir, targetVertexUndir);

    vector<string> worstWalk;

    double walk_dist = printShortestPath("undirected", source_nearest_node, dest_nearest_node, worstWalk);
    int walk_time = walk_dist * 18 / 20, optimalTime; // 4 kmph walk speed

    for(string &s: nodes){
        vector<string> busPath, beforeWalkPath, afterWalkPath;
        if(nodeMap.find(s) == nodeMap.end()) continue;
        // cout<<s<<"\n";
        double initialWalk = printShortestPath("undirected", source_nearest_node, s, beforeWalkPath);
        // if(!initialWalk) continue;
        string stop_id = nodeMap[s];
        int trip = lower_bound(stopArrivalMap[stop_id].begin(), stopArrivalMap[stop_id].end(), make_pair(dep_time, ""), comparePairs) - stopArrivalMap[stop_id].begin();
        if(trip == stopArrivalMap[stop_id].size()) continue;
        // cout<<stopArrivalMap[stop_id][trip].first<<" "<<stopArrivalMap[stop_id][trip].second<<"\n";
        int differenceInSeconds = timeDifferenceInSeconds(dep_time, stopArrivalMap[stop_id][trip].first), time;
        if(differenceInSeconds > walk_time) continue;
        string trip_id = stopArrivalMap[stop_id][trip].second, busTime = stopArrivalMap[stop_id][trip].first;
        int time_taken = differenceInSeconds;
        for(int i = 0; i < tripDetails[trip_id].size() - 1; i++){
            pair<string, string> cur = tripDetails[trip_id][i], next = tripDetails[trip_id][i + 1];
            if(cur.second < busTime) continue;
            if(haversineDistance(stopLatLonMap[cur.first], dest_lat, dest_lon) < haversineDistance(stopLatLonMap[next.first], dest_lat, dest_lon)){
                // cout<<"Bus "<<i<<"\n";
                if(!busPath.empty()) busPath.pop_back();
                int bus_distance = printShortestPath("directed", stopMap[cur.first], stopMap[next.first], busPath);
                time_taken += timeDifferenceInSeconds(cur.second, next.second);
                // cout<<next.first<<" ";
            }
            else{
                int finalWalk = printShortestPath("undirected", stopMap[cur.first], dest_nearest_node, afterWalkPath);
                time_taken += finalWalk * 18 / 20;
                // cout<<"\n";
            }
        }

        if(busPath.size() == 0 || time_taken > walk_time) continue;

        cout<<"Initial walk: \n";
        for(string &f: beforeWalkPath) cout<<f<<", ";
        cout<<"\nBus path:\n";
        for(string &f: busPath) cout<<f<<", ";
        cout<<"\nAfter walk:\n";
        for(string &f: afterWalkPath) cout<<f<<", ";
        cout<<"\nTime taken: "<<time_taken<<"\n";
        cout<<"\n\n\n";
        //get the trip and then all subsequent stops and check for every further stop whether the stop is closer to destination or not and use flags.
    }

    return 0;
}
