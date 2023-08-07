
#  MITACS GRI Project: Multimodal Transit Routing Algorithm

## Algorithm:

Currently, this algorithm finds all such paths wherein a person walks from his current location to nearby bus stops, takes a bus, and goes no further ahead than when he starts becoming farther from the destination and gets down at the bus stop nearest to the destination and then walks again to the target location.

Input: (Latitude, Longitude) of source and destination, preferredWalkDistance and departureTime

Output: A set of all paths that involve a walkable distance less than or equal to "preferredWalkDistance" and may or may not involve a bus in the travel path.

The basic concept is Dijkstra's algorithm, which has been further modified to accommodate combinations of multiple modes of travel.

Here assumption is that person walks at a speed of 4 kmph and bus travels at a speed of 20kmph or more. However, the time taken by the bus to travel from one bus stop to another is the difference between the arrival times of the bus at the 2 stops.

The optimization in this algorithm is happening in the time domain, meaning that the person might have to travel a longer distance than necessary in order to save time. For example, to travel 1 km, there might be a shorter walk path, but the bus might cover larger distance in shorter time due to higher speed, enabling the person to cover the distance in a shorter period of time.

### Terms: 

departureTime: This is the time that the person sets out to begin the travel

preferredWalkDistance: This refers to the maximum distance that the person is ready to walk in the multimodal transit path

worstPath: This is the path related to worst-case scenario wherein the person has to walk entirely between source and destination

printShortestPath (type of graph, source node, destination node, path): This is a function that finds the path between source and destination, with the undirected graph leading to walk path calculation and directed graph leading to bus path calculation

walk_time: This refers to the time taken to walk from source to destination entirely. This acts as an upper bound for all our optimisation in the time domain.

beforeWalkPath: This is the path of walk from source to the nearby bus stops.

afterWalkPath: This is the path of walk from bus stop to the target.

busPath: This refers to the intermediate path of bus between 2 bus stops closest to source and target seperately.

initialWalk: The distance traveled through "beforeWalkPath".

trip_id: Each bus follows a trajectory and and this trajectory itself is trip, which is essentially a vector of pair of stop_id and time of arrival of the bus at that stop.

nodes: vector of all nodes that are at a distance less than or equal to "preferredWalkDistance"

trip: Refers to the trip of the first bus that comes at the bus stop after person has walked up to the bus stop.

nodeMap: Map that maps the nearest node of a stop, to the stop_id of that stop

stopArrivalMap: Map that maps a stop_id to a vector of pair of arrival time and trip id of the buses that visit a bus stop

busTime: This refers to the first bus that arrives at the bus stop after departure time

time_taken: This refers to the total time taken to travel in a path, in seconds.

tripDetails: Map that maps the trip id to a vector of pair of stop id and arrival time of the bus at that stop

finalWalk: Distance travelled from the bus stop to the destination

### Pseudocode:

    vector<string> worstPath;

    double walk_dist = printShortestPath("undirected", source_nearest_node, dest_nearest_node, worstPath);
    int walk_time = walk_dist * 18 / 20; // 4 kmph walk speed 4 * 5 / 18 = 20 / 18 (20 / 18 = dist / time)

    for(string &s: nodes){
        vector<string> busPath, beforeWalkPath, afterWalkPath;

        //if the node is not a bus stop, then skip it.
        if(nodeMap.find(s) == nodeMap.end()) continue;

        //accounts for the distance travelled from source to the nearby bus stops
        double initialWalk = printShortestPath("undirected", source_nearest_node, s, beforeWalkPath);

        //skip the bus stop if distance to be travelled is 0
        if(!initialWalk) continue;
        
        string stop_id = nodeMap[s];

        //This basically finds the first bus to reach the bus stop after departure time
        int trip = lower_bound(stopArrivalMap[stop_id].begin(), stopArrivalMap[stop_id].end(), make_pair(dep_time, ""), comparePairs) - stopArrivalMap[stop_id].begin();

        //If there is no bus arriving at the bus stop after departure time, then skip this bus stop.
        if(trip == stopArrivalMap[stop_id].size()) continue;

        //calculation of time taken to walk from source node to the current bus stop
        int differenceInSeconds = timeDifferenceInSeconds(dep_time, stopArrivalMap[stop_id][trip].first), time;

        //if the time taken to walk upto the bus stop is greater than the time taken to walk upto the destination node, then skip such bus stops
        if(differenceInSeconds > walk_time) continue;

        
        string trip_id = stopArrivalMap[stop_id][trip].second, busTime = stopArrivalMap[stop_id][trip].first;
        int time_taken = differenceInSeconds;

        //track a bus's path and check how far the person can travel using the bus to be nearest to the destination
        for(int i = 0; i < tripDetails[trip_id].size() - 1; i++){

            //keeping track of current and next bus stops
            pair<string, string> cur = tripDetails[trip_id][i], next = tripDetails[trip_id][i + 1];

            //skipping all buses that visited the bus stop before the departure time
            if(cur.second < busTime) continue;

            //if displacement between current bus stop and destination is greater than displacement between next bus stop and destination, then go to the next bus stop, else get down from the bus.
            if(haversineDistance(stopLatLonMap[cur.first], dest_lat, dest_lon) < haversineDistance(stopLatLonMap[next.first], dest_lat, dest_lon)){
                if(!busPath.empty()) busPath.pop_back();
                printShortestPath("directed", stopMap[cur.first], stopMap[next.first], busPath);
                time_taken += timeDifferenceInSeconds(cur.second, next.second);
            }
            else{
                int finalWalk = printShortestPath("undirected", stopMap[cur.first], dest_nearest_node, afterWalkPath);
                time_taken += finalWalk * 18 / 20;
            }
        }

        //If time taken is greater than that spent in walking from source to destination, then skip it.
        if(time_taken > walk_time) continue;

        std::cout<<"Initial walk: \n";
        for(string &f: beforeWalkPath) std::cout<<f<<", ";
        std::cout<<"\nBus path:\n";
        for(string &f: busPath) std::cout<<f<<", ";
        std::cout<<"\nAfter walk:\n";
        for(string &f: afterWalkPath) std::cout<<f<<", ";
        std::cout<<"\nTime taken: "<<time_taken<<"\n";
    }
