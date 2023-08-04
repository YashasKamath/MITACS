### MITACS GRI Project on Developing Multimodal Transit Routing Algorithm, a combo of walk and bus modes of travel with a focus on optimising time taken to reach target from current location

# Algorithm:

About the algorithm: Currently, this algorithm finds all such paths wherein a person walks from his current location to nearby bus stops, takes a bus, and goes no further ahead than when he starts becoming farther from the destination and gets down at the bus stop nearest to the destination and then walks again to the target location.

Terms: 
worstPath: This is the path related to worst-case scenario wherein the person has to walk entirely between source and destination
printShortestPath (type of graph, source node, destination node, path): This is a function that finds the path between source and destination, with the undirected graph leading to walk path calculation and directed graph leading to bus path calculation
walk_time: This refers to the time taken to walk from source to destination entirely. This acts as an upper bound for all our optimisation in the time domain.
beforeWalkPath: This is the path of walk from source to the nearby bus stops.
afterWalkPath: This is the path of walk from bus stop to the target.
busPath: This refers to the intermediate path of bus between 2 bus stops closest to source and target seperately.
initialWalk: The distance traveled through "beforeWalkPath".
trip_id: Each bus follows a trajectory and and this trajectory itself is trip, which is essentially a vector of pair of stop_id and time of arrival of the bus at that stop.

    vector<string> worstPath;

    double walk_dist = printShortestPath("undirected", source_nearest_node, dest_nearest_node, worstPath);
    int walk_time = walk_dist * 18 / 20; // 4 kmph walk speed 4 * 5 / 18 = 20 / 18 (20 / 18 = dist / time)
    
    // for(string w: worstPath) cout<<w<<", ";

    for(string &s: nodes){
        vector<string> busPath, beforeWalkPath, afterWalkPath;
        if(nodeMap.find(s) == nodeMap.end()) continue;
        double initialWalk = printShortestPath("undirected", source_nearest_node, s, beforeWalkPath);
        if(!initialWalk) continue;
        // cout<<"Consider:\n";
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
                printShortestPath("directed", stopMap[cur.first], stopMap[next.first], busPath);
                time_taken += timeDifferenceInSeconds(cur.second, next.second);
                // cout<<next.first<<" ";
            }
            else{
                int finalWalk = printShortestPath("undirected", stopMap[cur.first], dest_nearest_node, afterWalkPath);
                time_taken += finalWalk * 18 / 20;
                std::cout<<"\n";
            }
        }
        // if(busPath.size() == 0 || time_taken > walk_time) continue;

        std::cout<<"Initial walk: \n";
        for(string &f: beforeWalkPath) std::cout<<f<<", ";
        std::cout<<"\nBus path:\n";
        for(string &f: busPath) std::cout<<f<<", ";
        std::cout<<"\nAfter walk:\n";
        for(string &f: afterWalkPath) std::cout<<f<<", ";
        std::cout<<"\nTime taken: "<<time_taken<<"\n";
        std::cout<<"\n\n\n";
        std::cout<<"\n\n\n";
        // //get the trip and then all subsequent stops and check for every further stop whether the stop is closer to destination or not and use flags.
    }
