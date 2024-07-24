#include <iostream>
#include <vector>
#include <cmath> // For ceil
#include <numeric> // For std::accumulate
#include <cassert> // For assert
#include<bits/stdc++.h>
#include <fstream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <limits>
#include <unordered_map>
#include <climits>

#include "aggregationTree.hpp"



AggregationTree::AggregationTree(){
    fullList = getContextInfo();
    CHList = fullList;
    if (!initializeGraph()) {
        std::cout << "Error when reading network topology" << std::endl;
    }
    //std::cout << "Finish initialization!" << std::endl;
}


void AggregationTree::init_labels()
{
    memset(lx, 0, sizeof(lx)); // initiate lx vector to 0
    memset(ly, 0, sizeof(ly)); // initiate ly vector to 0
    for (int x = 0; x < n; x++)
        for (int y = 0; y < n; y++)
            lx[x] = std::max(lx[x], cost[x][y]);
}


void AggregationTree::update_labels()
{
    int x, y;
    int delta = 99999999; //init delta as infinity
    for (y = 0; y < n; y++) //calculate delta using slack
        if (!T[y])
            delta = std::min(delta, slack[y]);
    for (x = 0; x < n; x++) //update X labels
        if (S[x])
            lx[x] -= delta;
    for (y = 0; y < n; y++) //update Y labels
        if (T[y])
            ly[y] += delta;
    for (y = 0; y < n; y++) //update slack array
        if (!T[y])
            slack[y] -= delta;
}


void AggregationTree::add_to_tree(int x, int prev_iousx)
//x - current vertex,prev_iousx - vertex from X before x in the alternating path, //so we add edges (prev_iousx, xy[x]), (xy[x], x)
{
    S[x] = true; //add x to S
    // prev_ious[x] = prev_iousx; //we need this when augmenting
    for (int y = 0; y < n; y++) //update slacks, because we add new vertex to S
        if (lx[x] + ly[y] - cost[x][y] < slack[y]) {
            slack[y] = lx[x] + ly[y] - cost[x][y];
            slackx[y] = x;
        }
}


void AggregationTree::augment(std::vector<int> &allocation) //main function of the algorithm
{
    if (max_match == n) return; //check whether matching is already perfect
    int x, y, root; //just counters and root vertex
    int q[31], wr = 0, rd = 0; //q - queue for bfs, wr,rd - write and read
    //pos in queue
    memset(S, false, sizeof(S)); //init set S to false
    memset(T, false, sizeof(T)); //init set T to false
    memset(prev_ious, -1, sizeof(prev_ious)); //init set prev_ious - for the alternating tree

    for (x = 0; x < n; x++) //finding root of the tree
    {
        if (xy[x] == -1) {
            q[wr++] = root = x;
            prev_ious[x] = -2;
            S[x] = true;
            break;
        }
    }

    for (y = 0; y < n; y++) //initializing slack array
    {
        slack[y] = lx[root] + ly[y] - cost[root][y];
        slackx[y] = root;
    }

    //second part of augment() function
    while (true) //main cycle
    {
        while (rd < wr) //building tree with bfs cycle
        {
            x = q[rd++]; //current vertex from X part
            for (y = 0; y < n; y++) //iterate through all edges in equality graph
                if (cost[x][y] == lx[x] + ly[y] && !T[y]) {
                    if (yx[y] == -1) break; //an exposed vertex in Y found, so
                    //augmenting path exists!
                    T[y] = true; //else just add y to T,
                    q[wr++] = yx[y]; //add vertex yx[y], which is matched
                    //with y, to the queue
                    add_to_tree(yx[y], x); //add edges (x,y) and (y,yx[y]) to the tree
                }
            if (y < n)
                break; //augmenting path found!
        }
        if (y < n)
            break; //augmenting path found!

        update_labels(); //augmenting path not found, so improve labeling
        wr = rd = 0;
        for (y = 0; y < n; y++)
            //in this cycle we add edges that were added to the equality graph as a
            //result of improving the labeling, we add edge (slackx[y], y) to the tree if
            //and only if !T[y] && slack[y] == 0, also with this edge we add another one
            //(y, yx[y]) or augment the matching, if y was exposed
            if (!T[y] && slack[y] == 0) {
                if (yx[y] == -1) //exposed vertex in Y found - augmenting path exists!
                {
                    x = slackx[y];
                    break;
                } else {
                    T[y] = true; //else just add y to T,
                    if (!S[yx[y]]) {
                        q[wr++] = yx[y]; //add vertex yx[y], which is matched with
                        //y, to the queue
                        add_to_tree(yx[y], slackx[y]); //and add edges (x,y) and (y,
                        //yx[y]) to the tree
                    }
                }
            }
        if (y < n) break; //augmenting path found!
    }

    if (y < n) //we found augmenting path!
    {
        max_match++; //increment matching
        //in this cycle we inverse edges along augmenting path
        for (int cx = x, cy = y, ty; cx != -2; cx = prev_ious[cx], cy = ty) {
            ty = xy[cx];
            yx[cy] = cx;
            xy[cx] = cy;
        }
        augment(allocation); //recall function, go to step 1 of the algorithm
    }
}//end of augment() function

int AggregationTree::hungarian(std::vector<int> &allocation) {

    int ret = 0; //weight of the optimal matching
    max_match = 0; //number of vertices in current matching
    memset(xy, -1, sizeof(xy));
    memset(yx, -1, sizeof(yx));
    init_labels(); //step 0
    augment(allocation); //steps 1-3

    for (int x = 0; x < n; x++) {
        ret += cost[x][xy[x]];
        allocation[x * n + xy[x]] = 1;
    }//forming answer there

    return ret;
}

std::vector<int> AggregationTree::assignmentProblem(int Arr[], int N)
{
    n = N;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cost[i][j] = -1 * Arr[i * n + j];

    std::vector<int> allocation(N * N, 0);
    hungarian(allocation);
    return allocation;
}
// End of Hungarian algorithm



// Start of node list management, used for aggregation tree construction
// Return nodes except "con" and "forwarder"
std::vector <std::string> AggregationTree::getContextInfo() {
    std::ifstream file(filename);
    std::vector <std::string> nodes;

    if (!file.is_open()) {
        std::cerr << "Fail to open file: " << filename << std::endl;
        return nodes;
    }

    std::string line;
    bool linkSection = false;

    while (getline(file, line)) {
        if (line == "router") {
            linkSection = true;
            continue;
        } else if (line == "link") {
            break;
        }

        if (linkSection && !line.empty() && line.find("forwarder") != 0 && line.find("con") != 0) {
            std::istringstream iss(line);
            std::string nodeName;
            iss >> nodeName;
            nodes.push_back(nodeName);
        }
    }

    return nodes;
}

std::vector <std::string> AggregationTree::deleteNodes(std::vector <std::string> deletedList, std::vector <std::string> oldList) {
    std::vector <std::string> newList = oldList;

    auto newEnd = std::remove_if(newList.begin(), newList.end(),
                                 [&deletedList](const std::string &element) {
                                     return std::find(deletedList.begin(), deletedList.end(), element) !=
                                            deletedList.end();
                                 });

    newList.erase(newEnd, newList.end());

    return newList;
}

// End of node list management



// Start of link cost computation, used for all cases
// Define a custom type for easier readability
//typedef std::pair<std::string, std::string> NodePair;

// Function to initialize the adjacency list from the file
bool AggregationTree::initializeGraph() {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    bool linkSection = false;

    while (getline(file, line)) {
        if (line == "link") {
            linkSection = true;
            continue;
        }

        if (linkSection && !line.empty()) {
            std::istringstream iss(line);
            std::string node1, node2;
            std::string speed; // Placeholder for speed
            int cost;
            std::string delay; // Placeholder for delay
            int priority; // Placeholder for priority

            iss >> node1 >> node2 >> speed >> cost >> delay >> priority;

            // Store the connections in an adjacency list
            graph[node1].push_back(std::make_pair(node2, cost));
            graph[node2].push_back(std::make_pair(node1, cost));
        }
    }

    file.close();
    return true;
}


// Function to find the minimum link cost between any two nodes using Dijkstra's algorithm
int AggregationTree::findLinkCost(const std::string& start, const std::string& end) {
    if (graph.find(start) == graph.end() || graph.find(end) == graph.end())
        return -1; // Nodes are not present in the graph

    std::priority_queue<std::pair<std::string, int>, std::vector<std::pair<std::string, int>>, Compare> pq;
    std::unordered_map<std::string, int> distances;

    // Initialize distances to maximum
    for (const auto& node : graph) {
        distances[node.first] = INT_MAX;
    }

    // Start from the source node
    pq.push({start, 0});
    distances[start] = 0;

    while (!pq.empty()) {
        auto current = pq.top();
        pq.pop();

        std::string currentNode = current.first;
        int currentCost = current.second;

        // Shortest path to end node found
        if (currentNode == end) {
            return currentCost;
        }

        // Traverse all adjacent nodes
        for (const auto& neighbor : graph[currentNode]) {
            std::string nextNode = neighbor.first;
            int nextCost = neighbor.second;
            int newCost = currentCost + nextCost;

            // Check if a cheaper path is found
            if (newCost < distances[nextNode]) {
                distances[nextNode] = newCost;
                pq.push({nextNode, newCost});
            }
        }
    }

    std::cout << "Error happened, no route is found!!!!!!!!!!!" << std::endl;
    return -1; // No path found
}

// End of link cost calculation

// Compare two clusters for balanced k-means algorithm
bool AggregationTree::compareClusters(const std::vector<std::vector<std::string>>& cluster1, const std::vector<std::vector<std::string>>& cluster2) {
    if (cluster1.size() != cluster2.size()) {
        return false; // Different number of clusters
    }

    // Iterate through each cluster
    for (size_t i = 0; i < cluster1.size(); ++i) {
        // Create copies of the inner vectors
        std::vector<std::string> cluster1_inner = cluster1[i];
        std::vector<std::string> cluster2_inner = cluster2[i];

        // Sort both inner vectors
        std::sort(cluster1_inner.begin(), cluster1_inner.end());
        std::sort(cluster2_inner.begin(), cluster2_inner.end());

        // Compare sorted vectors
        if (cluster1_inner != cluster2_inner) {
            return false; // Different elements in at least one corresponding pair of clusters
        }
    }
    return true; // All corresponding pairs of clusters contain the same elements
}


std::vector<std::vector<std::string>> AggregationTree::balancedKMeans(int N, int C, int numClusters, std::vector<int> clusterAssignment,
                                                                      std::vector<std::string> dataPointNames, std::vector<std::vector<std::string>> clusters) {
    // Vector to hold the output matrix in 1D format
    std::vector<int> output(N * N);

    // Calculate the matrix
    for (int i = 0; i < N; ++i) {
        int dataPoint = i;
        int clusterOfDataPoint = clusterAssignment[dataPoint];

        // "clusterNodes" stores the data points this cluster has
        const auto& clusterNodes = clusters[clusterOfDataPoint];

        for (int j = 0; j < N; ++j) {
            std::string otherDataPoint = dataPointNames[j];
            long long totalCost = 0;

            // Compute total cost between current data point and all nodes in this cluster
            for (const auto& node : clusterNodes) {
                totalCost += findLinkCost(otherDataPoint, node);
            }

            // Compute average cost and store it into output matrix
            int averageCost = static_cast<int>(totalCost / clusterNodes.size());
            output[i * N + j] = averageCost;
        }
    }

    // Output the matrix in 1D vector format
    std::cout << "\nCost Matrix:" << std::endl;
    for (int i = 0; i < N * N; ++i) {
        if (i % N == 0) std::cout << std::endl;  // New line for each row in visualization
        std::cout << output[i] << " ";
    }
    std::cout << std::endl;

    // Start of Hungarian algorithm
    // Cast data type
    const size_t MAX_SIZE = output.size();
    int Hungarian_input[MAX_SIZE];
    std::copy(output.begin(), output.end(), Hungarian_input);

    std::vector<int> allocation = assignmentProblem(Hungarian_input, N);
    std::cout << "\nAllocation Matrix:" << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << allocation[i * N + j] << " ";
        }
        std::cout << std::endl;
    }

    // Update the cluster for next iteration
    std::vector<std::vector<std::string>> newCluster(numClusters);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            if (allocation[i * N + j] == 1){
                newCluster[clusterAssignment[i]].push_back(dataPointNames[j]);
            }
        }
    }

    if (compareClusters(clusters, newCluster)){
        std::cout << "Hungarian algorithm converges." << std::endl;
        return newCluster;
    } else {

        // testing old cluster
        int i = 0;
        std::cout << "\nOld cluster" << std::endl;
        for (const auto& oldItemVector : clusters) {
            std::cout << "Cluster " << i << " " <<std::endl;
            for (const auto& oldStr : oldItemVector) {
                std::cout << oldStr << " ";
            }
            std::cout << std::endl;
            ++i;
        }

        // testing new cluster
        int j = 0;
        std::cout << "\nNew cluster" << std::endl;
        for (const auto& oldItemVector : newCluster) {
            std::cout << "Cluster " << j << " " <<std::endl;
            for (const auto& oldStr : oldItemVector) {
                std::cout << oldStr << " ";
            }
            std::cout << std::endl;
            ++j;
        }

    }
    balancedKMeans(N, C, numClusters, clusterAssignment, dataPointNames, clusters);
}



// Start of cluster head construction

std::string AggregationTree::findCH(std::vector<std::string> clusterNodes, std::vector<std::string> clusterHeadCandidate, std::string client) {

    std::string CH = client;
    // Initiate a large enough cost
    int leastCost = 1000;

    for (const auto& headCandidate : clusterHeadCandidate) {
        bool canBeCH = true;

        for (const auto& node : clusterNodes) {
            if (findLinkCost(node, client) < findLinkCost(node, headCandidate)) {
                canBeCH = false;
                break;
            }
        }

        // This candidate is closer to client
        long long totalCost = 0;
        if (canBeCH) {
            for (const auto& node : clusterNodes) {
                totalCost += findLinkCost(node, headCandidate);
            }
            int averageCost = static_cast<int>(totalCost / clusterNodes.size());

            if (averageCost < leastCost) {
                leastCost = averageCost;
                CH = headCandidate;
            }
        }
    }

    if (CH == client) {
        std::cerr << "No CH is found for current cluster!!!!!!!!!!!" << std::endl;
        return CH;
    } else {
        std::cout << "CH " << CH << " is chosen." << std::endl;
        return CH;
    }


}

// End of cluster head construction

bool AggregationTree::aggregationTreeConstruction(std::vector<std::string> dataPointNames, int C) {

    // Compute N
    int N = dataPointNames.size();

    // Compute the number of clusters k
    int numClusters = static_cast<int>(ceil(static_cast<double>(N) / C));

    // Create a vector to store cluster assignments
    std::vector<int> clusterAssignment(N);
    for (int i = 0; i < N; ++i) {
        clusterAssignment[i] = i % numClusters; // Cluster assignment based on modulus operation
    }

    // Output the cluster assignments
    std::cout << "Cluster initialization." << std::endl;
    std::cout << "There are " << numClusters << " clusters." << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "Data point " << dataPointNames[i] << " is in cluster " << clusterAssignment[i] << std::endl;
    }

    // Create a map of clusters to their data points, store data point's ID inside cluster's vector
    std::vector<std::vector<std::string>> clusters(numClusters);
    for (int i = 0; i < N; ++i) {
        clusters[clusterAssignment[i]].push_back(dataPointNames[i]);
    }

    // Start of balanced K-Means
    // Get the output at current layer (data point allocation for each cluster)
    std::vector<std::vector<std::string>> newCluster = balancedKMeans(N, C, numClusters, clusterAssignment, dataPointNames, clusters);

    // Construct the nodeList for CH allocation
    int i = 0;
    std::cout << "\nIterating new clusters." << std::endl;
    for (const auto& iteCluster: newCluster) {
        std::cout << "Cluster " << i << " contains the following nodes:" <<std::endl;
        for (const auto& iteNode: iteCluster) {
            CHList.erase(std::remove(CHList.begin(), CHList.end(), iteNode), CHList.end());
            std::cout << iteNode << " ";
        }
        std::cout << std::endl;
        ++i;
    }


    std::cout << "\nCurrent CHList for CH allocation later: " << std::endl;
    for (const auto& item : CHList) {
        std::cout << item << std::endl;
    }

    // Start CH allocation, currently ignore those can't find cluster head
    std::vector<std::string> newDataPoints;
    std::cout << "\nStarting CH allocation." << std::endl;
    for (const auto& clusterNodes : newCluster) {
        std::string clusterHead = findCH(clusterNodes, CHList, globalClient);
        if (clusterHead != globalClient) {
            CHList.erase(std::remove(CHList.begin(), CHList.end(), clusterHead), CHList.end());
            aggregationAllocation[clusterHead] = clusterNodes;
            newDataPoints.push_back(clusterHead);
        }
        else {
            std::cout << "Due to no cluster head found, combine these nodes into sub-tree." << std::endl;
            noCHTree.push_back(clusterNodes);
        }

    }

    std::cout << "\nCHList after CH allocation: " << std::endl;
    for (const auto& item : CHList) {
        std::cout << item << std::endl;
    }

    if (newDataPoints.size() < C){
        // If with an entire layer, no CH found, then extract one from subTree list
        if (newDataPoints.empty()){
            const auto& firstSubTree = noCHTree[0];
            aggregationAllocation[globalClient] = firstSubTree;
            noCHTree.erase(noCHTree.begin());
            return true;
        } else {
            aggregationAllocation[globalClient] = newDataPoints;
            return true;
        }
    } else {
        aggregationTreeConstruction(newDataPoints, C);
    }
}

// Initialize the nodes from the bottom for first iteration, i.e. get all producers
std::vector<std::string> AggregationTree::getProducers() {
    std::ifstream file(filename);
    std::vector<std::string> proNodes;
    std::string line;
    bool inRouterSection = false;

    if (!file) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return proNodes; // Return empty if file cannot be opened
    }

    while (std::getline(file, line)) {
        // Trim whitespace and check if the line is a section header
        line.erase(0, line.find_first_not_of(" \t\n\r\f\v")); // Trim leading whitespace
        if (line == "router") {
            inRouterSection = true;
        } else if (line == "link") {
            break; // Exit if we reach the "link" section
        } else if (inRouterSection) {
            // Check if the line starts with "pro"
            if (line.substr(0, 3) == "pro") {
                proNodes.push_back(line);
            }
        }
    }
    file.close();
    return proNodes;

}

// Compute the number of producers, used for compute model average on consumer
int AggregationTree::countProducers() {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return -1; // Return -1 or any specific error code to indicate failure
    }

    std::string line;
    bool inRouterSection = false;
    int producerCount = 0;

    // Read the file content line by line
    while (std::getline(file, line)) {
        // Trim the line from both sides (optional, depends on input cleanliness)
        line.erase(line.find_last_not_of(" \n\r\t") + 1);
        line.erase(0, line.find_first_not_of(" \n\r\t"));

        // Check for the starting point of the router section
        if (line == "router") {
            inRouterSection = true;
            continue;
        }

        // Check for the end of the router section (only "link" marks the end)
        if (line == "link") {
            break;  // Exit the loop as we are past the relevant router section
        }

        // Count the producer lines if we are within the router section
        if (inRouterSection && line.substr(0, 3) == "pro") {
            producerCount++;
        }
    }

    file.close(); // Close the file after reading
    return producerCount;
}


