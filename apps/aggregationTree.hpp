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
#include <set>


class AggregationTree {
public:

    AggregationTree();
    virtual ~AggregationTree(){};

    void init_labels();

    void update_labels();

    void add_to_tree(int x, int prev_iousx);

    void augment(std::vector<int> &allocation);

    int hungarian(std::vector<int> &allocation);

    std::vector<int> assignmentProblem(int Arr[], int N);

    // Comparator for priority queue
    struct Compare {
        bool operator()(const std::pair<std::string, int>& p1, const std::pair<std::string, int>& p2) {
            return p1.second > p2.second;
        }
    };

    std::vector <std::string> getContextInfo();

    std::vector <std::string> deleteNodes(std::vector <std::string> deletedList, std::vector <std::string> oldList);

    bool initializeGraph();

    int findLinkCost(const std::string& start, const std::string& end);

    bool compareClusters(const std::vector<std::vector<std::string>>& cluster1, const std::vector<std::vector<std::string>>& cluster2);

    std::vector<std::vector<std::string>> balancedKMeans(int N, int C, int numClusters, std::vector<int> clusterAssignment,
                                                         std::vector<std::string> dataPointNames, std::vector<std::vector<std::string>> clusters);

    std::string findCH(std::vector<std::string> clusterNodes, std::vector<std::string> clusterHeadCandidate, std::string client);

    bool aggregationTreeConstruction(std::vector<std::string> dataPointNames, int C);

    std::vector<std::string> getProducers();

    int countProducers();





    // Global variables
    std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> graph;
    std::string filename = "scratch/DataCenterTopology.txt";
    //std::string filename = "scratch/aggregationGrid.txt";
    std::vector<std::string> fullList;
    std::vector<std::string> CHList;
    std::string globalClient = "con0";
    std::map<std::string, std::vector<std::string>> aggregationAllocation;
    std::vector<std::vector<std::string>> noCHTree;

    // Hungarian variables
    int cost[1001][1001]; //cost matrix
    int n, max_match; //n workers and n jobs
    int lx[1001], ly[1001]; //labels of X and Y parts
    int xy[1001]; //xy[x] - vertex that is matched with x,
    int yx[1001]; //yx[y] - vertex that is matched with y
    bool S[1001], T[1001]; //sets S and T in algorithm
    int slack[1001]; //as in the algorithm description
    int slackx[1001]; //slackx[y] such a vertex, that
    int prev_ious[1001]; //array for memorizing alternating p

};