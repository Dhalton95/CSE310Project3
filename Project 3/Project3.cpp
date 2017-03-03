//
//  Project3.cpp
//  Project 3
//
//  Created by Dhalton Huber on 11/6/16.
//  Copyright © 2016 Dhalton Huber. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;

#define HEIGHT 63
#define WIDTH 63
#define DEPTH 852

struct GridNode
{
    float summation;
};


// Function to compute S-sub-x-x for Pearson Correlation Coefficient
float Sxx(float ***array, int x, int y)
{
    // All float viariables used in calculation
    float toReturn;
    float summation = 0.0;
    float xBar;
    float sum = 0.0;
    int divisor = DEPTH;
    
    
    // Calculate the sum for time series to compute xBar
    for (int a = 0; a <= DEPTH; a++)
    {
        if (array[x][y][a] != 168 && array[x][y][a] != 157)
        {
            sum += array[x][y][a];
        }
        else
        {
            divisor--;
        }

    }
    
    // xBar is the sum divided by the number of elements
    xBar = sum / divisor;
    
    // For loop to calcualte the summation
    for (int i = 0; i < DEPTH; i++)
    {
        if (array[x][y][i] != 168 && array[x][y][i] != 157)
        {
            // summation (x-sub-i - xBar)^2
            summation += ((array[x][y][i] - xBar) * (array[x][y][i] - xBar));
        }
        
    }
    // Set toReturn and return it
    toReturn = summation;
    return toReturn;
}


// Function to calculate S-sub-x-y for Pearson Correlation Coefficient
float Sxy(float ***array, int x, int y, int xPrime, int yPrime)
{
    // All float variables used in calculations
    float toReturn;
    float summation = 0.0;
    float xBar;
    float yBar;
    float xSum = 0.0;
    float ySum = 0.0;
    int xDivisor = DEPTH;
    int yDivisor = DEPTH;
    
    // Calculate the sum for time series to compute x an y bar
    for (int a = 0; a <= DEPTH; a++)
    {
        if (array[x][y][a] != 168 && array[x][y][a] != 157)
        {
            xSum += array[x][y][a];
        }
        else
        {
            xDivisor--;
        }
        
        if (array[xPrime][yPrime][a] != 168 && array[xPrime][yPrime][a] != 157)
        {
            ySum += array[xPrime][yPrime][a];
        }
        else
        {
            yDivisor--;
        }
        
    }
    
    // Compute x and y bar
    xBar = xSum / xDivisor;
    yBar = ySum / yDivisor;
    
    // For loop to calculate the summation
    for (int i = 0; i < DEPTH; i++)
    {
        if (array[x][y][i] != 168 && array[x][y][i] != 157 && array[xPrime][yPrime][i] != 168 && array[xPrime][yPrime][i] != 157)
        {
            // Summation (x-sub-i - xBar)(y-sub-i - yBar)
            summation += ((array[x][y][i] - xBar) * (array[xPrime][yPrime][i] - yBar));
        }

    }
    // Set toReturn and return it
    toReturn = summation;
    return toReturn;
}

// Function to compute r for Pearson Correlation Coefficient
float r(float ***array, int x, int y, int xPrime, int yPrime)
{
    float toReturn = 0.0;
    
    toReturn = Sxy(array, x, y, xPrime, yPrime) / sqrt(Sxx(array, x, y)*Sxx(array, xPrime, yPrime));
    
    return toReturn;
}





// Node structure to store edges
struct AListNode
{
    int x;
    int y;
    struct AListNode* next;
    int cellNumber;
};

// Node structure to create a grid to hold edges
struct AList
{
    struct AListNode *head;
    
    bool discovered;
    int discoveryTime;
    int finishTime;
    bool visited;
    int cellNumber;
};


class Graph
{
    private:
    int V;
    int time;
    int cellNumberCounter = 1;

    public:
    AList **graph2D;
    int connectedComponents;
    
        Graph(int V, float ***graph3D)
        {
            this->V = V;
            // Allocate space for our 2D graph of AList structures
            graph2D = new AList*[WIDTH];
            for (int i = 0; i < WIDTH; i++)
            {
                graph2D[i] = new AList[HEIGHT];
            }
            
            // Set all pointers in grid to NULL to begin with
            for (int i = 0; i < HEIGHT; i++)
            {
                for (int j = 0; j < WIDTH; j++)
                {
                    graph2D[j][i].head = NULL;
                    graph2D[j][i].visited = false;
                    graph2D[j][i].cellNumber = cellNumberCounter;
                    cellNumberCounter++;
                }
                
            }
        }
    
    ~Graph(void)
    {
    }
    
    // Function to add an edge to the graph between two cells (at x,y and xPrime, yPrime)
    void addEdge(int x, int y, int xPrime, int yPrime)
    {
        // Create 2 new nodes with x, y and xPrime, yPrime
        AListNode* newNode1 = newAListNode(xPrime, yPrime);
        AListNode* newNode2 = newAListNode(x, y);
        
        
        // Insert the nodes at their proper locations at the head of the linked list
        newNode1->next = graph2D[x][y].head;
        newNode2->next = graph2D[xPrime][yPrime].head;
        graph2D[x][y].head = newNode1;
        graph2D[xPrime][yPrime].head = newNode2;
    }
    
    // Function to create and initialize a new AListNode
    AListNode* newAListNode(int x, int y)
    {
        // Create a new node to return
        AListNode* newNode = new AListNode;
        // Set its x and y values
        newNode->x = x;
        newNode->y = y;
        // Set its next pointer to NULL
        newNode->next = NULL;
        newNode->cellNumber = graph2D[x][y].cellNumber;
        // Return
        return newNode;
    }
    
    // Function to get the degree for a cell
    int ComputeDegree(int x, int y)
    {
        int toReturn = 0;
        // Count used to count up the number of edges for a cell
        int count = 0;
        
        
        
            // Create an AListNode object and set it equal to the head of the LL at the cell
            AListNode *traverser = graph2D[x][y].head;
        
            // Loop through the linked list and count up all nodes
            while (traverser != NULL)
            {
                count++;
                traverser = traverser->next;
            }
            
        
        // Set toReturn to the count and return it
        toReturn = count;
        return toReturn;
    }
    
    
    // Function to perform a Depth-First Search on the graph
    void DFS()
    {
        // Set the time to zero
        this->time = 0;
        connectedComponents = 0;
        
        // Loop through the 2D grid
        for (int y = 0; y < HEIGHT; y++)
        {
            for (int x = 0; x < WIDTH; x++)
            {
                // Set all vertexes to undiscovered
                graph2D[x][y].discovered = false;
                
                
            }
        }
        
        for (int y = 0; y < HEIGHT; y++)
        {
            for (int x = 0; x < WIDTH; x++)
            {
               // If any vertex in the graph has not been discovered by the call to DFS_VISIT yet
                if (this->graph2D[x][y].discovered == false)
                {
                    // Call DFS_VISIT and increment the number of connected components
                    
                    this->graph2D = DFS_VISIT(x, y, graph2D);
                    connectedComponents++;
                }
            }
        }
        
    }
    
    // Recursive part of the depth-first-search algorithm
    AList** DFS_VISIT(int x, int y, AList **graph)
    {
        // Increment time and add it to the vertex's discovery time
        this->time = this->time + 1;
        graph[x][y].discoveryTime = this->time;
        // Set the vertex to discovered
        graph[x][y].discovered = true;
        
        // Create AListNode object to traverse the linked list
        AListNode *traverser = graph2D[x][y].head;
        
        while (traverser != NULL)
        {
            
            // if the graph at the index of the traverser's coordinates has not been discovered
            if (this->graph2D[traverser->x][traverser->y].discovered == false)
            {
                    
                // recursively call DFS_VISIT
                memcpy(DFS_VISIT(traverser->x, traverser->y, graph), graph, sizeof(AList) * WIDTH * HEIGHT);
            }
            
            // Move the traverser along the list
            traverser = traverser->next;
        }
        // Increment time and set it as the vertex's finish time
        this->time = this->time + 1;
        graph[x][y].finishTime = this->time;
        
        return graph;
        
    }
    
    
    double ClusteringCoefficientVertex(int x, int y)
    {
        // Variables to keep track of the number of edges and vertices in the neighborhood of the vertex
        int vertices = 0;
        int edges = 0;
        
        // Create AListNode object to traverse the linked list of edges
        AListNode *traverser = graph2D[x][y].head;
        
        // While loop the initial vertex's edges
        while (traverser != NULL)
        {
            // If the vertex in the graph has not yet been visited
            if (graph2D[traverser->x][traverser->y].visited == false)
            {
                // Increment the number of vertices and set visited to true
                vertices++;
                graph2D[traverser->x][traverser->y].visited = true;
                
                // Create a second traverser to look at all of the edges adjacent to the vertex which is adjacent to the source
                AListNode *traverser2 = graph2D[traverser->x][traverser->y].head;
            
                // While loop to scan the "edges of the edges"
                while (traverser2 != NULL)
                {
                
                    // Create a third traverser to look at all of the edges adjacent to the source again
                    AListNode *traverser3 = graph2D[x][y].head;
                
                    while (traverser3 != NULL)
                    {
                    
                        // If the vertex adjacent to the source is also adjacent to another vertex which is adjacent to the source
                        if (traverser3->x == traverser2->x && traverser3->y == traverser2->y)
                        {
                            // If the x's and y's are not the same as the source vertex (so we don't count the source)
                            if (traverser2->x != x && traverser2->y != y)
                            {
                                // If the vertex in the graph has not been visited
                                if (graph2D[traverser2->x][traverser2->y].visited == false)
                                {
                                    // Set the visited boolean to true and increment the number of edges and vertices
                                    graph2D[traverser2->x][traverser2->y].visited = true;
                                    vertices++;
                                    edges++;
                                }
                            }
                        }
                    
                        // Move the third traverser along the linked list
                        traverser3 = traverser3->next;
                    }
                
                    // Move the second traverser along the linked list
                    traverser2 = traverser2->next;
                }
            }
            
            // Move the traverser along the list
            traverser = traverser->next;
        }
        
        
        
        // Reset all visited booleans to false for the next time this function is called
        for (int y = 0; y < HEIGHT; y++)
        {
            for (int x = 0; x < WIDTH; x++)
            {
                graph2D[x][y].visited = false;
            }
        }
        
        // Calculate the Clustering Coefficient for this vertex and return it
        double toReturn = 0.0;
        if (vertices != 0 && edges != 0)
        {
            toReturn = 2 * edges;
            double divisor = vertices * (vertices - 1);
            toReturn = toReturn / divisor;
        }
            
        
        return toReturn;
        
    }
    
    
    double MeanClusteringCoefficient()
    {
        // Sum all of the clustering coefficients up from all of the vertices
        double sum = 0.0;
        for (int y = 0; y < HEIGHT; y++)
        {
            for (int x = 0; x < WIDTH; x++)
            {
                sum = sum + ClusteringCoefficientVertex(x, y);
            }
        }
        // Calculate the mean clustering coefficient by dividing by the number of vertices and return that
        double meanCC = sum / (HEIGHT * WIDTH);
        return meanCC;
    }
    
    
    double CharacteristicPathLength()
    {
        
        // Allocate memory for an array of distances
        // The size of the array will be (# of nodes) * (# of nodes)
        // I couldn't think of another way to do floyd warshall's algorithm with a matrix of distances unless I assigned each node in my graph a number (cellNumber) and then made a matrix whose length and width were the amount of nodes in the graph. This way, I can set up a distance matrix where each x and y represent a cell in the graph by its cellNumber (although it was a little tricky)
        int **distArray;
        
        distArray = new int*[V];
        
        for (int y = 0; y < (V); y++)
        {
            distArray[y] = new int[V];
        }
        
        // Initalize all values in the distance matrix to a large number
        for (int y = 0; y < HEIGHT; y++)
        {
            for (int x = 0; x < WIDTH; x++)
            {
                distArray[x][y] = 10000000;
            }
        }
        
        // 2 outer for loops to help fill out the distance matrix
        for (int y = 0; y < V; y++)
        {
            for (int x = 0; x < V; x++)
            {
                // If x and y are the same (same node)
                if (x == y)
                {
                    // Set the distance to 0
                    distArray[x][y] = 0;
                }
                // Otherwise,
                else if (distArray[x][y] != 1)
                {
                    // Loop through the graph
                    for (int y2 = 0; y2 < HEIGHT; y2++)
                    {
                        for (int x2 = 0; x2 < WIDTH; x2++)
                        {
                            // If the cellNumber is equal to y (it could be x or y; graph is undirected)
                            if (graph2D[x2][y2].cellNumber == y)
                            {
                                // Traverse the edges to see if it is connected to x
                                AListNode *traverser = graph2D[x2][y2].head;
                                
                                while (traverser != NULL)
                                {
                                    // If the cell number of the vertex that it is connected to is equal to our other cellNumber (x)
                                    if (traverser->cellNumber == x)
                                    {
                                        // Set the distance between them to 1
                                        distArray[x][y] = 1;
                                        distArray[y][x] = 1;
                                    }
                                    
                                    // Move the traverser along the linked list
                                    traverser = traverser->next;
                                }
                            }
                            
                        }
                    }
                }
            }
        }
        
        // Now we have our distance array initially set up to perform the algorithm
        
        // Use floyd warshall to calculate the minimum path length between all pairs of vertices
        for (int k = 0; k < V; k++)
        {
            for (int i = 0; i < V; i++)
            {
                for (int j = 0; j < V; j++)
                {
                    if (distArray[i][k] + distArray[k][j] < distArray[i][j])
                    {
                        distArray[i][j] = distArray[i][k] + distArray[k][j];
                    }
                }
            }
        }
        
        // integer used to sum al of the path lengths
        int sum = 0;
        // Traverse the matrix of distances
        for (int y = 0; y < V; y++)
        {
            for (int x = 0; x < V; x++)
            {
                if (distArray[x][y] < 10000000)
                {
                    // sum up all of the distances that were calculated
                    sum = sum + distArray[x][y];
                }
            }
        }
        
        // Calculate the divisor (n choose 2)
        double divisor = ((HEIGHT * WIDTH) * (HEIGHT * WIDTH - 1) / 2);
        
        // Calculate the characteristic path length by dividing the sum by the divisor
        double CPL = sum / divisor;
        
        for (int i = 0; i < V; i++)
        {
            delete[] distArray[i];
        }
        delete[] distArray;
        
        // Return characteristic path length
        return CPL;
        
        
    }
    
};




int main()
{
    cout.flush();
    float ***array3D;
    
    // Allocate memory for 3D array
    array3D = new float**[WIDTH];
    for (int i = 0; i < WIDTH; ++i)
    {
        array3D[i] = new float*[HEIGHT];
        
        for (int j = 0; j < HEIGHT; ++j)
        {
            array3D[i][j] = new float[DEPTH];
        }
    }
    
    
    // depthTracker keeps track of current "depth" i.e. grid that we are filling out
    int depthTracker = 0;
    
    // Loop to traverse years
    for (int i = 1990; i <= 2005; i++)
    {
        
        // Convert the year int to a string
        string yearString = to_string(i);
        
        // Loop to traverse the weeks
        for (int j = 1; j <= 52; j++)
        {
            // Create string to hold week value
            string weekString;
            // If we are still in the first 10 weeks
            if (j < 10)
            {
                // Add a '0' in front of the number of the week for the file formatting
                weekString = "0" + to_string(j);
            }
            else
            {
                weekString = to_string(j);
            }
            
            // Generate the file name for each grid
            string fileName = "Beaufort_Sea_diffw" + weekString + "y" + yearString + "+landmask";
            
            // Create and open a filestream object to read binary data
            fstream inputFile;
            inputFile.open("CS310_project_subregion/" + yearString + "/" + fileName, ios::binary | ios::in);
            
            // Loop to handle height of 2D array
            for (int k = 0; k < HEIGHT; k++)
            {
                // Loop to handle width of 2D array
                for (int l = 0; l < WIDTH; l++)
                {
                    // Float variable to initially read in each item
                    float aValue;
                    
                    
                    
                    // Read each item and store it in aValue
                    inputFile.read((char *) &aValue, sizeof(float));
                    // Place the value we read and stored into aValue into the array at it's correct location
                    array3D[l][k][depthTracker] = aValue;
                    
                    
                    
                }
            }
            
            
            // Close each file after we finish with it
            inputFile.close();
            // Increment depthTracker every time we are about to generate a new fileName
            depthTracker++;
        }
        
    }
    
    
    // Integers to hold the number of edges in each graph
    int edges95 = 0;
    int edges925 = 0;
    int edges9 = 0;
    
    // Create new graphs for each threshold
    Graph graph95(HEIGHT * WIDTH, array3D);
    Graph graph925(HEIGHT * WIDTH, array3D);
    Graph graph9(HEIGHT * WIDTH, array3D);
    
    /************* ONE WAY OF ADDING THE EDGES ***************/
    
    
    /*// The 2 outer for loops are for x and y
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
                // The 2 inner for loops are for xPrime and yPrime
                // Set xPrime = x and yPrime = y so we don't make comparisons we already have
                for (int yPrime = y; yPrime < HEIGHT; yPrime++)
                {
                    for (int xPrime = x; xPrime < WIDTH; xPrime++)
                    {
                        // If both traversers are looking at the same cell
                        if (x == xPrime && y == yPrime)
                        {
                            // Do nothing
                        }
                        // In all other cases
                        else
                        {
                                // Calculate the Pearson Correlation Coefficient for the two points
                                float PCC = r(array3D, x, y, xPrime, yPrime);
                                
                                //cout << "PCC between " << x << ", " << y << " and " << xPrime << ", " << yPrime << " is " << PCC << endl;
                                
                                // If |PCC| is greater than 0.95
                                if (abs(PCC) >= 0.95)
                                {
                                    //cout << "Edge between (" << x << ", " << y << ")" << " & " << "(" << xPrime << ", " << yPrime << ")" << endl;
                                    // Update the proper edges count and add an edge to the proper graph
                                    edges95++;
                                    graph95.addEdge(x, y, xPrime, yPrime);
                                }
                                
                                // If |PCC| is greater than 0.925
                                if (abs(PCC) >= 0.925)
                                {
                                    // Update the proper edges count and add an edge to the proper graph
                                    edges925++;
                                    graph925.addEdge(x, y, xPrime, yPrime);
                                }
                                
                                // If |PCC| is greater than 0.9
                                if (abs(PCC) >= 0.9)
                                {
                                    // Update the proper edges count and add an edge to the proper graph
                                    edges9++;
                                    graph9.addEdge(x, y, xPrime, yPrime);
                                }
                            
                        }
                    }
                }
            
        }
        // Helps keep track of where algorithm is when calculating the PCCs
        cout << "Done with row " << y + 1 << endl;
    }*/
    
    
    
    /****************** SECOND WAY OF ADDING THE EDGES************/
    
    // Create a grid of GridNode structures dynamically
    GridNode **grid;
    // Allocate space for our 2D graph of AList structures
    grid = new GridNode*[WIDTH];
    for (int i = 0; i < WIDTH; i++)
    {
        grid[i] = new GridNode[HEIGHT];
    }
    
    // Calculate the summation needed for the Pearson Coefficient for each vertex ( so we only have to do this once)
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            grid[x][y].summation = Sxx(array3D, x, y);
        }
    }
    
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            
            for (int yPrime = y; yPrime < HEIGHT; yPrime++)
            {
                for (int xPrime = x; xPrime < WIDTH; xPrime++)
                {
                    // If x's and y's are the same
                    if (x == xPrime && y == yPrime)
                    {
                       // Do nothing
                    }
                    else
                    {
                        // In all other cases, calculate the coefficient using the calculations we have already done
                        float r = Sxy(array3D, x, y, xPrime, yPrime) / sqrt(grid[x][y].summation * grid[xPrime][yPrime].summation);
                        
                        // Add edges where necessary
                        if (abs(r) >= 0.95)
                        {
                            edges95++;
                            graph95.addEdge(x, y, xPrime, yPrime);
                        }
                        
                        if (abs(r) >= 0.925)
                        {
                            // Update the proper edges count and add an edge to the proper graph
                            edges925++;
                            graph925.addEdge(x, y, xPrime, yPrime);
                        }
                        
                        if (abs(r) >= 0.9)
                        {
                            // Update the proper edges count and add an edge to the proper graph
                            edges9++;
                            graph9.addEdge(x, y, xPrime, yPrime);
                        }
                    }
                }
            }
        }
        //cout << "Done with row: " << y + 1 << endl;
    }
    
    // Delete the grid after we make the graph
    
    for (int i = 0; i < WIDTH; i ++)
    {
        delete[] grid[i];
    }
    delete[] grid;
    
    /******************************************************/
    cout << "\nNumber of edges:\n";
    
    cout << ".95 edges: " << edges95 << endl;
    cout << ".925 edges: " << edges925 << endl;
    cout << ".9 edges: " << edges9 << endl;
    
    
    /************INVESTIGATIONS**********/
    
    
    
    /*// Variables used to keep track of the counts of each degree in our graph
    int count0 = 0;
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int count4 = 0;
    int count5 = 0;
    int count6 = 0;
    int count7 = 0;
    int count8 = 0;
    int count9 = 0;
    int count10 = 0;
    int count11 = 0;
    int count12 = 0;
    int count13 = 0;
    int count14 = 0;
    int count15 = 0;
    int count16 = 0;
    int count17 = 0;
    int count18 = 0;
    int count19 = 0;
    int count20 = 0;
    int count21 = 0;
    int count22 = 0;
    int count23 = 0;
    int count24 = 0;
    int count25 = 0;
    int count26 = 0;
    int count27 = 0;
    int count28 = 0;
    int count29 = 0;
    int count30 = 0;
    int count31 = 0;
    int count32 = 0;
    int count33 = 0;
    int count34 = 0;
    int count35 = 0;
    int count36 = 0;
    int count37 = 0;
    int count38 = 0;
    int count39 = 0;
    int count40 = 0;
    int count41 = 0;
    int count42 = 0;
    
    
    // For loops to traverse graph
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            // Call ComputeDegree to get the degree for the cell (can change graph name to get degrees for another threshold graph)
            int aDegree = graph9.ComputeDegree(x, y);
            
            // Update the proper degree count variable
            if (aDegree == 0)
            {
                count0++;
            }
            else if (aDegree == 1)
            {
                count1++;
            }
            else if (aDegree == 2)
            {
                count2++;
            }
            else if (aDegree == 3)
            {
                count3++;
            }
            else if (aDegree == 4)
            {
                count4++;
            }
            else if (aDegree == 5)
            {
                count5++;
            }
            else if (aDegree == 6)
            {
                count6++;
            }
            else if (aDegree == 7)
            {
                count7++;
            }
            else if (aDegree == 8)
            {
                count8++;
            }
            else if (aDegree == 9)
            {
                count9++;
            }
            else if (aDegree == 10)
            {
                count10++;
            }
            else if (aDegree == 11)
            {
                count11++;
            }
            else if (aDegree == 12)
            {
                count12++;
            }
            else if (aDegree == 13)
            {
                count13++;
            }
            else if (aDegree == 14)
            {
                count14++;
            }
            else if (aDegree == 15)
            {
                count15++;
            }
            else if (aDegree == 16)
            {
                count16++;
            }
            else if (aDegree ==17)
            {
                count17++;
            }
            else if (aDegree == 18)
            {
                count18++;
            }
            else if (aDegree == 19)
            {
                count19++;
            }
            else if (aDegree == 20)
            {
                count20++;
            }
            else if (aDegree == 21)
            {
                count21++;
            }
            else if (aDegree == 22)
            {
                count22++;
            }
            else if (aDegree == 23)
            {
                count23++;
            }
            else if (aDegree == 24)
            {
                count24++;
            }
            else if (aDegree == 25)
            {
                count25++;
            }
            else if (aDegree == 26)
            {
                count26++;
            }
            else if (aDegree == 27)
            {
                count27++;
            }
            else if (aDegree == 28)
            {
                count28++;
            }
            else if (aDegree == 29)
            {
                count29++;
            }
            else if (aDegree == 30)
            {
                count30++;
            }
            else if (aDegree == 31)
            {
                count31++;
            }
            else if (aDegree == 32)
            {
                count32++;
            }
            else if (aDegree == 33)
            {
                count33++;
            }
            else if (aDegree == 34)
            {
                count34++;
            }
            else if (aDegree == 35)
            {
                count35++;
            }
            else if (aDegree == 36)
            {
                count36++;
            }
            else if (aDegree == 37)
            {
                count37++;
            }
            else if (aDegree == 38)
            {
                count38++;
            }
            else if (aDegree == 39)
            {
                count39++;
            }
            else if (aDegree == 40)
            {
                count40++;
            }
            else if (aDegree == 41)
            {
                count41++;
            }
            else if (aDegree == 42)
            {
                count42++;
            }
            else if (aDegree > 20)
            {
                
                cout << aDegree << endl;
                
            }
            
        }
    }
    
    cout << "\nDegree Distributions:\n";
    cout << "0: " << count0 << endl;
    cout << "1: " << count1 << endl;
    cout << "2: " << count2 << endl;
    cout << "3: " << count3 << endl;
    cout << "4: " << count4 << endl;
    cout << "5: " << count5 << endl;
    cout << "6: " << count6 << endl;
    cout << "7: " << count7 << endl;
    cout << "8: " << count8 << endl;
    cout << "9: " << count9 << endl;
    cout << "10: " << count10 << endl;
    cout << "11: " << count11 << endl;
    cout << "12: " << count12 << endl;
    cout << "13: " << count13 << endl;
    cout << "14: " << count14 << endl;
    cout << "15: " << count15 << endl;
    cout << "16: " << count16 << endl;
    cout << "17: " << count17 << endl;
    cout << "18: " << count18 << endl;
    cout << "19: " << count19 << endl;
    cout << "20: " << count20 << endl;
    cout << "21: " << count21 << endl;
    cout << "22: " << count22 << endl;
    cout << "23: " << count23 << endl;
    cout << "24: " << count24 << endl;
    cout << "25: " << count25 << endl;
    cout << "26: " << count26 << endl;
    cout << "27: " << count27 << endl;
    cout << "28: " << count28 << endl;
    cout << "29: " << count29 << endl;
    cout << "30: " << count30 << endl;
    cout << "31: " << count31 << endl;
    cout << "32: " << count32 << endl;
    cout << "33: " << count33 << endl;
    cout << "34: " << count34 << endl;
    cout << "35: " << count35 << endl;
    cout << "36: " << count36 << endl;
    cout << "37: " << count37 << endl;
    cout << "38: " << count38 << endl;
    cout << "39: " << count39 << endl;
    cout << "40: " << count40 << endl;
    cout << "41: " << count41 << endl;
    cout << "42: " << count42 << endl;
    
    
    graph95.DFS();
    graph925.DFS();
    graph9.DFS();
    
    cout << "\nConnected Components:\n";
    cout << ".95: " << graph95.connectedComponents << endl;
    cout << ".925: " << graph925.connectedComponents<< endl;
    cout << ".9: " << graph9.connectedComponents << endl;*/
    
    
    cout << "Mean Clustering Coefficient of Graph threshold .95 " << graph95.MeanClusteringCoefficient() << endl;
    cout << "Mean Clustering Coefficient of Graph threshold.925 " << graph925.MeanClusteringCoefficient() << endl;
    cout << "Mean Clustering Coefficient of Graph threshold .9 " << graph9.MeanClusteringCoefficient() << endl;
    
    
    double meanDegree95 = 0.0;
    double meanDegree925 = 0.0;
    double meanDegree9 = 0.0;
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            meanDegree95 = meanDegree95 + graph95.ComputeDegree(x, y);
            meanDegree925 = meanDegree925 + graph925.ComputeDegree(x, y);
            meanDegree9 = meanDegree9 + graph9.ComputeDegree(x, y);
        }
    }
    
    meanDegree95 = meanDegree95 / (WIDTH * HEIGHT);
    meanDegree925 = meanDegree925 / (WIDTH * HEIGHT);
    meanDegree9 = meanDegree9 / (WIDTH * HEIGHT);
    
    cout << "Mean Clustering Coefficient of Random Graph threshold .95: " << (meanDegree95 / (WIDTH * HEIGHT)) << endl;
    cout << "Mean Clustering Coefficient of Random Graph threshold .925: " << (meanDegree925 / (WIDTH * HEIGHT)) << endl;
    cout << "Mean Clustering Coefficient of Random Graph threshold .9: " << (meanDegree9 / (WIDTH * HEIGHT)) << endl;
    
    double CPL95 = graph95.CharacteristicPathLength();
    graph95.~Graph();
    double CPL925 = graph925.CharacteristicPathLength();
    graph925.~Graph();
    double CPL9 = graph9.CharacteristicPathLength();
    graph9.~Graph();
    
    double CPL95Random = (log(HEIGHT * WIDTH) / log(meanDegree95));
    double CPL925Random = (log(HEIGHT * WIDTH) / log(meanDegree925));
    double CPL9Random = (log(HEIGHT * WIDTH) / log(meanDegree9));
    
    cout << "CPL of Graph threshold .95: " << CPL95 << endl;
    cout << "CPL of Graph threshold .925: " << CPL925 << endl;
    cout << "CPL of Graph threshold .9: " << CPL9 << endl;
    
    cout << "CPL of random Graph threshold .95: " << CPL95Random << endl;
    cout << "CPL of random Graph threshold .925: " << CPL925Random << endl;
    cout << "CPL of random Graph threshold .9: " << CPL9Random << endl;
    
    
    
    // Delete all arrays before finishing program
    for (int i = 0; i < HEIGHT; ++i)
    {
        for (int j = 0; j < WIDTH; ++j)
        {
            delete [] array3D[i][j];
        }
        
        delete [] array3D[i];
    }
    delete [] array3D;
    
    
}
