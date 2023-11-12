#include <iostream>
#include <stdlib.h>
#include <utility>
#include <cmath> // square function
#include <algorithm> // find function
#include <vector>


int main () {
    int size;
    typedef std::pair<int,int> P;
    
    srand(time(nullptr));

    std::cout<<"Input graph size: "<<std::endl;
    std::cin >> size;

    std::vector<P> coords;
    std::vector<std::vector<float>> matrix;
     

    for ( int i = 0, x, y; i < size; i++) { // setting coordinates for points
        y = rand() % 10000;
        x = rand() % 10000;
        P coordinates(x, y);
        std::cout << coordinates.first << " ";
        std::cout << coordinates.second << "\n";
        coords.push_back(coordinates);
    }
    

    
    for ( int i = 0; i < size; i++) { // calculating distance between two points and puting it in matrix
        matrix.push_back(std::vector<float>(size, -1));
        for (int j = 0; j < size; j++){
            if (i != j){
                float distance = sqrt((coords[j].first - coords[i].first)*(coords[j].first - coords[i].first) + (coords[j].second - coords[i].second)*(coords[j].second - coords[i].second));
                matrix[i][j] = distance;
            }
        }
    }
    
    
    std::vector<int> path;
    std::vector<int> visited(size, 0);
    int lenght = 0;
    path.push_back(0);
    visited[0] = 1;
    int dest;
    float allDist = 0.0;

    while(path.size() < size) {
        int current = path.back();
        float currentValue = 100000000000000000;
        for (int j = 0; j < size; j ++) {
            if(current!=j && !visited[j] && matrix[current][j] >= 0 && matrix[current][j] <= currentValue){
                dest = j;
                currentValue = matrix[current][j];
            }
        }
        path.push_back(dest);
        visited[dest] = 1;
        allDist += currentValue;

    }
    path.push_back(0);
    allDist += matrix[0][dest];
    
    //showing results
    
    std::cout << "Lenght of path " << allDist << "\n Path:";
    
    for ( int i = 0; i < size; i++ ){
        std::cout << path[i] << " ";
    }
    return 0;
}