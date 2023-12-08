#include <iostream>
#include <stdlib.h>
#include <utility>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <math.h>

std::vector<std::pair<int, int>> ChooseInput(int argc, char* argv[]);
std::vector<std::pair<int, int>> GenerateRandomGraph(int size);
std::vector<std::pair<int, int>> GetGraphFromFile(std::string file);
std::vector<std::vector<double>> Matrix(std::vector<std::pair<int, int>> coords, int size);
std::pair<double, std::vector<int>> greedyAlgorithmTSP(int size, std::vector<std::vector<double>> matrix);
std::pair<double, std::vector<int>> AntsAlgorithm(int ants, int iterations, float alpha, float beta, float p, float Q, float c, std::vector<std::vector<double>> &Matrix);
std::pair<std::vector<float>, std::vector<int>> DataForAnts();
void OutputToTerminal(double allDistance, std::vector<int> path, int size);
void SaveOutput(double allDistance, std::vector<int> path);
void probabilityValue(int current, float alpha, float beta, std::vector<std::vector<double>> &Matrix, std::vector<long double> &probabilityAnts, std::vector<std::vector<float>> &trailIntensity, std::vector<int> &allowAnt);
int next(std::vector<long double> &probability, std::vector<std::vector<double>> &Matrix, std::vector<int> &allowAnt);
int ChooseAnAlgorithm(int choice);
double GetDistance(std::vector<std::pair<int, int>> coords, int i, int j);
long double randomLongDouble();

int main(int argc, char *argv[]){
    std::vector<std::pair<int, int>> coords;
    std::vector<std::vector<double>> matrix;
    std::pair<double, std::vector<int>> output;
    std::pair<std::vector<float>, std::vector<int>> dataAnts;
    std::vector<float> FloatDataForAnts;
    std::vector<int> IntDataForAnts;
    std::vector<int> path;
    double allDistance;
    int choice;
    
    
    coords = ChooseInput(argc, argv);
    int size = coords.size();
    matrix = Matrix(coords, size);

    choice = ChooseAnAlgorithm(choice);
    dataAnts = DataForAnts();
    FloatDataForAnts = dataAnts.first;
    IntDataForAnts = dataAnts.second;

    if(choice){
        output = AntsAlgorithm(IntDataForAnts[0], IntDataForAnts[1], FloatDataForAnts[0], FloatDataForAnts[1], FloatDataForAnts[2], FloatDataForAnts[3], FloatDataForAnts[4], matrix);
    } else {
        output = greedyAlgorithmTSP(size, matrix);
    }
    SaveOutput(output.first, output.second);
    OutputToTerminal(output.first, output.second, size);
    return 0;
}



int ChooseAnAlgorithm(int choice){
    std::cout<<"Choose a method:"<<std::endl;
    std::cout<<"(0) Greedy algorithm"<<std::endl;
    std::cout<<"(1) Ant-colony algorithm"<<std::endl;
    std::cin>>choice;
    return choice;
}

std::pair<std::vector<float>, std::vector<int>> DataForAnts() {
    std::vector<float> FloatData(5, 0);
    std::vector<int> IntData(2,0);
    float alpha, beta, p, Q, c;
    int ants, iterations;

    std::cout<<"Enter Alpha: "<<std::endl;
    std::cin>>alpha;
    FloatData.push_back(alpha);
    std::cout<<"Enter Beta:"<<std::endl;
    std::cin>>beta;
    FloatData.push_back(beta);
    std::cout<<"Enter p: "<<std::endl;
    std::cin>>p;
    FloatData.push_back(p);
    std::cout<<"Enter Q: "<<std::endl;
    std::cin>>Q;
    FloatData.push_back(Q);
    std::cout<<"Enter c: "<<std::endl;
    std::cin>>c;
    FloatData.push_back(c);

    std::cout<<"Enter number of ants: "<<std::endl;
    std::cin>>ants;
    IntData.push_back(ants);
    std::cout<<"Enter number of iterations: "<<std::endl;
    std::cin>>iterations;
    IntData.push_back(iterations);

    return std::pair<std::vector<float>, std::vector<int>> (FloatData, IntData);
}

std::vector<std::pair<int, int>> ChooseInput(int argc, char* argv[]){
    std::string file;
    std::vector<std::pair<int, int>> coords;
    int size, choice;

    if (argc==2){
        file = argv[1];
        coords = GetGraphFromFile(file);
    } else if (argc == 1) {
        std::cout<<"Input graph size: "<<std::endl;
        std::cin>>size;
        coords = GenerateRandomGraph(size);
    }

    return coords;
}

std::vector<std::pair<int, int>> GenerateRandomGraph(int size) {
    std::ofstream outputFile("random.txt");
    std::vector<std::pair<int, int>> coords;
    std::pair<int, int> newPair;
    int coord1, coord2;

    srand(time(nullptr));
    

    for (int i = 0; i < size; i++) { // setting coordinates for points
        coord1 = rand() % 10000;
        coord2 = rand() % 10000;
        newPair.first = coord1;
        newPair.second = coord2;
        std::cout << newPair.first << " ";
        std::cout << newPair.second << "\n";
        coords.push_back(newPair);
    }
    if (outputFile.is_open()) {
        outputFile << size << std::endl;
        for(int i = 0; i<size; i++){
            outputFile << i+1 << " " << coords[i].first <<" "<< coords[i].second << std::endl;
        }
    } else {
        std::cout<<"Opening file for a random graph input failed"<<std::endl;
    }
    return coords;
}
std::vector<std::pair<int, int>> GetGraphFromFile(std::string file) {
    std::ifstream inputFile(file);
    std::vector<std::pair<int, int>> coords;
    std::string coord1, coord2;
    int size, nr; 
    if (inputFile.is_open()) {
        inputFile >> size;
        while(inputFile >> nr >> coord1 >> coord2){
            coords.push_back(std::pair<int,int>(std::stoi(coord1), std::stoi(coord2)));
        }
        return coords;
    } else {
        std::cout<<"Opening file with coordinates failed"<<std::endl;
    }
    return coords; 

}
double GetDistance(std::vector<std::pair<int, int>> coords, int i, int j){
    return sqrt((coords[j].first - coords[i].first)*(coords[j].first - coords[i].first) + (coords[j].second - coords[i].second)*(coords[j].second - coords[i].second));
}
std::vector<std::vector<double>> Matrix(std::vector<std::pair<int, int>> coords, int size){
    std::vector<std::vector<double>> matrix;
    for ( int i = 0; i < size; i++) {
        matrix.push_back(std::vector<double>(size, -1));
        for (int j = 0; j < size; j++){
            if (i != j){
                matrix[i][j] = GetDistance(coords, i, j);
            }
        }
    }
    return matrix;
}
std::pair<double, std::vector<int>> greedyAlgorithmTSP(int size, std::vector<std::vector<double>> matrix){
    std::vector<int> path;
    std::vector<int> visited(size, 0);
    double allDistance = 0.0;
    int currentVertex, destination;
    double currentValue;

    visited[0] = 1;
    path.push_back(0);
    
    while(path.size()<size){
        currentValue = 100000000000000000;
        currentVertex = path.back();
        for (int i = 0; i < size; i++){
            if(currentVertex != i && !visited[i] && matrix[currentVertex][i] < currentValue) {
                currentValue = matrix[currentVertex][i];
                destination = i;
            }
        }
        path.push_back(destination);
        visited[destination] = 1;
        allDistance += currentValue;
    }
    path.push_back(0);
    allDistance += matrix[0][destination];

    return std::pair<double, std::vector<int>>(allDistance, path);
}
void SaveOutput(double allDistance, std::vector<int> path) {
    std::ofstream resultFile("output.txt");
    if (resultFile.is_open()) {
        resultFile << "Distance: "<< allDistance <<std::endl;
        for (int i = 0; i < path.size(); i++){
            resultFile << path[i] + 1 << "\n";
        }
        resultFile << std::endl;
    } else {
        std::cout<<"Could not open result file"<<std::endl;
    }
}
void OutputToTerminal(double allDistance, std::vector<int> path, int size) {
    std::cout << "Length of path: " << allDistance << std::endl << "Path: ";
    for ( int i = 0; i <= size; i++ ){
        std::cout << path[i] + 1 << " ";
    }
}

//Ants
long double randomLongDouble() {
    std::srand(time(nullptr));
    return (long double)(rand()) / (long double) (RAND_MAX);
}

int next(std::vector<long double> &probability, std::vector<std::vector<double>> &Matrix, std::vector<int> &allowAnt){
    long double choice;
    int size = Matrix.size();
    choice = randomLongDouble();

    for (int i = 0; i < size; i++) {
        if (choice < probability[i]) {
            return i;
        }
        choice -= probability[i];
    }
    for (int j = 0; j < size; j++) {
        if(allowAnt[j] == 1) {
            return j;
        }
    }
    return -1;
}
void probabilityValue(int current, float alpha, float beta, std::vector<std::vector<double>> &Matrix, std::vector<long double> &probabilityAnts, std::vector<std::vector<float>> &trailIntensity, std::vector<int> &allowAnt){
    
    int size = Matrix.size();
    long double numerator, value, valueEnd, denominator = 0.0;

    for (int i = 0; i < size; i++){
        if(allowAnt[i] == 0) {
            continue;
        }
        if(i = current) {
            allowAnt[i] = 0;
            continue;
        }
        value = pow(trailIntensity[current][i], alpha) * pow(1.0/Matrix[current][i], beta);
        if(value>=0){
            denominator+= value;
        }
    }
    for (int i = 0; i < size; i++){
        if(allowAnt[i] == 0){
            probabilityAnts[i] = 0;
        } else {
            if (denominator <= 0) {
                probabilityAnts[i] = 0;
                continue;
            }
            numerator = pow(trailIntensity[current][i], alpha)/pow(Matrix[current][i], beta);
            valueEnd = numerator/denominator;
            probabilityAnts[i] = numerator/denominator;
        }
    }
}
std::pair<double, std::vector<int>> AntsAlgorithm(int ants, int iterations, float alpha, float beta, float p, float Q, float c, std::vector<std::vector<double>> &Matrix){
    int size = Matrix.size();
    float dist = 0.0;
    std::vector<std::vector<float>> trailIntensity (size, std::vector<float>(size, c));
    std::vector<std::vector<float>> trailIntensityTmp (size, std::vector<float>(size, 0));
    std::vector<std::vector<int>> antsTrail (ants, std::vector<int>());
    std::vector<std::vector<int>> allowAnt (ants, std::vector<int>(size, 1));
    std::vector<std::vector<long double>> probabilityAnts (ants, std::vector<long double> (size, 0));
    std::vector<int> path;
    double allDist = 100000000000000000;


    for (int a = 0; a < ants; a++){
        antsTrail[a].push_back(rand() % (size - 0) + 0);
        allowAnt[a][antsTrail[a][0]] = 0;
    }
    for(int k = 0; k < size -1; k++) {
        for (int a = 0; a < ants; a++){
            int current = antsTrail[a].back();
            probabilityValue(current, alpha, beta, Matrix, probabilityAnts[a], trailIntensity, allowAnt[a]);
            int point = next(probabilityAnts[a], Matrix, allowAnt[a]);
            antsTrail[a].push_back(point);
            allowAnt[a][point] = 0;
        }
    }
    for (int a = 0; a < ants; a++) {
        dist = 0.0;
        antsTrail[a].push_back(antsTrail[a][0]);
        for (int k = 0; k < size; k++) {
            dist+=Matrix[antsTrail[a][k]][antsTrail[a][k+1]];
        }
        for(int k = 0; k < size; k++) {
            trailIntensityTmp[antsTrail[a][k+1]][antsTrail[a][k]] += Q/dist;
            trailIntensityTmp[antsTrail[a][k]][antsTrail[a][k+1]] += Q/dist;
        }
        if (dist < allDist){
            allDist = dist;
            path = antsTrail[a];
        }
        for (int k = 0; k < size; k++) {
            for (int l = k+1; l < size; l++) {
                trailIntensity[k][l] = (p*trailIntensity[k][l]) + trailIntensityTmp[k][l];
                trailIntensity[l][k] = trailIntensity[k][l];
                trailIntensityTmp[k][l] = 0.0;
                trailIntensityTmp[l][k] = 0.0;
            }
        }
        for (int a = 0; a < ants; a++) {
            antsTrail[a].clear();
            allowAnt[a] = std::vector<int> (size, 1);
        }
        
    }
    return std::pair<double, std::vector<int>>(allDist, path);
}

