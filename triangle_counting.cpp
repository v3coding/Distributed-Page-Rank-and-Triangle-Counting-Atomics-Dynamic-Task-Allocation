#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

static int numberOfThreads;

//array1 = the set of in neighbors for vertex in question u
//len1 = the number of in neighbors for vertex u
//array2 = the set of out neighbors of the neighbor vertex v which is a neighbor of u
//len 2 = the number of out neighbors of neighbor v which is a neighbor of u
//u = the vertex in question, u
//v = a neighbor of u, the neighbor being examined
long countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2,
                    uintV u, uintV v) {

  uintE i = 0, j = 0; // indexes for array1 and array2
  long count = 0;

  //if u == v, meaning, if the neighbor of a vertex v is actually itself
  if (u == v)
    return count;

  //iterate while condition 1, iterator for (u's set of neighbors, which is i, is less than the total neighbors) (i < len1)
  //iterate while condition 2, iterator for (the vertex in question, u, has a neighbor, v, and we are iterating for v's out neighbors, for the number of out neighbors) (j < len2)
  //check for triangle
  //we need to split this in to edges
  while ((i < len1) && (j < len2)) {
    //If array[i] (an in neighbor for u) is equal to (array2[j] (an out neighbor for v) then we have a triangle based on the fact that v is an out neighbor for u
    if (array1[i] == array2[j]) {
      //if array1[i] != u (the in neighbor of u is not u itself) and (array1[i] != v) the in neighbor for u is not its neighbor v
      if ((array1[i] != u) && (array1[i] != v)) {
        count++;
      } else {
        // triangle with self-referential edge -> ignore
      }
      i++;
      j++;
    } else if (array1[i] < array2[j]) {
      i++;
    } else {
      j++;
    }
  }
  return count;
}

void triangleCountParallel(uintV start, uintV finish, Graph &g, uintV n, uintV& triangle_count,double& time_taken){
  timer t2;
  t2.start();
    for (uintV u = start; u < finish; u++) {
    // For each outNeighbor v, find the intersection of inNeighbor(u) and
    // outNeighbor(v)
    uintE out_degree = g.vertices_[u].getOutDegree();
    for (uintE i = 0; i < out_degree; i++) {
      uintV v = g.vertices_[u].getOutNeighbor(i);
      triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                       g.vertices_[u].getInDegree(),
                                       g.vertices_[v].getOutNeighbors(),
                                       g.vertices_[v].getOutDegree(), u, v);
    }
  }
  time_taken = t2.stop();
}

void triangleCountVertexDriver(Graph &g) {
  uintV n = g.n_;
  long triangle_count = 0;
  double time_taken = 0.0;
  timer t1;
  uintV Counter[numberOfThreads];
  for(int i = 0; i < numberOfThreads; i++){
    Counter[i] = 0;
  }
  // The outNghs and inNghs for a given vertex are already sorted

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  double numThreads = (double) numberOfThreads;
  double numVertexes = (double) n;
  double start = 0;
  double finish = 0;
  uintV startCaster = 0;
  uintV finishCaster = 0;
  std::thread threads[numberOfThreads];
  double timers[numberOfThreads];
  for(int i = 0; i < numberOfThreads; i++){
    timers[i] = 0;
  }

  t1.start();
  // Process each edge <u,v>
    for (int i = 0; i < numberOfThreads; i++){
    if(((start + numVertexes) / numThreads) < n){
      finish = start + (numVertexes / numThreads);
    }
    else {
      finish = numVertexes;
    }
    startCaster = (uintV) start;
    finishCaster = (uintV) finish;
    threads[i] = std::thread(triangleCountParallel,startCaster,finishCaster,std::ref(g),n,std::ref(Counter[i]),std::ref(timers[i]));
    start = finish;
  }
  for(int i = 0; i < numberOfThreads; i++){
    threads[i].join();
  }

  time_taken = t1.stop();

  triangle_count = 0;
  for(int i = 0; i < numberOfThreads; i++){
    triangle_count += Counter[i];
  }
  // -------------------------------------------------------------------
  // Here, you can just print the number of non-unique triangles counted by each
  // thread std::cout << "thread_id, triangle_count, time_taken\n"; Print the
  // above statistics for each thread Example output for 2 threads: thread_id,
  // triangle_count, time_taken 1, 102, 0.12 0, 100, 0.12
  std::cout << "thread_id, triangle_count, time_taken" << std::endl;
  for(int i = numberOfThreads-1; i >= 0; i--){
    std::cout << i << ", " << Counter[i] << ", " << timers[i] << " " << std::endl; 
  }

  // Print the overall statistics
  std::cout << "Number of triangles : " << triangle_count << "\n";
  std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
}




void triangleCountEdgeDriver(Graph &g) {
  uintV n = g.n_;
  uintV edges = g.m_;
  long triangle_count = 0;
  double time_taken = 0.0;
  timer t1;

  int edgeStart[numberOfThreads];
  int edgeEnd[numberOfThreads];
  uintV vertexStart[numberOfThreads];
  uintV vertexEnd[numberOfThreads];

  for(int i = 0; i < numberOfThreads; i++){
    edgeStart[i] = 0;
    edgeEnd[i] = 0;
    vertexStart[i] = 0;
    vertexEnd[i] = 0;
  }


  double start = 0;
  double finish = 0;
  double numEdges = g.m_;
  double numThreads = numberOfThreads;


  for (int i = 0; i < numberOfThreads; i++){
    if(((start + numEdges) / numThreads) < edges){
      finish = start + (numEdges / numThreads);
    }
    else {
      finish = numEdges;
    }
    edgeStart[i] = start;
    edgeEnd[i] = finish;
    start = finish;
  }

  uintE currentEdge = 0;
  int edgeIterator = 0;
  uintV Counter[numberOfThreads];
  for(int i = 0; i < numberOfThreads; i++){
    Counter[i] = 0;
  }
  double timers[numberOfThreads];
  for(int i = 0; i < numberOfThreads; i++){
    timers[i] = 0;
  }


  edgeStart[0] = 0;
  vertexStart[0] = 0;
  std::thread threads[numberOfThreads];

  for(uintV x = 0; x < n; x++){
    currentEdge += g.vertices_[x].getOutDegree();
      if(currentEdge >= edgeEnd[edgeIterator]){
        vertexEnd[edgeIterator] = x;
        edgeIterator++;
        vertexStart[edgeIterator] = x;
      }
      if(edgeIterator == numberOfThreads-1){
        x = n;
        vertexEnd[numberOfThreads-1] = n;
      }
  }
  //std::cout << "Total number of verticies = " << g.n_ << std::endl;
  //for(int i = 0; i < numberOfThreads; i++){
  //  std::cout << "Thread " << i << " will beging iterate from vertex " << vertexStart[i] << " to vertex " << vertexEnd[i] << " covering edges from " << edgeStart[i] << " to " << edgeEnd[i] << std::endl;
  //}

  t1.start();

 // std::cout << "here" << std::endl;

  for(int i = 0; i < numberOfThreads; i++){
    threads[i] = std::thread(triangleCountParallel,vertexStart[i],vertexEnd[i],std::ref(g),n,std::ref(Counter[i]),std::ref(timers[i]));
  }

  for(int i = 0; i < numberOfThreads; i++){
    threads[i].join();
  }

  time_taken = t1.stop();

  triangle_count = 0;
  for(int i = 0; i < numberOfThreads; i++){
    triangle_count += Counter[i];
  }
  // -------------------------------------------------------------------
  // Here, you can just print the number of non-unique triangles counted by each
  // thread std::cout << "thread_id, triangle_count, time_taken\n"; Print the
  // above statistics for each thread Example output for 2 threads: thread_id,
  // triangle_count, time_taken 1, 102, 0.12 0, 100, 0.12
  std::cout << "thread_id, triangle_count, time_taken" << std::endl;
  for(int i = numberOfThreads-1; i >= 0; i--){
    std::cout << i << ", " << Counter[i] << ", " << timers[i] << " " << std::endl; 
  }

  // Print the overall statistics
  std::cout << "Number of triangles : " << triangle_count << "\n";
  std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
}

void triangleCountSerial(Graph &g) {
  uintV n = g.n_;
  long triangle_count = 0;
  double time_taken = 0.0;
  timer t1;
  t1.start();
  for (uintV u = 0; u < n; u++) {
    uintE out_degree = g.vertices_[u].getOutDegree();
    for (uintE i = 0; i < out_degree; i++) {
      uintV v = g.vertices_[u].getOutNeighbor(i);
      triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                       g.vertices_[u].getInDegree(),
                                       g.vertices_[v].getOutNeighbors(),
                                       g.vertices_[v].getOutDegree(), u, v);
    }
  }
  time_taken = t1.stop();
  std::cout << "Number of triangles : " << triangle_count << "\n";
  std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
}


int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "triangle_counting_serial",
      "Count the number of triangles using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
          {"strategy", "Strategy to be used",
           cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  uint strategy = cl_options["strategy"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";
  numberOfThreads = n_workers;
  std::cout << "Task decomposition strategy : " << strategy << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";

  switch (strategy) {
  case 0:
    std::cout << "\nSerial\n";
    triangleCountSerial(g);
    break;
  case 1:
    std::cout << "\nVertex-based work partitioning\n";
    triangleCountVertexDriver(g);
    break;
  case 2:
    std::cout << "\nEdge-based work partitioning\n";
    triangleCountEdgeDriver(g);
    break;
  case 3:
    std::cout << "\nDynamic task mapping\n";
    break;
  default:
    break;
  }

  return 0;
}
