#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <unistd.h>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
#endif

static int numberOfThreads;
static CustomBarrier* my_barrier;
static int iterating;
std::atomic<int> currentVertex;
static uintV granularityInt;

void pageRankParallel(uintV k, uintV n, int max_iters, std::atomic<PageRankType> *pr_curr, std::atomic<PageRankType> *pr_next, Graph &g, double& time_taken)
{
  timer t1;
  time_taken = 0.0;
  PageRankType tempRank;
  std::thread::id this_id = std::this_thread::get_id();

  t1.start();
  for (int iter = 0; iter < max_iters; iter++){
    for (uintV u = k; u < n; u++){
      uintE out_degree = g.vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++)
      {
        uintV v = g.vertices_[u].getOutNeighbor(i);
        tempRank = pr_next[v];
        while (!pr_next[v].compare_exchange_weak(tempRank, pr_next[v] + (pr_curr[u] / out_degree))){}
      }
    }
    my_barrier->wait();
    for (uintV v = k; v < n; v++){
      pr_next[v] = PAGE_RANK(pr_next[v]);
      pr_curr[v].store(pr_next[v]);
      pr_next[v] = 0.0;
    }
    my_barrier->wait();
  }
  time_taken = t1.stop();
}

void pageRankParallelDriver(Graph &g, int max_iters){
  uintV n = g.n_;

  std::atomic<PageRankType> *pr_curr = new std::atomic<PageRankType>[n];
  std::atomic<PageRankType> *pr_next = new std::atomic<PageRankType>[n];
  my_barrier = new CustomBarrier(numberOfThreads);
  std::thread threads[numberOfThreads];

  for (uintV i = 0; i < n; i++)
  {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

  double times[numberOfThreads];
  double time_taken = 0;

  double numThreads = (double) numberOfThreads;
  double numVertexes = (double) n;
  double start = 0;
  double finish = 0;
  uintV startCaster = 0;
  uintV finishCaster = 0;

  timer t2;
  t2.start();
  for (int i = 0; i < numberOfThreads; i++){
    if(((start + numVertexes) / numThreads) < n){
      finish = start + (numVertexes / numThreads);
    }
    else {
      finish = numVertexes;
    }
    startCaster = (uintV) start;
    finishCaster = (uintV) finish;
    threads[i] = std::thread(pageRankParallel, startCaster, finishCaster, max_iters, std::ref(pr_curr), std::ref(pr_next), std::ref(g), std::ref(times[i]));
    start = finish;
  }
  for(int i = 0; i < numberOfThreads; i++){
    threads[i].join();
  }
  time_taken = t2.stop();

  // -------------------------------------------------------------------
  // std::cout << "thread_id, time_taken\n";
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, time_taken
  // 0, 0.12
  // 1, 0.12
  std::cout << "thread_id, time taken" << std::endl;
  for (int i = 0; i < numberOfThreads; i++)
  {
    std::cout << i << " " << times[i] << std::endl;
  }
//
  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++)
  {
    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

void pageRankEdgeDriver(Graph &g, int max_iters){
  uintV n = g.n_;

  std::atomic<PageRankType> *pr_curr = new std::atomic<PageRankType>[n];
  std::atomic<PageRankType> *pr_next = new std::atomic<PageRankType>[n];
  my_barrier = new CustomBarrier(numberOfThreads);
  std::thread threads[numberOfThreads];
  int startVertex[numberOfThreads];
  int finishVertex[numberOfThreads];
  double timers[numberOfThreads];
  uintV verticesCounted[numberOfThreads];
  uintE edgesCounted[numberOfThreads];
  for(int i = 0; i < numberOfThreads; i++){
    timers[i] = 0;
    verticesCounted[i] = 0;
    edgesCounted[i] = 0;
  }

  for (uintV i = 0; i < n; i++)
  {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

  double times[numberOfThreads];
  double time_taken = 0;

  double numThreads = (double) numberOfThreads;
  double numVertexes = (double) n;
  double start = 0;
  double finish = 0;
  uintV startCaster = 0;
  uintV finishCaster = 0;

  double numEdges = g.m_;

  int edgeStart[numberOfThreads];
  int edgeEnd[numberOfThreads];
  uintV vertexStart[numberOfThreads];
  uintV vertexEnd[numberOfThreads];
  double edges = g.m_;

  for(int i = 0; i < numberOfThreads; i++){
    edgeStart[i] = 0;
    edgeEnd[i] = 0;
    vertexStart[i] = 0;
    vertexEnd[i] = 0;
  }

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

  int edgeIterator = 0;
  uintE currentEdge = 0;

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

  for(int i = 0; i < numberOfThreads; i++){
    std::cout << "Thread " << i << " will begin at " << vertexStart[i] << " and end at " << vertexEnd[i] << " covering edges from " << edgeStart[i] << " to " << edgeEnd[i] << std::endl;
  }

  timer t2;
  t2.start();
  for (int i = 0; i < numberOfThreads; i++){
    if(((start + numVertexes) / numThreads) < n){
      finish = start + (numVertexes / numThreads);
    }
    else {
      finish = numVertexes;
    }
    startCaster = (uintV) start;
    finishCaster = (uintV) finish;
    threads[i] = std::thread(pageRankParallel, vertexStart[i], vertexEnd[i], max_iters, std::ref(pr_curr), std::ref(pr_next), std::ref(g), std::ref(times[i]));
    start = finish;
  }
  for(int i = 0; i < numberOfThreads; i++){
    threads[i].join();
  }
  time_taken = t2.stop();

  // -------------------------------------------------------------------
  // std::cout << "thread_id, time_taken\n";
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, time_taken
  // 0, 0.12
  // 1, 0.12
  std::cout << "thread_id, time taken" << std::endl;
  for (int i = 0; i < numberOfThreads; i++)
  {
    std::cout << i << " " << times[i] << std::endl;
  }
//
  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++)
  {
    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

uintV getNextVertexToBeProcessed(Graph &g){
  uintV index;
  if(currentVertex < g.n_){
    index = std::atomic_fetch_add(&currentVertex,granularityInt);
  }else{
    return -1;
  }
  return index;
}

void pageRankDynamic(int max_iters, std::atomic<PageRankType> *pr_curr, std::atomic<PageRankType> *pr_next, Graph &g, double& time_taken, int& verticiesComputed){
  timer t1;
  uintV n = g.n_;
  time_taken = 0.0;
  PageRankType tempRank;
  uintV u = 1;
  int touched[n];
  for(int i = 0; i < n; i++){
    touched[i] = 0;
  }
  uintV end = 0;

  t1.start();
  for (int iter = 0; iter < max_iters; iter++){
    //for each iter, dynamically grab each edge?
    u = getNextVertexToBeProcessed(g);
   // std::cout << "Starting at .. " << u << std::endl;
    while(u >= 0){
      if(u + granularityInt < n){
        end = granularityInt;
       // std::cout << "end outer = " << end << std::endl;
      } else{
        end = n - u;
      //  std::cout << "end = " << end << std::endl;
      }
      for(int i = 0; i < end; i++){
        touched[u] = 1;
        verticiesComputed++;
        uintE out_degree = g.vertices_[u].getOutDegree();
        for (uintE i = 0; i < out_degree; i++){
          uintV v = g.vertices_[u].getOutNeighbor(i);
          tempRank = pr_next[v];
          while (!pr_next[v].compare_exchange_weak(tempRank, pr_next[v] + (pr_curr[u] / out_degree))){}
        }
        u++;
      }
      u = getNextVertexToBeProcessed(g);
  }
    my_barrier->wait();

    for(int i = 0; i < n; i++){
      if(touched[i] == 1){
        pr_next[i] = PAGE_RANK(pr_next[i]);
        pr_curr[i].store(pr_next[i]);
        pr_next[i] = 0.0;
        touched[i] = 0;
      }
    }
    my_barrier->wait();
    currentVertex = 0;
    my_barrier->wait();
//-------------------
  }
  time_taken = t1.stop();
}

void pageRankDynamicDriver(Graph &g, int max_iters){
  uintV n = g.n_;
  iterating = 1;
  currentVertex.store(0);

  std::atomic<PageRankType> *pr_curr = new std::atomic<PageRankType>[n];
  std::atomic<PageRankType> *pr_next = new std::atomic<PageRankType>[n];
  my_barrier = new CustomBarrier(numberOfThreads);
  std::thread threads[numberOfThreads];
  int verticesComputed[numberOfThreads];

  for(int i = 0; i < numberOfThreads; i++){
    verticesComputed[i] = 0;
  }

  for (uintV i = 0; i < n; i++)
  {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

  double times[numberOfThreads];
  double time_taken = 0;

  timer t2;
  t2.start();
  for (int i = 0; i < numberOfThreads; i++){
    threads[i] = std::thread(pageRankDynamic, max_iters, std::ref(pr_curr), std::ref(pr_next), std::ref(g), std::ref(times[i]),std::ref(verticesComputed[i]));
  }
  for(int i = 0; i < numberOfThreads; i++){
    threads[i].join();
  }

  int total_vertices = 0;
  for(int i = 0; i < numberOfThreads; i++){
    std::cout << "The number of vertices computed for thread " << i << " is " << verticesComputed[i] << std::endl;
    total_vertices += verticesComputed[i];
  }
  std::cout << "And the total number of vertices computed was " << total_vertices << std::endl;
  std::cout << "And the total number of vertices in the graph was " << g.n_ << std::endl;

  time_taken = t2.stop();

  // -------------------------------------------------------------------
  // std::cout << "thread_id, time_taken\n";
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, time_taken
  // 0, 0.12
  // 1, 0.12
  std::cout << "thread_id, time taken" << std::endl;
  for (int i = 0; i < numberOfThreads; i++)
  {
    std::cout << i << " " << times[i] << std::endl;
  }
//
  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++)
  {
    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}


void pageRankSerial(Graph &g, int max_iters) {
  uintV n = g.n_;
  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];
  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }
  // Push based pagerank
  timer t1;
  double time_taken = 0.0;
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = 0; u < n; u++) {
      uintE out_degree = g.vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = g.vertices_[u].getOutNeighbor(i);
        pr_next[v] += (pr_curr[u] / out_degree);
      }
    }
    for (uintV v = 0; v < n; v++) {
      pr_next[v] = PAGE_RANK(pr_next[v]);

      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v];
      pr_next[v] = 0.0;
    }
  }
  time_taken = t1.stop();
  // -------------------------------------------------------------------
  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++) {
    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_push",
      "Calculate page_rank using serial and parallel execution");
  options.add_options(
      "",
      {
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
          {"nIterations", "Maximum number of iterations",
           cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
          {"strategy", "Strategy to be used",
           cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
          {"granularity", "Granularity to be used",
           cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  uint strategy = cl_options["strategy"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  granularityInt = 1;
  uint granularity = cl_options["granularity"].as<uint>();
  granularityInt = granularity;
  std::cout << "Granularity int = " << granularityInt << std::endl;
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
//
#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using FLOAT\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";
  numberOfThreads = n_workers;
  std::cout << "Task decomposition strategy : " << strategy << "\n";
  std::cout << "Iterations : " << max_iterations << "\n";
  std::cout << "Granularity : " << granularity << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  switch (strategy) {
  case 0:
    std::cout << "\nSerial\n";
    pageRankSerial(g, max_iterations);
    break;
  case 1:
    std::cout << "\nVertex-based work partitioning\n";
    pageRankParallelDriver(g, max_iterations);
    break;
  case 2:
    std::cout << "\nEdge-based work partitioning\n";
    pageRankEdgeDriver(g, max_iterations);
    break;
  case 3:
    std::cout << "\nDynamic task mapping\n";
    pageRankDynamicDriver(g, max_iterations);
    break;
  default:
    break;
  }

  return 0;
}
