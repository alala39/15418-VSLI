/**
 * Parallel VLSI Wire Routing via OpenMP
 * Aditya Lala(alala), Name 2(andrew_id 2)
 */

#include "wireroute.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <omp.h>
#include <algorithm>
#include "mic.h"

#define BUFSIZE 1024

static int _argc;
static const char **_argv;

/* Starter code function, don't touch */
const char *get_option_string(const char *option_name,
			      const char *default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return _argv[i + 1];
  return default_value;
}

/* Starter code function, do not touch */
int get_option_int(const char *option_name, int default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return atoi(_argv[i + 1]);
  return default_value;
}

/* Starter code function, do not touch */
float get_option_float(const char *option_name, float default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return (float)atof(_argv[i + 1]);
  return default_value;
}

/* Starter code function, do not touch */
static void show_help(const char *program_path)
{
    printf("Usage: %s OPTIONS\n", program_path);
    printf("\n");
    printf("OPTIONS:\n");
    printf("\t-f <input_filename> (required)\n");
    printf("\t-n <num_of_threads> (required)\n");
    printf("\t-p <SA_prob>\n");
    printf("\t-i <SA_iters>\n");
}

//TODO NEED TO FIRST REMOVE THE CURRENT WIRE BEFORE FINDING THE COST FOR THE MOST OPTIMAL PATH 

void removeEqualXWire(int x,int y1,int y2,cost_t *costs,int dim_x,int dim_y)
{
    int start = std::min(y1,y2);
    int end = std::max(y1,y2);
    while(start<=end)
    {
        if(x==29 && start==9)
        printf("BeforeRX: %d\n",costs[start * dim_x + x]);
        if(costs[start * dim_x + x])
        costs[start * dim_x + x]-=1;
        if(x==29 && start==9)
        printf("AfterX: %d\n",costs[start * dim_x + x]);
        start+=1;
    }
}

void removeEqualYWire(int x1,int x2,int y,cost_t *costs,int dim_x,int dim_y)
{
    int start = std::min(x1,x2);
    int end = std::max(x1,x2);
    while(start <= end)
    {
        if(start==29 && y==9)
        printf("BeforeRY: %d\n",costs[y * dim_x + start]);
        if(costs[y * dim_x + start])
        costs[y * dim_x + start]-=1;
        if(start==29 && y==9)
        printf("AfterRY: %d\n",costs[y * dim_x + start]);
        start+=1;
    }
}

void placeXWire(int x,int y1,int y2,cost_t *costs,int dim_x,int dim_y)
{
    int start = std::min(y1,y2);
    int end = std::max(y1,y2);
    while(start<=end)
    {
        if(x==29 && start==9)
        printf("BeforePX: %d\n",costs[start * dim_x + x]);
        costs[start * dim_x + x]+=1;
        if(x==29 && start==9)
        printf("AfterPX: %d\n",costs[start * dim_x + x]);
        start+=1;
    }
}

void placeYWire(int x1,int x2,int y,cost_t *costs,int dim_x,int dim_y)
{
    int start = std::min(x1,x2);
    int end = std::max(x1,x2);
    while(start <= end)
    {
        if(start==29 && y==9)
        printf("BeforePY: %d\n",costs[y * dim_x + start]);
        costs[y * dim_x + start]+=1;
        if(start==29 && y==9)
        printf("AfterPY: %d\n",costs[y * dim_x + start]);
        start+=1;
    }
}

void placeWire(wire_t *wire,cost_t *costs,int dim_x,int dim_y)
{
    //Assumes that the wire being placed has at least one bend
    int x1 = wire->x1;
    int x2 = wire->x2;
    int y1 = wire->y1;
    int y2 = wire->y2;
    int bend1x = wire->bend1x;
    int bend2x = wire->bend2x;
    int bend1y = wire->bend1y;
    int bend2y = wire->bend2y;
    bool BADWIRE = (x1==16) && (y1==9) && (bend1x==29) && (bend1y==9) && (bend2x==29) && (bend2y==23) && (x2==30) && (y2==23);
    if(bend1x!=undef && bend1y!=undef)  //At least one bend
    {
        if(x1==bend1x)  //Segment 1
        {
            placeXWire(x1,y1,bend1y,costs,dim_x,dim_y);
        }
        else
        {
            assert(y1==bend1y);
            placeYWire(x1,bend1x,y1,costs,dim_x,dim_y);
        }
        if(bend2x!=undef && bend2y!=undef)  //Two bends
        {
            if(bend1x==bend2x)  //Segment 2
            {
                placeXWire(bend1x,bend1y,bend2y,costs,dim_x,dim_y);
            }
            else
            {
                assert(bend1y==bend2y && bend1x!=bend2x);
                placeYWire(bend1x,bend2x,bend2y,costs,dim_x,dim_y);
            }
            if(bend2x==x2)  //Segment 3
            {
                placeXWire(x2,bend2y,y2,costs,dim_x,dim_y);
            }
            else
            {
                assert(bend2y==y2);
                placeYWire(bend2x,x2,y2,costs,dim_x,dim_y);
            }
        }
        else //Only one bend
        {
            assert(bend2x==undef && bend2y==undef);
            if(bend1x==x2)
            {
                placeXWire(x2,bend1y,y2,costs,dim_x,dim_y);
            }
            else
            {
                assert(bend1y==y2);
                placeYWire(bend1x,x2,y2,costs,dim_x,dim_y);
            }
        }
    }
    else  //No bend
    {
        assert(bend1x==undef && bend1y==undef && bend2x==undef && bend2y==undef);
        if(x1==x2)
        {
            placeXWire(x1,y1,y2,costs,dim_x,dim_y);
        }
        else
        {
            assert(y1==y2);
            placeYWire(x1,x2,y1,costs,dim_x,dim_y);
        }
    }
    if(bend1x!=undef && bend1y!=undef)
    {
        costs[bend1y * dim_x + bend1x]-=1;
    }
    if(bend2x!=undef && bend2y!=undef)
    {
        costs[bend2y * dim_x + bend2x]-=1;
    }
}

//functions that get the cost along the path of a given wire
int costPerLine(int x1,int x2,int y1,int y2,cost_t *costs, int dim_x, int dim_y)
{
    int result = 0;
    if(x1==x2)
    {
        int start = std::min(y1,y2);
        int end = std::max(y1,y2);
        while(start <= end)
        {
            //result += costs[start * dim_x + x1];
            if(costs[start * dim_x + x1]>result)
            result = costs[start * dim_x + x1];
            start+=1;
        }
    }
    else
    {
        assert(y1==y2);
        int start = std::min(x1,x2);
        int end = std::max(x1,x2);
        while(start <= end)
        {
            //result += costs[y1 * dim_x + start];
            if(costs[y1 * dim_x + start])
            result = costs[y1 * dim_x + start];
            start+=1;
        }
    }
    return result;
}

int costEqualX(int x, int y1,int y2,cost_t *costs,int dim_x,int dim_y,int delta,wire_t *wire,int *deltaUsed)
{
    int result = 0;
    int start = std::min(y1,y2);
    int end = std::max(y1,y2);
    result = costPerLine(x,x,y1,y2,costs,dim_x,dim_y);
    wire_t w = {x,x,y1,y2,undef,undef,undef,undef};
    *wire = w;
    if(y1==y2)
    {
        return 0;
    }
    int curMin = -1;
    if(delta>0)
    {
        for(int i = 1;i<=delta/2;i++)
        {
            if(x+i >= dim_x || x-i < 0)
            {
                continue;
            }
            int res1 = std::max(costPerLine(x,x+i,start,start,costs,dim_x,dim_y),std::max(costPerLine(x+i,x+i,start,end,costs,dim_x,dim_y),costPerLine(x+i,x,end,end,costs,dim_x,dim_y)));
            int res2 = std::max(costPerLine(x,x-i,start,start,costs,dim_x,dim_y),std::max(costPerLine(x-i,x-i,start,end,costs,dim_x,dim_y),costPerLine(x-i,x,end,end,costs,dim_x,dim_y)));
            result = std::min(result,std::min(res1,res2));
            if(result==res1 && (result<=curMin || curMin==-1 || *deltaUsed==undef))
            {
                wire_t w = {x,x,y1,y2,x+i,x+i,y1,y2};
                *wire = w;
                curMin = res1;
                *deltaUsed = i;
            }
            else if(result==res2 && (result<=curMin || curMin==-1 || *deltaUsed==undef))
            {
                wire_t w = {x,x,y1,y2,x-i,x-i,y1,y2};
                *wire = w;
                curMin = res2;
                *deltaUsed = i;
            }
            else  //Else case is that the straight result is the min
            {
                curMin = result;
                *deltaUsed = 0; //No delta is used in this case
            }
        }
    }
    return result;
}

int checkWire(int x1,int x2,int y1,int y2,cost_t *costs,int dim_x,int dim_y,wire_t *wire)
{
    //wire_t w(x1,x2,y1,y2);
    //*wire = w;
    return std::min(costPerLine(x1,x1,y1,y2,costs,dim_x,dim_y)+costPerLine(x2,x2,y1,y2,costs,dim_x,dim_y),costPerLine(x1,x2,y1,y1,costs,dim_x,dim_y)+costPerLine(x1,x2,y2,y2,costs,dim_x,dim_y));
}

/////////////////////////////////////
int checkAllHorizontal(int x1,int x2,int y1,int y2,cost_t *costs,int dim_x,int dim_y,int delta,wire_t *wire)
{
    //TODO THIS FUNCTION WILL CHECK ALL HORIZONTAL (CHANGES IN X) VARIATIONS OF THE WIRE CONFIGURATION
    int result = 0;
    int startX = std::min(x1,x2);
    int endX = std::max(x1,x2);
    int lowerBound = startX - (delta/2);
    int upperBound = endX + (delta/2);
    int startY;
    int endY;
    int curMin = -1;
    if(startX==x1)
    {
        startY = y1;
        endY = y2;
    }
    else
    {
        assert(startX==x2);
        startY = y2;
        endY = y1;
    }
    for(int i = lowerBound;i <=upperBound;i++)
    {
        if(i >= dim_x || i < 0)
        {
            continue;
        }
        int rec = std::max(costPerLine(startX,i,startY,startY,costs,dim_x,dim_y),std::max(costPerLine(i,i,startY,endY,costs,dim_x,dim_y),costPerLine(i,endX,endY,endY,costs,dim_x,dim_y)));
        if(rec<=curMin || curMin==-1)
        {
            if(i!=startX && i!=endX)
            {
                curMin = rec;
                //TODO STORE THE NEW WIRE
                wire_t w = {startX,endX,startY,endY,i,i,startY,endY};
                *wire = w;
            }
            else if(i==startX)
            {
                curMin = rec;
                //TODO STORE THE NEW WIRE
                wire_t w = {startX,endX,startY,endY,i,undef,endY,undef};
                *wire = w;
            } 
            else if(i==endX)
            {
                curMin = rec;
                //TODO STORE THE NEW WIRE
                wire_t w = {startX,endX,startY,endY,i,undef,startY,undef};
                *wire = w;
            }
        }
    }
    return curMin;
}
/////////////////////////////////////

int checkAllVertical(int x1,int x2,int y1,int y2,cost_t *costs,int dim_x,int dim_y,int delta,wire_t *wire)
{
    int result = 0;
    int startY = std::min(y1,y2);
    int endY = std::max(y1,y2);
    int lowerBound = startY - (delta/2);
    int upperBound = endY + (delta/2);
    int startX;
    int endX;
    int curMin = -1;
    if(startY==y1)
    {
        startX = x1;
        endX = x2;
    }
    else
    {
        assert(startY==y2);
        startX = x2;
        endX = x1;
    }
    for(int i = lowerBound;i<upperBound;i++)
    {
        if(i>=dim_y || i<0)
        {
            continue;
        }
        int rec = std::max(costPerLine(startX,startX,startY,i,costs,dim_x,dim_y),std::max(costPerLine(startX,endX,i,i,costs,dim_x,dim_y),costPerLine(endX,endX,i,endY,costs,dim_x,dim_y)));
        if(rec<=curMin || curMin==-1)
        {
            if(i!=startY && i !=endY)
            {
                curMin = rec;
                wire_t w = {startX,endX,startY,endY,startX,endX,i,i};
                *wire = w;
            }
            else if(i == startY)
            {
                curMin = rec;
                wire_t w = {startX,endX,startY,endY,endX,undef,i,undef};
                *wire = w;
            }
            else if(i == endY)
            {
                curMin = rec;
                wire_t w = {startX,endX,startY,endY,startX,undef,i,undef};
                *wire = w;
            } 
        }
    }
    return curMin;
}
////////////////////////////////////

int costEqualY(int x1, int x2,int y,cost_t *costs,int dim_x,int dim_y,int delta,wire_t *wire,int *deltaUsed)
{
    int result = 0;
    int start = std::min(x1,x2);
    int end = std::max(x1,x2);
    int curMin = -1;
    result = costPerLine(x1,x2,y,y,costs,dim_x,dim_y);
    wire_t w = {x1,x2,y,y,undef,undef,undef,undef};
    *wire = w;
    if(x1==x2)
    {
        return 0;
    }
    if(delta>0)
    {
        for(int i = 1;i<=delta/2;i++)
        {
            if(y+i>=dim_y || y-i < 0)  //Boundary check
            {
                continue;
            }
            int res1 = std::max(costPerLine(start,start,y,y+i,costs,dim_x,dim_y),std::max(costPerLine(start,end,y+i,y+i,costs,dim_x,dim_y),costPerLine(end,end,y+i,y,costs,dim_x,dim_y)));
            int res2 = std::max(costPerLine(start,start,y,y-i,costs,dim_x,dim_y),std::max(costPerLine(start,end,y-i,y-i,costs,dim_x,dim_y),costPerLine(end,end,y-i,y,costs,dim_x,dim_y)));
            result = std::min(result,std::min(res1,res2));
            if(result==res1 && (res1<=curMin || curMin==-1 || *deltaUsed==undef))
            {
                wire_t w = {x1,x2,y,y,x1,x2,y+i,y+i};
                *wire = w;
                curMin = res1;
                *deltaUsed = i;
            }
            else if(result==res2 && (res1<=curMin || curMin==-1 || *deltaUsed==undef))
            {
                wire_t w = {x1,x2,y,y,x1,x2,y-i,y-i};
                *wire = w;
                curMin = res2;
                *deltaUsed = i;
            }
            else //Case where straight line is the min
            {
                curMin = result;
                *deltaUsed = 0;
            }
        }
    }
    return result;
} 

void serialAlg(wire_t *wires,cost_t *costs,int SA_iters,int dim_x,int dim_y,int num_of_wires,int delta)
{
  int counter = 0;
  while(counter < SA_iters)
  {
    for(int i = 0;i<num_of_wires;i++)
    {
      int x1 = wires[i].x1;
      int y1 = wires[i].y1;
      int x2 = wires[i].x2;
      int y2 = wires[i].y2;
      int bend1x = wires[i].bend1x;
      int bend1y = wires[i].bend1y;
      int bend2x = wires[i].bend2x;
      int bend2y = wires[i].bend2y;
      bool BADWIRE = (x1==16) && (y1==9) && (bend1x==29) && (bend1y==9) && (bend2x==29) && (bend2y==23) && (x2==30) && (y2==23);
      int curMin = 0;
      int minimum = 0;
      //First check to see if points are in a straight line.
      if(x1==x2 || y1==y2)  //Points are in a straight line
      {
        //TODO NEED TO CHECK TO SEE THAT THE MAX OF THE COSTS ALONG THIS WIRE ARE 
        //LESS THAN THE CURRENT MIN COST. JUST CHECKING DISTANCE IS NOT ENOUGH
        //curMin[i] = cost;
        //int curCostWire = costPerLine(x1,x2,y1,y2,costs,dim_x,dim_y);
        if(x1==x2)
        {
            removeEqualXWire(x1,y1,y2,costs,dim_x,dim_y);
            wire_t EqX = {x1,x2,y1,y2,undef,undef,undef,undef};
            int deltaUsed = undef;
            int curCostWire = costEqualX(x1,y1,y2,costs,dim_x,dim_y,delta,&EqX,&deltaUsed);
            placeWire(&EqX,costs,dim_x,dim_y);
            //TODO NEED TO PLACE EQX
        }
        else
        {
            assert(y1==y2);
            removeEqualYWire(x1,x2,y1,costs,dim_x,dim_y);
            wire_t EqY = {x1,x2,y1,y2,undef,undef,undef,undef};
            int deltaUsed = undef;
            int curCostWire = costEqualY(x1,x2,y1,costs,dim_x,dim_y,delta,&EqY,&deltaUsed);
            placeWire(&EqY,costs,dim_x,dim_y);
            //TODO NEED TO PLACE EQY
        }
      }
      else  //Have a bend
      {
        //TODO FIRST CHECK HORIZONTAL PATHS, THEN VERTICAL PATHS FOR NEW MIN.
        if(x1==wires[i].bend1x)
        {
            removeEqualXWire(x1,y1,wires[i].bend1y,costs,dim_x,dim_y);
        }
        else if(y1==wires[i].bend1y)
        { 
            //assert(y1==wires[i].bend1y);
            removeEqualYWire(x1,wires[i].bend1x,y1,costs,dim_x,dim_y);
        }
        if(wires[i].bend2x!=undef && wires[i].bend2y!=undef)
        {
            if(wires[i].bend1x==wires[i].bend2x)
            {  
                removeEqualXWire(wires[i].bend1x,wires[i].bend1y,wires[i].bend2y,costs,dim_x,dim_y);
            }
            else if(wires[i].bend1y==wires[i].bend2y)
            {
                //assert(wires[i].bend1y==wires[i].bend2y);
                removeEqualYWire(wires[i].bend1x,wires[i].bend2x,wires[i].bend1y,costs,dim_x,dim_y);
            }
            if(wires[i].bend2x==x2)
            {
                removeEqualXWire(x2,wires[i].bend2y,y2,costs,dim_x,dim_y);
            }
            else if(wires[i].bend2y==y2)
            { 
                //assert(wires[i].bend2y==y2);
                removeEqualYWire(wires[i].bend2x,x2,y2,costs,dim_x,dim_y);
            }
        }
        else if(wires[i].bend1x==x2)
        {
            removeEqualXWire(x2,wires[i].bend1y,y2,costs,dim_x,dim_y);
        }
        else if(wires[i].bend1y==y2)
        {
            //assert(wires[i].bend1y==y2);
            removeEqualYWire(wires[i].bend1x,x2,y2,costs,dim_x,dim_y);
        }
        //wire_t noB = {x1,x2,y1,y2,undef,undef,undef,undef};
        wire_t H = {x1,x2,y1,y2,undef,undef,undef,undef};
        wire_t V = {x1,x2,y1,y2,undef,undef,undef,undef};
        //int noBend = checkWire(x1,x2,y1,y2,costs,dim_x,dim_y,&noB); 
        int Horizontal = checkAllHorizontal(x1,x2,y1,y2,costs,dim_x,dim_y,delta,&H);
        int Vertical = checkAllVertical(x1,x2,y1,y2,costs,dim_x,dim_y,delta,&V);
        minimum = std::min(Horizontal,Vertical);
        /*if(minimum==noBend)
        {
            //TODO NEED TO PLACE NOB
            wires[i] = noB;
            placeWire(&noB,costs,dim_x,dim_y);
        }*/
        if(minimum==Horizontal)
        {
            //TODO NEED TO PLACE H
            wires[i] = H;
            placeWire(&H,costs,dim_x,dim_y);
        }
        else
        {
            assert(minimum==Vertical);
            //TODO NEED TO PLACE V
            wires[i] = V;
            placeWire(&V,costs,dim_x,dim_y);
        }
      }
    }
    counter++;
  }
}


int main(int argc, const char *argv[])
{
  using namespace std::chrono;
  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::duration<double> dsec;

  auto init_start = Clock::now();
  double init_time = 0;
 
  _argc = argc - 1;
  _argv = argv + 1;

  /* You'll want to use these parameters in your algorithm */
  const char *input_filename = get_option_string("-f", NULL);
  int num_of_threads = get_option_int("-n", 1);
  double SA_prob = get_option_float("-p", 0.1f);
  int SA_iters = get_option_int("-i", 5);

  int error = 0;

  if (input_filename == NULL) {
    printf("Error: You need to specify -f.\n");
    error = 1;
  }

  if (error) {
    show_help(argv[0]);
    return 1;
  }
  
  printf("Number of threads: %d\n", num_of_threads);
  printf("Probability parameter for simulated annealing: %lf.\n", SA_prob);
  printf("Number of simulated anneling iterations: %d\n", SA_iters);
  printf("Input file: %s\n", input_filename);

  FILE *input = fopen(input_filename, "r");

  if (!input) {
    printf("Unable to open file: %s.\n", input_filename);
    return -1;
  }
 
  int dim_x, dim_y;
  int delta;
  int num_of_wires;

  fscanf(input, "%d %d\n", &dim_x, &dim_y);
  fscanf(input, "%d\n", &delta);
  assert(delta >= 0 && delta%2 == 0);
  fscanf(input, "%d\n", &num_of_wires);

  wire_t *wires = (wire_t *)calloc(num_of_wires, sizeof(wire_t));
  (void)wires;
 
  int x1, y1, x2, y2;
  int index = 0;
  while (fscanf(input, "%d %d %d %d\n", &x1, &y1, &x2, &y2) != EOF) {
    /* PARSE THE INPUT FILE HERE.
     * Define wire_t in wireroute.h and store 
     * x1, x2, y1, and y2 into the wires array allocated above
     * based on your wire_t definition. */
     wires[index].x1 = x1;
     wires[index].x2 = x2;
     wires[index].y1 = y1;
     wires[index].y2 = y2;
     wires[index].bend1x = undef;
     wires[index].bend2x = undef;
     wires[index].bend1y = undef;
     wires[index].bend2y = undef;
     index++;
  }

  if (index != num_of_wires) {
    printf("Error: wire count mismatch");
    return -1;
  }

  cost_t *costs = (cost_t *)calloc(dim_x * dim_y, sizeof(cost_t));
  (void)costs;
  /* INITIALIZE YOUR COST MATRIX HERE */

  /*for(size_t x = 0;x<dim_x;x++)  //No wires have been placed yet
  {
    for(size_t y = 0;y<dim_y;y++)
    {
      costs[x][y] = 0;
    }
  }*/  //TODO DO NOT NEED THIS AS THE COST ARRAY IS INIT WITH CALLOC

  
  /* Initialize additional data structures needed in the algorithm 
   * here if you feel it's needed. */

  error = 0;

  init_time += duration_cast<dsec>(Clock::now() - init_start).count();
  printf("Initialization Time: %lf.\n", init_time);

  auto compute_start = Clock::now();
  double compute_time = 0;
#ifdef RUN_MIC /* Use RUN_MIC to distinguish between the target of compilation */

  /* This pragma means we want the code in the following block be executed in 
   * Xeon Phi.
   */
#pragma offload target(mic) \
  inout(wires: length(num_of_wires) INOUT)    \
  inout(costs: length(dim_x*dim_y) INOUT)
#endif
  {
    /* Implement the wire routing algorithm here
     * Feel free to structure the algorithm into different functions
     * Don't use global variables.
     * Use OpenMP to parallelize the algorithm. 
     * You should really implement as much of this (if not all of it) in
     * helper functions. */
      //TODO WILL IMPLEMENT ANNEALING HERE
      serialAlg(wires,costs,SA_iters,dim_x,dim_y,num_of_wires,delta);
  }

  compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
  printf("Computation Time: %lf.\n", compute_time);

  /* OUTPUT YOUR RESULTS TO FILES HERE 
   * When you're ready to output your data to files, uncommment this chunk of
   * code and fill in the specified blanks indicated by comments. More about
   * this in the README. */
  char input_filename_cpy[BUFSIZE];
  strcpy(input_filename_cpy, input_filename);
  char *filename = basename(input_filename_cpy);
  char output_filename[BUFSIZE];


  sprintf(output_filename, "costs_%s_%d.txt", filename, num_of_threads);
  FILE *output_costs_file = fopen(output_filename, "w");
  if (!output_costs_file) {
    printf("Error: couldn't output costs file");
    return -1;
  }

  fprintf(output_costs_file, "%d %d\n", dim_x, dim_y);
  for(int row = 0;row<dim_y;row++)
  {
      for(int col = 0;col<dim_x;col++)
      {
          fprintf(output_costs_file, "%d ",costs[row * dim_x + col]);
      }
      fprintf(output_costs_file,"\n");
  }
  
  // WRITE COSTS TO FILE HERE 

  fclose(output_costs_file);


  sprintf(output_filename, "output_%s_%d.txt", filename, num_of_threads);
  FILE *output_routes_file = fopen(output_filename, "w");
  if (!output_routes_file) {
    printf("Error: couldn't output routes file");
    return -1;
  }

  fprintf(output_routes_file, "%d %d\n", dim_x, dim_y);
  fprintf(output_routes_file, "%d\n",delta);
  fprintf(output_routes_file, "%d\n", num_of_wires);

  for(int i = 0;i<num_of_wires;i++)
  {
      if(wires[i].bend1x!=undef && wires[i].bend1y!=undef)
      {
          if(wires[i].bend2x!=undef && wires[i].bend2y!=undef)
          {
              fprintf(output_routes_file, "%d %d %d %d %d %d %d %d\n",wires[i].x1,wires[i].y1,wires[i].bend1x,wires[i].bend1y,wires[i].bend2x,wires[i].bend2y,wires[i].x2,wires[i].y2);
          }
          else
          {
              fprintf(output_routes_file, "%d %d %d %d %d %d\n",wires[i].x1,wires[i].y1,wires[i].bend1x,wires[i].bend1y,wires[i].x2,wires[i].y2);
          }
      }
      else
      {
          fprintf(output_routes_file, "%d %d %d %d\n",wires[i].x1,wires[i].y1,wires[i].x2,wires[i].y2);
      }
  }
  
  // WRITE WIRES TO FILE HERE

  fclose(output_routes_file);

  return 0;
}
