#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

class MeasureTime {

  long mtime, seconds, useconds;
  struct timeval start, end;

public:
  MeasureTime(){

  }

  void startMeasuring(){
    gettimeofday(&start, NULL);
  }

  long  endMeasuring(){
    gettimeofday(&end, NULL);
    
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

    // printf("Elapsed time: %ld milliseconds\n", mtime);
    return mtime;
  }


};

/*int main()
{
  MeasureTime mt;
  mt.startMeasuring();
  mt.endMeasuring();
  
    return 0;
    }*/


