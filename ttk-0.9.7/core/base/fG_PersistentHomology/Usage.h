#ifndef USAGE_H
#define USAGE_H

#include "sys/types.h"
#include "sys/resource.h"


#ifdef __APPLE__
    #define UNIT 1
#else
    #define UNIT 0
#endif

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

class MemoryUsage{

private:
    int who = RUSAGE_SELF;
    struct rusage usage;
    int ret;


public:

    inline float getValue_in_KB(bool output){

         ret=getrusage(who,&usage);
         if(output) cout << "Memory Usage: " << usage.ru_maxrss/pow(1024.0,UNIT) << " KB" << endl;

         return usage.ru_maxrss/(1024.0);
    }

    inline float getValue_in_MB(bool output){

         ret=getrusage(who,&usage);
         if(output) cout << "Memory Usage: " << usage.ru_maxrss/pow(1024.0,UNIT+1) << " MB" << endl;

         return usage.ru_maxrss/(1024.0*1024.0);
    }

    inline float getValue_in_GB(bool output){

         ret=getrusage(who,&usage);
         if(output) cout << "Memory Usage: " << usage.ru_maxrss/pow(1024.0,UNIT+2) << " GB" << endl;

         return usage.ru_maxrss/(1024.0*1024.0*1024.0);
    }
};



#endif // USAGE_H
