#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sharedspice.h>
#include <sstream>
#include <vector>

SUITE(Temperature) {
    bool no_bg = true;
    int vecgetnumber = 0;
    double v2dat;
    static bool has_break = false;
    int testnumber = 0;
    void alterp(int sig);
    static bool errorflag = false;

#ifndef _MSC_VER
    pthread_t mainthread;
#endif // _MSC_VER

    int ng_getchar(char* outputreturn, int ident, void* userdata);

    int ng_getstat(char* outputreturn, int ident, void* userdata);

    int ng_exit(int, bool, bool, int ident, void*);

    int ng_thread_runs(bool noruns, int ident, void* userdata);

    int ng_initdata(pvecinfoall intdata, int ident, void* userdata);

    int ng_data(pvecvaluesall vdata, int numvecs, int ident, void* userdata);

    int cieq(char* p, char* s);

    int ciprefix(const char* p, const char* s);

    /* Callback function called from bg thread in ngspice to transfer
    any string created by printf or puts. Output to stdout in ngspice is
    preceded by token stdout, same with stderr.*/
    int ng_getchar(char* outputreturn, int ident, void* userdata) {
        printf("%s\n", outputreturn);
        return 0;
    }

    int ng_getstat(char* outputreturn, int ident, void* userdata) {
        printf("%s\n", outputreturn);
        return 0;
    }

    int ng_thread_runs(bool noruns, int ident, void* userdata) {
        no_bg = noruns;
        if (noruns)
            printf("bg not running\n");
        else
            printf("bg running\n");

        return 0;
    }

    /* Callback function called from bg thread in ngspice once per accepted data point */
    int ng_data(pvecvaluesall vdata, int numvecs, int ident, void* userdata) {
        int* ret;

        v2dat = vdata->vecsa[vecgetnumber]->creal;
        if (!has_break && (v2dat > 0.5)) {
            //        printf("Data V(2) value: %f\n", v2dat);
            /* using signal SIGTERM by sending to main thread, alterp() then is run from the main thread,
              (not on Windows though!)  */
#ifndef _MSC_VER
            if (testnumber == 4)
                pthread_kill(mainthread, SIGTERM);
#endif
            has_break = true;
            printf("Pause requested, setpoint reached\n");
            /* leave bg thread for a while to allow halting it from main */
#if defined(__MINGW32__) || defined(_MSC_VER)
            Sleep(100);
#else
            usleep(100000);
#endif
            //        ret = ((int * (*)(char*)) ngSpice_Command_handle)("bg_halt");
        }
        return 0;
    }

    /* Callback function called from bg thread in ngspice once upon intialization
       of the simulation vectors)*/
    int ng_initdata(pvecinfoall intdata, int ident, void* userdata) {
        int i;
        int vn = intdata->veccount;
        for (i = 0; i < vn; i++) {
            printf("Vector: %s\n", intdata->vecs[i]->vecname);
            /* find the location of V(2) */
            if (cieq(intdata->vecs[i]->vecname, "V(2)"))
                vecgetnumber = i;
        }
        return 0;
    }

    /* Callback function called from bg thread in ngspice if fcn controlled_exit()
       is hit. Do not exit, but unload ngspice. */
    int ng_exit(int exitstatus, bool immediate, bool quitexit, int ident, void* userdata) {
        if (quitexit) {
            printf("DNote: Returned form quit with exit status %d\n", exitstatus);
            exit(exitstatus);
        }
        if (immediate) {
            printf("DNote: Unloading ngspice inmmediately is not possible\n");
            printf("DNote: Can we recover?\n");
        }

        else {
            printf("DNote: Unloading ngspice is not possible\n");
            printf("DNote: Can we recover? Send 'quit' command to ngspice.\n");
            errorflag = true;
            ngSpice_Command("quit 5");
            //        raise(SIGINT);
        }

        return exitstatus;
    }

    /* Funcion called from main thread upon receiving signal SIGTERM */
    void alterp(int sig) {
        ngSpice_Command("bg_halt");
    }

    /* Case insensitive str eq. */
    /* Like strcasecmp( ) XXX */
    int cieq(char* p, char* s) {
        while (*p) {
            if ((isupper(*p) ? tolower(*p) : *p) != (isupper(*s) ? tolower(*s) : *s))
                return (false);
            p++;
            s++;
        }
        return (*s ? false : true);
    }

    /* Case insensitive prefix. */
    int ciprefix(const char* p, const char* s) {
        while (*p) {
            if ((isupper(*p) ? tolower(*p) : *p) != (isupper(*s) ? tolower(*s) : *s))
                return (false);
            p++;
            s++;
        }
        return (true);
    }
    // TEST(Test_Ngspice){

    //     int ret, i;
    //     char *curplot, *vecname;
    //     char ** circarray;
    //     char **vecarray;

    //     ret = ngSpice_Init(ng_getchar, ng_getstat, ng_exit,  ng_data, ng_initdata, ng_thread_runs, NULL);
    //     testnumber = 2;

    //     /* create a circuit that fails due to missing include */
    //     ret = ngSpice_Command("circbyline fail test");
    //     ret = ngSpice_Command("circbyline V1 1 0 1");
    //     ret = ngSpice_Command("circbyline R1 1 0 1");
    //     ret = ngSpice_Command("circbyline .dc V1 0 1 0.1");
    //     ret = ngSpice_Command("circbyline .end");
    //     /* wait to catch error signal, if available */
    // #if defined(__MINGW32__) || defined(_MSC_VER)
    //     Sleep(100);
    // #else
    //     usleep(100000);
    // #endif

    //     ret = ngSpice_Command("bg_run");
    //     ret = ngSpice_Command("bg_halt");
    //     ret = ngSpice_Command("listing");
    //     ret = ngSpice_Command("alter R1=2");
    //     ret = ngSpice_Command("bg_resume");
    //     ret = ngSpice_Command("write test3.raw V(2)");
    //     printf("rawfile testout3.raw created\n");
    // }
}